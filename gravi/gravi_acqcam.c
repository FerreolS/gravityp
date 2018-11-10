/* $Id: gravi_vis.c,v 1.10 2014/11/12 15:10:40 nazouaoui Exp $
 *
 * This file is part of the GRAVI Pipeline
 * Copyright (C) 2002,2003 European Southern Observatory
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/**
 * @defgroup gravi_acqcam  Acquisition Camera
 *
 * This module process the data of the acquisition camera. The main function
 * called by the recipe gravity_vis is @c gravi_reduce_acqcam(). Two kind of acquisition
 * camera images are reduced as described in section Algorithms/Processing ACQ :
 * - the field images are reduced by the function @c gravi_acqcam_field()
 * - the pupil images are reduced by the function @c gravi_acqcam_pupil()
 *
 */
/**@{*/

/*-----------------------------------------------------------------------------
                                   Includes
 -----------------------------------------------------------------------------*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cpl.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include <hdrl_strehl.h>

#include "gravi_dfs.h"
#include "gravi_data.h"
#include "gravi_pfits.h"

#include "gravi_cpl.h"
#include "gravi_utils.h"

#include "gravi_acqcam.h"

/*-----------------------------------------------------------------------------
                               Private prototypes
 -----------------------------------------------------------------------------*/

double exp1 (double x);
double sin1 (double x);
int gravi_acqcam_xy_diode (const double v[], double *xd, double *yd);

cpl_error_code gravi_acqcam_get_pup_ref (cpl_propertylist * header,
                                         cpl_size tel,
                                         cpl_vector * pupref);

cpl_error_code gravi_acqcam_get_diode_ref (cpl_propertylist * header,
                                           cpl_size tel,
                                           cpl_vector * output);

static int gravi_acqcam_spot (const double x_in[], const double v[], double *result);
static int gravi_acqcam_xy_sub (const double v[], double *xsub, double *ysub);

cpl_error_code gravi_acqcam_spot_imprint (cpl_image * img, cpl_vector * a);

static int gravi_acqcam_spot_dfda (const double x_in[], const double v[], double result[]);

cpl_error_code gravi_acqcam_fit_spot (cpl_image * img, cpl_size ntry,
                                      cpl_vector * a, 
                                      int fitAll,
                                      int * nspot);

double gravi_acqcam_z2meter (double PositionPixels);

cpl_error_code gravi_acqcam_pupil (cpl_image * mean_img,
                                   cpl_imagelist * acqcam_imglist,
                                   cpl_propertylist * header,
                                   cpl_table * acqcam_table,
                                   cpl_propertylist * o_header);

cpl_error_code gravi_acqcam_field (cpl_image * mean_img,
                                   cpl_imagelist * acqcam_imglist,
                                   cpl_propertylist * header,
                                   cpl_table * acqcam_table,
                                   cpl_propertylist * o_header);

cpl_error_code gravi_acq_fit_gaussian (cpl_image * img, double *x, double *y,
                                       double *ex, double *ey, cpl_size size);

cpl_error_code gravi_acq_measure_strehl(cpl_image * img, double x, double y, 
                                        double pscale, double *SR, cpl_propertylist * header);

cpl_error_code gravi_acq_measure_max(cpl_image * img, double x, double y, double size, double * img_max);

/* This global variable optimises the computation
 * of partial derivative on fitted parameters */
const extern int * GRAVI_LVMQ_FREE;
const int * GRAVI_LVMQ_FREE = NULL;

/* Number of parameters in the model 'gravi_acqcam_spot' 
 * And position of parameters */
#define GRAVI_SPOT_NA    30
#define GRAVI_SPOT_SUB   0
#define GRAVI_SPOT_ANGLE 8
#define GRAVI_SPOT_SCALE 9
#define GRAVI_SPOT_DIODE 10
#define GRAVI_SPOT_FWHM  13
#define GRAVI_SPOT_FLUX  14

#define GRAVI_ACQ_PUP_FLUX 1e6

/*-----------------------------------------------------------------------------
                             Functions code
 -----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
 * @brief Preprocess the ACQ images: correct bad pixels, clean from
 *        pupil background via blinking, filter median bias.
 * 
 * @param output_data:    the output gravi_data where the cleaned imagelist
 *                        will be saved as IMAGING_DATA_ACQ.
 * @param input_data:     the input gravi_data with the raw imagelist
 * @param bad_map:        the gravi_data containing the bad pixel map
 *                        for ACQ (in extension IMAGING_DATA_ACQ).
 *
 * \exception CPL_ERROR_NULL_INPUT input data is missing
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_preproc_acqcam (gravi_data *output_data,
                                     gravi_data *input_data,
                                     gravi_data *bad_map)
{
    gravi_msg_function_start(1);
    cpl_ensure_code (output_data, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code (input_data,  CPL_ERROR_NULL_INPUT);

    /* Check if extension exist */
    if (!gravi_data_has_extension (input_data, GRAVI_IMAGING_DATA_ACQ_EXT)) {
        gravi_msg_warning (cpl_func,"Cannot preproc the ACQCAM, not data in file");
        return CPL_ERROR_NONE;
    }
    if (!gravi_data_has_extension (bad_map, GRAVI_IMAGING_DATA_ACQ_EXT)) {
        gravi_msg_warning (cpl_func,"Cannot preproc the ACQCAM, no badpixel in BAD");
        return CPL_ERROR_NONE;
    }

    /* Construct a mask of badpixels */
    cpl_image * badpix_img = gravi_data_get_img (bad_map, GRAVI_IMAGING_DATA_ACQ_EXT);
    cpl_mask * badpix_mask = cpl_mask_threshold_image_create (badpix_img, 0.5, 10000);
    CPLCHECK_MSG ("Cannot get BAD map for ACQ");

    /* Get the imagelist */
    cpl_imagelist * imglist;
    imglist = gravi_data_get_cube (input_data, GRAVI_IMAGING_DATA_ACQ_EXT);
    CPLCHECK_MSG ("Cannot get image for ACQ");
    
    /* Allocate new memory */
    imglist = cpl_imagelist_duplicate (imglist);
    
    /* 
     * Loop on images to cleanup-badpixels 
     */
    cpl_size nrow = cpl_imagelist_get_size (imglist);
    for (cpl_size row = 0; row < nrow; row++) {

        /* Get image */
        cpl_image * img = cpl_imagelist_get (imglist, row);

        /* Cleanup-badpixel */
        cpl_image_reject_from_mask (img, badpix_mask);
        cpl_detector_interpolate_rejected (img);
        CPLCHECK_MSG ("Cannot clean-up badpixel");
    }
    
    FREE (cpl_mask_delete, badpix_mask);

    /* Get the size */
    cpl_image * img = cpl_imagelist_get (imglist, 0);
    cpl_size nx = cpl_image_get_size_x (img);
    cpl_size ny = cpl_image_get_size_y (img);

    /* FIXME: Deal with the case full-frame in a better way */
    cpl_size ury = (ny>1100) ? 1200 : 750;        

    /* 
     * Remove the pupil background by the mean of blinking
     */

    if (nrow == 1) {
        gravi_msg_warning ("FIXME","Cannot remove blinked pupil (no blink)");
        
    } else {
        cpl_msg_info (cpl_func, "Remove the blinking");

        /* Get pupil diode flux for all images. Deal with the
         * change of format (FIXME: make it more clean) */
        double flux[nrow];
        for (cpl_size row = 0; row < nrow; row ++) {
            cpl_image * img = cpl_imagelist_get (imglist,row);
            if (nx < 1100)                
                flux[row] = cpl_image_get_flux_window (img, 1, 800, 1000, 1000) - 
                            cpl_image_get_median_window (img, 1, 800, 1000, 1000) * 1000 * 200;
            else
                flux[row] = cpl_image_get_flux_window (img, 1, 1200, 2048, 1536) -
                            cpl_image_get_median_window (img, 1, 1200, 2048, 1536) * 2048 * 336;

            cpl_msg_debug (cpl_func, "flux %lli = %.2f",row, flux[row]);
        }
            
        for (cpl_size row = 0; row < nrow; row ++) {

            /* Find the best possible blink, default is
             * current (hence frame will be zero) */
            cpl_size blink = row;
            if ( (flux[row] - flux[CPL_MAX(row-1,0)] ) > GRAVI_ACQ_PUP_FLUX) blink = row-1;
            if ( (flux[row] - flux[CPL_MIN(row+1,nrow-1)] ) > GRAVI_ACQ_PUP_FLUX) blink = row+1;

            cpl_msg_debug (cpl_func, "row %lli is debiased with %lli",row,blink);
            
            /* Remove the blink only for pupil part of the camera */
            gravi_image_subtract_window (cpl_imagelist_get (imglist, row),
                                         cpl_imagelist_get (imglist, blink),
                                         1, ury, nx, ny, 1, ury);
            CPLCHECK_MSG ("Cannot remove blinked pupil");
        }
    }

    /* 
     * Run a median filter for bias
     */
    cpl_msg_info (cpl_func, "Remove a median filter");
    
    cpl_mask * kernel = cpl_mask_new (11, 11);
    cpl_mask_not (kernel);
    
    for (cpl_size row = 0; row < nrow; row ++) {
        if (row %10 == 0)
            cpl_msg_info_overwritable (cpl_func, "Median filter of ACQ %lld over %lld", row+1, nrow);
        
        cpl_image * img = cpl_imagelist_get (imglist, row);
        cpl_image * unfiltered_img = cpl_image_extract (img, 1, ury, nx, ny);
        cpl_image * filtered_img   = cpl_image_duplicate (unfiltered_img);
        
        cpl_image_filter_mask (filtered_img, unfiltered_img, kernel,
                               CPL_FILTER_MEDIAN, CPL_BORDER_FILTER);
        gravi_image_subtract_window (img, filtered_img,
                                     1, ury, nx, ny, 1, 1);
        
        FREE (cpl_image_delete, unfiltered_img);
        FREE (cpl_image_delete, filtered_img);
        CPLCHECK_MSG ("Cannot run median filter");
    }

    FREE (cpl_mask_delete, kernel);

    /* 
     * Set in output 
     */
    gravi_data_add_cube (output_data, NULL, GRAVI_IMAGING_DATA_ACQ_EXT, imglist);
    
    gravi_msg_function_exit(1);
    return CPL_ERROR_NONE;   
}

/*----------------------------------------------------------------------------*/

/* Fast sin function */
double sin1 (double x)
{
    /* Within -pi +pi*/
    while (x >=  CPL_MATH_PI) x -= CPL_MATH_2PI;
    while (x <= -CPL_MATH_PI) x += CPL_MATH_2PI;

    /* Within -pi/2 +pi/2*/
    if (x >  CPL_MATH_PI_2) x =  CPL_MATH_PI - x;
    if (x < -CPL_MATH_PI_2) x = -CPL_MATH_PI - x;

    double x2 = x*x, x3 = x2*x, x5 = x3*x2;
    return 0.9996949 * x - 0.1656700 * x3 + 0.0075134 * x5;
}

/* Fast exp function */
double exp1 (double x) {
  x = 1.0 + x / 256.0;
  x *= x; x *= x; x *= x; x *= x;
  x *= x; x *= x; x *= x; x *= x;
  return x;
}

/* Compute the 4 diode positions from {angle, scaling, dx, dy} 
 * The model assume the 4 diodes form a rectangle centered on
 * the pupil. Hence this model has a 180deg degeneracy */
int gravi_acqcam_xy_diode (const double v[], double *xd, double *yd)
{
    double dx = v[GRAVI_SPOT_DIODE+0];
    double dy = v[GRAVI_SPOT_DIODE+1];
    
    /* Angle */
    double ang = (v[GRAVI_SPOT_ANGLE] - v[GRAVI_SPOT_DIODE+2])* CPL_MATH_RAD_DEG;
    double sang = sin1 (ang) * v[GRAVI_SPOT_SCALE];
    double cang = sin1 (ang + CPL_MATH_PI_2) * v[GRAVI_SPOT_SCALE];
    
    /* Diode arrangement */
    xd[0] = -cang * dx + sang * dy;
    xd[1] = -cang * dx - sang * dy;
    xd[2] =  cang * dx - sang * dy;
    xd[3] =  cang * dx + sang * dy;
    
    yd[0] =  cang * dy + sang * dx;
    yd[1] = -cang * dy + sang * dx;
    yd[2] = -cang * dy - sang * dx;
    yd[3] =  cang * dy - sang * dx;
    
    return 0;
}

/* Compute the 4 sub-aperture positions from the sub-aperture modes */
static int gravi_acqcam_xy_sub (const double v[], double *xsub, double *ysub)
{
    /* Sub-apperture arrangement */
    const double * vd = v + GRAVI_SPOT_SUB;
    xsub[0] = vd[0] + vd[1] + vd[2] + vd[3];
    xsub[1] = vd[0] - vd[1] + vd[2] - vd[3];
    xsub[2] = vd[0] + vd[1] - vd[2] - vd[3];
    xsub[3] = vd[0] - vd[1] - vd[2] + vd[3];
    
    ysub[0] = vd[4] + vd[5] + vd[6] + vd[7];
    ysub[1] = vd[4] - vd[5] + vd[6] - vd[7];
    ysub[2] = vd[4] + vd[5] - vd[6] - vd[7];
    ysub[3] = vd[4] - vd[5] - vd[6] + vd[7];

    return 0;
}

/*----------------------------------------------------------------------------*/

static int gravi_acqcam_spot (const double x_in[], const double v[], double *result)
{
    *result = 0.0;

    /* Static parameters */
    double xd[4], yd[4], xsub[4], ysub[4];
    gravi_acqcam_xy_diode (v, xd, yd);
    gravi_acqcam_xy_sub (v, xsub, ysub);

    double fwhm2 = v[GRAVI_SPOT_FWHM];
    
    /* Loop on diode and appertures.
     * The capture range is 2.FWHM */
    for (int diode = 0; diode < 4 ; diode++) {
        for (int sub = 0; sub < 4 ; sub++) {
            double xf = (x_in[0] - xsub[sub] - xd[diode]);
            double yf = (x_in[1] - ysub[sub] - yd[diode]);
            double dist = (xf*xf + yf*yf) / fwhm2;
            if (dist < 4.0) *result += v[GRAVI_SPOT_FLUX+sub*4+diode] * exp1 (-dist);
        }
    }

    return 0;
}

/*----------------------------------------------------------------------------*/

static int gravi_acqcam_spot_dfda (const double x_in[], const double v[], double result[])
{
    double next = 0.0, here = 0.0, epsilon = 1e-8;

    double vlocal[GRAVI_SPOT_NA];
    memcpy (vlocal, v, sizeof(double)*GRAVI_SPOT_NA);

    /* Compute value in-place */
    gravi_acqcam_spot (x_in, vlocal, &here);

    /* Fill with zeros */
    for (int a = 0; a < GRAVI_SPOT_NA; a++) result[a] = 0.0;
        
    /* Loop on parameters to compute finite differences 
     * FIXME: The analytical derivative may be faster, but
     * wasn't true in first tests, thus keep these. */
    for (int a = 0; a < GRAVI_SPOT_FLUX; a++) {
        if (GRAVI_LVMQ_FREE[a] != 0) {
            vlocal[a] += epsilon;
            gravi_acqcam_spot (x_in, vlocal, &next);
            vlocal[a] -= epsilon;
            result[a] = (next - here) / epsilon;
        }
    }

    /* The intensities are trivial analytic derivative */
    for (int a = GRAVI_SPOT_FLUX; a < GRAVI_SPOT_NA; a++) {
        if (GRAVI_LVMQ_FREE[a] != 0) {
            result[a] = here / v[a];
        }
    }

    return 0;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Imprint a cross (pixel=0) in the image of the pupil spot
 * 
 * @param img:    input/output img
 * @param a:      paramters for spot position
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_acqcam_spot_imprint (cpl_image * img, cpl_vector * a)
{
    gravi_msg_function_start(0);
    cpl_ensure_code (img, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code (a,   CPL_ERROR_NULL_INPUT);

    cpl_size nx = cpl_image_get_size_x (img);
    cpl_size ny = cpl_image_get_size_y (img);

    /* Static parameters */
    double xd[4], yd[4], xsub[4], ysub[4];
    double * v = cpl_vector_get_data (a);
    gravi_acqcam_xy_diode (v, xd, yd);
    gravi_acqcam_xy_sub (v, xsub, ysub);
    
    /* Loop on diode and appertures */
    for (int diode = 0; diode < 4 ; diode++) {
        for (int sub = 0; sub < 4 ; sub++) {
            cpl_size xf = roundl(xsub[sub] + xd[diode]) + 1;
            cpl_size yf = roundl(ysub[sub] + yd[diode]) + 1;
	    if (xf < 2 || xf > nx-2 || yf < 2 || yf > ny-2) continue;
            cpl_image_set (img, xf,   yf, 0);
            cpl_image_set (img, xf-1, yf, 0);
            cpl_image_set (img, xf+1, yf, 0);
            cpl_image_set (img, xf, yf+1, 0);
            cpl_image_set (img, xf, yf-1, 0);
            CPLCHECK ("Cannot imprint cross in image");
        }
    }

    gravi_msg_function_exit(0);
    return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Read the diode position from header into the vector output
 * 
 * @param header:   input header
 * @param tel:      requested beam (0..3)
 * @param output:   output vector, shall be already allocated
 *
 * \exception CPL_ERROR_NULL_INPUT input data is missing
 * \exception CPL_ERROR_ILLEGAL_INPUT tel outside limits
 *
 * Read the diode position from header into the vector output. Assume
 *        the four diodes form a rectangle centered on the pupil center.
 *
 *        - output[8]   = rotation  [deg], set to 0.0
 *        - output[9]   = scale  [pix/m]
 *        - output[10]  = dx [m]
 *        - output[11]  = dy [m]
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_acqcam_get_diode_ref (cpl_propertylist * header,
                                           cpl_size tel,
                                           cpl_vector * output)
{
    gravi_msg_function_start(0);
    cpl_ensure_code (header,          CPL_ERROR_NULL_INPUT);
    cpl_ensure_code (tel>=0 && tel<4, CPL_ERROR_ILLEGAL_INPUT);
    cpl_ensure_code (output,          CPL_ERROR_ILLEGAL_INPUT);

    /* Get the telescope name and ID */
    const char * telname = gravi_conf_get_telname (tel, header);

    /* Check telescope name */
    if (!telname) cpl_msg_error (cpl_func, "Cannot read the telescope name");
    cpl_ensure_code (telname, CPL_ERROR_ILLEGAL_INPUT);
    
    /* Hardcoded theoretical positions in mm */

    /* If UTs or ATs, select scaling, rotation, and spacing */
    if (telname[0] == 'U') {
        cpl_vector_set (output, GRAVI_SPOT_ANGLE,  0.0);
        cpl_vector_set (output, GRAVI_SPOT_SCALE, 16.225);
        cpl_vector_set (output, GRAVI_SPOT_DIODE+0, 0.363);
        cpl_vector_set (output, GRAVI_SPOT_DIODE+1, 0.823);
        cpl_vector_set (output, GRAVI_SPOT_DIODE+2, 0);
    } else if (telname[0] == 'A') {
        cpl_vector_set (output, GRAVI_SPOT_ANGLE,  0.0);
        cpl_vector_set (output, GRAVI_SPOT_SCALE, 73.0154);
        cpl_vector_set (output, GRAVI_SPOT_DIODE+0, 0.122);
        cpl_vector_set (output, GRAVI_SPOT_DIODE+1, 0.158);
        cpl_vector_set (output, GRAVI_SPOT_DIODE+2, 0);
    } else 
        return cpl_error_set_message (cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                      "Cannot get telescope name");
        
    gravi_msg_function_exit(0);
    return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Re-arrange the sub-aperture position into the output vector
 *
 * @param header:   input header
 * @param tel:      requested beam (0..3)
 * @param output:   output vector, shall be already allocated
 *
 * \exception CPL_ERROR_NULL_INPUT input data is missing
 * \exception CPL_ERROR_ILLEGAL_INPUT tel outside limits
 *
 *
 * Read the sub-aperture position for the pupil sensor,
 *        and re-arrange them into the output vector
 *
 *        - output[0] = x0+x1+x2+x3   (center of sub-apertures)
 *        - output[1] = x0-x1+x2-x3
 *        - output[2] = x0+x1-x2-x3
 *        - output[3] = x0-x1-x2+x3
 *        - output[4..7] = same for y
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_acqcam_get_pup_ref (cpl_propertylist * header,
                                         cpl_size tel,
                                         cpl_vector * output)
{
    gravi_msg_function_start(0);
    cpl_ensure_code (header,          CPL_ERROR_NULL_INPUT);
    cpl_ensure_code (tel>=0 && tel<4, CPL_ERROR_ILLEGAL_INPUT);
    cpl_ensure_code (output,          CPL_ERROR_ILLEGAL_INPUT);

    cpl_size ntel = 4, nsub = 4;
    cpl_size nsx = 0, nsy = 0, sx = 0, sy = 0;

    /* If sub-windowing, we read the sub-window size and
     * the sub-window start for pupil */
    if ( cpl_propertylist_has (header, "ESO DET1 FRAMES NX") ) {
        char name[90];
        
        nsx = cpl_propertylist_get_int (header, "ESO DET1 FRAMES NX");
        sprintf (name, "ESO DET1 FRAM%lld STRX", 3*ntel + tel + 1);
        sx = cpl_propertylist_get_int (header, name);
        
        nsy = cpl_propertylist_get_int (header, "ESO DET1 FRAMES NY");
        sprintf (name, "ESO DET1 FRAM%lld STRY", 3*ntel + tel + 1);
        sy = cpl_propertylist_get_int (header, name);
        
        CPLCHECK_MSG ("Cannot get sub-windowing parameters");
    }
    
    cpl_msg_debug (cpl_func,"sub-window pupil %lli sx= %lld sy = %lld", tel, sx, sy);
        
    /* Read the sub-apperture reference positions 
     * Converted to accound for sub-windowing 
     * In vector convention (start at 0,0) */
    double xsub[4], ysub[4];
    for (int sub = 0; sub < nsub ; sub++) {
        double xv = gravi_pfits_get_ptfc_acqcam (header, sub*ntel + tel + 1);
        double yv = gravi_pfits_get_ptfc_acqcam (header, sub*ntel + tel + 17);
        xsub[sub] = xv - (sx - tel*nsx) - 1;
        ysub[sub] = yv - (sy - 3*nsy) - 1;
        cpl_msg_debug (cpl_func,"pupil %lli subC %i = %10.4f,%10.4f",
                       tel, sub, xsub[sub], ysub[sub]);
        CPLCHECK_MSG ("Cannot get pupil reference position");
    }

    /* All linear combination of sub-appertures center */
    cpl_vector_set (output, GRAVI_SPOT_SUB+0, 0.25 * (xsub[0] + xsub[1] + xsub[2] + xsub[3]));
    cpl_vector_set (output, GRAVI_SPOT_SUB+1, 0.25 * (xsub[0] - xsub[1] + xsub[2] - xsub[3]));
    cpl_vector_set (output, GRAVI_SPOT_SUB+2, 0.25 * (xsub[0] + xsub[1] - xsub[2] - xsub[3]));
    cpl_vector_set (output, GRAVI_SPOT_SUB+3, 0.25 * (xsub[0] - xsub[1] - xsub[2] + xsub[3]));
    
    cpl_vector_set (output, GRAVI_SPOT_SUB+4, 0.25 * (ysub[0] + ysub[1] + ysub[2] + ysub[3]));
    cpl_vector_set (output, GRAVI_SPOT_SUB+5, 0.25 * (ysub[0] - ysub[1] + ysub[2] - ysub[3]));
    cpl_vector_set (output, GRAVI_SPOT_SUB+6, 0.25 * (ysub[0] + ysub[1] - ysub[2] - ysub[3]));
    cpl_vector_set (output, GRAVI_SPOT_SUB+7, 0.25 * (ysub[0] - ysub[1] - ysub[2] + ysub[3]));
    
    gravi_msg_function_exit(0);
    return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Fit a pupil spot pattern into an image.
 * 
 * @param img:    input image
 * @param ntry:   number of random starting point
 * @param a:      vector of parameter, modified in-place
 * @param fitAll: if set to 1, the high order of sub-aperture position (more
 *                than center and focus), the diode scaling, and the diode
 *                intensities are let free.
 * @param nspot:  is filled with the number of detected spots.
 * 
 * Fit a pupil spot pattern into an image. The global minimum is found
 * first with a fit with a large capture range in the whole image (after
 * a 5x5 binning to minimize the computation time).
 * The number of random starting point is given by
 * the ntry parameter. If ntry==1, the starting is kept unmodified.
 * Then 10x10 pixels around each expected spot are extracted and fit
 * with true Gaussian (free FWHM and amplitude).
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_acqcam_fit_spot (cpl_image * img,
                                      cpl_size ntry,
                                      cpl_vector * a,
                                      int fitAll,
                                      int * nspot)
{
    gravi_msg_function_start(0);
    cpl_ensure_code (img, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code (a,   CPL_ERROR_NULL_INPUT);
    cpl_ensure_code (ntry>0, CPL_ERROR_ILLEGAL_INPUT);
    cpl_ensure_code (nspot, CPL_ERROR_NULL_INPUT);

    int nv = 0;
    cpl_size x0 = cpl_vector_get (a, GRAVI_SPOT_SUB+0);
    cpl_size y0 = cpl_vector_get (a, GRAVI_SPOT_SUB+4);
    CPLCHECK_MSG ("Cannot get values valid");

    /* Compute RMS in the central region */
    double RMS = gravi_image_get_noise_window (img, x0-25, y0-25, x0+25, y0+25);

    /* The image is surely empty */
    if (RMS == 0) { *nspot = 0; return CPL_ERROR_NONE;}

    /*
     * Coarse: fit with a re-bin image
     */

    /* To lower the number of point, we extract a window of 175x175
     * around the center of sub-apertures, rebin with 5x5 pixels */
    cpl_size nw = 175, n_mean = 5;
    cpl_size nc = nw/n_mean, nint = n_mean*n_mean;
    
    /* Allocate vector and matrix for the fit */
    cpl_matrix * x_matrix  = cpl_matrix_new (nc*nc, 2);
    cpl_vector * y_vector  = cpl_vector_new (nc*nc);
    cpl_vector * sy_vector = cpl_vector_new (nc*nc);

    /* Loop on pixel in the considered window */
    for (cpl_size x = 0; x < nc; x++) {
      for (cpl_size y = 0; y < nc; y++) {
	/* Average 3x3 pixels */
	double z_mean = 0.0, x_mean = 0.0, y_mean = 0.0;
	for (int i=-2; i<=2; i++) {
	  for (int j=-2; j<=2; j++) {
	    cpl_size x1 = x*n_mean+i+x0-(n_mean*nc)/2;
	    cpl_size y1 = y*n_mean+j+y0-(n_mean*nc)/2;
	    z_mean += cpl_image_get (img, x1+1, y1+1, &nv);
	    x_mean += x1; y_mean += y1;
	  }
	}
	/* Set in matrix and vectors */
	cpl_matrix_set (x_matrix, x*nc+y, 0, x_mean/nint);
	cpl_matrix_set (x_matrix, x*nc+y, 1, y_mean/nint);
	cpl_vector_set (y_vector, x*nc+y, z_mean/nint);
	cpl_vector_set (sy_vector, x*nc+y, RMS);
	CPLCHECK_MSG ("Cannot fill matrix/vector");
      }
    } /* End loop on re-sampled pixels*/

    /* Output for global minimisation */
    cpl_vector * a_start = cpl_vector_duplicate (a);
    cpl_vector * a_tmp = cpl_vector_duplicate (a);
    double chisq_final = 1e10;
    srand(1);

    /* Set the fwhm to 6 to force a large capture range */
    cpl_vector_set (a_start, GRAVI_SPOT_FWHM, 6.*6.);

    /* If needed, define a realistic value for the amplitude (for numerical stability) */
    if (fitAll) for (int d=0;d<16;d++) cpl_vector_set (a_start, GRAVI_SPOT_FLUX+d, RMS);

    /* Fit sub-aperture mean position; and diode rotation */
    const int ia_global[] = {1,0,0,0, 1,0,0,0, 1,0, 0,0,0, 0,
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    GRAVI_LVMQ_FREE = ia_global;

    /* Loop on various starting points */
    for (cpl_size try = 0; try < ntry; try++) {

        /* Move starting point in position (+-10pix) and angle (entire circle) */
        cpl_vector_copy (a_tmp, a_start);
        if (try > 0) {
            cpl_vector_set (a_tmp, GRAVI_SPOT_SUB+0,
                            cpl_vector_get (a_tmp, GRAVI_SPOT_SUB+0) +
                            (rand()%20) - 10);
            cpl_vector_set (a_tmp, GRAVI_SPOT_SUB+4,
                            cpl_vector_get (a_tmp, GRAVI_SPOT_SUB+4) +
                            (rand()%20) - 10);
            cpl_vector_set (a_tmp, GRAVI_SPOT_ANGLE,
                            cpl_vector_get (a_tmp, GRAVI_SPOT_ANGLE) +
                            (rand()%180));
        }

        /* Fit from this starting point */
        double chisq;
        cpl_fit_lvmq (x_matrix, NULL, y_vector, sy_vector,
                      a_tmp, GRAVI_LVMQ_FREE,
                      gravi_acqcam_spot,
                      gravi_acqcam_spot_dfda,
                      CPL_FIT_LVMQ_TOLERANCE,
                      CPL_FIT_LVMQ_COUNT,
                      CPL_FIT_LVMQ_MAXITER,
                      NULL, &chisq, NULL);
        CPLCHECK_MSG ("Cannot fit global minimum");

        /* Check chi2 and keep if lowest so far */
        if (chisq < chisq_final) {
            cpl_msg_debug (cpl_func, "chisq_final = %f -> %f", chisq_final, chisq);
            chisq_final = chisq;
            cpl_vector_copy (a, a_tmp);
        }

    } /* End loop on try starting points */

    FREE (cpl_matrix_delete, x_matrix);
    FREE (cpl_vector_delete, y_vector);
    FREE (cpl_vector_delete, sy_vector);
    FREE (cpl_vector_delete, a_tmp);
    FREE (cpl_vector_delete, a_start);

    
    /* 
     * Fine: fit 10 pixel around each spot with true Gaussian
     */
    
    cpl_size nf = 10, ndiode = 4, nsub = 4;
    x_matrix  = cpl_matrix_new (nf*nf*ndiode*nsub, 2);
    y_vector  = cpl_vector_new (nf*nf*ndiode*nsub);
    sy_vector = cpl_vector_new (nf*nf*ndiode*nsub);
    
    double xd[4], yd[4], xsub[4], ysub[4];
    gravi_acqcam_xy_diode (cpl_vector_get_data (a), xd, yd);
    gravi_acqcam_xy_sub (cpl_vector_get_data (a), xsub, ysub);
    
    /* Loop on diode and appertures to fill the matrix and
     * vector for the fine fit */
    for (cpl_size v = 0, diode = 0; diode < ndiode ; diode++) {
        for (int sub = 0; sub < nsub ; sub++) {
            cpl_size xf = roundl(xsub[sub] + xd[diode]);
            cpl_size yf = roundl(ysub[sub] + yd[diode]);
            
            /* Extract 10 pixels around each spot */
            double mf = 0.0;
            for (cpl_size x = xf-nf/2; x < xf+nf/2; x++) {
                for (cpl_size y = yf-nf/2; y < yf+nf/2; y++) {
                    double value = cpl_image_get (img, x+1, y+1, &nv);
                    if (value > mf) mf = value;
                    cpl_matrix_set (x_matrix, v, 0, x);
                    cpl_matrix_set (x_matrix, v, 1, y);
                    cpl_vector_set (y_vector, v, value);
                    cpl_vector_set (sy_vector, v, RMS);
                    v++;
                    CPLCHECK_MSG ("Cannot fill matrix");
                }
            }
            
            /* Save the local maximum of this spot as
             * the starting point for its amplitude */
            cpl_vector_set (a, GRAVI_SPOT_FLUX+sub*4+diode, mf);
        }
    }

    /* Set FWHM to a realist value in [pix**2] */
    cpl_vector_set (a, GRAVI_SPOT_FWHM, 2.3*2.3);

    /* Fit all sub-aperture position; rotation and scaling of diodes;
     * and individual intensities of spots. (high order and scaling
     * are only fitted if fitAll != 0) */
    int F = fitAll;
    const int ia_fine[] = {1,F,1,F, 1,1,F,F, 1,F, 0,0,0, 0,
                           1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
    GRAVI_LVMQ_FREE = ia_fine;

    /* Fit from this starting point */
    double chisq_fine = 0.0;
    cpl_errorstate prestate = cpl_errorstate_get();
    cpl_fit_lvmq (x_matrix, NULL, y_vector, sy_vector,
                  a, GRAVI_LVMQ_FREE,
                  gravi_acqcam_spot,
                  gravi_acqcam_spot_dfda,
                  CPL_FIT_LVMQ_TOLERANCE,
                  CPL_FIT_LVMQ_COUNT,
                  CPL_FIT_LVMQ_MAXITER,
                  NULL, &chisq_fine, NULL);

    /* Recover on error but return NULL */
    if (cpl_error_get_code() == CPL_ERROR_CONTINUE) {
      cpl_errorstate_set (prestate);
      cpl_msg_warning (cpl_func, "Cannot fit pupil spots...");
      for (int d=0; d<16; d++) cpl_vector_set (a, GRAVI_SPOT_FLUX+d, 0);
    }
    
    CPLCHECK_MSG ("Cannot fit fine");

    cpl_msg_debug (cpl_func, "chisq_final = %.2f -> fine = %.2f",
                   chisq_final, chisq_fine);
    
    FREE (cpl_matrix_delete, x_matrix);
    FREE (cpl_vector_delete, y_vector);
    FREE (cpl_vector_delete, sy_vector);

    /* Count the number of valid spots based on their amplitude */
    *nspot = 0;
    for (int d=0; d<16; d++) if (cpl_vector_get (a, GRAVI_SPOT_FLUX+d) > 3.*RMS) (*nspot)++;

    gravi_msg_function_exit(0);
    return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/*
 * @brief measure Strehl Ratio of the source at the given location
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_acq_measure_strehl(cpl_image * img, double x, double y, 
                                        double pscale, double *SR, cpl_propertylist * header)
{
    gravi_msg_function_start(0);
    cpl_ensure_code (img, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code (SR, CPL_ERROR_NULL_INPUT);

    const char * telname = gravi_conf_get_telname (0, header);
    CPLCHECK ("Cannot get telescope name");
    hdrl_parameter * strehl_params;

    /* Hardcoded Mirror diameters */
    if (telname[0] == 'U') {
        strehl_params = 
            hdrl_strehl_parameter_create (1.8e-6, 8.0/2.0, 1.126/2.0,
                                          pscale*1e-3, pscale*1e-3,
                                          0.8, 0.8, 1.0);
    } else if (telname[0] == 'A') {
        strehl_params = 
            hdrl_strehl_parameter_create (1.8e-6, 1.8/2.0, 0.138/2.0,
                                          pscale*1e-3, pscale*1e-3,
                                          0.8, 0.8, 1.0);
    } else {
        cpl_error_set_message (cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                               "Cannot recognise the telescope");
        return CPL_ERROR_ILLEGAL_INPUT;
    }

    hdrl_strehl_result strehl;
    cpl_image * sub_image = cpl_image_extract (img, x-50, y-50, x+50, y+50);
    hdrl_image * sub_hdrl = hdrl_image_create (sub_image, NULL);

    strehl = hdrl_strehl_compute (sub_hdrl, strehl_params);
    *SR = (double) strehl.strehl_value.data;
    FREE (hdrl_image_delete, sub_hdrl);
    FREE (cpl_image_delete, sub_image);
    hdrl_parameter_delete(strehl_params);

    gravi_msg_function_exit(0);
    return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/*
 * @brief measure maximum of the source at the given location
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_acq_measure_max(cpl_image * img, double x, double y, double size, double * img_max)
{
    gravi_msg_function_start(0);
    cpl_ensure_code (img, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code (img_max, CPL_ERROR_NULL_INPUT);
    
    cpl_image * small_img = cpl_image_extract(img, x-size, y-size, x+size, y+size);
    *img_max = cpl_image_get_max(small_img);
    cpl_image_delete(small_img);
    
    gravi_msg_function_exit(0);
    return CPL_ERROR_NONE;
}


/*----------------------------------------------------------------------------*/
/**
 * @brief Fit a Gaussian into an image, and mark the position.
 * 
 * @param img:    input image
 * @param x,y:    input/output position (guess and best fit)
 * @param ex,ey:  output uncertainties on x and y
 * @param size:   size of box to consider
 * 
 * \exception CPL_ERROR_NULL_INPUT input data is missing
 *
 * The function use cpl_fit_image_gaussian to fit a Gaussian into the image
 * The best-fit position is then fill with 0 in the image.
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_acq_fit_gaussian (cpl_image * img, double *x, double *y,
                                       double *ex, double *ey, cpl_size size)
{
    gravi_msg_function_start(0);
    cpl_ensure_code (img, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code (x,  CPL_ERROR_NULL_INPUT);
    cpl_ensure_code (y,  CPL_ERROR_NULL_INPUT);
    cpl_ensure_code (ex, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code (ey, CPL_ERROR_NULL_INPUT);

    /* Fill first guess */
    cpl_array * parameters = cpl_array_new (7, CPL_TYPE_DOUBLE);
    cpl_array_fill_window_invalid (parameters, 0, 7);
    cpl_array_set (parameters, 2, 0);
    cpl_array_set (parameters, 3, *x);
    cpl_array_set (parameters, 4, *y);
    cpl_array_set (parameters, 5, 3);
    cpl_array_set (parameters, 6, 3);
    CPLCHECK_MSG ("Error creating parameter table");
    
    double med = cpl_image_get_median_window (img,
                     (cpl_size)(*x)-size, (cpl_size)(*y)-size,
                     (cpl_size)(*x)+size, (cpl_size)(*y)+size);
    cpl_array_set (parameters, 0, med);
    CPLCHECK_MSG ("Error getting median");
    
    /* Fit Gaussian */
    double rms=0.;
    cpl_fit_image_gaussian (img, NULL, (cpl_size)(*x), (cpl_size)(*y),
                            size, size, parameters,
                            NULL, NULL, &rms, NULL, NULL,
                            NULL, NULL, NULL, NULL);
    CPLCHECK_MSG ("Error fitting Gaussian");

    /* Set back */
    *x = cpl_array_get (parameters,3,NULL);
    *y = cpl_array_get (parameters,4,NULL);
    CPLCHECK_MSG ("Error reading fit result");

    /* Check errors */
    /* reject result if either:       */
    /*          * peak is below 3*rms */
    /*          * sigma_x > 5 pixels  */
    /*          * sigma_y > 5 pixels  */
    double A   = cpl_array_get (parameters, 1, NULL);
    double rho = cpl_array_get (parameters, 2, NULL);
    double sx  = cpl_array_get (parameters, 5, NULL);
    double sy  = cpl_array_get (parameters, 6, NULL);

    if ( A < 0. ) {
      // detection is just not significant
      cpl_msg_info (cpl_func, "rejecting fit: x=%g, y=%g, SNR=%g, sx=%g, sy=%g",
                    *x, *y, A/(rms * 2.*M_PI*sx*sy*sqrt(1-rho*rho)), sx, sy);
      *x = 0.;
      *y = 0.;
      *ex = -1.;
      *ey = -1.;
    } else {
      // cf. Condon 1996, PASP 109:166
      double cst = 2. * sqrt(2. * M_PI * (1. - rho*rho) * sx * sy) * rms / A;
      *ex=cst*sx;
      *ey=cst*sy;
    }
    CPLCHECK_MSG("Error checking significance of fit result");

    /* Fill image with zero at the detected position */
    //if (*x > 0. && *y > 0.) {
    //  cpl_image_set (img, (cpl_size)(*x), (cpl_size)(*y), 0.0);
    //}
    //CPLCHECK_MSG("Error setting peak to zero");
    
    /* Delete */
    FREE (cpl_array_delete, parameters);
    
    gravi_msg_function_exit(0);
    return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Reduce the images of pupil from ACQ
 * 
 * @param mean_img:          input mean image
 * @param acqcam_imglist:    input image list
 * @param header:            input header
 * @param acqcam_table:      output table
 * @param o_header:          output header
 * 
 * \exception CPL_ERROR_NULL_INPUT input data is missing
 *
 * The routine analyse the pupil from ACQ and create QC parameters in the
 * header, as well as columns in the acqcam_table.
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_acqcam_pupil (cpl_image * mean_img,
                                   cpl_imagelist * acqcam_imglist,
                                   cpl_propertylist * header,
                                   cpl_table * acqcam_table,
                                   cpl_propertylist * o_header)
{
    gravi_msg_function_start(1);
    cpl_ensure_code (mean_img,       CPL_ERROR_NULL_INPUT);
    cpl_ensure_code (acqcam_imglist, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code (header,         CPL_ERROR_NULL_INPUT);
    cpl_ensure_code (acqcam_table,   CPL_ERROR_NULL_INPUT);
    cpl_ensure_code (o_header,       CPL_ERROR_NULL_INPUT);

    /* Number of row */
    cpl_size nrow = cpl_imagelist_get_size (acqcam_imglist);

    /* Compute separation */
    double sobj_x = gravi_pfits_get_sobj_x (header);
    double sobj_y = gravi_pfits_get_sobj_y (header);
    double sobj_sep = sqrt(sobj_x*sobj_x + sobj_y*sobj_y);

    /* Pupil positions (or we use array of 3)  */
    gravi_table_new_column (acqcam_table, "PUPIL_NSPOT", NULL, CPL_TYPE_INT);
    gravi_table_new_column (acqcam_table, "PUPIL_X", "pix", CPL_TYPE_DOUBLE);
    gravi_table_new_column (acqcam_table, "PUPIL_Y", "pix", CPL_TYPE_DOUBLE);
    gravi_table_new_column (acqcam_table, "PUPIL_Z", "pix", CPL_TYPE_DOUBLE);
    gravi_table_new_column (acqcam_table, "PUPIL_R", "deg", CPL_TYPE_DOUBLE);
    gravi_table_new_column (acqcam_table, "PUPIL_U", "m", CPL_TYPE_DOUBLE);
    gravi_table_new_column (acqcam_table, "PUPIL_V", "m", CPL_TYPE_DOUBLE);
    gravi_table_new_column (acqcam_table, "PUPIL_W", "m", CPL_TYPE_DOUBLE);
    gravi_table_new_column (acqcam_table, "OPD_PUPIL", "m", CPL_TYPE_DOUBLE);
    
    int nspot = 0, ntel = 4;
    
    /* Loop on tel */
    for (int tel = 0; tel < ntel; tel++) {
        cpl_msg_info (cpl_func, "Compute pupil position for beam %i", tel+1);

        /* Get the conversion angle xy to uv in [rad] */
        double fangle = gravi_pfits_get_fangle_acqcam (header, tel);
        double cfangle = cos(fangle * CPL_MATH_RAD_DEG);
        double sfangle = sin(fangle * CPL_MATH_RAD_DEG);
        CPLCHECK_MSG ("Cannot read ESO INS DROTOFF#");
        
        /* Get the orientation of star */
        double drotoff = gravi_pfits_get_drotoff (header, tel);
        double cdrotoff = cos(drotoff * CPL_MATH_RAD_DEG);
        double sdrotoff = sin(drotoff * CPL_MATH_RAD_DEG);
        CPLCHECK_MSG ("Cannot read ESO INS DROTOFF#");

        /* Allocate memory */
        cpl_vector * a_start = cpl_vector_new (GRAVI_SPOT_NA);
        cpl_vector_fill (a_start, 0.0);
        
        /* Read the sub-apperture reference positions 
         * Converted to accound for sub-windowing 
         * In vector convention (start at 0,0) */
        gravi_acqcam_get_pup_ref (header, tel, a_start);
        gravi_acqcam_get_diode_ref (header, tel, a_start);
        CPLCHECK_MSG ("Cannot read ACQ data for pupil in header");

        cpl_vector * a_final = cpl_vector_duplicate (a_start);
        CPLCHECK_MSG ("Cannot prepare parameters");
                
        /* Fit pupil sensor spots in this image, with various
         * starting points to converge to the global minimum */
        gravi_acqcam_fit_spot (mean_img, 30, a_final, 1, &nspot);
        CPLCHECK_MSG ("Cannot fit rotation and center");

        /* Offsets in pixels */
        double xpos = cpl_vector_get (a_final,GRAVI_SPOT_SUB+0) -
                      cpl_vector_get (a_start,GRAVI_SPOT_SUB+0);
        double ypos = cpl_vector_get (a_final,GRAVI_SPOT_SUB+4) -
                      cpl_vector_get (a_start,GRAVI_SPOT_SUB+4);

        /* Scale and rotation in deg (rectangle, thus 180deg symetrie) */
        double scale = cpl_vector_get (a_final,GRAVI_SPOT_SCALE);
        double angle = cpl_vector_get (a_final,GRAVI_SPOT_ANGLE);
        if (angle < 0)   angle += 180;
        if (angle > 180) angle -= 180;
        cpl_vector_set (a_final,GRAVI_SPOT_ANGLE,angle);

        /* In UV [m] */
        double upos = (cfangle * xpos - sfangle * ypos) / scale;
        double vpos = (sfangle * xpos + cfangle * ypos) / scale;
        
        /* Add best position as a cross in image */
        gravi_acqcam_spot_imprint (mean_img, a_final);
        
        /* Add QC parameters */
        char qc_name[100];
        
        sprintf (qc_name, "ESO QC ACQ FIELD%i NORTH_ANGLE", tel+1);
        cpl_msg_info (cpl_func, "%s = %f", qc_name, fangle);
        cpl_propertylist_update_double (o_header, qc_name, fangle);
        cpl_propertylist_set_comment (o_header, qc_name, "[deg] y->x, predicted North direction on ACQ");
        
        sprintf (qc_name, "ESO QC ACQ PUP%i NSPOT", tel+1);
        cpl_msg_info (cpl_func, "%s = %i", qc_name, nspot);
        cpl_propertylist_update_int (o_header, qc_name, nspot);
        cpl_propertylist_set_comment (o_header, qc_name, "nb. of pupil spot in ACQ");

        sprintf (qc_name, "ESO QC ACQ PUP%i ANGLE", tel+1);
        cpl_msg_info (cpl_func, "%s = %f", qc_name, angle);
        cpl_propertylist_update_double (o_header, qc_name, angle);
        cpl_propertylist_set_comment (o_header, qc_name, "[deg] y->x, diode angle on ACQ");

        sprintf (qc_name, "ESO QC ACQ PUP%i SCALE", tel+1);
        cpl_msg_info (cpl_func, "%s = %f", qc_name, scale);
        cpl_propertylist_update_double (o_header, qc_name, scale);
        cpl_propertylist_set_comment (o_header, qc_name, "[pix/m] diode scale on ACQ");

        sprintf (qc_name, "ESO QC ACQ PUP%i XPOS", tel+1);
        cpl_msg_info (cpl_func, "%s = %f", qc_name, xpos);
        cpl_propertylist_update_double (o_header, qc_name, xpos);
        cpl_propertylist_set_comment (o_header, qc_name, "[pix] pupil x-shift in ACQ");
        
        sprintf (qc_name, "ESO QC ACQ PUP%i YPOS", tel+1);
        cpl_msg_info (cpl_func, "%s = %f", qc_name, ypos);
        cpl_propertylist_update_double (o_header, qc_name, ypos);
        cpl_propertylist_set_comment (o_header, qc_name, "[pix] pupil y-shift in ACQ");

        sprintf (qc_name, "ESO QC ACQ PUP%i UPOS", tel+1);
        cpl_msg_info (cpl_func, "%s = %f", qc_name, upos);
        cpl_propertylist_update_double (o_header, qc_name, upos);
        cpl_propertylist_set_comment (o_header, qc_name, "[m] pupil u-shift in ACQ");
        
        sprintf (qc_name, "ESO QC ACQ PUP%i VPOS", tel+1);
        cpl_msg_info (cpl_func, "%s = %f", qc_name, vpos);
        cpl_propertylist_update_double (o_header, qc_name, vpos);
        cpl_propertylist_set_comment (o_header, qc_name, "[m] pupil v-shift in ACQ");

        /* Loop on all images */
        for (cpl_size row = 0; row < nrow; row++) {
            if (row %10 == 0 || row == (nrow-1))
                cpl_msg_info_overwritable (cpl_func, "Fit image %lld over %lld", row+1, nrow);
            
            /* Get data and first guess */
            cpl_image * img = cpl_imagelist_get (acqcam_imglist, row);
            cpl_vector * a_row = cpl_vector_duplicate (a_final);

            /* Fit pupil sensor spots in this image. */
            gravi_acqcam_fit_spot (img, 1, a_row, 0, &nspot);
            CPLCHECK_MSG ("Cannot fit sub-appertures of image");

            /* If spot detected */
            if (nspot < 4) {
                cpl_table_set (acqcam_table, "PUPIL_NSPOT", row*ntel+tel, 0);
            }
            else {
                /* Add best position as a cross in image */
                gravi_acqcam_spot_imprint (img, a_row);

                /* Remove reference */
                cpl_vector_subtract (a_row, a_start);
                
                /* Get rotation [deg] and shift [pix]*/
                double r_shift = cpl_vector_get (a_row, GRAVI_SPOT_ANGLE);
                double x_shift = cpl_vector_get (a_row, GRAVI_SPOT_SUB+0);
                double y_shift = cpl_vector_get (a_row, GRAVI_SPOT_SUB+4);
                double z_shift = -0.5 * ( cpl_vector_get (a_row, GRAVI_SPOT_SUB+2) +
                                          cpl_vector_get (a_row, GRAVI_SPOT_SUB+5));
                
                /* In UV [m] */
                double u_shift = (cfangle * x_shift - sfangle * y_shift) / scale;
                double v_shift = (sfangle * x_shift + cfangle * y_shift) / scale;               
                double w_shift = gravi_acqcam_z2meter (z_shift);
                
                cpl_table_set (acqcam_table, "PUPIL_NSPOT", row*ntel+tel, nspot);
                cpl_table_set (acqcam_table, "PUPIL_X", row*ntel+tel, x_shift);
                cpl_table_set (acqcam_table, "PUPIL_Y", row*ntel+tel, y_shift);
                cpl_table_set (acqcam_table, "PUPIL_Z", row*ntel+tel, z_shift);
                cpl_table_set (acqcam_table, "PUPIL_R", row*ntel+tel, r_shift);
                cpl_table_set (acqcam_table, "PUPIL_U", row*ntel+tel, u_shift);
                cpl_table_set (acqcam_table, "PUPIL_V", row*ntel+tel, v_shift);
                cpl_table_set (acqcam_table, "PUPIL_W", row*ntel+tel, w_shift);

                /* Compute the OPD_PUPIL */
                double opd_pupil = sobj_sep * GRAVI_MATH_RAD_MAS / scale * 
                                   (x_shift * cdrotoff + y_shift * sdrotoff);
                cpl_table_set (acqcam_table, "OPD_PUPIL", row*ntel+tel, opd_pupil);
            }
            
            FREE (cpl_vector_delete, a_row);
        }

        FREE (cpl_vector_delete, a_start);
        FREE (cpl_vector_delete, a_final);
    } /* End loop on tel */
        
    gravi_msg_function_exit(1);
    return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Reduce the images of field from ACQ
 * 
 * @param mean_img:          input mean image
 * @param acqcam_imglist:    input image list
 * @param header:            input header
 * @param acqcam_table:      output table
 * @param o_header:          output header
 * 
 * \exception CPL_ERROR_NULL_INPUT input data is missing
 * \exception CPL_ERROR_ILLEGAL_INPUT Cannot find the data of the expected
 * telescope
 *
 * The routine analyse the field from ACQ and create QC parameters in the
 * header, as well as columns in the acqcam_table.
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_acqcam_field (cpl_image * mean_img,
                                   cpl_imagelist * acqcam_imglist,
                                   cpl_propertylist * header,
                                   cpl_table * acqcam_table,
                                   cpl_propertylist * o_header)
{
    gravi_msg_function_start(1);
    cpl_ensure_code (mean_img,       CPL_ERROR_NULL_INPUT);
    cpl_ensure_code (acqcam_imglist, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code (header,         CPL_ERROR_NULL_INPUT);
    cpl_ensure_code (acqcam_table,   CPL_ERROR_NULL_INPUT);
    cpl_ensure_code (o_header,       CPL_ERROR_NULL_INPUT);
    
    char qc_name[100];
    int ntel = 4;
    
    /* Number of row */
    cpl_size nrow = cpl_imagelist_get_size (acqcam_imglist);
    
    /* Create columns */
    gravi_table_new_column (acqcam_table, "FIELD_SC_X", "pix", CPL_TYPE_DOUBLE);
    gravi_table_new_column (acqcam_table, "FIELD_SC_Y", "pix", CPL_TYPE_DOUBLE);
    gravi_table_new_column (acqcam_table, "FIELD_SC_XERR", "pix", CPL_TYPE_DOUBLE);
    gravi_table_new_column (acqcam_table, "FIELD_SC_YERR", "pix", CPL_TYPE_DOUBLE);
    gravi_table_new_column (acqcam_table, "FIELD_FT_X", "pix", CPL_TYPE_DOUBLE);
    gravi_table_new_column (acqcam_table, "FIELD_FT_Y", "pix", CPL_TYPE_DOUBLE);
    gravi_table_new_column (acqcam_table, "FIELD_FT_XERR", "pix", CPL_TYPE_DOUBLE);
    gravi_table_new_column (acqcam_table, "FIELD_FT_YERR", "pix", CPL_TYPE_DOUBLE);
    gravi_table_new_column (acqcam_table, "FIELD_SCALE", "mas/pix", CPL_TYPE_DOUBLE);
    gravi_table_new_column (acqcam_table, "FIELD_SCALEERR", "mas/pix", CPL_TYPE_DOUBLE);
    gravi_table_new_column (acqcam_table, "FIELD_FIBER_DX", "pix", CPL_TYPE_DOUBLE);
    gravi_table_new_column (acqcam_table, "FIELD_FIBER_DY", "pix", CPL_TYPE_DOUBLE);
    gravi_table_new_column (acqcam_table, "FIELD_FIBER_DXERR", "pix", CPL_TYPE_DOUBLE);
    gravi_table_new_column (acqcam_table, "FIELD_FIBER_DYERR", "pix", CPL_TYPE_DOUBLE);
    gravi_table_new_column (acqcam_table, "FIELD_STREHL", "ratio", CPL_TYPE_DOUBLE);

    /* Position of roof center on full frame */
    double roof_x[] = {274.4, 787.1, 1236.1, 1673.4};
    double roof_y[] = {242.3, 247.7, 225.8, 235.6};
    
    /* Position of single-field spot on full frame */
    double spot_x[] = {289. , 798.2, 1245.5, 1696.};
    double spot_y[] = {186.5, 187.5,  172.5,  178.};
    
    /* Default position angle of roof */
    double roof_pos[] = {38.49, 38.54, 38.76, 39.80};
    
    /* Static acqcam field calibration error */
    /* calibration from Frank, March 2017, using Marcel      	
    double calib_dx[] = {-0.3483, -1.0251, -0.5432, -0.2024} ;
    double calib_dy[] = { 0.3089, -0.5221, -0.2686, -0.3843} ; */
    /* calibration from Julien, using GJ65 observation with Sylvestre tracking on 2018-01-03, email 2018-01-10 
    double calib_dx[] = {-0.61, -1.69, -0.97,  0.00} ;
    double calib_dy[] = { 0.53, -1.07, -0.49, -0.47} ; */
    /* calibration from Oli, using GJ65 observation with Sylvestre tracking on 2018-01-03, email 2018-01-19/
    double calib_dx[] = {-0.63, -1.67, -0.97,  0.02} ;
    double calib_dy[] = { 0.55, -1.07, -0.49, -0.45} ; */
    /* FE: removing obsolete calib_dx, because now handled by ICS, or by 
       correcting fitsheader retrospectively */
    double calib_dx[] = {0,0,0,0} ;
    double calib_dy[] = {0,0,0,0} ;
    
    /* If sub-windowing, we read the sub-window size */
    cpl_size nsx = 512;
    cpl_size nsy = 512;
    if ( cpl_propertylist_has (header, "ESO DET1 FRAMES NX") ) {
      nsx = cpl_propertylist_get_int (header, "ESO DET1 FRAMES NX");
      nsy = cpl_propertylist_get_int (header, "ESO DET1 FRAMES NY");
    }
    
    /* Compute separation */
    double sobj_x = gravi_pfits_get_sobj_x (header);
    double sobj_y = gravi_pfits_get_sobj_y (header);
    
    /* Accumulated mapping offsets */
    double sobj_offx=0., sobj_offy=0, sobj_drho=0., sobj_dth=0.;
    if (cpl_propertylist_has(header, "ESO INS SOBJ OFFX")) {
        /* mapping / mosaicing blind offsets from otiginal acquisition
         * to current position */
        sobj_offx = cpl_propertylist_get_double (header, "ESO INS SOBJ OFFX");
        sobj_offy = cpl_propertylist_get_double (header, "ESO INS SOBJ OFFY");
        
        /* distance from acquired SC position to current SC position in mas */
        sobj_drho = sqrt(sobj_offx*sobj_offx+sobj_offy*sobj_offy);
        
        /* position angle on sky from acquired SC position to     */
        /* current SC position, neglecting anamorphism variations */
        sobj_dth = atan2(sobj_offx, sobj_offy)*180./M_PI;
    }
    CPLCHECK_MSG ("Cannot get separation");
    
    /* Recover position of originally acquired star, before any blind offset */
    double dx_in = sobj_x - sobj_offx;
    double dy_in = sobj_y - sobj_offy;
    double rho_in = sqrt(dx_in*dx_in + dy_in*dy_in);
    
    /* Force zero separation for SINGLE mode */
    if ( gravi_pfits_get_mode (header) == MODE_SINGLE) rho_in = 0.;
    CPLCHECK_MSG ("Error reading header information");
    
    /* Loop on tel */
    for (int tel = 0; tel < ntel; tel++) {
        char name[90];
        double rp=roof_pos[tel]; // default value of roof position angle
        cpl_size sx=tel*512+1, sy=1;
        cpl_msg_info (cpl_func, "Compute field position for beam %i", tel+1);
        
        /* Read roof position angle */
        sprintf (name, "ESO INS DROTOFF%d", tel + 1);
        if ( cpl_propertylist_has (header, name) ) {
            rp = cpl_propertylist_get_double(header, name);
            CPLCHECK ("Cannot get rotation");
        }
        
        /* Approx. position angle of the binary, left from top */ 
        double approx_PA = 270.-rp;
        
        /* Get the telescope name and ID */
        const char * telname = gravi_conf_get_telname (tel, header);
        
        /* Check telescope name */
        if (!telname) cpl_msg_error (cpl_func, "Cannot read the telescope name");
        cpl_ensure_code (telname, CPL_ERROR_ILLEGAL_INPUT);
        
        /* Hardcoded approx. plate-scale in mas/pix. */
        double scale = 0.0;
        if (telname[0] == 'U') {
            scale = 18.;
        } else if (telname[0] == 'A') {
            if (tel == 0) {
                scale = 76.8;
            } else if (tel == 1) {
                scale = 78.0;
            } else if (tel == 2) {
                scale = 77.0;
            } else if (tel == 3) {
                scale = 84.6;
            }
        }
        
        /* Position of the fibres */
        double fiber_xft=0.;
        double fiber_yft=0.;
        double fiber_xsc=0.;
        double fiber_ysc=0.;
        sprintf (name, "ESO ACQ FIBER FT%dX", tel + 1);
        if (cpl_propertylist_has(header, name)) {
            fiber_xft=cpl_propertylist_get_double(header, name);
            sprintf (name, "ESO ACQ FIBER FT%dY", tel + 1);
            fiber_yft=cpl_propertylist_get_double(header, name);
            sprintf (name, "ESO ACQ FIBER SC%dX", tel + 1);
            fiber_xsc=cpl_propertylist_get_double(header, name);
            sprintf (name, "ESO ACQ FIBER SC%dY", tel + 1);
            fiber_ysc=cpl_propertylist_get_double(header, name);
        }
        double fiber_ft_sc_x=fiber_xsc-fiber_xft;
        double fiber_ft_sc_y=fiber_ysc-fiber_yft;
        
        /* Get the North position angle on the camera */
        double fangle = gravi_pfits_get_fangle_acqcam (header, tel);
        CPLCHECK ("Cannot get rotation");
        
        /* Mapping/mosaicing offset on acq cam axes, in mas, */
        /* neglecting amnamorphism variations */
        double sobj_offx_cam = sobj_drho * sin((fangle+sobj_dth)/180.*M_PI);
        double sobj_offy_cam = sobj_drho * cos((fangle+sobj_dth)/180.*M_PI);
        
        /* If sub-windowing, we read the sub-window start for field */
        if ( nsx != 512 ) {
            sprintf (name, "ESO DET1 FRAM%d STRX", tel + 1);
            sx = cpl_propertylist_get_int (header, name);
            
            sprintf (name, "ESO DET1 FRAM%d STRY", tel + 1);
            sy = cpl_propertylist_get_int (header, name);
            
            CPLCHECK_MSG ("Cannot get sub-windowing parameters");
        }
        
        cpl_msg_debug (cpl_func,"sub-window field %i sx= %lld sy = %lld", tel, sx, sy);
        
        /*--------------------------------------------------*/
        /* Guess position of SC and FT targes in mean image */
        /*--------------------------------------------------*/
        
        double xFT, yFT, xSC, ySC;
        
        if (rho_in == 0.) {
            /* TODO: close dual-field */
            /* Single-field case */
            /* Simply shift the best spot from full frame to cut-out */
            xFT = spot_x[tel] - sx + nsx*tel + 1;
            yFT = spot_y[tel] - sy + 1;
            xSC = xFT;
            ySC = yFT;
        } else {
            /* Pixel position of roof center on cut-out frame */
            /* Shift from full frame to cut-out */
            double cutout_roof_x = roof_x[tel] - sx + nsx*tel + 1;
            double cutout_roof_y = roof_y[tel] - sy + 1;
            
            cpl_msg_info (cpl_func, "    ROOF_X = %f", cutout_roof_x);
            cpl_msg_info (cpl_func, "    ROOF_Y = %f", cutout_roof_y);
            
            /* Approx pixel offset from SC to FT, divided by 2 */
            double approx_dx=0.5*rho_in*sin(approx_PA*M_PI/180.)/scale;
            double approx_dy=0.5*rho_in*cos(approx_PA*M_PI/180.)/scale;
            
            /* Expected position of the two stars */
            xFT = cutout_roof_x - approx_dx ;
            yFT = cutout_roof_y - approx_dy ;
            xSC = cutout_roof_x + approx_dx ;
            ySC = cutout_roof_y + approx_dy ;
        }
        
        cpl_msg_info (cpl_func, "guess SC_X = %f", xSC);
        cpl_msg_info (cpl_func, "guess SC_Y = %f", ySC);
        cpl_msg_info (cpl_func, "guess FT_X = %f", xFT);
        cpl_msg_info (cpl_func, "guess FT_Y = %f", yFT);
        
        /*------------------------------------------*/
        /* Measure SC target position in mean image */
        /*------------------------------------------*/
        
        /* Box size */
        /* Optimal size has been determined empirically. */
        /* Too small value (15) will sometimes miss even a bright target. */
        /* Too large value (25) will often pick a wrong star (or backround noise) */
        cpl_size size = 20;
        
        double xSCguess=xSC, ySCguess=ySC, xFTguess=xFT, yFTguess=yFT;
        double ex, ey;
        double qc_val=0.;
        
        gravi_acq_fit_gaussian (mean_img, &xSC, &ySC, &ex, &ey, size);
        CPLCHECK_MSG("Error fitting SC");
        
        /*--------------------------------------------------------*/
        /* Add QC parameters for SC target position in mean image */
        /*--------------------------------------------------------*/
        
        sprintf (qc_name, "ESO QC ACQ FIELD%i SC_X", tel+1);
        if (xSC==0.) {
            /* Fitting failed: put QC to 0., reset xSC to gues value */
            qc_val = 0.;
            xSC = xSCguess;
        } else {
            /* Fiting succeeded: shift into full frame */
            qc_val = xSC + sx - 1 - nsx*tel;
        }
        cpl_msg_info (cpl_func, "%s = %f", qc_name, qc_val);
        cpl_propertylist_update_double (o_header, qc_name, qc_val);
        cpl_propertylist_set_comment (o_header, qc_name, "[pixel] position in mean image");
        
        sprintf (qc_name, "ESO QC ACQ FIELD%i SC_Y", tel+1);
        if (ySC==0.) {
            /* Fitting failed: put QC to 0., reset ySC to gues value */
            qc_val = 0.;
            ySC = ySCguess;
        } else {
            /* Fiting succeeded: shift into full frame */
            qc_val =  ySC + sy -1;
        }
        cpl_msg_info (cpl_func, "%s = %f", qc_name, qc_val);
        cpl_propertylist_update_double (o_header, qc_name, qc_val);
        cpl_propertylist_set_comment (o_header, qc_name, "[pixel] position in mean image");
        
        /*------------------------------------------*/
        /* Measure FT target position in mean image */
        /*------------------------------------------*/
        
        if (rho_in != 0.) {
            /* Detec FT */
            gravi_acq_fit_gaussian (mean_img, &xFT, &yFT, &ex, &ey, size);
            CPLCHECK_MSG("Error fitting FT");
        } else {
            xFT=xSC;
            yFT=ySC;
        }
        
        /*--------------------------------------------------------*/
        /* Add QC parameters for SC target position in mean image */
        /*--------------------------------------------------------*/
        
        sprintf (qc_name, "ESO QC ACQ FIELD%i FT_X", tel+1);
        if (xFT==0.) {
            /* Fitting failed: put QC to 0., reset xFT to gues value */
            qc_val = 0.;
            xFT = xFTguess;
        } else {
            /* Fiting succeeded: shift into full frame */
            qc_val = xFT + sx - 1 - nsx*tel;
        }
        cpl_msg_info (cpl_func, "%s = %f", qc_name, qc_val);
        cpl_propertylist_update_double (o_header, qc_name, qc_val);
        cpl_propertylist_set_comment (o_header, qc_name, "[pixel] position in mean image");
        
        sprintf (qc_name, "ESO QC ACQ FIELD%i FT_Y", tel+1);
        if (yFT==0.) {
            /* Fitting failed: put QC to 0., reset yFT to gues value */
            qc_val = 0.;
            yFT = yFTguess;
        } else {
            /* Fiting succeeded: shift into full frame */
            qc_val =  yFT + sy -1;
        }
        cpl_msg_info (cpl_func, "%s = %f", qc_name, qc_val);
        cpl_propertylist_update_double (o_header, qc_name, qc_val);
        cpl_propertylist_set_comment (o_header, qc_name, "[pixel] position in mean image");
        
        /*---------------------------------------*/
        /* If in dual field, measure plate scale */
        /*---------------------------------------*/
        
        if (rho_in != 0.) {
            sprintf (qc_name, "ESO QC ACQ FIELD%i SCALE", tel+1);
            double sep = sqrt((ySC-yFT)*(ySC-yFT)+(xSC-xFT)*(xSC-xFT));
            double pscale = sep ? rho_in/sep : 0.;
            qc_val=pscale;
            cpl_msg_info (cpl_func, "%s = %f", qc_name, qc_val);
            cpl_propertylist_update_double (o_header, qc_name, qc_val);
            cpl_propertylist_set_comment (o_header, qc_name,
                                        "[mas/pixel] plate-scale in the "
                                        "FT-SC direction");
        }
        
        /*------------------------------------------------*/
        /* If in dual field, measure fibre position error */
        /*------------------------------------------------*/
        
        if (rho_in != 0.) {
            /* Error in SC fibre positioning */
            /* The three terms are */
            /*  - offset from FT target as detected to original SC */
            /*    target as detected; */
            /*  - blind offset command from mapping template projected */
            /*    on acqisition camera; */
            /*  - offset from FT fiber to SC fiber. */
            sprintf (qc_name, "ESO QC ACQ FIELD%i SC_FIBER_DX", tel+1);
            qc_val = 0;
            if (scale) {
                qc_val = (xSC-xFT) + sobj_offx_cam/scale - fiber_ft_sc_x - calib_dx[tel];
            }
            cpl_msg_info (cpl_func, "%s = %f", qc_name, qc_val);
            cpl_propertylist_update_double (o_header, qc_name, qc_val);
            cpl_propertylist_set_comment (o_header, qc_name,
                                        "[pixel] dx from SC fiber to "
                                        "SC object");
            
            sprintf (qc_name, "ESO QC ACQ FIELD%i SC_FIBER_DY", tel+1);
            qc_val = 0;
            if (scale) {
                qc_val = (ySC-yFT) + sobj_offy_cam/scale - fiber_ft_sc_y - calib_dy[tel];
            }
            cpl_msg_info (cpl_func, "%s = %f", qc_name, qc_val);
            cpl_propertylist_update_double (o_header, qc_name, qc_val);
            cpl_propertylist_set_comment (o_header, qc_name,
                                        "[pixel] dx from SC fiber to "
                                        "SC object");
        }
        
        /*-----------------------------*/
        /* Average Strehl computation  */
        /*-----------------------------*/
        
        double max_on_average, strehl_on_average;
        
        /* Measure strehl and maximum on averaged image */
        gravi_acq_measure_strehl(mean_img, xFT, yFT, scale, &strehl_on_average, header);
        gravi_acq_measure_max(mean_img, xFT, yFT, 15, &max_on_average);
        
        /* Update Strehl QC */
        sprintf (qc_name, "ESO QC ACQ FIELD%i STREHL", tel+1);
        cpl_msg_info (cpl_func, "%s = %f", qc_name, strehl_on_average);
        cpl_propertylist_update_double (o_header, qc_name, strehl_on_average);
        cpl_propertylist_set_comment (o_header, qc_name, "Average Strehl from stacked AcqCam images");
        
        /*----------------------------*/
        /* Now process frame by frame */
        /*----------------------------*/
        
        /* Loop on all images */
        for (cpl_size row = 0; row < nrow; row++) {
            if (row %10 == 0 || row == (nrow-1))
                cpl_msg_info_overwritable (cpl_func, "Fit image %lld over %lld", row+1, nrow);
            
            /* Get data */
            cpl_image * img = cpl_imagelist_get (acqcam_imglist, row);
            CPLCHECK_MSG("Error getting image");
            
            /*-------------------------------------------------*/
            /* SC target position computation of current frame */
            /*-------------------------------------------------*/
            double xsc = xSC, ysc = ySC, exsc=0., eysc=0.;
            gravi_acq_fit_gaussian (img, &xsc, &ysc, &exsc, &eysc, size);
            CPLCHECK_MSG("Error fitting SC");          
            /* Shift back positions to full frame */
            if (xsc != 0.) xsc += sx - 1 - nsx*tel;
            if (ysc != 0.) ysc += sy - 1;
            cpl_table_set (acqcam_table, "FIELD_SC_X", row*ntel+tel, xsc);
            cpl_table_set (acqcam_table, "FIELD_SC_Y", row*ntel+tel, ysc);
            cpl_table_set (acqcam_table, "FIELD_SC_XERR", row*ntel+tel, exsc);
            cpl_table_set (acqcam_table, "FIELD_SC_YERR", row*ntel+tel, eysc);
            CPLCHECK_MSG("Error setting SC columns");
            
            /*-------------------------------------------------*/
            /* FT target position computation of current frame */
            /*-------------------------------------------------*/
            double xft = xFT, yft = yFT, exft=0., eyft=0.;
            if (rho_in != 0.) {
                gravi_acq_fit_gaussian (img, &xft, &yft, &exft, &eyft, size);
                CPLCHECK_MSG("Error fitting FT");
                /* Shift back positions to full frame */
                if (xft != 0.) xft += sx - 1 - nsx*tel;
                if (yft != 0.) yft += sy - 1;
            } else {
                xft=xsc;
                yft=ysc;
                exft=exsc;
                eyft=eysc;
            }
            cpl_table_set (acqcam_table, "FIELD_FT_X", row*ntel+tel, xft);
            cpl_table_set (acqcam_table, "FIELD_FT_Y", row*ntel+tel, yft);
            cpl_table_set (acqcam_table, "FIELD_FT_XERR", row*ntel+tel, exft);
            cpl_table_set (acqcam_table, "FIELD_FT_YERR", row*ntel+tel, eyft);
            CPLCHECK_MSG("Error setting FT column");

            /*------------------------------------------*/
            /* Plate scale computation of current frame */
            /*------------------------------------------*/
            double ft_sc_x = xsc - xft;
            double ft_sc_y = ysc - yft;
            double eft_sc_x = sqrt(exsc*exsc+exft*exft);
            double eft_sc_y = sqrt(eysc*eysc+eyft*eyft);
            double sep = sqrt(ft_sc_x*ft_sc_x+ft_sc_y*ft_sc_y);
            double pscale = sep ? rho_in/sep : 0.;
            double escale = 0.;
            if (sep) {
                escale = rho_in/(sep*sep*sep)*
                (ft_sc_x*eft_sc_x+ft_sc_y*eft_sc_y) ;
            }
            cpl_table_set (acqcam_table, "FIELD_SCALE", row*ntel+tel, pscale);
            cpl_table_set (acqcam_table, "FIELD_SCALEERR", row*ntel+tel, escale);
            
            /*---------------------------------------------*/
            /* Error in fibre positioning of current frame */
            /*---------------------------------------------*/
            /* The three terms are */
            /*  - offset from FT target as detected to original SC */
            /*  target as detected; */
            /*  - blind offset command from mapping template projected */
            /*  on acqisition camera; */
            /*  - offset from FT fiber to SC fiber. */
            double corrx=0, corry=0., ecorrx=0., ecorry=0.;
            if (pscale) {
                corrx = ft_sc_x + sobj_offx_cam/pscale - fiber_ft_sc_x - calib_dx[tel];
                corry = ft_sc_y + sobj_offy_cam/pscale - fiber_ft_sc_y - calib_dy[tel];
                double tmp = escale/(pscale*pscale);
                tmp *= tmp;
                ecorrx = sqrt(eft_sc_x*eft_sc_x + eft_sc_y*eft_sc_y + sobj_offx_cam*sobj_offx_cam*tmp);
                ecorry = sqrt(eft_sc_x*eft_sc_x + eft_sc_y*eft_sc_y + sobj_offy_cam*sobj_offy_cam*tmp);
            }
            cpl_table_set (acqcam_table, "FIELD_FIBER_DX", row*ntel+tel, corrx);
            cpl_table_set (acqcam_table, "FIELD_FIBER_DY", row*ntel+tel, corry);
            cpl_table_set (acqcam_table, "FIELD_FIBER_DXERR", row*ntel+tel, ecorrx);
            cpl_table_set (acqcam_table, "FIELD_FIBER_DYERR", row*ntel+tel, ecorry);
            
            /*-------------------------------------*/
            /* Strehl computation of current frame */
            /*-------------------------------------*/
            double max_on_frame;
            gravi_acq_measure_max(img, xFT, yFT, 15, &max_on_frame);
            cpl_table_set (acqcam_table, "FIELD_STREHL", row*ntel+tel, strehl_on_average*(max_on_frame/max_on_average) );
            
        } /* End loop on images */

        /* Add some QC */
        double sc_std_x = gravi_table_get_column_std (acqcam_table, "FIELD_SC_X", tel, ntel);
        sprintf (qc_name, "ESO QC ACQ FIELD%i SC_X STD", tel+1);
        cpl_propertylist_update_double (o_header, qc_name, sc_std_x);
        cpl_propertylist_set_comment (o_header, qc_name, "[pix] Std of field position of SC");

        double sc_std_y = gravi_table_get_column_std (acqcam_table, "FIELD_SC_Y", tel, ntel);
        sprintf (qc_name, "ESO QC ACQ FIELD%i SC_Y STD", tel+1);
        cpl_propertylist_update_double (o_header, qc_name, sc_std_y);
        cpl_propertylist_set_comment (o_header, qc_name, "[pix] Std of field position of SC");

        double ft_std_x = gravi_table_get_column_std (acqcam_table, "FIELD_FT_X", tel, ntel);
        sprintf (qc_name, "ESO QC ACQ FIELD%i FT_X STD", tel+1);
        cpl_propertylist_update_double (o_header, qc_name, ft_std_x);
        cpl_propertylist_set_comment (o_header, qc_name, "[pix] Std of field position of FT");

        double ft_std_y = gravi_table_get_column_std (acqcam_table, "FIELD_FT_Y", tel, ntel);
        sprintf (qc_name, "ESO QC ACQ FIELD%i FT_Y STD", tel+1);
        cpl_propertylist_update_double (o_header, qc_name, ft_std_y);
        cpl_propertylist_set_comment (o_header, qc_name, "[pix] Std of field position of FT");
        
    } /* End loop on tel */
    
    gravi_msg_function_exit(1);
    return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Reduce the ACQ camera images
 *  
 * @param output_data:  The output gravi_data where the OI_VIS_ACQ table
 *                      will be created, with ndit * ntel rows.
 * @param input_data:   The input gravi_data here the ACQ imagelist is
 *                      read.
 * \exception CPL_ERROR_NULL_INPUT input data is missing
 * 
 * The routine process the PUPIL sensor. It creates a table with
 * the columns PUPIL_NSPOT (number of detected spots), PUPIL_R (rotation angle)
 * of telescope diode, PUPIL_X, PUPIL_Y, PUPIL_Z (shifts in [m]) (@c gravi_acqcam_pupil).
 * The routine process also the FIELD sensor. It creates a table with the position
 * of the FT and SC target, the fiber shift DX, DY, SCALE and the STREHL (@c gravi_acqcam_field).
 * The TIME in [us] is also stored.
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_reduce_acqcam (gravi_data * output_data,
                                    gravi_data * input_data)
{
    gravi_msg_function_start(1);
    cpl_ensure_code (output_data, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code (input_data,  CPL_ERROR_NULL_INPUT);
    
    /* Check if extension exist */
    if (!gravi_data_has_extension (input_data, GRAVI_IMAGING_DATA_ACQ_EXT)) {
        gravi_msg_warning (cpl_func, "Cannot reduce the ACQCAM, not data");
        return CPL_ERROR_NONE;
    }

    /* Get the data and header */
    cpl_propertylist * header, * o_header;
    header = gravi_data_get_header (input_data);
    o_header = gravi_data_get_header (output_data);
    
    cpl_imagelist * acqcam_imglist;
    acqcam_imglist = gravi_data_get_cube (input_data, GRAVI_IMAGING_DATA_ACQ_EXT);
    CPLCHECK_MSG ("Cannot get data or header");

    /* Build the table */
    cpl_size ntel = 4;
    cpl_size nrow = cpl_imagelist_get_size (acqcam_imglist);
    
    cpl_table * acqcam_table;
    acqcam_table = cpl_table_new (nrow * ntel);

    /* Time column shall contain the time from PCR.ACQ.START in [us] */
    cpl_table_new_column (acqcam_table, "TIME", CPL_TYPE_INT);
    cpl_table_set_column_unit (acqcam_table, "TIME", "us");

    /* Loop on DIT in cube to fill the TIME column 
     * same value for all 4 beams*/
    for (cpl_size row = 0; row < nrow; row++) {
        double time = gravi_pfits_get_time_acqcam (header, row);
        for (int tel = 0; tel < ntel; tel ++)
            cpl_table_set (acqcam_table, "TIME", row*ntel+tel, time);
    }
    
    /* Compute mean image */
    cpl_image * mean_img = cpl_imagelist_collapse_create (acqcam_imglist);

    /* Compute FIELD columns */
    gravi_acqcam_field (mean_img, acqcam_imglist, header,
                        acqcam_table, o_header);
	CPLCHECK_MSG ("Cannot reduce field images");
    
    
    /* Compute PUPIL columns */
    gravi_acqcam_pupil (mean_img, acqcam_imglist, header,
                        acqcam_table, o_header);
	CPLCHECK_MSG ("Cannot reduce pupil images");
    
    /* Add this output table in the gravi_data */
	gravi_data_add_img (output_data, NULL, GRAVI_IMAGING_DATA_ACQ_EXT, mean_img);
	CPLCHECK_MSG ("Cannot add acqcam_table in data");
    
	gravi_data_add_table (output_data, NULL, GRAVI_OI_VIS_ACQ_EXT, acqcam_table);
	CPLCHECK_MSG ("Cannot add acqcam_table in data");
    
    gravi_msg_function_exit(1);
    return CPL_ERROR_NONE;   
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Convert z_shift from [pixel] to [meters]
 *        Formula extracted from gvacqPupilTracker.c
 */
/*----------------------------------------------------------------------------*/

double gravi_acqcam_z2meter (double PositionPixels)
{
    double f_PT      = 14e-3;    /* pupil tracker lenslet FL*/
    double f_lens    = 467e-3;   /* folding optics lens FL */
    double Llambda   = 1.2e-6;   /* laser diode wavelength */
    double D_beam    = 18e-3;    /* meter */
    double D_pixel   = 18e-6;
    double D_AT      = 1.8;      /* m */
    double D_lenslet = 2 * 1.015e-3;
    
    double longDef;
    longDef = 8 * (f_PT / D_lenslet) * (f_PT / D_lenslet) * 3.5 * D_pixel *
        D_beam / (f_PT * D_lenslet) * Llambda / CPL_MATH_2PI * PositionPixels;
    
    return f_lens * f_lens * longDef / (f_PT + longDef) / f_PT * (D_AT / D_lenslet);
}

/**@}*/
