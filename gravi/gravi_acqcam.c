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
 * @defgroup gravity_acqcam  TBD
 */
/**@{*/

/*-----------------------------------------------------------------------------
                                   Includes
 -----------------------------------------------------------------------------*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#define _XOPEN_SOURCE
#include <cpl.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "gravi_dfs.h"
#include "gravi_data.h"
#include "gravi_pfits.h"

#include "gravi_cpl.h"
#include "gravi_utils.h"

#include "gravi_acqcam.h"

/*-----------------------------------------------------------------------------
                               Private prototypes
 -----------------------------------------------------------------------------*/

int gravi_acqcam_isblink (cpl_imagelist * imglist, cpl_size pos);

cpl_error_code gravi_acqcam_get_pup_ref (cpl_propertylist * header,
                                         cpl_size tel,
                                         cpl_vector * pupref);

cpl_error_code gravi_acqcam_get_diode_ref (cpl_propertylist * header,
                                           cpl_size tel,
                                           cpl_vector * output);

static int gravi_acqcam_spot (const double x_in[], const double v[], double *result);

cpl_error_code gravi_acqcam_spot_imprint (cpl_image * img, cpl_vector * a);

static int gravi_acqcam_spot_dfda (const double x_in[], const double v[], double result[]);

cpl_error_code gravi_acqcam_fit_spot (cpl_image * img, cpl_size ntry,
                                      cpl_vector * a, 
				      int fitAll,
				      int * nspot);

double gravi_acqcam_z2meter (double PositionPixels);

/* This global variable optimises the computation
 * of partial derivative on fitted parameters */
const extern int * GRAVI_LVMQ_FREE;
const int * GRAVI_LVMQ_FREE = NULL;

/* Number of parameters in the model 'gravi_acqcam_spot' 
 * And position of parameters */
#define GRAVI_SPOT_NA    29
#define GRAVI_SPOT_SUB   0
#define GRAVI_SPOT_ANGLE 8
#define GRAVI_SPOT_SCALE 9
#define GRAVI_SPOT_DIODE 10
#define GRAVI_SPOT_FWHM  12
#define GRAVI_SPOT_FLUX  13

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

    cpl_propertylist * data_header;
    data_header = gravi_data_get_header (input_data);
    CPLCHECK_MSG ("Cannot get data or header");

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
    
        int blink = gravi_acqcam_isblink (imglist, 0) == 1 ? 0 : 1;
        
        /* Remove the blink only for pupil */
        for (cpl_size row = 0; row < nrow; row ++) {
            gravi_image_subtract_window (cpl_imagelist_get (imglist, row),
                                         cpl_imagelist_get (imglist, blink),
                                         1, ury, nx, ny, 1, ury);
            CPLCHECK_MSG ("Cannot remove blinked pupil");
            if (row == blink && row < nrow-2) blink +=2;
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
inline double sin1 (double x)
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
inline double exp1 (double x) {
  x = 1.0 + x / 256.0;
  x *= x; x *= x; x *= x; x *= x;
  x *= x; x *= x; x *= x; x *= x;
  return x;
}

/* Compute the 4 diode positions from {angle, scaling, dx, dy} 
 * The model assume the 4 diodes form a rectangle centered on
 * the pupil. Hence this model has a 180deg degeneracy */
inline int gravi_acqcam_xy_diode (const double v[], double *xd, double *yd)
{
    /* Angle */
    double ang = v[GRAVI_SPOT_ANGLE] * CPL_MATH_RAD_DEG;
    double sang = sin1 (ang) * v[GRAVI_SPOT_SCALE];
    double cang = sin1 (ang + CPL_MATH_PI_2) * v[GRAVI_SPOT_SCALE];

    double dx = v[GRAVI_SPOT_DIODE+0];
    double dy = v[GRAVI_SPOT_DIODE+1];
    
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
inline int gravi_acqcam_xy_sub (const double v[], double *xsub, double *ysub)
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
 * @brief Read the diode position from header into the vector output. Assume
 *        the four diodes form a rectangle centered on the pupil center.
 * 
 *        output[8]   = rotation  [deg], set to 0.0
 *        output[9]   = scale  [pix/m]
 *        output[10]  = dx [m]
 *        output[11]  = dy [m]
 * 
 * @param header:   input header
 * @param tel:      requested beam (0..3)
 * @param output:   output vector, shall be already allocated
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
    // int telid = atoi (telname+2) - 1;
    CPLCHECK ("Cannot get telescope name");

    /* Hardcoded theoretical positions in mm */

    /* If UTs or ATs, select scaling and position */
    if (telname[0] == 'U') {
        cpl_vector_set (output, GRAVI_SPOT_ANGLE,  0.0);
        cpl_vector_set (output, GRAVI_SPOT_SCALE, 16.225);
        cpl_vector_set (output, GRAVI_SPOT_DIODE+0, 0.363);
        cpl_vector_set (output, GRAVI_SPOT_DIODE+1, 0.823);
    } else if (telname[0] == 'A') {
        cpl_vector_set (output, GRAVI_SPOT_ANGLE,  0.0);
        cpl_vector_set (output, GRAVI_SPOT_SCALE, 73.0154);
        cpl_vector_set (output, GRAVI_SPOT_DIODE+0, 0.122);
        cpl_vector_set (output, GRAVI_SPOT_DIODE+1, 0.158);
    } else 
        return cpl_error_set_message (cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                      "Cannot get telescope name");
        
    gravi_msg_function_exit(0);
    return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Read the sub-aperture position for the pupil sensor,
 *        and re-arrange them into the output vector
 *
 *        output[0] = x0+x1+x2+x3   (center of sub-apertures)
 *        output[1] = x0-x1+x2-x3
 *        output[2] = x0+x1-x2-x3
 *        output[3] = x0-x1-x2+x3
 *        output[4..7] = same for y
 *
 * @param header:   input header
 * @param tel:      requested beam (0..3)
 * @param output:   output vector, shall be already allocated
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
 * @brief Fit a pupil spot pattern into an image. The global minimum is found
 *        first with a fit which is, by-design, a correlation with the brightest
 *        pixels of the image. The number of random starting point is given by
 *        the ntry parameter. If ntry==1, the starting is kept unmodified.
 *        Then 10x10 pixels around each expected spot are extracted and fit
 *        with true Gaussian (free FWHM and amplitude).
 * 
 * @param img:    input image
 * @param ntry:   number of random starting point
 * @param a:      vector of parameter, modified in-place
 * @param nspot:  is filled with the number of detected spots.
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
    cpl_size nx = cpl_image_get_size_x (img);
    cpl_size ny = cpl_image_get_size_y (img);

    cpl_size x0 = cpl_vector_get (a, GRAVI_SPOT_SUB+0);
    cpl_size y0 = cpl_vector_get (a, GRAVI_SPOT_SUB+4);
    CPLCHECK_MSG ("Cannot get values valid");

    /* Compute RMS in the central region */
    double RMS = gravi_image_get_noise_window (img, x0-25, y0-25, x0+25, y0+25);

    /* The image is surely empty */
    if (RMS == 0) { *nspot = 0; return CPL_ERROR_NONE;}

    /*
     * Coarse: correlation with a re-bin image
     */

    /* To lower the number of point, we extract a window of 100x100
     * around the center of sub-apertures, rebin with 3x3 pixels */
    cpl_size nw = 100, n_mean = 3;
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
	for (int i=-1; i<=1; i++) {
	  for (int j=-1; j<=1; j++) {
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
	cpl_vector_set (sy_vector, x*nc+y, 1.0);
	CPLCHECK_MSG ("Cannot fill matrix/vector");
      }
    } /* End loop on re-sampled pixels*/

    /* Normalize for numerical stability. */
    cpl_vector_divide_scalar (y_vector, RMS);

    /* Output for global minimisation */
    cpl_vector * a_start = cpl_vector_duplicate (a);
    cpl_vector * a_tmp = cpl_vector_duplicate (a);
    double chisq_final = 1e10;
    srand(1);

    /* Fit sub-aperture mean position; and diode rotation */
    const int ia_global[] = {1,0,0,0, 1,0,0,0, 1,0,0,0, 0,
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    GRAVI_LVMQ_FREE = ia_global;

    /* Loop on various starting points */
    for (cpl_size try = 0; try < ntry; try++) {

        /* Move starting point in position (+-10pix) and angle (entire circle) */
        cpl_vector_copy (a_tmp, a_start);
        if (try > 0) {
            cpl_vector_set (a_tmp, GRAVI_SPOT_SUB+0, cpl_vector_get (a_tmp, GRAVI_SPOT_SUB+0) + (rand()%20) - 10);
            cpl_vector_set (a_tmp, GRAVI_SPOT_SUB+4, cpl_vector_get (a_tmp, GRAVI_SPOT_SUB+4) + (rand()%20) - 10);
            cpl_vector_set (a_tmp, GRAVI_SPOT_ANGLE, cpl_vector_get (a_tmp, GRAVI_SPOT_ANGLE) + (rand()%180));
        }

        /* Set the fwhm to 6 and amplitude to 1.0, to force
         * a pseudo-correlation with large capture range */
        cpl_vector_set (a_tmp, GRAVI_SPOT_FWHM, 6.*6.);
        for (int d=0;d<16;d++) cpl_vector_set (a_tmp, GRAVI_SPOT_FLUX+d, 1.0);

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

    /* Impose FWHM to a realist value in [pix**2] */
    cpl_vector_set (a, GRAVI_SPOT_FWHM, 2.3*2.3);

    /* Fit all sub-aperture position; rotation and scaling of diodes;
     * and individual intensities of spots */
    int F = fitAll;
    const int ia_fine[] = {1,F,1,F, 1,1,F,F, 1,F,0,0, 0,
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
/**
 * @brief Reduce the ACQ camera images
 *  
 * @param output_data:  The output gravi_data where the OI_VIS_ACQ table
 *                      will be created, with ndit * ntel rows.
 * @param input_data:   The input gravi_data here the ACQ imagelist is
 *                      read.
 * 
 * The routine only process the PUPIL sensor so far. It creates a table with
 * the columns PUPIL_NSPOT (number of detected spots), PUPIL_R (rotation angle)
 * of telescope diode, PUPIL_X, PUPIL_Y, PUPIL_Z (shifts in [m]).
 * The TIME in [us] is also stored.
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_reduce_acqcam (gravi_data * output_data,
                                    gravi_data * input_data)
{
    gravi_msg_function_start(1);
    cpl_ensure_code (output_data, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code (input_data,  CPL_ERROR_NULL_INPUT);
    
    char qc_name[100];

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
    
    /* 
     * Compute TIME column
     */

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
    
    /* 
     * Compute PUPIL columns
     */

    /* Pupil positions (or we use array of 3)  */
    cpl_table_new_column (acqcam_table, "PUPIL_NSPOT", CPL_TYPE_INT);
    cpl_table_new_column (acqcam_table, "PUPIL_X", CPL_TYPE_DOUBLE);
    cpl_table_new_column (acqcam_table, "PUPIL_Y", CPL_TYPE_DOUBLE);
    cpl_table_new_column (acqcam_table, "PUPIL_Z", CPL_TYPE_DOUBLE);
    cpl_table_new_column (acqcam_table, "PUPIL_R", CPL_TYPE_DOUBLE);
    cpl_table_set_column_unit (acqcam_table, "PUPIL_X", "m");
    cpl_table_set_column_unit (acqcam_table, "PUPIL_Y", "m");
    cpl_table_set_column_unit (acqcam_table, "PUPIL_Z", "m");
    cpl_table_set_column_unit (acqcam_table, "PUPIL_R", "deg");
    
    /* Compute mean image */
    cpl_image * mean_img = cpl_imagelist_collapse_create (acqcam_imglist);
    int nspot = 0;
    
    /* Loop on tel */
    for (int tel = 0; tel < ntel; tel++) {
        cpl_msg_info (cpl_func, "Compute pupil position for beam %i", tel+1);

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
                
        double scale = cpl_vector_get (a_final,GRAVI_SPOT_SCALE);
        double xpos = cpl_vector_get (a_final,GRAVI_SPOT_SUB+0) -
                      cpl_vector_get (a_start,GRAVI_SPOT_SUB+0);
        double ypos = cpl_vector_get (a_final,GRAVI_SPOT_SUB+4) -
                      cpl_vector_get (a_start,GRAVI_SPOT_SUB+4);

        /* Add best position as a cross in image */
        gravi_acqcam_spot_imprint (mean_img, a_final);
        
        /* Add QC parameters */
        sprintf (qc_name, "ESO QC ACQ PUP%i NSPOT", tel+1);
        cpl_msg_info (cpl_func, "%s = %i", qc_name, nspot);
        cpl_propertylist_update_int (o_header, qc_name, nspot);
        cpl_propertylist_set_comment (o_header, qc_name, "nb. of pupil spot in ACQ");

        sprintf (qc_name, "ESO QC ACQ PUP%i ANGLE", tel+1);
        cpl_msg_info (cpl_func, "%s = %f", qc_name, cpl_vector_get (a_final,GRAVI_SPOT_ANGLE));
        cpl_propertylist_update_double (o_header, qc_name, cpl_vector_get (a_final,GRAVI_SPOT_ANGLE));
        cpl_propertylist_set_comment (o_header, qc_name, "[deg] diode angle on ACQ");

        sprintf (qc_name, "ESO QC ACQ PUP%i SCALE", tel+1);
        cpl_msg_info (cpl_func, "%s = %f", qc_name, scale);
        cpl_propertylist_update_double (o_header, qc_name, scale);
        cpl_propertylist_set_comment (o_header, qc_name, "[pix/m] diode scale on ACQ");

        sprintf (qc_name, "ESO QC ACQ PUP%i FWHM", tel+1);
        cpl_msg_info (cpl_func, "%s = %f", qc_name, sqrt (cpl_vector_get (a_final,GRAVI_SPOT_FWHM)));
        cpl_propertylist_update_double (o_header, qc_name, sqrt (cpl_vector_get (a_final,GRAVI_SPOT_FWHM)));
        cpl_propertylist_set_comment (o_header, qc_name, "[pix] spot fwhm in ACQ");

        sprintf (qc_name, "ESO QC ACQ PUP%i XPOS", tel+1);
        cpl_msg_info (cpl_func, "%s = %f", qc_name, xpos);
        cpl_propertylist_update_double (o_header, qc_name, xpos);
        cpl_propertylist_set_comment (o_header, qc_name, "[pix] pupil x-shift in ACQ");
        
        sprintf (qc_name, "ESO QC ACQ PUP%i YPOS", tel+1);
        cpl_msg_info (cpl_func, "%s = %f", qc_name, ypos);
        cpl_propertylist_update_double (o_header, qc_name, ypos);
        cpl_propertylist_set_comment (o_header, qc_name, "[pix] pupil y-shift in ACQ");
        
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
                
                /* Compute latteral shift [m] and rotation [deg] */
                double x_shift = cpl_vector_get (a_row, GRAVI_SPOT_SUB+0) / scale;
                double y_shift = cpl_vector_get (a_row, GRAVI_SPOT_SUB+4) / scale;
                double r_shift = cpl_vector_get (a_row, GRAVI_SPOT_ANGLE);
                
                /* Compute longitudinal shift [m] */
                double z_shift = -0.5 * ( cpl_vector_get (a_row, GRAVI_SPOT_SUB+2) +
                                          cpl_vector_get (a_row, GRAVI_SPOT_SUB+5));
                z_shift = gravi_acqcam_z2meter (z_shift);
                
                cpl_table_set (acqcam_table, "PUPIL_NSPOT", row*ntel+tel, nspot);
                cpl_table_set (acqcam_table, "PUPIL_X", row*ntel+tel, x_shift);
                cpl_table_set (acqcam_table, "PUPIL_Y", row*ntel+tel, y_shift);
                cpl_table_set (acqcam_table, "PUPIL_Z", row*ntel+tel, z_shift);
                cpl_table_set (acqcam_table, "PUPIL_R", row*ntel+tel, r_shift);
            }
            
            FREE (cpl_vector_delete, a_row);
        }

        FREE (cpl_vector_delete, a_start);
        FREE (cpl_vector_delete, a_final);
    } /* End loop on tel */
    
    /* 
     * Add this output table in the gravi_data 
     */
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
    
    double LongitudinalDefocusShift = 8 * (f_PT / D_lenslet) * (f_PT / D_lenslet) * 3.5 * D_pixel *
                                      D_beam / (f_PT * D_lenslet) * Llambda / CPL_MATH_2PI * PositionPixels;
    
    return f_lens * f_lens * LongitudinalDefocusShift / (f_PT + LongitudinalDefocusShift) / f_PT * (D_AT / D_lenslet);
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Return 1 if flux(pos) < flux(pos+1) * 
 */
/*----------------------------------------------------------------------------*/

int gravi_acqcam_isblink (cpl_imagelist * imglist, cpl_size pos)
{
    gravi_msg_function_start(1);
    cpl_ensure (imglist, CPL_ERROR_NULL_INPUT, -1);
    
    int blinkOFF = 1;

    /* Ensure two images */
    cpl_size nrow = cpl_imagelist_get_size (imglist);
    if (nrow == 1) return blinkOFF;

    /* Get size */
    cpl_image * img = cpl_imagelist_get (imglist,0);
    cpl_size nx = cpl_image_get_size_x (img);

    /* Part of the image to analyse */
    cpl_size llx, urx, lly, ury;
    if (nx < 1100) {
        llx =  1;
        urx =  1000;
        lly =  800;
        ury =  1000;
    } else {
        llx =  1;
        urx =  2048;
        lly =  1200;
        ury =  1536;
    }

    /* tell if the first image is blink on or off */
    double flux0 = cpl_image_get_flux_window (cpl_imagelist_get (imglist,pos), llx, lly, urx, ury);
    double flux1 = cpl_image_get_flux_window (cpl_imagelist_get (imglist,pos+1), llx, lly, urx, ury);
    if (flux0 > flux1) blinkOFF = 0;

    gravi_msg_function_exit(1);
    return blinkOFF;
}



/**@}*/

