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

cpl_error_code gravi_acqcam_get_pup_ref (cpl_propertylist * header, cpl_size tel,
                                         cpl_vector * pupref);

static int gravi_acqcam_spot (const double x_in[], const double v[], double *result);
cpl_error_code gravi_acqcam_spot_imprint (cpl_image * img, cpl_vector * a);
int gravi_acqcam_spot_count (cpl_image * img, cpl_vector * a, double threshold);

static int gravi_acqcam_spot_dfda (const double x_in[], const double v[], double result[]);

cpl_error_code gravi_acqcam_fit_spot (cpl_image * img, cpl_size nrand,
                                      cpl_vector * a, const int ia[],
                                      int * nspot);

/*-----------------------------------------------------------------------------
                             Functions code
 -----------------------------------------------------------------------------*/

cpl_error_code gravi_preproc_acqcam (gravi_data *output_data,
                                     gravi_data *input_data,
                                     gravi_data *bad_map)
{
    gravi_msg_function_start(1);
    cpl_ensure_code (output_data, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code (input_data,  CPL_ERROR_NULL_INPUT);

    /* Check if extension exist */
    if (!gravi_data_has_extension (input_data, GRAVI_IMAGING_DATA_ACQ_EXT)) {
        cpl_msg_warning (cpl_func, "Cannot preproc the ACQCAM, not data");
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

/* Compute the 4 diode positions */
inline int gravi_acqcam_sub_diode (const double v[], double *xd, double *yd)
{
    /* Angle */
    double ang = v[8] * CPL_MATH_RAD_DEG;
    // double sang = sin (ang) * 0.5;
    // double cang = cos (ang) * 0.5;
    double sang = sin1 (ang) * 0.5;
    double cang = sin1 (ang + CPL_MATH_PI_2) * 0.5;

    /* Diode arrangement */
    xd[0] = -cang * v[9] + sang * v[10];
    xd[1] = -cang * v[9] - sang * v[10];
    xd[2] =  cang * v[9] - sang * v[10];
    xd[3] =  cang * v[9] + sang * v[10];
    
    yd[0] =  cang * v[10] + sang * v[9];
    yd[1] = -cang * v[10] + sang * v[9];
    yd[2] = -cang * v[10] - sang * v[9];
    yd[3] =  cang * v[10] - sang * v[9];
    
    return 0;
}

/* Compute the 4 sub-aperture positions */
inline int gravi_acqcam_sub_pos (const double v[], double *xsub, double *ysub)
{
    /* Sub-apperture arrangement */
    xsub[0] = v[0] + v[1] + v[2] + v[3];
    xsub[1] = v[0] - v[1] + v[2] - v[3];
    xsub[2] = v[0] + v[1] - v[2] - v[3];
    xsub[3] = v[0] - v[1] - v[2] + v[3];
    
    ysub[0] = v[4] + v[5] + v[6] + v[7];
    ysub[1] = v[4] - v[5] + v[6] - v[7];
    ysub[2] = v[4] + v[5] - v[6] - v[7];
    ysub[3] = v[4] - v[5] - v[6] + v[7];

    return 0;
}

/*----------------------------------------------------------------------------*/

static int gravi_acqcam_spot (const double x_in[], const double v[], double *result)
{
    double fwhm2 = 6.0*6.0;
    *result = 0.0;

    /* Static parameters */
    double xd[4], yd[4], xsub[4], ysub[4];
    gravi_acqcam_sub_diode (v, xd, yd);
    gravi_acqcam_sub_pos (v, xsub, ysub);

    /* Loop on diode and appertures.
     * The capture range is 2.FWHM */
    for (int diode = 0; diode < 4 ; diode++) {
        for (int sub = 0; sub < 4 ; sub++) {
            double xf = (x_in[0] - xsub[sub] - xd[diode]);
            double yf = (x_in[1] - ysub[sub] - yd[diode]);
            double dist = (xf*xf + yf*yf) / fwhm2;
            if (dist < 4.0) *result += exp1 (-dist);
        }
    }

    return 0;
}

/*----------------------------------------------------------------------------*/

static int gravi_acqcam_spot_dfda (const double x_in[], const double v[], double result[])
{
    cpl_size na = 11;
    double next = 0.0, here = 0.0, epsilon = 1e-8;

    double vlocal[na];
    memcpy (vlocal, v, sizeof(double)*na);

    /* Loop on parameters to compute finite differences */
    for (int a = 0; a < na; a++) {

        vlocal[a] += epsilon;
        gravi_acqcam_spot (x_in, vlocal, &next);
        
        vlocal[a] -= 2.*epsilon;
        gravi_acqcam_spot (x_in, vlocal, &here);
        
        result[a] = (next - here) / (2.*epsilon);
        vlocal[a] += epsilon;
    }

    return 0;
}

/*----------------------------------------------------------------------------*/

int gravi_acqcam_spot_count (cpl_image * img, cpl_vector * a, double threshold)
{
    gravi_msg_function_start(0);
    cpl_ensure (img, CPL_ERROR_NULL_INPUT, -1);
    cpl_ensure (a,   CPL_ERROR_NULL_INPUT, -1);
    
    /* Static parameters */
    double xd[4], yd[4], xsub[4], ysub[4];
    double * v = cpl_vector_get_data (a);
    gravi_acqcam_sub_diode (v, xd, yd);
    gravi_acqcam_sub_pos (v, xsub, ysub);
    
    /* Loop on diode and appertures */
    int nspot = 0, nv = 0;
    for (int diode = 0; diode < 4 ; diode++) {
        for (int sub = 0; sub < 4 ; sub++) {
            cpl_size xf = roundl(xsub[sub] + xd[diode]) + 1;
            cpl_size yf = roundl(ysub[sub] + yd[diode]) + 1;
            if (cpl_image_get (img, xf, yf, &nv) > threshold) nspot ++;
        }
    }
    
    gravi_msg_function_exit(0);
    return nspot;
}

cpl_error_code gravi_acqcam_spot_imprint (cpl_image * img, cpl_vector * a)
{
    gravi_msg_function_start(0);
    cpl_ensure_code (img, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code (a,   CPL_ERROR_NULL_INPUT);

    /* Static parameters */
    double xd[4], yd[4], xsub[4], ysub[4];
    double * v = cpl_vector_get_data (a);
    gravi_acqcam_sub_diode (v, xd, yd);
    gravi_acqcam_sub_pos (v, xsub, ysub);
    
    /* Loop on diode and appertures */
    for (int diode = 0; diode < 4 ; diode++) {
        for (int sub = 0; sub < 4 ; sub++) {
            cpl_size xf = roundl(xsub[sub] + xd[diode]) + 1;
            cpl_size yf = roundl(ysub[sub] + yd[diode]) + 1;
            cpl_image_set (img, xf,   yf, 0);
            cpl_image_set (img, xf-1, yf, 0);
            cpl_image_set (img, xf+1, yf, 0);
            cpl_image_set (img, xf, yf+1, 0);
            cpl_image_set (img, xf, yf-1, 0);
        }
    }

    gravi_msg_function_exit(0);
    return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Read the reference position for the pupil sensor, and re-arrange them
 *
 * @param header:   input header
 * @param tel:      requested beam (0..3)
 * @param pupref:   output vector, shall be already allocated. Output are 
 *                  set in the first 8 elements: 4 sub-appertures x {X,Y}
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_acqcam_get_pup_ref (cpl_propertylist * header, cpl_size tel,
                                         cpl_vector * pupref)
{
    gravi_msg_function_start(0);
    cpl_ensure_code (header,          CPL_ERROR_NULL_INPUT);
    cpl_ensure_code (tel>=0 && tel<4, CPL_ERROR_ILLEGAL_INPUT);
    cpl_ensure_code (pupref,          CPL_ERROR_ILLEGAL_INPUT);

    char name[90];
    cpl_size ntel = 4;
    cpl_size nsub = 4;

    cpl_vector * xsub = cpl_vector_new (nsub);
    cpl_vector * ysub = cpl_vector_new (nsub);

    /* Read sub-windows size */
    cpl_size nsx = cpl_propertylist_get_int (header, "ESO DET1 FRAMES NX");
    cpl_size nsy = cpl_propertylist_get_int (header, "ESO DET1 FRAMES NY");
        
    /* Read sub-windows start for pupil */
    sprintf (name, "ESO DET1 FRAM%lld STRX", 3*ntel + tel + 1);
    cpl_size sx = cpl_propertylist_get_int (header, name);
    
    sprintf (name, "ESO DET1 FRAM%lld STRY", 3*ntel + tel + 1);
    cpl_size sy = cpl_propertylist_get_int (header, name);
    
    cpl_msg_debug (cpl_func,"pupil %lli sx= %lld sy = %lld", tel, sx, sy);
        
    /* Read the sub-apperture reference positions 
     * Converted to accound for sub-windowing 
     * In vector convention (start at 0,0) */
    for (int sub = 0; sub < nsub ; sub++) {
        double xv = gravi_pfits_get_ptfc_acqcam (header, sub*ntel + tel + 1);
        double yv = gravi_pfits_get_ptfc_acqcam (header, sub*ntel + tel + 17);
        cpl_vector_set (xsub, sub, xv - (sx - tel*nsx) - 1);
        cpl_vector_set (ysub, sub, yv - (sy - 3*nsy) - 1);
        cpl_msg_debug (cpl_func,"pupil %lli subC %i = %10.4f,%10.4f", tel, sub,
                       cpl_vector_get (xsub, sub),cpl_vector_get (ysub, sub));
        CPLCHECK_MSG ("Cannot get pupil reference position");
    }

    /* All linear combination of sub-appertures - X */
    cpl_vector_set (pupref, 0, 0.25 * 
                    (cpl_vector_get (xsub,0) + cpl_vector_get (xsub,1) +
                     cpl_vector_get (xsub,2) + cpl_vector_get (xsub,3)));
    cpl_vector_set (pupref, 1, 0.25 * 
                    (cpl_vector_get (xsub,0) - cpl_vector_get (xsub,1) +
                     cpl_vector_get (xsub,2) - cpl_vector_get (xsub,3)));
    cpl_vector_set (pupref, 2, 0.25 * 
                    (cpl_vector_get (xsub,0) - cpl_vector_get (xsub,2) +
                     cpl_vector_get (xsub,1) - cpl_vector_get (xsub,3)));
    cpl_vector_set (pupref, 3, 0.25 * 
                    (cpl_vector_get (xsub,0) - cpl_vector_get (xsub,1) +
                     cpl_vector_get (xsub,3) - cpl_vector_get (xsub,2)));
    
    /* All linear combination of sub-appertures - Y */
    cpl_vector_set (pupref, 4, 0.25 * 
                    (cpl_vector_get (ysub,0) + cpl_vector_get (ysub,1) +
                     cpl_vector_get (ysub,2) + cpl_vector_get (ysub,3)));
    cpl_vector_set (pupref, 5, 0.25 * 
                    (cpl_vector_get (ysub,0) - cpl_vector_get (ysub,1) +
                     cpl_vector_get (ysub,2) - cpl_vector_get (ysub,3)));
    cpl_vector_set (pupref, 6, 0.25 * 
                    (cpl_vector_get (ysub,0) - cpl_vector_get (ysub,2) +
                     cpl_vector_get (ysub,1) - cpl_vector_get (ysub,3)));
    cpl_vector_set (pupref, 7, 0.25 * 
                    (cpl_vector_get (ysub,0) - cpl_vector_get (ysub,1) +
                     cpl_vector_get (ysub,3) - cpl_vector_get (ysub,2)));
    
    /* Delete */
    FREE (cpl_vector_delete, xsub);
    FREE (cpl_vector_delete, ysub);

    gravi_msg_function_exit(0);
    return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Fit a pupil spot pattern into an image
 *
 * @param img:    input image
 * @param nrand:  number of random starting point (around specified parameters)
 * @param a:      vector of parameter, modified in-place, in the form:
 *                (x0+x1+x2+x3, x0-x1+x2-x3, x0+x1-x2-x3, x0-x1-x2+x3, 
 *                 y0+y1+y2+y3, y0-y1+y2-y3, y0+y1-y2-y3, y0-y1-y2+y3, 
 *                 diode_angle, delta_diode_x, delta_diode_y)
 * @param ia:     specifies which parameter is to be fitted, and which is to 
 *                be kept unmodified.
 * @param nspot:  is filled with the number of detected spots.
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_acqcam_fit_spot (cpl_image * img,
                                      cpl_size nrand,
                                      cpl_vector * a,
                                      const int ia[],
                                      int * nspot)
{
    gravi_msg_function_start(0);
    cpl_ensure_code (img, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code (a,   CPL_ERROR_NULL_INPUT);
    cpl_ensure_code (nrand>0, CPL_ERROR_ILLEGAL_INPUT);
    cpl_ensure_code (nspot, CPL_ERROR_NULL_INPUT);

    *nspot = 0;
    int nv = 0;
    cpl_size nx = cpl_image_get_size_x (img);
    cpl_size ny = cpl_image_get_size_y (img);

    cpl_size xsub = cpl_vector_get (a, 0);
    cpl_size ysub = cpl_vector_get (a, 4);
    CPLCHECK_MSG ("Cannot get values valid");

    /* Get only pixels around the center,
     * to estimate the noise */
    cpl_size ns = 50;
    cpl_vector * flux = cpl_vector_new (ns*ns);
    for (cpl_size v = 0, x = xsub-ns/2; x < xsub+ns/2; x++) {
        for (cpl_size y = ysub-ns/2; y < ysub+ns/2; y++) {
            cpl_vector_set (flux, v, cpl_image_get (img, x+1, y+1, &nv));
            v++;
            CPLCHECK_MSG ("Cannot fill vector");
        }
    }

    /* Compute typical error as the
     * median of spatial variation */
    cpl_vector_multiply (flux, flux);
    double RMS = sqrt (cpl_vector_get_median (flux));
    FREE (cpl_vector_delete, flux);

    double threshold = 10 * RMS;

    /* Get only valid pixels for this beam */
    cpl_size nw = 200;
    cpl_size nvalid = 0;
    cpl_vector * is_valid = cpl_vector_new (nx * ny);
    for (cpl_size x = xsub-nw/2; x < xsub+nw/2; x++) {
        for (cpl_size y = ysub-nw/2; y < ysub+nw/2; y++) {
            if (cpl_image_get (img, x+1, y+1, &nv) > threshold) {
                cpl_vector_set (is_valid, x*ny + y, 1.0);
                nvalid ++;
            } else {
                cpl_vector_set (is_valid, x*ny + y, 0.0);
                cpl_image_set (img, x+1, y+1, 0);
            }
                
            CPLCHECK_MSG ("Cannot compute valid");
        }
    }
    cpl_msg_debug (cpl_func, "nvalid = %lld with >%.2f (RMS = %.2f)",
                   nvalid, threshold, RMS);

    /* Verify enough pixel, otherwise return */
    if (nvalid < 100) {
        *nspot = 0;
        FREE (cpl_vector_delete, is_valid);
        return CPL_ERROR_NONE;
    }

    /* Fill matrix */
    cpl_matrix * x_matrix  = cpl_matrix_new (nvalid, 2);
    cpl_vector * y_vector  = cpl_vector_new (nvalid);
    cpl_vector * sy_vector = cpl_vector_new (nvalid);
    for (cpl_size v = 0, x = xsub-nw/2; x < xsub+nw/2; x++) {
        for (cpl_size y = ysub-nw/2; y < ysub+nw/2; y++) {
            if (cpl_vector_get (is_valid, x*ny + y)) {
                cpl_matrix_set (x_matrix, v, 0, x);
                cpl_matrix_set (x_matrix, v, 1, y);
                cpl_vector_set (y_vector, v, cpl_image_get (img, x+1, y+1, &nv));
                cpl_vector_set (sy_vector, v, 1.0);
                v ++;
                CPLCHECK_MSG ("Cannot fill matrix/vector");
            }
        }
    }
    FREE (cpl_vector_delete, is_valid);

    /* Normalize y so that the minimum of y is 1.0, which is the maximum
     * of the model function. So that the fit ressemble a correlation.
     * The sqrt is to avoid putting too much weight on bright spots */
    cpl_vector_divide_scalar (y_vector, cpl_vector_get_min (y_vector));
    cpl_vector_sqrt (y_vector);

    /* Output for global minimisation */
    cpl_vector * a_start = cpl_vector_duplicate (a);
    cpl_vector * a_tmp = cpl_vector_duplicate (a);
    double chisq_final = 1e10;
    srand(1);

    /* Loop on various starting points */
    for (cpl_size try = 0; try < nrand; try++) {

        /* Move starting point in position (+-10pix) and angle (entire circle) */
        cpl_vector_copy (a_tmp, a_start);
        if (try > 0) {
            if (ia[0]) cpl_vector_set (a_tmp, 0, cpl_vector_get (a_tmp, 0) + (rand()%20) - 10);
            if (ia[4]) cpl_vector_set (a_tmp, 4, cpl_vector_get (a_tmp, 4) + (rand()%20) - 10);
            if (ia[8]) cpl_vector_set (a_tmp, 8, cpl_vector_get (a_tmp, 8) + (rand()%180));
        }
        
        /* Fit from this starting point */
        double chisq;
        cpl_fit_lvmq (x_matrix, NULL, y_vector, sy_vector,
                      a_tmp, ia,
                      gravi_acqcam_spot,
                      gravi_acqcam_spot_dfda,
                      CPL_FIT_LVMQ_TOLERANCE,
                      CPL_FIT_LVMQ_COUNT,
                      CPL_FIT_LVMQ_MAXITER,
                      NULL, &chisq, NULL);
        CPLCHECK_MSG ("Cannot fit");

        /* Check chi2 and keep if lowest so far */
        if (chisq < chisq_final) {
            cpl_msg_debug (cpl_func, "chisq_final = %f -> %f", chisq_final, chisq);
            chisq_final = chisq;
            cpl_vector_copy (a, a_tmp);
        }

    } /* End loop on try starting points */

    /* Get the image value at the position of spots, and
     * count the number of spots */
    *nspot = gravi_acqcam_spot_count (img, a, threshold);

    FREE (cpl_vector_delete, a_tmp);
    FREE (cpl_vector_delete, a_start);
    FREE (cpl_matrix_delete, x_matrix);
    FREE (cpl_vector_delete, y_vector);
    FREE (cpl_vector_delete, sy_vector);

    gravi_msg_function_exit(0);
    return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Reduce the ACQ camera images
 *  
 * @param data        The gravi_data input/output
 * 
 * The routine read the images from the ACQ camera, reduce these images
 * with the XXXXX algorithm. The resulting tables is append into
 * data.
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
        cpl_msg_warning (cpl_func, "Cannot reduce the ACQCAM, not data");
        return CPL_ERROR_NONE;
    }

    /* Get the header */
    cpl_propertylist * header;
    header = gravi_data_get_header (input_data);
    CPLCHECK_MSG ("Cannot get data or header");

    /* Get the data and header */
    cpl_imagelist * acqcam_imglist;
    acqcam_imglist = gravi_data_get_cube (input_data, GRAVI_IMAGING_DATA_ACQ_EXT);

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
    cpl_table_set_column_unit (acqcam_table, "PUPIL_X", "pix");
    cpl_table_set_column_unit (acqcam_table, "PUPIL_Y", "pix");
    cpl_table_set_column_unit (acqcam_table, "PUPIL_Z", "pix");
    cpl_table_set_column_unit (acqcam_table, "PUPIL_R", "deg");
    
    /* Compute mean image */
    cpl_image * mean_img = cpl_imagelist_collapse_create (acqcam_imglist);
    int nspot = 0;

    /* Loop on tel */
    for (int tel = 0; tel < ntel; tel++) {
        cpl_msg_info (cpl_func, "Compute pupil position for beam %i", tel+1);

        /* Allocate memory */
        cpl_vector * a_start = cpl_vector_new (11);
        
        /* Read the sub-apperture reference positions 
         * Converted to accound for sub-windowing 
         * In vector convention (start at 0,0) */
        gravi_acqcam_get_pup_ref (header, tel, a_start);
        CPLCHECK_MSG ("Cannot read ACQ data for pupil in header");

        /* Diode rotation [deg], and spacing in x and y [pix] */
        cpl_vector_set (a_start, 8,  0.0); // deg
        cpl_vector_set (a_start, 9, 18.0); // pix
        cpl_vector_set (a_start,10, 24.0); // pix
        
        cpl_vector * a_final = cpl_vector_duplicate (a_start);
        CPLCHECK_MSG ("Cannot prepare parameters");
        
        /* First fit: global center and rotation only */
        const int ia_global[] = {1,0,0,0, 1,0,0,0, 1,0,0};
        gravi_acqcam_fit_spot (mean_img, 30, a_final, ia_global, &nspot);
        CPLCHECK_MSG ("Cannot fit rotation and center");

        /* Second fit: independ sub-appertures
         * and free diode spacing */
        const int ia_all[] = {1,1,1,1, 1,1,1,1, 1,1,1};
        gravi_acqcam_fit_spot (mean_img, 1, a_final, ia_all, &nspot);
        CPLCHECK_MSG ("Cannot fit sub-appertures");

        cpl_msg_info (cpl_func, "Found %i spots in mean img of tel %i", nspot, tel+1);
        
        /* Add best position as a cross in image */
        gravi_acqcam_spot_imprint (mean_img, a_final);

        /* Loop on all images */
        for (cpl_size row = 0; row < nrow; row++) {
            if (row %10 == 0)
                cpl_msg_info_overwritable (cpl_func, "Fit image %lld over %lld", row+1, nrow);
            
            /* Get data and first guess */
            cpl_image * img = cpl_imagelist_get (acqcam_imglist, row);
            cpl_vector * a_row = cpl_vector_duplicate (a_final);

            /* Fit all */
            const int ia_row[] = {1,1,1,1, 1,1,1,1, 1,0,0};
            gravi_acqcam_fit_spot (img, 1, a_row, ia_row, &nspot);
            CPLCHECK_MSG ("Cannot fit sub-appertures of image");

            /* Add best position as a cross in image */
            gravi_acqcam_spot_imprint (img, a_row);

            /* Remove reference */
            cpl_vector_subtract (a_row, a_start);

            /* Compute latteral shift */
            double x_shift = cpl_vector_get (a_row, 0);
            double y_shift = cpl_vector_get (a_row, 4);
            double r_shift = cpl_vector_get (a_row, 8);
            
            /* Compute longitudinal shift */
            double z_shift = -0.5 * ( cpl_vector_get (a_row, 2) +
                                      cpl_vector_get (a_row, 5));

            /* Fill table */
            cpl_table_set (acqcam_table, "PUPIL_NSPOT", row*ntel+tel, nspot);
            if (nspot > 8) {
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
 * @brief Return 1 if flux(pos) < flux(pos+1)
 * 
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

