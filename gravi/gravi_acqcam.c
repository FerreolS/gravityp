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

/*
 * History :
 *    26/11/2018 Read parameter roof_x, roof_y, sot_x, spot_y and roof_pos from calibration file
 *    28/11/2018 Read parameter plate scale from calibration file
 *    03/12/2018 Read gravi_acqcam_pupil parameter from calibration file
 *    10/01/2019 Fix compiler Warnings : inititialize with NULL : roof_x, roof_y, sot_x, spot_y and roof_pos
 *                                           unused parameter nsy
 *    12/09/2019 EKW Flip parameter from cpl_vector_set (output, GRAVI_SPOT_DIODE+0, 0.122);
 *                                       cpl_vector_set (output, GRAVI_SPOT_DIODE+1, 0.158)
 *                                  to   cpl_vector_set (output, GRAVI_SPOT_DIODE+0, 0.158);
 *                                       cpl_vector_set (output, GRAVI_SPOT_DIODE+1, 0.122)
 *    18/09/2019 EkW Double spot fitting as Ollis analysis showed problems
 *    16/10/2019 EkW Remove debug statements
 *    06/11/2019 EKW Add Franks double pscale = sep ? rho_in/sep : scale (instead 0) ;
 *    06/11/2019 EKW Commit with right PIPE-8675 as despite the line abovem, all changes are related to it
 *    15/07/2020 EKW JIRA PIPE-9123 : Pipeline faile with Strehl=NaN
 */
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
                               Defines & Macros
 -----------------------------------------------------------------------------*/

/* Number of parameters in the model 'gravi_acqcam_spot'
 * And position of parameters */
#define GRAVI_SPOT_NTEL 4
#define GRAVI_SPOT_NLENS 4
#define GRAVI_SPOT_NSPOT 4
#define GRAVI_SPOT_NFOCUS 21
#define GRAVI_SPOT_NSEARCH 45
#define GRAVI_SPOT_SWINDOW 28

#define GRAVI_SPOT_NA    30
#define GRAVI_SPOT_SUB   0
#define GRAVI_SPOT_ANGLE 8
#define GRAVI_SPOT_SCALE 9
#define GRAVI_SPOT_DIODE 10
#define GRAVI_SPOT_FWHM  13
#define GRAVI_SPOT_FLUX  14

#define GRAVI_ACQ_PUP_FLUX 1e6

#define CHECK_BIT(var,pos) ((var) & (1<<(pos)))

/*-----------------------------------------------------------------------------
                               Private prototypes
 -----------------------------------------------------------------------------*/

double exp1 (double x);
double sin1 (double x);
int gravi_acqcam_xy_diode (const double v[], double *xd, double *yd);

static int gravi_acqcam_spot (const double x_in[], const double v[], double *result);
static int gravi_acqcam_xy_sub (const double v[], double *xsub, double *ysub);

double gravi_acqcam_z2meter (double PositionPixels, gravi_data *static_param_data);

cpl_error_code gravi_acqcam_pupil (cpl_image * mean_img,
                                   cpl_imagelist * acqcam_imglist,
                                   cpl_propertylist * header,
                                   cpl_table * acqcam_table,
                                   cpl_propertylist * o_header,
								   gravi_data *static_param_data);

cpl_error_code gravi_acqcam_field (cpl_image * mean_img,
                                   cpl_imagelist * acqcam_imglist,
                                   cpl_propertylist * header,
                                   cpl_table * acqcam_table,
                                   cpl_propertylist * o_header,
                                   gravi_data *static_param_data);

cpl_error_code gravi_acq_fit_gaussian (cpl_image * img, double *x, double *y,
                                       double *ex, double *ey, cpl_size size);

cpl_error_code gravi_acq_measure_strehl(cpl_image * img, double x, double y, 
                                        double pscale, double *SR, cpl_propertylist * header);

cpl_error_code gravi_acq_measure_max(cpl_image * img, double x, double y, double size, double * img_max);

cpl_error_code gravi_image_fft_correlate (cpl_image *ia, cpl_image *ib, cpl_size *xd, cpl_size *yd);

cpl_error_code gravi_acqcam_pupil_v2 (cpl_image * mean_img,
                                      cpl_imagelist * acqcam_imglist,
                                   cpl_propertylist * header,
                                   cpl_table * acqcam_table,
                                   cpl_propertylist * o_header,
                                   gravi_data *static_param_data);

cpl_error_code gravi_acqcam_clean_pupil_v2(cpl_imagelist * acqcam_imglist, cpl_imagelist * pupilImage_filtered, const cpl_size ury);

cpl_error_code gravi_acqcam_select_good_frames_v2(cpl_imagelist * acqcam_imglist, cpl_imagelist * pupilImage_onFrames, cpl_array * good_frames);

cpl_error_code gravi_acqcam_get_pup_ref_v2 (cpl_propertylist * header,
                                            cpl_bivector *  diode_pos_subwindow);

cpl_error_code gravi_acqcam_get_diode_ref_v2 (cpl_propertylist * header,
                                              cpl_array * good_frames,
                                      cpl_vector * scale_vector,
                                      cpl_bivector **  diode_pos_telescope,
                                              int nrow_on);

cpl_error_code gravi_acqcam_get_diode_theoretical_v2(cpl_bivector *  diode_pos_subwindow,
                                                     cpl_bivector **  diode_pos_telescope,
                                                     cpl_bivector **  diode_pos_theoretical,
                                                     cpl_size nrow_on, int ury);

cpl_error_code gravi_acqcam_spot_imprint_v2(cpl_image *mean_img, cpl_bivector ** diode_pos_offset, cpl_bivector **  diode_pos_theoretical, int ury);


cpl_error_code gravi_acqcam_perform_shiftandadd_v2(cpl_imagelist * pupilImage_onFrames, cpl_imagelist ** pupilImage_shiftandadd, cpl_array * good_frames,
                                                   cpl_vector * focus_value,
                                                     cpl_bivector **  diode_pos_theoretical,
                                                     cpl_bivector **  diode_pos_offset, cpl_size nrow_on);

cpl_error_code gravi_acqcam_get_pupil_offset_v2(cpl_imagelist ** pupilImage_shiftandadd,
                                                     cpl_array * good_frames,
                                                cpl_bivector **  diode_pos_offset, cpl_propertylist * o_header,
                                                     cpl_size nrow_on);

cpl_error_code gravi_acqcam_set_pupil_table_v2(cpl_table * acqcam_table, cpl_propertylist * header, cpl_propertylist * o_header,  cpl_array * good_frames,
                                               cpl_vector* scale_vector, cpl_array * bad_frames_short, cpl_bivector **  diode_pos_offset , cpl_vector * focus_value, gravi_data *static_param_data);

cpl_image * gravi_image_extract(cpl_image * image_in, cpl_size llx, cpl_size lly, cpl_size urx, cpl_size ury);

double gravi_acqcam_defocus_scaling(int focus);

/* This global variable optimises the computation
 * of partial derivative on fitted parameters */
const int * GRAVI_LVMQ_FREE;
const int * GRAVI_LVMQ_FREE = NULL;


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
        gravi_msg_warning (cpl_func,"Cannot preproc the ACQCAM, no data in file");
        return CPL_ERROR_NONE;
    }
    if (!gravi_data_has_extension (bad_map, GRAVI_IMAGING_DATA_ACQ_EXT)) {
        gravi_msg_warning (cpl_func,"Cannot preproc the ACQCAM, no badpixel in BAD");
        return CPL_ERROR_NONE;
    }

    /* get image of badpixels */
    cpl_image * badpix_img = gravi_data_get_img (bad_map, GRAVI_IMAGING_DATA_ACQ_EXT);
    CPLCHECK_MSG ("Cannot get BAD map for ACQ");

    /* Get the imagelist */
    cpl_imagelist * imglist;
    imglist = gravi_data_get_cube (input_data, GRAVI_IMAGING_DATA_ACQ_EXT);
    CPLCHECK_MSG ("Cannot get image for ACQ");
    
    /* check for image size */
    cpl_image *  image0 =   cpl_imagelist_get (imglist, 0);
    if ( (cpl_image_get_size_x (image0) != cpl_image_get_size_x (badpix_img) ) ||
         (cpl_image_get_size_y (image0) != cpl_image_get_size_y (badpix_img) ) ) {
        gravi_msg_warning (cpl_func,"Cannot preproc the ACQCAM. Bad pixel mask does not have correct size.");
        return CPL_ERROR_NONE;
    }
    
    /* Construct a mask of badpixels */
    cpl_mask * badpix_mask = cpl_mask_threshold_image_create (badpix_img, 0.5, 10000);
    
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
#if 0
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
#endif

/*----------------------------------------------------------------------------*/
/**
 * @brief measure Strehl Ratio of the source at the given location
 *
 * @param img      input image
 * @param x,y      position in image in FITS convention (1,1 is lower,left)
 * @param pscale   plate scale
 * @param SR       returned strehl value
 * @param header   property list to read the telescope (UTs or ATs)
 *
 * The function extrat a sub-image of +-50pixel around the x,y position
 * and run the hdrl_strehl_compute function.
 **/
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_acq_measure_strehl(cpl_image * img, double x, double y, 
                                        double pscale, double *SR,
                                        cpl_propertylist * header)
{
    gravi_msg_function_start(0);
    cpl_ensure_code (img,        CPL_ERROR_NULL_INPUT);
    cpl_ensure_code (SR,         CPL_ERROR_NULL_INPUT);
    cpl_ensure_code (header,     CPL_ERROR_NULL_INPUT);
    cpl_ensure_code (pscale > 0, CPL_ERROR_ILLEGAL_INPUT);

    cpl_size hw = 50;
    cpl_ensure_code (x > hw, CPL_ERROR_ILLEGAL_INPUT);
    cpl_ensure_code (y > hw, CPL_ERROR_ILLEGAL_INPUT);

    /* Structure to define the parameters
       for fitting strehl in image */
    hdrl_parameter * strehl_params;

    /* Hardcoded Mirror diameters */
    const char * telname = gravi_conf_get_telname (0, header);
    CPLCHECK ("Cannot get telescope name");
    
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

    /* Extract sub-image */
    cpl_image * sub_image = cpl_image_extract (img, x-hw, y-hw, x+hw, y+hw);
    hdrl_image * sub_hdrl = hdrl_image_create (sub_image, NULL);

    /* Run Strehl algorithm from HDRL package */
    hdrl_strehl_result strehl;
    strehl = hdrl_strehl_compute (sub_hdrl, strehl_params);
    *SR = (double) strehl.strehl_value.data;

    /* Delete */
    FREE (cpl_image_delete, sub_image);
    FREE (hdrl_image_delete, sub_hdrl);
    FREE (hdrl_parameter_delete, strehl_params);

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
    
    cpl_size x_med_min=(cpl_size)(*x)-size;
    cpl_size x_med_max=(cpl_size)(*x)+size;
    cpl_size y_med_min=(cpl_size)(*y)-size;
    cpl_size y_med_max=(cpl_size)(*y)+size;
    cpl_size  nx = cpl_image_get_size_x (img);
    cpl_size  ny = cpl_image_get_size_y (img);
    if (x_med_min <1) x_med_min=1;
    if (y_med_min <1) y_med_min=1;
    if (x_med_max > nx) x_med_max=nx;
    if (y_med_max > ny) y_med_max=ny;
    
    double med = cpl_image_get_median_window (img,
                                              x_med_min, y_med_min,
                                              x_med_max, y_med_max);
    cpl_array_set (parameters, 0, med);
    CPLCHECK_MSG ("Error getting median");
    
    /* Fit Gaussian */
    double rms=0.;
    
    cpl_error_code fit_converged =  cpl_fit_image_gaussian (img, NULL, (cpl_size)(*x), (cpl_size)(*y),
                            size, size, parameters,
                            NULL, NULL, &rms, NULL, NULL,
                            NULL, NULL, NULL, NULL);
    
    if (fit_converged == CPL_ERROR_NONE)
    {
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
            
        // cf. Condon 1996, PASP 109:166
        double cst = 2. * sqrt(2. * M_PI * (1. - rho*rho) * sx * sy) * rms / A;
        *ex=cst*sx;
        *ey=cst*sy;
        
        if ( A < 0. ) {
          // detection is just not significant
          cpl_msg_info (cpl_func, "rejecting fit: x=%g, y=%g, SNR=%g, sx=%g, sy=%g",
                        *x, *y, A/(rms * 2.*M_PI*sx*sy*sqrt(1-rho*rho)), sx, sy);
          *x = 0.;
          *y = 0.;
          *ex = -1.;
          *ey = -1.;
        }
        CPLCHECK_MSG("Error checking significance of fit result");
    } else {
        cpl_msg_warning (cpl_func, "fit of pupil beacon did not converge");
        *x = 0.;
        *y = 0.;
        *ex = -1.;
        *ey = -1.;
        cpl_error_reset();
    }
    CPLCHECK_MSG ("Pupil beacon fit failed");

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
 * @brief Reduce the images of pupil from ACQ V2.0
 *
 * @param mean_img:  input image
 * @param acqcam_imglist:    input image list
 * @param header:            input header
 * @param acqcam_table:      output table
 * @param o_header:          output header
 * @param static_param_data:          gravity static dataset
 *
 * \exception CPL_ERROR_NULL_INPUT input data is missing
 *
 * The routine analyse the pupil from ACQ and create QC parameters in the
 * header, as well as columns in the acqcam_table.
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_acqcam_pupil_v2(cpl_image * mean_img, cpl_imagelist * acqcam_imglist, cpl_propertylist * header,
        cpl_table * acqcam_table, cpl_propertylist * o_header,
        gravi_data *static_param_data)
{
gravi_msg_function_start(1);
cpl_ensure_code(acqcam_imglist, CPL_ERROR_NULL_INPUT);
cpl_ensure_code(header, CPL_ERROR_NULL_INPUT);
cpl_ensure_code(acqcam_table, CPL_ERROR_NULL_INPUT);
cpl_ensure_code(o_header, CPL_ERROR_NULL_INPUT);
    
cpl_size  nx = cpl_image_get_size_x (cpl_imagelist_get (acqcam_imglist, 0));
cpl_size  ny = cpl_image_get_size_y (cpl_imagelist_get (acqcam_imglist, 0));

/*
 * constant ury is the beginning of the upper part of the acquisition camera
 * It delimitates where the pupil beacon images starts
 * Not sure it works in full frame case (ny>1100)
 */
const cpl_size ury = (ny>1100) ? 1200 : 745;
    
/*
 * First step. The goal is to get the images pupil, cutted to keep only pupil frames (using ury).
 * The output is the pupil, after kernel filtering with a Gaussian pdf the size of the pupil beacon
 */
    cpl_imagelist * pupilImage_filtered  =  cpl_imagelist_new ();
    gravi_acqcam_clean_pupil_v2(acqcam_imglist, pupilImage_filtered, ury);
    CPLCHECK_MSG("Cannot clean pupil");
    
/*
 * Second step. The goal is to select the frames with the blinking beacons on
 * The same routines gets the frames where the beacons are off
 * It uses theses frames to remove the background
 */
    cpl_array * good_frames= cpl_array_new(cpl_imagelist_get_size(pupilImage_filtered),CPL_TYPE_INT);
    cpl_imagelist * pupilImage_onFrames = cpl_imagelist_new ();
    gravi_acqcam_select_good_frames_v2(pupilImage_filtered,pupilImage_onFrames,good_frames);
    CPLCHECK_MSG("Cannot find blinking pupil frames");
    
/* Third step. The goal is to initialize the vectors which will be used to know where the diodes
 * are on the detectors (hence the use of bivectors, which contains x and y positions
 * the size of the arrays are defined using several fixed parameters:
 * GRAVI_SPOT_NTEL = 4 corresponds to the 4 telescopes
 * GRAVI_SPOT_NSPOT = 4 corresponds to the 4 diodes on each telescopes
 * GRAVI_SPOT_NLENS = 4 corresponds to the 4 lenslets on the acquisition camera
 * nrow_on= 122 (for example) which corresponds to the N frames where the diodes are ON */

    cpl_size nrow    = cpl_imagelist_get_size(pupilImage_filtered);
    cpl_size nrow_on = cpl_imagelist_get_size(pupilImage_onFrames);
    
    if ((nrow_on<1)||(nrow<4))
    {
        /* exiting with empty pupil table */
    cpl_array_fill_window_int (good_frames, 0, nrow_on, 0);
       gravi_acqcam_set_pupil_table_v2(acqcam_table, header, o_header, good_frames, NULL,  NULL, NULL, NULL, NULL);
       CPLCHECK_MSG("Cannot stor pupil offset values in OI_ACQ table");
        
        cpl_imagelist_delete(pupilImage_onFrames);
        cpl_imagelist_delete(pupilImage_filtered);
        cpl_array_delete(good_frames);
        cpl_msg_warning (cpl_func, "Cannot reduce the %lli frames (not enough frames)",nrow);
        gravi_msg_function_exit(1);
        return CPL_ERROR_NONE;
    }
    
    cpl_vector * scale_vector = cpl_vector_new(GRAVI_SPOT_NTEL);
    cpl_vector * focus_value = cpl_vector_new(GRAVI_SPOT_NTEL);
    cpl_array * bad_frames_short= cpl_array_new(nrow_on,CPL_TYPE_INT);
    cpl_array_fill_window_int (bad_frames_short, 0, nrow_on, 0);
    cpl_bivector *  diode_pos_subwindow =  cpl_bivector_new (GRAVI_SPOT_NLENS*GRAVI_SPOT_NTEL);
    cpl_bivector **  diode_pos_telescope = cpl_malloc (nrow_on * sizeof(cpl_bivector *)) ;
    cpl_bivector **  diode_pos_theoretical = cpl_malloc (nrow_on * sizeof(cpl_bivector *)) ;
    cpl_bivector **  diode_pos_offset = cpl_malloc (nrow_on * sizeof(cpl_bivector *)) ;
    cpl_imagelist ** pupilImage_shiftandadd = cpl_malloc (GRAVI_SPOT_NTEL * sizeof(cpl_imagelist *));
    
    for (int n = 0 ; n < nrow_on; n++)
    {
        diode_pos_telescope[n] = cpl_bivector_new (GRAVI_SPOT_NSPOT*GRAVI_SPOT_NTEL);
        diode_pos_theoretical[n] = cpl_bivector_new (GRAVI_SPOT_NFOCUS*GRAVI_SPOT_NSPOT*GRAVI_SPOT_NLENS*GRAVI_SPOT_NTEL);
        diode_pos_offset[n] = cpl_bivector_new (GRAVI_SPOT_NTEL);
    }
    
    for (int tel = 0 ; tel < GRAVI_SPOT_NTEL; tel++)
        pupilImage_shiftandadd[tel] = cpl_imagelist_new();
    CPLCHECK_MSG("Cannot initialized vectors");

/*
 * Fourth step. Here we get the position of the 4 subwindows (correspond to each lenslet), in pixel/detector space
 */
    gravi_acqcam_get_pup_ref_v2(header, diode_pos_subwindow);
    CPLCHECK_MSG("Cannot get pupil reference values");
    
/*
 * Fifth step. Here we get the position of the 4 diodes of each 4 telescopes, in pixel/detector space
 */
    gravi_acqcam_get_diode_ref_v2(header, good_frames, scale_vector, diode_pos_telescope, nrow_on);
    CPLCHECK_MSG("Cannot get telescope beacons position");
    
/*
 * Seventh step. Get a large array of diode vector position (diode_pos_theoretical)
 * this for each telescope, diode, lenslet, and focal number.
 */
    gravi_acqcam_get_diode_theoretical_v2(diode_pos_subwindow ,diode_pos_telescope,
                                          diode_pos_theoretical, nrow_on, ury);
    CPLCHECK_MSG("Cannot get theoretical position of beacons");
    
/*
 * Eigth step. Test each focus position to see which one is best. Then, collapse all spot images according to that focus value.
 * The result pupilImage_shiftandadd correspond to an image for each blinking frame and each telescope
 */
    gravi_acqcam_perform_shiftandadd_v2(pupilImage_onFrames, pupilImage_shiftandadd, good_frames, focus_value, diode_pos_theoretical,
                                                                                diode_pos_offset, nrow_on);
    CPLCHECK_MSG("Cannot find correct focus to shift and add pupil images");
    
/*
 * Nineth step. For each one of the pupilImage_shiftandadd, fit a gaussian to know where the pupil is
 * The measurement is done as an offset with respect to the theroetical position of the beacons
 */
    gravi_acqcam_get_pupil_offset_v2(pupilImage_shiftandadd, bad_frames_short, diode_pos_offset , o_header, nrow_on);
    CPLCHECK_MSG("Cannot get pupil offset values");
    
/*
 * Last step. Storing the value of the pupil position in a table
 */
    gravi_acqcam_set_pupil_table_v2(acqcam_table, header, o_header, good_frames, scale_vector, bad_frames_short, diode_pos_offset, focus_value, static_param_data);
    CPLCHECK_MSG("Cannot stor pupil offset values in OI_ACQ table");
    
 /*
  * Bonus. Add collapsed image to mean image and imprint best position as a cross in image
  */
    
    cpl_image *  pupilImage_onFrames_collapse=cpl_imagelist_collapse_create (pupilImage_onFrames);
    cpl_image *  pupilImage_onFrames_1frame=cpl_imagelist_get (pupilImage_onFrames,0);
    gravi_image_replace_window (mean_img, pupilImage_onFrames_1frame, 1, ury, nx, ny, 1, 1);
    gravi_acqcam_spot_imprint_v2(mean_img, diode_pos_offset, diode_pos_theoretical, ury);
    CPLCHECK_MSG("Cannot modify mean_image");
    
/* cleaning all arrays ( free memory ) */
    
for (int tel = 0 ; tel < GRAVI_SPOT_NTEL; tel++)
{
    if (pupilImage_shiftandadd[tel] != NULL)
        cpl_imagelist_delete( pupilImage_shiftandadd[tel]);
}
FREE (cpl_free, pupilImage_shiftandadd);

    
    for (int n = 0 ; n < nrow_on; n++)
    {
        cpl_bivector_delete( diode_pos_telescope[n]);
        cpl_bivector_delete( diode_pos_theoretical[n]);
        cpl_bivector_delete( diode_pos_offset[n]);
    }
    FREE (cpl_free, diode_pos_telescope);
    FREE (cpl_free, diode_pos_theoretical);
    FREE (cpl_free, diode_pos_offset);

    cpl_bivector_delete( diode_pos_subwindow);
    cpl_image_delete(pupilImage_onFrames_collapse);
    cpl_imagelist_delete(pupilImage_onFrames);
    cpl_imagelist_delete(pupilImage_filtered);
    cpl_array_delete(good_frames);
    cpl_array_delete(bad_frames_short);
    cpl_vector_delete(scale_vector);
    cpl_vector_delete(focus_value);

CPLCHECK_MSG("Failed at freeing memory");
gravi_msg_function_exit(1);
return CPL_ERROR_NONE;
}
    
/*----------------------------------------------------------------------------*/
/**
 * @brief Cleaning pupil images by cross-correlation with gaussian function
 *
 * @param acqcam_imglist:   input acqcam_imglist
 * @param pupilImage_filtered:   output pupilImage_filtered
 * @param ury:      y limit of pupil beacon camera
 *
 * \exception CPL_ERROR_NULL_INPUT input data is missing
 *
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_acqcam_clean_pupil_v2(cpl_imagelist * acqcam_imglist, cpl_imagelist * pupilImage_filtered, const cpl_size ury)
{
    gravi_msg_function_start(1);
    
    cpl_ensure_code(acqcam_imglist, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code(pupilImage_filtered, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code(ury, CPL_ERROR_NULL_INPUT);
    
    /* Number of row */
    cpl_size nrow = cpl_imagelist_get_size(acqcam_imglist);

    cpl_size  nx = cpl_image_get_size_x (cpl_imagelist_get (acqcam_imglist, 0));
    cpl_size  ny = cpl_image_get_size_y (cpl_imagelist_get (acqcam_imglist, 0));
    
    cpl_imagelist * pupilImage  =  cpl_imagelist_new ();

    for (cpl_size n = 0; n < nrow ; n++)
    {
    cpl_image * small_img_tmp = cpl_imagelist_get (acqcam_imglist, n);
    cpl_image * small_img = cpl_image_extract( small_img_tmp, 1, ury, nx, ny);
    cpl_imagelist_set(pupilImage, small_img,n);
    }
    
    nx=cpl_image_get_size_x (cpl_imagelist_get (pupilImage, 0));
    ny=cpl_image_get_size_y (cpl_imagelist_get (pupilImage, 0));
    
    cpl_image*   pupilImageBackground = cpl_imagelist_collapse_minmax_create    (pupilImage,(cpl_size)(nrow/8),(cpl_size)(5*nrow/8));
    
    cpl_imagelist_subtract_image (pupilImage,pupilImageBackground);
    CPLCHECK_MSG("Failure to subtract background");
    
    cpl_size Npref = 21;
    cpl_matrix * kernel1 = cpl_matrix_new (Npref,Npref);
    cpl_matrix * kernel2 = cpl_matrix_new (Npref,Npref);
    for (cpl_size x =0;x< Npref; x++)
    for (cpl_size y =0;y< Npref; y++)
    {
        double radius_square=(x-(Npref-1)/2)*(x-(Npref-1)/2)+(y-(Npref-1)/2)*(y-(Npref-1)/2);
        cpl_matrix_set(kernel1,x,y,exp(- radius_square /30));
        cpl_matrix_set(kernel2,x,y,exp(- radius_square /40));
    }
    
    cpl_matrix_divide_scalar (kernel1, cpl_matrix_get_mean (kernel1));
    cpl_matrix_divide_scalar (kernel2, cpl_matrix_get_mean (kernel2));
    cpl_matrix_subtract    (kernel1,kernel2);
    cpl_matrix_divide_scalar (kernel1, cpl_matrix_get_stdev (kernel1));
    cpl_matrix_divide_scalar (kernel1, cpl_matrix_get_max (kernel1));
    
    CPLCHECK_MSG("Error computing gaussian Kernel");
    
    for (cpl_size n = 0; n < nrow ; n++)
    {
    cpl_image * small_img_tmp = cpl_imagelist_get (pupilImage, n);
    cpl_image * filtered_img = cpl_image_new (nx, ny, CPL_TYPE_DOUBLE);
    cpl_image_filter    (    filtered_img,
                         small_img_tmp, kernel1,CPL_FILTER_LINEAR,CPL_BORDER_FILTER);
    cpl_imagelist_set(pupilImage_filtered, filtered_img,n);
        if (n %10 == 0 || n == (nrow-1))
        cpl_msg_info_overwritable (cpl_func, "Convolution of pupil frame %lli over %lli frames",n, nrow);
        
    }
    
    cpl_matrix_delete(kernel1);
    cpl_matrix_delete(kernel2);
    
    cpl_image_delete(pupilImageBackground);
    cpl_imagelist_delete(pupilImage);
    
    
    
    
    CPLCHECK_MSG("Pupil Fitting V2 does not work");
    gravi_msg_function_exit(1);
    return CPL_ERROR_NONE;
    
}
    

/*----------------------------------------------------------------------------*/
/**
* @brief select pupil frames with pupil beacon on. Clean pupil frames by substraction of images with pupil beacon off.
*
* @param pupilImage_filtered:   input imagelist with cleaned (filtered) images
 * @param pupilImage_onFrames:   output pupilImage with beacons on
 * @param good_frames:   array of integer. It tells which frames have beacon on or off
* @param ury:      y limit of pupil beacon camera
*
* \exception CPL_ERROR_NULL_INPUT input data is missing
*
 * Not that the output imagelist is a list of n_on images, with n_on the number of good_frames at off:
 *   n_on= sum(good_frames)
*/
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_acqcam_select_good_frames_v2(cpl_imagelist * pupilImage_filtered, cpl_imagelist * pupilImage_onFrames, cpl_array * good_frames)
{
    gravi_msg_function_start(1);
    int nv =0;
    
    cpl_ensure_code(pupilImage_filtered, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code(good_frames, CPL_ERROR_NULL_INPUT);
    
    /* Number of row , size of image */
    cpl_size nrow = cpl_imagelist_get_size(pupilImage_filtered);
    cpl_size  nx = cpl_image_get_size_x (cpl_imagelist_get (pupilImage_filtered, 0));
    cpl_size  ny = cpl_image_get_size_y (cpl_imagelist_get (pupilImage_filtered, 0));
    
    cpl_vector * pupil_max = cpl_vector_new(4);
    cpl_vector * pupil_mediam = cpl_vector_new(nrow);
    cpl_array * back_frames= cpl_array_new(nrow,CPL_TYPE_INT);
    
    for (cpl_size n = 0; n < nrow; n++)
    {
        
        cpl_image * image = cpl_imagelist_get (pupilImage_filtered,n);
        cpl_vector_set(pupil_max, 0, cpl_image_get_max_window (image, 1, 1, 250,  ny));
        cpl_vector_set(pupil_max, 1, cpl_image_get_max_window (image, 251, 1, 500,  ny));
        cpl_vector_set(pupil_max, 2, cpl_image_get_max_window (image, 501, 1, 750,  ny));
        cpl_vector_set(pupil_max, 3, cpl_image_get_max_window (image, 751, 1, nx,  ny));
        cpl_vector_set(pupil_mediam, n, cpl_vector_get_median (pupil_max));
        /*
        cpl_msg_info (cpl_func, "tototo image valeur %lli = %.2f",n, cpl_image_get(image, 101, 151, &nv));
        cpl_msg_info (cpl_func, "tototo %.2f,  %.2f,  %.2f,  %.2f", cpl_vector_get(pupil_max, 0),cpl_vector_get(pupil_max, 1),cpl_vector_get(pupil_max, 2),cpl_vector_get(pupil_max, 3));
        cpl_msg_info (cpl_func, "Finding min brightness level for beacon ON %lli: %.5f ADU",n, cpl_vector_get(pupil_mediam, n));*/
        
    }
    
    cpl_vector *   pupil_mediam_sort = cpl_vector_duplicate (pupil_mediam);
    cpl_vector_sort (pupil_mediam_sort,CPL_SORT_ASCENDING);
    cpl_vector * pupil_median_down = cpl_vector_extract (pupil_mediam_sort, 0, (nrow-1)/2, 1);
    cpl_vector * pupil_median_up = cpl_vector_extract (pupil_mediam_sort, nrow/2, nrow-1, 1);
    double threshold_up = cpl_vector_get_median(pupil_median_up);
    double threshold_down = cpl_vector_get_median(pupil_median_down);
    double threshold=(threshold_up+threshold_down)/2.0;
    
    CPLCHECK_MSG("Could not get median maximum of beacon flux ");
    cpl_msg_info (cpl_func, "Found threshold brightness level for beacon : %.5f ADU", threshold);
    
    for (cpl_size n = 0; n < nrow; n++)
    {
        if (cpl_vector_get(pupil_mediam, n) >threshold *1.1)
            cpl_array_set_int (good_frames,n,1);
        else
            cpl_array_set_int (good_frames,n,0);
        if (cpl_vector_get(pupil_mediam, n) <threshold * 0.9)
            cpl_array_set_int (back_frames,n,1);
        else
            cpl_array_set_int (back_frames,n,0);
    }
    
    cpl_size frames_background = 4;
    cpl_size n_goodFrames = 0;
    
    CPLCHECK_MSG("Failed to find blincking pupil files");
    
    for (cpl_size n = 0; n < nrow; n++)
        if (cpl_array_get_int(good_frames,n,&nv) == 1)
        {
            cpl_image * image = cpl_imagelist_get (pupilImage_filtered,n);
            /*cpl_msg_info (cpl_func, "Good frame image valeur %lli = %.2f",n, cpl_image_get(image, 101, 151, &nv));*/
            
            cpl_size n_frames_background = 0;
            cpl_imagelist * pupilImage_offFrames  =  cpl_imagelist_new ();
                for (cpl_size b = 1; b < nrow; b++)
                {
                    if (n+b < nrow)
                    if ((cpl_array_get_int(back_frames,n+b,&nv) == 1)&(n_frames_background<frames_background))
                    {
                        cpl_image * image_background = cpl_imagelist_get (pupilImage_filtered,n+b);
                        cpl_image * image_background_copy  =  cpl_image_duplicate (image_background);
                        cpl_imagelist_set(pupilImage_offFrames, image_background_copy, n_frames_background);
                        n_frames_background = n_frames_background + 1;
                    }
                    if (n-b >= 0)
                    if ((cpl_array_get_int(back_frames,n-b,&nv) == 1)&(n_frames_background<frames_background))
                    {
                        cpl_image * image_background = cpl_imagelist_get (pupilImage_filtered,n-b);
                        cpl_image * image_background_copy  =  cpl_image_duplicate (image_background);
                        cpl_imagelist_set(pupilImage_offFrames, image_background_copy, n_frames_background);
                        n_frames_background = n_frames_background + 1;
                    }
                }
            
            cpl_image * image_background_mean = cpl_imagelist_collapse_create (pupilImage_offFrames);
            cpl_image * image_background_subtracted = cpl_image_subtract_create (image, image_background_mean);
            cpl_imagelist_set(pupilImage_onFrames, image_background_subtracted, n_goodFrames);
            
            cpl_image_delete(image_background_mean);
            cpl_imagelist_delete(pupilImage_offFrames);
             
            n_goodFrames ++;
        }
    cpl_msg_info(cpl_func, "Found %lli frames with pupil beacon ON (over %lli frames)",n_goodFrames, nrow);
    
    cpl_vector_delete(pupil_max);
    cpl_vector_delete(pupil_mediam);
    cpl_vector_delete(pupil_median_up);
    cpl_vector_delete(pupil_median_down);
    cpl_vector_delete(pupil_mediam_sort);
    cpl_array_delete(back_frames);
    
    CPLCHECK_MSG("Pupil Selection of good files failed");
    gravi_msg_function_exit(1);
    return CPL_ERROR_NONE;
}
    
/*----------------------------------------------------------------------------*/
/**
 * @brief Get the reference pixels for the pupil guiding on the acquisition camera
 *
 * @param header:   input header
 * @param diode_pos_subwindow:   output bivector, positions of the subwindows in pupil camera
 *
 * \exception CPL_ERROR_NULL_INPUT input data is missing
 * \exception CPL_ERROR_ILLEGAL_INPUT tel outside limits
 *
 * The output bi-vector is of size 4x4=16 (4 telescopes, 4 lenslets on the acquisition camera).
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_acqcam_get_pup_ref_v2 (cpl_propertylist * header, cpl_bivector *  diode_pos_subwindow)
{
    gravi_msg_function_start(0);
    cpl_ensure_code (header,          CPL_ERROR_NULL_INPUT);
    cpl_ensure_code (diode_pos_subwindow,          CPL_ERROR_NULL_INPUT);

    cpl_size nsx = 0, nsy = 0, sx = 0, sy = 0;
    
    cpl_vector * x_pos_subwindow = cpl_bivector_get_x (diode_pos_subwindow);
    cpl_vector * y_pos_subwindow = cpl_bivector_get_y (diode_pos_subwindow);

    
    for (int tel = 0 ; tel < GRAVI_SPOT_NTEL; tel++)
    {
        /* If sub-windowing, we read the sub-window size and
         * the sub-window start for pupil */
        if ( cpl_propertylist_has (header, "ESO DET1 FRAMES NX") ) {
            char name[90];
            
            nsx = cpl_propertylist_get_int (header, "ESO DET1 FRAMES NX");
            sprintf (name, "ESO DET1 FRAM%d STRX", 3*GRAVI_SPOT_NTEL + tel + 1);
            sx = cpl_propertylist_get_int (header, name);
            
            nsy = cpl_propertylist_get_int (header, "ESO DET1 FRAMES NY");
            sprintf (name, "ESO DET1 FRAM%d STRY", 3*GRAVI_SPOT_NTEL + tel + 1);
            sy = cpl_propertylist_get_int (header, name);
            
            CPLCHECK_MSG ("Cannot get sub-windowing parameters");
        }
        
        cpl_msg_info (cpl_func,"sub-window pupil %i sx= %lld sy = %lld", tel, sx, sy);
            
        /* Read the sub-apperture reference positions
         * Converted to accound for sub-windowing
         * In vector convention (start at 0,0) */
        
        for (int lens = 0; lens < GRAVI_SPOT_NLENS ; lens++) {
            double xv = gravi_pfits_get_ptfc_acqcam (header, lens*GRAVI_SPOT_NTEL + tel + 1);
            double yv = gravi_pfits_get_ptfc_acqcam (header, lens*GRAVI_SPOT_NTEL + tel + 17);
            cpl_vector_set( x_pos_subwindow, lens*GRAVI_SPOT_NTEL+tel , xv - (sx - tel*nsx) - 1);
            cpl_vector_set( y_pos_subwindow, lens*GRAVI_SPOT_NTEL+tel , yv - (sy - 3*nsy)   - 1);
    /*        cpl_msg_info (cpl_func,"pupil %lli subC %i = %10.4f,%10.4f",
                           tel, lens, x_diode_subwindow[lens][tel], y_diode_subwindow[lens][tel]);*/
            CPLCHECK_MSG ("Cannot get pupil reference position");
        }
    }

    gravi_msg_function_exit(0);
    return CPL_ERROR_NONE;
}


/*----------------------------------------------------------------------------*/
/**
* @brief Get the position of the telescope diodes in pixels
*
* @param header:   input header
 * @param good_frames:   array of 0 and 1, which tells when the pupil beacons are on
 * @param scale_vector:  output vector with the scales used for each telescope
 * @param diode_pos_telescope:   output bivector, positions of the diode in pixel unit (and camera orientation)
*
* \exception CPL_ERROR_NULL_INPUT input data is missing
*
* The output bi-vector is of size n_onx4x4=16 (4 telescopes, 4 diodes).
 * To note, because of the rotation of the parralactic angle (which can be large close to zenith), the output bivector is
 * actually an array of n_on bivector, where n_on is the number of frames with the beacons light on.
*/
/*----------------------------------------------------------------------------*/


cpl_error_code gravi_acqcam_get_diode_ref_v2 (cpl_propertylist * header,
                                              cpl_array * good_frames,
                                      cpl_vector * scale_vector,
                                      cpl_bivector **  diode_pos_telescope,
                                              int nrow_on)
{
    gravi_msg_function_start(0);
    cpl_ensure_code (header,          CPL_ERROR_NULL_INPUT);
    cpl_ensure_code (good_frames,          CPL_ERROR_NULL_INPUT);

    int nv = 0;
    double dx,dy,scale;
    double parang1 = cpl_propertylist_get_double(header, "ESO ISS PARANG START");
    double parang2 = cpl_propertylist_get_double(header, "ESO ISS PARANG END");
    CPLCHECK_MSG("Cannot determine parallactic angle");
    
    cpl_size nrow = cpl_array_get_size(good_frames);
    if (fabs(cpl_array_get_mean(good_frames)*nrow -nrow_on) > 1e-2)
        cpl_msg_error (cpl_func, "Ratio of blinking frames different %f from %i",cpl_array_get_mean(good_frames)*nrow,nrow_on);
    
    for (int tel = 0 ; tel < GRAVI_SPOT_NTEL; tel++)
    {
        
        /* Get the telescope name and ID */
        const char * telname = gravi_conf_get_telname (tel, header);

        /* Check telescope name */
        if (!telname) cpl_msg_error (cpl_func, "Cannot read the telescope name");
        cpl_ensure_code (telname, CPL_ERROR_ILLEGAL_INPUT);
        CPLCHECK_MSG("Cannot get telescope names");
        
        /* Hardcoded theoretical positions in mm */

        /* If UTs or ATs, select scaling, rotation, and spacing */
        if (telname[0] == 'U') {
        // FE 2019-05-15 replaced with median measured for the whole 2017-2018
            // Galactic Center data set, which should be the best information at hand
            // cpl_vector_set (output, GRAVI_SPOT_SCALE, 16.225);
        if ( tel == 0 ) scale = 16.83;
        if ( tel == 1 ) scale = 17.42;
        if ( tel == 2 ) scale = 16.85;
        if ( tel == 3 ) scale = 17.41;
        // below information could be calculated from diode position in header
            // this would also give the 45 offset angle introduced above to calculate
        // the spot angle
            dx = 0.363;
            dy = 0.823;
        } else if (telname[0] == 'A') {
        // FE 2019-05-15 maybe we should also update the AT numbers
            scale = 73.0154;
            // EW, FE 2019-09-11: short and long side of AT beacons are
            // flipped when compared to UTs
            dx = 0.158;
            dy = 0.122;
        } else
            return cpl_error_set_message (cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                          "Cannot get telescope name");
        CPLCHECK_MSG("Cannot get telescope scale");
        cpl_vector_set(scale_vector,tel,scale);
    
        /* We see that fitting the spot angle fails if one spot is fainter than reflection from brightest spot
         resulting in wrong pupil_x,y, therefore we calculate angle from header information
         by design of telescope spiders. This model angle should be
         angle = north_angle + paralactic angle + 45
         Remark: checking the angles for the GC easter flare night, we get an offset of 45.75 degrees */
        double northangle = gravi_pfits_get_northangle_acqcam(header, tel);
        CPLCHECK_MSG("Cannot determine field angle");

            
        cpl_size n_on= 0;
        for (cpl_size n =0; n < nrow; n++)
        {
            /*cpl_msg_info (cpl_func, "Totos: %lli, %.2f",cpl_array_get_int(good_frames,n,&nv),cpl_array_get_mean(good_frames));*/
            if (cpl_array_get_int(good_frames,n,&nv) ==1)
            {
                double padif = parang2 - parang1;
                if (padif > 180)
                    padif -= 360;
                if (padif < -180)
                    padif += 360;
                double parang= parang1 + padif*n/(nrow-1);
                
                // TODO: FE handle case of calibration unit //
                double angle = northangle + parang + 45.;
                if (angle < 0)   angle += 180;
                if (angle > 180) angle -= 180;
                double cang = cos(angle * CPL_MATH_RAD_DEG) * scale;
                double sang = sin(angle * CPL_MATH_RAD_DEG) * scale;
                
                /* Diode arrangement */
                
                cpl_vector * x_pos_telescope = cpl_bivector_get_x (diode_pos_telescope[n_on]);
                cpl_vector * y_pos_telescope = cpl_bivector_get_y (diode_pos_telescope[n_on]);
                
                cpl_vector_set(x_pos_telescope, 0*GRAVI_SPOT_NTEL+tel, -cang * dx + sang * dy);
                cpl_vector_set(x_pos_telescope, 1*GRAVI_SPOT_NTEL+tel, -cang * dx - sang * dy);
                cpl_vector_set(x_pos_telescope, 2*GRAVI_SPOT_NTEL+tel,  cang * dx - sang * dy);
                cpl_vector_set(x_pos_telescope, 3*GRAVI_SPOT_NTEL+tel,  cang * dx + sang * dy);
                
                cpl_vector_set(y_pos_telescope, 0*GRAVI_SPOT_NTEL+tel,  cang * dy + sang * dx);
                cpl_vector_set(y_pos_telescope, 1*GRAVI_SPOT_NTEL+tel, -cang * dy + sang * dx);
                cpl_vector_set(y_pos_telescope, 2*GRAVI_SPOT_NTEL+tel, -cang * dy - sang * dx);
                cpl_vector_set(y_pos_telescope, 3*GRAVI_SPOT_NTEL+tel,  cang * dy - sang * dx);
                
                n_on ++;
                
                if (n_on == 10)
                {
                    cpl_msg_info (cpl_func, "angle, parang %.2f, %.2f",angle,parang);
                    cpl_msg_info (cpl_func, "dx and dy %.2f, %.2f",cang,sang);
                }
                    
                
            }
        }
    
    
    CPLCHECK_MSG("Cannot determine diode position in telescope space");
    }
    
    gravi_msg_function_exit(0);
    return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
* @brief Get the position of the telescope diodes in pixels
*
 * @param diode_pos_subwindow:   input bivector, positions of the subwindows in pupil camera
 * @param diode_pos_telescope:   input bivector, positions of the diode in pixel unit (and camera orientation)
 * @param diode_pos_theoretical:   output bivectori, expected position of the beacons on the acq camera
 .
*
* \exception CPL_ERROR_NULL_INPUT input data is missing
*
* The output bi-vector is of size n_onx11x4x4x4=16 (4 telescopes, 4 diodes, 4 lenslet).
 * To note, because of the rotation of the parralactic angle (which can be large close to zenith), the output bivector is
 * actually an array of n_on bivector, where n_on is the number of frames with the beacons light on.
 * Also To note, the position depends on the focus value. Therefore, we will propose 11 different positions for 11 values of defocs.
 * The number 11 is actually a parameter defined by GRAVI_SPOT_NFOCUS.
 * The scaling factor is defined by the function gravi_acqcam_defocus_scaling (A thrid order polynomial).
*/
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_acqcam_get_diode_theoretical_v2(cpl_bivector *  diode_pos_subwindow,
                                                     cpl_bivector **  diode_pos_telescope,
                                                     cpl_bivector **  diode_pos_theoretical,
                                                     cpl_size nrow_on, int ury)
{
    
    gravi_msg_function_start(1);
    cpl_ensure_code (diode_pos_subwindow,          CPL_ERROR_NULL_INPUT);
    cpl_ensure_code (diode_pos_telescope,          CPL_ERROR_NULL_INPUT);
    cpl_ensure_code (diode_pos_theoretical,          CPL_ERROR_NULL_INPUT);
    
    cpl_vector * x_pos_subwindow = cpl_bivector_get_x (diode_pos_subwindow);
    cpl_vector * y_pos_subwindow = cpl_bivector_get_y (diode_pos_subwindow);
    
    for (int tel = 0 ; tel < GRAVI_SPOT_NTEL; tel++)
    {
        double x_lenslet_mean=0.0;
        double y_lenslet_mean=0.0;
        for (int lens = 0 ; lens < GRAVI_SPOT_NLENS; lens++)
        {
            x_lenslet_mean+=cpl_vector_get(x_pos_subwindow,lens*GRAVI_SPOT_NTEL+tel)/GRAVI_SPOT_NLENS;
            y_lenslet_mean+=cpl_vector_get(y_pos_subwindow,lens*GRAVI_SPOT_NTEL+tel)/GRAVI_SPOT_NLENS;
        }
        
        cpl_msg_info (cpl_func, "Reference pixel position for tel %d : X = %.2f, Y= %.2f", tel, x_lenslet_mean, y_lenslet_mean);
        
        for (int lens = 0 ; lens < GRAVI_SPOT_NLENS; lens++)
        for (int focus = 0 ; focus < GRAVI_SPOT_NFOCUS; focus++)
        {
            double xlenslet=cpl_vector_get(x_pos_subwindow,lens*GRAVI_SPOT_NTEL+tel)+
            (cpl_vector_get(x_pos_subwindow,lens*GRAVI_SPOT_NTEL+tel)-x_lenslet_mean)*gravi_acqcam_defocus_scaling(focus);
            double ylenslet=cpl_vector_get(y_pos_subwindow,lens*GRAVI_SPOT_NTEL+tel)+
            (cpl_vector_get(y_pos_subwindow,lens*GRAVI_SPOT_NTEL+tel)-y_lenslet_mean)*gravi_acqcam_defocus_scaling(focus);
            
            for (int spot = 0 ; spot < GRAVI_SPOT_NSPOT; spot++)
            for (cpl_size n_on = 0 ; n_on < nrow_on; n_on++)
            {
                
                cpl_vector * x_pos_telescope = cpl_bivector_get_x (diode_pos_telescope[n_on]);
                cpl_vector * y_pos_telescope = cpl_bivector_get_y (diode_pos_telescope[n_on]);
                cpl_vector * x_pos_theoretical = cpl_bivector_get_x (diode_pos_theoretical[n_on]);
                cpl_vector * y_pos_theoretical = cpl_bivector_get_y (diode_pos_theoretical[n_on]);
                
                double x_theo = xlenslet + cpl_vector_get(x_pos_telescope,spot*GRAVI_SPOT_NTEL+tel);
                double y_theo = ylenslet + cpl_vector_get(y_pos_telescope,spot*GRAVI_SPOT_NTEL+tel)- ury +1;
                
                cpl_vector_set(x_pos_theoretical,((focus*GRAVI_SPOT_NSPOT+spot)*GRAVI_SPOT_NLENS+lens)*GRAVI_SPOT_NTEL+tel,x_theo);
                cpl_vector_set(y_pos_theoretical,((focus*GRAVI_SPOT_NSPOT+spot)*GRAVI_SPOT_NLENS+lens)*GRAVI_SPOT_NTEL+tel,y_theo);
            }
        }
    }
    
CPLCHECK_MSG("Cannot compute the theoretical coordinates of spots");
gravi_msg_function_exit(1);
return CPL_ERROR_NONE;
    
}

/*----------------------------------------------------------------------------*/
/**
* @brief imprint zeros at the postion of the detected pupil spots
*
 * @param mean_img:   input / output image,
 * @param diode_pos_offset:   input bivectori, the initial pupil shift caused by the shift and add algorithm
 * @param diode_pos_theoretical:   input bivectori, expected position of the beacons on the acq camera
 * @param ury:  separation between pupil images and rest of acq camera (typically 750)
 .
*/
/*----------------------------------------------------------------------------*/


cpl_error_code gravi_acqcam_spot_imprint_v2(cpl_image *mean_img, cpl_bivector ** diode_pos_offset, cpl_bivector **  diode_pos_theoretical, int ury)
{
    gravi_msg_function_start(1);
    
    cpl_ensure_code(mean_img, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code(diode_pos_offset, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code(diode_pos_theoretical, CPL_ERROR_NULL_INPUT);
    
    cpl_size nx = cpl_image_get_size_x (mean_img);
    cpl_size ny = cpl_image_get_size_y (mean_img);
    
    /* imprint position for first image */
    int n_on = 0;
    cpl_vector * x_pos_offset = cpl_bivector_get_x (diode_pos_offset[n_on]);
    cpl_vector * y_pos_offset = cpl_bivector_get_y (diode_pos_offset[n_on]);
    cpl_vector * x_pos_theoretical = cpl_bivector_get_x (diode_pos_theoretical[n_on]);
    cpl_vector * y_pos_theoretical = cpl_bivector_get_y (diode_pos_theoretical[n_on]);
    
    /* imprint position with no defocus */
    int focus = GRAVI_SPOT_NFOCUS/2;
    
    /* Loop on diode and appertures */
    for (int tel = 0 ; tel < GRAVI_SPOT_NTEL; tel++)
    for (int lens = 0 ; lens < GRAVI_SPOT_NLENS; lens++)
    for (int spot = 0 ; spot < GRAVI_SPOT_NSPOT; spot++)
    {
        double x_off = cpl_vector_get(x_pos_offset,tel);
        double y_off = cpl_vector_get(y_pos_offset,tel)+ury;
        double x_diode = cpl_vector_get(x_pos_theoretical,((focus*GRAVI_SPOT_NSPOT+spot)*GRAVI_SPOT_NLENS+lens)*GRAVI_SPOT_NTEL+tel);
        double y_diode = cpl_vector_get(y_pos_theoretical,((focus*GRAVI_SPOT_NSPOT+spot)*GRAVI_SPOT_NLENS+lens)*GRAVI_SPOT_NTEL+tel);
        cpl_size xf = roundl(x_diode + x_off) + 1;
        cpl_size yf = roundl(y_diode + y_off) + 1;
        if (xf < 2 || xf > nx-2 || yf < 2 || yf > ny-2) continue;
        cpl_image_set (mean_img, xf,   yf, 0);
        cpl_image_set (mean_img, xf-1, yf, 0);
        cpl_image_set (mean_img, xf+1, yf, 0);
        cpl_image_set (mean_img, xf, yf+1, 0);
        cpl_image_set (mean_img, xf, yf-1, 0);
        CPLCHECK ("Cannot imprint cross in image");
    }
    
    CPLCHECK ("Cannot imprint cross in image");
    gravi_msg_function_exit(1);
    return CPL_ERROR_NONE;
}


/*----------------------------------------------------------------------------*/
/**
* @brief perform shift and add of the pupil image according to the theoretical positions
*
 * @param pupilImage_onFrames:   input imagelist, pupil images where the beacons are on
 * @param pupilImage_shiftandadd:   output imagelist, pupil images after shift and add
 * @param good_frames:   input array of 0 and 1, which tells when the pupil beacons are on
 * @param focus_value:   output vector, store the measured focus value.
 * @param diode_pos_theoretical:   input bivectori, expected position of the beacons on the acq camera
 * @param diode_pos_offset:   output bivectori, the initial pupil shift caused by the shift and add algorithm
 * @param nrow_on:   the number of input frames
 .
*
* \exception CPL_ERROR_NULL_INPUT input data is missing
*
 *
 * The first step of this routine is to find the focus/defocus value.
 * To do so, it use brut force to find the focus wich gives the maximum flux (over a 55x55 add and shifted image)
 *
 * The output images list are the pupil plane images, but shift and added according to the theoretical positions,
 * and the maximum focus.
 * In other words, it co-adds all the beacon images on a single image, which has a smaller size.
 * the size of the shigt and added images is 59*2+1 (that is fixed by the parameter GRAVI_SPOT_NSEARCH)
 * It means that if the beacons moved by more tha 59 pixels, we will not find them.
 * But the margin should be big enough to cover all cases.
 *
* The output bi-vector is of size n_onx4  (4 telescopes, n_on frames).
 * They correspond to the x and y pupil shift observe at each 'beacon ON' frames
 * Not the the diode_pos_offset bivector is just an initialisation (different from zero) resulting from the cut of the shift
 * and add algorithm
*/
/*----------------------------------------------------------------------------*/


cpl_error_code gravi_acqcam_perform_shiftandadd_v2(cpl_imagelist * pupilImage_onFrames, cpl_imagelist ** pupilImage_shiftandadd, cpl_array * good_frames,
                                                   cpl_vector * focus_value,
                                                     cpl_bivector **  diode_pos_theoretical ,
                                                     cpl_bivector **  diode_pos_offset, cpl_size nrow_on)
{
    gravi_msg_function_start(1);
    
    cpl_ensure_code(pupilImage_onFrames, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code(pupilImage_shiftandadd, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code(good_frames, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code(focus_value, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code(diode_pos_theoretical, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code(diode_pos_offset, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code(nrow_on, CPL_ERROR_NULL_INPUT);
    
    if (cpl_imagelist_get_size(pupilImage_onFrames) != nrow_on)
        cpl_msg_error (cpl_func, "Problem with number of blinking frames");
    
    for (int tel = 0 ; tel < GRAVI_SPOT_NTEL; tel++)
    {
        cpl_vector * focus_max = cpl_vector_new(GRAVI_SPOT_NFOCUS);
        for (cpl_size focus = 0 ; focus < GRAVI_SPOT_NFOCUS; focus++)
        {
            cpl_imagelist * pupilImage_shifted= cpl_imagelist_new();
            cpl_size n_images = 0;
            for (cpl_size n = 0 ; n < nrow_on; n++)
            for (int lens = 0 ; lens < GRAVI_SPOT_NLENS; lens++)
            for (int spot = 0 ; spot < GRAVI_SPOT_NSPOT; spot++)
            {
                cpl_vector * x_pos_theoretical = cpl_bivector_get_x (diode_pos_theoretical[n]);
                cpl_vector * y_pos_theoretical = cpl_bivector_get_y (diode_pos_theoretical[n]);
                
                double x_theo=cpl_vector_get(x_pos_theoretical,((focus*GRAVI_SPOT_NSPOT+spot)*GRAVI_SPOT_NLENS+lens)*GRAVI_SPOT_NTEL+tel);
                double y_theo=cpl_vector_get(y_pos_theoretical,((focus*GRAVI_SPOT_NSPOT+spot)*GRAVI_SPOT_NLENS+lens)*GRAVI_SPOT_NTEL+tel);
                
                cpl_size x=(int) (x_theo + 0.5);
                cpl_size y=(int) (y_theo+ 0.5);
                cpl_image * verysmall_img_tmp = cpl_imagelist_get (pupilImage_onFrames, n);
                cpl_image * verysmall_img = gravi_image_extract( verysmall_img_tmp, x+1-GRAVI_SPOT_NSEARCH, y+1-GRAVI_SPOT_NSEARCH, x+1+GRAVI_SPOT_NSEARCH, y+1+GRAVI_SPOT_NSEARCH);
                cpl_imagelist_set(pupilImage_shifted, verysmall_img, n_images);
                n_images ++;
            }
            
            cpl_image * image_mean = cpl_imagelist_collapse_create (pupilImage_shifted);
            cpl_vector_set(focus_max, focus, cpl_image_get_max(image_mean));
            cpl_imagelist_delete(pupilImage_shifted);
            cpl_image_delete(image_mean);
        }
    
        int focus_max_pos = cpl_vector_get_maxpos (focus_max);
        cpl_msg_info (cpl_func, "Focus value for telescope %d : F = %.2f pixels", tel, GRAVI_SPOT_SWINDOW * gravi_acqcam_defocus_scaling(focus_max_pos));
        cpl_vector_set(focus_value,tel,gravi_acqcam_defocus_scaling(focus_max_pos)*100.);
        CPLCHECK_MSG("Cannot find optimum focus position");
        
        for (cpl_size n = 0 ; n < nrow_on; n++)
        {
            cpl_vector * x_pos_theoretical = cpl_bivector_get_x (diode_pos_theoretical[n]);
            cpl_vector * y_pos_theoretical = cpl_bivector_get_y (diode_pos_theoretical[n]);
            cpl_imagelist * pupilImage_shifted= cpl_imagelist_new();
            cpl_vector *  offs_x =   cpl_vector_new (GRAVI_SPOT_NLENS*GRAVI_SPOT_NSPOT);
            cpl_vector *  offs_y =   cpl_vector_new (GRAVI_SPOT_NLENS*GRAVI_SPOT_NSPOT);
            cpl_size n_images = 0;
            for (int lens = 0 ; lens < GRAVI_SPOT_NLENS; lens++)
            for (int spot = 0 ; spot < GRAVI_SPOT_NSPOT; spot++)
            {
                double x_theo=cpl_vector_get(x_pos_theoretical,((focus_max_pos*GRAVI_SPOT_NSPOT+spot)*GRAVI_SPOT_NLENS+lens)*GRAVI_SPOT_NTEL+tel);
                double y_theo=cpl_vector_get(y_pos_theoretical,((focus_max_pos*GRAVI_SPOT_NSPOT+spot)*GRAVI_SPOT_NLENS+lens)*GRAVI_SPOT_NTEL+tel);
                int x= (int) (x_theo+ 0.5);
                int y= (int) (y_theo + 0.5);
                cpl_image * verysmall_img_tmp = cpl_imagelist_get (pupilImage_onFrames, n);
                cpl_image * verysmall_img = gravi_image_extract( verysmall_img_tmp, x+1-GRAVI_SPOT_NSEARCH, y+1-GRAVI_SPOT_NSEARCH, x+1+GRAVI_SPOT_NSEARCH, y+1+GRAVI_SPOT_NSEARCH);
                /*if (n==10)
                {
                    cpl_msg_info (cpl_func, "Cutting at position : X/Y = %lli/%lli, [%lli,%lli,%lli]", x,y,lens,spot,tel);
                    cpl_msg_info (cpl_func, "image_tmp %lli shifted : X/Y = %lli/%lli/%lli 15/10 %.2f", n,tel,lens,spot,cpl_image_get(verysmall_img_tmp, 11, 16, &nv));
                    cpl_msg_info (cpl_func, "image %lli shifted : X/Y = %lli/%lli/%lli 15/10 %.2f", nrow_on,tel,lens,spot,cpl_image_get(verysmall_img, 11, 16, &nv));
                }*/
                /*if (n == 10)
                cpl_msg_info (cpl_func, "image %lli shifted : X/Y = %lli/%lli/%lli 15/10 %.2f -- off %.2f", n,tel,lens,spot,cpl_image_get(verysmall_img, 11, 16, &nv),x_diode_theo_fine[n][focus_max_pos][spot][lens][tel]);*/
                cpl_imagelist_set(pupilImage_shifted, verysmall_img, n_images);
                cpl_vector_set (offs_x, n_images, (x-x_theo));
                cpl_vector_set (offs_y, n_images, (y-y_theo));
                
                n_images ++;
            }
            
            cpl_bivector *  offs =   cpl_bivector_wrap_vectors (offs_x, offs_y);
            
            double ppos_x, ppos_y;
            cpl_image** cpl_image_combined =  cpl_geom_img_offset_saa    (pupilImage_shifted,offs,CPL_KERNEL_DEFAULT,
                                                                          0,0,CPL_GEOM_INTERSECT,&ppos_x,&ppos_y);
            CPLCHECK_MSG("Cannot do fine shift and add");
            
            cpl_vector * x_pos_offset = cpl_bivector_get_x (diode_pos_offset[n]);
            cpl_vector * y_pos_offset = cpl_bivector_get_y (diode_pos_offset[n]);
            
            cpl_image * image_mean = cpl_image_extract (cpl_image_combined[0], 4, 4, GRAVI_SPOT_NSEARCH*2-3, GRAVI_SPOT_NSEARCH*2-3);
            
            cpl_vector_set(x_pos_offset, tel, ppos_x+3);
            cpl_vector_set(y_pos_offset, tel, ppos_y+3);
                
            cpl_imagelist_set(pupilImage_shiftandadd[tel], image_mean, n);
            CPLCHECK_MSG("Cannot add shift and added image to imagelist");
            
            if (cpl_image_combined[0] != NULL) cpl_image_delete(cpl_image_combined[0]);
            if (cpl_image_combined[1] != NULL) cpl_image_delete(cpl_image_combined[1]);
            cpl_free(cpl_image_combined);
            
            cpl_bivector_delete(offs);
            cpl_imagelist_delete(pupilImage_shifted);
        }
        cpl_vector_delete(focus_max);
    }
    
    CPLCHECK_MSG("Fail at running the shift and add ESO algorithm");
    gravi_msg_function_exit(1);
    
    return CPL_ERROR_NONE;
}


/*----------------------------------------------------------------------------*/
/**
* @brief on the images shift and added, find maximum and perform gaussian fit
*
 * @param pupilImage_shiftandadd:   input imagelist, pupil images after shift and add
 * @param bad_frames_short:  output array of (0/1) int which tells if the pupil offset was found.
 * @param o_header:          output header
 * @param diode_pos_offset:   output bivectori, the initial pupil shift caused by the shift and add algorithm
 * @param nrow_on:   the number of input frames
 .
*
* \exception CPL_ERROR_NULL_INPUT input data is missing
*
 *
* The output bi-vector is of size n_onx4  (4 telescopes, n_on frames).
 * They correspond to the x and y pupil shift observe at each 'beacon ON' frames
 * The diode_pos_offset bivector correspond to the offset of the pupil with respect to the theoretical position
*/
/*----------------------------------------------------------------------------*/


cpl_error_code gravi_acqcam_get_pupil_offset_v2(cpl_imagelist ** pupilImage_shiftandadd,
                                                     cpl_array * bad_frames_short,
                                                cpl_bivector **  diode_pos_offset, cpl_propertylist * o_header,
                                                cpl_size nrow_on)
{
    gravi_msg_function_start(1);
    
    cpl_ensure_code(pupilImage_shiftandadd, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code(bad_frames_short, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code(diode_pos_offset, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code(o_header, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code(nrow_on, CPL_ERROR_NULL_INPUT);
    
    int nv =0;
    
    if (cpl_array_get_size(bad_frames_short)!= nrow_on)
        cpl_msg_error (cpl_func, "Problem with number of blinking frames");
    
    for (int tel = 0 ; tel < GRAVI_SPOT_NTEL; tel++)
    {
        double previous_xsc=0.0;
        double previous_ysc=0.0;
        double flux_sum=0.0;
        int n_sum=0;
    
        for (cpl_size n = 0 ; n < nrow_on; n++)
        {
            cpl_image * image = cpl_imagelist_get (pupilImage_shiftandadd[tel], n);
            CPLCHECK_MSG("Failure at getting saa image");
            
            cpl_size px, py;
            cpl_image_get_maxpos (image, &px, &py);
            CPLCHECK_MSG("Failure to read maximum position");
            
            double flux_max= cpl_image_get(image,px,py,&nv);
            double std = cpl_image_get_stdev  (image);
            
            
            double xsc = (double) px;
            double ysc = (double) py;
            double exsc=3., eysc=3.;
            cpl_size size = 15;
            
            gravi_acq_fit_gaussian (image, &xsc, &ysc, &exsc, &eysc, size);
            CPLCHECK_MSG("Failure at fitting pupil spot");
            
            cpl_vector * x_pos_offset = cpl_bivector_get_x (diode_pos_offset[n]);
            cpl_vector * y_pos_offset = cpl_bivector_get_y (diode_pos_offset[n]);
            CPLCHECK_MSG("Failure at reading bivector offsets");
            
            double x_offset_final= cpl_vector_get(x_pos_offset, tel) + xsc-GRAVI_SPOT_NSEARCH-1;
            double y_offset_final= cpl_vector_get(y_pos_offset, tel) + ysc-GRAVI_SPOT_NSEARCH-1;
            
            cpl_vector_set(x_pos_offset, tel, x_offset_final);
            cpl_vector_set(y_pos_offset, tel, y_offset_final);
            
            int previous;
            if ((flux_max/std < 4)||(exsc<0))
            {
                previous = cpl_array_get_int (bad_frames_short,n,&nv);
                cpl_array_set_int (bad_frames_short,n,previous|(1<<tel));
            } else {
            if (n>0)
                {
                    double distance=(xsc-previous_xsc)*(xsc-previous_xsc)+(ysc-previous_ysc)*(ysc-previous_ysc);
                    if (distance > 10*10)
                    {
                        previous = cpl_array_get_int (bad_frames_short,n-1,&nv);
                        cpl_array_set_int (bad_frames_short,n-1,previous|(1<<(tel+4)));
                        previous = cpl_array_get_int (bad_frames_short,n,&nv);
                        cpl_array_set_int (bad_frames_short,n,previous|(1<<(tel+4)));
                    }
                }
                flux_sum+=flux_max;
                n_sum+=1;
            }
            previous_xsc=xsc;
            previous_ysc=ysc;
            CPLCHECK_MSG("Failure at testing bad pupil measurments");
        }
        
        char qc_name[100];
        double flux_mean=flux_sum/(1e-7+n_sum);
        /* add QC parameters with pupil beacon flux */
            sprintf(qc_name, "ESO QC ACQ PUP%i SPOT_FLUX", tel + 1);
            cpl_msg_info(cpl_func, "%s = %f", qc_name,flux_mean);
            cpl_propertylist_update_double(o_header, qc_name, flux_mean);
            cpl_propertylist_set_comment(o_header, qc_name,"[ADU] mean spot flux");
        
        CPLCHECK_MSG("Failing to store QC value");
    }
    
    CPLCHECK_MSG("Failure at fitting gaussian on diode position");
    
    gravi_msg_function_exit(1);
    return CPL_ERROR_NONE;
}


/*----------------------------------------------------------------------------*/
/**
* @brief store the pupil offset in the table, and create QC parameters.
*
 * @param acqcam_table:   input /output table
 * @param header:          input header
 * @param o_header:          output header
 * @param good_frames:  input array of (0/1) int which tells if the pupil beacons are "on"
 * @param scale_vector: optional --   input  bivectori, the initial pupil shift caused by the shift and add algorithm
 * @param bad_frames_short:  optional --  input array of (0/1) int which tells if the pupil offset was found.
 * @param diode_pos_offset:  optional --   input bivectori, the initial pupil shift caused by the shift and add algorithm
 * @param focus_value:    optional -- input bivectori, the initial pupil shift caused by the shift and add algorithm
 * @param static_param_data:  optional --  static parameter table
 .
*
* \exception CPL_ERROR_NULL_INPUT input data is missing
*
 *
* Populate the table of in pupil offset (and u-v plane shift).
 * if no good_frames (sum(good_frames==0), than all the other parameters can be put to NULL
 * In that case, the table will be filled with zeros.
*/
/*----------------------------------------------------------------------------*/



cpl_error_code    gravi_acqcam_set_pupil_table_v2(cpl_table * acqcam_table, cpl_propertylist * header, cpl_propertylist * o_header, cpl_array * good_frames,
                                                  cpl_vector* scale_vector, cpl_array * bad_frames_short, cpl_bivector **  diode_pos_offset, cpl_vector * focus_value, gravi_data *static_param_data)
{

gravi_msg_function_start(1);
int nv =0;
    
cpl_ensure_code(acqcam_table, CPL_ERROR_NULL_INPUT);
cpl_ensure_code(header, CPL_ERROR_NULL_INPUT);
cpl_ensure_code(o_header, CPL_ERROR_NULL_INPUT);
cpl_ensure_code(good_frames, CPL_ERROR_NULL_INPUT);
cpl_size nrow = cpl_array_get_size(good_frames);
    
/* Pupil positions array  */
gravi_table_new_column(acqcam_table, "PUPIL_NSPOT", NULL, CPL_TYPE_INT);
gravi_table_new_column(acqcam_table, "PUPIL_X", "pix", CPL_TYPE_DOUBLE);
gravi_table_new_column(acqcam_table, "PUPIL_Y", "pix", CPL_TYPE_DOUBLE);
gravi_table_new_column(acqcam_table, "PUPIL_Z", "pix", CPL_TYPE_DOUBLE);
gravi_table_new_column(acqcam_table, "PUPIL_R", "deg", CPL_TYPE_DOUBLE);
gravi_table_new_column(acqcam_table, "PUPIL_U", "m", CPL_TYPE_DOUBLE);
gravi_table_new_column(acqcam_table, "PUPIL_V", "m", CPL_TYPE_DOUBLE);
gravi_table_new_column(acqcam_table, "PUPIL_W", "m", CPL_TYPE_DOUBLE);
gravi_table_new_column(acqcam_table, "OPD_PUPIL", "m", CPL_TYPE_DOUBLE);
    
double sobj_x = gravi_pfits_get_sobj_x(header);
double sobj_y = gravi_pfits_get_sobj_y(header);
CPLCHECK_MSG("Cannot determine SOBJ X and Y");
    
if (cpl_array_get_mean(good_frames)*nrow<.5)
{
    /* no good frames */
    /* filling the table with zeros */
    cpl_msg_info(cpl_func,"No good pupil images ==> no good pupil position available");
    
    for (int tel = 0; tel < GRAVI_SPOT_NTEL; tel++)
    for (cpl_size row = 0; row < nrow; row++)
    {
        cpl_table_set(acqcam_table, "PUPIL_NSPOT", row * GRAVI_SPOT_NTEL + tel, 0);
        CPLCHECK_MSG("Cannot put 0 in NSPOT data in the ACQ PUPIL");
    }
    
    gravi_msg_function_exit(1);
    return CPL_ERROR_NONE;
}
    
    /* looks like the table will not be empty. Checking additional parameters*/
cpl_ensure_code(scale_vector, CPL_ERROR_NULL_INPUT);
cpl_ensure_code(bad_frames_short, CPL_ERROR_NULL_INPUT);
cpl_ensure_code(diode_pos_offset, CPL_ERROR_NULL_INPUT);
cpl_ensure_code(focus_value, CPL_ERROR_NULL_INPUT);
cpl_ensure_code(static_param_data, CPL_ERROR_NULL_INPUT);

cpl_size nrow_on = cpl_array_get_size(bad_frames_short);

if (fabs(cpl_array_get_mean(good_frames)*nrow -nrow_on) > 1e-2)
    cpl_msg_error (cpl_func, "Ratio of blinking frames different %f from %lli",cpl_array_get_mean(good_frames)*nrow,nrow_on);

    

for (int tel = 0; tel < GRAVI_SPOT_NTEL; tel++)
{
    int n_on=0;
    
    /* Get the conversion angle xy to uv in [rad] */
    double northangle = gravi_pfits_get_northangle_acqcam(header, tel);
    double cfangle = cos(northangle * CPL_MATH_RAD_DEG);
    double sfangle = sin(northangle * CPL_MATH_RAD_DEG);
    double scale = cpl_vector_get(scale_vector,tel);
    CPLCHECK_MSG("Cannot determine field angle");
    
    /* values which are valids for all frames */
    double r_shift = 0.0;
    double z_shift = cpl_vector_get(focus_value,tel)/GRAVI_SPOT_SWINDOW;
    double w_shift = gravi_acqcam_z2meter(z_shift, static_param_data);
    
    /*values that does change */
    double x_shift = 0.0;
    double y_shift = 0.0;
    double u_shift = 0.0;
    double v_shift = 0.0;
    double opd_pupil = 0.0;
    
    /* buffers to average values */
    double x_shift_sum = 0.0;
    double y_shift_sum = 0.0;
    double u_shift_sum = 0.0;
    double v_shift_sum = 0.0;
    int n_sum = 0;
    int n_bad_snr = 0;
    int n_bad_distance = 0;
    
    for (cpl_size row = 0; row < nrow; row++)
    {
        /*cpl_msg_info (cpl_func, "Row number %lli n %lli",row,n_on);*/
                    
        if (cpl_array_get_int(good_frames,row,&nv) != 1)
        {
            cpl_table_set(acqcam_table, "PUPIL_NSPOT", row * GRAVI_SPOT_NTEL + tel, 0);
            CPLCHECK_MSG("Cannot put 0 in NSPOT data in the ACQ PUPIL");
        }
        else
        {
            
            int is_it_bad=cpl_array_get_int (bad_frames_short,n_on,&nv);
            
            cpl_vector * x_pos_offset = cpl_bivector_get_x (diode_pos_offset[n_on]);
            cpl_vector * y_pos_offset = cpl_bivector_get_y (diode_pos_offset[n_on]);
            CPLCHECK_MSG("Cannot put data from the offset bivector");

            x_shift = cpl_vector_get(x_pos_offset,tel);
            y_shift = cpl_vector_get(y_pos_offset,tel);

            /* In UV [m] */
            u_shift = (cfangle * x_shift - sfangle * y_shift)
                    / scale;
            v_shift = (sfangle * x_shift + cfangle * y_shift)
                    / scale;
            opd_pupil = -(u_shift * sobj_x + v_shift * sobj_y)
                    * GRAVI_MATH_RAD_MAS;
            
            CPLCHECK_MSG("Cannot prepare data to be put in the ACQ PUPIL table");
            
            
            if (CHECK_BIT(is_it_bad,tel)||CHECK_BIT(is_it_bad,tel+4))
            {
                cpl_table_set(acqcam_table, "PUPIL_NSPOT", row * GRAVI_SPOT_NTEL + tel, 0);
                
                /* increasing counters */
                if (CHECK_BIT(is_it_bad,tel  )) n_bad_snr += 1;
                if (CHECK_BIT(is_it_bad,tel+4)) n_bad_distance += 1;
                u_shift = 0.0;
                v_shift = 0.0;
                opd_pupil = 0.0;
                cpl_msg_info(cpl_func,"Pupil image number%4lli is un-usable for tel %i",row,tel);
            }
            else
            {
                cpl_table_set(acqcam_table, "PUPIL_NSPOT", row * GRAVI_SPOT_NTEL + tel,
                        16);
                
                /* averaging positions */
                x_shift_sum += x_shift;
                y_shift_sum += y_shift;
                u_shift_sum += u_shift;
                v_shift_sum += v_shift;
                
                /* increasing counter */
                n_sum += 1;
            }
            
            cpl_table_set(acqcam_table, "PUPIL_X", row * GRAVI_SPOT_NTEL + tel,
                    x_shift);
            cpl_table_set(acqcam_table, "PUPIL_Y", row * GRAVI_SPOT_NTEL + tel,
                    y_shift);
            cpl_table_set(acqcam_table, "PUPIL_Z", row * GRAVI_SPOT_NTEL + tel,
                    z_shift);
            cpl_table_set(acqcam_table, "PUPIL_R", row * GRAVI_SPOT_NTEL + tel,
                    r_shift);
            cpl_table_set(acqcam_table, "PUPIL_U", row * GRAVI_SPOT_NTEL + tel,
                    u_shift);
            cpl_table_set(acqcam_table, "PUPIL_V", row * GRAVI_SPOT_NTEL + tel,
                    v_shift);
            cpl_table_set(acqcam_table, "PUPIL_W", row * GRAVI_SPOT_NTEL + tel,
                    w_shift);
            cpl_table_set(acqcam_table, "OPD_PUPIL", row * GRAVI_SPOT_NTEL + tel,
                    opd_pupil);
            CPLCHECK_MSG("Cannot put data in the ACQ PUPIL table");
                
            n_on++;
        }
        CPLCHECK_MSG("Cannot put data in the ACQ PUPIL table");
    }
        
        /* Add QC parameters */
        char qc_name[100];

        sprintf(qc_name, "ESO QC ACQ FIELD%i NORTH_ANGLE", tel + 1);
        cpl_msg_info(cpl_func, "%s = %f", qc_name, northangle);
        cpl_propertylist_update_double(o_header, qc_name, northangle);
        cpl_propertylist_set_comment(o_header, qc_name,
                "[deg] y->x, predicted North direction on ACQ");
    
        sprintf(qc_name, "ESO QC ACQ PUP%i FRAMES", tel + 1);
        cpl_msg_info(cpl_func, "%s = %i", qc_name, n_sum);
        cpl_propertylist_update_double(o_header, qc_name, n_sum);
        cpl_propertylist_set_comment(o_header, qc_name,
                "Good ACQ pupil frames");
                    
        sprintf(qc_name, "ESO QC ACQ PUP%i BADSNR", tel + 1);
        cpl_msg_info(cpl_func, "%s = %i", qc_name, n_bad_snr);
        cpl_propertylist_update_double(o_header, qc_name, n_bad_snr);
        cpl_propertylist_set_comment(o_header, qc_name,
                "Frames with low SNR");
                    
        sprintf(qc_name, "ESO QC ACQ PUP%i JUMP", tel + 1);
        cpl_msg_info(cpl_func, "%s = %i", qc_name, n_bad_distance);
        cpl_propertylist_update_double(o_header, qc_name, n_bad_distance);
        cpl_propertylist_set_comment(o_header, qc_name,
                "Frames with position jump");
    
        sprintf(qc_name, "ESO QC ACQ PUP%i SCALE", tel + 1);
        cpl_msg_info(cpl_func, "%s = %f", qc_name, scale);
        cpl_propertylist_update_double(o_header, qc_name, scale);
        cpl_propertylist_set_comment(o_header, qc_name,
                "[pix/m] diode scale on ACQ");
    
        sprintf(qc_name, "ESO QC ACQ PUP%i FOCUS", tel + 1);
        cpl_msg_info(cpl_func, "%s = %f", qc_name, z_shift);
        cpl_propertylist_update_int(o_header, qc_name, z_shift);
        cpl_propertylist_set_comment(o_header, qc_name,
                "[pix] defocus of pupil plane");
                
        /* to avoid division by zero */
        if (n_sum ==0) n_sum+=1;

        sprintf(qc_name, "ESO QC ACQ PUP%i XPOS", tel + 1);
        cpl_msg_info(cpl_func, "%s = %f", qc_name, x_shift_sum/n_sum);
        cpl_propertylist_update_double(o_header, qc_name, x_shift_sum/n_sum);
        cpl_propertylist_set_comment(o_header, qc_name,
                "[pix] pupil x-shift in ACQ");

        sprintf(qc_name, "ESO QC ACQ PUP%i YPOS", tel + 1);
        cpl_msg_info(cpl_func, "%s = %f", qc_name, y_shift_sum/n_sum);
        cpl_propertylist_update_double(o_header, qc_name, y_shift_sum/n_sum);
        cpl_propertylist_set_comment(o_header, qc_name,
                "[pix] pupil y-shift in ACQ");

        sprintf(qc_name, "ESO QC ACQ PUP%i UPOS", tel + 1);
        cpl_msg_info(cpl_func, "%s = %f", qc_name, u_shift_sum/n_sum);
        cpl_propertylist_update_double(o_header, qc_name, u_shift_sum/n_sum);
        cpl_propertylist_set_comment(o_header, qc_name,
                "[m] pupil u-shift in ACQ");

        sprintf(qc_name, "ESO QC ACQ PUP%i VPOS", tel + 1);
        cpl_msg_info(cpl_func, "%s = %f", qc_name, v_shift_sum/n_sum);
        cpl_propertylist_update_double(o_header, qc_name, v_shift_sum/n_sum);
        cpl_propertylist_set_comment(o_header, qc_name,
                "[m] pupil v-shift in ACQ");
        
    
    CPLCHECK_MSG("Cannot put data in the ACQ PUPIL table");
    
    
    
}

gravi_msg_function_exit(1);
return CPL_ERROR_NONE;
}



/*----------------------------------------------------------------------------*/
/**
* @brief extract sub window of image (similar to cpl_image_extract)
*
 * @param image_in:   input image
 * @param header:          input header
 * @param o_header:          output header
 * @param good_frames:  input array of (0/1) int which tells if the pupil beacons are "on"
 * @param scale_vector: optional --   input  bivectori, the initial pupil shift caused by the shift and add algorithm
 * @param bad_frames_short:  optional --  input array of (0/1) int which tells if the pupil offset was found.
 * @param diode_pos_offset:  optional --   input bivectori, the initial pupil shift caused by the shift and add algorithm
 * @param focus_value:    optional -- input bivectori, the initial pupil shift caused by the shift and add algorithm
 * @param static_param_data:  optional --  static parameter table
 .
 *
* Pupilate output image with zero if corrdinates outside range (instead of giving error).
*/
/*----------------------------------------------------------------------------*/



cpl_image * gravi_image_extract(cpl_image * image_in, cpl_size llx, cpl_size lly, cpl_size urx, cpl_size ury)
    {
        
        cpl_size  nx = cpl_image_get_size_x (image_in);
        cpl_size  ny = cpl_image_get_size_y (image_in);
        cpl_size llx_new,lly_new,urx_new,ury_new;
        cpl_size xpos, ypos;
        
        cpl_size nx_new=urx+1-llx;
        cpl_size ny_new=ury+1-lly;
        cpl_image * image_out= cpl_image_new (nx_new, ny_new, cpl_image_get_type(image_in));
        cpl_image_fill_window (image_out, 1, 1, nx_new, ny_new, 0.0);
        
        
        if ((llx >= nx)||(lly >= ny)||(urx <= 0)||(ury <= 0))
        {
            cpl_msg_warning (cpl_func, "Cutting at x=(%lli,%lli) y=(%lli,%lli) is outside the window bondaries", llx, urx, lly, ury);
        }
        else
        {
            if (llx < 1)
            {
                llx_new=1;
                xpos=1-llx;
            } else
            {
                llx_new=llx;
                xpos=1;
            }
            if (lly < 1)
            {
                lly_new=1;
                ypos=1-lly;
            } else
            {
                lly_new=lly;
                ypos=1;
            }
            
            if (urx > nx)
                urx_new=nx;
            else
                urx_new=urx;
            if (ury > ny)
                ury_new=ny;
            else
                ury_new=ury;
            
            cpl_image * image_in_cut = cpl_image_extract(image_in,llx_new,lly_new,urx_new,ury_new);
            
            cpl_image_copy (image_out, image_in_cut, xpos, ypos);
            
            cpl_image_delete(image_in_cut);
        }
        
        return image_out;
    }
    
/*----------------------------------------------------------------------------*/
/**
* @brief gives focus value as in a list of value
*
*/
/*----------------------------------------------------------------------------*/


double gravi_acqcam_defocus_scaling(int focus)
{
    double defocus;
    defocus = 0.4*(2.0*focus/(GRAVI_SPOT_NFOCUS-1) - 1.0)*(2.0*focus/(GRAVI_SPOT_NFOCUS-1) - 1.0)*(2.0*focus/(GRAVI_SPOT_NFOCUS-1) - 1.0);
    
    return defocus;
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
                                   cpl_propertylist * o_header,
                                   gravi_data *static_param_data)
{
    gravi_msg_function_start(1);
    cpl_ensure_code (mean_img,       CPL_ERROR_NULL_INPUT);
    cpl_ensure_code (acqcam_imglist, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code (header,         CPL_ERROR_NULL_INPUT);
    cpl_ensure_code (acqcam_table,   CPL_ERROR_NULL_INPUT);
    cpl_ensure_code (o_header,       CPL_ERROR_NULL_INPUT);
    
    char qc_name[100];
    int ntel = 4;

    int sts_mode = 0;

    /* check the feed mode */
    if ( !strcmp(gravi_pfits_get_feed (header), "DUAL_STS") ) {
        sts_mode = DUAL_STS;
    }
    else{
        sts_mode = SINGLE_STS;
    }

    /* Number of row */
    cpl_size nrow = cpl_imagelist_get_size (acqcam_imglist);

    /* MODE_ONAXIS or MODE_OFFAIS. We query once
       at the begining of the function */
    int axis_mode = gravi_pfits_get_axis (header);
    
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

    /* ----------------- START EKW 26/11/2018 read constant parameter from calibration file */
    /* Position of roof center on full frame */
    //double roof_x[] = {274.4, 787.1, 1236.1, 1673.4};
    //double roof_y[] = {242.3, 247.7, 225.8, 235.6};
    double *roof_x = NULL;
    double *roof_y = NULL;
    
    /* Position of single-field spot on full frame */
    //double spot_x[] = {289. , 798.2, 1245.5, 1696.};
    //double spot_y[] = {186.5, 187.5,  172.5,  178.};
    double *spot_x = NULL;
    double *spot_y = NULL;
    
    /* Default position angle of roof */
    //double roof_pos[] = {38.49, 38.54, 38.76, 39.80};
    double *roof_pos = NULL;

    cpl_table * roof_pos_table = gravi_data_get_table (static_param_data, "ROOFPOS");
    CPLCHECK_MSG ("STATIC_PARAM not available in the SOF. It is mandatory for acqcam reduction.");

     if ( cpl_table_has_column(roof_pos_table , "roof_x") ) {
         roof_x= cpl_table_get_data_double (roof_pos_table, "roof_x");
         cpl_msg_info(cpl_func,"roof_x [0] : %e", roof_x[0] );
         cpl_msg_info(cpl_func,"roof_x [1] : %e", roof_x[1] );
         cpl_msg_info(cpl_func,"roof_x [2] : %e", roof_x[2] );
         cpl_msg_info(cpl_func,"roof_x [3] : %e", roof_x[3] );
     }
     else {
       cpl_msg_warning(cpl_func,"Cannot get the default values for roof_x ");
     }

     if ( cpl_table_has_column(roof_pos_table , "roof_y") ) {
         roof_y= cpl_table_get_data_double (roof_pos_table, "roof_y");
         cpl_msg_info(cpl_func,"roof_y [0] : %e", roof_y[0] );
         cpl_msg_info(cpl_func,"roof_y [1] : %e", roof_y[1] );
         cpl_msg_info(cpl_func,"roof_y [2] : %e", roof_y[2] );
         cpl_msg_info(cpl_func,"roof_y [3] : %e", roof_y[3] );
     }
     else {
       cpl_msg_warning(cpl_func,"Cannot get the default values for roof_y ");
     }

     if ( cpl_table_has_column(roof_pos_table , "spot_x") ) {
         spot_x= cpl_table_get_data_double (roof_pos_table, "spot_x");
         cpl_msg_info(cpl_func,"spot_x [0] : %e", spot_x[0] );
         cpl_msg_info(cpl_func,"spot_x [1] : %e", spot_x[1] );
         cpl_msg_info(cpl_func,"spot_x [2] : %e", spot_x[2] );
         cpl_msg_info(cpl_func,"spot_x [3] : %e", spot_x[3] );
     }
     else {
       cpl_msg_warning(cpl_func,"Cannot get the default values for spot_x ");
     }

     if ( cpl_table_has_column(roof_pos_table , "spot_y") ) {
         spot_y= cpl_table_get_data_double (roof_pos_table, "spot_y");
         cpl_msg_info(cpl_func,"spot_y [0] : %e", spot_y[0] );
         cpl_msg_info(cpl_func,"spot_y [1] : %e", spot_y[1] );
         cpl_msg_info(cpl_func,"spot_y [2] : %e", spot_y[2] );
         cpl_msg_info(cpl_func,"spot_y [3] : %e", spot_y[3] );
     }
     else {
       cpl_msg_warning(cpl_func,"Cannot get the default values for spot_y ");
     }

     if ( cpl_table_has_column(roof_pos_table , "roof_pos") ) {
         roof_pos= cpl_table_get_data_double (roof_pos_table, "roof_pos");
         cpl_msg_info(cpl_func,"roof_pos [0] : %e", roof_pos[0] );
         cpl_msg_info(cpl_func,"roof_pos [1] : %e", roof_pos[1] );
         cpl_msg_info(cpl_func,"roof_pos [2] : %e", roof_pos[2] );
         cpl_msg_info(cpl_func,"roof_pos [3] : %e", roof_pos[3] );
     }
     else {
       cpl_msg_warning(cpl_func,"Cannot get the default values for roof_pos ");
     }
    /* ------------------ END EKW 26/11/2018 read constant parameter from calibration file */



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
    /* EKW 10/01/2019 cpl_size nsy = 512; */
    if ( cpl_propertylist_has (header, "ESO DET1 FRAMES NX") ) {
      nsx = cpl_propertylist_get_int (header, "ESO DET1 FRAMES NX");
      /* EKW 10/01/2019  nsy = cpl_propertylist_get_int (header, "ESO DET1 FRAMES NY"); */
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
        else {
            cpl_msg_info (cpl_func, "%s not in header, use %f", name, rp);
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
        char telid[128];
        if (telname[0] == 'U') {
            /* scale = 18.; */
            scale = cpl_propertylist_get_double (gravi_data_get_plist(static_param_data,GRAVI_PRIMARY_HDR_EXT), "ESO PLATE SCALE UT");
            cpl_msg_info (cpl_func,"PLATE SCALE UT is  : %e",scale);
        } else if (telname[0] == 'A') {
            sprintf(telid, "ESO PLATE SCALE AT%d", tel+1);
            scale = cpl_propertylist_get_double (gravi_data_get_plist(static_param_data,GRAVI_PRIMARY_HDR_EXT), telid);
            cpl_msg_info (cpl_func,"PLATE SCALE AT%d is : %e",tel+1, scale);
            /*
            if (tel == 0) {
                scale = 76.8;
            } else if (tel == 1) {
                scale = 78.0;
            } else if (tel == 2) {
                scale = 77.0;
            } else if (tel == 3) {
                scale = 84.6;
            }
            */
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
        double northangle = gravi_pfits_get_northangle_acqcam (header, tel);
        CPLCHECK ("Cannot get rotation");
        
        /* Mapping/mosaicing offset on acq cam axes, in mas, */
        /* neglecting amnamorphism variations */
        double sobj_offx_cam = sobj_drho * sin((northangle+sobj_dth)/180.*M_PI);
        double sobj_offy_cam = sobj_drho * cos((northangle+sobj_dth)/180.*M_PI);
        
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
        
        if (axis_mode == MODE_ONAXIS) {
            /* TODO: close dual-field */
            /* Single-field case */
            /* Simply shift the best spot from full frame to cut-out */
            xFT = spot_x[tel] - sx + nsx*tel + 1;
            yFT = spot_y[tel] - sy + 1;
            xSC = xFT;
            ySC = yFT;
        } else if (sts_mode == DUAL_STS) {
            xFT = fiber_xft - sx + nsx*tel + 1;
            yFT = fiber_yft - sy + 1;
            xSC = fiber_xsc - sx + nsx*tel + 1;
            ySC = fiber_ysc - sy + 1;
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
        
        if (axis_mode != MODE_ONAXIS) {
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
        
        if (axis_mode != MODE_ONAXIS) {
            sprintf (qc_name, "ESO QC ACQ FIELD%i SCALE", tel+1);
            double sep = sqrt((ySC-yFT)*(ySC-yFT)+(xSC-xFT)*(xSC-xFT));
	    /* FE 2019-10-31 in case separation could not  be measured, i.e. 0, then use default scale */
	    double pscale = scale;
	    if (sts_mode != DUAL_STS) pscale = sep ? rho_in/sep : scale;
	    /*            double pscale = sep ? rho_in/sep : 0.; */
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
        
        if (axis_mode != MODE_ONAXIS) {
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
	
        /*
         * JIRA PIPE-9123 : Pipeline faile with Strehl=NaN
         * For 2017 AT data the spots jump up by 30 pixel and therefore are not
         * catched in the subwindow around roof_pos , cut_out +/- 50 pixel
        */
        if (cpl_error_get_code() != CPL_ERROR_NONE  ) {
           cpl_msg_info (cpl_func, "WARNING Filling STREHL due to NaN");
           cpl_error_reset();
           strehl_on_average = 0.0;
        }
	
        gravi_acq_measure_max(mean_img, xFT, yFT, 15, &max_on_average);
        
        /* Update Strehl QC */
        sprintf (qc_name, "ESO QC ACQ FIELD%i STREHL", tel+1);
        cpl_msg_info (cpl_func, "%s = %f", qc_name, strehl_on_average);
        cpl_propertylist_update_double (o_header, qc_name, strehl_on_average);
        cpl_propertylist_set_comment (o_header, qc_name, "Average Strehl from stacked AcqCam images");


        /* Adding Strehl for SC channel in wide mode (PIPE-9913) */
        if (sts_mode == DUAL_STS) {
            gravi_acq_measure_strehl(mean_img, xSC, ySC, scale, &strehl_on_average, header);
	
            /*
             * JIRA PIPE-9123 : Pipeline faile with Strehl=NaN
             * For 2017 AT data the spots jump up by 30 pixel and therefore are not
             * catched in the subwindow around roof_pos , cut_out +/- 50 pixel
            */
            if (cpl_error_get_code() != CPL_ERROR_NONE  ) {
               cpl_msg_info (cpl_func, "WARNING Filling STREHL SC due to NaN");
               cpl_error_reset();
               strehl_on_average = 0.0;
            }
        
            /* Update Strehl QC */
            sprintf (qc_name, "ESO QC ACQ FIELD%i STREHLSC", tel+1);
            cpl_msg_info (cpl_func, "%s = %f", qc_name, strehl_on_average);
            cpl_propertylist_update_double (o_header, qc_name, strehl_on_average);
            cpl_propertylist_set_comment (o_header, qc_name, "Average Strehl from stacked AcqCam images");
        }
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
            if (axis_mode != MODE_ONAXIS) {
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
            double pscale = scale;
	    if (sts_mode != DUAL_STS) pscale = sep ? rho_in/sep : 0.;
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

        double strehl_std = gravi_table_get_column_std (acqcam_table, "FIELD_STREHL", tel, ntel);
        sprintf (qc_name, "ESO QC ACQ FIELD%i STREHL STD", tel+1);
        cpl_propertylist_update_double (o_header, qc_name, strehl_std);
        cpl_propertylist_set_comment (o_header, qc_name, "Std of FT strehl");
        
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
                                    gravi_data * input_data,
                                    gravi_data * static_param_data)
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
    
    /* Check if ISS data present (to avoid crash) */
    if (!gravi_conf_get_telname (0, header)) {
        gravi_msg_warning (cpl_func, "Cannot reduce the ACQCAM, no ISS keywords");
        return CPL_ERROR_NONE;
    }
    
    cpl_imagelist * acqcam_imglist;
    cpl_imagelist * acqcam_imglist_v2;
    acqcam_imglist = gravi_data_get_cube (input_data, GRAVI_IMAGING_DATA_ACQ_EXT);
    acqcam_imglist_v2 = gravi_data_get_cube (input_data, GRAVI_IMAGING_DATA_ACQ_EXT);
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
                        acqcam_table, o_header, static_param_data);

	CPLCHECK_MSG ("Cannot reduce field images");
    
    /* Compute PUPIL columns with algorithm 2.0*/
    gravi_acqcam_pupil_v2 (mean_img, acqcam_imglist_v2, header,
                        acqcam_table, o_header, static_param_data);
    CPLCHECK_MSG ("Cannot reduce pupil images");
    
    /* Add this output table in the gravi_data */
    cpl_propertylist * plist_acq_cam = cpl_propertylist_new ();
    cpl_propertylist_update_string (plist_acq_cam, "INSNAME", INSNAME_ACQ);
	gravi_data_add_img (output_data, plist_acq_cam, GRAVI_IMAGING_DATA_ACQ_EXT, mean_img);
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

double gravi_acqcam_z2meter (double PositionPixels, gravi_data * static_param_data)
{
    double f_PT      ; /*  = 14e-3;        */ /* pupil tracker lenslet FL*/
    double f_lens    ; /*  = 467e-3;       */ /* folding optics lens FL */
    double Llambda   ; /*  = 1.2e-6;       */  /* laser diode wavelength */
    double D_beam    ; /*  = 18e-3;        */  /* meter */
    double D_pixel   ; /*  = 18e-6;        */
    double D_AT      ; /*  = 1.8;          */  /* m */
    double D_lenslet ; /*  = 2 * 1.015e-3; */

    f_PT = cpl_propertylist_get_double (gravi_data_get_plist(static_param_data,GRAVI_PRIMARY_HDR_EXT), "ESO LENS fPT");
    cpl_msg_debug (cpl_func,"ESO LENS fPT is   : %e", f_PT);

    f_lens = cpl_propertylist_get_double (gravi_data_get_plist(static_param_data,GRAVI_PRIMARY_HDR_EXT), "ESO LENS flens");
    cpl_msg_debug (cpl_func,"ESO LENS f_lens is  : %e", f_lens);

    Llambda = cpl_propertylist_get_double (gravi_data_get_plist(static_param_data,GRAVI_PRIMARY_HDR_EXT), "ESO LENS Llambda");
    cpl_msg_debug (cpl_func,"ESO LENS Llambda is  : %e", Llambda);

    D_beam = cpl_propertylist_get_double (gravi_data_get_plist(static_param_data,GRAVI_PRIMARY_HDR_EXT), "ESO LENS Dbeam");
    cpl_msg_debug (cpl_func,"ESO LENS D_beam is  : %e", D_beam);

    D_pixel = cpl_propertylist_get_double (gravi_data_get_plist(static_param_data,GRAVI_PRIMARY_HDR_EXT), "ESO LENS Dpixel");
    cpl_msg_debug (cpl_func,"ESO LENS D_pixel is  : %e", D_pixel);

    D_AT = cpl_propertylist_get_double (gravi_data_get_plist(static_param_data,GRAVI_PRIMARY_HDR_EXT), "ESO LENS DAT");
    cpl_msg_debug (cpl_func,"ESO LENS D_AT is  : %e", D_AT);

    D_lenslet = cpl_propertylist_get_double (gravi_data_get_plist(static_param_data,GRAVI_PRIMARY_HDR_EXT), "ESO LENS Dlenslet");
    cpl_msg_debug (cpl_func,"ESO LENS D_lenslet is  : %e", D_lenslet);


    double longDef;
    longDef = 8 * (f_PT / D_lenslet) * (f_PT / D_lenslet) * 3.5 * D_pixel *
        D_beam / (f_PT * D_lenslet) * Llambda / CPL_MATH_2PI * PositionPixels;
    
    return f_lens * f_lens * longDef / (f_PT + longDef) / f_PT * (D_AT / D_lenslet);
}


/*----------------------------------------------------------------------------*/
/**
 * @brief Correlate two images using FFT.
 *  
 * @param ia Input image cut down to one telescope
 * @param ib Input Model
 * @param xd Output x shift
 * @param yd Output y shift
 *
 * Input images are not destroyed
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_image_fft_correlate (cpl_image *ia, cpl_image *ib, cpl_size *xd, cpl_size *yd)
{
  gravi_msg_function_start(0);                                                                                                    
  cpl_ensure_code (ia, CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (ib, CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (xd, CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (yd, CPL_ERROR_NULL_INPUT);
  
  cpl_error_code error = CPL_ERROR_NONE;
  cpl_type type  = cpl_image_get_type(ia);
  cpl_size nx    = cpl_image_get_size_x(ia);
  cpl_size ny    = cpl_image_get_size_y(ia);

  cpl_image * ic = cpl_image_new(nx, ny, type);
  cpl_image * fa = cpl_image_new(nx, ny, type | CPL_TYPE_COMPLEX);
  cpl_image * fb = cpl_image_new(nx, ny, type | CPL_TYPE_COMPLEX);
  cpl_image * fc = cpl_image_new(nx, ny, type | CPL_TYPE_COMPLEX);

  cpl_imagelist * iab = cpl_imagelist_new();
  cpl_imagelist * fab = cpl_imagelist_new();

  /* Put image ia into imagelist iab */
  error = cpl_imagelist_set(iab, ia, 0);

  /* Put image ib into imagelist iab */
  error = cpl_imagelist_set(iab, ib, 1);

  /* Put empty Fourier image fa in imagelist fab */
  error = cpl_imagelist_set(fab, fa, 0);

  /* Put empty Fourier image fb in imagelist fab*/
  error = cpl_imagelist_set(fab, fb, 1);

  /* FFT on the image list iab and fill FFT image list fab */
  error = cpl_fft_imagelist(fab, iab, CPL_FFT_FORWARD);

  /* conjugate fourier image list fb to fc */
  error = cpl_image_conjugate(fc, fb);

  /* Multiply fourier image list fa and fc*/
  error = cpl_image_multiply(fc, fa);

  /* Reverse FFT */
  error = cpl_fft_image(ic, fc, CPL_FFT_BACKWARD | CPL_FFT_NOSCALE);

  /* Find max position -> shift*/
  error = cpl_image_get_maxpos(ic, xd, yd);

  /* Unwrap for "negative" maximum position */
  if (*xd > nx/2) *xd = *xd - nx;
  if (*yd > ny/2) *yd = *yd - ny;

  /* subtract one for come from fits convention (lower left pixel = 1,1) to shift vector */
  (*xd)--;
  (*yd)--;

  //  cpl_image_save (ia, "ia.fits", CPL_TYPE_DOUBLE,NULL, CPL_IO_DEFAULT);
  //  cpl_image_save (ib, "ib.fits", CPL_TYPE_DOUBLE,NULL, CPL_IO_DEFAULT);
  //  cpl_image_save (ic, "ic.fits", CPL_TYPE_DOUBLE,NULL, CPL_IO_DEFAULT);

  /* Free memory allocated by this routine */
  FREE (cpl_imagelist_unwrap, iab);
  FREE (cpl_imagelist_unwrap, fab);
  FREE (cpl_image_delete, ic);
  FREE (cpl_image_delete, fa);
  FREE (cpl_image_delete, fb);
  FREE (cpl_image_delete, fc);

  gravi_msg_function_exit(0);
  return (error);
}

/**@}*/
