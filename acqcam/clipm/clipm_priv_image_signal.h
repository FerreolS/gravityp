
/*********************************************************************
 * E.S.O. - VLT project
 *
 * "@(#) $Id: clipm_priv_image_signal.h 218576 2011-08-23 18:25:12Z cgarcia $"
 *
 * Private functions for image signal processing
 *
 * who       when        what
 * --------  ----------  ----------------------------------------------
 * hlorch    2008-08-28  created
 */
 
#ifndef CLIPM_PRIV_IMAGE_SIGNAL_H
#define CLIPM_PRIV_IMAGE_SIGNAL_H

/*-----------------------------------------------------------------------------
    Includes
 -----------------------------------------------------------------------------*/

#include <cpl.h>

/*-----------------------------------------------------------------------------
    Declaration Block
 -----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#endif

/*-----------------------------------------------------------------------------
 *     Defines
 *      -----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
 * @hideinitializer
 * @brief   Loop over a pixel buffer, only regarding good pixels.
 * @param   TYPE        C data type (e.g. int, float, double)
 * @param   IMSIZE_XY   Int buffer of size 2, containing the image size [x, y]
 * @param   WDWSIZE_XY  Int buffer of size 2, containing the window size
 * @param   DATA_WDW    Pointer to bad pixel mask data (window)
 * @param   BADP_WDW    Pointer to bad pixel mask data (window)
 * @param   YACTION     Action to perform first in each row
 * @param   XACTION     Action to always perform per pixel
 * @param   XACTION_ONLY_WITH_BPM   Action to perform only if a bad pixel mask
 *                                  is provided (for speed reasons)
 */
/*----------------------------------------------------------------------------*/
#define clipm_priv_image_LOOP_GOOD_DATA_CONST( \
                                            TYPE, \
                                            IMSIZE_XY, \
                                            WDWSIZE_XY, \
                                            DATA_WDW, \
                                            BADP_WDW, \
                                            YACTION, \
                                            XACTION, \
                                            XACTION_ONLY_WITH_BPM) \
do { \
    int         x, \
                y; \
    const TYPE  *data; \
    \
    if (BADP_WDW != NULL) \
    { \
        const cpl_binary    *badp; \
        for (   y = 0, data = DATA_WDW, badp = BADP_WDW; \
                y < WDWSIZE_XY[1]; \
                y++, data += IMSIZE_XY[0], badp += IMSIZE_XY[0]) \
        { \
            YACTION \
            for (x = 0; x < WDWSIZE_XY[0]; x++) \
                if (! badp[x]) \
                { \
                    XACTION \
                    XACTION_ONLY_WITH_BPM \
                } \
        } \
    } \
    else \
    { \
        for (   y = 0, data = DATA_WDW; \
                y < WDWSIZE_XY[1]; \
                y++, data += IMSIZE_XY[0]) \
        { \
            YACTION \
            for (x = 0; x < WDWSIZE_XY[0]; x++) \
            { \
                XACTION \
            } \
        } \
    } \
} while (0)


/*-----------------------------------------------------------------------------
    Prototypes
 -----------------------------------------------------------------------------*/

double          clipm_priv_image_get_mean_perimeter(
                                            const void  *data,
                                            cpl_type    type,
                                            double      xc,
                                            double      yc,
                                            double      r,
                                            int         xsize,
                                            int         ysize,
                                            const cpl_binary
                                                        *badp,
                                            int         *nrpix);

cpl_error_code  clipm_priv_image_estimate_fwhm_xy(
                                            const cpl_image *img,
                                            const double    *xy_peakpos,
                                            double          bg_level,
                                            double          *out_xy_fwhm,
                                            double          *out_xy_middle,
                                            double          *out_xy_edge_sigma);

cpl_error_code  clipm_priv_image_get_kappa_sigma(
                                            const cpl_image *img,
                                            const int       window_xxyy[4],
                                            double          kappa,
                                            double          initial_limits[],
                                            int             nmax_iterations,
                                            double          *out_mean,
                                            double          *out_sigma,
                                            double          *out_kappasigma,
                                            int             *out_nused,
                                            int             *out_niterations);

double          clipm_priv_image_estimate_low_kappa_sigma(
                                            const cpl_image *image,
                                            const int       window_xxyy[4],
                                            double          kappa,
                                            double          *out_sigma,
                                            int             *out_nused);

cpl_error_code  clipm_priv_image_get_psf_sigma(
                                            const cpl_image *img,
                                            const int       window_xxyy[4],
                                            const double    *centre_xy,
                                            double          bg_level,
                                            double          cut_lower,
                                            double          cut_upper,
                                            double          *sigma_xy,
                                            double          gain,
                                            double          *centre_err_xy);

cpl_error_code  clipm_priv_image_collapse(  const cpl_image *image,
                                            const int       window_xxyy[4],
                                            cpl_array       **horizontal,
                                            cpl_array       **vertical,
                                            cpl_array       **x_weight_map,
                                            cpl_array       **y_weight_map);

int             clipm_priv_image_pixel_count_below(
                                            const cpl_image *img,
                                            const int       window_xxyy[4],
                                            double          limit,
                                            int             *out_nbad);

cpl_image       *clipm_priv_image_get_sat(  const cpl_image *image,
                                            const int       window_xxyy[4],
                                            cpl_image       **contrib_sat);

cpl_image       *clipm_priv_image_conv_matrix(
                                            const cpl_image     *input,
                                            const cpl_matrix    *kernel,
                                            const int           window_xxyy[4],
                                            int                 extend_bpm,
                                            int                 int2double);

cpl_image       *clipm_priv_image_filter_lowpass(
                                            const cpl_image *image,
                                            const int       window_xxyy[4],
                                            double          sigma);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}   /* extern "C" */
#endif

#endif /* CLIPM_PRIV_IMAGE_SIGNAL_H */
