
/*********************************************************************
 * E.S.O. - VLT project
 *
 * "@(#) $Id: clipm_image_signal.h 218576 2011-08-23 18:25:12Z cgarcia $"
 *
 * Functions for image signal processing
 *
 * who       when        what
 * --------  ----------  ----------------------------------------------
 * hlorch    2008-08-28  created
 */
 
#ifndef CLIPM_IMAGE_SIGNAL_H
#define CLIPM_IMAGE_SIGNAL_H

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
    Prototypes
 -----------------------------------------------------------------------------*/

double          clipm_image_signal_estimate_bg_in_region(
                                            const cpl_image *img,
                                            const int       window_xxyy[4],
                                            double          *out_sigma,
                                            int             *out_nused);

double          clipm_image_signal_estimate_fwhm_round(
                                            const cpl_image *img,
                                            double          x_peakpos,
                                            double          y_peakpos,
                                            double          bg_level);

cpl_error_code  clipm_image_signal_get_barycentre(
                                            const cpl_image *img,
                                            const int       window_xxyy[4],
                                            double          bg_level,
                                            double          lower_cutlevel,
                                            double          *out_xy_centre,
                                            double          *out_weight,
                                            int             *out_nused);


/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}   /* extern "C" */
#endif

#endif /* CLIPM_IMAGE_SIGNAL_H */
