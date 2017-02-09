
/*********************************************************************
 * E.S.O. - VLT project
 *
 * "@(#) $Id: clipm_priv_checks.h 177213 2008-12-05 15:33:51Z hlorch $"
 *
 * Private functions for checking parameters
 *
 * who       when        what
 * --------  ----------  ----------------------------------------------
 * hlorch    2006-08-12  created
 */
 
#ifndef CLIPM_PRIV_CHECKS_H
#define CLIPM_PRIV_CHECKS_H

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

int             clipm_priv_checks_is_window_full_image(
                                            const int       window_xxyy[4],
                                            const cpl_image *image);
cpl_error_code  clipm_priv_checks_window_image(
                                            const int       window_xxyy[4],
                                            const cpl_image *image,
                                            int             allow_window_NULL,
                                            int             *img_size_xy,
                                            int             *window_size_xy,
                                            int             *buffer_start_xy);
cpl_error_code  clipm_priv_checks_window_minmax(
                                            const int       window_xxyy[4],
                                            int             ndims,
                                            int             allow_window_NULL);
cpl_error_code  clipm_priv_checks_window_guarantee(
                                            int             window_xxyy[4],
                                            int             xsize,
                                            int             ysize,
                                            int             min_wsize);
cpl_error_code  clipm_priv_checks_window_guarantee_image(
                                            int             window_xxyy[4],
                                            const cpl_image *image,
                                            int             min_wsize);
cpl_error_code  clipm_priv_checks_window_guarantee_window(
                                            int             window_xxyy[4],
                                            const int       *ref_window,
                                            int             min_wsize);
cpl_error_code  clipm_priv_checks_images_match(
                                            const cpl_image *img1,
                                            const cpl_image *img2,
                                            cpl_type        *type,
                                            int             *xsize,
                                            int             *ysize);
cpl_error_code  clipm_priv_checks_imtype_any(
                                            const cpl_image *image,
                                            cpl_type        *out_type);
cpl_error_code  clipm_priv_checks_imtype_float(
                                            const cpl_image *image,
                                            cpl_type        *out_type);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}   /* extern "C" */
#endif

#endif /* CLIPM_PRIV_CHECKS_H */
