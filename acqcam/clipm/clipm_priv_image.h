
/*********************************************************************
 * E.S.O. - VLT project
 *
 * "@(#) $Id: clipm_priv_image.h 177213 2008-12-05 15:33:51Z hlorch $"
 *
 * Private functions for handling images
 *
 * who       when        what
 * --------  ----------  ----------------------------------------------
 * hlorch    2006-08-12  created
 */
 
#ifndef CLIPM_PRIV_IMAGE_H
#define CLIPM_PRIV_IMAGE_H

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

cpl_image       *clipm_priv_image_wrap(     const int   *size,
                                            cpl_type    type,
                                            void        *data);

double          clipm_priv_image_get_nearest_good(
                                            const cpl_image *image,
                                            double          x,
                                            double          y,
                                            int             *out_x,
                                            int             *out_y);

const cpl_binary    *clipm_priv_image_bpm_get_if_exist(
                                            const cpl_image *image);

cpl_error_code  clipm_priv_image_bpm_or(    cpl_image       *mod_image,
                                            const cpl_image *ref_image);

cpl_error_code  clipm_priv_image_bpm_border_bad(
                                            cpl_image   *image,
                                            int         left,
                                            int         right,
                                            int         lower,
                                            int         upper);

cpl_error_code  clipm_priv_image_bpm_reject_above(
                                            cpl_image       *image,
                                            const int       *window_xxyy,
                                            double          limit);

cpl_vector      *clipm_priv_image_copy_good_data_vector(
                                            const cpl_image *image,
                                            const int       *window_xxyy,
                                            int             *out_ngood);

cpl_error_code  clipm_priv_image_get_data_const(
                                            const cpl_image *image,
                                            const int       *window_xxyy,
                                            int             allow_window_NULL,
                                            int             *img_size_xy,
                                            int             *window_size_xy,
                                            int             *buffer_start_xy,
                                            cpl_type        *type,
                                            const void      **data_start,
                                            const cpl_binary
                                                            **badp_start);

cpl_image       *clipm_priv_image_extract_cast(
                                            const cpl_image *input,
                                            const int       *window_xxyy,
                                            cpl_type        outtype);

cpl_image       *clipm_priv_image_extract_round(
                                            const cpl_image *input,
                                            const int       *window_xxyy);

cpl_error_code  clipm_priv_image_fill(      cpl_image       *image,
                                            const int       *window_xxyy,
                                            double          value);

cpl_error_code  clipm_priv_image_fill_circle(
                                            cpl_image       *img_modified,
                                            double          centre_xy[2],
                                            double          radius,
                                            double          brightness,
                                            int             anti_alias,
                                            int             additive);

cpl_error_code  clipm_priv_image_fill_polygon(
                                            cpl_image           *image,
                                            const cpl_matrix    *points,
                                            double              brightness,
                                            int                 anti_alias,
                                            int                 additive);

cpl_error_code  clipm_priv_image_fill_rectangle(
                                            cpl_image           *image,
                                            const double        centre_xy[2],
                                            double              angle,
                                            const double        size_lw[2],
                                            double              brightness,
                                            int                 anti_alias,
                                            int                 additive);

cpl_error_code  clipm_priv_image_save_debug(const cpl_image *img,
                                            const char      *filename);

void            clipm_priv_image_null(      cpl_image       **i);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}   /* extern "C" */
#endif

#endif /* CLIPM_PRIV_IMAGE_H */
