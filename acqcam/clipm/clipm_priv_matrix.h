
/*********************************************************************
 * E.S.O. - VLT project
 *
 * "@(#) $Id: clipm_priv_matrix.h 177213 2008-12-05 15:33:51Z hlorch $"
 *
 * PRIVATE functions for matrix handling
 *
 * who       when        what
 * --------  ----------  ----------------------------------------------
 * hlorch    2007-02-07  created
 */
 
#ifndef CLIPM_PRIV_MATRIX_H
#define CLIPM_PRIV_MATRIX_H

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

cpl_matrix      *clipm_priv_matrix_new_from_image_window(
                                            const cpl_image     *image,
                                            const int           window_xxyy[4],
                                            const int           vflip);

cpl_error_code  clipm_priv_matrix_round(    cpl_matrix          *matrix);

cpl_matrix      *clipm_priv_matrix_get_median_rows(
                                            const cpl_matrix    *matrix);

cpl_matrix      *clipm_priv_matrix_get_mean_rows(
                                            const cpl_matrix    *matrix);

cpl_matrix      *clipm_priv_matrix_select_cols(
                                            const cpl_matrix    *matrix,
                                            const int           *selection,
                                            int                 select_nonzero);

cpl_error_code  clipm_priv_matrix_copy_col_vector(
                                            cpl_matrix          *matrix,
                                            const cpl_vector    *in,
                                            int                 col,
                                            int                 start_row);
cpl_error_code  clipm_priv_matrix_copy_vector_col(
                                            cpl_vector          *in,
                                            const cpl_matrix    *matrix,
                                            int                 col,
                                            int                 start_row);

cpl_error_code  clipm_priv_matrix_transform_create_rot2d(
                                            double              angle,
                                            const double        *centre_xy,
                                            cpl_matrix          **transform,
                                            cpl_matrix          **shift);

cpl_error_code  clipm_priv_matrix_transform_invert(
                                            const cpl_matrix    *transf,
                                            const cpl_matrix    *shift,
                                            cpl_matrix          **inv_transf,
                                            cpl_matrix          **inv_shift);

cpl_matrix      *clipm_priv_matrix_transform_points(
                                            const cpl_matrix    *points,
                                            const cpl_matrix    *transformation,
                                            const cpl_matrix    *shift);

cpl_matrix      *clipm_priv_matrix_create_corners_rectangle(
                                            const double        centre_xy[2],
                                            double              angle,
                                            const double        size_lw[2]);

void            clipm_priv_matrix_null(     cpl_matrix          **m);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}   /* extern "C" */
#endif

#endif /* CLIPM_PRIV_MATRIX_H */
