
/*********************************************************************
 * E.S.O. - VLT project
 *
 * "@(#) $Id: clipm_aperture.h 177213 2008-12-05 15:33:51Z hlorch $"
 *
 * Functions for aperture characterisation
 *
 * who       when        what
 * --------  ----------  ----------------------------------------------
 * hlorch    2006-12-08  created
 */

#ifndef CLIPM_APERTURE_H
#define CLIPM_APERTURE_H

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

cpl_error_code  clipm_aperture_detect_circles(
                                            const cpl_image *input,
                                            const int       window_xxyy[4],
                                            cpl_matrix      **out_circles);

cpl_error_code  clipm_aperture_characterise_circular(
                                            const cpl_image *input,
                                            const int       window_xxyy[4],
                                            double          *out_centre_x,
                                            double          *out_centre_y,
                                            double          *out_radius,
                                            double          *out_r_sigma,
                                            cpl_matrix      **statistics);

cpl_error_code  clipm_aperture_characterise_rectangular(
                                            const cpl_image *input,
                                            const int       window_xxyy[4],
                                            double          *centre_x,
                                            double          *centre_y,
                                            double          *angle,
                                            double          *length,
                                            double          *width,
                                            cpl_matrix      **statistics);

/*cpl_error_code  clipm_aperture_characterise_square(
                                            const cpl_image *input,
                                            const int   window_xxyy[4],
                                            cpl_vector  **out_centre_coord,
                                            double      *out_size,
                                            double      *out_angle,
                                            cpl_mask    **out_mask,
                                            cpl_stats   **out_statistics,
                                            unsigned    stats_bitmask);*/

cpl_error_code  clipm_aperture_slitpos(     cpl_image *input,
                                            int     max_width,
                                            double  *out_x_centre,
                                            double  *out_y_centre,
                                            double  *out_angle,
                                            int     *out_y_size);

cpl_matrix      *clipm_aperture_get_rectangle_corners(
                                            double  centre_x,
                                            double  centre_y,
                                            double  angle,
                                            double  length,
                                            double  width);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}   /* extern "C" */
#endif

#endif  /* CLIPM_APERTURE_H */
