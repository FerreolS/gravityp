
/*********************************************************************
 * E.S.O. - VLT project
 *
 * "@(#) $Id: clipm_centroiding.h 177213 2008-12-05 15:33:51Z hlorch $"
 *
 * Functions for centroiding inside image windows
 *
 * who       when        what
 * --------  ----------  ----------------------------------------------
 * hlorch    2006-08-18  created
 */

#ifndef CLIPM_CENTROIDING_H
#define CLIPM_CENTROIDING_H

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

cpl_error_code  clipm_centroiding_gauss(    const cpl_image *image,
                                            const int       window_xxyy[4],
                                            double          *xy_centre,
                                            double          *xy_centre_err,
                                            double          *xy_sigma,
                                            double          *xy_sigma_err,
                                            double          *xy_fwhm,
                                            double          *xy_fwhm_err,
                                            double          *centre_intensity,
                                            int             robustness);

cpl_error_code  clipm_centroiding_moment(   const cpl_image *image,
                                            const int       window_xxyy[4],
                                            int             allow_wdw_enlarge,
                                            double          gain,
                                            double  *out_xy_centre,
                                            double  *out_xy_centre_err,
                                            double  *out_xy_sigma,
                                            double  *out_xy_fwhm,
                                            double  *centre_intensity);

cpl_error_code  clipm_centroiding_multi_gauss(
                                            const cpl_image     *image,
                                            const cpl_matrix    *locations,
                                            unsigned int         areasize,
                                            cpl_matrix  **xy_centre,
                                            cpl_matrix  **xy_centre_err,
                                            cpl_matrix  **xy_sigma,
                                            cpl_matrix  **xy_sigma_err,
                                            cpl_matrix  **xy_fwhm,
                                            cpl_matrix  **xy_fwhm_err,
                                            cpl_matrix  **centre_intensities,
                                            cpl_array   **all_error_codes,
                                            int         robustness);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}   /* extern "C" */
#endif

#endif  /* CLIPM_CENTROIDING_H */
