
/*********************************************************************
 * E.S.O. - VLT project
 *
 * "@(#) $Id: clipm_priv_array.h 176658 2008-11-26 18:02:24Z hlorch $"
 *
 * PRIVATE functions for array handling
 *
 * who       when        what
 * --------  ----------  ----------------------------------------------
 * hlorch    2007-04-27  created
 */
 
#ifndef CLIPM_PRIV_ARRAY_H
#define CLIPM_PRIV_ARRAY_H

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

cpl_array       *clipm_priv_array_new_from_image_row(
                                            const cpl_image *image,
                                            int             row_ndx);

cpl_array       *clipm_priv_array_new_from_image_col(
                                            const cpl_image *image,
                                            int             col_ndx);

double          clipm_priv_array_estimate_fwhm(
                                            const cpl_array     *input,
                                            double              centre,
                                            double              bg_level,
                                            int                 *out_maxindex,
                                            double              *out_middlepos,
                                            double              *out_edgeslope);

void            clipm_priv_array_null(      cpl_array       **a);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}   /* extern "C" */
#endif

#endif /* CLIPM_PRIV_ARRAY_H */
