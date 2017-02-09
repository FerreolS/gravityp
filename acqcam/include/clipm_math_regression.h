
/*********************************************************************
 * E.S.O. - VLT project
 *
 * "@(#) $Id: clipm_math_regression.h 169734 2008-07-02 16:20:41Z hlorch $"
 *
 * Functions for doing regression
 *
 * who       when        what
 * --------  ----------  ----------------------------------------------
 * hlorch    2007-03-05  created
 */

#ifndef CLIPM_MATH_REGRESSION_H
#define CLIPM_MATH_REGRESSION_H

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

cpl_polynomial  *clipm_math_regression_linear_1d(
                                            const double *x,
                                            const double *y,
                                            int N,
                                            double *out_mse);
                                            
cpl_polynomial  *clipm_math_regression_linear(
                                            const cpl_matrix *Xpos,
                                            const cpl_vector *Y,
                                            double *out_mse);
                                            
cpl_matrix      *clipm_math_regression_linear_series(
                                            const cpl_matrix *Xpos,
                                            const cpl_matrix *Y,
                                            double *out_mse);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}   /* extern "C" */
#endif

#endif  /* CLIPM_MATH_REGRESSION_H */
