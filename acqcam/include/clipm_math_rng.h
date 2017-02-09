
/*********************************************************************
 * E.S.O. - VLT project
 *
 * "@(#) $Id: clipm_math_rng.h 152480 2007-06-12 11:16:32Z hlorch $"
 *
 * Functions for generation of random values
 *
 * who       when        what
 * --------  ----------  ----------------------------------------------
 * hlorch    2006-08-08  created
 */

#ifndef CLIPM_MATH_RNG_H
#define CLIPM_MATH_RNG_H

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

double          clipm_math_rng_uniform50(   void);
double          clipm_math_rng_gaussian(    void);
void            *clipm_math_rng_poisson_pointpattern_2d(
                                            double  xmin,
                                            double  ymin,
                                            double  xmax,
                                            double  ymax,
                                            int     homogeneity,
                                            int     N,
                                            cpl_type type);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}   /* extern "C" */
#endif

#endif  /* CLIPM_MATH_RNG_H */
