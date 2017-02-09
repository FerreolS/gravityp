
/*********************************************************************
 * E.S.O. - VLT project
 *
 * "@(#) $Id: clipm_math.h 166884 2008-04-29 15:22:18Z hlorch $"
 *
 * Functions for basic mathematical operations
 *
 * who       when        what
 * --------  ----------  ----------------------------------------------
 * hlorch    2006-05-11  created
 */

#ifndef CLIPM_MATH_H
#define CLIPM_MATH_H

/*-----------------------------------------------------------------------------
    Includes
 -----------------------------------------------------------------------------*/

#include <math.h>
/* ensure to have CPL mathematical constants */
#include "clipm_compatibility_replacements.h"

/*-----------------------------------------------------------------------------
    Declaration Block
 -----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#endif

/*-----------------------------------------------------------------------------
    Prototypes
 -----------------------------------------------------------------------------*/

float           clipm_math_arctan_0to2pi(   double x,
                                            double y);
int             clipm_math_round_d2i(       double in);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}   /* extern "C" */
#endif

#endif  /* CLIPM_MATH_H */
