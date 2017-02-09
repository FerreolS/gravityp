/*********************************************************************
 * E.S.O. - VLT project
 *
 * "@(#) $Id: clipm_math.c 166884 2008-04-29 15:22:18Z hlorch $"
 *
 * Functions for basic mathematical operations
 *
 * who       when        what
 * --------  ----------  ----------------------------------------------
 * hlorch    2006-05-11  created
 */

/**
 * @defgroup clipm_math Basic Functions and Constants
 * @ingroup math_analysis
 *
 * This module provides functions for basic mathematical operations, and
 * definitions of mathematical constants. These definitions are made
 * for the case that they are not available in <math.h>.
 *
 * @par Synopsis:
 * @code
#include "clipm_math.h"
 * @endcode
 */

/*-----------------------------------------------------------------------------
    Includes
 -----------------------------------------------------------------------------*/

#include "clipm_math.h"
#include "clipm_priv_math.h"
#include "clipm_priv_error.h"
 
/*-----------------------------------------------------------------------------
    Implementation
 -----------------------------------------------------------------------------*/
/**@{*/

/*----------------------------------------------------------------------------*/
/**
 * @brief   Compute the arcus tangens for the whole angle range.
 * @param   x  Adjacent leg
 * @param   y  Opposite leg
 * @return  The angle in the range \f$[0\ldots 2\pi]\f$
 */
/*----------------------------------------------------------------------------*/
float           clipm_math_arctan_0to2pi(   double x,
                                            double y) {
    float angle;
    
    if (x > 0)
        angle = atan(y / x);
    else if (x < 0)
        angle = atan(y / x) + CPL_MATH_PI;
    else
        angle = clipm_priv_math_sign(y) * CPL_MATH_PI / 2;
    if (angle < 0)
        angle += 2 * CPL_MATH_PI;

    return angle;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief   Round a double value to the nearest integer.
 * @param   in  Input value
 * @return  The rounded value
 * 
 * This function is intended to replace the non-ANSI-C function rint().
 */
/*----------------------------------------------------------------------------*/
int             clipm_math_round_d2i(       double in) {
    double offset = in < 0 ? -0.5 : 0.5;
    return (int)(in + offset);
}

/**@}*/



