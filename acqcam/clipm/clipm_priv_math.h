
/*********************************************************************
 * E.S.O. - VLT project
 *
 * "@(#) $Id: clipm_priv_math.h 166626 2008-04-24 15:49:12Z hlorch $"
 *
 * Inline, non-public functions for basic mathematical operations
 *
 * who       when        what
 * --------  ----------  ----------------------------------------------
 * hlorch    2006-08-21  created
 */

/**
 * @internal
 * @defgroup clipm_priv_math Math Functions
 * @ingroup internal_docs
 *
 * This header file provides private math functions.
 * 
 * @par Synopsis:
 * @code
 *   #include "clipm_priv_math.h"
 * @endcode
 */
/** @{ */

#ifndef CLIPM_PRIV_MATH_H
#define CLIPM_PRIV_MATH_H

/*-----------------------------------------------------------------------------
    Includes
 -----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------
    Function Macros
 -----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
 * @brief   Get the maximum
 * @param   a   Number a
 * @param   b   Number b
 * @return  The maximum
 */
/*----------------------------------------------------------------------------*/
#define clipm_priv_math_max(a, b) \
    ((a) > (b) ? (a) : (b))

/*----------------------------------------------------------------------------*/
/**
 * @brief   Get the minimum
 * @param   a   Number a
 * @param   b   Number b
 * @return  The minimum
 */
/*----------------------------------------------------------------------------*/
#define clipm_priv_math_min(a, b) \
    ((a) < (b) ? (a) : (b))

/*----------------------------------------------------------------------------*/
/**
 * @brief   Get the sign of a number
 * @param   in  Input value
 * @return  +1 if (in >= 0), else -1
 */
/*----------------------------------------------------------------------------*/
#define clipm_priv_math_sign(in) \
    ((in) >= 0 ? 1 : -1)

/*-----------------------------------------------------------------------------
    Declaration Block
 -----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#endif

/*----------------------------------------------------------------------------*/
/**
 * @brief    Round a value to the nearest integer.
 * @param    in  Input value
 * @return   The rounded value
 * 
 * This macro is intended to replace the non-ANSI-C function rint().
 */
/*----------------------------------------------------------------------------*/
#define clipm_priv_math_round_d2i(in) \
    ((in < 0) ? (int)(in - 0.5) : (int)(in + 0.5))

/** @} */
/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}   /* extern "C" */
#endif

#endif  /* CLIPM_PRIV_MATH_H */
