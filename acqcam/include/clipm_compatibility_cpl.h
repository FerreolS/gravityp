
/*********************************************************************
 * E.S.O. - VLT project
 *
 * "@(#) $Id: clipm_compatibility_cpl.h 166626 2008-04-24 15:49:12Z hlorch $"
 *
 * Define the compatible CPL version
 *
 * who       when        what
 * --------  ----------  ----------------------------------------------
 * hlorch    2007-08-16  created
 */

/**
 * @defgroup clipm_compatibility_cpl    Get Installed CPL Version
 * @ingroup cpl_compat_macros
 * 
 * Below, the compatible CPL version is defined. It influences compilation.
 * 
 * @par Synopsis:
 * @code
 *   #include "clipm_compatibility_cpl.h"
 * @endcode
 */
/** @{ */

#ifndef CLIPM_COMPATIBILITY_CPL_H
#define CLIPM_COMPATIBILITY_CPL_H

/*-----------------------------------------------------------------------------
    Includes
 -----------------------------------------------------------------------------*/

#include <cpl.h>

/*-----------------------------------------------------------------------------
    Defines
 -----------------------------------------------------------------------------*/

#ifndef CPL_VERSION_CODE
    #error "CPL_VERSION_CODE not defined in cpl_version.h"
#endif

/**
 * @brief   Determine the major CPL version during compile time.
 * @hideinitializer
 */
#define CLIPM_GET_INSTALLED_CPL_VERSION   CPL_VERSION_CODE/65536

#if CLIPM_GET_INSTALLED_CPL_VERSION < 3
    #error "Incompatible CPL version found (too old)."
#endif

/*#if (CPL_VERSION_CODE/65536) >= 4
    #define CLIPM_GET_INSTALLED_CPL_VERSION    4
    *#warning "Compiling for CPL version 4."*
#elif (CPL_VERSION_CODE/65536) >= 3
    #define CLIPM_GET_INSTALLED_CPL_VERSION    3
    *#warning "Compiling for CPL version 3."*
#else
    #error "Incompatible CPL version found (too old)."
#endif
*/

/*----------------------------------------------------------------------------*/
/** @} */
#endif /* CLIPM_COMPATIBILITY_CPL_H */
