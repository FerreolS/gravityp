
/*********************************************************************
 * E.S.O. - VLT project
 *
 * "@(#) $Id: clipm_math_correlation.h 169734 2008-07-02 16:20:41Z hlorch $"
 *
 * Functions for signal correlations
 *
 * who       when        what
 * --------  ----------  ----------------------------------------------
 * hlorch    2007-03-05  created
 */

#ifndef CLIPM_MATH_CORRELATION_H
#define CLIPM_MATH_CORRELATION_H

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
    Types
 -----------------------------------------------------------------------------*/

/**
 * @defgroup    clipm_math_correlation_modes  Coverage Modes
 * @ingroup     clipm_math_correlation
 * @brief       Coverage modes for cross-signal operations.
 * 
 * These modes define the size of the output of operations, which cross-process
 * two different input signals with each other, e.g. like a convolution or a
 * cross-correlation.
 * 
 * Details on the options are documented in the respective functions.
 */
/** @{ */

/** @brief Interpolation option type */
typedef unsigned int clipm_coverage_mode;

/** @brief  Return the portion of the cross-operation that
 *          is computed without the zero-padded edges. */
extern const clipm_coverage_mode CLIPM_COVERAGE_VALID;
/** @brief  Returns the central part of the result that is
 *          the same size as the first input signal. */
extern const clipm_coverage_mode CLIPM_COVERAGE_SAME;
/** @brief  Compute the full result. */
extern const clipm_coverage_mode CLIPM_COVERAGE_FULL;
/** @brief  Use user-defined output size. */
extern const clipm_coverage_mode CLIPM_COVERAGE_CUSTOM;

/** @} */

/*-----------------------------------------------------------------------------
    Prototypes
 -----------------------------------------------------------------------------*/

int             clipm_math_get_coverage_size_1d(
                                            int size1,
                                            int size2,
                                            int custom_xy_outsize,
                                            clipm_coverage_mode cov);

cpl_error_code  clipm_math_get_coverage_size(
                                            const int *size1,
                                            const int *size2,
                                            const int *custom_xy_outsize,
                                            int *out_size,
                                            int ndims,
                                            clipm_coverage_mode cov);

cpl_matrix      *clipm_math_xcorr_image(    const cpl_image *image1,
                                            const cpl_image *image2,
                                            const int       *window1,
                                            const int       *window2,
                                            clipm_coverage_mode cov,
                                            const int   *custom_xy_outsize,
                                            cpl_matrix  **overlap_map);

cpl_matrix      *clipm_math_xcorr_matrix(   const cpl_matrix    *m1,
                                            const cpl_matrix    *m2,
                                            clipm_coverage_mode cov,
                                            const int   *custom_xy_outsize,
                                            cpl_matrix  **overlap_map);

cpl_matrix      *clipm_math_normxcorr_image(const cpl_image *image1,
                                            const cpl_image *image2,
                                            const int       *window1,
                                            const int       *window2,
                                            clipm_coverage_mode cov,
                                            const int       *custom_xy_outsize);

cpl_matrix      *clipm_math_normxcorr_matrix(
                                            const cpl_matrix    *m1,
                                            const cpl_matrix    *m2,
                                            clipm_coverage_mode cov,
                                            const int   *custom_xy_outsize);

cpl_image       *clipm_math_conv_image_matrix(
                                            const cpl_image     *image,
                                            const cpl_matrix    *kernelmat,
                                            const int           *im_window,
                                            clipm_coverage_mode cov,
                                            const int   *custom_xy_outsize,
                                            int         norm_to_kernel,
                                            cpl_matrix  **overlap_map);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}   /* extern "C" */
#endif

#endif  /* CLIPM_MATH_CORRELATION_H */
