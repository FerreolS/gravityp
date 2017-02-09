/*********************************************************************
 * E.S.O. - VLT project
 *
 * "@(#) $Id: clipm_centroiding.c 218576 2011-08-23 18:25:12Z cgarcia $"
 *
 * Functions for centroiding inside image windows
 *
 * who       when        what
 * --------  ----------  ----------------------------------------------
 * hlorch    2006-08-18  created
 */

/**
 * @defgroup clipm_centroiding Peak Centroiding
 * @ingroup image_analysis
 *
 * This module provides functions for centroiding inside image windows.
 *
 * @par Synopsis:
 * @code
#include "clipm_centroiding.h"
 * @endcode
 */
/** @{ */
/*-----------------------------------------------------------------------------
    Includes
 -----------------------------------------------------------------------------*/

#include <math.h>
#include <string.h>

#include "clipm_centroiding.h"
#include "clipm_math.h"
#include "clipm_image_signal.h"
#include "clipm_priv_image.h"
#include "clipm_priv_image_signal.h"
#include "clipm_compatibility_replacements.h"
#include "clipm_priv_array.h"
#include "clipm_priv_error.h"
#include "clipm_priv_checks.h"
#include "clipm_priv_math.h"
#include "clipm_priv_matrix.h"
#include "clipm_priv_optimize.h"
#include "clipm_priv_vector.h"

/*-----------------------------------------------------------------------------
    Defines
 -----------------------------------------------------------------------------*/

/******************************************************************************/
/* Doxygen: exclude the following macro(s) due to bug with @internal keyword */
/** @} */

/**
 * @internal
 * @brief   Minimum window size for centroiding.
 * 
 * This constant defines the minimum required windowsize, respecting the
 * constraint of the function cpl_vector_fit_gaussian().
 */
#define CLIPM_CENTROIDING_MIN_WINDOWSIZE    5
/** @addtogroup clipm_centroiding */
/** @{ */
/******************************************************************************/

/*-----------------------------------------------------------------------------
    Private Prototypes
 -----------------------------------------------------------------------------*/

static
cpl_error_code  _clipm_centroiding_fit_window_to_range(
                                            int             *window_to_fit,
                                            const cpl_image *ref_image,
                                            const int       *ref_window,
                                            int             allow_exceed_refwdw,
                                            int             min_windowsize);
static
int             _clipm_centroiding_window_compare_size(
                                            const int       *wdw_A,
                                            const int       *wdw_B);
static
cpl_error_code  _clipm_centroiding_array_extract_common_valid(
                                            const cpl_array *a1,
                                            const cpl_array *a2,
                                            cpl_vector      **out_v1,
                                            cpl_vector      **out_v2,
                                            cpl_vector      **out_pos);

/*-----------------------------------------------------------------------------
    Private Implementation
 -----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
 * @brief   Guarantee that a window fits into an image range or into another
 *          window, cut it if necessary.
 * @param   window_to_fit       Coordinate buffer of the form
 *                              {xa, xb, ya, yb},
 *                              minimum/maximum order is irrelevant
 * @param   ref_image           Reference image
 * @param   ref_window          (Optional) reference window coordinate buffer,
 *                              can be NULL
 * @param   allow_exceed_refwdw Flag to allow @a window_to_fit to exceed
 *                              @a ref_window, in other words to ignore
 *                              @a ref_window
 * @param   min_windowsize      Minimum window size, @a window_to_fit is
 *                              enlarged accordingly if necessary
 * @return  CPL error code
 * 
 * @par Error Handling:
 * The following errors can occur:
 * - CPL_ERROR_NULL_INPUT: @a ref_image is NULL
 * - CPL_ERROR_ACCESS_OUT_OF_RANGE:
 *   - @a ref_window coordinates are outside the image range
 *   - @a min_windowsize exceeds the size of @a ref_image
 *   - @a min_windowsize exceeds the size of @a ref_window (if not NULL)
 * - CPL_ERROR_ILLEGAL_INPUT: @a min_windowsize < 1
 */
/*----------------------------------------------------------------------------*/
static
cpl_error_code  _clipm_centroiding_fit_window_to_range(
                                            int             *window_to_fit,
                                            const cpl_image *ref_image,
                                            const int       *ref_window,
                                            int             allow_exceed_refwdw,
                                            int             min_windowsize)
{
    CLIPM_TRY
    {
        clipm_priv_checks_window_image(     ref_window,
                                            ref_image,
                                            1,
                                            NULL,
                                            NULL,
                                            NULL);
        CLIPM_TRY_CHECK_ERROR_STATE();
        
        if (ref_window == NULL || allow_exceed_refwdw != 0)
            clipm_priv_checks_window_guarantee_image(
                                            window_to_fit,
                                            ref_image,
                                            min_windowsize);
        else
            clipm_priv_checks_window_guarantee_window(
                                            window_to_fit,
                                            ref_window,
                                            min_windowsize);
    }
    CLIPM_CATCH
    {
    }
    
    return CLIPM_ERROR_GET_NEW_SINCE_TRY();
}

/*----------------------------------------------------------------------------*/
/**
 * @brief   Compare the sizes of two window areas.
 * @param   wdw_A   Coordinate buffer of the form
 *                  {xa, xb, ya, yb}, must NOT be NULL,
 *                  minimum/maximum order is irrelevant
 * @param   wdw_B   Coordinate buffer of the form
 *                  {xa, xb, ya, yb}, must NOT be NULL,
 *                  minimum/maximum order is irrelevant
 * @return  1 if (A > B), 0 if (A == B), -1 if (A < B)
 * 
 * Provision of NULL pointers will crash the application.
 */
/*----------------------------------------------------------------------------*/
static
int             _clipm_centroiding_window_compare_size(
                                            const int       *wdw_A,
                                            const int       *wdw_B)
{
    int size_A[2],
        size_B[2];
    int dim;
    int A,
        B;
    for (dim = 0; dim < 2; dim++)
    {
        size_A[dim] = abs(wdw_A[2*dim] - wdw_A[2*dim+1]) + 1;
        size_B[dim] = abs(wdw_B[2*dim] - wdw_B[2*dim+1]) + 1;
    }
    A = size_A[0]*size_A[1];
    B = size_B[0]*size_B[1];
    if (A > B)
        return 1;
    else if (A == B)
        return 0;
    else
        return -1;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief   Extract array entries that are valid in two arrays.
 * @param   a1      Input array 1
 * @param   a2      Input array 2
 * @param   out_v1  (Output) vector 1
 * @param   out_v2  (Output) vector 2
 * @param   out_pos (Output) index positions
 * @return  CPL error code
 * 
 * 
 * 
 * 
 * 
 */
/*----------------------------------------------------------------------------*/
static
cpl_error_code  _clipm_centroiding_array_extract_common_valid(
                                            const cpl_array *a1,
                                            const cpl_array *a2,
                                            cpl_vector      **out_v1,
                                            cpl_vector      **out_v2,
                                            cpl_vector      **out_pos)
{
    CLIPM_TRY
    {
        int         size,
                    outsize;
        cpl_type    type1,
                    type2;
        
        /* init output */
        clipm_priv_vector_null(out_v1);
        clipm_priv_vector_null(out_v2);
        clipm_priv_vector_null(out_pos);
        
        /* check input */
        CLIPM_TRY_CHECK_AUTOMSG(            a1 != NULL, CPL_ERROR_NULL_INPUT);
        CLIPM_TRY_CHECK_AUTOMSG(            a2 != NULL, CPL_ERROR_NULL_INPUT);
        CLIPM_TRY_CHECK_AUTOMSG(            out_v1 != NULL,
                                            CPL_ERROR_NULL_INPUT);
        CLIPM_TRY_CHECK_AUTOMSG(            out_v2 != NULL,
                                            CPL_ERROR_NULL_INPUT);
        CLIPM_TRY_CHECK_AUTOMSG(            out_pos != NULL,
                                            CPL_ERROR_NULL_INPUT);
        
        size = cpl_array_get_size(a1);
        CLIPM_TRY_CHECK_AUTOMSG(            size == cpl_array_get_size(a2),
                                            CPL_ERROR_INCOMPATIBLE_INPUT);
        CLIPM_TRY_CHECK_AUTOMSG(            size > 0,
                                            CPL_ERROR_DATA_NOT_FOUND);
        
        type1 = cpl_array_get_type(a1);
        CLIPM_TRY_CHECK_AUTOMSG(            type1 == CPL_TYPE_INT
                                            || type1 == CPL_TYPE_FLOAT
                                            || type1 == CPL_TYPE_DOUBLE,
                                            CPL_ERROR_INVALID_TYPE);
        type2 = cpl_array_get_type(a2);
        CLIPM_TRY_CHECK_AUTOMSG(            type2 == CPL_TYPE_INT
                                            || type2 == CPL_TYPE_FLOAT
                                            || type2 == CPL_TYPE_DOUBLE,
                                            CPL_ERROR_INVALID_TYPE);
        
        /* determine output size */
        if ((! cpl_array_has_invalid(a1))
            && (! cpl_array_has_invalid(a2)))
        {
            outsize = size;
        }
        else
        {
            int n,
                valid1,
                valid2;
            outsize = 0;
            for (n = 0; n < size; n++)
            {
                cpl_array_get(a1, n, &valid1);
                cpl_array_get(a2, n, &valid2);
                if (valid1 == 0 && valid2 == 0)
                    outsize++;
            }
            CLIPM_TRY_ASSERT_ERROR_STATE();
        }
        CLIPM_TRY_ASSERT_ERROR_STATE();
        
        CLIPM_TRY_CHECK(                    outsize > 0,
                                            CPL_ERROR_DATA_NOT_FOUND,
                                            "a1, a2",
                                            "no common valid entries");
        
        *out_v1 = cpl_vector_new(outsize);
        *out_v2 = cpl_vector_new(outsize);
        *out_pos = cpl_vector_new(outsize);
        CLIPM_TRY_ASSERT_ERROR_STATE();
        
        /* copy data */
        if (outsize == size)
        {
            int n;
            for (n = 0; n < outsize; n++)
            {
                cpl_vector_set(*out_v1, n, cpl_array_get(a1, n, NULL));
                cpl_vector_set(*out_v2, n, cpl_array_get(a2, n, NULL));
                cpl_vector_set(*out_pos, n, n);
            }
            CLIPM_TRY_ASSERT_ERROR_STATE();
        }
        else
        {
            int nin,
                nout,
                valid1,
                valid2;
            nout = 0;
            for (nin = 0; nin < size; nin++)
            {
                double  val1,
                        val2;
                val1 = cpl_array_get(a1, nin, &valid1);
                val2 = cpl_array_get(a2, nin, &valid2);
                if (valid1 == 0 && valid2 == 0)
                {
                    cpl_vector_set(*out_v1, nout, val1);
                    cpl_vector_set(*out_v2, nout, val2);
                    cpl_vector_set(*out_pos, nout, nin);
                    nout++;
                }
            }
            CLIPM_TRY_ASSERT_ERROR_STATE();
        }
    }
    CLIPM_CATCH
    {
        clipm_priv_vector_null(out_v1);
        clipm_priv_vector_null(out_v2);
        clipm_priv_vector_null(out_pos);
    }
    
    return CLIPM_ERROR_GET_NEW_SINCE_TRY();
}

/*-----------------------------------------------------------------------------
    Implementation
 -----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
 * @brief   Determine the position of an object in an image window, for x and
 *          y separately, by collapsing the image and fitting a gaussian.
 * @param   image               Input image (FITS convention)
 * @param   window_xxyy         Coordinate buffer of the form
 *                              {xa, xb, ya, yb}, can be NULL,
 *                              minimum/maximum order is irrelevant
 * @param   xy_centre           Output centre coordinate buffer of size 2
 *                              ([x, y]), can be NULL
 * @param   xy_centre_err       Output centre error buffer (1 sigma) of size 2
 *                              ([x_err, y_err]), can be NULL
 * @param   xy_sigma            Output sigma buffer of size 2
 *                              ([sx, sy]), can be NULL
 * @param   xy_sigma_err        Output sigma error buffer (1 sigma) of size 2
 *                              ([sx_err, sy_err]), can be NULL
 * @param   xy_fwhm             Output full-width-half-maximum buffer of size 2
 *                              ([fwhm_x, fwhm_y]), can be NULL
 * @param   xy_fwhm_err         Output FWHM error buffer (1 sigma) of size 2
 *                              ([fwhm_x_err, fwhm_y_err]), can be NULL
 * @param   centre_intensity    Output intensity of the (closest) centre pixel,
 *                              can be NULL
 * @param   robustness          Maximum number of retries (>= 0)
 * @return  CPL error code
 * 
 * @par Overview:
 * The central position of an enclosed object is computed by fitting a gaussian
 * to the marginal distributions respectively in both dimensions.
 * \n\n
 * The centre (in x and y) is returned in @a xy_centre, the sigma of the
 * gaussian in @a xy_sigma, and the full-width-half-maximum in @a xy_fwhm.
 * \n\n
 * The measurement errors of each of these values are returned as 1 sigma in
 * @a xy_centre_err, @a xy_sigma_err, and @a xy_fwhm_err respectively. To get
 * the error variance, these values must be squared by the user.
 * \n\n
 * If the fitting was successful, then the intensity of the pixel which is
 * closest to the center is returned in @a centre_intensity.
 * \n\n
 * @a window contains at the end the last used window coordinates, if
 * @a robustness > 0 (see below).
 * 
 * @par Robustness:
 * If @a robustness > 0 and the fitting fails in the first try, then the window
 * is resized automatically (around a centre guess, see below) and the fitting
 * is retried. This is repeated a maximum of @a robustness + 1 times.\n
 * In this case, @a window contains at the end the last used window
 * coordinates.
 * \n\n
 * The fitting is generally first repeated maximal @a robustness times until it
 * is successful. Then (if @a robustness > 0):
 * - if it was a failure, then just one other try is added,
 * - if it was successful,
 *   one final fitting process is again performed in the +- 3 sigma
 *   region in x and y around the found centroid for improving the result.
 * .
 * Reasonable values for @a robustness are:
 * - 0: for one pure gaussian fit (with exactly the given window, if provided),
 *   and no retries
 * - 1..2: for good signal/noise ratio
 * - < ca. 10: for weak signals, eventually side peaks
 * - >= ca. 10: other strange cases? None known so far.
 * .
 * If there are other peaks than the main object to be centered, it might be
 * possible that the fitting locks onto them. With enough difference in weight,
 * this should seldom happen. For the case that there are two peaks with
 * diagonal offset, the x and y fit could lock onto the different peaks. This
 * would be detected by a not successful final fitting process or strange
 * results.
 * 
 * @par Fitting Constraints:
 * - It is possible that a successful fit is returned. This does not
 *   necessarily mean that the correct centroid has been found. Therefore,
 *   always the errors (uncertainties) should be evaluated and checked,
 *   whether these are in a tolerable range<b>(!!!)</b>.
 * - If the sigma-error (sigma-uncertainty) is greater than half the sigma,
 *   then the fitting is interpreted as a failure (and evt. retried
 *   depending on @a robustness).
 * 
 * @par What It Does Not:
 * This function was designed to fit a gaussian object with a flat background
 * and limited noise.
 * Searching gaussian peaks in difficult patterns is not what it does and can.
 * 
 * @note
 * As a rule of thumb, to get low measurement errors, the @a window or @a image
 * should enclose at least the 3 sigma region of the fitted gaussian, i.e.
 * for example the "whole" airy disk of a star.
 * 
 * @par Bad Pixel Handling:
 * Bad pixel maps are supported. This means that bad pixels are omitted
 * during computation of the marginal distributions, this means during
 * averaging. If this is not desired, bad pixels should be interpolated before
 * calling this function.
 * 
 * @par Error Handling:
 * - The following error codes can be set and returned:
 *   - CPL_ERROR_NULL_INPUT if @a image is NULL
 *   - CPL_ERROR_ILLEGAL_INPUT if:
 *     - the size of @a window or @a image in any dimension is less then 5
 *     - @a robustness < 0
 *   - CPL_ERROR_ACCESS_OUT_OF_RANGE if @a window contains outlying coordinates
 *   - CPL_ERROR_INVALID_TYPE if @a image is not of type int, float, double
 *   - If the fitting fails finally, one of the following codes is returned:
 *     - CPL_ERROR_ILLEGAL_OUTPUT
 *       if the centre is outside @a window or @a image (this function is
 *       intended for centroiding on existing objects, otherwise have a look at
 *       pure fitting functions like e.g. cpl_vector_fit_gaussian())
 *     - CPL_ERROR_CONTINUE:
 *       - if obviously a local but not global minimum was found, or
 *       - if the fitting algorithm could not converge
 *     - CPL_ERROR_SINGULAR_MATRIX if the covariance matrix could not be
 *       computed
 * - If a fitting of at least one dimension fails to converge (i.e. the returned
 *   error is CPL_ERROR_CONTINUE), then the following
 *   values are returned for the respective dimension (x or y):
 *   - centre:    the median position after estimating the signal background
 *   - sigma:     the median of the absolute residuals multiplied by 1.4828
 *   - fwhm:      sigma * 2.35482 (\f$\sigma\cdot 2\cdot\sqrt{2 \cdot \ln(2)}\f$)
 *   - centre_err:        -1
 *   - sigma_err:         -1
 *   - fwhm_err:          -1
 *   - centre_intensity:  not computed
 * - In all other error cases (wrong input parameters), the output values are
 *   undefined.
 * 
 * @todo
 * - flag to allow enlargement of window
 * 
 * @par Example:
@code
#include <cpl.h>
#include <clipm.h>

cpl_error_code  myfitting(  cpl_image *image,
                            int x1,
                            int y1,
                            int x2,
                            int y2)
{ 
    double  icent;
    int     window[4] = {x1, x2, y1, y2};
    double  center[2],
            sigma[2],
            center_errorsigma[2],
            sigma_errorsigma[2];
    cpl_error_code
            error;
    const double
            max_center_error = 0.5,
            max_sigma_error = 3;
    
    error = clipm_centroiding_gauss(
                            image,
                            window,
                            center,
                            center_errorsigma,
                            sigma,
                            sigma_errorsigma,
                            NULL,
                            NULL,
                            &icent,
                            7);
    if (error != CPL_ERROR_NONE ||
        center_errorsigma[0] > max_center_error ||
        center_errorsigma[1] > max_center_error ||
        sigma_errorsigma[0] > max_sigma_error ||
        sigma_errorsigma[1] > max_sigma_error)
        printf("Failure\n");
    else {
        printf("Success\n");
        printf("Pixel %f at [%.2f, %.2f]\n",
                            icent,
                            center[0],
                            center[1]);
    }
    
    return error;
}
@endcode
 */
/*----------------------------------------------------------------------------*/
cpl_error_code  clipm_centroiding_gauss(    const cpl_image *image,
                                            const int       window_xxyy[4],
                                            double          *xy_centre,
                                            double          *xy_centre_err,
                                            double          *xy_sigma,
                                            double          *xy_sigma_err,
                                            double          *xy_fwhm,
                                            double          *xy_fwhm_err,
                                            double          *centre_intensity,
                                            int             robustness)
{
    /*cpl_vector  *collapsed_errors[2] = {NULL, NULL};*/
    cpl_vector  *collapsed_signal_copy = NULL,
                *collapsed_contrib_copy = NULL,
                *collapsed_positions = NULL;
    cpl_array   *collapsed_signal[2] = {NULL, NULL},
                *collapsed_contrib[2] = {NULL, NULL};
    int         dim;

    CLIPM_TRY
    {
        int         imsize[2], /* [x, y] respectively */
                    wdwsize[2],
                    start[2];
        int         iteration_window[2][2],
                    old_window[2][2];
        int         iteration,
                    refine = 0;
        cpl_error_code
                    fitting_error[2];
        double      centre_backup[2],
                    sigma_backup[2],
                    fwhm_backup[2],
                    centre_err_backup[2],
                    sigma_err_backup[2],
                    fwhm_err_backup[2]; /* if arg is NULL */
        
        /* we need centre and sigma and centre_err and sigma_err, so
         * redirect if NULL is provided */
        if (xy_centre == NULL)
            xy_centre = centre_backup;
        if (xy_sigma == NULL)
            xy_sigma = sigma_backup;
        if (xy_fwhm == NULL)
            xy_fwhm = fwhm_backup;
        if (xy_centre_err == NULL)
            xy_centre_err = centre_err_backup;
        if (xy_sigma_err == NULL)
            xy_sigma_err = sigma_err_backup;
        if (xy_fwhm_err == NULL)
            xy_fwhm_err = fwhm_err_backup;
        
        /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         * Preparation
         * - check parameters
         * - set up some variables
         *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
        clipm_priv_checks_window_image(     window_xxyy,
                                            image,
                                            1,
                                            imsize,
                                            wdwsize,
                                            start);
        CLIPM_TRY_CHECK_ERROR_STATE();
        
        CLIPM_TRY_CHECK_AUTOMSG(            wdwsize[0] >=
                                              CLIPM_CENTROIDING_MIN_WINDOWSIZE,
                                            CPL_ERROR_ILLEGAL_INPUT);
        CLIPM_TRY_CHECK_AUTOMSG(            wdwsize[1] >=
                                              CLIPM_CENTROIDING_MIN_WINDOWSIZE,
                                            CPL_ERROR_ILLEGAL_INPUT);
        CLIPM_TRY_CHECK_AUTOMSG(            robustness >= 0,
                                            CPL_ERROR_ILLEGAL_INPUT);

        /* init iteration window */
        for (dim = 0; dim < 2; dim++)
        {
            iteration_window[dim][0] =      start[dim] + 1;
            iteration_window[dim][1] =      start[dim] + wdwsize[dim];
        }
        
        /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         * - Collapse the image
         * - Do the fitting, and iterate in the case that a fitting fails
         *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
        iteration = 0;
        while (1)
        {
            cpl_error_code  errc;
            
            CLIPM_ERROR_RECOVER_TRYSTATE();
            
            /* backup iteration window */
            memcpy(*old_window, *iteration_window, 4*sizeof(int));
            
            /* start the show */
            clipm_priv_image_collapse(      image,
                                            *iteration_window,
                                            &collapsed_signal[0],
                                            &collapsed_signal[1],
                                            &collapsed_contrib[0],
                                            &collapsed_contrib[1]);
            CLIPM_TRY_CHECK_ERROR_STATE();
            
            /* DON'T INTERRUPT this loop, try both dimensions before
             * returning an error, only ASSERTs allowed */
            for (dim = 0; dim < 2; dim++)
            {
                double      area,
                            offset,
                            next_wdw_halfsize;
                
                /* first I used collapsed_contrib to compute the sigmas
                 * and took 1/sqrt of it, inserting 1e12 where the
                 * contribution was zero to prevent the corresponding signal
                 * value from contributing to cpl_vector_fit_gaussian().
                 * Unfortunately cpl_vector_fit_gaussian does not use the
                 * sigmas when finding first guesses, so this didn't work
                 * reliably, and now here I use a function to extract only
                 * the valid values from the array collapsed_signal. */
                errc = _clipm_centroiding_array_extract_common_valid(
                                            collapsed_signal[dim],
                                            collapsed_contrib[dim],
                                            &collapsed_signal_copy,
                                            &collapsed_contrib_copy,
                                            &collapsed_positions);
                CLIPM_TRY_ASSERT(           errc == CPL_ERROR_NONE);
                
                /* convert contribution to relative errors by taking 1/sqrt() */
                errc = cpl_vector_power(    collapsed_contrib_copy,
                                            -0.5);
                CLIPM_TRY_ASSERT(           errc == CPL_ERROR_NONE);
                
                fitting_error[dim] = clipm_priv_optimize_gaussian(
                                            collapsed_positions,
                                            NULL,
                                            collapsed_signal_copy,
                                            collapsed_contrib_copy,
                                            CPL_FIT_CENTROID |
                                                CPL_FIT_STDEV |
                                                CPL_FIT_AREA |
                                                CPL_FIT_OFFSET,
                                            0, /* no optimisation, do it here */
                                            &xy_centre[dim],
                                            &xy_centre_err[dim],
                                            &xy_sigma[dim],
                                            &xy_sigma_err[dim],
                                            &xy_fwhm[dim],
                                            &xy_fwhm_err[dim],
                                            &area,
                                            &offset,
                                            NULL, NULL, NULL, NULL, NULL);
                
                if (fitting_error[dim] == CPL_ERROR_NONE)
                {
                    /* if this was successful, restrict to 3-sigma region */
                    next_wdw_halfsize = 3 * xy_sigma[dim];
                    /* make it a failure if values are not acceptable */
                    if (xy_centre_err[dim] > wdwsize[dim])
                    {
                        /* catch the case that the wrong local fitting
                         * minimum is found, which mostly should result here
                         * in a very high uncertainty. so set an error
                         * code and continue below setting centre_err etc. */
                        fitting_error[dim] = CPL_ERROR_CONTINUE;
                        if (! CLIPM_ERROR_IS_SET())
                            CLIPM_ERROR_SET_MSG(
                                fitting_error[dim],
                                "", "stuck in a local fitting solution");
                    }
                    else if (xy_sigma_err[dim] > xy_sigma[dim]/2)
                    {
                        /* catch the case that the gauss width error sigma
                         * is higher than half the gauss sigma */
                        fitting_error[dim] = CPL_ERROR_CONTINUE;
                        if (! CLIPM_ERROR_IS_SET())
                            CLIPM_ERROR_SET_MSG(
                                fitting_error[dim],
                                "", "uncertainty of result too high");
                    }
                }
                else
                {
                    /* if this was a failure, restrict to 1-sigma region,
                     * because then sigma was determined using some median
                     * and is usually too high */
                    next_wdw_halfsize = xy_sigma[dim];
                }
                xy_centre[dim] += iteration_window[dim][0]; /* window offset */
                
                /* define new iteration window, cut to valid img range later */
                iteration_window[dim][0] =  floor(  xy_centre[dim] -
                                                    next_wdw_halfsize);
                iteration_window[dim][1] =  ceil(   xy_centre[dim] +
                                                    next_wdw_halfsize);
            }
            
            /* stop if this was the refinement step,
             * or if no iteration is wanted */
            if (    robustness == 0 ||
                    refine)
                break;
            
            errc = _clipm_centroiding_fit_window_to_range(
                                            *iteration_window,
                                            image,
                                            window_xxyy,
                                            0,
                                            CLIPM_CENTROIDING_MIN_WINDOWSIZE);
            CLIPM_TRY_ASSERT(errc == CPL_ERROR_NONE);
            
            /* stop if window didn't change */
            if (0 == memcmp(*old_window, *iteration_window, 4*sizeof(int)))
                break;
            
            
            if (refine)
                break;
            
            /* if successful for the first time, or if at iteration limit,
             * then add one refinement step,
             * this refinement step reports also if x and y did not find
             * the same centroid */
            if (    iteration >= robustness ||
                    (   fitting_error[0] == CPL_ERROR_NONE &&
                        fitting_error[1] == CPL_ERROR_NONE))
                refine = 1;
            
            iteration++;
        }
        
        /* we operated using error code variables, now get back to CPL system
         * and set the respective CPL error */
        CLIPM_TRY_CHECK_ERROR_STATE();
        
        /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         * Last checks
         *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
        
        /* check if fitted centroid position is inside image, but only after
         * fitting in both x/y, maybe a user wants the values
         * anyway so don't exit the loop above */
        for (dim = 0; dim < 2; dim++)
            CLIPM_TRY_CHECK(                xy_centre[dim] >= 1 &&
                                            xy_centre[dim] <= imsize[dim],
                                            CPL_ERROR_ILLEGAL_OUTPUT,
                                            "",
                    "the gaussian fitting returned a center outside the"
                    "regarded data range, what is not appropriately modeled");
        
        /* get intensity of the central pixel */
        if (centre_intensity != NULL)
            *centre_intensity = clipm_priv_image_get_nearest_good(
                                            image,
                                            clipm_math_round_d2i(xy_centre[0]),
                                            clipm_math_round_d2i(xy_centre[1]),
                                            NULL,
                                            NULL);
        
    }
    CLIPM_CATCH
    {
    }
    
    for (dim = 0; dim < 2; dim++)
    {
        cpl_array_delete(collapsed_signal[dim]);
        cpl_array_delete(collapsed_contrib[dim]);
    }
    
    cpl_vector_delete(collapsed_signal_copy);
    cpl_vector_delete(collapsed_contrib_copy);
    cpl_vector_delete(collapsed_positions);

    return CLIPM_ERROR_GET_NEW_SINCE_TRY();
}

/*----------------------------------------------------------------------------*/
/**
 * @brief   Determine the barycenter of an object in an image window.
 * @param   image               Input image (FITS convention)
 * @param   window_xxyy         Coordinate buffer of the form
 *                              {xa, xb, ya, yb}, can be NULL,
 *                              minimum/maximum order is irrelevant
 * @param   allow_wdw_enlarge   Flag to allow the function to increase the
 *                              window size when it detects that the window
 *                              cuts the wings of the centered object,
 *                              default should be 0
 * @param   gain                (Optional) detector gain input in unit
 *                              [1/electrons], required for computation of
 *                              @a out_xy_centre_err, ignored if
 *                              @a out_xy_centre_err is NULL
 * @param   out_xy_centre       (Output) centre coordinate buffer of size 2
 *                              ([x, y])
 * @param   out_xy_centre_err   (Optional output) centre uncertainty buffer of
 *                              size 2 ([x, y]),
 *                              can be NULL
 * @param   out_xy_sigma        (Optional output) sigma buffer of size 2
 *                              ([x, y]),
 *                              can be NULL
 * @param   out_xy_fwhm         (Optional Output) FWHM (buffer of size 2
 *                              [x, y]), measured (not derived from sigma)
 * @param   centre_intensity    (Optional output) intensity of the (closest)
 *                              centre pixel,
 *                              can be NULL
 * @return  CPL error code
 * 
 * @par Principle:
 * - The background is subtracted, and the barycentre position of an
 *   enclosed object is computed, and returned in @a out_xy_centre.
 * - If @a out_xy_centre_err is not NULL, the uncertainty of the centre position
 *   is computed. This is done by statistical analysis, and so the @a gain is
 *   required therefore.
 * - If @a out_xy_sigma is not NULL, the standard deviation is computed in a
 *   region of +/- 1.1*FWHM.
 * - If @a out_xy_fwhm is not NULL, the measured full-width-half-maximum is
 *   returned.
 * 
 * @par Bad Pixel Handling:
 * Bad pixel maps are supported. This means that bad pixels are ignored during
 * the computation. But be aware that this may lead to a shift of
 * the barycenter, so a proper way of interpolation should better be applied
 * before centroiding.
 * 
 * @par Constraints:
 * - The peak must be positive, i.e. brighter than the background.
 * - @a image must be of size 3x3 or greater.
 * 
 * @par Algorithm
 * - The following process is iterated, until the @a window is stable:
 *   - The background and its noise are computed using iterative kappa-sigma
 *     clipping, using values inside +/- 3 sigma,
 *   - the barycentre is computed for values above background + 3*sigma.
 *   - in the case of failure, the @a window is enlarged,
 *   - the FWHM is measured,
 *   - and the @a window is refined to cover the +/- 2*FWHM region.
 * - If required, the PSF sigma @a out_xy_sigma and the centre uncertainty
 *   @a out_xy_centre_err are computed in the +/- 1.1*FWHM region.
 * - If required, the intensity of the pixel closest to the centre (and not
 *   flagged as bad) is returned.
 * 
 * @par Error Handling:
 * In the case of error, all output values are set to -1.0.@n
 * The following errors can be set and returned:
 * - CPL_ERROR_NULL_INPUT: if @a image or @a xy_centre is NULL
 * - CPL_ERROR_INVALID_TYPE: @a image is not of type CPL_TYPE_INT,
 *   CPL_TYPE_FLOAT or CPL_TYPE_DOUBLE
 * - CPL_ERROR_ACCESS_OUT_OF_RANGE: @a window specifies coordinates outside
 *   the @a image plane
 * - CPL_ERROR_ILLEGAL_INPUT:
 *   - the @a image (window) size is not >= 3x3
 * - CPL_ERROR_DATA_NOT_FOUND:
 *   - less than 2 pixels are good, or
 *   - around the centre of the estimated peak are (too many) bad pixels
 * - CPL_ERROR_CONTINUE: the barycentre computation could not succeed.
 * 
 * @todo
 * - CLIPM_CENTROIDING_MIN_WINDOWSIZE
 * - replace +/- 2 fwhm by _clipm_priv_image_get_optimal_box_width() (or
 *   better, use the output from clipm_priv_image_estimate_fwhm_xy())
 */
/*----------------------------------------------------------------------------*/
cpl_error_code  clipm_centroiding_moment(   const cpl_image *image,
                                            const int       window_xxyy[4],
                                            int             allow_wdw_enlarge,
                                            double          gain,
                                            double  *out_xy_centre,
                                            double  *out_xy_centre_err,
                                            double  *out_xy_sigma,
                                            double  *out_xy_fwhm,
                                            double  *centre_intensity)
{
    CLIPM_TRY
    {
        int         dim,
                    iter,
                    threshold_too_high = 0;
        const int   max_iterations = 10;    /* arbitrary limit */
        int         wdwsize[2],             /* [x, y] */
                    start[2];
        int         iteration_window[2][2],
                    old_window[2][2];
        double      xy_centre[2],
                    xy_sigma_tmp[2],
                    xy_fwhm[2];
        double      bg_mean,
                    bg_sigma,
                    threshold = 0.0;
        cpl_error_code
                    errc;
        
        /* check the input */
        clipm_priv_checks_window_image(     window_xxyy,
                                            image,
                                            1,
                                            NULL,
                                            wdwsize,
                                            start);
        CLIPM_TRY_CHECK_ERROR_STATE();
        CLIPM_TRY_CHECK_AUTOMSG(            out_xy_centre != NULL,
                                            CPL_ERROR_NULL_INPUT);
        
        CLIPM_TRY_CHECK_AUTOMSG(            gain > 0.0 ||
                                                out_xy_centre_err == NULL,
                                            CPL_ERROR_ILLEGAL_INPUT);
        CLIPM_TRY_CHECK(                    wdwsize[0] >=3 &&
                                            wdwsize[1] >=3,
                                            CPL_ERROR_ILLEGAL_INPUT,
                                            "image (window) size",
                                            "must be >= 3x3");
        
        /* init iteration window */
        for (dim = 0; dim < 2; dim++)
        {
            iteration_window[dim][0] =      start[dim] + 1;
            iteration_window[dim][1] =      start[dim] + wdwsize[dim];
        }
        
        /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         * Process:
         * - Once, estimate background and sigma,
         * Iterate:
         * 1. compute barycentre for values above mean + 3 sigma
         *    (on error, refine background and sigma),
         * 2. estimate lateral size (fwhm),
         * 3. refine window around +/- 2*fwhm region,
         * 4. exit if window did not change
         * 
         * Generally:
         * - use a threshold >= background, so that the barycentre is
         *   always inside the region
         *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
        for (iter = 0; iter < max_iterations; iter++)
        {
            int             n,
                            repeat_fwhm;
            /*printf("window [%d,%d]-[%d,%d]\n",
                                            iteration_window[0][0],
                                            iteration_window[1][0],
                                            iteration_window[0][1],
                                            iteration_window[1][1]);*/

            if (    iter == 0
                    || _clipm_centroiding_window_compare_size(
                                            *iteration_window,
                                            *old_window) > 0
                )
            {
                errc = clipm_priv_image_get_kappa_sigma(
                                            image,
                                            *iteration_window,
                                            3.0,    /* kappa */
                                            NULL,
                                            15,     /* nr of iterations */
                                            &bg_mean,
                                            &bg_sigma,
                                            NULL,
                                            NULL,
                                            NULL);
                CLIPM_TRY_ASSERT(           errc == CPL_ERROR_NONE ||
                                            errc == CPL_ERROR_CONTINUE ||
                                            errc == CPL_ERROR_DATA_NOT_FOUND ||
                                            errc == CPL_ERROR_INVALID_TYPE);
                /* complains if nr(good pixels) < 2 */
                CLIPM_ERROR_SET_MSG_IF_CODE(CPL_ERROR_CONTINUE,
                                            "no peak could be found",
                                            "");
                CLIPM_TRY_CHECK_ERROR_STATE();
                /*printf("bg %.2f, sigma %.2f\n", bg_mean, bg_sigma);*/
                threshold = bg_mean + 3 * bg_sigma;
            }

            threshold_too_high = 0;

            for (n = 0; n <= threshold_too_high; n++) /* do once or twice */
            {
                errc = clipm_image_signal_get_barycentre(
                                            image,
                                            *iteration_window,
                                            bg_mean,
                                            threshold,
                                            xy_centre,
                                            NULL,
                                            NULL);
                CLIPM_TRY_ASSERT(           errc == CPL_ERROR_NONE ||
                                            errc == CPL_ERROR_DATA_NOT_FOUND ||
                                            errc == CPL_ERROR_DIVISION_BY_ZERO);
                if (errc == CPL_ERROR_DIVISION_BY_ZERO)
                {
                    if (! threshold_too_high && allow_wdw_enlarge != 0)
                    {
                        /* maybe the window around a star is too small, so that
                         * the background and sigma values are too high.
                         * in this case compute the barycentre of all values
                         * above background level.
                         * So:
                         * - decrease threshold (for this time)
                         * - repeat barycentre computation ONCE
                         * - hope for window size increase after fwhm estimation
                         * - in next iteration, refine background mean and sigma
                         * - since the window needs to be increased now maybe
                         *   multiple times, refine bg and sigma forever!!
                         */
                        threshold_too_high = 1;
                        threshold = bg_mean;
                        CLIPM_ERROR_RECOVER_TRYSTATE();
                    }
                    else
                    {
                        CLIPM_TRY_EXIT_WITH_ERROR_MSG(
                                            CPL_ERROR_CONTINUE,
                                            "no peak could be found",
                                            "");
                    }
                }
                CLIPM_TRY_CHECK_ERROR_STATE();
            }
        
            repeat_fwhm = 0;
            for (n = 0; n <= repeat_fwhm; n++)
            {
                errc = clipm_priv_image_estimate_fwhm_xy(
                                            image,
                                            xy_centre,
                                            bg_mean,
                                            xy_fwhm,
                                            NULL,
                                            NULL);
                /* if it didn't work due to really bad data or too small
                 * window size, try once again... */
                if (errc == CPL_ERROR_CONTINUE && !repeat_fwhm)
                {
                    CLIPM_ERROR_RECOVER_TRYSTATE();
                    bg_mean = cpl_image_get_min_window(
                                            image,
                                            iteration_window[0][0],
                                            iteration_window[1][0],
                                            iteration_window[0][1],
                                            iteration_window[1][1]);
                    repeat_fwhm = 1;
                }
            }
            /* if the barycentre is outside... (yes that's possible) */
            if (    errc == CPL_ERROR_ACCESS_OUT_OF_RANGE ||
                    errc == CPL_ERROR_CONTINUE)
                CLIPM_TRY_EXIT_WITH_ERROR_MSG(
                                            CPL_ERROR_CONTINUE,
                                            "no peak could be found",
                                            "");
            CLIPM_TRY_CHECK_ERROR_STATE();
            
            /*printf("barycentre at [%.2f, %.2f], fwhm (%.2f, %.2f)\n",
                                            xy_centre[0],
                                            xy_centre[1],
                                            xy_fwhm[0],
                                            xy_fwhm[1]);*/
            
            /* backup iteration window */
            memcpy(*old_window, *iteration_window, sizeof(int)*4);
            
            /* define new iteration window, and cut to valid image range */
            for (dim = 0; dim < 2; dim++)
            {
                iteration_window[dim][0] =  floor(  xy_centre[dim] -
                                                    2*xy_fwhm[dim]);
                iteration_window[dim][1] =  ceil(   xy_centre[dim] +
                                                    2*xy_fwhm[dim]);
            }
            _clipm_centroiding_fit_window_to_range(
                                            *iteration_window,
                                            image,
                                            window_xxyy,
                                            allow_wdw_enlarge,
                                            3);
            CLIPM_TRY_ASSERT_ERROR_STATE(); /* the above code should not fail */
            
            /* stop if window didn't change */
            if (0 == memcmp(*old_window, *iteration_window, sizeof(int)*4))
                break;
        }
        CLIPM_TRY_CHECK_ERROR_STATE();
        
        for (dim = 0; dim < 2; dim++)
        {
            out_xy_centre[dim] = xy_centre[dim];
            if (out_xy_fwhm != NULL)
                out_xy_fwhm[dim] = xy_fwhm[dim];
        }

        /* clipm_priv_image_get_psf_sigma() needs xy_sigma */
        if (    out_xy_centre_err != NULL &&
                out_xy_sigma == NULL)
            out_xy_sigma = xy_sigma_tmp;
        
        if (out_xy_sigma != NULL)
        {
            double max;
            /* define new iteration window for better numerical results */
            for (dim = 0; dim < 2; dim++)
            {
                iteration_window[dim][0] =  floor(  xy_centre[dim] -
                                                    1.1*xy_fwhm[dim]);
                iteration_window[dim][1] =  ceil(   xy_centre[dim] +
                                                    1.1*xy_fwhm[dim]);
            }
            clipm_priv_checks_window_guarantee_image(
                                            *iteration_window,
                                            image,
                                            3);     /* min size */
            CLIPM_TRY_ASSERT_ERROR_STATE(); /* the above code should not fail */
            max = cpl_image_get_max_window( image,
                                            iteration_window[0][0],
                                            iteration_window[1][0],
                                            iteration_window[0][1],
                                            iteration_window[1][1]);
            errc = clipm_priv_image_get_psf_sigma(
                                            image,
                                            *iteration_window,
                                            xy_centre,
                                            bg_mean,
                                            bg_mean,
                                            max,
                                            out_xy_sigma,
                                            gain,
                                            out_xy_centre_err);
            if (errc == CPL_ERROR_DIVISION_BY_ZERO)
                CLIPM_TRY_EXIT_WITH_ERROR_MSG(
                                            CPL_ERROR_CONTINUE,
                                            "no peak could be found",
                                            "");
            CLIPM_TRY_CHECK_ERROR_STATE();
        }
        
        if (centre_intensity != NULL)
        {
            *centre_intensity = clipm_priv_image_get_nearest_good(
                                            image,
                                            xy_centre[0],
                                            xy_centre[1],
                                            NULL,
                                            NULL);
            CLIPM_TRY_CHECK_ERROR_STATE();
        }
    }
    CLIPM_CATCH
    {
        int dim;
        for (dim = 0; dim < 2; dim++)
        {
            if (out_xy_centre != NULL)
                out_xy_centre[dim] = -1.0;
            if (out_xy_centre_err != NULL)
                out_xy_centre_err[dim] = -1.0;
            if (out_xy_sigma != NULL)
                out_xy_sigma[dim] = -1.0;
            if (out_xy_fwhm != NULL)
                out_xy_fwhm[dim] = -1.0;
        }
        if (centre_intensity != NULL)
            *centre_intensity = -1.0;
    }
    
    return CLIPM_ERROR_GET_NEW_SINCE_TRY();
}

/*----------------------------------------------------------------------------*/
/**
 * @brief   Invoke clipm_centroiding_gauss() at several locations in an image.
 * @param   image               Input image (FITS convention)
 * @param   locations           Center coordinates of sub-regions (windows)
 *                              as 2 x N matrix
 * @param   areasize            Edge-length of the square sub-regions (windows),
 *                              at least 5
 * @param   xy_centre           Output centre coordinate matrix of size 2 x N
 *                              ([x, y]), can be NULL
 * @param   xy_centre_err       Output centre error matrix (1 sigma) of size
 *                              2 x N, can be NULL
 * @param   xy_sigma            Output sigma matrix of size
 *                              2 x N, can be NULL
 * @param   xy_sigma_err        Output sigma error matrix (1 sigma) of size
 *                              2 x N, can be NULL
 * @param   xy_fwhm             Output full-width-half-maximum matrix of size
 *                              2 x N, can be NULL
 * @param   xy_fwhm_err         Output FWHM error matrix (1 sigma) of size
 *                              2 x N, can be NULL
 * @param   centre_intensities  Output intensity matrix of the (closest) centre
 *                              pixels of size 1 x N,
 *                              can be NULL
 * @param   all_error_codes     Output array of type CPL_TYPE_INT and size N,
 *                              can be NULL
 * @param   robustness          Maximum number of respective retries (>= 0)
 * @return  CPL error code
 * 
 * @par Description:
 * - The function clipm_centroiding_gauss() is invoked @a N times, for @a N
 *   different square windows of @a image. The windows are defined by
 *   their centers, given by @a locations, which must be of the form:
 *     \f[ locations = \left( \begin{array}{lllc}
       x_0 & x_1 & \cdots & x_{N-1} \\
       y_0 & y_1 & \cdots & y_{N-1}
       \end{array} \right) \f]
 * - The edge-length of the square windows is @a given by @a areasize.
 * - @a all_error_codes returns all @a N returned error codes of
 *   clipm_centroiding_gauss().
 * - @a centre_intensities is returned as an 1 x @a N matrix containing the
 *   @a N intensities of the pixels closest to the @a N detected centres.
 * - All other output parameters are returned as 2 x @a N matrices, containing
 *   in each column the corresponding values, for @a x in row 0, and for @a y
 *   in row 1; For further details, please refer to clipm_centroiding_gauss().
 * - Any output parameter can be NULL.
 * - The window corner indices are rounded.
 * 
 * @par Error Handling:
 * The error handling is divided into two sections, the return value of this
 * function, and the returned centroiding error states in @a all_error_codes.
 * -# The following error codes are set and returned by this function:
 *    - CPL_ERROR_NULL_INPUT: @a image or @a locations are NULL
 *    - CPL_ERROR_ILLEGAL_INPUT:
 *      - the number of rows in @a ref_locations is not 2
 *      - @a area_size < 5
 *      - @a robustness < 0
 *    - CPL_ERROR_INVALID_TYPE: @a image is not of type CPL_TYPE_INT,
 *      CPL_TYPE_FLOAT, or CPL_TYPE_DOUBLE
 *      @n @n
 * -# These error codes are stored in @a all_error_codes if either the
 *    horizontal or vertical centroiding failed:
 *    - 0: success
 *    - CPL_ERROR_ACCESS_OUT_OF_RANGE: the specified area exceeds the
 *      boundaries of the respective image
 *    - CPL_ERROR_ILLEGAL_OUTPUT: the determined centroid seems to be outside
 *      the window
 *    - CPL_ERROR_CONTINUE: the centroiding failed
 *    - CPL_ERROR_SINGULAR_MATRIX: the uncertainties could not be computed
 *    - -1: this value just catches potential library version conflicts
 */
/*----------------------------------------------------------------------------*/
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
                                            int         robustness)
{
    cpl_matrix  **out_xy_params[6];
    
    CLIPM_TRY
    {
        int         n,
                    loc,
                    nlocations;
        int         window[4];
        double      centre_intens;
        double      xy_centre_params[6][2];
        
        out_xy_params[0] = xy_centre;
        out_xy_params[1] = xy_centre_err;
        out_xy_params[2] = xy_sigma;
        out_xy_params[3] = xy_sigma_err;
        out_xy_params[4] = xy_fwhm;
        out_xy_params[5] = xy_fwhm_err;
        
        /* init output params */
        for (n = 0; n < 6; n++)
            clipm_priv_matrix_null(out_xy_params[n]);
        clipm_priv_matrix_null(centre_intensities);
        clipm_priv_array_null(all_error_codes);
        
        /* check input params */
        CLIPM_TRY_EXIT_IFN(
            nlocations = cpl_matrix_get_ncol(locations));
        if (image == NULL)
            CLIPM_TRY_EXIT_WITH_ERROR(      CPL_ERROR_NULL_INPUT);
        CLIPM_TRY_CHECK_AUTOMSG(            areasize >= 5,
                                            CPL_ERROR_ILLEGAL_INPUT);
        if (cpl_matrix_get_nrow(locations) != 2)
            CLIPM_TRY_EXIT_WITH_ERROR_MSG(  CPL_ERROR_ILLEGAL_INPUT,
                                            "locations",
                                            CLIPM_MSG_ERR_2ROWXY);

        /* init output params */
        for (n = 0; n < 6; n++)
            if (out_xy_params[n] != NULL)
                *out_xy_params[n] = cpl_matrix_new(2, nlocations);
        if (centre_intensities != NULL)
            *centre_intensities = cpl_matrix_new(1, nlocations); /* one row */;
        if (all_error_codes != NULL)
            *all_error_codes = cpl_array_new(nlocations, CPL_TYPE_INT);
        CLIPM_TRY_CHECK_ERROR_STATE();
        
        for (loc = 0; loc < nlocations; loc++)
        {
            cpl_error_code  local_error;
            int             dim;
            
            window[0] = clipm_math_round_d2i(
                                            cpl_matrix_get(locations, 0, loc) -
                                                (areasize-1)/2);
            window[1] = window[0] + areasize - 1;
            window[2] = clipm_math_round_d2i(
                                            cpl_matrix_get(locations, 1, loc) -
                                                (areasize-1)/2);
            window[3] = window[2] + areasize - 1;
            
            local_error = clipm_centroiding_gauss(
                                            image,
                                            window,
                                            xy_centre_params[0], /* centre */
                                            xy_centre_params[1], /* c_err */
                                            xy_centre_params[2], /* sigma */
                                            xy_centre_params[3], /* s_err */
                                            xy_centre_params[4], /* fwhm */
                                            xy_centre_params[5], /* fwhm_err */
                                            &centre_intens,
                                            robustness);
            /* if the fitting itself failed */
            if (local_error == CPL_ERROR_ILLEGAL_OUTPUT ||
                local_error == CPL_ERROR_CONTINUE ||
                local_error == CPL_ERROR_SINGULAR_MATRIX ||
                local_error == CPL_ERROR_ACCESS_OUT_OF_RANGE)
            {
                CLIPM_ERROR_RECOVER_TRYSTATE();
            }
            /* if there was a fatal error... */
            CLIPM_TRY_CHECK_ERROR_STATE();
            
            for (n = 0; n < 6; n++)
                if (out_xy_params[n] != NULL)
                    for (dim = 0; dim < 2; dim++)
                        cpl_matrix_set(     *out_xy_params[n],
                                            dim,
                                            loc,
                                            xy_centre_params[n][dim]);
            if (centre_intensities != NULL)
                cpl_matrix_set(*centre_intensities, 0, loc, centre_intens);
            
            if (all_error_codes != NULL)
            {
                /* CPL_ERROR_NONE is currently 0, but to be on the safe side for
                 * all times, ensure that the returned value in this case is
                 * really 0. For the very unlikely case in future that NONE is
                 * not 0, return any other value than 0 in local_error. */
                if (local_error == CPL_ERROR_NONE)
                    local_error = 0;
                else if (local_error == 0) /* if not CPL_ERROR_NONE set to not 0 */
                    local_error = -1;
                cpl_array_set_int(*all_error_codes, loc, local_error);
            }
            CLIPM_TRY_CHECK_ERROR_STATE();
        }
    }
    CLIPM_CATCH
    {
        int n;
        for (n = 0; n < 6; n++)
            clipm_priv_matrix_null(out_xy_params[n]);
        clipm_priv_matrix_null(centre_intensities);
        clipm_priv_array_null(all_error_codes);
    }
    
    return CLIPM_ERROR_GET_NEW_SINCE_TRY();
}


/** @} */
