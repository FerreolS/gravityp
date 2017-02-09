/*********************************************************************
 * E.S.O. - VLT project
 *
 * "@(#) $Id: clipm_priv_image_signal.c 218576 2011-08-23 18:25:12Z cgarcia $"
 *
 * Private functions for image signal processing
 *
 * who       when        what
 * --------  ----------  ----------------------------------------------
 * hlorch    2008-08-28  created
 */

/**
 * @internal
 * @defgroup clipm_priv_image_signal    Image Signal Processing
 * @ingroup internal_image
 *
 * This module provides private image signal processing functions.
 *
 * @par Synopsis:
 * @code
#include "clipm_priv_image.h"
 * @endcode
 */
/**@{*/

/*-----------------------------------------------------------------------------
    Includes
 -----------------------------------------------------------------------------*/

#include "clipm_priv_image_signal.h"

#include "clipm_math.h"

#include "clipm_priv_checks.h"
#include "clipm_compatibility_replacements.h"
#include "clipm_priv_array.h"
#include "clipm_priv_error.h"
#include "clipm_priv_image.h"
#include "clipm_priv_math.h"

#include <float.h>
#include <limits.h>


/*-----------------------------------------------------------------------------
    Private Prototypes
 -----------------------------------------------------------------------------*/

static
double          _clipm_priv_image_get_optimal_box_width(
                                            double  signal_height,
                                            double  fwhm,
                                            double  slope_at_halfmax,
                                            double  wanted_sigma_coverage,
                                            double  *out_sigma_estimate,
                                            int     *out_is_edge_sigma);

/*-----------------------------------------------------------------------------
    Private Implementation
 -----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
 * @brief   Average values in a buffer along a perimeter.
 * @param   data    Data buffer (unchecked)
 * @param   type    CPL data type
 * @param   xc      X centre
 * @param   yc      Y centre
 * @param   r       Radius
 * @param   xsize   X buffer size
 * @param   ysize   Y buffer size
 * @param   badp    (Optional) bad pixel buffer of type cpl_binary,
 *                  can be NULL
 * @param   nrpix   (Optional output) number of averaged pixels
 * @return  Integral value
 * 
 * @par Error Handling:
 * The following error cases can happen:
 * - CPL_ERROR_DATA_NOT_FOUND: no pixels were integrated, e.g. due to bad
 *   pixels or out of range
 * - CPL_ERROR_INVALID_TYPE: the type is not int, float or double
 */
/*----------------------------------------------------------------------------*/
double          clipm_priv_image_get_mean_perimeter(
                                            const void  *data,
                                            cpl_type    type,
                                            double      xc,
                                            double      yc,
                                            double      r,
                                            int         xsize,
                                            int         ysize,
                                            const cpl_binary
                                                        *badp,
                                            int         *nrpix)
{
    double  p = 0.0;
    
    CLIPM_TRY
    {
        if (nrpix != NULL)
            *nrpix = 0;
        
        /* truncating to integer shall be the same as rounding */
        xc += 0.5;
        yc += 0.5;
        
        r = fabs(r);
        
        if (r < 0.1)
            r = 0.1;
    
/* --- throw in a macro to ease life ---------------------------------------- */
#define         clipm_priv_image_get_mean_perimeter_LOOP( \
                                            TYPE, \
                                            ACTION) \
do { \
    double  a; \
    int     x, \
            y; \
    \
    for (a = 0; a < CPL_MATH_2PI; a += 1/r) \
    { \
        x = (int)(xc + r * cos(a)); \
        y = (int)(yc + r * sin(a)); \
        if (    x >= 0 && x < xsize && \
                y >= 0 && y < ysize) \
        { \
            ACTION \
        } \
    } \
} while (0)

#define         clipm_priv_image_get_mean_perimeter_BODY( \
                                            TYPE) \
do { \
    int     w = 0; \
    TYPE    Tperi = 0; \
    const TYPE \
            *Td; \
    Td = (const TYPE*)data; \
    \
    if (badp != NULL) \
    { \
        clipm_priv_image_get_mean_perimeter_LOOP( \
                                            TYPE, \
            if (! badp[y*xsize + x]) \
            { \
                Tperi += Td[y*xsize + x]; \
                w++; \
            } \
        ); \
    } \
    else \
    { \
        clipm_priv_image_get_mean_perimeter_LOOP( \
                                            TYPE, \
            Tperi += Td[y*xsize + x]; \
            w++; \
        ); \
    } \
    if (nrpix != NULL) \
        *nrpix = w; \
    if (w > 0) \
        p = (double)Tperi / (double)w; \
    else \
        CLIPM_TRY_EXIT_WITH_ERROR(CPL_ERROR_DATA_NOT_FOUND); \
} while (0)
/*----------------------------------------------------------------------------*/

        switch (type)
        {
            case CPL_TYPE_INT: {
                clipm_priv_image_get_mean_perimeter_BODY(int);
                break;
            }
            case CPL_TYPE_FLOAT: {
                clipm_priv_image_get_mean_perimeter_BODY(float);
                break;
            }
            case CPL_TYPE_DOUBLE: {
                clipm_priv_image_get_mean_perimeter_BODY(double);
                break;
            }
            default:
                CLIPM_TRY_EXIT_WITH_ERROR_MSG(
                                            CPL_ERROR_INVALID_TYPE,
                                            "data",
                                            "must be int, float or double");
        }
    }
    CLIPM_CATCH
    {
    }
    
    return p;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief   Force the iterative kappa-sigma clipping to converge against the dark.
 * @param   image       Input image
 * @param   window_xxyy Coordinate buffer of the form
 *                      {x1a, x1b, y1a, y1b},
 *                      minimum/maximum order is irrelevant
 * @param   kappa       Kappa
 * @param   out_sigma   (Optional output) RMS of the estimated background
 *                      -1.0 in the case of error
 * @param   out_nused   (Optional output) number of used pixels,
 *                      0 in the case of error
 * @return  The estimated background level
 * 
 * @par Principle:
 * - This function is intended to find the background of a (small) image
 *   region, which might contain other (brighter) objects. It is not intended
 *   to estimate the background of a big picture.
 * - The minimum, median and maximum pixel values are used as initial limits
 *   for two runs of iterative kappa-sigma clipping. While these both converge
 *   to different values, the minimum, median and maximum are now taken from
 *   the lower part and the iterations of kappa-sigma clipping are re-run,
 *   until both converge to the same values. This is then an indication that
 *   both belong to the same statistics. It requires of course, that the
 *   kappa-sigma clipping routine can leave its initial limits.
 * 
 * @par Constraints:
 * The background is expected to be darker than other objects in the image
 * region.
 * 
 * @par Bad Pixel Handling:
 * Bad pixels are left out.
 * 
 * @par Error Handling:
 * The following error codes can be set and returned:
 * - CPL_ERROR_NULL_INPUT: @a image is NULL
 * - CPL_ERROR_INVALID_TYPE: @a image is not of type CPL_TYPE_INT,
 *   CPL_TYPE_FLOAT or CPL_TYPE_DOUBLE
 * - CPL_ERROR_ACCESS_OUT_OF_RANGE: @a window_xxyy specifies coordinates outside
 *   the @a image plane
 * - CPL_ERROR_ILLEGAL_INPUT: @a kappa < 1.0
 * - CPL_ERROR_DATA_NOT_FOUND: less than 5 pixels are good
 */
/*----------------------------------------------------------------------------*/
double          clipm_priv_image_estimate_low_kappa_sigma(
                                            const cpl_image *image,
                                            const int       window_xxyy[4],
                                            double          kappa,
                                            double          *out_sigma,
                                            int             *out_nused)
{
    cpl_vector  *img_sorted = NULL;
    double      result_mean = -1.0,
                result_sigma = -1.0;
    int         nused = 0;
    
    CLIPM_TRY
    {
        int             imsize[2],
                        wdwsize[2],
                        buffstart[2];
        int             ngood,
                        npix;
        
        /* check input */
        clipm_priv_checks_window_image(     window_xxyy,
                                            image,
                                            1,
                                            imsize,
                                            wdwsize,
                                            buffstart);
        CLIPM_TRY_CHECK_ERROR_STATE();
        clipm_priv_checks_imtype_any(       image, NULL);
        CLIPM_TRY_CHECK_ERROR_STATE();
        CLIPM_TRY_CHECK_AUTOMSG(            kappa >= 1.0,
                                            CPL_ERROR_ILLEGAL_INPUT);
        
        /* create a sorted vector containing all image pixel values */
        img_sorted = clipm_priv_image_copy_good_data_vector(
                                            image,
                                            window_xxyy,
                                            &ngood);
        CLIPM_TRY_CHECK_ERROR_STATE();
        
        CLIPM_TRY_CHECK(                    ngood >= 5,
                                            CPL_ERROR_DATA_NOT_FOUND,
                                            "image",
                                            "has less than 5 good pixels");
        
        cpl_vector_sort(img_sorted, +1);
        CLIPM_TRY_ASSERT_ERROR_STATE();
        
        npix = ngood;
        while(1)
        {
            double  initial_limits[2];
            double  min,
                    max,
                    median,
                    mean1,
                    mean2,
                    sigma1,
                    sigma2;
            
            min = cpl_vector_get(img_sorted, 0);
            max = cpl_vector_get(img_sorted, npix-1);
            median = cpl_vector_get(img_sorted, (npix-1)/2);
        
            initial_limits[0] = min;
            initial_limits[1] = median;
            clipm_priv_image_get_kappa_sigma(
                                            image,
                                            window_xxyy,
                                            kappa,
                                            initial_limits,
                                            10,
                                            &mean1,
                                            &sigma1,
                                            NULL,
                                            &nused,
                                            NULL);
            CLIPM_TRY_CHECK_ERROR_STATE();
            /*cpl_msg_info(__func__, "npix: %d", npix);*/
            /*cpl_msg_info(__func__, "bg: %.2f, sigma: %.2f", mean1, sigma1);*/
            initial_limits[0] = median;
            initial_limits[1] = max;
            clipm_priv_image_get_kappa_sigma(
                                            image,
                                            window_xxyy,
                                            kappa,
                                            initial_limits,
                                            10,
                                            &mean2,
                                            &sigma2,
                                            NULL,
                                            NULL,
                                            NULL);
            CLIPM_TRY_CHECK_ERROR_STATE();
            /*cpl_msg_info(__func__, "bg: %.2f, sigma: %.2f", mean2, sigma2);*/
            
            /* if the same within 1 sigma, then treat as converged */
            if (    npix > ngood / 20 && npix > 10
                    && (fabs(mean1 - mean2) > sigma1
                        || fabs(mean1 - mean2) > sigma2))
                npix /= 2;
            else
            {
                result_mean = mean1;
                result_sigma = sigma1;
                break;
            }
        }
    }
    CLIPM_CATCH
    {
        result_mean = -1.0;
        result_sigma = -1.0;
        nused = 0;
    }
    
    cpl_vector_delete(img_sorted);
    
    if (out_sigma != NULL)
        *out_sigma = result_sigma;
    if (out_nused != NULL)
        *out_nused = nused;
    
    return result_mean;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief   Compute sigma width by the edge slope, and return coverage width.
 * @param   signal_height           Difference between maximum and background
 * @param   fwhm                    Full-width half-maximum, measured
 * @param   slope_at_halfmax        Steepness measured at half maximum
 * @param   wanted_sigma_coverage   The wanted coverage of the distribution
 *                                  wings, in sigma, recommended: 2.58 (99%)
 * @param   out_sigma_estimate      (Optional output) the estimated sigma
 * @param   out_is_edge_sigma       (Optional output) see below
 * @return  The proposed window width
 * 
 * @par Principle:
 * - First, it is tested whether the object has a gaussian shape:
 *   - If the computation of sigma using the FWHM (and assuming a Gaussian)
 *     leads to a similar value as the computation using the steepness and
 *     the signal height, then the existence of a Gaussian is guessed.
 *   - In this case, the resulting width is
 *     (2 * wanted_sigma_coverage * sigma_computed_by_steepness)
 *   - @a out_is_edge_sigma, if provided, will be set to 0.
 * - Otherwise, a shape is assumed which was generated by convolving a flat
 *   area (step function) with a Gaussian.
 *   - For this case, the sigma is accordingly computed using the steepness
 *     and the signal height.
 *   - The resulting width is
 *     (FWHM + 2 * wanted_sigma_coverage * sigma)
 *   - @a out_is_edge_sigma, if provided, will be set to 1.
 * 
 * @par Constraints:
 * It must be ensured that:
 * - @a signal_height > 0,
 * - @a fwhm > 0, and
 * - @a slope_at_halfmax != 0.
 * 
 * @par Error Handling:
 * No other errors can occur.
 * 
 * @todo
 * - Look at if (result > 2*fwhm) result = 2*fwhm;
 */
/*----------------------------------------------------------------------------*/
static
double          _clipm_priv_image_get_optimal_box_width(
                                            double  signal_height,
                                            double  fwhm,
                                            double  slope_at_halfmax,
                                            double  wanted_sigma_coverage,
                                            double  *out_sigma_estimate,
                                            int     *out_is_edge_sigma)
{
    double  sigma_by_slope = -1.0,
            result = -1.0;
    
    CLIPM_TRY
    {
        double  sigma_by_fwhm;
        
        if (out_is_edge_sigma != NULL)
            *out_is_edge_sigma = 0;
        
        CLIPM_TRY_ASSERT(                   signal_height > 0);
        CLIPM_TRY_ASSERT(                   fwhm > 0);
        CLIPM_TRY_ASSERT(                   slope_at_halfmax != 0.0);
        
        /* for the case of a gaussian:
         * sigma = height / slope * sqrt(- (log(0.5) / 2.0)) */
        sigma_by_slope =    0.58870501125773733175
                            * signal_height / fabs(slope_at_halfmax);
        sigma_by_fwhm = fwhm / CPL_MATH_FWHM_SIG;
        result = 2 * wanted_sigma_coverage * sigma_by_slope;
        
        /* if it is rather a step function convolved with a gaussian */
        if (sigma_by_fwhm > 1.5 * sigma_by_slope)
        {
            /* sigma = height / slope / sqrt(2*Pi) */
            sigma_by_slope =    signal_height / fabs(slope_at_halfmax)
                                / CPL_MATH_SQRT2PI;
            result = fwhm + 2 * wanted_sigma_coverage * sigma_by_slope;
            
            if (out_is_edge_sigma != NULL)
                *out_is_edge_sigma = 1;
        }
        
        if (result > 2*fwhm)
            result = 2*fwhm;
        
        if (out_sigma_estimate != NULL)
            *out_sigma_estimate = sigma_by_slope;
    }
    CLIPM_CATCH
    {
        result = -1.0;
        if (out_sigma_estimate != NULL)
            *out_sigma_estimate = -1.0;
    }
    
    return result;
}

/*-----------------------------------------------------------------------------
    Implementation
 -----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
 * @brief   Estimate the background in a small image region.
 * @param   img             Input image
 * @param   window_xxyy     Coordinate buffer of the form
 *                          {xa, xb, ya, yb}, can be NULL,
 *                          minimum/maximum order is irrelevant
 * @param   out_sigma       (Optional output) background sigma,
 *                          can be NULL, returns -1.0 in the case of error
 * @param   out_nused       (Optional output) number of found background pixels,
 *                          can be NULL, returns 0 in the case of error
 * @return                  Background value, -1.0 in the case of error
 * 
 * @par Introduction:
 * - This function is intended to find the background of a (small) image
 *   region, which might contain other (brighter) objects. It is not intended
 *   to estimate the background of a big picture.
 * - It solely operates on pixel statistics, i.e. it does not assume any shape
 *   of non-background objects.
 * 
 * @par Principle:
 * - The background is found initially by a strategy using iterative
 *   kappa-sigma-clipping:
 *   - the pixel value range is divided into an upper and a lower half range,
 *   - for both ranges an iterative kappa-sigma-clipping is performed, but it
 *     is allowed to leave the initial range (this will happen if the range
 *     limit is within a major gaussian value distribution),
 *   - both results are compared, if they are equal (within 1 sigma), then
 *     the result is regarded as a solution, else it means that different
 *     value distributions do exist, and this procedure is repeated with the
 *     lower half range.
 * - If 10% of the pixels inside the image region are below the found
 *   background distribution (this means if 11% are below the 99% limit which is
 *   mean - 2.58*sigma), then the background detection (using iterative
 *   kappa-sigma-clipping) is repeated for these 10%. If the result is good
 *   (judged by whether more than 5% of the pixels contribute to this solution),
 *   then this means that not the background has been found before, and
 *   the new result is accepted, otherwise the old result is restored. If
 *   the new result is accepted, then the check whether 10% are below the found
 *   background distribution is done again and this step is repeated if
 *   necessary, until this is not the case anymore.
 * 
 * @par Constraints:
 * - The background is expected to be darker than other objects in the image
 *   region.
 * - The image region must contain at least 5 non-bad pixels.
 * 
 * @par Error Handling:
 * The following error codes can be set and returned:
 * - CPL_ERROR_NULL_INPUT: @a image is NULL
 * - CPL_ERROR_INVALID_TYPE: @a image is not of type CPL_TYPE_INT,
 *   CPL_TYPE_FLOAT or CPL_TYPE_DOUBLE
 * - CPL_ERROR_ACCESS_OUT_OF_RANGE: @a window_xxyy specifies coordinates outside
 *   the @a image plane
 * - CPL_ERROR_DATA_NOT_FOUND: less than 5 pixels are good
 */
/*----------------------------------------------------------------------------*/
double          clipm_priv_image_estimate_bg_in_region(
                                            const cpl_image *img,
                                            const int       window_xxyy[4],
                                            double          *out_sigma,
                                            int             *out_nused)
{
    cpl_image   *tmp_img = NULL;
    double      bg_result = -1.0,
                sigma_result = -1.0;
    
    CLIPM_TRY
    {
        int             imsize[2],
                        wdwsize[2],
                        buffstart[2];
        int             nr_pix_below,
                        nused,
                        nbad,
                        ngood;
        double          backup_result_bg,
                        backup_result_sigma;
        int             backup_nused = 0;
        const double    Kappa = 2.58,               /* 99% */
                        Limit_lower_pixels = 0.10,  /* if 10% below distrib.*/
                        Limit_not_enough = 0.05;
        
        clipm_priv_checks_window_image(     window_xxyy,
                                            img,
                                            1,
                                            imsize,
                                            wdwsize,
                                            buffstart);
        CLIPM_TRY_CHECK_ERROR_STATE();
        clipm_priv_checks_imtype_any(       img, NULL);
        CLIPM_TRY_CHECK_ERROR_STATE();
        
        bg_result = clipm_priv_image_estimate_low_kappa_sigma(
                                            img,
                                            window_xxyy,
                                            Kappa,
                                            &sigma_result,
                                            &nused);
        CLIPM_TRY_CHECK_ERROR_STATE();
        /*cpl_msg_info(__func__, "bg: %.2f, sigma: %.2f, nused: %d", bg_result, sigma_result, nused);*/
        
        nr_pix_below = clipm_priv_image_pixel_count_below(
                                            img,
                                            window_xxyy,
                                            bg_result - Kappa * sigma_result,
                                            &nbad);
        CLIPM_TRY_ASSERT_ERROR_STATE();
            
        ngood = wdwsize[0]*wdwsize[1] - nbad;
        backup_result_bg = bg_result;
        backup_result_sigma = sigma_result;
        backup_nused = nused;
        
        /* if 10% below distribution (i.e. 11% below 2.58 sigma) */
        while ( (double)nr_pix_below / (double)ngood > Limit_lower_pixels + 0.01
                && nr_pix_below > 5)
        {
            if (tmp_img == NULL)
                tmp_img = cpl_image_extract(img,
                                            buffstart[0] + 1,
                                            buffstart[1] + 1,
                                            buffstart[0] + wdwsize[0],
                                            buffstart[1] + wdwsize[1]);
            CLIPM_TRY_ASSERT_ERROR_STATE();
            
            clipm_priv_image_bpm_reject_above(
                                            tmp_img,
                                            NULL,
                                            bg_result - Kappa * sigma_result);
            CLIPM_TRY_ASSERT_ERROR_STATE();
            
            /* if the last computation was fine, then backup the results
             * before trying again */
            if ((double)nused / (double)ngood >= Limit_not_enough)
            {
                backup_result_bg = bg_result;
                backup_result_sigma = sigma_result;
                backup_nused = nused;
            }
            
            /*cpl_msg_info(__func__, "iterating for lower pixels");*/
            bg_result = clipm_priv_image_estimate_low_kappa_sigma(
                                            tmp_img,
                                            NULL,
                                            Kappa,
                                            &sigma_result,
                                            &nused);
            if (CLIPM_ERROR_GET_NEW_SINCE_TRY() == CPL_ERROR_DATA_NOT_FOUND)
            {
                CLIPM_ERROR_RECOVER_TRYSTATE();
                nused = 0;
                break;
            }
            CLIPM_TRY_ASSERT_ERROR_STATE();
            /*cpl_msg_info(__func__, "bg: %.2f, sigma: %.2f, nused: %d", bg_result, sigma_result, nused);*/
            
            nr_pix_below = clipm_priv_image_pixel_count_below(
                                            tmp_img,
                                            window_xxyy,
                                            bg_result - Kappa * sigma_result,
                                            &nbad);
            CLIPM_TRY_ASSERT_ERROR_STATE();
        }
        
        if (tmp_img != NULL
            && (double)nused / (double)ngood < Limit_not_enough)
        {
            /* if really not enough pixels create a background, then
             * revert to last good values */
            bg_result = backup_result_bg;
            sigma_result = backup_result_sigma;
            nused = backup_nused;
            /*cpl_msg_info(__func__, "reverting to: bg: %.2f, sigma: %.2f, nused: %d", bg_result, sigma_result, nused);*/
        }
        
        if (out_sigma != NULL)
            *out_sigma = sigma_result;
        if (out_nused != NULL)
            *out_nused = nused;
    }
    CLIPM_CATCH
    {
        bg_result = -1.0;
        
        if (out_sigma != NULL)
            *out_sigma = -1.0;
        if (out_nused != NULL)
            *out_nused = 0;
    }
    
    cpl_image_delete(tmp_img);
    
    return bg_result;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief   Estimate the FWHM of a round object.
 * @param   img         Input image
 * @param   x_peakpos   Horizontal peak position
 * @param   y_peakpos   Vertical peak position
 * @param   bg_level    Input background level to subtract
 * @return  The full-widt-half-maximum, -1.0 in the case of error
 * 
 * @par Overview:
 * Beginning at the provided peak position, the diameter is searched, at which
 * the mean along the perimeter has fallen to half of the peak position.
 * \n\n
 * The applied principle is just a fast and rough estimation, with the purpose
 * of a quick estimation of a peak's size. But the radius-dependent brightness
 * evaluation uses the whole perimeter, so that noise in outer radiuses is
 * decreased. The peak is therefore expected to be round.
 * 
 * @par Bad Pixel Handling:
 * Bad pixel maps are supported. This means that bad pixels are ignored during
 * the computation.
 * 
 * @par Error Handling:
 * This function can set the following error codes:
 * - CPL_ERROR_ACCESS_OUT_OF_RANGE: if @a peakpos is outside the image range
 * - CPL_ERROR_NULL_INPUT: if @a img is NULL
 * - CPL_ERROR_INVALID_TYPE: @a img is not of type CPL_TYPE_INT,
 *   CPL_TYPE_FLOAT or CPL_TYPE_DOUBLE
 * - CPL_ERROR_CONTINUE:
 *   - the peak is darker than @a bg_level
 *   - the FWHM could not be found
 * - CPL_ERROR_DATA_NOT_FOUND: around the centre, only bad pixels could be
 *   found
 */
/*----------------------------------------------------------------------------*/
double          clipm_priv_image_estimate_fwhm_round(
                                            const cpl_image *img,
                                            double          x_peakpos,
                                            double          y_peakpos,
                                            double          bg_level)
{
    double      fwhm = -1.0;
    const cpl_binary
                *badp = NULL;
    
    CLIPM_TRY
    {
        cpl_type    type;
        int         w,
                    imsize[2];
        const void  *data;
        double      d,
                    peak,
                    halfpeak,
                    p,
                    plast,
                    r,
                    rlast,
                    r_cross = 0.0;
        
        clipm_priv_image_get_data_const(    img,
                                            NULL,
                                            1,
                                            imsize,
                                            NULL,
                                            NULL,
                                            &type,
                                            &data,
                                            &badp);
        CLIPM_TRY_CHECK_ERROR_STATE();
        
        CLIPM_TRY_CHECK_AUTOMSG(
                x_peakpos >= 1 && x_peakpos <= imsize[0] &&
                y_peakpos >= 1 && y_peakpos <= imsize[1],
            CPL_ERROR_ACCESS_OUT_OF_RANGE);
        
        /* make buffer position */
        x_peakpos--;
        y_peakpos--;
        
        /* diagonal image size as limit */
        d = sqrt(imsize[0]*imsize[0] + imsize[1]*imsize[1]);
        
        /* Find the central peak value. If it happens to be a bad pixel,
         * increase the radius and integrate along the perimeter, until
         * the computation was successful. */
        r = - 1.0;
        do
        {
            r += 1.0;
            CLIPM_ERROR_RECOVER_TRYSTATE();
            peak = clipm_priv_image_get_mean_perimeter(
                                            data, type,
                                            x_peakpos, y_peakpos,
                                            r,
                                            imsize[0], imsize[1],
                                            badp, NULL)
                    - bg_level;
        } while (   CLIPM_ERROR_GET_NEW_SINCE_TRY() ==
                    CPL_ERROR_DATA_NOT_FOUND &&
                    r < d);
        if (badp != NULL)
            CLIPM_ERROR_SET_MSG_IF_CODE(    CPL_ERROR_DATA_NOT_FOUND,
                                            "",
                                            "possibly too many bad pixels");
        CLIPM_TRY_CHECK_ERROR_STATE();
        
        CLIPM_TRY_CHECK(                    peak > 0,
                                            CPL_ERROR_CONTINUE,
                                            "",
                                            "peak is not brighter"
                                            " than background");
        
        halfpeak = peak / 2.0;
        
        /* Find the point where the mean along the perimeter has fallen
         * below half the peak. */
        p = peak;
        plast = p;
        rlast = r;
        w = 1;
        while (r < d)
        {
            double  tmp;
            r += 1.0;
            
            tmp = clipm_priv_image_get_mean_perimeter(
                                            data, type,
                                            x_peakpos, y_peakpos,
                                            r,
                                            imsize[0], imsize[1],
                                            badp,
                                            &w);
            if (w > 0)
            {
                p = tmp - bg_level;
                
                if (p <= halfpeak)
                {
                    if (p-plast != 0)
                        r_cross = (halfpeak-plast)*(r-rlast)/(p-plast) + rlast;
                    else
                        r_cross = (r + rlast)/2.0;
                    break;
                }
                else
                {
                    plast = p;
                    rlast = r;
                }
            }
            else
            {
                /* only bad pixels found, so jump over this radius step */
                CLIPM_ERROR_RECOVER_TRYSTATE();
            }
        }
        
        /* check output */
        if (r_cross >= d || r >= d || r_cross < 0.5 || w == 0)
            CLIPM_TRY_EXIT_WITH_ERROR_MSG(  CPL_ERROR_CONTINUE,
                                            "",
                                            "the fwhm could not be found");
        
        fwhm = 2 * r_cross;
    }
    CLIPM_CATCH
    {
        fwhm = -1.0;
    }
    
    return fwhm;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief   Estimate the FWHM of an object separable in x and y.
 * @param   img                 Input image
 * @param   xy_peakpos          Peak position (double buffer of length 2 [x, y])
 * @param   bg_level            Input background level to subtract
 * @param   out_xy_fwhm         Output full-width-half-maximum
 *                              (double buffer of length 2 [x, y])
 * @param   out_xy_middle       (Optional output) middle position between
 *                              left and right FWHM edges
 *                              (double buffer of length 2 [x, y]),
 *                              can be NULL
 * @param   out_xy_edge_sigma   (Optional output) sigma of the smoothing of the
 *                              edge, if not measurable -1.0, see also below,
 *                              (double buffer of length 2 [x, y]),
 *                              can be NULL
 * @return  CPL error code
 * 
 * @par Principle:
 * - Beginning at the provided peak position, the full-width-half-maximum is
 *   first searched along the crossing row and column. Thenafter, the region is
 *   collapsed in x and y in a region around (+/- edge + 2.6 times the underlying edge
 *   sigma, if this could not be determined, then in the +/- FWHM region),
 *   and the FWHM is searched in the collapsed (marginal) signals.
 * - The peak must be brighter than the background.
 * - The image size must be equal or greater than 3x3 pixels.
 * - If the object appears to be a flat object, then @a out_xy_edge_sigma will
 *   return the estimated sigma of the smoothing function in x and y
 *   respectively. If not measurable, or if the object does not appear to be
 *   a smoothed step function, then it will contain -1.0 respectively. 
 * 
 * @par Constraints:
 * Collapsing the signal in x and y means, that the result is only accurate
 * for separable functions, like squares or a gaussian.
 * 
 * @note
 * It is allowed to provide the same pointer for @a xy_peakpos and
 * @a out_xy_middle, for example for refinement purposes.
 * 
 * @par Bad Pixel Handling:
 * Bad pixel maps are supported. This means that bad pixels are ignored during
 * the computation.
 * 
 * @par Error Handling:
 * This function can set the following error codes:
 * - CPL_ERROR_NULL_INPUT: if @a img or @a xy_peakpos or @a out_xy_fwhm
 *   is NULL
 * - CPL_ERROR_ACCESS_OUT_OF_RANGE: if @a xy_peakpos is outside the image
 *   range
 * - CPL_ERROR_ILLEGAL_INPUT: the image size is less than 3x3
 * - CPL_ERROR_INVALID_TYPE: @a img is not of type CPL_TYPE_INT,
 *   CPL_TYPE_FLOAT or CPL_TYPE_DOUBLE
 * - CPL_ERROR_CONTINUE: the FWHM could not be found
 * - CPL_ERROR_DATA_NOT_FOUND: too many bad pixels prevent the computation of
 *   the FWHM
 * .
 * In the case of error, {-1.0, -1.0} is returned.
 * 
 * @todo
 * - test with flat object in unit test
 */
/*----------------------------------------------------------------------------*/
cpl_error_code  clipm_priv_image_estimate_fwhm_xy(
                                            const cpl_image *img,
                                            const double    *xy_peakpos,
                                            double          bg_level,
                                            double          *out_xy_fwhm,
                                            double          *out_xy_middle,
                                            double          *out_xy_edge_sigma)
{
    cpl_array   *xy_signals[2] = { NULL, NULL };
    int         dim;
    
    CLIPM_TRY
    {
        int         imsize[2],
                    cent[2],
                    maxpos[2];
        int         iteration_window[2][2];
        double      fwhm[2];
        
        /* check input */
        CLIPM_TRY_CHECK_AUTOMSG(            img != NULL,
                                            CPL_ERROR_NULL_INPUT);
        clipm_priv_checks_imtype_any(       img, NULL);
        CLIPM_TRY_CHECK_ERROR_STATE();
        
        CLIPM_TRY_CHECK_AUTOMSG(            xy_peakpos != NULL,
                                            CPL_ERROR_NULL_INPUT);
        CLIPM_TRY_CHECK_AUTOMSG(            out_xy_fwhm != NULL,
                                            CPL_ERROR_NULL_INPUT);
        
        imsize[0] = cpl_image_get_size_x(   img);
        imsize[1] = cpl_image_get_size_y(   img);
        CLIPM_TRY_CHECK_AUTOMSG(
                xy_peakpos[0] >= 1 && xy_peakpos[0] <= imsize[0] &&
                xy_peakpos[1] >= 1 && xy_peakpos[1] <= imsize[1],
            CPL_ERROR_ACCESS_OUT_OF_RANGE);

        CLIPM_TRY_CHECK(                    imsize[0] >=3 &&
                                            imsize[1] >=3,
                                            CPL_ERROR_ILLEGAL_INPUT,
                                            "image size",
                                            "must be >= 3x3");
        
        /* extract input */
        cent[0] = clipm_math_round_d2i(xy_peakpos[0]);
        cent[1] = clipm_math_round_d2i(xy_peakpos[1]);
        
        xy_signals[0] = clipm_priv_array_new_from_image_row(img, cent[1]);
        xy_signals[1] = clipm_priv_array_new_from_image_col(img, cent[0]);
        CLIPM_TRY_ASSERT_ERROR_STATE();
        
        /* get FWHM, and restrict regarded area to +/- 2.6*sigma edge */
        for (dim = 0; dim < 2; dim++)
        {
            double          middle,
                            slope;
            cpl_error_code  errc;
            fwhm[dim] = clipm_priv_array_estimate_fwhm(
                                            xy_signals[dim],
                                            xy_peakpos[dim] - 1.0,
                                            bg_level,
                                            &(maxpos[dim]),
                                            &middle,
                                            &slope);
            errc  = CLIPM_ERROR_GET_NEW_SINCE_TRY();
            CLIPM_TRY_ASSERT(               errc == CPL_ERROR_NONE ||
                                            errc == CPL_ERROR_CONTINUE ||
                                            errc == CPL_ERROR_DATA_NOT_FOUND);
            CLIPM_TRY_CHECK_ERROR_STATE();
/*printf("bg_level %.2f, height %.2f, slope %.2f\n",
bg_level,
cpl_array_get(
    xy_signals[dim],
    maxpos[dim], NULL) - bg_level,
slope);*/
            
            /* if determination of edge slope was successful */
            if (slope > 0)
            {
                double  new_width;
                new_width = _clipm_priv_image_get_optimal_box_width(
                                            cpl_array_get(
                                                xy_signals[dim],
                                                maxpos[dim], NULL) - bg_level,
                                            fwhm[dim],
                                            slope,
                                            /* 2.6 sigma * correcting factor
                                             * for the measured slope between
                                             * 30% and 70% in
                                             * clipm_priv_array_estimate_fwhm
                                             */
                                            2.6 * 0.9607,
                                            NULL,
                                            NULL);
                CLIPM_TRY_ASSERT_ERROR_STATE();
                /* FITS */
                iteration_window[dim][0] = floor(middle - new_width/2.0) + 1;
                iteration_window[dim][1] = ceil( middle + new_width/2.0) + 1;
/*printf("!!!! estimating window width, fwhm = %.2f, additional+/-: %.2f, new width: %.2f\n",
fwhm[dim], (new_width-fwhm[dim])/2.0, new_width);*/
            }
            else
            {   /* FITS */
                iteration_window[dim][0] = floor(middle - fwhm[dim]) + 1;
                iteration_window[dim][1] = ceil( middle + fwhm[dim]) + 1;
            }
            maxpos[dim] += 1;   /* FITS */

            clipm_priv_array_null(&xy_signals[dim]);
        }

        /* collapse image in iteration window */
        clipm_priv_checks_window_guarantee_image(
                                            *iteration_window,
                                            img,
                                            3);

        clipm_priv_image_collapse(          img,
                                            *iteration_window,
                                            &xy_signals[0],
                                            &xy_signals[1],
                                            NULL,
                                            NULL);
        CLIPM_TRY_CHECK_ERROR_STATE();

        /* and do it again */
        for (dim = 0; dim < 2; dim++)
        {
            double          middlepos,
                            max,
                            slope;
            int             max_ndx;
            cpl_error_code  errc;

            fwhm[dim] = clipm_priv_array_estimate_fwhm(
                                            xy_signals[dim],
                                            maxpos[dim]
                                                - iteration_window[dim][0],
                                            bg_level,
                                            &max_ndx,
                                            &middlepos,
                                            &slope);
            errc  = CLIPM_ERROR_GET_NEW_SINCE_TRY();
            CLIPM_TRY_ASSERT(               errc == CPL_ERROR_NONE ||
                                            errc == CPL_ERROR_CONTINUE ||
                                            errc == CPL_ERROR_DATA_NOT_FOUND);
            CLIPM_TRY_CHECK_ERROR_STATE();
            
            out_xy_fwhm[dim] = fwhm[dim];
            if (out_xy_middle != NULL)
                out_xy_middle[dim] =        middlepos
                                            + (double)iteration_window[dim][0];
            
            if (out_xy_edge_sigma != NULL)
            {
                double  edge_sigma_estimate;
                int     is_edge_sigma;
                max = cpl_array_get(xy_signals[dim], max_ndx, NULL);
                _clipm_priv_image_get_optimal_box_width(
                                            max - bg_level,
                                            fwhm[dim],
                                            slope,
                                            2.56,
                                            &edge_sigma_estimate,
                                            &is_edge_sigma);
                if (CLIPM_ERROR_GET_NEW_SINCE_TRY() != CPL_ERROR_NONE)
                {
                    CLIPM_ERROR_RECOVER_TRYSTATE();
                    is_edge_sigma = 0;
                }
                if (is_edge_sigma)
                    out_xy_edge_sigma[dim] = edge_sigma_estimate;
                else
                    out_xy_edge_sigma[dim] = -1.0;
            }
        }
    }
    CLIPM_CATCH
    {
        for (dim = 0; dim < 2; dim++)
        {
            if (out_xy_fwhm != NULL)
                out_xy_fwhm[dim] = -1.0;
            if (out_xy_middle != NULL)
                out_xy_middle[dim] = -1.0;
            if (out_xy_edge_sigma != NULL)
                out_xy_edge_sigma[dim] = -1.0;
        }
    }
    
    for (dim = 0; dim < 2; dim++)
        clipm_priv_array_null(&xy_signals[dim]);
    
    return CLIPM_ERROR_GET_NEW_SINCE_TRY();
}

/*----------------------------------------------------------------------------*/
/**
 * @brief   Determine mean and sigma by iteratively ignoring outliers.
 * @param   img             Input image
 * @param   window_xxyy     Coordinate buffer of the form
 *                          {xa, xb, ya, yb}, can be NULL,
 *                          minimum/maximum order is irrelevant
 * @param   kappa           Ignore (abs(value - mean) > kappa * sigma),
 *                          must be >= 1.0
 * @param   initial_limits  (Optional) double buffer of the form {l1, l2},
 *                          can be NULL,
 *                          minimum/maximum order is irrelevant
 * @param   nmax_iterations Maximum number of iterations, must be >= 0
 * @param   out_mean        (Optional) mean,
 *                          can be NULL
 * @param   out_sigma       (Optional) sigma,
 *                          can be NULL
 * @param   out_kappasigma  (Optional) kappa * sigma,
 *                          can be NULL
 * @param   out_nused       (Optional) number of remaining (used) pixels,
 *                          can be NULL
 * @param   out_niterations (Optional) number of performed iterations,
 *                          can be NULL
 * @return  CPL error code
 * 
 * @par Principle:
 * - If @a initial_limits is given, only values in this intial interval are
 *   used for the inital mean and sigma computation. This does not mean
 *   however, that the solution diverges out of this interval and later
 *   converges to a new range.
 * 
 * @par Bad Pixel Handling:
 * Bad pixel maps are supported. This means that bad pixels are ignored during
 * the computation.
 * 
 * @par Error Handling:
 * The following error codes can be set and returned:
 * - CPL_ERROR_NULL_INPUT: @a image is NULL
 * - CPL_ERROR_INVALID_TYPE: @a image is not of type CPL_TYPE_INT,
 *   CPL_TYPE_FLOAT or CPL_TYPE_DOUBLE
 * - CPL_ERROR_ACCESS_OUT_OF_RANGE: @a window_xxyy specifies coordinates outside
 *   the @a image plane
 * - CPL_ERROR_ILLEGAL_INPUT:
 *   - @a kappa < 1.0, or
 *   - @a nmax_iterations < 0, or
 *   - @a img contains less than 2 pixels
 * - CPL_ERROR_DATA_NOT_FOUND: less than 2 pixels are good
 * - CPL_ERROR_CONTINUE: all pixels were removed during iteration, please try
 *     a higher @a kappa (since @a kappa is limited to >= 1.0, this should
 *     normally not happen)
 * .
 * In the case of error, the following values are returned:
 * - @a out_mean: -1.0
 * - @a out_sigma: -1.0
 * - @a out_kappasigma: -1.0
 * - @a out_nused: 0
 * - @a out_niterations: 0
 * 
 * @todo
 * - test initial_limits in unit test, and that they can be left
 */
/*----------------------------------------------------------------------------*/
cpl_error_code  clipm_priv_image_get_kappa_sigma(
                                            const cpl_image *img,
                                            const int       window_xxyy[4],
                                            double          kappa,
                                            double          initial_limits[],
                                            int             nmax_iterations,
                                            double          *out_mean,
                                            double          *out_sigma,
                                            double          *out_kappasigma,
                                            int             *out_nused,
                                            int             *out_niterations)
{
    CLIPM_TRY
    {
        int         imsize[2],    /* x, y */
                    wdwsize[2],
                    start[2];
        cpl_type    type;
        const void  *data_start = NULL;
        const cpl_binary
                    *badp_start = NULL;
        
        /* init output */
        if (out_mean != NULL)
            *out_mean = -1.0;
        if (out_sigma != NULL)
            *out_sigma = -1.0;
        if (out_kappasigma != NULL)
            *out_kappasigma = -1.0;
        if (out_nused != NULL)
            *out_nused = 0;
        if (out_niterations != NULL)
            *out_niterations = 0;
        
        clipm_priv_image_get_data_const(    img,
                                            window_xxyy,
                                            1,
                                            imsize,
                                            wdwsize,
                                            start,
                                            &type,
                                            &data_start,
                                            &badp_start);
        CLIPM_TRY_CHECK_ERROR_STATE();
        
        CLIPM_TRY_CHECK_AUTOMSG(            kappa >= 1.0,
                                            CPL_ERROR_ILLEGAL_INPUT);
        CLIPM_TRY_CHECK_AUTOMSG(            nmax_iterations >= 0,
                                            CPL_ERROR_ILLEGAL_INPUT);
        
/* --- throw in a macro to ease life ---------------------------------------- */
#define clipm_priv_image_get_kappa_sigma_BODY( \
                                            TYPE) \
do { \
    TYPE    Ttemp, \
            l_limit = 0, \
            u_limit = 0; \
    int     N, \
            iter, \
            nused, \
            nused_last; \
    double  dtemp, \
            mean, \
            sigma; \
    \
    if (initial_limits != NULL) \
    { \
        l_limit = clipm_priv_math_min(initial_limits[0], initial_limits[1]); \
        u_limit = clipm_priv_math_max(initial_limits[0], initial_limits[1]); \
    } \
    \
    N = wdwsize[0] * wdwsize[1]; \
    CLIPM_TRY_CHECK(                        N >= 2, \
                                            CPL_ERROR_ILLEGAL_INPUT, \
                                            "img", \
                                            "must contain at least 2 pixels"); \
    \
    /* init mean */ \
    Ttemp = 0; \
    N = 0; \
    if (initial_limits == NULL) \
    { \
        clipm_priv_image_LOOP_GOOD_DATA_CONST(  TYPE, imsize, wdwsize, \
                                            data_start, badp_start, \
            \
            , \
            Ttemp += data[x]; \
            N++; \
            , \
            \
        ); \
    } \
    else \
    { \
        clipm_priv_image_LOOP_GOOD_DATA_CONST(  TYPE, imsize, wdwsize, \
                                            data_start, badp_start, \
            \
            , \
            if (data[x] >= l_limit && data[x] <= u_limit) \
            { \
                Ttemp += data[x]; \
                N++; \
            } \
            , \
            \
        ); \
    } \
    CLIPM_TRY_CHECK(                        N >= 2, \
                                            CPL_ERROR_DATA_NOT_FOUND, \
                                            "img", \
                                            "most pixels are bad"); \
    mean = (double)Ttemp / N; \
    \
    /* init sigma */ \
    sigma = 0; \
    if (initial_limits == NULL) \
    { \
        clipm_priv_image_LOOP_GOOD_DATA_CONST(  TYPE, imsize, wdwsize, \
                                            data_start, badp_start, \
            \
            , \
            dtemp = (double)data[x] - mean; \
            sigma += dtemp * dtemp; \
            , \
            \
        ); \
    } \
    else \
    { \
        clipm_priv_image_LOOP_GOOD_DATA_CONST(  TYPE, imsize, wdwsize, \
                                            data_start, badp_start, \
            \
            , \
            if (data[x] >= l_limit && data[x] <= u_limit) \
            { \
                dtemp = (double)data[x] - mean; \
                sigma += dtemp * dtemp; \
            } \
            , \
            \
        ); \
    } \
    sigma = sqrt(sigma / (N-1)); \
    \
    /* iterate */ \
    nused_last = 0; \
    nused = N; \
    for (iter = 0; iter < nmax_iterations; iter++) \
    { \
        double  clip, \
                dsqs = 0.0; \
        TYPE    sum = 0; \
        \
        clip = kappa * sigma; \
        nused = 0; \
        \
        clipm_priv_image_LOOP_GOOD_DATA_CONST(TYPE, imsize, wdwsize, \
                                            data_start, badp_start, \
            \
            , \
            dtemp = (double)data[x] - mean; \
            if (fabs(dtemp) <= clip) \
            { \
                nused++; \
                sum += data[x]; \
                dsqs += dtemp * dtemp; \
            } \
            , \
            \
        ); \
        \
        /* if nused == nused_last, we still need to compute the new sigma */ \
        sigma = sqrt(dsqs / (nused-1)); \
        if (nused <= 2 || nused == nused_last) \
            break; \
        \
        mean = (double)sum / nused; \
        nused_last = nused; \
    } \
    \
    /* this should not happen with a kappa >= 1.0 */ \
    CLIPM_TRY_CHECK(                        nused > 0, \
                                            CPL_ERROR_CONTINUE, \
                                            "", \
                                            "kappa too small?"); \
    \
    if (out_mean != NULL) \
        *out_mean = mean; \
    if (out_sigma != NULL) \
        *out_sigma = sigma; \
    if (out_kappasigma != NULL) \
        *out_kappasigma = kappa * sigma; \
    if (out_nused != NULL) \
        *out_nused = nused; \
    if (out_niterations != NULL) \
        *out_niterations = iter; \
    \
} while (0)
/* -------------------------------------------------------------------------- */

        switch(type)
        {
            case    CPL_TYPE_INT:
                clipm_priv_image_get_kappa_sigma_BODY(int);
                break;
            case    CPL_TYPE_FLOAT:
                clipm_priv_image_get_kappa_sigma_BODY(float);
                break;
            case    CPL_TYPE_DOUBLE:
                clipm_priv_image_get_kappa_sigma_BODY(double);
                break;
            default:
                CLIPM_TRY_EXIT_WITH_ERROR_MSG(
                                            CPL_ERROR_INVALID_TYPE,
                                            "image",
                                            "must be int, float or double");
        }
    }
    CLIPM_CATCH
    {
    }
    
    return CLIPM_ERROR_GET_NEW_SINCE_TRY();
}

/*----------------------------------------------------------------------------*/
/**
 * @brief   Compute the barycentre of an object.
 * @param   img             Input image
 * @param   window_xxyy     Coordinate buffer of the form
 *                          {xa, xb, ya, yb}, can be NULL,
 *                          minimum/maximum order is irrelevant
 * @param   bg_level        Input background level to subtract
 * @param   lower_cutlevel  Lower limit of values to use (separately from
 *                          @a bg_level) (values can be equal to or greater)
 * @param   out_xy_centre   Mandatory (double) buffer of size 2 to contain
 *                          the barycentre position after success
 * @param   out_weight      (Optional) output weight of the peak,
 *                          can be NULL
 * @param   out_nused       (Optional) output number of used values,
 *                          can be NULL
 * @return  CPL error code
 * 
 * @par Bad Pixel Handling:
 * Bad pixel maps are supported. This means that bad pixels are ignored during
 * the computation.
 * 
 * @note
 * With a combination of positive and negative values in the image, it is
 * possible that the returned barycentre is outside the image. If not tolerated,
 * this should be checked afterwards. This can be avoided by ensuring
 * lower_cutlevel >= bg_level.
 * 
 * @par Error Handling:
 * The following errors can be set and returned:
 * - CPL_ERROR_NULL_INPUT: @a image is NULL
 * - CPL_ERROR_INVALID_TYPE: @a image is not of type CPL_TYPE_INT,
 *   CPL_TYPE_FLOAT or CPL_TYPE_DOUBLE
 * - CPL_ERROR_ACCESS_OUT_OF_RANGE: @a window_xxyy specifies coordinates outside
 *   the @a image plane
 * - CPL_ERROR_DATA_NOT_FOUND: all pixels are bad
 * - CPL_ERROR_DIVISION_BY_ZERO:
 *   - no pixel was above @a lower_cutlevel
 *     (not that a division by zero has occurred, just one error has to be used
 *     for it)
 *   - the sum over the signal is close to zero
 * .
 * In the case of error, the following values are returned:
 * - @a out_xy_centre: {-1.0, -1.0}
 * - @a out_weight: 0
 * - @a out_nused: 0
 * 
 * @todo:
 * - test bad pixel handling in unit test
 */
/*----------------------------------------------------------------------------*/
cpl_error_code  clipm_priv_image_get_barycentre(
                                            const cpl_image *img,
                                            const int       window_xxyy[4],
                                            double          bg_level,
                                            double          lower_cutlevel,
                                            double          *out_xy_centre,
                                            double          *out_weight,
                                            int             *out_nused)
{
    CLIPM_TRY
    {
        int         imsize[2],    /* x, y */
                    wdwsize[2],
                    start[2];
        cpl_type    type;
        const void  *data_start = NULL;
        const cpl_binary
                    *badp_start = NULL;
        
        /* init output */
        if (out_xy_centre != NULL)
        {
            out_xy_centre[0] = -1.0;
            out_xy_centre[1] = -1.0;
        }
        if (out_weight != NULL)
            *out_weight = 0.0;
        if (out_nused != NULL)
            *out_nused = 0;
        
        clipm_priv_image_get_data_const(    img,
                                            window_xxyy,
                                            1,
                                            imsize,
                                            wdwsize,
                                            start,
                                            &type,
                                            &data_start,
                                            &badp_start);
        CLIPM_TRY_CHECK_ERROR_STATE();
        
        CLIPM_TRY_CHECK_AUTOMSG(            out_xy_centre != NULL,
                                            CPL_ERROR_NULL_INPUT);
        
        {
            int     nused = 0,
                    ngood = 0,
                    dim;
            double  val,
                    sum = 0.0;
            double  moment[2] = {0.0, 0.0};

/* --- throw in a macro to ease life ---------------------------------------- */
#define clipm_priv_image_get_barycentre_BODY( \
                                            TYPE) \
do { \
    clipm_priv_image_LOOP_GOOD_DATA_CONST(  TYPE, imsize, wdwsize, \
                                            data_start, badp_start, \
        \
        , \
        if ((double)data[x] >= lower_cutlevel) \
        { \
            val = (double)data[x] - bg_level; \
            sum += val; \
            moment[0] += val * x; \
            moment[1] += val * y; \
            nused++; \
        } \
        ngood++; \
        , \
        \
    ); \
} while (0)
/* -------------------------------------------------------------------------- */
            switch(type)
            {
                case    CPL_TYPE_INT:
                    clipm_priv_image_get_barycentre_BODY(int);
                    break;
                case    CPL_TYPE_FLOAT:
                    clipm_priv_image_get_barycentre_BODY(float);
                    break;
                case    CPL_TYPE_DOUBLE:
                    clipm_priv_image_get_barycentre_BODY(double);
                    break;
                default:
                    CLIPM_TRY_EXIT_WITH_ERROR_MSG(
                                            CPL_ERROR_INVALID_TYPE,
                                            "image",
                                            "must be int, float or double");
            }
            
            CLIPM_TRY_CHECK(                ngood >= 1,
                                            CPL_ERROR_DATA_NOT_FOUND,
                                            "img",
                                            "all pixels are bad");
            
            CLIPM_TRY_CHECK(                nused >= 1,
                                            CPL_ERROR_DIVISION_BY_ZERO,
                                            "img",
                                            "no pixel value above threshold");
            
            CLIPM_TRY_CHECK(                fabs(sum) > DBL_MIN,
                                            CPL_ERROR_DIVISION_BY_ZERO,
                                            "img",
                                            "sum of signal is close to zero");
            
            if (out_xy_centre != NULL)
                for (dim = 0; dim < 2; dim++)
                    out_xy_centre[dim] = (moment[dim] / sum)
                                        + start[dim] + 1.0; /* FITS */
            if (out_weight != NULL)
                *out_weight = sum;
            if (out_nused != NULL)
                *out_nused = nused;
        }
    }
    CLIPM_CATCH
    {
    }
    
    return CLIPM_ERROR_GET_NEW_SINCE_TRY();
}

/*----------------------------------------------------------------------------*/
/**
 * @brief   Compute the sigma (RMS) of an object (i.e. the PSF).
 * @param   img             Input image
 * @param   window_xxyy     Coordinate buffer of the form
 *                          {xa, xb, ya, yb}, can be NULL,
 *                          minimum/maximum order is irrelevant
 * @param   centre_xy       Peak position (double buffer of length 2 [x, y])
 * @param   bg_level        Input background level to subtract
 * @param   cut_lower       Lower limit for values to be taken into account
 * @param   cut_upper       Upper limit for values to be taken into account
 * @param   sigma_xy        (Output) Mandatory (double) buffer of size 2 to
 *                          contain the barycentre position after success
 * @param   gain            (Optional) detector gain value (input),
 *                          in unit "ADU / electron",
 *                          required for the computation of @a centre_err_xy,
 *                          will be ignored otherwise
 * @param   centre_err_xy   (Optional output) contains the centre error after
 *                          success (double buffer of length 2 [x, y]),
 *                          needs the @a gain to be set (!),
 *                          can be NULL
 * @return  CPL error code
 * 
 * @par Principle:
 * The lateral sigma (or PSF sigma) on an object is computed (or in better
 * words, the root mean square (RMS)). It assumes, that the centre has already
 * been well determined.
 * 
 * @par Bad Pixel Handling:
 * Bad pixel maps are supported. This means that bad pixels are ignored during
 * the computation.
 * 
 * @par Constraints:
 * - The object peak must be positive (brighter than the background).
 * - The image or image window size must be >= 3x3 pixels.
 * - If @a centre_err_xy is requested, @a gain must be > 0.
 * - Generally, at least 2 pixels must be not bad and fulfill
 *   @a cut_lower <= pixel <= @a cut_upper.
 * - Theoretically it is possible (although rarely the case) that extreme
 *   negative noise values can mess up the sigma computation.
 * - The integral of the signal divided by the gain must be positive
 * 
 * @par Error Handling:
 * The following errors can be set and returned:
 * - CPL_ERROR_NULL_INPUT:
 *   - @a image is NULL
 *   - @a centre_xy is NULL
 *   - @a sigma_xy is NULL
 * - CPL_ERROR_ILLEGAL_INPUT:
 *   - @a centre_err_xy is requested, but gain == 0.0
 *   - the image (window) size is not >= 3x3
 *   - cut_upper <= cut_lower
 * - CPL_ERROR_INVALID_TYPE: @a image is not of type CPL_TYPE_INT,
 *   CPL_TYPE_FLOAT or CPL_TYPE_DOUBLE
 * - CPL_ERROR_ACCESS_OUT_OF_RANGE: @a window_xxyy specifies coordinates outside
 *   the @a image plane
 * - CPL_ERROR_DATA_NOT_FOUND: all pixels (except none or 1) are bad
 * - CPL_ERROR_DIVISION_BY_ZERO: not at least 2 pixels fulfilled
 *   @a cut_lower <= pixel <= @a cut_upper
 *   (not that a division by zero has occurred, just one error has to be used
 *   for it)
 * - CPL_ERROR_CONTINUE:
 *   - extreme noise outliers led to invalid sigma computation
 *   - @a centre_err_xy is requested and was tried to be computed,
 *     but the flux is not positive or not sufficient
 * .
 * In the case of error, the following values are returned:
 * - @a sigma_xy: {-1.0, -1.0}
 * - @a centre_err_xy: {-1.0, -1.0}
 * 
 * @todo:
 * - implement unit test
 */
/*----------------------------------------------------------------------------*/
cpl_error_code  clipm_priv_image_get_psf_sigma(
                                            const cpl_image *img,
                                            const int       window_xxyy[4],
                                            const double    *centre_xy,
                                            double          bg_level,
                                            double          cut_lower,
                                            double          cut_upper,
                                            double          *sigma_xy,
                                            double          gain,
                                            double          *centre_err_xy)
{
    CLIPM_TRY
    {
        int         imsize[2],              /* x, y */
                    wdwsize[2],
                    start[2];
        int         dim;
        cpl_type    type;
        const void  *data_start = NULL;
        const cpl_binary
                    *badp_start = NULL;
        
        /* init output */
        for (dim = 0; dim < 2; dim++)
        {
            if (sigma_xy != NULL)
                sigma_xy[dim] = -1.0;
            if (centre_err_xy != NULL)
                centre_err_xy[dim] = -1.0;
        }
        
        CLIPM_TRY_CHECK_AUTOMSG(            img != NULL,
                                            CPL_ERROR_NULL_INPUT);
        CLIPM_TRY_CHECK_AUTOMSG(            centre_xy != NULL,
                                            CPL_ERROR_NULL_INPUT);
        CLIPM_TRY_CHECK_AUTOMSG(            sigma_xy != NULL,
                                            CPL_ERROR_NULL_INPUT);
        
        clipm_priv_image_get_data_const(    img,
                                            window_xxyy,
                                            1,
                                            imsize,
                                            wdwsize,
                                            start,
                                            &type,
                                            &data_start,
                                            &badp_start);
        CLIPM_TRY_CHECK_ERROR_STATE();
        
        CLIPM_TRY_CHECK_AUTOMSG(            gain != 0.0
                                            || centre_err_xy == NULL,
                                            CPL_ERROR_ILLEGAL_INPUT);
        CLIPM_TRY_CHECK_AUTOMSG(            cut_lower < cut_upper,
                                            CPL_ERROR_ILLEGAL_INPUT);
        CLIPM_TRY_CHECK        (            wdwsize[0] >=3 &&
                                            wdwsize[1] >=3,
                                            CPL_ERROR_ILLEGAL_INPUT,
                                            "image (window) size",
                                            "must be >= 3x3");
        
        {
            int     nused = 0,
                    ngood = 0;
            double  val,
                    dx,
                    dy,
                    sum = 0.0;
            double  moment[2] = {0.0, 0.0},
                    buffercentre[2];
            
            for (dim = 0; dim < 2; dim++)
                buffercentre[dim] = centre_xy[dim] - ((double)start[dim] + 1.0);
            
/* --- throw in a macro to ease life ---------------------------------------- */
#define clipm_priv_image_get_psf_sigma_BODY(TYPE) \
do { \
    clipm_priv_image_LOOP_GOOD_DATA_CONST(  TYPE, imsize, wdwsize, \
                                            data_start, badp_start, \
        dy = (double)y - buffercentre[1]; \
        , \
        dx = (double)x - buffercentre[0]; \
        if ((double)data[x] >= cut_lower && \
            (double)data[x] <= cut_upper) \
        { \
            val = (double)data[x] - bg_level; \
            sum += val; \
            moment[0] += val * dx * dx; \
            moment[1] += val * dy * dy; \
            nused++; \
        } \
        ngood++; \
        , \
        \
    ); \
} while (0)
/* -------------------------------------------------------------------------- */
            switch(type)
            {
                case    CPL_TYPE_INT:
                    clipm_priv_image_get_psf_sigma_BODY(int);
                    break;
                case    CPL_TYPE_FLOAT:
                    clipm_priv_image_get_psf_sigma_BODY(float);
                    break;
                case    CPL_TYPE_DOUBLE:
                    clipm_priv_image_get_psf_sigma_BODY(double);
                    break;
                default:
                    CLIPM_TRY_EXIT_WITH_ERROR_MSG(
                                            CPL_ERROR_INVALID_TYPE,
                                            "image",
                                            "must be int, float or double");
            }
            
            CLIPM_TRY_CHECK(                ngood >= 2,
                                            CPL_ERROR_DATA_NOT_FOUND,
                                            "img",
                                            "most pixels are bad");
            
            CLIPM_TRY_CHECK(                nused >= 2,
                                            CPL_ERROR_DIVISION_BY_ZERO,
                                            "img",
                                            "almost no pixel values inside"
                                            " cutlevels");
            
            if (sigma_xy != NULL)
                for (dim = 0; dim < 2; dim++)
                {
                    CLIPM_TRY_CHECK(        sum * moment[dim] > 0.0,
                                            CPL_ERROR_CONTINUE,
                                            "",
                                            "extreme noise or distortions");
                    sigma_xy[dim] = sqrt(moment[dim] / sum);
                }
            
            if (centre_err_xy != NULL)
            {
                double  nr_electrons;
                CLIPM_TRY_CHECK(            (nr_electrons = sum/gain) > 0.0,
                                            CPL_ERROR_CONTINUE,
                                            "centre error computation",
                                            "flux is negative");
                for (dim = 0; dim < 2; dim++)
                    centre_err_xy[dim] = sigma_xy[dim] / sqrt(nr_electrons);
            }
        }
            
    }
    CLIPM_CATCH
    {
        int dim;
        for (dim = 0; dim < 2; dim++)
        {
            if (sigma_xy != NULL)
                sigma_xy[dim] = -1.0;
            if (centre_err_xy != NULL)
                centre_err_xy[dim] = -1.0;
        }
    }
    
    return CLIPM_ERROR_GET_NEW_SINCE_TRY();
}

/*----------------------------------------------------------------------------*/
/**
 * @brief   Collapse an image (by averaging) in both dimensions.
 * @param   image           Input image
 * @param   window_xxyy     Window coordinate buffer of the form
 *                          {xa, xb, ya, yb}, can be NULL
 * @param   horizontal      (Output) horizontal array of type CPL_TYPE_DOUBLE
 * @param   vertical        (Output) vertical array of type CPL_TYPE_DOUBLE
 * @param   x_weight_map    (Optional output) horizontal weighting map
 *                          of type CPL_TYPE_INT
 * @param   y_weight_map    (Optional output) vertical weighting map
 *                          of type CPL_TYPE_INT
 * @return  CPL error code
 * 
 * @par Principle:
 * - The image is averaged horizontally (all columns) and stored in
 *   @a *vertical, and vertically (all rows) and stored in @a *horizontal.
 *   In other words, the array names reflect the remaining signal direction.
 * - If provided, the arrays @a x_weight_map and @a y_weight_map contain at
 *   the end the number of non-bad pixels that had contributed to an average
 *   value in the above vectors. NULL is returned in the case of error.
 * 
 * @par Bad Pixel Handling:
 * Bad pixel maps are supported. This means that bad pixels are omitted
 * during averaging. Completely bad rows or columns will lead to invalid
 * entries in @a horizontal or @a vertical, respectively.
 * 
 * @par Error Handling:
 * In the case of error, @a horizontal and @a vertical are set to NULL.@n
 * The following error codes are reported:
 * - CPL_ERROR_NULL_INPUT: if @a image, @a horizontal or @a vertical is NULL
 * - CPL_ERROR_INVALID_TYPE if @a image is not of type int, float or double
 * - CPL_ERROR_ACCESS_OUT_OF_RANGE if @a window contains coordinates outside the
 *   image range
 * 
 * @todo
 * - implement unit test
 */
/*----------------------------------------------------------------------------*/
cpl_error_code  clipm_priv_image_collapse(  const cpl_image *image,
                                            const int       window_xxyy[4],
                                            cpl_array       **horizontal,
                                            cpl_array       **vertical,
                                            cpl_array       **x_weight_map,
                                            cpl_array       **y_weight_map)
{
    cpl_array   *xweight = NULL,
                *yweight = NULL;
    
    CLIPM_TRY
    {
        double      *xmeandata, /* x signal (collapsed in y) */
                    *ymeandata; /* y signal (collapsed in x) */
        int         *xweightdata = NULL,
                    *yweightdata = NULL;
        int         imsize[2],
                    windowsize[2];
        int         x,
                    y;
        cpl_type    type;
        const void  *data_start = NULL;
        const cpl_binary
                    *badp_start = NULL;

        /* initialize output */
        clipm_priv_array_null(horizontal);
        clipm_priv_array_null(vertical);
        clipm_priv_array_null(x_weight_map);
        clipm_priv_array_null(y_weight_map);
        
        /* check input */
        CLIPM_TRY_CHECK_AUTOMSG(image != NULL,      CPL_ERROR_NULL_INPUT);
        CLIPM_TRY_CHECK_AUTOMSG(horizontal != NULL, CPL_ERROR_NULL_INPUT);
        CLIPM_TRY_CHECK_AUTOMSG(vertical != NULL,   CPL_ERROR_NULL_INPUT);
        
        clipm_priv_image_get_data_const(    image,
                                            window_xxyy,
                                            1,
                                            imsize,
                                            windowsize,
                                            NULL,
                                            &type,
                                            &data_start,
                                            &badp_start);
        CLIPM_TRY_CHECK_ERROR_STATE();
        
        (*horizontal) = cpl_array_new(windowsize[0], CPL_TYPE_DOUBLE);
        cpl_array_fill_window_double(*horizontal, 0, windowsize[0], 0.0);
        xmeandata = cpl_array_get_data_double(*horizontal);
        (*vertical) = cpl_array_new(windowsize[1], CPL_TYPE_DOUBLE);
        cpl_array_fill_window_double(*vertical, 0, windowsize[1], 0.0);
        ymeandata = cpl_array_get_data_double(*vertical);
        CLIPM_TRY_ASSERT_ERROR_STATE();
        
        if (badp_start != NULL)
        {
            xweight = cpl_array_new(windowsize[0], CPL_TYPE_INT);
            cpl_array_fill_window_int(xweight, 0, windowsize[0], 0);
            xweightdata = cpl_array_get_data_int(xweight);
            yweight = cpl_array_new(windowsize[1], CPL_TYPE_INT);
            cpl_array_fill_window_int(yweight, 0, windowsize[1], 0);
            yweightdata = cpl_array_get_data_int(yweight);
            CLIPM_TRY_ASSERT_ERROR_STATE();
        }
        
/* -------------------------------------------------------------------------- */
/* Macro for collapsing the image at the same time in both dimensions */
#define         clipm_COLLAPSE_IMAGE_TO_ARRAY_BODY(  TYPE) \
do { \
    clipm_priv_image_LOOP_GOOD_DATA_CONST(  TYPE, imsize, windowsize, \
                                            data_start, badp_start, \
        \
        , \
        xmeandata[x] += data[x]; \
        ymeandata[y] += data[x]; \
        , \
        xweightdata[x]++; \
        yweightdata[y]++; \
    ); \
} while (0)
/* -------------------------------------------------------------------------- */

        /* collapse it */
        switch (type) {
            case CPL_TYPE_INT: {
                clipm_COLLAPSE_IMAGE_TO_ARRAY_BODY(  int);
                break;
            }
            case CPL_TYPE_FLOAT: {
                clipm_COLLAPSE_IMAGE_TO_ARRAY_BODY(  float);
                break;
            }
            case CPL_TYPE_DOUBLE: {
                clipm_COLLAPSE_IMAGE_TO_ARRAY_BODY(  double);
                break;
            }
            default:
                CLIPM_TRY_EXIT_WITH_ERROR_MSG(
                                            CPL_ERROR_INVALID_TYPE,
                                            "input images", 
                                            "must be int, float or double");
        }
        
        /* normalize */
        if (badp_start != NULL)
        {
            for (x = 0; x < windowsize[0]; x++)
            {
                if (xweightdata[x] > 0)
                    xmeandata[x] /= xweightdata[x];
                else
                    cpl_array_set_invalid(*horizontal, x);
            }
            for (y = 0; y < windowsize[1]; y++)
            {
                if (yweightdata[y] > 0)
                    ymeandata[y] /= yweightdata[y];
                else
                    cpl_array_set_invalid(*vertical, y);
            }
            
            if (x_weight_map != NULL)
            {
                *x_weight_map = xweight;
                xweight = NULL;
            }
            if (y_weight_map != NULL)
            {
                *y_weight_map = yweight;
                yweight = NULL;
            }
        }
        else
        {
            for (x = 0; x < windowsize[0]; x++)
                xmeandata[x] /= windowsize[1];
            for (y = 0; y < windowsize[1]; y++)
                ymeandata[y] /= windowsize[0];
            
            /* check if weight maps were requested */
            if (x_weight_map != NULL && *x_weight_map == NULL)
            {
                *x_weight_map = cpl_array_new(windowsize[0], CPL_TYPE_INT);
                cpl_array_fill_window_int(  *x_weight_map,
                                            0,
                                            windowsize[0],
                                            windowsize[1]); /* weighting */
            }
            if (y_weight_map != NULL && *y_weight_map == NULL)
            {
                *y_weight_map = cpl_array_new(windowsize[1], CPL_TYPE_INT);
                cpl_array_fill_window_int(  *y_weight_map,
                                            0,
                                            windowsize[1],
                                            windowsize[0]); /* weighting */
            }
            CLIPM_TRY_ASSERT_ERROR_STATE();
        }
        
    }
    CLIPM_CATCH
    {
        clipm_priv_array_null(horizontal);
        clipm_priv_array_null(vertical);
        clipm_priv_array_null(x_weight_map);
        clipm_priv_array_null(y_weight_map);
    }
    
    cpl_array_delete(xweight);
    cpl_array_delete(yweight);
    
    return CLIPM_ERROR_GET_NEW_SINCE_TRY();
}

/*----------------------------------------------------------------------------*/
/**
 * @brief   Count all the pixels (strictly) below a given limit.
 * @param   img             Input image
 * @param   window_xxyy     Coordinate buffer of the form
 *                          {xa, xb, ya, yb}, can be NULL,
 *                          minimum/maximum order is irrelevant
 * @param   limit           Limit
 * @param   out_nbad        (Optional output) number of bad pixels in
 *                          image (window), returns -1 in the case of error
 * @return  Number of pixels below @a limit, 0 in the case of error
 * 
 * @par Bad Pixel Handling:
 * Bad pixels are not counted.
 * 
 * @par Error Handling:
 * The following errors can occur:
 * - CPL_ERROR_NULL_INPUT: @a img is NULL
 * - CPL_ERROR_INVALID_TYPE: @a image is not of type CPL_TYPE_INT,
 *   CPL_TYPE_FLOAT or CPL_TYPE_DOUBLE
 * - CPL_ERROR_ACCESS_OUT_OF_RANGE: @a window_xxyy specifies coordinates outside
 *   the image range
 */
/*----------------------------------------------------------------------------*/
int             clipm_priv_image_pixel_count_below(
                                            const cpl_image *img,
                                            const int       window_xxyy[4],
                                            double          limit,
                                            int             *out_nbad)
{
    int result = 0,
        ngood = 0;
    
    CLIPM_TRY
    {
        int         imsize[2],
                    wdwsize[2];
        cpl_type    type;
        const void  *data_start = NULL;
        const cpl_binary
                    *badp_start = NULL;
        
        clipm_priv_image_get_data_const(    img,
                                            window_xxyy,
                                            1,
                                            imsize,
                                            wdwsize,
                                            NULL,
                                            &type,
                                            &data_start,
                                            &badp_start);
        CLIPM_TRY_CHECK_ERROR_STATE();
        
/* --- throw in a macro to ease life ---------------------------------------- */
#define clipm_priv_image_pixel_count_below_BODY(  TYPE) \
do { \
    TYPE    tlim = (TYPE) limit; \
    clipm_priv_image_LOOP_GOOD_DATA_CONST(  TYPE, imsize, wdwsize, \
                                            data_start, badp_start, \
        \
        , \
        if (data[x] < tlim) \
            result++; \
        , \
        ngood++; \
    ); \
} while (0)
/* -------------------------------------------------------------------------- */
        switch(type)
        {
            case    CPL_TYPE_INT:
                clipm_priv_image_pixel_count_below_BODY(int);
                break;
            case    CPL_TYPE_FLOAT:
                clipm_priv_image_pixel_count_below_BODY(float);
                break;
            case    CPL_TYPE_DOUBLE:
                clipm_priv_image_pixel_count_below_BODY(double);
                break;
            default:
                CLIPM_TRY_EXIT_WITH_ERROR_MSG(
                                            CPL_ERROR_INVALID_TYPE,
                                            "image",
                                            "must be int, float or double");
        }
        
        if (out_nbad != NULL)
        {
            if (badp_start != NULL)
                *out_nbad = wdwsize[0]*wdwsize[1] - ngood;
            else
                *out_nbad = 0;
        }
    }
    CLIPM_CATCH
    {
        if (out_nbad != NULL)
            *out_nbad = -1;
    }
    
    return result;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*int             clipm_priv_image_pixel_healing(
                                            const cpl_image *img,
                                            const int       window_xxyy[4],
                                            double          deviation_limit,
                                            int             heal_dark,
                                            int             heal_bright,
                                            int             direction,
                                            int             *out_nbad)
{
    int result = 0,
        ngood = 0;
    
    CLIPM_TRY
    {
        int         imsize[2],
                    wdwsize[2],
                    start[2];
        cpl_type    type;
        const void  *data_start = NULL;
        const cpl_binary
                    *badp_start = NULL;
        
        clipm_priv_image_get_data_const(    img,
                                            window_xxyy,
                                            1,
                                            imsize,
                                            wdwsize,
                                            start,
                                            &type,
                                            &data_start,
                                            &badp_start);
        CLIPM_TRY_CHECK_ERROR_STATE();
        
* --- throw in a macro to ease life ---------------------------------------- *
#define clipm_priv_image_pixel_count_below_BODY(  TYPE) \
do { \
    TYPE    tlim = (TYPE) limit; \
    clipm_priv_image_LOOP_GOOD_DATA_CONST(  TYPE, imsize, wdwsize, \
                                            data_start, badp_start, \
        \
        , \
        if (data[x] < tlim) \
            result++; \
        , \
        ngood++; \
    ); \
} while (0)
* -------------------------------------------------------------------------- *
        switch(type)
        {
            case    CPL_TYPE_INT:
                clipm_priv_image_pixel_count_below_BODY(int);
                break;
            case    CPL_TYPE_FLOAT:
                clipm_priv_image_pixel_count_below_BODY(float);
                break;
            case    CPL_TYPE_DOUBLE:
                clipm_priv_image_pixel_count_below_BODY(double);
                break;
            default:
                CLIPM_TRY_EXIT_WITH_ERROR_MSG(
                                            CPL_ERROR_INVALID_TYPE,
                                            "image",
                                            "must be int, float or double");
        }
        
        if (out_nbad != NULL)
        {
            if (badp_start != NULL)
                *out_nbad = wdwsize[0]*wdwsize[1] - ngood;
            else
                *out_nbad = 0;
        }
    }
    CLIPM_CATCH
    {
        if (out_nbad != NULL)
            *out_nbad = -1;
    }
    
    return result;
}
*/

/*----------------------------------------------------------------------------*/
/**
 * @brief   Convolve an image with a kernel stored in a matrix.
 * @param   input       Input image
 * @param   kernel      Kernel
 * @param   window_xxyy Coordinate buffer of the form
 *                      {xa, xb, ya, yb}, can be NULL,
 *                      minimum/maximum order is irrelevant
 * @param   extend_bpm  Flag whether to extend the bad pixel flag to the
 *                      range of coverage of the kernel (incl. the image
 *                      border), otherwise normalise to covered kernel weight
 * @param   int2double  Flag whether to convert integer to double (see below)
 * @return  Convolved image, NULL in the case of error
 * 
 * @par Principle:
 * - An image is convolved with a kernel stored in a matrix.
 * - The matrix is used as literally defined and displayed on the screen
 *   (this means, the data buffer content is treated as vertically mirrored
 *   compared to an image).
 * - If @a window_xxyy is given, an output image with the corresponding size
 *   is created. The surrounding area in the input image, which is covered by
 *   the kernel, is used.
 * - If @a extend_bpm is true, then every pixel seeing a bad pixel or the border
 *   of the input image within the range of coverage of the kernel, is flagged
 *   as bad.
 * - If @a extend_bpm is false, then every output pixel is weighted with
 *   (kernel_sum / covered_kernel_sum), to compensate for gaps by bad pixels or
 *   input image borders. In this case, the sum over the kernel must be
 *   different from zero! (For example, a gradient operator, usually summing up
 *   to zero makes no sense here, this should be used with
 *   @a extend_bpm = true). If (kernel_sum / covered_kernel_sum) >= 1e3, then
 *   the respective pixel is flagged as bad (this is particularly useful for
 *   kernels containing zeros, like a shift operator).
 * - The centre of the kernel matrix is interpreted as the center of the
 *   impulse response.
 * - This function uses float or double intermediate values. If the input image
 *   is of type CPL_TYPE_INT, then double intermediate values are used. Per
 *   default, in this case the result is rounded back to integer. If
 *   @a int2double is true, then the result will be left as it is, this means
 *   the output for an input of type integer will be of type double.
 * 
 * @par Constraints:
 * - The kernel size must be odd.
 * - If @a extend_bpm is false, then the sum over the kernel must be different
 *   from zero.
 * 
 * @par Bad Pixel Handling:
 * - See above.
 * 
 * @par Error Handling:
 * The following errors may occur:
 * - CPL_ERROR_NULL_INPUT: @a input or @a kernel is NULL
 * - CPL_ERROR_ACCESS_OUT_OF_RANGE: a window coordinate exceeds the image range
 * - CPL_ERROR_ILLEGAL_INPUT:
 *   - the kernel size is not odd
 *   - @a extend_bpm == false and the sum over the kernel is close to zero
 *     (< 1e-3)
 * - CPL_ERROR_INVALID_TYPE: @a input is not of type CPL_TYPE_INT,
 *   CPL_TYPE_FLOAT or CPL_TYPE_DOUBLE
 * 
 * @todo:
 * - implement unit test
 */
/*----------------------------------------------------------------------------*/
cpl_image       *clipm_priv_image_conv_matrix(
                                            const cpl_image     *input,
                                            const cpl_matrix    *kernel,
                                            const int           window_xxyy[4],
                                            int                 extend_bpm,
                                            int                 int2double)
{
    cpl_image   *result = NULL,
                *norm = NULL,
                *tmpimg = NULL;
    
    CLIPM_TRY
    {
        int             imsize[2],
                        outsize[2],
                        wdwbufferstart[2],
                        kernelsize[2],
                        koffset[2];
        const void      *indata;
        void            *outdata;
        double          *normdata = NULL;
        double          kernelnorm = 0.0;
        int             xk,
                        yk,
                        dim;
        cpl_type        type;
        const cpl_binary
                        *inbpmdata;
        cpl_mask        *outbpm = NULL;
        cpl_binary      *outbpmdata = NULL;
        
        clipm_priv_checks_window_image(     window_xxyy,
                                            input,
                                            1,
                                            imsize,
                                            outsize,
                                            wdwbufferstart);
        CLIPM_TRY_CHECK_ERROR_STATE();
        clipm_priv_checks_imtype_any(       input, &type);
        CLIPM_TRY_CHECK_ERROR_STATE();
        
        CLIPM_TRY_CHECK_AUTOMSG(            kernel != NULL,
                                            CPL_ERROR_NULL_INPUT);
        kernelsize[0] = cpl_matrix_get_ncol(kernel);
        kernelsize[1] = cpl_matrix_get_nrow(kernel);
        
        CLIPM_TRY_CHECK_AUTOMSG(            kernelsize[0] % 2,
                                            CPL_ERROR_ILLEGAL_INPUT);
        CLIPM_TRY_CHECK_AUTOMSG(            kernelsize[1] % 2,
                                            CPL_ERROR_ILLEGAL_INPUT);
        
        indata = cpl_image_get_data_const(input);
        inbpmdata = clipm_priv_image_bpm_get_if_exist(input);
        
        /* create output image, if input is integer, then create double
         * as intermediate */
        if (type == CPL_TYPE_INT)
            result = cpl_image_new(outsize[0], outsize[1], CPL_TYPE_DOUBLE);
        else
            result = cpl_image_new(outsize[0], outsize[1], type);
        CLIPM_TRY_ASSERT_ERROR_STATE();
        
        outdata = cpl_image_get_data(result);
        if (inbpmdata != NULL || (extend_bpm))
        {
            outbpm = cpl_image_get_bpm(result);
            outbpmdata = cpl_mask_get_data(outbpm);
        }
        
        for (dim = 0; dim < 2; dim++)
        {
            koffset[dim] = (kernelsize[dim] - 1) / 2;
        }
        
        if (!extend_bpm)
        {
            /* partial kernels may be applied, so normalize at the end */
            kernelnorm =    cpl_matrix_get_mean(kernel)
                            * kernelsize[0]
                            * kernelsize[1];
            CLIPM_TRY_CHECK(                fabs(kernelnorm) > 1e-3,
                                            CPL_ERROR_ILLEGAL_INPUT,
                                            "kernel",
                                            "must have a mean different from "
                                            "zero, if (extend_bpm==false)");
            
            norm = cpl_image_new(outsize[0], outsize[1], CPL_TYPE_DOUBLE);
            normdata = cpl_image_get_data_double(norm);
            CLIPM_TRY_ASSERT_ERROR_STATE();
            /* set all pixels bad, to set them good again later */
            if (inbpmdata != NULL)
                cpl_mask_not(outbpm);
        }
        
        /* loop over kernel */
        for (yk = -koffset[1]; yk < kernelsize[1]-koffset[1]; yk++)
            for (xk = -koffset[0]; xk < kernelsize[0]-koffset[0]; xk++)
            {
                double      kval;
                int         outbuffstart[2] = {0, 0},
                            inbuffstart[2],
                            inbuffend[2];
                int         x,
                            y,
                            rowlength,
                            indata_offset,
                            outdata_offset;
                
                kval = cpl_matrix_get(kernel, yk+koffset[1], xk+koffset[0]);
                
                inbuffstart[0] = wdwbufferstart[0] - xk;
                inbuffstart[1] = wdwbufferstart[1] + yk; /* inverted y index */
                
                for (dim = 0; dim < 2; dim++)
                {
                    inbuffend[dim] = inbuffstart[dim] + outsize[dim];
                    outbuffstart[dim] = 0;
                    /* shrink buffer window if exceeds input image range */
                    if (inbuffstart[dim] < 0)
                    {
                        outbuffstart[dim] = - inbuffstart[dim];
                        inbuffstart[dim] = 0;
                    }
                    inbuffend[dim] =    inbuffend[dim] < imsize[dim] ?
                                        inbuffend[dim] : imsize[dim];
                }
                rowlength = inbuffend[0] - inbuffstart[0];
                
                indata_offset = inbuffstart[1]*imsize[0] + inbuffstart[0];
                outdata_offset = outbuffstart[1]*outsize[0] + outbuffstart[0];
                
                
/* --- throw in a macro to ease life ---------------------------------------- */
#define     clipm_priv_image_conv_BODY(     INTYPE, \
                                            OUTTYPE) \
do { \
                const INTYPE    *inrowdata; \
                OUTTYPE         *outrowdata; \
                inrowdata = (const INTYPE*)indata + indata_offset; \
                outrowdata = (OUTTYPE*)outdata + outdata_offset; \
                \
                if (inbpmdata != NULL) \
                { \
                    const cpl_binary    *inrowbpm; \
                    cpl_binary          *outrowbpm; \
                    inrowbpm = inbpmdata + indata_offset; \
                    outrowbpm = outbpmdata + outdata_offset; \
                    \
                    if (normdata != NULL) /* normalize to used kernel part */ \
                    { \
                        double  *normrowdata; \
                        normrowdata = normdata + outdata_offset; \
                        \
                        for (y = inbuffstart[1]; y < inbuffend[1]; y++) \
                        { \
                            for (x = 0; x < rowlength; x++) \
                            { \
                                if (! inrowbpm[x]) \
                                { \
                                    outrowdata[x] += kval * inrowdata[x]; \
                                    normrowdata[x] += kval; \
                                    /* flag good (there was contribution) */ \
                                    outrowbpm[x] = CPL_BINARY_0; \
                                } \
                            } \
                            inrowdata += imsize[0]; \
                            normrowdata += outsize[0]; \
                            outrowdata += outsize[0]; \
                            inrowbpm += imsize[0]; \
                            outrowbpm += outsize[0];\
                        } \
                    } \
                    else /* flag bad border later */ \
                    { \
                        for (y = inbuffstart[1]; y < inbuffend[1]; y++) \
                        { \
                            for (x = 0; x < rowlength; x++) \
                            { \
                                if (! inrowbpm[x]) \
                                    outrowdata[x] += kval * inrowdata[x]; \
                                else \
                                    /* flag as bad, a pixel was missing */ \
                                    outrowbpm[x] = CPL_BINARY_1; \
                            } \
                            inrowdata += imsize[0]; \
                            outrowdata += outsize[0]; \
                            inrowbpm += imsize[0]; \
                            outrowbpm += outsize[0];\
                        } \
                    } \
                } \
                else    /* if (inbpmdata == NULL) */ \
                { \
                    if (normdata != NULL) /* normalize to used kernel part */ \
                    { \
                        double  *normrowdata; \
                        normrowdata = normdata + outdata_offset; \
                        \
                        for (y = inbuffstart[1]; y < inbuffend[1]; y++) \
                        { \
                            for (x = 0; x < rowlength; x++) \
                            { \
                                outrowdata[x] += kval * inrowdata[x]; \
                                normrowdata[x] += kval; \
                            } \
                            inrowdata += imsize[0]; \
                            normrowdata += outsize[0]; \
                            outrowdata += outsize[0]; \
                        } \
                    } \
                    else /* flag bad border later */ \
                    { \
                        for (y = inbuffstart[1]; y < inbuffend[1]; y++) \
                        { \
                            for (x = 0; x < rowlength; x++) \
                            { \
                                outrowdata[x] += kval * inrowdata[x]; \
                            } \
                            inrowdata += imsize[0]; \
                            outrowdata += outsize[0]; \
                        } \
                    } \
                } \
} while (0)
/* -------------------------------------------------------------------------- */
            
                switch (type)
                {
                    case CPL_TYPE_INT:
                        clipm_priv_image_conv_BODY(int, double);
                        break;
                    case CPL_TYPE_FLOAT:
                        clipm_priv_image_conv_BODY(float, float);
                        break;
                    case CPL_TYPE_DOUBLE:
                        clipm_priv_image_conv_BODY(double, double);
                        break;
                    default:
                        CLIPM_TRY_EXIT_WITH_ERROR_MSG(
                                            CPL_ERROR_INVALID_TYPE,
                                            "input image", 
                                            "must be int, float or double");
                }
            }
        
        if (normdata == NULL) /* if extend bad pixels, then also flag borders */
        {
            int llwidth[2],
                urwidth[2];
            for (dim = 0; dim < 2; dim++)
            {
                llwidth[dim] =  koffset[dim] - wdwbufferstart[dim];
                llwidth[dim] = llwidth[dim] < 0 ? 0 : llwidth[dim];
                urwidth[dim] =  wdwbufferstart[dim] + outsize[dim]
                                + koffset[dim]
                                - imsize[dim];
                urwidth[dim] = urwidth[dim] < 0 ? 0 : urwidth[dim];
            }
            
            clipm_priv_image_bpm_border_bad(result,
                                            llwidth[0],
                                            urwidth[0],
                                            llwidth[1],
                                            urwidth[1]);
            CLIPM_TRY_ASSERT_ERROR_STATE();
        }
        else    /* if not extend bad pixels, then normalize */
        {
            /* get output bpm to also flag pixels with unsufficient kernel
             * contribution as bad */
            if (outbpmdata == NULL)
            {
                outbpm = cpl_image_get_bpm(result);
                outbpmdata = cpl_mask_get_data(outbpm);
                CLIPM_TRY_ASSERT(outbpmdata != NULL);
            }
            
/* --- throw in a macro to ease life ---------------------------------------- */
#define     clipm_priv_image_conv_NORM(     TYPE) \
do { \
                TYPE  *outalldata; \
                int tot, \
                    n; \
                \
                outalldata = (TYPE*)outdata; \
                tot = outsize[0] * outsize[1]; \
                \
                for (n = 0; n < tot; n++) \
                { \
                    if (*outbpmdata == CPL_BINARY_0) \
                    { \
                        double  factor; \
                        factor = (*normdata) / kernelnorm; \
                        if (factor > 1e-3) \
                            *outalldata /= factor; \
                        else \
                            *outbpmdata = CPL_BINARY_1; /* bad */ \
                    } \
                    normdata++; \
                    outalldata++; \
                    outbpmdata++; \
                } \
} while (0)
/* -------------------------------------------------------------------------- */
                switch (type)
                {
                    case CPL_TYPE_INT: /* double intermediate */
                        clipm_priv_image_conv_NORM(double);
                        break;
                    case CPL_TYPE_FLOAT:
                        clipm_priv_image_conv_NORM(float);
                        break;
                    case CPL_TYPE_DOUBLE:
                        clipm_priv_image_conv_NORM(double);
                        break;
                    default:
                        break;
                }
        }
        
        /* remove bad pixel mask, if all are good */
        if (cpl_image_count_rejected(result) == 0)
            cpl_image_accept_all(result);
        
        /* convert back to integer if input was int */
        if (type == CPL_TYPE_INT && !int2double)
        {
            tmpimg = result;
            result = clipm_priv_image_extract_round(
                                            tmpimg,
                                            NULL);
            CLIPM_TRY_CHECK_ERROR_STATE();
            
            clipm_priv_image_null(&tmpimg);
        }
    }
    CLIPM_CATCH
    {
        clipm_priv_image_null(&result);
    }
    
    cpl_image_delete(norm);
    cpl_image_delete(tmpimg);
    
    return result;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief   Convolve an image with a gaussian bell curve.
 * @param   image           Input image
 * @param   window_xxyy     Coordinate buffer of the form
 *                          {xa, xb, ya, yb}, can be NULL,
 *                          minimum/maximum order is irrelevant
 * @param   sigma           Sigma of the gaussian
 * @return  The new softened image, NULL on error
 * 
 * @par Principle:
 * - Since a 2d-Gaussian is separable, actually two convolutions are performed,
 *   one in x and y respectively. This leads to slightly different results if
 *   bad pixels are present, which must be considered neglegible if using
 *   this function.
 * - If @a window_xxyy is specified, the surrounding signal is also regarded
 *   in the convolution (to an extent of 4 sigma).
 * 
 * @par Constraints:
 * - Images can be of type CPL_TYPE_INT, CPL_TYPE_FLOAT or CPL_TYPE_DOUBLE.
 * - @a sigma must be greater than or equal to 0.
 * 
 * @par Error Handling:
 * Possible error codes set in this function:
 * - CPL_ERROR_NULL_INPUT if @a image is NULL
 * - CPL_ERROR_ACCESS_OUT_OF_RANGE: a window coordinate is outside the image
 * - CPL_ERROR_ILLEGAL_INPUT if @a sigma <= 0
 * - CPL_ERROR_INVALID_TYPE if the passed image type is not supported
 * 
 * @todo
 * - implement unit test
 */
/*----------------------------------------------------------------------------*/
cpl_image       *clipm_priv_image_filter_lowpass(
                                            const cpl_image *image,
                                            const int       window_xxyy[4],
                                            double          sigma)
{
    int         kernelsize=0, n, center;
    cpl_matrix  *kernelh = NULL,
                *kernelv = NULL;
    cpl_image   *output = NULL,
                *temp_im = NULL;
    double      *kdata, w;

    CLIPM_TRY
    {
        int         windowbuffer[4];
        int         *localwindow = NULL;
        int         imsize[2];
        int         halfsize,
                    original_lower = 0,
                    original_upper = 0;
        cpl_type    intype;
        
        /* check input */
        clipm_priv_checks_window_image(     window_xxyy,
                                            image,
                                            1,
                                            imsize,
                                            NULL, NULL);
        CLIPM_TRY_CHECK_ERROR_STATE();
        clipm_priv_checks_imtype_any(       image, &intype);
        CLIPM_TRY_CHECK_ERROR_STATE();
        
        CLIPM_TRY_CHECK_AUTOMSG(            sigma > 0,
                                            CPL_ERROR_ILLEGAL_INPUT);
        
        halfsize = ceil(4*sigma);
        kernelsize = 2*halfsize + 1;
        center = halfsize;
        
        kernelh = cpl_matrix_new(1, kernelsize);
        CLIPM_TRY_ASSERT_ERROR_STATE();
        
        kdata = cpl_matrix_get_data(kernelh);
        kernelv = cpl_matrix_wrap(kernelsize, 1, kdata);
        CLIPM_TRY_ASSERT_ERROR_STATE();

        w = 1.0/(sigma * CPL_MATH_SQRT2PI);
        for (n = 0; n < kernelsize; n++)
            kdata[n] = w * exp(- ((double)(n-center)/sigma)*
                                            ((double)(n-center)/sigma) / 2);
        
        /* enlarge window vertically, so that the surrounding signal is also
         * useable for the second convolution */
        if (window_xxyy != NULL)
        {
            localwindow = windowbuffer;
            memcpy(localwindow, window_xxyy, 4*sizeof(int));
            clipm_priv_checks_window_guarantee( /* ensure min/max order */
                                            localwindow,
                                            imsize[0],
                                            imsize[1],
                                            1);
            original_lower = localwindow[2];
            original_upper = localwindow[3];
            localwindow[2] -= halfsize;     /* extend down */
            localwindow[3] += halfsize;     /* extend up */
            clipm_priv_checks_window_guarantee( /* ensure fit image (cut) */
                                            localwindow,
                                            imsize[0],
                                            imsize[1],
                                            1);
        }
        else
            localwindow = NULL;
        
        /* convolve horizontally */
        temp_im = clipm_priv_image_conv_matrix(
                                            image,
                                            kernelh,
                                            localwindow,
                                            0,  /* no bad pixel dilation */
                                            1); /* int 2 double */
        CLIPM_TRY_ASSERT_ERROR_STATE();
        
        /* provide a window which cuts the lower/upper extended stripes */
        if (localwindow != NULL)
        {
            localwindow[1] -= localwindow[0] - 1;   /* right */
            localwindow[0] = 1; /* left */
            localwindow[2] = original_lower - localwindow[2] + 1;
            localwindow[3] = original_upper - original_lower + localwindow[2];
        }
        
        output = clipm_priv_image_conv_matrix(
                                            temp_im,
                                            kernelv,
                                            localwindow,
                                            0,  /* no bad pixel dilation */
                                            1); /* don't get int here anyway */
        CLIPM_TRY_ASSERT_ERROR_STATE();
        
        if (intype == CPL_TYPE_INT)
        {
            cpl_image_delete(temp_im);
            temp_im = clipm_priv_image_extract_round(output, NULL);
            CLIPM_TRY_ASSERT_ERROR_STATE();
            
            cpl_image_delete(output);
            output = temp_im;
            temp_im = NULL;
        }
    }
    CLIPM_CATCH
    {
        clipm_priv_image_null(&output);
    }

    cpl_matrix_unwrap(kernelv);
    cpl_matrix_delete(kernelh);
    cpl_image_delete(temp_im);
    
    return output;
}

/**@}*/
