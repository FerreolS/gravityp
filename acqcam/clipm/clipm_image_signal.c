/*********************************************************************
 * E.S.O. - VLT project
 *
 * "@(#) $Id: clipm_image_signal.c 253013 2014-03-19 17:18:08Z cgarcia $"
 *
 * Functions for image signal processing
 *
 * who       when        what
 * --------  ----------  ----------------------------------------------
 * hlorch    2008-08-28  created
 */

/**
 * @defgroup clipm_image_signal    Image Signal Processing
 *
 * This module provides image signal processing functions.
 *
 * @par Synopsis:
 * @code
#include "clipm_image_signal.h"
 * @endcode
 */
/**@{*/

/*-----------------------------------------------------------------------------
    Includes
 -----------------------------------------------------------------------------*/

#include "clipm_image_signal.h"

#include "clipm_math.h"

#include "clipm_priv_image_signal.h"
#include "clipm_priv_checks.h"
#include "clipm_compatibility_replacements.h"
#include "clipm_priv_array.h"
#include "clipm_priv_error.h"
#include "clipm_priv_image.h"
#include "clipm_priv_math.h"

#include <float.h>
#include <limits.h>

/*-----------------------------------------------------------------------------
    Defines
 -----------------------------------------------------------------------------*/


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
double          clipm_image_signal_estimate_bg_in_region(
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
                                            img,
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
double          clipm_image_signal_estimate_fwhm_round(
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
cpl_error_code  clipm_image_signal_get_barycentre(
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

/**@}*/
