/*********************************************************************
 * E.S.O. - VLT project
 *
 * "@(#) $Id: clipm_priv_array.c 177213 2008-12-05 15:33:51Z hlorch $"
 *
 * Private functions for handling arrays
 *
 * who       when        what
 * --------  ----------  ----------------------------------------------
 * hlorch    2007-04-27  created
 */

/**
 * @internal
 * @defgroup clipm_priv_array Array Handling
 * @ingroup internal_docs
 *
 * This module provides private array handling functions.
 *
 * @par Synopsis:
 * @code
#include "clipm_priv_vector.h"
 * @endcode
 */
/**@{*/

/*-----------------------------------------------------------------------------
    Includes
 -----------------------------------------------------------------------------*/

#include "clipm_priv_array.h"

#include "clipm_math.h"
#include "clipm_compatibility_replacements.h"
#include "clipm_priv_checks.h"
#include "clipm_priv_error.h"
#include "clipm_priv_image.h"

#include <string.h> /* memcpy() */

/*-----------------------------------------------------------------------------
    Private Prototypes
 -----------------------------------------------------------------------------*/

static
cpl_error_code  _clipm_priv_array_find_limits(
                                            const double    *data,
                                            const int       *badflags,
                                            int             length,
                                            int             startpos,
                                            double          limit,
                                            double          *left,
                                            double          *right);

/*-----------------------------------------------------------------------------
    Private Implementation
 -----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
 * @brief   In an array, find the left and right positions falling under limit.
 * @param   data        Data array
 * @param   badflags    (Optional) array containing bad flags
 * @param   length      Array size
 * @param   startpos    Starting position
 * @param   limit       Value limit
 * @param   left        Output left position
 * @param   right       Output right position
 * @return  CPL error code 
 * 
 * @par Error Handling:
 * - In the case of error, all output values are set to -1.
 * - The following error codes can be set by this function:
 *   - CPL_ERROR_CONTINUE: the FWHM could not be found, or the peak is not
 *     brighter than @a bg_level
 *   - CPL_ERROR_DATA_NOT_FOUND: too many bad pixels prevent the computation of
 *     the FWHM
 */
/*----------------------------------------------------------------------------*/
static
cpl_error_code  _clipm_priv_array_find_limits(
                                            const double    *data,
                                            const int       *badflags,
                                            int             length,
                                            int             startpos,
                                            double          limit,
                                            double          *left,
                                            double          *right)
{
    CLIPM_TRY
    {
        double          lpos,
                        rpos;
        int             id,
                        iu,
                        last_id,
                        last_iu,
                        goodcountl,
                        goodcountr;

        /* find limit indices */
        last_id = last_iu = -1;
        goodcountl = goodcountr = 0;
        for (id = startpos; id >= 0; id--)
            if (badflags == NULL || !badflags[id])
            {
                goodcountl++;
                if (data[id] <= limit)
                    break;
                last_id = id;
            }
        for (iu = startpos; iu < length; iu++)
            if (badflags == NULL || !badflags[iu])
            {
                goodcountr++;
                if (data[iu] <= limit)
                    break;
                last_iu = iu;
            }
        CLIPM_TRY_CHECK(                    goodcountl >= 2 &&
                                            goodcountr >= 2,
                                            CPL_ERROR_DATA_NOT_FOUND,
                                            "",
                                            "too many bad pixels");
        CLIPM_TRY_CHECK(                    id >= 0 &&
                                            iu < length,
                                            CPL_ERROR_CONTINUE,
                                            "",
                                            "fwhm could not be found");
        
        /* compute left and right limit positions */
        if (data[id]-data[last_id] != 0)
            lpos =  (limit-data[last_id])
                    * ((double)(id-last_id))
                    / (data[id]-data[last_id])
                    + last_id;
        else
            lpos = ((double)(last_id + id)) / 2.0;
        
        if (data[iu]-data[last_iu] != 0)
            rpos =  (limit-data[last_iu])
                    * ((double)(iu-last_iu))
                    / (data[iu]-data[last_iu])
                    + last_iu;
        else
            rpos = ((double)(last_iu + iu)) / 2.0;
        
        *left = lpos;
        *right = rpos;
    }
    CLIPM_CATCH
    {
        *left = -1.0;
        *right = -1.0;
    }
    
    return CLIPM_ERROR_GET_NEW_SINCE_TRY();
}

/*-----------------------------------------------------------------------------
    Implementation
 -----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
 * @brief   Extract an image row and store it in an array.
 * @param   image   Input image
 * @param   row_ndx Row index (FITS convention)
 * @return  The new array
 * 
 * @par Bad Pixel Handling:
 * - Bad pixels will lead to invalid entried in the output array.
 * 
 * @par Error Handling:
 * The following errors can occur:
 * - CPL_ERROR_NULL_INPUT: @a image is NULL
 * - CPL_ERROR_INVALID_TYPE: @a image is not of type CPL_TYPE_INT,
 *   CPL_TYPE_FLOAT or CPL_TYPE_DOUBLE
 * - CPL_ERROR_ACCESS_OUT_OF_RANGE: @a row_ndx is outside the vertical
 *   image range
 */
/*----------------------------------------------------------------------------*/
cpl_array       *clipm_priv_array_new_from_image_row(
                                            const cpl_image *image,
                                            int             row_ndx)
{
    cpl_array   *A = NULL;
    
    CLIPM_TRY
    {
        int         imsize[2],
                    wdwsize[2];
        int         window[2][2];
        cpl_type    type;
        const void  *data_start = NULL;
        const cpl_binary
                    *badp_start = NULL;
        
        CLIPM_TRY_CHECK_AUTOMSG(            image != NULL,
                                            CPL_ERROR_NULL_INPUT);
        clipm_priv_checks_imtype_any(       image, NULL);
        CLIPM_TRY_CHECK_ERROR_STATE();
        
        window[0][0] = 1;
        window[0][1] = cpl_image_get_size_x(image);
        window[1][0] = row_ndx;
        window[1][1] = row_ndx;
        
        clipm_priv_image_get_data_const(    image,
                                            *window,
                                            1,
                                            imsize,
                                            wdwsize,
                                            NULL,
                                            &type,
                                            &data_start,
                                            &badp_start);
        CLIPM_TRY_CHECK_ERROR_STATE();
        
        A = cpl_array_new(wdwsize[0], type);
        CLIPM_TRY_ASSERT_ERROR_STATE();
        
        switch (type)
        {
            case CPL_TYPE_INT:
                cpl_array_copy_data_int(A, data_start);
                break;
            case CPL_TYPE_FLOAT:
                cpl_array_copy_data_float(A, data_start);
                break;
            case CPL_TYPE_DOUBLE:
                cpl_array_copy_data_double(A, data_start);
                break;
            default:
                CLIPM_ERROR_SET_MSG(        CPL_ERROR_INVALID_TYPE,
                                            "image",
                                            "must be int, float or double");
        }
        CLIPM_TRY_ASSERT_ERROR_STATE();
        
        if (badp_start != NULL)
        {
            int x;
            for (x = 0; x < wdwsize[0]; x++)
                if (badp_start[x] != CPL_BINARY_0)
                    cpl_array_set_invalid(A, x);
        }
        CLIPM_TRY_ASSERT_ERROR_STATE();
    }
    CLIPM_CATCH
    {
        clipm_priv_array_null(&A);
    }
    
    return A;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief   Extract an image column and store it in an array.
 * @param   image   Input image
 * @param   col_ndx Column index (FITS convention)
 * @return  The new array
 * 
 * @par Bad Pixel Handling:
 * - Bad pixels will lead to invalid entried in the output array.
 * 
 * @par Error Handling:
 * The following errors can occur:
 * - CPL_ERROR_NULL_INPUT: @a image is NULL
 * - CPL_ERROR_INVALID_TYPE: @a image is not of type CPL_TYPE_INT,
 *   CPL_TYPE_FLOAT or CPL_TYPE_DOUBLE
 * - CPL_ERROR_ACCESS_OUT_OF_RANGE: @a col_ndx is outside the horizontal
 *   image range
 */
/*----------------------------------------------------------------------------*/
cpl_array       *clipm_priv_array_new_from_image_col(
                                            const cpl_image *image,
                                            int             col_ndx)
{
    cpl_array   *A = NULL;
    
    CLIPM_TRY
    {
        int         imsize[2],
                    wdwsize[2];
        int         window[2][2];
        cpl_type    type;
        const void  *data_start = NULL;
        const cpl_binary
                    *badp_start = NULL;
        
        CLIPM_TRY_CHECK_AUTOMSG(            image != NULL,
                                            CPL_ERROR_NULL_INPUT);
        window[0][0] = col_ndx;
        window[0][1] = col_ndx;
        window[1][0] = 1;
        window[1][1] = cpl_image_get_size_y(image);
        
        clipm_priv_image_get_data_const(    image,
                                            *window,
                                            1,
                                            imsize,
                                            wdwsize,
                                            NULL,
                                            &type,
                                            &data_start,
                                            &badp_start);
        CLIPM_TRY_CHECK_ERROR_STATE();
        clipm_priv_checks_imtype_any(       image, NULL);
        CLIPM_TRY_CHECK_ERROR_STATE();
        
        A = cpl_array_new(wdwsize[1], type);
        CLIPM_TRY_ASSERT_ERROR_STATE();
        
        switch (type)
        {
            int y;
            case CPL_TYPE_INT:
            {
                const int       *data = data_start;
                for (y = 0; y < wdwsize[1]; y++, data += imsize[0])
                    cpl_array_set_int(A, y, *data);
                break;
            }
            case CPL_TYPE_FLOAT:
            {
                const float     *data = data_start;
                for (y = 0; y < wdwsize[1]; y++, data += imsize[0])
                    cpl_array_set_float(A, y, *data);
                break;
            }
            case CPL_TYPE_DOUBLE:
            {
                const double    *data = data_start;
                for (y = 0; y < wdwsize[1]; y++, data += imsize[0])
                    cpl_array_set_double(A, y, *data);
                break;
            }
            default:
                CLIPM_ERROR_SET_MSG(        CPL_ERROR_INVALID_TYPE,
                                            "image",
                                            "must be int, float or double");
        }
        CLIPM_TRY_ASSERT_ERROR_STATE();
        
        if (badp_start != NULL)
        {
            int y;
            for (y = 0; y < wdwsize[1]; y++, badp_start += imsize[0])
                if (*badp_start != CPL_BINARY_0)
                    cpl_array_set_invalid(A, y);
        }
        CLIPM_TRY_ASSERT_ERROR_STATE();
    }
    CLIPM_CATCH
    {
        clipm_priv_array_null(&A);
    }
    
    return A;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief   Search maximum from given centre position and find FWHM.
 * @param   input           Input vector
 * @param   centre          Centre position
 * @param   bg_level        Background level
 * @param   out_maxindex    (Optional output) index of found maximum value,
 *                          can be NULL
 * @param   out_middlepos   (Optional output) middle position between
 *                          left and right FWHM edges,
 *                          can be NULL
 * @param   out_edgeslope   (Optional output) steepness of the edges,
 *                          returns -1.0 if unsuccessful,
 *                          can be NULL
 * @return  Full-width-half-maximum, -1.0 in the case of error
 * 
 * @par Principle:
 * - The maximum is searched in ascending direction.
 * - From the maximum's position, the borders where the signal falls below
 *   the mean of the maximum and @a bg_level are searched.
 * - If @a out_edgeslope is given, then the steepness is measured between
 *   the 30% and 70% percent levels. This is done for both end, if one fails,
 *   then the other value is returned, otherwise the mean of both.
 * 
 * @note
 * - @a out_maxindex does not return the position of the maximum value of the
 *   whole @a input vector, but of the maximum found by starting at @a centre.
 * - This function jumps over "bad" flagges regions.
 * - The FWHM is interpolated using non-bad values.
 * - If @a out_edgeslope can not be computed, NO error is set
 * 
 * @par Error Handling:
 * - In the case of error, all output values are set to -1.
 * - The following error codes can be set by this function:
 *   - CPL_ERROR_NULL_INPUT: @a input vector is NULL
 *   - CPL_ERROR_ACCESS_OUT_OF_RANGE: @a centre is outside the @a input
 *     vector range
 *   - CPL_ERROR_INVALID_TYPE: @a input is not of type CPL_TYPE_INT,
 *     CPL_TYPE_FLOAT or CPL_TYPE_DOUBLE
 *   - CPL_ERROR_CONTINUE: the FWHM could not be found, or the peak is not
 *     brighter than @a bg_level
 *   - CPL_ERROR_DATA_NOT_FOUND: too many bad pixels prevent the computation of
 *     the FWHM
 * 
 * @todo
 * - Force _clipm_priv_vector_find_limits() to fail in unit test
 * - Check out_edgeslope in unit test, also with flat object
 */
/*----------------------------------------------------------------------------*/
double          clipm_priv_array_estimate_fwhm(
                                            const cpl_array     *input,
                                            double              centre,
                                            double              bg_level,
                                            int                 *out_maxindex,
                                            double              *out_middlepos,
                                            double              *out_edgeslope)
{
    double      fwhm = -1.0;
    double      *data = NULL;
    int         *bad = NULL;
    
    CLIPM_TRY
    {
        double          lpos,
                        rpos,
                        max,
                        halfmax;
        int             length,
                        maxpos,
                        id,
                        iu,
                        goodcount,
                        n;
        cpl_type        type;
        
        CLIPM_TRY_CHECK_AUTOMSG(            input != NULL,
                                            CPL_ERROR_NULL_INPUT);
        length = cpl_array_get_size(input);
        
        CLIPM_TRY_CHECK_AUTOMSG(            centre >= 0,
                                            CPL_ERROR_ACCESS_OUT_OF_RANGE);
        CLIPM_TRY_CHECK_AUTOMSG(            centre <= length - 1.0,
                                            CPL_ERROR_ACCESS_OUT_OF_RANGE);
        
        type = cpl_array_get_type(input);
        CLIPM_TRY_CHECK(                    type == CPL_TYPE_INT
                                            || type == CPL_TYPE_FLOAT
                                            || type == CPL_TYPE_DOUBLE,
                                            CPL_ERROR_INVALID_TYPE,
                                            "input",
                                            "must be int, float or double");
        
        data = cpl_malloc(length * sizeof(*data));
        bad  = cpl_malloc(length * sizeof(*bad));
        CLIPM_TRY_CHECK_ERROR_STATE();
        
        /* copy input data into "data" and valid flags into "bad" */
        for (n = 0; n < length; n++)
            data[n] = cpl_array_get(input, n, bad + n);
        CLIPM_TRY_ASSERT_ERROR_STATE();
        
        /* find real maximum */
        maxpos = clipm_math_round_d2i(centre);
        max = halfmax = bg_level;
        goodcount = 0;
        for (id = iu = maxpos; id >= 0 || iu < length; id--, iu++)
        {
            if (id >= 0 && !bad[id])
            {
                goodcount++;
                if (data[id] > max)
                {
                    max = data[id];
                    halfmax = (max + bg_level) / 2.0;
                    maxpos = id;
                }
            }
            if (iu < length && !bad[iu])
            {
                goodcount++;
                if (data[iu] > max)
                {
                    max = data[iu];
                    halfmax = (max + bg_level) / 2.0;
                    maxpos = iu;
                }
                else if (data[iu] < halfmax)
                    iu = length;    /* break searching upwards */
            }
            if (    id >= 0 && !bad[id]
                    && data[id] < halfmax)
                id  = -1;           /* break searching downwards */
        }
        
        CLIPM_TRY_CHECK(                    goodcount >= 1,
                                            CPL_ERROR_DATA_NOT_FOUND,
                                            "",
                                            "too many bad pixels");
        /* the check below also prevents negative fwhm */
        CLIPM_TRY_CHECK(                    max > bg_level,
                                            CPL_ERROR_CONTINUE,
                                            "",
                                            "peak is not brighter"
                                            " than background");
        /* we cannot find the fwhm if the peak is at the border */
        CLIPM_TRY_CHECK(                    maxpos > 0 && maxpos < length -1,
                                            CPL_ERROR_CONTINUE,
                                            "",
                                            "fwhm could not be found");
        
        _clipm_priv_array_find_limits(      data,
                                            bad,
                                            length,
                                            maxpos,
                                            halfmax,
                                            &lpos,
                                            &rpos);
        CLIPM_TRY_CHECK_ERROR_STATE();
        
        fwhm = rpos - lpos;
        
        CLIPM_TRY_ASSERT(                   fwhm >= 0.0);
        
        if (out_maxindex != NULL)
            *out_maxindex = maxpos;
        if (out_middlepos != NULL)
            *out_middlepos = (lpos + rpos) / 2.0;
        
        if (out_edgeslope != NULL)
        {
            double  lup,
                    rup,
                    ldown,
                    rdown;
            cpl_error_code  errc;
            
            errc = _clipm_priv_array_find_limits(
                                            data,
                                            bad,
                                            length,
                                            maxpos,
                                            bg_level + (max-bg_level)*0.7,
                                            &lup,
                                            &rup);
            if (errc == CPL_ERROR_NONE)
                errc = _clipm_priv_array_find_limits(
                                            data,
                                            bad,
                                            length,
                                            maxpos,
                                            bg_level + (max-bg_level)*0.3,
                                            &ldown,
                                            &rdown);
            if (errc == CPL_ERROR_NONE)
            {
                double  lslope,
                        rslope,
                        ldev,
                        rdev;
                
                lslope = (max-bg_level)*0.4 / fabs(lup - ldown);
                rslope = (max-bg_level)*0.4 / fabs(rup - rdown);
                
                /* check both slopes, and if necessary choose the best,
                 * check if the fwhm point is out of the center of the slope */
                ldev = fabs((lup + ldown)/2.0 - lpos);
                rdev = fabs((rup + rdown)/2.0 - rpos);
                if (ldev > fabs(lup + ldown) / 2.0 / 1.5)
                    lslope = -1;
                if (rdev > fabs(rup + rdown) / 2.0 / 1.5)
                    rslope = -1;
                
                if (lslope >= 0 && rslope >= 0)
                    *out_edgeslope = (lslope + rslope) / 2.0;
                else if (lslope >= 0 && rslope < 0)
                    *out_edgeslope = lslope;
                else if (lslope < 0 && rslope >= 0)
                    *out_edgeslope = rslope;
                else 
                    *out_edgeslope = -1.0;
            }
            else
            {
                CLIPM_ERROR_RECOVER_TRYSTATE();
                *out_edgeslope = -1.0;
            }
        }
    }
    CLIPM_CATCH
    {
        fwhm = -1.0;
        if (out_maxindex != NULL)
            *out_maxindex = -1;
        if (out_middlepos != NULL)
            *out_middlepos = -1.0;
        if (out_edgeslope != NULL)
            *out_edgeslope = -1.0;
    }
    
    if (data != NULL)
        cpl_free(data);
    if (bad != NULL)
        cpl_free(bad);
    
    return fwhm;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief   Delete a CPL array object and set the pointer to NULL.
 * @param   a   Pointer to array pointer
 * @return  Nothing
 * 
 * The following code is executed:
 * @code
    if (a != NULL)
    {
        cpl_array_delete(*a);  // checks for NULL pointer
        *a = NULL;
    }
 * @endcode
 * 
 * @par Error Handling:
 * No error can occur here.
 */
/*----------------------------------------------------------------------------*/
void            clipm_priv_array_null(      cpl_array       **a)
{
    if (a != NULL)
    {
        cpl_array_delete(*a);  /* checks for NULL pointer */
        *a = NULL;
    }
}

/**@}*/
