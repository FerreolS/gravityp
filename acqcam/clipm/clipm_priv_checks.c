/*********************************************************************
 * E.S.O. - VLT project
 *
 * "@(#) $Id: clipm_priv_checks.c 177213 2008-12-05 15:33:51Z hlorch $"
 *
 * Private functions for checking parameters
 *
 * who       when        what
 * --------  ----------  ----------------------------------------------
 * hlorch    2006-08-12  created
 */

/**
 * @internal
 * @defgroup clipm_priv_checks Parameter Checks
 * @ingroup internal_docs
 *
 * This module provides parameter checks.
 *
 * @par Synopsis:
 * @code
#include "clipm_priv_checks.h"
 * @endcode
 */
/**@{*/

/*-----------------------------------------------------------------------------
    Includes
 -----------------------------------------------------------------------------*/

#include "clipm_priv_checks.h"

#include "clipm_priv_error.h"
#include "clipm_priv_math.h"

/*-----------------------------------------------------------------------------
    Defines
 -----------------------------------------------------------------------------*/

#define max(a, b) ((a) > (b) ? (a) : (b))
#define min(a, b) ((a) < (b) ? (a) : (b))

/*-----------------------------------------------------------------------------
    Implementation
 -----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
 * @brief   Check if window coordinates represent the whole image
 * @param   window_xxyy         Coordinate buffer of the form
 *                              {x1a, x1b, y1a, y1b},
 *                              minimum/maximum order is irrelevant
 * @param   image               Image (FITS convention)
 * @return  CPL error code
 * 
 * @par Principle:
 * If @a window_xxyy is != NULL and the coordinates in @a window_xxyy
 * represent exactly the whole image plane, then 1 is returned.
 * Otherwise 0 is returned.
 * 
 * @par Error Handling:
 * The following codes are set, if the check fails:
 * - CPL_ERROR_NULL_INPUT: @a image is NULL
 * .
 * No error is set if @a window_xxyy is NULL or if the coordinates in
 * @a window_xxyy are outside the image range.
 */
/*----------------------------------------------------------------------------*/
int             clipm_priv_checks_is_window_full_image(
                                            const int       window_xxyy[4],
                                            const cpl_image *image)
{
    int result = 0;
    
    CLIPM_TRY
    {
        int xsize,
            ysize;
        
        xsize = cpl_image_get_size_x(image);
        ysize = cpl_image_get_size_y(image);
        CLIPM_TRY_CHECK_ERROR_STATE();
        
        if (    window_xxyy != NULL &&
                (   (   window_xxyy[0] == 1 &&
                        window_xxyy[1] == xsize) ||
                    (   window_xxyy[0] == xsize &&
                        window_xxyy[1] == 1)
                ) &&
                (   (   window_xxyy[2] == 1 &&
                        window_xxyy[3] == ysize) ||
                    (   window_xxyy[2] == ysize &&
                        window_xxyy[3] == 1)
                )
            )
        {
            result = 1;
        }
    }
    CLIPM_CATCH
    {
    }
    
    return result;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief   Check window coordinates with the corresponding image
 * @param   window_xxyy         Coordinate buffer of the form
 *                              {x1a, x1b, y1a, y1b},
 *                              minimum/maximum order is irrelevant
 * @param   image               Image (FITS convention)
 * @param   allow_window_NULL   If !0, then no error is set if (window==NULL)
 * @param   img_size_xy         Buffer of size 2, to which the image size [x, y]
 *                              is put out, can be NULL
 * @param   window_size_xy      Buffer of size 2, to which the window size
 *                              is determined and written, can be NULL
 * @param   buffer_start_xy     Buffer of size 2, to which the lower left window
 *                              indices are written, starting at 0 (not FITS),
 *                              can be NULL
 * @return  CPL error code
 * 
 * @par Overview:
 * - The window coordinates are checked for being inside the image range.
 * - If an input is NULL, or if the coordinates are out of range, then the
 *   appropriate CPL error is set and returned.
 * - If @a out_window_size != NULL, then the window size is calculated and
 *   written to this buffer. If @a window_xxyy == NULL, then @a out_window_size
 *   is set to the image size!
 * 
 * @par Error Handling:
 * The following codes are set, if the check fails:
 * - CPL_ERROR_NULL_INPUT: @a image and/or @a window_xxyy is NULL
 * - CPL_ERROR_ACCESS_OUT_OF_RANGE: window coordinates are outside the image
 * In the case of error, the output parameter values are undefined.
 * 
 * @note
 * - This function does not check min/max conditions. Use
 *   clipm_priv_checks_window_minmax() for this purpose.
 * - Specifying buffers of too small size
 *   (sizeof(window) < 4 || sizeof(out_img_size) < 2 ||
 *   sizeof(out_window_size) < 2 || sizeof(out_buffer_start_xy) < 2)
 *   will crash the application!
 */
/*----------------------------------------------------------------------------*/
cpl_error_code  clipm_priv_checks_window_image(
                                            const int       window_xxyy[4],
                                            const cpl_image *image,
                                            int             allow_window_NULL,
                                            int             *img_size_xy,
                                            int             *window_size_xy,
                                            int             *buffer_start_xy)
{
    CLIPM_TRY
    {
        int n,
            d;
        int temp_size[2];
        if (img_size_xy == NULL)
            img_size_xy = (int*)temp_size;
        
        CLIPM_TRY_CHECK_AUTOMSG(            image != NULL,
                                            CPL_ERROR_NULL_INPUT);
        CLIPM_TRY_CHECK_AUTOMSG(            allow_window_NULL
                                                || window_xxyy != NULL,
                                            CPL_ERROR_NULL_INPUT);
        
        img_size_xy[0] = cpl_image_get_size_x(image);
        img_size_xy[1] = cpl_image_get_size_y(image);
        CLIPM_TRY_CHECK_ERROR_STATE(); /* complains if image == NULL */
        
        if (window_xxyy != NULL) {
            for (d = 0; d < 2; d++)
            {
                for (n = 0; n < 2; n++)
                {
                    if (!(window_xxyy[2*d+n] >= 1))
                    {
                        char    msg[50];
                        sprintf(            msg,
                                            "!((window_xxyy[%d]=%d) >= 1)",
                                            2*d+n,
                                            window_xxyy[2*d+n]);
                        CLIPM_TRY_EXIT_WITH_ERROR_MSG(
                                            CPL_ERROR_ACCESS_OUT_OF_RANGE,
                                            "",
                                            msg);
                    }
                    if (!(window_xxyy[2*d+n] <= img_size_xy[d]))
                    {
                        char    msg[50];
                        sprintf(            msg,
                                            "!((window_xxyy[%d]=%d) <="
                                                " (img_size_xy[%d]=%d))",
                                            2*d+n,
                                            window_xxyy[2*d+n],
                                            d,
                                            img_size_xy[d]);
                        CLIPM_TRY_EXIT_WITH_ERROR_MSG(
                                            CPL_ERROR_ACCESS_OUT_OF_RANGE,
                                            "",
                                            msg);
                    }
                    /*CLIPM_TRY_CHECK_AUTOMSG(window_xxyy[2*d+n] >= 1,
                                            CPL_ERROR_ACCESS_OUT_OF_RANGE);
                    CLIPM_TRY_CHECK_AUTOMSG(window_xxyy[2*d+n]
                                                <= img_size_xy[d],
                                            CPL_ERROR_ACCESS_OUT_OF_RANGE);*/
                }
            }
            if (window_size_xy != NULL) {
                window_size_xy[0] = abs(window_xxyy[1] - window_xxyy[0]) +1;
                window_size_xy[1] = abs(window_xxyy[3] - window_xxyy[2]) +1;
            }
            if (buffer_start_xy != NULL) {
                buffer_start_xy[0] = min(window_xxyy[0], window_xxyy[1]) -1;
                buffer_start_xy[1] = min(window_xxyy[2], window_xxyy[3]) -1;
            }
        } else {
            if (window_size_xy != NULL) {
                window_size_xy[0] = img_size_xy[0];
                window_size_xy[1] = img_size_xy[1];
            }
            if (buffer_start_xy != NULL) {
                buffer_start_xy[0] = 0;
                buffer_start_xy[1] = 0;
            }
        }
    }
    CLIPM_CATCH
    {
    }
    
    return CLIPM_ERROR_GET_NEW_SINCE_TRY();
}

/*----------------------------------------------------------------------------*/
/**
 * @brief   Check coordinates for min/max condition.
 * @param   window_xxyy Coordinate buffer containing [xmin, xmax, ymin, ymax,
 *                      zmin, zmax,...], its size must be 2*ndims
 * @param   ndims       Number of dimensions
 * @param   allow_window_NULL   If > 0, then no error is set if (window==NULL)
 * @return  CPL error code
 * 
 * @par Overview:
 * The min coordinates are checked for being <= the max coordinates. If they
 * violate this, then CPL_ERROR_ACCESS_OUT_OF_RANGE is set and returned.
 * 
 * @par Error Handling:
 * The following codes are set, if the check fails:
 * - CPL_ERROR_NULL_INPUT: @a image is NULL
 * 
 * @note
 * Specifying a buffer of too small size
 * (sizeof(window) < 4) will crash the application!
 */
/*----------------------------------------------------------------------------*/
cpl_error_code  clipm_priv_checks_window_minmax(
                                            const int       window_xxyy[4],
                                            int             ndims,
                                            int             allow_window_NULL)
{
    int dim;

    CLIPM_TRY
    {
        CLIPM_TRY_CHECK_AUTOMSG(
            allow_window_NULL || window_xxyy != NULL,
            CPL_ERROR_NULL_INPUT);
        
        if (window_xxyy != NULL)
            for (dim = 0; dim < ndims; dim++)
                CLIPM_TRY_CHECK_AUTOMSG(
                    window_xxyy[2*dim] <= window_xxyy[2*dim+1],
                    CPL_ERROR_ACCESS_OUT_OF_RANGE);
    }
    CLIPM_CATCH
    {
    }
    
    return CLIPM_ERROR_GET_NEW_SINCE_TRY();
}

/*----------------------------------------------------------------------------*/
/**
 * @brief   Cut window coordinates to a FITS image range if outside, expand
 *          if necessary.
 * @param   window_xxyy Coordinate buffer of the form
 *                      {x1a, x1b, y1a, y1b} in FITS convention,
 *                      minimum/maximum order is irrelevant
 * @param   xsize       Horizontal image size
 * @param   ysize       Vertical image size
 * @param   min_windowsize   Minimum window size
 * @return  CPL error code
 * 
 * @par Principle:
 * This function cuts window coordinates so that they fit into a FITS image
 * range:
 * -# The window entries are ordered to {x_min, y_min, x_max, y_max}.
 * -# If @a window is already smaller than @a min_windowsize, then it is
 *    expanded in the respective dimension.
 * -# If @a window exceeds the FITS image range, it is cut.
 * -# If, by the cutting, @a window gets smaller than @a min_windowsize in x or
 *    y, it is expanded into the opposite direction. This is to ensure a minimum
 *    window size for certain purposes, e.g. like gaussian fitting etc.
 * 
 * @par Error Handling:
 * The following error codes can be set and returned:
 * - CPL_ERROR_NULL_INPUT: @a window is NULL
 * - CPL_ERROR_ACCESS_OUT_OF_RANGE: @a min_windowsize > @a xsize or @a ysize
 * - CPL_ERROR_ILLEGAL_INPUT: @a min_windowsize < 1
 * .
 * In the case of error, @a window is not modified.
 */
/*----------------------------------------------------------------------------*/
cpl_error_code  clipm_priv_checks_window_guarantee(
                                            int             window_xxyy[4],
                                            int             xsize,
                                            int             ysize,
                                            int             min_windowsize)
{
    CLIPM_TRY
    {
        int dim;
        int *wp;
        int xysize[2];
        
        if (window_xxyy == NULL)
            CLIPM_TRY_EXIT_WITH_ERROR(      CPL_ERROR_NULL_INPUT);

        CLIPM_TRY_CHECK_AUTOMSG(
            min_windowsize <= xsize,
            CPL_ERROR_ACCESS_OUT_OF_RANGE);
        CLIPM_TRY_CHECK_AUTOMSG(
            min_windowsize <= ysize,
            CPL_ERROR_ACCESS_OUT_OF_RANGE);
        
        CLIPM_TRY_CHECK_AUTOMSG(
            min_windowsize >= 1,
            CPL_ERROR_ILLEGAL_INPUT);
        
        xysize[0] = xsize;
        xysize[1] = ysize;
        
        for (dim = 0; dim < 2; dim++)
        {
            int temp,
                s;
            
            wp = window_xxyy + 2*dim;
            
            /* ensure min/max order */
            if (wp[1] < wp[0])
            {
                temp = wp[1];
                wp[1] = wp[0];
                wp[0] = temp;
            }

            /* expand if necessary */
            if ((s = wp[1] - wp[0] + 1) < min_windowsize)
            {
                wp[0] -= (min_windowsize - s) / 2;
                wp[1] = wp[0] + min_windowsize - 1;
            }

            /* cut to FITS range */
            if (wp[0] < 1)
            {
                wp[0] = 1;
                /* check min_windowsize */
                if (wp[1] < min_windowsize)
                    wp[1] = min_windowsize;
            }
            if (wp[1] > xysize[dim])
            {
                wp[1] = xysize[dim];
                if (wp[0] > (temp = xysize[dim] - min_windowsize + 1))
                    wp[0] = temp;
            }
        }
        
    }
    CLIPM_CATCH
    {
    }
    
    return CLIPM_ERROR_GET_NEW_SINCE_TRY();
}

/*----------------------------------------------------------------------------*/
/**
 * @brief   Cut window coordinates to a FITS image range if outside, expand
 *          if necessary.
 * @param   window_xxyy Coordinate buffer of the form
 *                      {x1a, x1b, y1a, y1b} in FITS convention,
 *                      minimum/maximum order is irrelevant
 * @param   image       The reference image
 * @param   min_windowsize   Minimum window size
 * @return  CPL error code
 * 
 * @see
 * For the complete documentation, please refer to
 * clipm_priv_checks_window_guarantee(). If @a image is NULL, then also
 * CPL_ERROR_NULL_INPUT is set and returned.
 */
/*----------------------------------------------------------------------------*/
cpl_error_code  clipm_priv_checks_window_guarantee_image(
                                            int             window_xxyy[4],
                                            const cpl_image *image,
                                            int             min_windowsize)
{
    CLIPM_TRY
    {
        int xsize,
            ysize;
        
        xsize = cpl_image_get_size_x(       image);
        ysize = cpl_image_get_size_y(       image);
        CLIPM_TRY_CHECK_ERROR_STATE();
        
        clipm_priv_checks_window_guarantee( window_xxyy,
                                            xsize,
                                            ysize,
                                            min_windowsize);
    }
    CLIPM_CATCH
    {
    }
    
    return CLIPM_ERROR_GET_NEW_SINCE_TRY();
}

/*----------------------------------------------------------------------------*/
/**
 * @brief   Cut window coordinates to an existing reference window range,
 *          expand if necessary.
 * @param   window_xxyy     Coordinate buffer of the form
 *                          {x1a, x1b, y1a, y1b} in FITS convention,
 *                          minimum/maximum order is irrelevant
 * @param   ref_window      Coordinate buffer of the form
 *                          {x1a, x1b, y1a, y1b} in FITS convention,
 *                          minimum/maximum order is irrelevant
 * @param   min_windowsize  Minimum window size
 * @return  CPL error code
 * 
 * @see
 * For the complete documentation, please refer to
 * clipm_priv_checks_window_guarantee(). If @a image is NULL, then also
 * CPL_ERROR_NULL_INPUT is set and returned.
 */
/*----------------------------------------------------------------------------*/
cpl_error_code  clipm_priv_checks_window_guarantee_window(
                                            int             window_xxyy[4],
                                            const int       *ref_window,
                                            int             min_windowsize)
{
    CLIPM_TRY
    {
        int size[2],
            lld[2];
        int dim;
        
        CLIPM_TRY_CHECK_AUTOMSG(            window_xxyy != NULL,
                                            CPL_ERROR_NULL_INPUT);
        CLIPM_TRY_CHECK_AUTOMSG(            ref_window != NULL,
                                            CPL_ERROR_NULL_INPUT);
        
        for (dim = 0; dim < 2; dim++)
        {
            /* delta from (1,1) to lower left */
            lld[dim]                = clipm_priv_math_min(
                                            ref_window[2*dim],
                                            ref_window[2*dim+1])
                                        - 1;
            size[dim]               = abs(  ref_window[2*dim]
                                            - ref_window[2*dim+1])
                                        + 1;
            window_xxyy[2*dim]      -= lld[dim];
            window_xxyy[2*dim+1]    -= lld[dim];
        }
        
        clipm_priv_checks_window_guarantee( window_xxyy,
                                            size[0],
                                            size[1],
                                            min_windowsize);

        for (dim = 0; dim < 2; dim++)
        {
            window_xxyy[2*dim]      += lld[dim];
            window_xxyy[2*dim+1]    += lld[dim];
        }
    }
    CLIPM_CATCH
    {
    }
    
    return CLIPM_ERROR_GET_NEW_SINCE_TRY();
}

/*----------------------------------------------------------------------------*/
/**
 * @brief   Check whether two images match in the requested features.
 * @param   img1    Image 1
 * @param   img2    Image 2
 * @param   type    (output) CPL type, can be NULL
 * @param   xsize   (output) Horizontal size, can be NULL
 * @param   ysize   (output) Vertical size, can be NULL
 * @return  CPL error code
 * 
 * @par Error Handling:
 * The following error codes can be set:
 * - CPL_ERROR_NULL_INPUT: @a img1 or @a img2 is NULL
 * - CPL_ERROR_TYPE_MISMATCH:
 *   @a type != NULL and @a img1 and @a img2 don't have the same type
 * - CPL_ERROR_INCOMPATIBLE_INPUT:
 *   - @a xsize != NULL and @a img1 and @a img2 differ in their horizontal size,
 *     or
 *   - @a ysize != NULL and @a img1 and @a img2 differ in their vertical size
 */
/*----------------------------------------------------------------------------*/
cpl_error_code  clipm_priv_checks_images_match(
                                            const cpl_image *img1,
                                            const cpl_image *img2,
                                            cpl_type        *type,
                                            int             *xsize,
                                            int             *ysize)
{
    CLIPM_TRY
    {
        if (type)
            *type = CPL_TYPE_INVALID;
        if (xsize)
            *xsize = 0;
        if (ysize)
            *ysize = 0;
        
        if (!img1 || !img2)
            CLIPM_TRY_EXIT_WITH_ERROR(      CPL_ERROR_NULL_INPUT);
        
        if (type)
        {
            *type = cpl_image_get_type(img1);
            if (*type != cpl_image_get_type(img2))
            CLIPM_TRY_EXIT_WITH_ERROR_MSG(  CPL_ERROR_TYPE_MISMATCH,
                                            "images",
                                            CLIPM_MSG_ERR_DIFFTYPES);
        }
        
        if (xsize)
            *xsize = cpl_image_get_size_x(img1);
        if (ysize)
            *ysize = cpl_image_get_size_y(img1);

        if (xsize && ysize)
        {
            if (    *xsize != cpl_image_get_size_x(img2)
                ||  *ysize != cpl_image_get_size_y(img2))
            CLIPM_TRY_EXIT_WITH_ERROR_MSG(  CPL_ERROR_INCOMPATIBLE_INPUT,
                                            "images",
                                            CLIPM_MSG_ERR_DIFFSIZES);
        }
        else if (xsize)
        {
            if (*xsize != cpl_image_get_size_x(img2))
            CLIPM_TRY_EXIT_WITH_ERROR_MSG(  CPL_ERROR_INCOMPATIBLE_INPUT,
                                            "images",
                                            "differ in width");
        }
        else if (ysize)
        {
            if (*ysize != cpl_image_get_size_y(img2))
            CLIPM_TRY_EXIT_WITH_ERROR_MSG(  CPL_ERROR_INCOMPATIBLE_INPUT,
                                            "images",
                                            "differ in height");
        }
    }
    CLIPM_CATCH
    {
        if (type)
            *type = CPL_TYPE_INVALID;
        if (xsize)
            *xsize = 0;
        if (ysize)
            *ysize = 0;
    }
    
    return CLIPM_ERROR_GET_NEW_SINCE_TRY();
}

/*----------------------------------------------------------------------------*/
/**
 * @brief   Check whether an image is of an allowed type.
 * @param   image       Input image
 * @param   out_type    (Optional output) image type
 * @return  CPL error code
 * 
 * @par Error Handling:
 * The following error codes can be returned:
 * - CPL_ERROR_NULL_INPUT: @a image is NULL
 * - CPL_ERROR_INVALID_TYPE: @a image is neither of type CPL_TYPE_INT,
 *   CPL_TYPE_FLOAT, nor CPL_TYPE_DOUBLE
 */
/*----------------------------------------------------------------------------*/
cpl_error_code  clipm_priv_checks_imtype_any(
                                            const cpl_image *image,
                                            cpl_type        *out_type)
{
    CLIPM_TRY
    {
        cpl_type    type;
        
        CLIPM_TRY_CHECK_AUTOMSG(            image != NULL,
                                            CPL_ERROR_NULL_INPUT);
        
        type = cpl_image_get_type(image);
        
        CLIPM_TRY_CHECK(                    type == CPL_TYPE_INT
                                            || type == CPL_TYPE_FLOAT
                                            || type == CPL_TYPE_DOUBLE,
                                            CPL_ERROR_INVALID_TYPE,
                                            "image",
                                            "must be int, float or double");
        
        if (out_type != NULL)
            *out_type = type;
    }
    CLIPM_CATCH
    {
        if (out_type != NULL)
            *out_type = CPL_TYPE_INVALID;
    }
    
    return CLIPM_ERROR_GET_NEW_SINCE_TRY();
}

/*----------------------------------------------------------------------------*/
/**
 * @brief   Check whether an image is of a floating-point type.
 * @param   image       Input image
 * @param   out_type    (Optional output) image type
 * @return  CPL error code
 * 
 * @par Error Handling:
 * The following error codes can be returned:
 * - CPL_ERROR_NULL_INPUT: @a image is NULL
 * - CPL_ERROR_INVALID_TYPE: @a image is neither of type CPL_TYPE_FLOAT
 *   nor CPL_TYPE_DOUBLE
 */
/*----------------------------------------------------------------------------*/
cpl_error_code  clipm_priv_checks_imtype_float(
                                            const cpl_image *image,
                                            cpl_type        *out_type)
{
    CLIPM_TRY
    {
        cpl_type    type;
        
        CLIPM_TRY_CHECK_AUTOMSG(            image != NULL,
                                            CPL_ERROR_NULL_INPUT);
        
        type = cpl_image_get_type(image);
        
        CLIPM_TRY_CHECK(                    type == CPL_TYPE_FLOAT
                                            || type == CPL_TYPE_DOUBLE,
                                            CPL_ERROR_INVALID_TYPE,
                                            "image",
                                            "must be float or double");
        
        if (out_type != NULL)
            *out_type = type;
    }
    CLIPM_CATCH
    {
        if (out_type != NULL)
            *out_type = CPL_TYPE_INVALID;
    }
    
    return CLIPM_ERROR_GET_NEW_SINCE_TRY();
}

/**@}*/
