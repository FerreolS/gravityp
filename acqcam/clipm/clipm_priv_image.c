/*********************************************************************
 * E.S.O. - VLT project
 *
 * "@(#) $Id: clipm_priv_image.c 177716 2008-12-15 15:10:00Z hlorch $"
 *
 * Private functions for handling images
 *
 * who       when        what
 * --------  ----------  ----------------------------------------------
 * hlorch    2007-01-11  created
 */

/**
 * @internal
 * @defgroup clipm_priv_image Image Handling Basics
 * @ingroup internal_image
 *
 * This module provides private image handling functions.
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

#include "clipm_priv_image.h"

#include "clipm_math.h"
#include "clipm_math_correlation.h"

#include "clipm_priv_checks.h"
#include "clipm_compatibility_replacements.h"
#include "clipm_priv_error.h"
#include "clipm_priv_math.h"
#include "clipm_priv_matrix.h"
#include "clipm_priv_vector.h"

/*-----------------------------------------------------------------------------
    Private Prototypes
 -----------------------------------------------------------------------------*/

static
int             _clipm_priv_image_is_point_in_polygon(
                                            int             N,
                                            const double    *xp,
                                            const double    *yp,
                                            double          x,
                                            double          y);

static
int             _clipm_priv_image_segment_cuts_polygon(
                                            int             N,
                                            const double    *xp,
                                            const double    *yp,
                                            double          x1,
                                            double          y1,
                                            double          x2,
                                            double          y2);

static
double          _clipm_priv_image_get_pixel_area_in_polygon(
                                            int             N,
                                            const double    *xp,
                                            const double    *yp,
                                            double          x,
                                            double          y);

static
cpl_mask        *_clipm_priv_image_mark_polygon_edges(
                                            const int       window_xxyy[4],
                                            int             N,
                                            const double    *xp,
                                            const double    *yp);

/*-----------------------------------------------------------------------------
    Private Implementation
 -----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
 * @brief Return whether a point lies inside a polygon.
 * @param   N   Number of polygon points
 * @param   xp  X polygon corner coordinates
 * @param   yp  Y polygon corner coordinates
 * @param   x   X point position
 * @param   y   Y point position
 * @return  1 or 0
 * 
 * @see
 * - http://local.wasp.uwa.edu.au/~pbourke/geometry/insidepoly/
 * 
 * Code by Randolph Franklin.
 * 
 * @note
 * - Specifying NULL pointers will crash the application.
 */
/*----------------------------------------------------------------------------*/
static
int             _clipm_priv_image_is_point_in_polygon(
                                            int             N,
                                            const double    *xp,
                                            const double    *yp,
                                            double          x,
                                            double          y)
{
    int i,
        j,
        c = 0;
    
    for (i = 0, j = N-1; i < N; j = i, i++)
    {
        if ((((yp[i] <= y) && (y < yp[j])) ||
             ((yp[j] <= y) && (y < yp[i]))) &&
            (x < (xp[j] - xp[i]) * (y - yp[i]) / (yp[j] - yp[i]) + xp[i]))
            c = !c;
    }
    return c;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief   Return whether a line segment touches/intersects with at least one
 *          polygon segment.
 * @param   N   Number of polygon points
 * @param   xp  X polygon corner coordinates
 * @param   yp  Y polygon corner coordinates
 * @param   x1  X segment start
 * @param   y1  Y segment start
 * @param   x2  X segment end
 * @param   y2  Y segment end
 * @return  1 or 0
 * 
 * @note
 * - Specifying NULL pointers will crash the application.
 */
/*----------------------------------------------------------------------------*/
static
int             _clipm_priv_image_segment_cuts_polygon(
                                            int             N,
                                            const double    *xp,
                                            const double    *yp,
                                            double          x1,
                                            double          y1,
                                            double          x2,
                                            double          y2)
{
    int     i,
            j,
            c = 0;
    double  x3, y3, x4, y4,
            x43, x31, x21, x13,
            y43, y31, y21, y13,
            s, t, denom;
    /* xA = x1 + s(x2 - x1),  xB = x3 + t(x4 - x3)
     * yA = y1 + s(y2 - y1),  yB = y3 + t(y4 - y3)
     * xA := xB,  yA := yB
     */
    for (i = 0, j = N-1; i < N; j = i, i++)
    {
        x3 = xp[i];
        y3 = yp[i];
        x4 = xp[j];
        y4 = yp[j];
        x43 = x4 - x3;
        x31 = x3 - x1;
        x21 = x2 - x1;
        x13 = x1 - x3;
        y43 = y4 - y3;
        y31 = y3 - y1;
        y21 = y2 - y1;
        y13 = y1 - y3;
        denom = x43*y21 - y43*x21;
        if (fabs(N) > 1e-9)
        {
            t = (y31*x21 - x31*y21) / denom;
            s = (y13*x43 - x13*y43) / (-denom);
            c |= (s >= 0.0 && s <= 1.0 && t >= 0.0 && t <= 1.0);
        }
    }
    return c;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief   Return pixel fraction which is inside a polygon.
 * @param   N   Number of polygon points
 * @param   xp  X polygon corner coordinates
 * @param   yp  Y polygon corner coordinates
 * @param   x   X point position
 * @param   y   Y point position
 * @return  Value in the range [0...1], quantization step is 0.01;
 * 
 * @note
 * - Specifying NULL pointers will crash the application.
 */
/*----------------------------------------------------------------------------*/
static
double          _clipm_priv_image_get_pixel_area_in_polygon(
                                            int             N,
                                            const double    *xp,
                                            const double    *yp,
                                            double          x,
                                            double          y)
{
    int     xi,
            yi;
    double  xs,
            ys,
            w = 0.0;
    
    ys = y - 0.45;
    for (yi = 0; yi < 10; yi++, ys += 0.1)
    {
        xs = x - 0.45;
        for (xi = 0; xi < 10; xi++, xs += 0.1)
            w += _clipm_priv_image_is_point_in_polygon(N, xp, yp, xs, ys);
    }
    w /= 100;
    
    return w;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief   Return a mask according to a window flagging polygon edge pixels.
 * @param   window_xxyy     Coordinate buffer of the form
 *                          {xa, xb, ya, yb}, can be NULL,
 *                          minimum/maximum order is irrelevant
 * @param   N               Number of polygon points
 * @param   xp              X polygon corner coordinates
 * @param   yp              Y polygon corner coordinates
 * @return  The mask, NULL in the case of error
 * 
 * @par Error Handling:
 * The following errors can be set and returned:
 * - CPL_ERROR_NULL_INPUT: any input pointer is NULL
 * 
 * @todo
 * - implement unit test (through a wrapper func)
 * - the implementation is inefficient, since all pixels are checked. it could
 *   be improved by walking along the lines and marking the pixels (wherever
 *   a polygon segment crosses a horizontal or vertical pixel border).
 */
/*----------------------------------------------------------------------------*/
static
cpl_mask        *_clipm_priv_image_mark_polygon_edges(
                                            const int       window_xxyy[4],
                                            int             N,
                                            const double    *xp,
                                            const double    *yp)
{
    cpl_mask    *mlow = NULL,
                *mleft = NULL;
    
    CLIPM_TRY
    {
        cpl_binary  *mlowdat,
                    *mleftdat;
        int         wdwsize[2],
                    ll[2],
                    ur[2];
        int         x,
                    y,
                    dim;
        double      xpos,
                    ypos;
        
        CLIPM_TRY_CHECK_AUTOMSG(            window_xxyy != NULL,
                                            CPL_ERROR_NULL_INPUT);
        CLIPM_TRY_CHECK_AUTOMSG(            xp != NULL,
                                            CPL_ERROR_NULL_INPUT);
        CLIPM_TRY_CHECK_AUTOMSG(            yp != NULL,
                                            CPL_ERROR_NULL_INPUT);
        
        for (dim = 0; dim < 2; dim++)
        {
            ll[dim] = clipm_priv_math_min(  window_xxyy[2*dim],
                                            window_xxyy[2*dim+1]);
            ur[dim] = clipm_priv_math_max(  window_xxyy[2*dim],
                                            window_xxyy[2*dim+1]);
            wdwsize[dim] = ur[dim] - ll[dim] + 1;
        }
        
        mlow = cpl_mask_new(wdwsize[0], wdwsize[1]);
        mleft = cpl_mask_new(wdwsize[0], wdwsize[1]);
        CLIPM_TRY_CHECK_ERROR_STATE();
        
        /* fill two masks, representing whether the (lower|left) pixel
         * edge cuts a segment of the polygon */
        mlowdat = cpl_mask_get_data(mlow);
        mleftdat = cpl_mask_get_data(mleft);
        CLIPM_TRY_ASSERT_ERROR_STATE();
        ypos = (double)(ll[1]) - 0.5;
        for (y = 0; y < wdwsize[1]; y++, ypos += 1.0)
        {
            xpos = (double)(ll[0]) - 0.5;
            for (x = 0; x < wdwsize[0]; x++, xpos += 1.0)
            {
                *mlowdat++ = _clipm_priv_image_segment_cuts_polygon(
                                            N,
                                            xp, yp,
                                            xpos, ypos,
                                            xpos + 1.0, ypos);
                *mleftdat++ = _clipm_priv_image_segment_cuts_polygon(
                                            N,
                                            xp, yp,
                                            xpos, ypos,
                                            xpos, ypos + 1.0);
            }
        }
        
        /* propagate: if the (lower|left) pixel edge touches the polygon,
         * then the (upper|right) edge of the (lower|left) pixel touches
         * it as well */
        mlowdat = cpl_mask_get_data(mlow);
        mleftdat = cpl_mask_get_data(mleft);
        CLIPM_TRY_ASSERT_ERROR_STATE();
        /* propagate down */
        for (y = 0; y < wdwsize[1] - 1; y++)
        {
            for (x = 0; x < wdwsize[0]; x++)
            {
                *mlowdat |= mlowdat[wdwsize[0]];
                mlowdat++;
            }
        }
        ypos = (double)ur[1] + 0.5;
        xpos = (double)ll[0] - 0.5;
        for (x = 0; x < wdwsize[0]; x++, xpos += 1.0)
        {
            *mlowdat++ |= _clipm_priv_image_segment_cuts_polygon(
                                            N,
                                            xp, yp,
                                            xpos, ypos,
                                            xpos + 1.0, ypos);;
        }
        /* propagate left */
        ypos = (double)ll[1] - 0.5;
        xpos = (double)ur[0] + 0.5;
        for (y = 0; y < wdwsize[1]; y++, ypos += 1.0)
        {
            for (x = 0; x < wdwsize[0] - 1; x++)
            {
                *mleftdat |= mleftdat[1];
                mleftdat++;
            }
            *mleftdat++ |= _clipm_priv_image_segment_cuts_polygon(
                                            N,
                                            xp, yp,
                                            xpos, ypos,
                                            xpos, ypos + 1.0);;
        }
        
        /* 
         * create mask indicating whether at least one of a pixel's edges
         * touches the polygon edge
         */
        cpl_mask_or(mlow, mleft);
        CLIPM_TRY_ASSERT_ERROR_STATE();
    }
    CLIPM_CATCH
    {
        cpl_mask_delete(mlow);
        mlow = NULL;
    }

    cpl_mask_delete(mleft);
    
    return mlow;
}

/*-----------------------------------------------------------------------------
    Implementation
 -----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
 * @brief   Create an image by wrapping it around an existing data buffer.
 * @param   size                Buffer of size 2, containing the image size
 *                              [x, y], can be NULL
 * @param   type                Either CPL_TYPE_INT, CPL_TYPE_FLOAT, or
 *                              CPL_TYPE_DOUBLE
 * @param   data                The data buffer of the specified type and size
 * @return  Newly allocated image, NULL in the case of error
 * 
 * @par Error Handling:
 * Possible error codes set in this function:
 * - CPL_ERROR_NULL_INPUT:
 *   - @a size is NULL
 *   - @a data is NULL
 * - CPL_ERROR_ILLEGAL_INPUT: an entry of size is <= 0
 * - CPL_ERROR_INVALID_TYPE: if the passed image type is not supported
 * 
 * @note
 * - Specifying a @a size buffer with less then 2 entries
 *   (sizeof(size) < 2)
 *   will crash the application!
 */
/*----------------------------------------------------------------------------*/
cpl_image       *clipm_priv_image_wrap(     const int   *size,
                                            cpl_type    type,
                                            void        *data)
{
    cpl_image   *out = NULL;
    
    CLIPM_TRY
    {
        CLIPM_TRY_CHECK_AUTOMSG(            size != NULL,
                                            CPL_ERROR_NULL_INPUT);
        CLIPM_TRY_CHECK_AUTOMSG(            data != NULL,
                                            CPL_ERROR_NULL_INPUT);
        switch(type)
        {
            case    CPL_TYPE_INT:
                out = cpl_image_wrap_int(   size[0],
                                            size[1],
                                            data);
                break;
            case    CPL_TYPE_FLOAT:
                out = cpl_image_wrap_float( size[0],
                                            size[1],
                                            data);
                break;
            case    CPL_TYPE_DOUBLE:
                out = cpl_image_wrap_double(size[0],
                                            size[1],
                                            data);
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
        clipm_priv_image_null(&out);
    }
    
    return out;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief   Return closest pixel which is not flagged as bad.
 * @param   image   Input image
 * @param   x       X position (FITS convention)
 * @param   y       Y position (FITS convention)
 * @param   out_x   (Optional output) x index of found pixel (FITS convention)
 * @param   out_y   (Optional output) y index of found pixel (FITS convention)
 * @return  The pixel value, -1.0 on error
 * 
 * @par Principle:
 * Returns the pixel value at the respective rounded position. If the pixel is
 * bad, the next good pixel in minimum distance is searched (in circles around
 * the provided position).
 * 
 * @par Constraints:
 * Images can be of type CPL_TYPE_INT, CPL_TYPE_FLOAT or CPL_TYPE_DOUBLE.
 * 
 * @par Error Handling:
 * Possible error codes set in this function:
 * - CPL_ERROR_NULL_INPUT: if @a image is NULL
 * - CPL_ERROR_INVALID_TYPE: if the passed image type is not supported
 * - CPL_ERROR_ACCESS_OUT_OF_RANGE: if x or y are out of the image range
 * - CPL_ERROR_DATA_NOT_FOUND: if no good pixel was found
 * .
 * In the case of error, the return values are:
 * - Function: -1.0
 * - @a out_x: 0
 * - @a out_y: 0
 */
/*----------------------------------------------------------------------------*/
double          clipm_priv_image_get_nearest_good(
                                            const cpl_image *image,
                                            double          x,
                                            double          y,
                                            int             *out_x,
                                            int             *out_y)
{
    double val = -1.0;
    
    CLIPM_TRY
    {
        cpl_type    type;
        int         imsize[2];
        const void  *data;
        const cpl_binary
                    *badp;
        double      dmax;
        
        clipm_priv_image_get_data_const(    image,
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
                x >= 1 && x <= imsize[0] &&
                y >= 1 && y <= imsize[1],
            CPL_ERROR_ACCESS_OUT_OF_RANGE);
        
        dmax = sqrt(imsize[0]*imsize[0] + imsize[1]*imsize[1]);
        
        /* FITS */
        x -= 1.0;
        y -= 1.0;
        
/* --- throw in a macro to ease life ---------------------------------------- */
#define         clipm_priv_image_get_nearest_good_BODY(TYPE) \
do { \
    const TYPE  *Tdata = data; \
    int         xi, \
                yi, \
                i; \
    xi = (int)(x + 0.5); \
    yi = (int)(y + 0.5); \
    \
    i = yi*imsize[0] + xi; \
    if (badp == NULL || ! badp[i]) \
    { \
        val = Tdata[i]; \
        if (out_x != NULL) \
            *out_x = xi + 1; \
        if (out_y != NULL) \
            *out_y = yi + 1; \
    } \
    else \
    { \
        double  r; \
        int     found = 0; \
        \
        for (   r = CPL_MATH_SQRT1_2; \
                r < dmax && ! found; \
                r += CPL_MATH_SQRT1_2) \
        { \
            double  a, \
                    dmin = r + 1, \
                    dx, \
                    dy, \
                    d; \
            for (a = 0; a < CPL_MATH_2PI; a += 1.0/r) \
            { \
                xi = (int)(x + 0.5 + r * cos(a)); \
                yi = (int)(y + 0.5 + r * sin(a)); \
                if (xi >= 0 && xi < imsize[0] && \
                    yi >= 0 && yi < imsize[1]) \
                { \
                    i = yi*imsize[0] + xi; \
                    if (! badp[i]) \
                    { \
                        dx = x - (double)xi; \
                        dy = y - (double)yi; \
                        d = sqrt(dx*dx + dy*dy); \
                        if (d < dmin) \
                        { \
                            found = 1; \
                            val = Tdata[i]; \
                            dmin = d; \
                            if (out_x != NULL) \
                                *out_x = xi + 1; \
                            if (out_y != NULL) \
                                *out_y = yi + 1; \
                        } \
                    } \
                } \
            } \
        } \
        \
        if (! found) \
            CLIPM_TRY_EXIT_WITH_ERROR_MSG(  CPL_ERROR_DATA_NOT_FOUND, \
                                            "", "no good pixel found"); \
    } \
} while (0)
/* -------------------------------------------------------------------------- */
        
        switch (type)
        {
            case CPL_TYPE_INT: {
                clipm_priv_image_get_nearest_good_BODY(int);
                break;
            }
            case CPL_TYPE_FLOAT: {
                clipm_priv_image_get_nearest_good_BODY(float);
                break;
            }
            case CPL_TYPE_DOUBLE: {
                clipm_priv_image_get_nearest_good_BODY(double);
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
        val = - 1.0;
        if (out_x != NULL)
            *out_x = 0;
        if (out_y != NULL)
            *out_y = 0;
    }
    
    return val;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief   Return bad pixel mask data, if there are bad pixels, otherwise NULL.
 * @param   image   Input image
 * @return  BPM pointer, NULL in case there are no bad pixels, or in the case of
 *          error
 * 
 * @par Error Handling:
 * The following error codes can be set:
 * - CPL_ERROR_NULL_INPUT: @a image is NULL
 */
/*----------------------------------------------------------------------------*/
const cpl_binary    *clipm_priv_image_bpm_get_if_exist(
                                            const cpl_image *image)
{
    const cpl_binary    *badp_data = NULL;
    
    CLIPM_TRY
    {
        CLIPM_TRY_CHECK_AUTOMSG(            image != NULL,
                                            CPL_ERROR_NULL_INPUT);
        
        if (cpl_image_count_rejected(image) > 0)
        {
            const cpl_mask  *bpm;
            bpm = cpl_image_get_bpm_const(image);
            CLIPM_TRY_CHECK_ERROR_STATE();
            badp_data = cpl_mask_get_data_const(bpm);
        }
    }
    CLIPM_CATCH
    {
        badp_data = NULL;
    }
    
    return badp_data;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief   Perform a logical OR between the bad pixel masks of 2 images.
 * @param   mod_image   Image to modify
 * @param   ref_image   Reference image
 * @return  CPL error code
 * 
 * @par Error Handling:
 * The following error codes can be set:
 * - CPL_ERROR_NULL_INPUT: @a mod_image or @a ref_image is NULL
 * - CPL_ERROR_INCOMPATIBLE_INPUT: @a mod_image and @a ref_image have different
 *   sizes
 */
/*----------------------------------------------------------------------------*/
cpl_error_code  clipm_priv_image_bpm_or(    cpl_image       *mod_image,
                                            const cpl_image *ref_image)
{
    CLIPM_TRY
    {
        int                 xsize,
                            ysize;
        
        clipm_priv_checks_images_match(     mod_image,
                                            ref_image,
                                            NULL,
                                            &xsize,
                                            &ysize);
        
        if (cpl_image_count_rejected(ref_image) > 0)
        {
            const cpl_mask  *ref_bpm;
            cpl_mask        *mod_bpm;
            
            ref_bpm = cpl_image_get_bpm_const(ref_image);
            mod_bpm = cpl_image_get_bpm(mod_image);
            CLIPM_TRY_CHECK_ERROR_STATE();
            
            cpl_mask_or(mod_bpm, ref_bpm);
        }        
    }
    CLIPM_CATCH
    {
    }
    
    return CLIPM_ERROR_GET_NEW_SINCE_TRY();
}

/*----------------------------------------------------------------------------*/
/**
 * @brief   Mark the border of an image as bad.
 * @param   image   Image to modify
 * @param   left    Left border width
 * @param   right   Right border width
 * @param   lower   Lower border width
 * @param   upper   Upper border width
 * @return  CPL error code
 * 
 * @par Principle:
 * - The width parameters must be >= 0, no maximum limit is required (this
 *   means, some width value big enough will just flag the whole image as
 *   bad).
 * 
 * @par Error Handling:
 * The following error codes can be set:
 * - CPL_ERROR_NULL_INPUT: @a image is NULL
 * - CPL_ERROR_ILLEGAL_INPUT: a width parameter is < 0
 */
/*----------------------------------------------------------------------------*/
cpl_error_code  clipm_priv_image_bpm_border_bad(
                                            cpl_image   *image,
                                            int         left,
                                            int         right,
                                            int         lower,
                                            int         upper)
{
    CLIPM_TRY
    {
        int         xsize,
                    ysize;
        cpl_mask    *bpm;
        cpl_binary  *bpmdata;
        
        CLIPM_TRY_CHECK_AUTOMSG(            image != NULL,
                                            CPL_ERROR_NULL_INPUT);
        CLIPM_TRY_CHECK_AUTOMSG(            left >= 0,
                                            CPL_ERROR_ILLEGAL_INPUT);
        CLIPM_TRY_CHECK_AUTOMSG(            right >= 0,
                                            CPL_ERROR_ILLEGAL_INPUT);
        CLIPM_TRY_CHECK_AUTOMSG(            lower >= 0,
                                            CPL_ERROR_ILLEGAL_INPUT);
        CLIPM_TRY_CHECK_AUTOMSG(            upper >= 0,
                                            CPL_ERROR_ILLEGAL_INPUT);
        
        xsize = cpl_image_get_size_x(image);
        ysize = cpl_image_get_size_y(image);
        CLIPM_TRY_ASSERT_ERROR_STATE();
        
        if (left == 0 && right == 0
            && lower == 0 && upper == 0)
        CLIPM_TRY_EXIT();
        
        bpm = cpl_image_get_bpm(image);
        bpmdata = cpl_mask_get_data(bpm);
        CLIPM_TRY_ASSERT_ERROR_STATE();
        
        if (left + right >= xsize || lower + upper >= ysize)
        {
            int n,
                tot;
            tot = xsize * ysize;
            for (n = 0; n < tot; n++)
                *bpmdata++ = CPL_BINARY_1;
        }
        else
        {
            int         x,
                        y,
                        xbetween,
                        ybetween;
            
            for (y = 0; y < lower; y++)
                for (x = 0; x < xsize; x++)
                    *bpmdata++ = CPL_BINARY_1;
            xbetween = xsize - left - right;
            ybetween = ysize - lower - upper;
            for (y = 0; y < ybetween; y++)
            {
                for (x = 0; x < left; x++)
                    *bpmdata++ = CPL_BINARY_1;
                bpmdata += xbetween;
                for (x = 0; x < right; x++)
                    *bpmdata++ = CPL_BINARY_1;
            }
            for (y = 0; y < upper; y++)
                for (x = 0; x < xsize; x++)
                    *bpmdata++ = CPL_BINARY_1;
        }
    }
    CLIPM_CATCH
    {
    }
    
    return CLIPM_ERROR_GET_NEW_SINCE_TRY();
}

/*----------------------------------------------------------------------------*/
/**
 * @brief   Reject all pixels equal to or greater than a limit.
 * @param   image           Image to modify
 * @param   window_xxyy     Coordinate buffer of the form
 *                          {xa, xb, ya, yb}, can be NULL,
 *                          minimum/maximum order is irrelevant
 * @param   limit           Limit
 * @return  CPL error code
 * 
 * @par Bad Pixel Handling:
 * Bad pixels are kept as bad.
 * 
 * @par Error Handling:
 * The following errors can be set and returned:
 * - CPL_ERROR_NULL_INPUT: @a image is NULL
 * - CPL_ERROR_INVALID_TYPE: @a image is not of type CPL_TYPE_INT,
 *   CPL_TYPE_FLOAT or CPL_TYPE_DOUBLE
 * - CPL_ERROR_ACCESS_OUT_OF_RANGE: @a window_xxyy specifies coordinates outside
 *   the image range
 */
/*----------------------------------------------------------------------------*/
cpl_error_code  clipm_priv_image_bpm_reject_above(
                                            cpl_image       *image,
                                            const int       *window_xxyy,
                                            double          limit)
{
    CLIPM_TRY
    {
        int         rejected = 0,
                    had_bad_pix;
        int         imsize[2],
                    wdwsize[2],
                    start[2];
        cpl_type    type;
        const void  *data_start = NULL;
        const cpl_binary
                    *badp_start = NULL;
        
        clipm_priv_image_get_data_const(    image,
                                            window_xxyy,
                                            1,
                                            imsize,
                                            wdwsize,
                                            start,
                                            &type,
                                            &data_start,
                                            &badp_start);
        CLIPM_TRY_CHECK_ERROR_STATE();
        
        had_bad_pix = (badp_start != NULL);
        if (!had_bad_pix)
        {
            cpl_mask    *bpm;
            bpm = cpl_image_get_bpm(image);
            badp_start = cpl_mask_get_data(bpm);
            badp_start += start[1]*imsize[0] + start[0];
        }
        
/* --- throw in a macro to ease life ---------------------------------------- */
#define clipm_priv_image_bpm_reject_above_BODY(  TYPE) \
do { \
    TYPE        *data = (TYPE*)data_start; \
    cpl_binary  *badp = (cpl_binary*)badp_start; \
    TYPE        tlim = (TYPE)limit; \
    int     x, \
            y; \
    \
    if (had_bad_pix) \
    { \
        for (y = 0; y < wdwsize[1]; y++) \
        { \
            for (x = 0; x < wdwsize[0]; x++) \
            { \
                if (!badp[x] && data[x] >= tlim) \
                { \
                    badp[x] = CPL_BINARY_1; \
                    rejected++; \
                } \
            } \
            data += imsize[0]; \
            badp += imsize[0]; \
        } \
    } \
    else \
    { \
        for (y = 0; y < wdwsize[1]; y++) \
        { \
            for (x = 0; x < wdwsize[0]; x++) \
            { \
                if (data[x] >= tlim) \
                { \
                    badp[x] = CPL_BINARY_1; \
                    rejected++; \
                } \
            } \
            data += imsize[0]; \
            badp += imsize[0]; \
        } \
    } \
} while (0)
/* -------------------------------------------------------------------------- */
        switch(type)
        {
            case    CPL_TYPE_INT:
                clipm_priv_image_bpm_reject_above_BODY(int);
                break;
            case    CPL_TYPE_FLOAT:
                clipm_priv_image_bpm_reject_above_BODY(float);
                break;
            case    CPL_TYPE_DOUBLE:
                clipm_priv_image_bpm_reject_above_BODY(double);
                break;
            default:
                CLIPM_TRY_EXIT_WITH_ERROR_MSG(
                                            CPL_ERROR_INVALID_TYPE,
                                            "image",
                                            "must be int, float or double");
        }
        
        /* erase bpm */
        if (rejected == 0 && !had_bad_pix)
            cpl_image_accept_all(image);
    }
    CLIPM_CATCH
    {
    }
    
    return CLIPM_ERROR_GET_NEW_SINCE_TRY();
}

/*----------------------------------------------------------------------------*/
/**
 * @brief   Copy all non-bad pixels into a new vector of appropriate size.
 * @param   image       Input image
 * @param   window_xxyy Coordinate buffer of the form
 *                      {x1a, x1b, y1a, y1b},
 *                      minimum/maximum order is irrelevant
 * @param   out_ngood   (Optional output) number of good pixels,
 *                      contains 0 in the case of error
 * @return  The new vector, NULL in the case of error
 * 
 * @par Error Handling:
 * The following errors may occur:
 * - CPL_ERROR_NULL_INPUT: @a image is NULL
 * - CPL_ERROR_INVALID_TYPE: @a image is not of type CPL_TYPE_INT,
 *   CPL_TYPE_FLOAT or CPL_TYPE_DOUBLE
 * - CPL_ERROR_ACCESS_OUT_OF_RANGE: @a window_xxyy specifies coordinates outside
 *   the image range
 * - CPL_ERROR_DATA_NOT_FOUND: @a image contains only bad pixels
 */
/*----------------------------------------------------------------------------*/
cpl_vector      *clipm_priv_image_copy_good_data_vector(
                                            const cpl_image *image,
                                            const int       *window_xxyy,
                                            int             *out_ngood)
{
    cpl_vector  *out = NULL;
    
    CLIPM_TRY
    {
        int         imsize[2],              /* x, y */
                    wdwsize[2],
                    start[2];
        int         ngood;
        cpl_type    type;
        const void  *data_start = NULL;
        const cpl_binary
                    *badp_start = NULL;
        double      *vdat;
        
        clipm_priv_image_get_data_const(    image,
                                            window_xxyy,
                                            1,
                                            imsize,
                                            wdwsize,
                                            start,
                                            &type,
                                            &data_start,
                                            &badp_start);
        CLIPM_TRY_CHECK_ERROR_STATE();
        
        ngood = wdwsize[0]*wdwsize[1];
        if (badp_start != NULL)
        {
            int x,
                y;
            const cpl_binary    *badp = badp_start;
            for (y = 0; y < wdwsize[1]; y++, badp += imsize[0])
                for (x = 0; x < wdwsize[0]; x++)
                    if (badp[x])
                        ngood--;
        }
        
        if (out_ngood != NULL)
            *out_ngood = ngood;
        
        CLIPM_TRY_CHECK(                    ngood > 0,
                                            CPL_ERROR_DATA_NOT_FOUND,
                                            "image",
                                            "has only bad pixels");
        
        out = cpl_vector_new(ngood);
        vdat = cpl_vector_get_data(out);
        CLIPM_TRY_ASSERT_ERROR_STATE();
        
/* --- throw in a macro to ease life ---------------------------------------- */
#define clipm_priv_image_copy_good_data_vector_BODY(  TYPE) \
do { \
    const TYPE          *data = (TYPE*)data_start; \
    const cpl_binary    *badp = badp_start; \
    int     x, \
            y; \
    \
    if (badp_start != NULL) \
    { \
        for (y = 0; y < wdwsize[1]; y++) \
        { \
            for (x = 0; x < wdwsize[0]; x++) \
            { \
                if (!badp[x]) \
                { \
                    *vdat++ = data[x]; \
                } \
            } \
            data += imsize[0]; \
            badp += imsize[0]; \
        } \
    } \
    else \
    { \
        for (y = 0; y < wdwsize[1]; y++) \
        { \
            for (x = 0; x < wdwsize[0]; x++) \
                *vdat++ = data[x]; \
            \
            data += imsize[0]; \
        } \
    } \
} while (0)
/* -------------------------------------------------------------------------- */
        switch(type)
        {
            case    CPL_TYPE_INT:
                clipm_priv_image_copy_good_data_vector_BODY(int);
                break;
            case    CPL_TYPE_FLOAT:
                clipm_priv_image_copy_good_data_vector_BODY(float);
                break;
            case    CPL_TYPE_DOUBLE:
                clipm_priv_image_copy_good_data_vector_BODY(double);
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
        clipm_priv_vector_null(&out);
        
        if (out_ngood != NULL)
            *out_ngood = 0;
    }
    
    return out;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief   Get image and window sizes, and data pointers to image (window).
 * @param   image               Image (FITS convention)
 * @param   window_xxyy         Coordinate buffer of the form
 *                              {x1a, x1b, y1a, y1b},
 *                              minimum/maximum order is irrelevant
 * @param   allow_window_NULL   If !0, then no error is set if (window==NULL)
 * @param   img_size_xy         (Optional) buffer of size 2, to which the image
 *                              size [x, y] is put out, can be NULL
 * @param   window_size_xy      (Optional) buffer of size 2, to which the window
 *                              size is determined and written, can be NULL
 * @param   buffer_start_xy     (Optional) buffer of size 2, to which the lower
 *                              left window indices are written (buffer,
 *                              not FITS), starting at 0,
 *                              can be NULL
 * @param   type                (Optional) output CPL-type of image
 * @param   data_start          (Optional) output pointer to image data, if
 *                              @a window_xxyy is given, then it points to the
 *                              lower left corner of the data
 * @param   badp_start          (Optional) output bad pixel map pointer,
 *                              returns NULL if no bad pixels are set, if
 *                              @a window_xxyy is given, then it points to the
 *                              lower left corner of the data
 * @return  CPL error code
 * 
 * @par Overview:
 * - The window coordinates are checked for being inside the image range.
 * - If an input is NULL, or if the coordinates are out of range, then the
 *   appropriate CPL error is set and returned.
 * - If @a out_window_size != NULL, then the window size is calculated and
 *   written to this buffer. If @a window_xxyy == NULL, then @a out_window_size
 *   is set to the image size!
 * - @a data_start and @a badp_start are (if provided) set to the respective
 *   buffers, if @a window_xxyy is provided, they point to the position in the
 *   data buffers corresponding to the lower left window corner.
 * 
 * @par Error Handling:
 * The following codes are set, if the check fails:
 * - CPL_ERROR_NULL_INPUT: @a image and/or @a window_xxyy is NULL
 * - CPL_ERROR_ACCESS_OUT_OF_RANGE: @a window_xxyy coordinates are outside the
 *   image
 * In the case of error, the output parameter values are undefined.
 * 
 * @note
 * - This function does not check min/max conditions of @a window_xxyy. Use
 *   clipm_priv_checks_window_minmax() for this purpose.
 * - Specifying buffers of too small size
 *   (sizeof(window) < 4 || sizeof(out_img_size) < 2 ||
 *   sizeof(out_window_size) < 2 || sizeof(out_buffer_start_xy) < 2)
 *   will crash the application!
 */
/*----------------------------------------------------------------------------*/
cpl_error_code  clipm_priv_image_get_data_const(
                                            const cpl_image *image,
                                            const int       *window_xxyy,
                                            int             allow_window_NULL,
                                            int             *img_size_xy,
                                            int             *window_size_xy,
                                            int             *buffer_start_xy,
                                            cpl_type        *type,
                                            const void      **data_start,
                                            const cpl_binary
                                                            **badp_start)
{
    CLIPM_TRY
    {
        int         tmp_imsize[2],
                    tmp_bufferstart[2];
        int         pointer_offset = 0;
        cpl_type    _type;
        
        /* init output */
        if (data_start != NULL)
            *data_start = NULL;
        if (badp_start != NULL)
            *badp_start = NULL;
        if (type != NULL)
            *type = CPL_TYPE_INVALID;
        
        /* redirect output pointers which are needed */
        if (img_size_xy == NULL)
            img_size_xy = (int*)tmp_imsize;
        if (buffer_start_xy == NULL)
            buffer_start_xy = (int*)tmp_bufferstart;
        if (type == NULL)
            type = &_type;
        
        clipm_priv_checks_window_image(     window_xxyy,
                                            image,
                                            allow_window_NULL,
                                            img_size_xy,
                                            window_size_xy,
                                            buffer_start_xy);
        CLIPM_TRY_CHECK_ERROR_STATE();
        
        *type = cpl_image_get_type(image);
        CLIPM_TRY_CHECK_ERROR_STATE();
        
        pointer_offset =    buffer_start_xy[1] * img_size_xy[0] +
                            buffer_start_xy[0];
        
        if (data_start != NULL)
        {
            *data_start = (const char*)cpl_image_get_data_const(image)
                            + cpl_type_get_sizeof(*type) * pointer_offset;
            CLIPM_TRY_CHECK_ERROR_STATE();
        }
        
        if (badp_start != NULL)
        {
            *badp_start = clipm_priv_image_bpm_get_if_exist(image);
            if (*badp_start != NULL)
                *badp_start += pointer_offset;
        }
    }
    CLIPM_CATCH
    {
        if (data_start != NULL)
            *data_start = NULL;
        if (badp_start != NULL)
            *badp_start = NULL;
        *type = CPL_TYPE_INVALID;
    }
    
    return CLIPM_ERROR_GET_NEW_SINCE_TRY();
}

/*----------------------------------------------------------------------------*/
/**
 * @brief   Extract an image window and cast into a new image.
 * @param   input       Input image (FITS convention)
 * @param   window_xxyy Coordinate buffer of the form
 *                      {xa, xb, ya, yb}, can be NULL,
 *                      minimum/maximum order is irrelevant
 * @param   outtype     Output image type
 * @return  The new extracted and casted image, NULL in the case of error
 * 
 * @par Principle:
 * The image data (from a window) is copied and casted to the new type.
 * Valid types are CPL_TYPE_INT, CPL_TYPE_FLOAT and CPL_TYPE_DOUBLE (for
 * both input and output images).
 * 
 * @par Error Handling:
 * Possible error codes set in this function:
 * - CPL_ERROR_NULL_INPUT if @a image is NULL
 * - CPL_ERROR_INVALID_TYPE:
 *   - @a out_type is not supported
 *   - the type of the @a input image is not supported
 * - CPL_ERROR_ACCESS_OUT_OF_RANGE if a coordinate is out of the image range
 * 
 * @todo
 * - implement unit test
 */
/*----------------------------------------------------------------------------*/
cpl_image       *clipm_priv_image_extract_cast(
                                            const cpl_image *input,
                                            const int       *window_xxyy,
                                            cpl_type        outtype)
{
    cpl_image   *output = NULL;

    CLIPM_TRY
    {
        int         imsize[2],
                    windowsize[2],
                    buffer_start[2];
        cpl_type    intype;
        const void  *data_start = NULL;
        const cpl_binary
                    *badp_start = NULL;
        
        clipm_priv_image_get_data_const(    input,
                                            window_xxyy,
                                            1,
                                            imsize,
                                            windowsize,
                                            buffer_start,
                                            &intype,
                                            &data_start,
                                            &badp_start);
        CLIPM_TRY_CHECK_ERROR_STATE();
        
/* --- throw in a macro to ease life ---------------------------------------- */
/* copy the image window into a double buffer, and subtract the mean */
#define clipm_COPY_WINDOW_TO_BUFFER(        OUT_CTYPE, \
                                            INTYPE) \
do { \
    const INTYPE    *imdata; \
    OUT_CTYPE       *out_data; \
    int             x, y; \
    \
    imdata = data_start; \
    out_data = cpl_malloc(                  windowsize[0] \
                                            * windowsize[1] \
                                            * sizeof(*out_data)); \
    output = clipm_priv_image_wrap(         windowsize, \
                                            outtype, \
                                            out_data); \
    if (output == NULL) \
    { \
        cpl_free(out_data); \
        out_data = NULL; \
    } \
    CLIPM_TRY_ASSERT_ERROR_STATE(); \
    \
    if (badp_start == NULL) \
    { \
        for (y = 0; y < windowsize[1]; y++, imdata += imsize[0]) \
            for (x = 0; x < windowsize[0]; x++) \
                *(out_data++) = imdata[x]; \
    } \
    else \
    { \
        cpl_mask    *out_bpm = NULL; \
        cpl_binary  *out_badp = NULL; \
        const cpl_binary \
                    *in_badp; \
        out_bpm = cpl_image_get_bpm(output); \
        CLIPM_TRY_CHECK_ERROR_STATE(); \
        out_badp = cpl_mask_get_data(out_bpm); \
        CLIPM_TRY_CHECK_ERROR_STATE(); \
        in_badp = badp_start; \
        \
        for (   y = 0; \
                y < windowsize[1]; \
                y++, imdata += imsize[0], in_badp += imsize[0]) \
            for (x = 0; x < windowsize[0]; x++) \
            { \
                *(out_data++) = imdata[x]; \
                *(out_badp++) = in_badp[x]; \
            } \
    } \
} while (0)
/* -------------------------------------------------------------------------- */

/* --- throw in another macro to ease life ---------------------------------- */
/* call the above macro multiple times */
#define clipm_priv_image_extract_SWITCH_INTYPE( \
                                            OUT_CTYPE) \
do { \
                switch(intype) { \
                    case CPL_TYPE_INT: \
                        clipm_COPY_WINDOW_TO_BUFFER( \
                                            OUT_CTYPE, \
                                            int); \
                        break; \
                    case CPL_TYPE_FLOAT: \
                        clipm_COPY_WINDOW_TO_BUFFER( \
                                            OUT_CTYPE, \
                                            float); \
                        break; \
                    case CPL_TYPE_DOUBLE: \
                        clipm_COPY_WINDOW_TO_BUFFER( \
                                            OUT_CTYPE, \
                                            double); \
                        break; \
                    default: /* don't exit */ \
                        CLIPM_ERROR_SET_MSG(CPL_ERROR_INVALID_TYPE, \
                                            "input", \
                                            "must be int, float or double"); \
                } \
} while (0)
/* -------------------------------------------------------------------------- */

        switch(outtype) {
            case CPL_TYPE_INT:
                clipm_priv_image_extract_SWITCH_INTYPE(
                                            int);
                break;
            case CPL_TYPE_FLOAT:
                clipm_priv_image_extract_SWITCH_INTYPE(
                                            float);
                break;
            case CPL_TYPE_DOUBLE:
                clipm_priv_image_extract_SWITCH_INTYPE(
                                            double);
                break;
            default:
                CLIPM_ERROR_SET_MSG(        CPL_ERROR_INVALID_TYPE,
                                            "outtype",
                                            "must be int, float or double");
        }
    }
    CLIPM_CATCH
    {
        clipm_priv_image_null(&output);
    }
    
    return output;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief   Extract an image window and round to integer.
 * @param   input       Input image (FITS convention)
 * @param   window_xxyy Coordinate buffer of the form
 *                      {xa, xb, ya, yb}, can be NULL,
 *                      minimum/maximum order is irrelevant
 * @return  The new rounded image, NULL in the case of error
 * 
 * @par Principle:
 * The image data (from a window) is copied and rounded to integer.
 * 
 * @par Error Handling:
 * Possible error codes set in this function:
 * - CPL_ERROR_NULL_INPUT if @a image is NULL
 * - CPL_ERROR_INVALID_TYPE if the input image type is not CPL_TYPE_INT,
 *   CPL_TYPE_FLOAT or CPL_TYPE_DOUBLE
 * - CPL_ERROR_ACCESS_OUT_OF_RANGE if a coordinate is out of the image range
 */
/*----------------------------------------------------------------------------*/
cpl_image       *clipm_priv_image_extract_round(
                                            const cpl_image *input,
                                            const int       *window_xxyy)
{
    cpl_image   *output = NULL;

    CLIPM_TRY
    {
        int         imsize[2],
                    windowsize[2],
                    buffer_start[2];
        cpl_type    intype;
        const void  *data_start = NULL;
        const cpl_binary
                    *badp_start = NULL;
        
        clipm_priv_image_get_data_const(    input,
                                            window_xxyy,
                                            1,
                                            imsize,
                                            windowsize,
                                            buffer_start,
                                            &intype,
                                            &data_start,
                                            &badp_start);
        CLIPM_TRY_CHECK_ERROR_STATE();
        
/* --- throw in a macro to ease life ---------------------------------------- */
/* copy the image window into a double buffer, and subtract the mean */
#define clipm_priv_image_extract_round_BODY(INTYPE) \
do { \
    const INTYPE    *imdata; \
    int             *out_data; \
    int             x, y; \
    \
    imdata = data_start; \
    out_data = cpl_malloc(                  windowsize[0] \
                                            * windowsize[1] \
                                            * sizeof(*out_data)); \
    output = clipm_priv_image_wrap(         windowsize, \
                                            CPL_TYPE_INT, \
                                            out_data); \
    if (output == NULL) \
    { \
        cpl_free(out_data); \
        out_data = NULL; \
    } \
    CLIPM_TRY_ASSERT_ERROR_STATE(); \
    \
    if (badp_start == NULL) \
    { \
        if (intype == CPL_TYPE_INT) \
        { \
            for (y = 0; y < windowsize[1]; y++, imdata += imsize[0]) \
                for (x = 0; x < windowsize[0]; x++) \
                    *(out_data++) = imdata[x]; \
        } \
        else \
        { \
            for (y = 0; y < windowsize[1]; y++, imdata += imsize[0]) \
                for (x = 0; x < windowsize[0]; x++) \
                    *(out_data++) = clipm_math_round_d2i(imdata[x]); \
        } \
    } \
    else \
    { \
        cpl_mask    *out_bpm = NULL; \
        cpl_binary  *out_badp = NULL; \
        const cpl_binary \
                    *in_badp; \
        out_bpm = cpl_image_get_bpm(output); \
        CLIPM_TRY_CHECK_ERROR_STATE(); \
        out_badp = cpl_mask_get_data(out_bpm); \
        CLIPM_TRY_CHECK_ERROR_STATE(); \
        in_badp = badp_start; \
        \
        if (intype == CPL_TYPE_INT) \
        { \
            for (   y = 0; \
                    y < windowsize[1]; \
                    y++, imdata += imsize[0], in_badp += imsize[0]) \
                for (x = 0; x < windowsize[0]; x++) \
                { \
                    *(out_data++) = imdata[x]; \
                    *(out_badp++) = in_badp[x]; \
                } \
        } \
        else \
        { \
            for (   y = 0; \
                    y < windowsize[1]; \
                    y++, imdata += imsize[0], in_badp += imsize[0]) \
                for (x = 0; x < windowsize[0]; x++) \
                { \
                    *(out_data++) = clipm_math_round_d2i(imdata[x]); \
                    *(out_badp++) = in_badp[x]; \
                } \
        } \
    } \
} while (0)
/* -------------------------------------------------------------------------- */

        switch(intype) {
            case CPL_TYPE_INT:
                clipm_priv_image_extract_round_BODY(int);
                break;
            case CPL_TYPE_FLOAT:
                clipm_priv_image_extract_round_BODY(float);
                break;
            case CPL_TYPE_DOUBLE:
                clipm_priv_image_extract_round_BODY(double);
                break;
            default:
                CLIPM_ERROR_SET_MSG(CPL_ERROR_INVALID_TYPE,
                                    "input",
                                    "must be int, float or double");
        }
    }
    CLIPM_CATCH
    {
        clipm_priv_image_null(&output);
    }
    
    return output;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief   Apply a value to all pixels inside the image or a window of the
 *          image.
 * @param   image       The image to draw into (FITS convention)
 * @param   window_xxyy Coordinate buffer of the form
 *                      {xa, xb, ya, yb}, can be NULL,
 *                      minimum/maximum order is irrelevant
 * @param   value       The value to be set (true or false)
 * @return  CPL error code
 * 
 * @par Overview:
 * - If @a window is given, all pixels inside the window are set to the
 *   provided @a value.
 * - If @a window is NULL, all pixels in the image are set
 *   to the provided @a value.
 * - If the type of the @a image is CPL_TYPE_INT, then @a value is rounded to
 *   the nearest integer.
 * 
 * @par Bad Pixel Handling:
 * - Bad pixels are flagged as good.
 * 
 * @par Error Handling:
 * The following errors can be set and returned:
 * - CPL_ERROR_NULL_INPUT: @a image is NULL
 * - CPL_ERROR_INVALID_TYPE: @a image is not of type CPL_TYPE_INT,
 *   CPL_TYPE_FLOAT or CPL_TYPE_DOUBLE
 * - CPL_ERROR_ACCESS_OUT_OF_RANGE: @a window specifies coordinates outside
 *   the @a image plane
 * .
 * In the case of error, the @a image is not changed.
 * 
 * @todo
 * - implement unit test
 */
/*----------------------------------------------------------------------------*/
cpl_error_code  clipm_priv_image_fill(      cpl_image       *image,
                                            const int       *window_xxyy,
                                            double          value)
{
    CLIPM_TRY
    {
        int         imsize[2],
                    wdwsize[2];
        cpl_type    type;
        const void  *data_start = NULL;
        const cpl_binary
                    *badp_start = NULL;
        
        clipm_priv_image_get_data_const(    image,
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
#define clipm_priv_image_fill_BODY(  TYPE) \
do { \
    int         x, \
                y; \
    TYPE        *data; \
    TYPE        setvalue; \
    cpl_binary  *badp; \
    \
    if (type == CPL_TYPE_INT) \
        setvalue = clipm_math_round_d2i(value); \
    else \
        setvalue = value; \
    \
    for (   y = 0, data = (TYPE*)data_start, badp = (cpl_binary*)badp_start; \
            y < wdwsize[1]; \
            y++, data += imsize[0], badp += imsize[0]) \
    { \
        if (badp_start != NULL) \
            for (x = 0; x < wdwsize[0]; x++) \
            { \
                data[x] = setvalue; \
                badp[x] = CPL_BINARY_0; \
            } \
        else \
            for (x = 0; x < wdwsize[0]; x++) \
            { \
                data[x] = setvalue; \
            } \
    } \
} while (0)
/* -------------------------------------------------------------------------- */

        switch(type)
        {
            case    CPL_TYPE_INT:
                clipm_priv_image_fill_BODY( int);
                break;
            case    CPL_TYPE_FLOAT:
                clipm_priv_image_fill_BODY( float);
                break;
            case    CPL_TYPE_DOUBLE:
                clipm_priv_image_fill_BODY( double);
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
 * @brief   Fill a circle.
 * @param   img_modified    The image to draw into
 * @param   centre_xy       Centre position (double buffer of length 2 [x, y])
 * @param   radius          The circle's radius
 * @param   brightness      The circle's content value
 * @param   anti_alias      Flag whether to soften the edge
 * @param   additive        Flag whether to add the circle on top of the
 *                          existing image content (@a additive != 0),
 *                          or whether to set the area to the @a brightness
 *                          value
 * @return  CPL error code
 * 
 * @par Principle:
 * - If @a additive is true, then the filled circle will be added to the
 *   current image content, otherwise the circle is filled with @a brightness.
 * - If @a anti_alias is true, then the edges are softened.
 * - If the type of @a img_modified is CPL_TYPE_INT, then the resulting values
 *   are rounded.
 * - @a radius can be negative, in this case no pixel is modified, and no error
 *   is set. If the circle is outside the range, no error is set either.
 * 
 * @note
 * Be aware, that a very small circle can also fit between pixels, so that none
 * might be set (if @a anti_alias is false).
 * 
 * @par Bad Pixel Handling:
 * - If @a additive is false, then bad pixels inside the circle are flagged as
 *   good.
 * - If @a additive is true, then only good pixels inside the circle are
 *   modified.
 * - If @a anti_alias is true, then only good pixels in the range of the
 *   softened edge are modified (regardless of @a additive).
 * 
 * @par Error Handling:
 * The following errors can be set and returned:
 * - CPL_ERROR_NULL_INPUT: @a img_modified or @a centre_xy is NULL
 * - CPL_ERROR_INVALID_TYPE: @a img_modified is not of type CPL_TYPE_INT,
 *   CPL_TYPE_FLOAT or CPL_TYPE_DOUBLE
 * .
 * In the case of error, @a img_modified is not changed.
 */
/*----------------------------------------------------------------------------*/
cpl_error_code  clipm_priv_image_fill_circle(
                                            cpl_image       *img_modified,
                                            double          centre_xy[2],
                                            double          radius,
                                            double          brightness,
                                            int             anti_alias,
                                            int             additive)
{
    CLIPM_TRY
    {
        int                 wdw[2][2];
        cpl_type            type;
        int                 imsize[2],
                            wdwsize[2],
                            start[2];
        const void          *data_start;
        const cpl_binary    *badp_start;
        double              r_in,
                            r_out,
                            edge_width;
        int                 dim;
        
        CLIPM_TRY_CHECK_AUTOMSG(            img_modified != NULL,
                                            CPL_ERROR_NULL_INPUT);
        CLIPM_TRY_CHECK_AUTOMSG(            centre_xy != NULL,
                                            CPL_ERROR_NULL_INPUT);
        
        imsize[0] = cpl_image_get_size_x(img_modified);
        imsize[1] = cpl_image_get_size_y(img_modified);
        CLIPM_TRY_ASSERT_ERROR_STATE();
        
        if (radius < 0)
            CLIPM_TRY_EXIT();
        
        if (anti_alias)
        {
            edge_width = CPL_MATH_SQRT2;
            r_in = radius - edge_width / 2.0;
            r_out = radius + edge_width / 2.0;
        }
        else
        {
            edge_width = 0;
            r_out = r_in = radius;
        }
        
        for (dim = 0; dim < 2; dim++)
        {
            wdw[dim][0] = ceil(centre_xy[dim] - r_out);
            wdw[dim][1] = floor(centre_xy[dim] + r_out);
        }
        
        clipm_priv_checks_window_guarantee( *wdw,
                                            imsize[0],
                                            imsize[1],
                                            1);
        CLIPM_TRY_ASSERT_ERROR_STATE();
        
        clipm_priv_image_get_data_const(    img_modified,
                                            *wdw,
                                            0,
                                            imsize,
                                            wdwsize,
                                            start,
                                            &type,
                                            &data_start,
                                            &badp_start);
        CLIPM_TRY_ASSERT_ERROR_STATE();
        
/* --- throw in a macro to ease life ---------------------------------------- */
#define clipm_image_fill_circle_BODY(  TYPE) \
do { \
    int         x, \
                y; \
    TYPE        *data = (TYPE*)data_start; \
    TYPE        maxbright; \
    cpl_binary  *badp = (cpl_binary*)badp_start; \
    double      dx, \
                dy, \
                d, \
                w; \
    \
    if (type == CPL_TYPE_INT) \
        maxbright = clipm_math_round_d2i(brightness); \
    else \
        maxbright = brightness; \
    \
    for (   y = 0; \
            y < wdwsize[1]; \
            y++, data += imsize[0], badp += imsize[0]) \
    { \
        /* possible because we know this is lower-left */ \
        dy = (double)(y + start[1]+1) - centre_xy[1]; \
        dx = (double)(start[0]+1) - centre_xy[0]; \
        \
        if (additive) \
        { \
            for (   x = 0; x < wdwsize[0]; x++, dx += 1.0) \
                if (badp_start == NULL || badp[x] == CPL_BINARY_0) \
                { \
                    d = sqrt(dx*dx + dy*dy); \
                    if (d < r_out) \
                    { \
                        if (d <= r_in) \
                            data[x] += maxbright; \
                        else \
                        { \
                            w = (r_out - d)/edge_width; \
                            if (type == CPL_TYPE_INT) \
                                data[x] += clipm_math_round_d2i(w*brightness); \
                            else \
                                data[x] += w*brightness; \
                        } \
                    } \
                } \
        } \
        else \
        { \
            for (   x = 0; x < wdwsize[0]; x++, dx += 1.0) \
            { \
                d = sqrt(dx*dx + dy*dy); \
                if (d < r_out) \
                { \
                    if (d <= r_in) \
                    { \
                        data[x] = maxbright; \
                        if (badp_start != NULL) \
                            badp[x] = CPL_BINARY_0; \
                    } \
                    else if (badp_start == NULL || badp[x] == CPL_BINARY_0) \
                    { \
                        w = (r_out - d)/edge_width; \
                        if (type == CPL_TYPE_INT) \
                            data[x] = clipm_math_round_d2i( \
                                        w*brightness + (1.0-w)*data[x]); \
                        else \
                            data[x] = w*brightness + (1.0-w)*data[x]; \
                    } \
                    /* else do nothing */ \
                } \
            } \
        } \
    } \
} while (0)
/* -------------------------------------------------------------------------- */
        
        switch(type)
        {
            case    CPL_TYPE_INT:
                clipm_image_fill_circle_BODY(int);
                break;
            case    CPL_TYPE_FLOAT:
                clipm_image_fill_circle_BODY(float);
                break;
            case    CPL_TYPE_DOUBLE:
                clipm_image_fill_circle_BODY(double);
                break;
            default:
                CLIPM_TRY_EXIT_WITH_ERROR_MSG(
                                            CPL_ERROR_INVALID_TYPE,
                                            "img_modified",
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
 * @brief   Fill a polygon.
 * @param   image           The image to draw into
 * @param   points          2xN corner point matrix
 * @param   brightness      The polygon's content value
 * @param   anti_alias      Flag whether to soften the edge
 * @param   additive        Flag whether to add the circle on top of the
 *                          existing image content (@a additive != 0),
 *                          or whether to set the area to the @a brightness
 *                          value
 * @return  CPL error code
 * 
 * @par Principle:
 * - If @a additive is true, then the filled polygon will be added to the
 *   current image content, otherwise the polygon is filled with @a brightness.
 * - If @a anti_alias is true, then the edges are softened.
 * - Anti-alias values are computed numerically.
 * - If the type of @a image is CPL_TYPE_INT, then the resulting values
 *   are rounded.
 * 
 * @note
 * - Be aware, that a very narrow structures can also fit between pixels, so
 *   that none might be set, if @a anti_alias is false.
 * - The number of points @a N must be equal or greater than 3.
 * 
 * @par Bad Pixel Handling:
 * - If @a additive is false, then bad pixels inside the polygon are flagged as
 *   good.
 * - If @a additive is true, then only good pixels inside the polygon are
 *   modified.
 * - If @a anti_alias is true, then only good pixels in the range of the
 *   softened edge are modified (regardless of @a additive).
 * 
 * @par Error Handling:
 * The following errors can be set and returned:
 * - CPL_ERROR_NULL_INPUT: @a image or @a points is NULL
 * - CPL_ERROR_ILLEGAL_INPUT: the @a points matrix has an improper size
 * - CPL_ERROR_INVALID_TYPE: @a img_modified is not of type CPL_TYPE_INT,
 *   CPL_TYPE_FLOAT or CPL_TYPE_DOUBLE
 * .
 * In the case of error, @a image is not changed.
 */
/*----------------------------------------------------------------------------*/
cpl_error_code  clipm_priv_image_fill_polygon(
                                            cpl_image           *image,
                                            const cpl_matrix    *points,
                                            double              brightness,
                                            int                 anti_alias,
                                            int                 additive)
{
    cpl_matrix  *dimpos = NULL;
    cpl_mask    *edgemask = NULL;
    
    CLIPM_TRY
    {
        int         imsize[2],
                    wdwsize[2],
                    start[2];
        int         wdw[2][2];
        int         dim,
                    Np;
        cpl_type    type;
        const void  *data_start = NULL;
        const cpl_binary
                    *badp_start = NULL;
        const double
                    *xpdata,
                    *ypdata;
        
        CLIPM_TRY_CHECK_AUTOMSG(            image != NULL,
                                            CPL_ERROR_NULL_INPUT);
        CLIPM_TRY_CHECK_AUTOMSG(            points != NULL,
                                            CPL_ERROR_NULL_INPUT);
        CLIPM_TRY_CHECK_AUTOMSG(            cpl_matrix_get_nrow(points) == 2,
                                            CPL_ERROR_ILLEGAL_INPUT);
        CLIPM_TRY_CHECK_AUTOMSG(            cpl_matrix_get_ncol(points) >= 3,
                                            CPL_ERROR_ILLEGAL_INPUT);
        
        Np = cpl_matrix_get_ncol(points);
        
        for (dim = 0; dim < 2; dim++)
        {
            dimpos = cpl_matrix_extract(points, dim, 0, 1, 1, 1, Np);
            CLIPM_TRY_ASSERT_ERROR_STATE();
            
            wdw[dim][0] = floor(cpl_matrix_get_min(dimpos));
            wdw[dim][1] = ceil(cpl_matrix_get_max(dimpos));
            cpl_matrix_delete(dimpos);
            dimpos = NULL;
        }
        
        xpdata = cpl_matrix_get_data_const(points);
        ypdata = xpdata + Np;
        CLIPM_TRY_ASSERT_ERROR_STATE();
        
        clipm_priv_checks_window_guarantee_image(
                                            *wdw,
                                            image,
                                            1);
        CLIPM_TRY_ASSERT_ERROR_STATE();
        
        clipm_priv_image_get_data_const(    image,
                                            *wdw,
                                            1,
                                            imsize,
                                            wdwsize,
                                            start,
                                            &type,
                                            &data_start,
                                            &badp_start);
        CLIPM_TRY_ASSERT_ERROR_STATE();
        
        if (anti_alias)
        {
            edgemask = _clipm_priv_image_mark_polygon_edges(
                                            *wdw,
                                            Np,
                                            xpdata,
                                            ypdata);
            CLIPM_TRY_ASSERT_ERROR_STATE();
        }
        
/* --- throw in a macro to ease life ---------------------------------------- */
#define clipm_priv_image_fill_polygon_BODY(  TYPE) \
do { \
    int         xi, \
                yi; \
    TYPE        *data; \
    TYPE        maxbright; \
    cpl_binary  *badp; \
    const cpl_binary \
                *edgedata = NULL; \
    double      xpos, \
                ypos; \
    \
    if (type == CPL_TYPE_INT) \
        maxbright = clipm_math_round_d2i(brightness); \
    else \
        maxbright = brightness; \
    \
    data = (TYPE*)data_start; \
    badp = (cpl_binary*)badp_start; \
    if (edgemask != NULL) \
        edgedata = cpl_mask_get_data_const(edgemask); \
    CLIPM_TRY_ASSERT_ERROR_STATE(); \
    \
    ypos = start[1] + 1; \
    for (   yi = 0; yi < wdwsize[1]; \
            yi++, data += imsize[0], badp += imsize[0], ypos += 1.0) \
    { \
        xpos = start[0] + 1; \
        for (xi = 0; xi < wdwsize[0]; xi++, xpos += 1.0, edgedata++) \
        { \
            if (edgemask == NULL || *edgedata == CPL_BINARY_0) \
            { \
                if (_clipm_priv_image_is_point_in_polygon( \
                                            Np, \
                                            xpdata, ypdata, \
                                            xpos, ypos)) \
                { \
                    if (additive) \
                    { \
                        if (badp_start == NULL || badp[xi] == CPL_BINARY_0)\
                            data[xi] += maxbright; \
                    } \
                    else \
                    { \
                        data[xi] = maxbright; \
                        if (badp_start != NULL) \
                            badp[xi] = CPL_BINARY_0; \
                    } \
                } \
            } \
            else if (badp_start == NULL || badp[xi] == CPL_BINARY_0) \
            { \
                double  w; \
                w = _clipm_priv_image_get_pixel_area_in_polygon( \
                                            Np, \
                                            xpdata, ypdata, \
                                            xpos, ypos); \
                if (additive) \
                    w = w*brightness + (double)data[xi]; \
                else \
                    w = w*brightness + (1.0-w)*data[xi]; \
                if (type == CPL_TYPE_INT) \
                    data[xi] = clipm_math_round_d2i(w); \
                else \
                    data[xi] = w; \
            } \
        } \
    } \
} while (0)
/* -------------------------------------------------------------------------- */
        
        switch(type)
        {
            case    CPL_TYPE_INT:
                clipm_priv_image_fill_polygon_BODY(int);
                break;
            case    CPL_TYPE_FLOAT:
                clipm_priv_image_fill_polygon_BODY(float);
                break;
            case    CPL_TYPE_DOUBLE:
                clipm_priv_image_fill_polygon_BODY(double);
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
    
    cpl_matrix_delete(dimpos);
    cpl_mask_delete(edgemask);
    
    return CLIPM_ERROR_GET_NEW_SINCE_TRY();
}

/*----------------------------------------------------------------------------*/
/**
 * @brief   Fill a rectangle.
 * @param   image           The image to draw into
 * @param   centre_xy       Centre coordinate as (x,y) tuple
 * @param   angle           Orientation
 * @param   size_lw         Length and width as a data tuple, the first entry
 *                          (size_lw[0]) will point into the direction of
 *                          @a angle
 * @param   brightness      The polygon's content value
 * @param   anti_alias      Flag whether to soften the edge
 * @param   additive        Flag whether to add the circle on top of the
 *                          existing image content (@a additive != 0),
 *                          or whether to set the area to the @a brightness
 *                          value
 * @return  CPL error code
 * 
 * @par Principle:
 * - If @a additive is true, then the filled polygon will be added to the
 *   current image content, otherwise the polygon is filled with @a brightness.
 * - If @a anti_alias is true, then the edges are softened.
 * - Anti-alias values are computed numerically.
 * - If the type of @a image is CPL_TYPE_INT, then the resulting values
 *   are rounded.
 * 
 * @par Bad Pixel Handling:
 * - If @a additive is false, then bad pixels inside the rectangle are flagged
 *   as good.
 * - If @a additive is true, then only good pixels inside the rectangle are
 *   modified.
 * - If @a anti_alias is true, then only good pixels in the range of the
 *   softened edge are modified (regardless of @a additive).
 * 
 * @par Error Handling:
 * The following errors can be set and returned:
 * - CPL_ERROR_NULL_INPUT: any input pointer is NULL
 * - CPL_ERROR_INVALID_TYPE: @a img_modified is not of type CPL_TYPE_INT,
 *   CPL_TYPE_FLOAT or CPL_TYPE_DOUBLE
 * .
 * In the case of error, @a image is not changed.
 * 
 * @todo
 * - implement unit test
 */
/*----------------------------------------------------------------------------*/
cpl_error_code  clipm_priv_image_fill_rectangle(
                                            cpl_image           *image,
                                            const double        centre_xy[2],
                                            double              angle,
                                            const double        size_lw[2],
                                            double              brightness,
                                            int                 anti_alias,
                                            int                 additive)
{
    cpl_matrix  *rectpoly = NULL;
    
    CLIPM_TRY
    {
        CLIPM_TRY_CHECK_AUTOMSG(            image != NULL,
                                            CPL_ERROR_NULL_INPUT);
        CLIPM_TRY_CHECK_AUTOMSG(            centre_xy != NULL,
                                            CPL_ERROR_NULL_INPUT);
        CLIPM_TRY_CHECK_AUTOMSG(            size_lw != NULL,
                                            CPL_ERROR_NULL_INPUT);
        clipm_priv_checks_imtype_any(image, NULL);
        CLIPM_TRY_CHECK_ERROR_STATE();
        
        rectpoly = clipm_priv_matrix_create_corners_rectangle(
                                            centre_xy,
                                            angle,
                                            size_lw);
        CLIPM_TRY_ASSERT_ERROR_STATE();
        
        clipm_priv_image_fill_polygon(      image,
                                            rectpoly,
                                            brightness,
                                            anti_alias,
                                            additive);
        CLIPM_TRY_ASSERT_ERROR_STATE();
    }
    CLIPM_CATCH
    {
    }
    
    cpl_matrix_delete(rectpoly);
    
    return CLIPM_ERROR_GET_NEW_SINCE_TRY();
}

/*----------------------------------------------------------------------------*/
/**
 * @brief   Write an image to disk and print a message.
 * @param   img
 * @param   filename
 * @return  CPL error code
 * 
 * @par Principle:
 * - The image is written to disk.
 * - A message is printed.
 * - Bad pixels are filled with min-(max-min) for visualisation.
 * 
 * @note
 * This function is intended for debugging/hacking and unit testing purposes
 * only.
 * 
 * @par Error Handling:
 * The following error codes can be set and returned:
 * - CPL_ERROR_NULL_INPUT: @a img or @a filename is NULL
 * - CPL_ERROR_TYPE_MISMATCH: the image type is not supported
 * - CPL_ERROR_FILE_NOT_CREATED: the file could not be created
 * 
 * @todo
 * - implement unit test
 */
/*----------------------------------------------------------------------------*/
cpl_error_code  clipm_priv_image_save_debug(const cpl_image *img,
                                            const char      *filename)
{
    cpl_image   *tmpim = NULL;
    
    CLIPM_TRY
    {
        int         nbad = 0;
        
        CLIPM_TRY_CHECK_AUTOMSG(            img != NULL,
                                            CPL_ERROR_NULL_INPUT);
        CLIPM_TRY_CHECK_AUTOMSG(            filename != NULL,
                                            CPL_ERROR_NULL_INPUT);
        
        nbad = cpl_image_count_rejected(img);
        CLIPM_TRY_ASSERT_ERROR_STATE();
        
        if (nbad > 0)
        {
            int     sizex,
                    sizey;
            double  bval;
            
            tmpim = cpl_image_duplicate(img);
            CLIPM_TRY_ASSERT_ERROR_STATE();
            
            sizex = cpl_image_get_size_x(img);
            sizey = cpl_image_get_size_y(img);
            CLIPM_TRY_ASSERT_ERROR_STATE();
            
            if (nbad < sizex * sizey)
            {
                double  min,
                        max;
                min = cpl_image_get_min(tmpim);
                max = cpl_image_get_max(tmpim);
                
                if (min < max)
                    bval = 2 * min - max;
                else
                    bval = -1.0;
            }
            else
            {
                bval = -1.0;
            }
            
            cpl_image_fill_rejected(tmpim, bval);
            CLIPM_TRY_ASSERT_ERROR_STATE();
        }
        else
            tmpim = (cpl_image*)img;
        
        cpl_msg_info("", "> Writing debug image \'%s\'.", filename);
        
        cpl_image_save(                     tmpim,
                                            filename,
                                            CPL_BPP_IEEE_FLOAT,
                                            NULL,
                                            CPL_IO_DEFAULT);
        CLIPM_TRY_CHECK_ERROR_STATE();
    }
    CLIPM_CATCH
    {
    }
    
    if (tmpim != img)
        cpl_image_delete(tmpim);
    
    return CLIPM_ERROR_GET_NEW_SINCE_TRY();
}

/*----------------------------------------------------------------------------*/
/**
 * @brief   Delete a CPL image object and set the pointer to NULL.
 * @param   i   Pointer to image pointer
 * @return  Nothing
 * 
 * The following code is executed:
 * @code
    if (i != NULL)
    {
        cpl_image_delete(*i);  // checks for NULL pointer
        *i = NULL;
    }
 * @endcode
 * 
 * @par Error Handling:
 * No error can occur here.
 */
/*----------------------------------------------------------------------------*/
void            clipm_priv_image_null(      cpl_image       **i)
{
    if (i != NULL)
    {
        cpl_image_delete(*i);  /* checks for NULL pointer */
        *i = NULL;
    }
}

/**@}*/
