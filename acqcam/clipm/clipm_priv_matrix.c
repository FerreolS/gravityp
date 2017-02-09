/*********************************************************************
 * E.S.O. - VLT project
 *
 * "@(#) $Id: clipm_priv_matrix.c 177213 2008-12-05 15:33:51Z hlorch $"
 *
 * Private functions for handling matrices
 *
 * who       when        what
 * --------  ----------  ----------------------------------------------
 * hlorch    2007-02-07  created
 */

/**
 * @internal
 * @defgroup clipm_priv_matrix Matrix Handling
 * @ingroup internal_docs
 *
 * This module provides private matrix handling functions.
 *
 * @par Synopsis:
 * @code
#include "clipm_priv_matrix.h"
 * @endcode
 */
/**@{*/

/*-----------------------------------------------------------------------------
    Includes
 -----------------------------------------------------------------------------*/

#include "clipm_priv_matrix.h"

#include "clipm_math.h"
#include "clipm_priv_checks.h"
#include "clipm_compatibility_replacements.h"
#include "clipm_priv_error.h"

#include <math.h>

/*-----------------------------------------------------------------------------
    Implementation
 -----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
 * @brief   Extract an image window into a matrix, flip vertically.
 * @param   image       Input image (FITS convention)
 * @param   window_xxyy Coordinate buffer (FITS) of the form
 *                      {xa, xb, ya, yb}, can be NULL,
 *                      minimum/maximum order is irrelevant
 * @param   vflip       Flag whether to vertically flip the data buffer 
 * @return  The matrix
 * 
 * @par Principle:
 * - The image data (from a window) is copied into a matrix.
 * - If @a vflip is not zero, then the data are flipped
 *   vertically (in terms of indexing) to respect the FITS convention.
 * - Images can be of type CPL_TYPE_INT, CPL_TYPE_FLOAT or CPL_TYPE_DOUBLE.
 * 
 * @par Error Handling:
 * Possible error codes set in this function:
 * - CPL_ERROR_NULL_INPUT if image is NULL
 * - CPL_ERROR_INVALID_TYPE if the passed image type is not supported
 * - CPL_ERROR_ACCESS_OUT_OF_RANGE if a coordinate is out of the image range
 * 
 * @par Example:
 * With @a vflip != 0, the following image content
 * @code
    column  1   2   3   4
    row   __________________
    4   |   A   B   C   D
    3   |   E   F   G   H
    2   |   I   J   K   L
    1   |   M   N   O   P
 * @endcode
 * together with the window {2, 4, 2, 4} will result in the matrix:
 * @code
    column  0   1   2
    row   __________________
    0   |   B   C   D
    1   |   F   G   H
    2   |   J   K   L
 * @endcode
 * 
 * @todo
 * Remove when cpl_matrix_new_from_image_window() is available.
 */
/*----------------------------------------------------------------------------*/
cpl_matrix      *clipm_priv_matrix_new_from_image_window(
                                            const cpl_image     *image,
                                            const int           window_xxyy[4],
                                            const int           vflip)
{
    int         imsize[2], windowsize[2], windowstart[2]; /* x and y */
    cpl_type    type;
    cpl_matrix  *matrix = NULL;
    double      *mdata = NULL;

    CLIPM_TRY
    {
        /* check input, check ranges, get image/window size */
        clipm_priv_checks_window_image(     window_xxyy,
                                            image,
                                            1,
                                            imsize,
                                            windowsize,
                                            windowstart);
        CLIPM_TRY_CHECK_ERROR_STATE();

        mdata = cpl_malloc(windowsize[0] * windowsize[1] * sizeof(double));
        CLIPM_TRY_CHECK_ERROR_STATE();
        
/* --- throw in a macro to ease life ---------------------------------------- */
/* copy the image window into a double buffer, and subtract the mean */
#define clipm_COPY_WINDOW_TO_MATRIX( \
                                            TYPE, \
                                            image, imsize, windowsize, \
                                            buffer_start, \
                                            out_buffer, \
                                            flip) \
do { \
    const TYPE  *imdata; \
    double      *out_data; \
    int         x, \
                y, \
                vstep = 1; \
    imdata = cpl_image_get_data_const(image); \
    /* start at top image window border*/ \
    imdata += buffer_start[1]*imsize[0] + buffer_start[0]; \
    if (flip) \
    { \
        imdata += (windowsize[1]-1)*imsize[0]; \
        vstep *= -1; \
    } \
    out_data = out_buffer; \
    for (y = 0; y < windowsize[1]; y++) { \
        for (x = 0; x < windowsize[0]; x++) \
            *(out_data++) = imdata[x]; \
        imdata += vstep * imsize[0]; /* go up/down one image row */ \
    } \
} while (0)
/* -------------------------------------------------------------------------- */

        switch(type = cpl_image_get_type(image)) {
            case CPL_TYPE_INT: {
                clipm_COPY_WINDOW_TO_MATRIX(int,
                                            image, imsize, windowsize,
                                            windowstart,
                                            mdata,
                                            vflip);
                break;
            }
            case CPL_TYPE_FLOAT: {
                clipm_COPY_WINDOW_TO_MATRIX(float,
                                            image, imsize, windowsize,
                                            windowstart,
                                            mdata,
                                            vflip);
                break;
            }
            case CPL_TYPE_DOUBLE: {
                clipm_COPY_WINDOW_TO_MATRIX(double,
                                            image, imsize, windowsize,
                                            windowstart,
                                            mdata,
                                            vflip);
                break;
            }
            default:
                CLIPM_TRY_EXIT_WITH_ERROR_MSG(
                                            CPL_ERROR_INVALID_TYPE,
                                            "image",
                                            "must be int, float or double");
        }
        
        matrix = cpl_matrix_wrap(windowsize[1], windowsize[0], mdata);
        CLIPM_TRY_CHECK_ERROR_STATE();
    }
    CLIPM_CATCH
    {
        cpl_matrix_unwrap(matrix);
        matrix = NULL;
        if (mdata != NULL)
            cpl_free(mdata);
    }
    
    return matrix;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief   Round the elements of a matrix to the nearest integer.
 * @param   matrix  The matrix
 * @return  CPL error code
 * 
 * @par Principle:
 * Element by element, the elements are rounded.
 * 
 * @par Error Handling:
 * If @a matrix is NULL, CPL_ERROR_NULL_INPUT is returned.
 * 
 * @todo:
 * Remove when cpl_matrix_round() is available.
 */
/*----------------------------------------------------------------------------*/
cpl_error_code  clipm_priv_matrix_round(    cpl_matrix          *matrix)
{
    CLIPM_TRY
    {
        double *ddata;
        int n, nrows, ncols;

        CLIPM_TRY_EXIT_IFN(
            ddata = cpl_matrix_get_data(matrix));

        nrows = cpl_matrix_get_nrow(matrix);
        ncols = cpl_matrix_get_ncol(matrix);
        CLIPM_TRY_CHECK_ERROR_STATE();

        for (n = 0; n < ncols * nrows; n++)
        {
            *ddata = clipm_math_round_d2i(*ddata);
            ddata++;
        }
    }
    CLIPM_CATCH
    {
    }
    
    return CLIPM_ERROR_GET_NEW_SINCE_TRY();
}

/*----------------------------------------------------------------------------*/
/**
 * @brief   Get the medians of the values in the matrix rows respectively.
 * @param   matrix  N x M input matrix
 * @return  N x 1 output matrix, NULL on error
 * 
 * @par Principle:
 * The output matrix has one column and the same number of rows as the input
 * matrix, with each column containing the median of all values in the
 * same row in the input matrix.
 * 
 * @par Error Handling:
 * If @a matrix is NULL, CPL_ERROR_NULL_INPUT is set and NULL is returned.
 */
/*----------------------------------------------------------------------------*/
cpl_matrix      *clipm_priv_matrix_get_median_rows(
                                            const cpl_matrix    *matrix)
{
    cpl_matrix      *out = NULL;
    const cpl_vector
                    *v = NULL;
    int             nrows, ncols;
    const double    *indata;
    double          *outdata;
    
    CLIPM_TRY
    {
        int n;
        
        CLIPM_TRY_EXIT_IFN(
            indata = cpl_matrix_get_data_const(matrix));
        
        nrows = cpl_matrix_get_nrow(matrix);
        ncols = cpl_matrix_get_ncol(matrix);
        
        CLIPM_TRY_EXIT_IFN(
            out = cpl_matrix_new(nrows, 1));
        outdata = cpl_matrix_get_data(out);
        
        for (n = 0; n < nrows; n++)
        {
            CLIPM_TRY_EXIT_IFN(
                v = cpl_vector_wrap(ncols, (double*)indata + n*ncols));
            *(outdata++) = cpl_vector_get_median_const(v);
            cpl_vector_unwrap((cpl_vector*)v);
        }
    }
    CLIPM_CATCH
    {
        cpl_matrix_delete(out);
        out = NULL;
        cpl_vector_unwrap((cpl_vector*)v);
    }
    
    return out;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief   Get the means of the values in the matrix rows respectively.
 * @param   matrix  N x M input matrix
 * @return  N x 1 output matrix, NULL on error
 * 
 * @par Principle:
 * The output matrix has one column and the same number of rows as the input
 * matrix, with each column containing the mean of all values in the
 * same row in the input matrix.
 * 
 * @par Error Handling:
 * If @a matrix is NULL, CPL_ERROR_NULL_INPUT is set and NULL is returned.
 */
/*----------------------------------------------------------------------------*/
cpl_matrix      *clipm_priv_matrix_get_mean_rows(
                                            const cpl_matrix    *matrix)
{
    cpl_matrix      *out = NULL;
    const cpl_vector
                    *v = NULL;
    int             nrows, ncols;
    const double    *indata;
    double          *outdata;
    
    CLIPM_TRY
    {
        int n;
        
        CLIPM_TRY_EXIT_IFN(
            indata = cpl_matrix_get_data_const(matrix));
        
        nrows = cpl_matrix_get_nrow(matrix);
        ncols = cpl_matrix_get_ncol(matrix);
        
        CLIPM_TRY_EXIT_IFN(
            out = cpl_matrix_new(nrows, 1));
        outdata = cpl_matrix_get_data(out);
        
        for (n = 0; n < nrows; n++)
        {
            CLIPM_TRY_EXIT_IFN(
                v = cpl_vector_wrap(ncols, (double*)indata + n*ncols));
            *(outdata++) = cpl_vector_get_mean(v);
            cpl_vector_unwrap((cpl_vector*)v);
        }
    }
    CLIPM_CATCH
    {
        cpl_matrix_delete(out);
        out = NULL;
        cpl_vector_unwrap((cpl_vector*)v);
    }
    
    return out;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief   Return selected columns in a new matrix.
 * @param   matrix          N x M input matrix
 * @param   selection       Integer buffer of size M containing the selection
 *                          flags
 * @param   select_nonzero  Flag to either select columns with a selection flag
 *                          equal to zero or different from zero
 * @return  N x nselected output matrix, NULL on error
 * 
 * @par Principle:
 * If @a select_nonzeros is 0, then all columns with @a selection[column]==0
 * are selected, otherwise all columns with @a selection[column]!=0 are
 * selected.
 * 
 * @par Error Handling:
 * The following error cases may occur, and the corresponding error codes set:
 * - CPL_ERROR_NULL_INPUT: @a matrix is NULL
 * - CPL_ERROR_DATA_NOT_FOUND: no columns are selected
 * .
 * In the case of error, NULL is returned.
 */
/*----------------------------------------------------------------------------*/
cpl_matrix      *clipm_priv_matrix_select_cols(
                                            const cpl_matrix    *matrix,
                                            const int           *selection,
                                            int                 select_nonzero)
{
    int             nrows, ncols, nselected = 0, n, r;
    cpl_matrix      *result = NULL;
    const double    *indata;
    double          *rdata;
    
    CLIPM_TRY
    {
        nrows = cpl_matrix_get_nrow(matrix);
        ncols = cpl_matrix_get_ncol(matrix);
        CLIPM_TRY_CHECK_ERROR_STATE();
        if (selection == NULL)
            CLIPM_TRY_EXIT_WITH_ERROR(  CPL_ERROR_NULL_INPUT);
        
        for (n = 0; n < ncols; n++)
        {
            if ((selection[n] && select_nonzero) ||
                    (!selection[n] && !select_nonzero))
                nselected++;
        }
        
        if (nselected == 0)
            CLIPM_TRY_EXIT_WITH_ERROR_MSG(  CPL_ERROR_DATA_NOT_FOUND,
                                            "",
                                            "selection is empty");
        
        CLIPM_TRY_EXIT_IFN(
            result = cpl_matrix_new(nrows, nselected));
        
        indata = cpl_matrix_get_data_const(matrix);
        rdata = cpl_matrix_get_data(result);
        CLIPM_TRY_CHECK_ERROR_STATE();
        
        for (n = 0; n < ncols; n++)
        {
            if ((selection[n] && select_nonzero) ||
                (!selection[n] && !select_nonzero))
            {
                for (r = 0; r < nrows; r++)
                    rdata[r*nselected] = indata[r*ncols];
                
                rdata++;
            }
            indata++;
        }

    }
    CLIPM_CATCH
    {
        cpl_matrix_delete(result);
        result = NULL;
    }
    
    return result;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief   Copy the content of a vector into a matrix column (or a part of it).
 * @param   matrix      Modified matrix
 * @param   vector      Input vector
 * @param   col         Matrix column index
 * @param   start_row   Index of the starting matrix row
 * @return  CPL error code
 * 
 * @par Principle:
 * A @a vector is copied into a column of the provided @a matrix, starting at
 * the matrix row @a start_row. Only the number of vector elements is copied.
 * @a start_row plus the @a vector's size must be less than or equal to the
 * number of @a matrix rows.
 * 
 * @par Error Handling:
 * This function can return the following error codes:
 * - CPL_ERROR_NULL_INPUT: @a matrix or @a vector is NULL
 * - CPL_ERROR_ACCESS_OUT_OF_RANGE:
 *   - @a col < 0, or @a col >= the number of
 *     columns of the @a matrix, or
 *   - @a start_row < 0, or
 *   - @a start_row + the @a vector's size > the number of @a matrix rows
 * .
 * In the case of error, the @a matrix is not modified.
 */
/*----------------------------------------------------------------------------*/
cpl_error_code  clipm_priv_matrix_copy_col_vector(
                                            cpl_matrix          *matrix,
                                            const cpl_vector    *vector,
                                            int                 col,
                                            int                 start_row)
{
    CLIPM_TRY
    {
        int     nrows,
                ncols,
                vectorsize,
                n;
        double  *mdata;
        const double
                *vdata;
        
        mdata = cpl_matrix_get_data(matrix);
        vdata = cpl_vector_get_data_const(vector);
        CLIPM_TRY_CHECK_ERROR_STATE();
        
        nrows = cpl_matrix_get_nrow(matrix);
        ncols = cpl_matrix_get_ncol(matrix);
        vectorsize = cpl_vector_get_size(vector);
        
        CLIPM_TRY_CHECK_AUTOMSG(
            col >= 0 && col < ncols,
            CPL_ERROR_ACCESS_OUT_OF_RANGE);
        CLIPM_TRY_CHECK_AUTOMSG(
            start_row >= 0 && start_row <= nrows - vectorsize,
            CPL_ERROR_ACCESS_OUT_OF_RANGE);
        
        mdata += col + start_row*ncols;
        for (n = 0; n < vectorsize; n++)
        {
            *mdata = *(vdata++);
            mdata += ncols;
        }
    }
    CLIPM_CATCH
    {
    }
    
    return CLIPM_ERROR_GET_NEW_SINCE_TRY();
}

/*----------------------------------------------------------------------------*/
/**
 * @brief   Copy the content of a matrix column (or a part of it) into a vector.
 * @param   vector      Modified vector
 * @param   matrix      Input matrix
 * @param   col         Matrix column index
 * @param   start_row   Index of the starting matrix row
 * @return  CPL error code
 * 
 * @par Principle:
 * A matrix column, starting at row @a start_row, is copied into the provided
 * @a vector. Only the number of vector elements is regarded. @a start_row
 * plus the @a vector's size must be less than or equal to the number of
 * @a matrix rows.
 * 
 * @par Error Handling:
 * This function can return the following error codes:
 * - CPL_ERROR_NULL_INPUT: @a vector or @a matrix is NULL
 * - CPL_ERROR_ACCESS_OUT_OF_RANGE:
 *   - @a col < 0, or @a col >= the number of
 *     columns of the @a matrix, or
 *   - @a start_row < 0, or
 *   - @a start_row + the @a vector's size > the number of @a matrix rows
 * .
 * In the case of error, the @a vector is not modified.
 */
/*----------------------------------------------------------------------------*/
cpl_error_code  clipm_priv_matrix_copy_vector_col(
                                            cpl_vector          *vector,
                                            const cpl_matrix    *matrix,
                                            int                 col,
                                            int                 start_row)
{
    CLIPM_TRY
    {
        int     nrows,
                ncols,
                vectorsize,
                n;
        const double
                *mdata;
        double  *vdata;
        
        mdata = cpl_matrix_get_data_const(matrix);
        vdata = cpl_vector_get_data(vector);
        CLIPM_TRY_CHECK_ERROR_STATE();
        
        nrows = cpl_matrix_get_nrow(matrix);
        ncols = cpl_matrix_get_ncol(matrix);
        vectorsize = cpl_vector_get_size(vector);
        
        CLIPM_TRY_CHECK_AUTOMSG(
            col >= 0 && col < ncols,
            CPL_ERROR_ACCESS_OUT_OF_RANGE);
        CLIPM_TRY_CHECK_AUTOMSG(
            start_row >= 0 && start_row <= nrows - vectorsize,
            CPL_ERROR_ACCESS_OUT_OF_RANGE);
        
        mdata += col + start_row*ncols;
        for (n = 0; n < vectorsize; n++)
        {
            *(vdata++) = *mdata;
            mdata += ncols;
        }
    }
    CLIPM_CATCH
    {
    }
    
    return CLIPM_ERROR_GET_NEW_SINCE_TRY();
}

/*----------------------------------------------------------------------------*/
/**
 * @brief   Create a 2d rotation matrix.
 * @param   angle       Rotation angle
 * @param   centre_xy   (Optional) rotation centre (double buffer of size 2),
 *                      if it is NULL, then (0, 0) is assumed as centre
 * @param   transform   (Optional output) transformation matrix
 * @param   shift       (Optional output) shift vector
 * @return  CPL error code
 * 
 * @par Error Handling:
 * No error should occur in this function.
 */
/*----------------------------------------------------------------------------*/
cpl_error_code  clipm_priv_matrix_transform_create_rot2d(
                                            double              angle,
                                            const double        centre_xy[2],
                                            cpl_matrix          **transform,
                                            cpl_matrix          **shift)
{
    cpl_matrix  *tmp_transform = NULL,
                *tmp = NULL;
    
    CLIPM_TRY
    {
        /* init output */
        clipm_priv_matrix_null(transform);
        clipm_priv_matrix_null(shift);
        
        if (transform == NULL && shift == NULL)
            CLIPM_TRY_EXIT();
        if (transform == NULL)
            transform = &tmp_transform;
        
        *transform = cpl_matrix_new(2, 2);
        
        cpl_matrix_set(*transform, 0, 0, cos(angle));
        cpl_matrix_set(*transform, 0, 1, -sin(angle));
        cpl_matrix_set(*transform, 1, 0, sin(angle));
        cpl_matrix_set(*transform, 1, 1, cos(angle));
        CLIPM_TRY_ASSERT_ERROR_STATE();
        
        if (shift == NULL)
            CLIPM_TRY_EXIT();
        
        *shift = cpl_matrix_new(2, 1);
        if (centre_xy != NULL)
        {
            cpl_matrix_set(*shift, 0, 0, centre_xy[0]);
            cpl_matrix_set(*shift, 1, 0, centre_xy[1]);
        }
        CLIPM_TRY_ASSERT_ERROR_STATE();
        
        tmp = cpl_matrix_product_create(*transform, *shift);
        CLIPM_TRY_ASSERT_ERROR_STATE();
        
        cpl_matrix_subtract(*shift, tmp);
        CLIPM_TRY_ASSERT_ERROR_STATE();
    }
    CLIPM_CATCH
    {
        clipm_priv_matrix_null(transform);
        clipm_priv_matrix_null(shift);
    }
    
    cpl_matrix_delete(tmp_transform);
    cpl_matrix_delete(tmp);
    
    return CLIPM_ERROR_GET_NEW_SINCE_TRY();
}

/*----------------------------------------------------------------------------*/
/**
 * @brief   Compute the inverse transform.
 * @param   transf      Input transformation matrix
 * @param   shift       Input shift vector
 * @param   inv_transf  Output transformation matrix
 * @param   inv_shift   Output shift vector
 * @return  CPL error code
 * 
 * @par Principle:
 * The following expression "out = transf * in + shift" is inverted so that
 * "in = inv_transf * out + inv_shift".
 * 
 * @par Error Handling:
 * The following error codes can be set and returned:
 * - CPL_ERROR_NULL_INPUT: any input is NULL
 * - CPL_ERROR_ILLEGAL_INPUT: @a transf is not a square matrix
 * - CPL_ERROR_SINGULAR_MATRIX: @a transf cannot be inverted
 * - CPL_ERROR_INCOMPATIBLE_INPUT: the number of columns of @a transf is not
 *   equal to the number of rows of @a shift
 * 
 * @todo
 * - implement unit test
 */
/*----------------------------------------------------------------------------*/
cpl_error_code  clipm_priv_matrix_transform_invert(
                                            const cpl_matrix    *transf,
                                            const cpl_matrix    *shift,
                                            cpl_matrix          **inv_transf,
                                            cpl_matrix          **inv_shift)
{
    CLIPM_TRY
    {
        /* init output */
        clipm_priv_matrix_null(inv_transf);
        clipm_priv_matrix_null(inv_shift);
        
        CLIPM_TRY_CHECK_AUTOMSG(            transf != NULL,
                                            CPL_ERROR_NULL_INPUT);
        CLIPM_TRY_CHECK_AUTOMSG(            shift != NULL,
                                            CPL_ERROR_NULL_INPUT);
        CLIPM_TRY_CHECK_AUTOMSG(            inv_transf != NULL,
                                            CPL_ERROR_NULL_INPUT);
        CLIPM_TRY_CHECK_AUTOMSG(            inv_shift != NULL,
                                            CPL_ERROR_NULL_INPUT);
        
        *inv_transf = cpl_matrix_invert_create(transf);
        CLIPM_ERROR_SET_MSG_IF_CODE(        CPL_ERROR_SINGULAR_MATRIX,
                                            "",
                                            "transform not invertible");
        CLIPM_TRY_CHECK_ERROR_STATE();
        
        *inv_shift = cpl_matrix_product_create(
                                            *inv_transf,
                                            shift);
        CLIPM_TRY_CHECK_ERROR_STATE();
        
        cpl_matrix_multiply_scalar(         *inv_shift,
                                            -1.0);
        CLIPM_TRY_ASSERT_ERROR_STATE();
    }
    CLIPM_CATCH
    {
        clipm_priv_matrix_null(inv_transf);
        clipm_priv_matrix_null(inv_shift);
    }
    
    return CLIPM_ERROR_GET_NEW_SINCE_TRY();
}

/*----------------------------------------------------------------------------*/
/**
 * @brief   Transform a pointlist.
 * @param   points          ND x NP (number of dimensions, times number of
 *                          points) matrix containing the point coordinates
 * @param   transformation  The ND x ND square transformation matrix,
 *                          can be NULL
 * @param   shift           The ND x 1 shift matrix (vector),
 *                          can be NULL
 * @return  The transformed points, NULL on error
 * 
 * @par Transformation:
 * The transformation is done by the following formula (where K denotes the
 * column of @a points):
 * \f[ point_K' = transformation \times point_K + shift \f]
 * 
 * @par Parameter Details:
 * - The matrices @a points and @a shift are expected in the following forms:\n
 *   \f$ points = \left( \begin{array}{lllc}
      x_0 & x_1 & \cdots & x_{NP-1} \\
      y_0 & y_1 & \cdots & y_{NP-1} \\
      z_0 & z_1 & \cdots & z_{NP-1} \\
      \vdots & \vdots & & \vdots
     \end{array} \right) \f$, and
     \f$shift = \left( \begin{array}{c}
     s_x\\
     s_y\\
     s_z\\
     \vdots
     \end{array} \right)\f$.
 * 
 * @par Constraints:
 * - The number of dimensions (ND) can be any > 0.
 * - The number of points (NP) can be any > 0.
 * - The @a transformation matrix must be square (ND x ND)
 * - The @a shift matrix represents a vector and must be of size ND x 1
 * - If @a transformation and @a shift are both NULL, then @a points is just
 *   copied.
 * 
 * @par Error Handling:
 * - CPL_ERROR_NULL_INPUT if @a points is NULL
 * - CPL_ERROR_INCOMPATIBLE_INPUT if the sizes of the input parameters do not
 *   fit together
 * - CPL_ERROR_ILLEGAL_INPUT:
 *   - if @a shift has more than one column
 *   - if @a transformation is not square
 */
/*----------------------------------------------------------------------------*/
cpl_matrix      *clipm_priv_matrix_transform_points(
                                            const cpl_matrix    *points,
                                            const cpl_matrix    *transformation,
                                            const cpl_matrix    *shift)
{
    cpl_matrix  *result = NULL;
    
    CLIPM_TRY
    {
        if (transformation != NULL)
        {
            CLIPM_TRY_CHECK(
                cpl_matrix_get_ncol(transformation) ==
                    cpl_matrix_get_nrow(transformation),
                CPL_ERROR_ILLEGAL_INPUT,
                "transformation matrix",
                "must be square");
            result = cpl_matrix_product_create(transformation, points);
            CLIPM_ERROR_SET_MSG_IF_CODE(
                CPL_ERROR_INCOMPATIBLE_INPUT,
                "points matrix",
                "nr of rows/dimensions must equal "
                "nr of columns of transf. matrix");
        }
        else
            result = cpl_matrix_duplicate(points);
        CLIPM_TRY_CHECK_ERROR_STATE();
        
        if (shift != NULL)
        {
            double          *pdata;
            const double    *sdata;
            int             p, dim, ndims, npoints;
            
            CLIPM_TRY_EXIT_IFN(
                ndims = cpl_matrix_get_nrow(points));
            npoints = cpl_matrix_get_ncol(points);
            CLIPM_TRY_CHECK(
                cpl_matrix_get_nrow(shift) == ndims,
                CPL_ERROR_INCOMPATIBLE_INPUT,
                "points, shift",
                "nr of dimensions (matrix rows) different");
            CLIPM_TRY_CHECK(
                cpl_matrix_get_ncol(shift) == 1,
                CPL_ERROR_ILLEGAL_INPUT,
                "shift",
                "must have 1 column");
            
            pdata = cpl_matrix_get_data(result);
            sdata = cpl_matrix_get_data_const(shift);
            
            for (dim = 0; dim < ndims; dim++)
            {
                for (p = 0; p < npoints; p++)
                    *(pdata++) += *sdata;
                sdata++;
            }
        }
    }
    CLIPM_CATCH
    {
        clipm_priv_matrix_null(&result);
    }
    
    return result;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief   Create a matrix containing the corner points of a rectangle.
 * @param   centre_xy   Centre coordinate as (x,y) tuple
 * @param   angle       Orientation
 * @param   size_lw     Length and width as a data tuple, the first entry
 *                      (size_lw[0]) will point into the direction of @a angle
 * @return  The matrix of size 2x4, NULL in the case of error
 * 
 * The output matrix will contain (x,y) tuples in 4 columns, respectively.
 * 
 * @par Error Handling:
 * - CPL_ERROR_NULL_INPUT: if any input pointer is NULL
 * 
 * @todo
 * - implement unit test
 */
/*----------------------------------------------------------------------------*/
cpl_matrix      *clipm_priv_matrix_create_corners_rectangle(
                                            const double        centre_xy[2],
                                            double              angle,
                                            const double        size_lw[2])
{
    cpl_matrix  *rectpoly_flat = NULL,
                *out_rectpoly_rot = NULL,
                *transform = NULL,
                *shift = NULL;
    
    CLIPM_TRY
    {
        CLIPM_TRY_CHECK_AUTOMSG(            centre_xy != NULL,
                                            CPL_ERROR_NULL_INPUT);
        CLIPM_TRY_CHECK_AUTOMSG(            size_lw != NULL,
                                            CPL_ERROR_NULL_INPUT);
        
        rectpoly_flat = cpl_matrix_new(2, 4);
        CLIPM_TRY_ASSERT_ERROR_STATE();
        
        shift = cpl_matrix_wrap(2, 1, (double*)centre_xy);
        CLIPM_TRY_ASSERT_ERROR_STATE();
        
        cpl_matrix_set(rectpoly_flat, 0, 0, -size_lw[0]/2.0);
        cpl_matrix_set(rectpoly_flat, 1, 0, -size_lw[1]/2.0);
        cpl_matrix_set(rectpoly_flat, 0, 1, -size_lw[0]/2.0);
        cpl_matrix_set(rectpoly_flat, 1, 1, size_lw[1]/2.0);
        cpl_matrix_set(rectpoly_flat, 0, 2, size_lw[0]/2.0);
        cpl_matrix_set(rectpoly_flat, 1, 2, size_lw[1]/2.0);
        cpl_matrix_set(rectpoly_flat, 0, 3, size_lw[0]/2.0);
        cpl_matrix_set(rectpoly_flat, 1, 3, -size_lw[1]/2.0);
        CLIPM_TRY_ASSERT_ERROR_STATE();
        
        clipm_priv_matrix_transform_create_rot2d(
                                            angle,
                                            NULL,
                                            &transform,
                                            NULL);
        CLIPM_TRY_ASSERT_ERROR_STATE();
        
        out_rectpoly_rot = clipm_priv_matrix_transform_points(
                                            rectpoly_flat,
                                            transform,
                                            shift);
        CLIPM_TRY_ASSERT_ERROR_STATE();
    }
    CLIPM_CATCH
    {
        clipm_priv_matrix_null(&out_rectpoly_rot);
    }
    
    cpl_matrix_delete(rectpoly_flat);
    cpl_matrix_delete(transform);
    
    cpl_matrix_unwrap(shift);
    
    return out_rectpoly_rot;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief   Delete a CPL matrix object and set the pointer to NULL.
 * @param   m   Pointer to matrix pointer
 * @return  Nothing
 * 
 * The following code is executed:
 * @code
    if (m != NULL)
    {
        cpl_matrix_delete(*m);  // checks for NULL pointer
        *m = NULL;
    }
 * @endcode
 * 
 * @par Error Handling:
 * No error can occur here.
 */
/*----------------------------------------------------------------------------*/
void            clipm_priv_matrix_null(     cpl_matrix          **m)
{
    if (m != NULL)
    {
        cpl_matrix_delete(*m);  /* checks for NULL pointer */
        *m = NULL;
    }
}

/**@}*/
