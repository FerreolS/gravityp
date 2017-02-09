/*********************************************************************
 * E.S.O. - VLT project
 *
 * "@(#) $Id: clipm_priv_vector.c 176658 2008-11-26 18:02:24Z hlorch $"
 *
 * Private functions for handling vectors
 *
 * who       when        what
 * --------  ----------  ----------------------------------------------
 * hlorch    2007-04-27  created
 */

/**
 * @internal
 * @defgroup clipm_priv_vector Vector Handling
 * @ingroup internal_docs
 *
 * This module provides private vector handling functions.
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

#include "clipm_priv_vector.h"

#include "clipm_math.h"
#include "clipm_compatibility_replacements.h"
#include "clipm_priv_error.h"

#include <string.h> /* memcpy() */

/*-----------------------------------------------------------------------------
    Implementation
 -----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
 * @brief   Return minimum and corresponding index of a vector.
 * @param   v       Input vector
 * @param   index   Output index, where the minimum has been found,
 *                  can be NULL, returns -1 in the case of error
 * @return  The minimum value, 0 in the case of error
 * 
 * @par Error Handling:
 * Possible error code set by this function:
 * - CPL_ERROR_NULL_INPUT: @a v is NULL
 */
/*----------------------------------------------------------------------------*/
double          clipm_priv_vector_get_min(  const cpl_vector    *v,
                                            int                 *index)
{
    double  min = 0;
    
    CLIPM_TRY
    {
        const double
                *data;
        int     size,
                n,
                ndx = 0;
        
        if (index != NULL)
            *index = -1;
        
        CLIPM_TRY_EXIT_IFN(
            data = cpl_vector_get_data_const(v));
        size = cpl_vector_get_size(v);
        
        min = *(data++);
        for (n = 1; n < size; n++)
        {
            if (min > *data)
            {
                min = *data;
                ndx = n;
            }
            data++;
        }
        
        if (index != NULL)
            *index = ndx;
    }
    CLIPM_CATCH
    {
    }
    
    return min;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief   Return maximum and corresponding index of a vector.
 * @param   v       Input vector
 * @param   index   Output index, where the maximum has been found,
 *                  can be NULL, returns -1 in the case of error
 * @return  The maximum value, 0 in the case of error
 * 
 * @par Error Handling:
 * Possible error code set by this function:
 * - CPL_ERROR_NULL_INPUT: @a v is NULL
 */
/*----------------------------------------------------------------------------*/
double          clipm_priv_vector_get_max(  const cpl_vector    *v,
                                            int                 *index)
{
    double  max = 0;
    
    CLIPM_TRY
    {
        const double
                *data;
        int     size,
                n,
                ndx = 0;
        
        if (index != NULL)
            *index = -1;
        
        CLIPM_TRY_EXIT_IFN(
            data = cpl_vector_get_data_const(v));
        size = cpl_vector_get_size(v);
        
        max = *(data++);
        for (n = 1; n < size; n++)
        {
            if (max < *data)
            {
                max = *data;
                ndx = n;
            }
            data++;
        }
        
        if (index != NULL)
            *index = ndx;
    }
    CLIPM_CATCH
    {
    }
    
    return max;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief   Prepend and append @em borderwidth 0-elements to @a v,
 *          respectively.
 * @param   v           Input vector
 * @param   prepend_nr  Number of elements that are going to be prepended at
 *                      the beginning of @a v
 * @param   append_nr   Number of elements that are going to be appended to
 *                      the end of @a v
 * @return  Resulting vector, NULL in the case of error
 * 
 * @par Error Handling:
 * The following error codes can be set:
 * - CPL_ERROR_NULL_INPUT: @a v is NULL
 * - CPL_ERROR_ILLEGAL_INPUT: @a prepend_nr < 0 or @a append_nr < 0
 */
/*----------------------------------------------------------------------------*/
cpl_vector      *clipm_priv_vector_expand(  const cpl_vector    *v,
                                            int                 prepend_nr,
                                            int                 append_nr)
{
    cpl_vector  *result = NULL;
    
    CLIPM_TRY
    {
        int     size,
                n;
        const double
                *vdata;
        double  *rdata;
        
        size = cpl_vector_get_size(v);
        CLIPM_TRY_CHECK_ERROR_STATE();
        
        CLIPM_TRY_CHECK_AUTOMSG(
            prepend_nr >= 0,
            CPL_ERROR_ILLEGAL_INPUT);
        CLIPM_TRY_CHECK_AUTOMSG(
            append_nr >= 0,
            CPL_ERROR_ILLEGAL_INPUT);
        
        result = cpl_vector_new(size + prepend_nr + append_nr);
        CLIPM_TRY_CHECK_ERROR_STATE();
        
        vdata = cpl_vector_get_data_const(v);
        rdata = cpl_vector_get_data(result);
        CLIPM_TRY_CHECK_ERROR_STATE();
        
        for (n = 0; n < prepend_nr; n++)
            *(rdata++) = 0;
        memcpy(                             rdata,
                                            vdata,
                                            size * sizeof(*vdata));
        rdata += size;
        for (n = 0; n < append_nr; n++)
            *(rdata++) = 0;
    }
    CLIPM_CATCH
    {
        clipm_priv_vector_null(&result);
    }
    
    return result;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief   Integrate a signal stored in a vector.
 * @param   v       Input vector
 * @return  The minimum value, 0 in the case of error
 * 
 * @par Error Handling:
 * Possible error code set by this function:
 * - CPL_ERROR_NULL_INPUT: @a v is NULL
 */
/*----------------------------------------------------------------------------*/
cpl_vector      *clipm_priv_vector_integrate(
                                            const cpl_vector    *v)
{
    cpl_vector  *result = NULL;
    
    CLIPM_TRY
    {
        int     n,
                size;
        const double
                *indata;
        double  *outdata;
        
        CLIPM_TRY_EXIT_IFN(
            indata = cpl_vector_get_data_const(v));
        
        size = cpl_vector_get_size(v);
        CLIPM_TRY_EXIT_IFN(
            result = cpl_vector_new(size));
        outdata = cpl_vector_get_data(result);
        
        *(outdata) = *(indata++);
        for (n = 1; n < size; n++)
        {
            outdata[1] = *outdata + *(indata++);
            outdata++;
        }
    }
    CLIPM_CATCH
    {
        clipm_priv_vector_null(&result);
    }
    
    return result;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief   Differentiate a signal stored in a vector.
 * @param   v       Input vector
 * @return  The minimum value, 0 in the case of error
 * 
 * @par Error Handling:
 * Possible error code set by this function:
 * - CPL_ERROR_NULL_INPUT: @a v is NULL
 */
/*----------------------------------------------------------------------------*/
cpl_vector      *clipm_priv_vector_differentiate(
                                            const cpl_vector    *v)
{
    cpl_vector  *result = NULL;
    
    CLIPM_TRY
    {
        int     n,
                size;
        const double
                *indata;
        double  *outdata;
        
        CLIPM_TRY_EXIT_IFN(
            indata = cpl_vector_get_data_const(v));
        
        size = cpl_vector_get_size(v);
        CLIPM_TRY_EXIT_IFN(
            result = cpl_vector_new(size));
        outdata = cpl_vector_get_data(result);
        
        *(outdata++) = 0;
        indata++;
        for (n = 1; n < size; n++)
        {
            *(outdata++) = indata[0] - indata[-1];
            indata++;
        }
        /**outdata = indata[0] - indata[-1];*/
    }
    CLIPM_CATCH
    {
        clipm_priv_vector_null(&result);
    }
    
    return result;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief   Get the n-dimensional normal vector to n-1 direction vectors.
 * @param   dv  Array of n-1 direction vectors
 * @return  Normal vector, NULL in the case of error
 * 
 * @par Overview:
 * - The number of vectors in @a dv must be n-1, where n is the length of the
 *   vectors. All vectors must have the same length.
 * - The number of dimensions must be at least 2.
 * - The number of dimensions @em n is determined from the first vector in
 *   @a dv.
 * 
 * @par Mathematical Meaning:
 * - The result vector is also referred to in literature as
 *   \f$n = dv_1 \times dv_2 \times \ldots \times dv_{ndims-1}\f$.
 * - The norm of the result vector equals the (n-1)-dimensional volume of the
 *   parallelepiped spanned by the input vectors. This means, parallel input
 *   vectors will cause a result vector of zero length.
 * 
 * @note
 * If the size of @a dv is less than n-1, the function will crash the
 * application!
 * 
 * @par Error Handling:
 * The following error codes can be set by this function:
 * - CPL_ERROR_NULL_INPUT: @a dv == NULL, or any of the entries of @a dv
 *                          is NULL
 * - CPL_ERROR_INCOMPATIBLE_INPUT: not all entries (vectors) of @a dv have
 *                          the same size
 * - CPL_ERROR_ILLEGAL_INPUT: the number of dimensions is less than 2
 */
/*----------------------------------------------------------------------------*/
cpl_vector      *clipm_priv_vector_get_normal(
                                            const cpl_vector    **dv)
{
    cpl_vector  *resultv = NULL;
    cpl_matrix  *altm = NULL;
    const double
                **dv_data = NULL;
    
    CLIPM_TRY
    {
        int     vectorsize,
                nvectors,
                i,
                alternate;
        double  *mdata = NULL;
        
        if (dv == NULL)
            CLIPM_TRY_EXIT_WITH_ERROR(      CPL_ERROR_NULL_INPUT);
        if ((*dv) == NULL)
            CLIPM_TRY_EXIT_WITH_ERROR(      CPL_ERROR_NULL_INPUT);
        
        vectorsize = cpl_vector_get_size(*dv);
        nvectors = vectorsize - 1;
        
        CLIPM_TRY_CHECK_AUTOMSG(
            vectorsize >= 2,
            CPL_ERROR_ILLEGAL_INPUT);
        
        dv_data = cpl_malloc(nvectors * sizeof(*dv_data));
        CLIPM_TRY_CHECK_ERROR_STATE();
        
        for (i = 0; i < nvectors; i++)
        {
            if (dv[i] == NULL)
                CLIPM_TRY_EXIT_WITH_ERROR(  CPL_ERROR_NULL_INPUT);
            if (cpl_vector_get_size(dv[i]) != vectorsize)
                CLIPM_TRY_EXIT_WITH_ERROR_MSG(
                                            CPL_ERROR_INCOMPATIBLE_INPUT,
                                            "vectors",
                                            CLIPM_MSG_ERR_DIFFSIZES);
            dv_data[i] = cpl_vector_get_data_const(dv[i]);
        }
        
        altm = cpl_matrix_new(nvectors, nvectors);
        CLIPM_TRY_CHECK_ERROR_STATE();
        
        resultv = cpl_vector_new(vectorsize);
        CLIPM_TRY_CHECK_ERROR_STATE();
        
        alternate = 1;
        for (i = 0; i < vectorsize; i++)
        {
            /* copy vectors into matrix, except i-th row */
            int     row, col;
            double  det;

            mdata = cpl_matrix_get_data(altm);
            CLIPM_TRY_CHECK_ERROR_STATE();
            
            for (row = 0; row < i; row++)
            {
                for (col = 0; col < nvectors; col++)
                    *(mdata++) = dv_data[col][row];
            }
            row++;
            for ( ; row < vectorsize; row++)
            {
                for (col = 0; col < nvectors; col++)
                    *(mdata++) = dv_data[col][row];
            }
            
            det = cpl_matrix_get_determinant(altm);
            CLIPM_TRY_CHECK_ERROR_STATE();
            
            cpl_vector_set(resultv, i, alternate * det);
            CLIPM_TRY_CHECK_ERROR_STATE();
            
            alternate = -alternate;
        }
    }
    CLIPM_CATCH
    {
        clipm_priv_vector_null(&resultv);
        
        if (cpl_error_get_code() == CPL_ERROR_UNSPECIFIED)
            CLIPM_ERROR_SET_MSG(            CPL_ERROR_SINGULAR_MATRIX,
                                            "",
                                            "this error was mathematically "
                                            "not thought to be possible");
    }
    
    cpl_matrix_delete(altm);
    if (dv_data != NULL)
        cpl_free(dv_data);
    
    return resultv;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief   Delete a CPL vector object and set the pointer to NULL.
 * @param   v   Pointer to vector pointer
 * @return  Nothing
 * 
 * The following code is executed:
 * @code
    if (v != NULL)
    {
        cpl_vector_delete(*v);  // checks for NULL pointer
        *v = NULL;
    }
 * @endcode
 * 
 * @par Error Handling:
 * No error can occur here.
 */
/*----------------------------------------------------------------------------*/
void            clipm_priv_vector_null(     cpl_vector          **v)
{
    if (v != NULL)
    {
        cpl_vector_delete(*v);  /* checks for NULL pointer */
        *v = NULL;
    }
}

/**@}*/
