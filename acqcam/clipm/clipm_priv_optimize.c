/*********************************************************************
 * E.S.O. - VLT project
 *
 * "@(#) $Id: clipm_priv_optimize.c 176658 2008-11-26 18:02:24Z hlorch $"
 *
 * Private functions for optimization
 *
 * who       when        what
 * --------  ----------  ----------------------------------------------
 * hlorch    2007-04-30  created
 */

/**
 * @internal
 * @defgroup clipm_priv_optimize Optimization Procedures
 * @ingroup internal_docs
 *
 * This module provides private optimization procedures.
 *
 * @par Synopsis:
 * @code
#include "clipm_priv_optimize.h"
 * @endcode
 */
/**@{*/

/*-----------------------------------------------------------------------------
    Includes
 -----------------------------------------------------------------------------*/

#include "clipm_priv_optimize.h"

#include "clipm_math.h"
#include "clipm_compatibility_replacements.h"
#include "clipm_priv_error.h"
#include "clipm_priv_matrix.h"
#include "clipm_priv_vector.h"

/*-----------------------------------------------------------------------------
    Implementation
 -----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
 * @brief   Apply a Gaussian fitting until it is successful.
 * @param   xvalues         Positions to fit
 * @param   xsigmas         Uncertainties of positions to fit (NOT supported)
 * @param   yvalues         Values to fit
 * @param   ysigmas         Uncertainties of values to fit
 * @param   fit_opts        Fitting options: bitwise combination of
 *                          CPL_FIT_CENTROID,
 *                          CPL_FIT_STDEV,
 *                          CPL_FIT_AREA,
 *                          CPL_FIT_OFFSET or
 *                          CPL_FIT_ALL (for all above)
 * @param   max_iterations  Maximum number of iterations
 * @param   x0              (output) Center of best fit gaussian. If
 *                          CPL_FIT_CENTROID is not set, this is also an input
 *                          parameter.
 * @param   x0_uncertainty  (output) Optional uncertainty of x0.
 * @param   sigma           (output) Width of best fit gaussian. A positive
 *                          number on success. If CPL_FIT_STDEV is not set, this
 *                          is also an input parameter.
 * @param   sigma_uncertainty   (output) Optional uncertainty of sigma.
 * @param   fwhm            (output) FWHM of best fit gaussian. A positive
 *                          number on success.
 * @param   fwhm_uncertainty    (output) Optional uncertainty of fwhm.
 * @param   area            (output) Area of gaussian. A positive number on
 *                          success. If CPL_FIT_AREA is not set, this is also an
 *                          input parameter.
 * @param   offset          (output) Fitted background level. If CPL_FIT_OFFSET
 *                          is not set, this is also an input parameter.
 * @param   mse             (output) If non-NULL, the mean squared error of the
 *                          best fit is returned.
 * @param   red_chisq       (output) If non-NULL, the reduced chi square of the
 *                          best fit is returned. This requires the @a ysigmas
 *                          to be specified.
 * @param   covariance      (output) If non-NULL, the formal covariance matrix
 *                          of the best fit is returned. This requires
 *                          @a ysigmas to be specified. The order of fit
 *                          parameters in the covariance matrix is defined as
 *                          (x0, sigma, area, offset), for example the (3,3)
 *                          element of the matrix (counting from zero) is the
 *                          variance of the fitted offset. The matrix must be
 *                          deallocated by calling cpl_matrix_delete() . On
 *                          error, NULL is returned.
 * @param   lower_ndx       Output lower index of iteration window, can be NULL
 * @param   upper_ndx       Output upper index of iteration window, can be NULL
 * @return  CPL error code
 * 
 * @par Principle:
 * A gaussian fitting is performed. If @a max_iterations > 0, a narrow window
 * around the found peak is defined, and the fitting repeated in this data
 * window. This is intended to exclude noise influences far away from the peak.
 * 
 * @par Constraints:
 * The vector @a xvalues must be sorted. If it is not, the iteration might fail.
 * 
 * @par Iteration Details:
 * - Iteration window definition:
 *   - If a gaussian fitting was successful, the next window is defined to cover
 *     the +/- 5 sigma region around the peak (using the determined sigma).
 *   - If a gaussian fitting was not successful, the +/- 1 sigma region is used,
 *     using the centre and sigma guess values from cpl_vector_fit_gaussian().
 *     This guessed sigma is usually much wider than the real one, that's why.
 *   - The iteration window is always at least 5 data points wide.
 * - The iteration is continued until either the window is stable or until
 *   @a max_iteration steps are performed.
 * 
 * @par Output Parameter Details:
 * The output parameters and returned error codes are defined to be the same
 * as from function @a cpl_vector_fit_gaussian(), with the following
 * extensions:
 * -# @a cpl_fit_gaussian() is called, and on failure a fitting window is
 *    defined to constrain the data and try another fitting, as long as
 *    the maximum number of trials @a max_iterations is not reached.
 * -# If @a max_iterations allows, on success a fitting is tried once more
 *    in a window defined as the +/- 5 sigma region, to improve the result.
 * -# @a xvalues does not need to be provided. If it is NULL, the indizes are
 *    taken as positions, starting with 0.
 * 
 * @par Error Handling:
 * The following error cases might happen, and the corresponding error code set:
 * - CPL_ERROR_NULL_INPUT: @a yvalues or @a x0 or @a sigma or @a area or
 *   @a offset is NULL
 * - CPL_ERROR_UNSUPPORTED_MODE:
 *   - @a xsigmas != NULL
 *   - the specified @a fit_opts is not a bitwise
 *     combination of the allowed values (e.g. 0 or 1)
 * - CPL_ERROR_ILLEGAL_INPUT:
 *   - @a max_iterations < 0
 *   - the length of @a yvalues < 5
 *   - any input noise values, @a sigma or @a area (if input) is non-positive
 * - CPL_ERROR_INCOMPATIBLE_INPUT: any input vectors have different sizes
 * - CPL_ERROR_CONTINUE: the optimization could not converge
 * - CPL_ERROR_SINGULAR_MATRIX: the @a covariance matrix could not be calculated
 * 
 * @todo
 * - Enable cpl_vector_fit_gaussian() error history again when fixed
 *   [DFS05410].
 * - Return area and offset uncertainties.
 * - Fix iteration window for non-equidistant xvalues.
 */
/*----------------------------------------------------------------------------*/
cpl_error_code  clipm_priv_optimize_gaussian(
                                            const cpl_vector    *xvalues,
                                            const cpl_vector    *xsigmas,
                                            const cpl_vector    *yvalues,
                                            const cpl_vector    *ysigmas,
                                            cpl_fit_mode        fit_opts,
                                            int                 max_iterations,
                                            double      *x0,
                                            double      *x0_uncertainty,
                                            double      *sigma,
                                            double      *sigma_uncertainty,
                                            double      *fwhm,
                                            double      *fwhm_uncertainty,
                                            double      *area,
                                            double      *offset,
                                            double      *mse,
                                            double      *red_chisq,
                                            cpl_matrix  **covariance,
                                            int         *lower_ndx,
                                            int         *upper_ndx)
{
    cpl_vector  *used_xvalues = NULL,
                *used_ysigmas = NULL;
    cpl_matrix  *intern_cov_mat = NULL,
                **used_cov_mat = &intern_cov_mat;
    double      intern_red_chisq;
    
    CLIPM_TRY
    {
        cpl_error_code  local_error = CPL_ERROR_CONTINUE;
        int     iteration = 0,
                iterated_on_success = 0,
                vectorsize,
                n,
                imin,
                imax;
        const double
                *xval_data = NULL,
                *ysig_data = NULL,
                *yval_data = NULL,
                *xsig_data = NULL;
        
        /* initialize */
        clipm_priv_matrix_null(covariance);
        
        /* check input */
        CLIPM_TRY_CHECK_AUTOMSG(yvalues != NULL,    CPL_ERROR_NULL_INPUT);
        CLIPM_TRY_CHECK_AUTOMSG(x0 != NULL,         CPL_ERROR_NULL_INPUT);
        CLIPM_TRY_CHECK_AUTOMSG(sigma != NULL,      CPL_ERROR_NULL_INPUT);
        CLIPM_TRY_CHECK_AUTOMSG(area != NULL,       CPL_ERROR_NULL_INPUT);
        CLIPM_TRY_CHECK_AUTOMSG(offset != NULL,     CPL_ERROR_NULL_INPUT);
        
        /* provide covariance only if needed */
        if (covariance != NULL)
            used_cov_mat = covariance;
        else if (   x0_uncertainty == NULL &&
                    sigma_uncertainty == NULL &&
                    fwhm_uncertainty == NULL)
            used_cov_mat = NULL;
        
        /* xsigmas is still unsupported. nevertheless make the code below
         * independent from this check */
        CLIPM_TRY_CHECK(                    xsigmas == NULL,
                                            CPL_ERROR_UNSUPPORTED_MODE,
                                            "", "x sigmas not implemented");
        
        /* get common vectorsize */
        vectorsize = cpl_vector_get_size(yvalues);
        
        /* parameter checks */
        CLIPM_TRY_CHECK_AUTOMSG(            max_iterations >= 0,
                                            CPL_ERROR_ILLEGAL_INPUT);
        CLIPM_TRY_CHECK_AUTOMSG(            vectorsize >= 5,
                                            CPL_ERROR_ILLEGAL_INPUT);
        
        /* we need xvalues, so if not provided, create equidistant values */
        if (xvalues == NULL)
        {
            used_xvalues = cpl_vector_new(vectorsize);
            for (n = 0; n < vectorsize; n++)
                cpl_vector_set(used_xvalues, n, n);
        }
        else
            used_xvalues = (cpl_vector*)xvalues;
        
        /* we might need ysigmas */
        if (ysigmas == NULL)
        {
            /* if covariances are needed, ysigmas are also required */
            if (used_cov_mat != NULL)
            {
                used_ysigmas = cpl_vector_new(vectorsize);
                cpl_vector_fill(used_ysigmas, 1.0);
            }
        }
        else
            used_ysigmas = (cpl_vector*)ysigmas;

        if (    cpl_vector_get_size(used_xvalues) != vectorsize ||
                cpl_vector_get_size(used_ysigmas) != vectorsize)
            CLIPM_TRY_EXIT_WITH_ERROR_MSG(  CPL_ERROR_INCOMPATIBLE_INPUT,
                                            "input vectors",
                                            CLIPM_MSG_ERR_DIFFSIZES);
        
        xval_data = cpl_vector_get_data(used_xvalues);
        yval_data = cpl_vector_get_data_const(yvalues);
        if (xsigmas != NULL)
            xsig_data = cpl_vector_get_data_const(xsigmas);
        if (used_ysigmas != NULL)
            ysig_data = cpl_vector_get_data_const(used_ysigmas);

        /* start the iteration of the gaussian fitting */
        imin = 0;
        imax = vectorsize-1;
        while   (   iteration < max_iterations+1 &&
                    (   local_error == CPL_ERROR_CONTINUE ||
                        (   local_error == CPL_ERROR_NONE &&
                            iterated_on_success == 0
                )   )   )
        {
            cpl_vector  *xval_iter = NULL,
                        *ysig_iter = NULL;
            const cpl_vector
                        *xsig_iter = NULL,
                        *yval_iter = NULL;
            int         old_imin = -1,
                        old_imax = -1;
            
            /* do one refinement step after fitting was successful */
            if (local_error == CPL_ERROR_NONE)
                iterated_on_success++;
            
            /* adjust iteration window */
            if (iteration > 0)  /* if not the first try */
            {
                double  wdwhw;    /* window half width in sigmas */
                
                /* if an error happened, then the sigma was determined
                 * using some median computation, which normally
                 * results in high values, so restrict to 1-sigma
                 * region, otherwise to 5-sigma region */
                wdwhw = (local_error == CPL_ERROR_NONE) ? 5.0 : 1.0;
                
                if (xvalues != NULL)
                {
                    imin = cpl_vector_find(xvalues, *x0 - wdwhw*(*sigma));
                    imax = cpl_vector_find(xvalues, *x0 + wdwhw*(*sigma));
                }
                else
                {
                    imin = clipm_math_round_d2i((*x0) - wdwhw*(*sigma));
                    imax = clipm_math_round_d2i((*x0) + wdwhw*(*sigma));
                }
                imin = imin < 0 ? 0 : imin;
                imax = imax < 0 ? 0 : imax;
                imin = imin < vectorsize ? imin : vectorsize - 1;
                imax = imax < vectorsize ? imax : vectorsize - 1;
                
                /* at least 5 data points to fit */
                if (imax - imin < 4)
                {
                    imin -= (5 - (imax - imin + 1))/2;
                    
                    imin = imin < 0 ? 0 : imin;
                    imin = imin > vectorsize - 5 ? vectorsize - 5 : imin;

                    imax = imin + 4;
                }
                
                if (lower_ndx != NULL)
                    *lower_ndx = imin;
                if (upper_ndx != NULL)
                    *upper_ndx = imax;
            }
            
            /* if window is the same as before, then exit */
            if (    imin == old_imin &&
                    imax == old_imax)
                break;
            old_imin = imin;
            old_imax = imax;
            
            /* delete a previously calculated covariance matrix */
            clipm_priv_matrix_null(used_cov_mat);
            
            xval_iter = cpl_vector_wrap(imax - imin + 1,
                                            (double*)xval_data + imin);
            yval_iter = cpl_vector_wrap(imax - imin + 1,
                                            (double*)yval_data + imin);
            if (xsig_data != NULL)
                xsig_iter = cpl_vector_wrap(imax - imin + 1,
                                            (double*)xsig_data + imin);
            if (ysig_data != NULL)
                ysig_iter = cpl_vector_wrap(imax - imin + 1,
                                            (double*)ysig_data + imin);

            CLIPM_ERROR_RECOVER_TRYSTATE();
            
            local_error = cpl_vector_fit_gaussian(
                                            xval_iter, xsig_iter,
                                            yval_iter, ysig_iter,
                                            fit_opts,
                                            x0,
                                            sigma,
                                            area,
                                            offset,
                                            mse,
                                            &intern_red_chisq,
                                            used_cov_mat);

            /* calculate other values */
            if (fwhm != NULL)
                *fwhm = (*sigma) * CPL_MATH_FWHM_SIG;

            cpl_vector_unwrap(xval_iter);
            cpl_vector_unwrap((cpl_vector*)yval_iter);
            cpl_vector_unwrap((cpl_vector*)xsig_iter);
            cpl_vector_unwrap(ysig_iter);
            
            iteration++;
        }
        if (CLIPM_ERROR_GET_NEW_SINCE_TRY() == CPL_ERROR_INVALID_TYPE)
            CLIPM_TRY_EXIT_WITH_ERROR_MSG(  CPL_ERROR_UNSUPPORTED_MODE,
                                            "fit_opts", "");
        /* cpl_vector_fit_gaussian() floods the error history when failing */
        if (CLIPM_ERROR_GET_NEW_SINCE_TRY() == CPL_ERROR_CONTINUE)
        {
            CLIPM_ERROR_RECOVER_TRYSTATE();
            CLIPM_ERROR_SET_MSG(
                CPL_ERROR_CONTINUE,
                "", "gaussian fitting failed");
        }
        /*CLIPM_ERROR_SET_MSG_IF_CODE(
            CPL_ERROR_CONTINUE,
            "", "gaussian fitting failed");*/
        CLIPM_ERROR_SET_MSG_IF_CODE(
            CPL_ERROR_SINGULAR_MATRIX,
            "", "the covariances could not be computed");
        CLIPM_ERROR_SET_MSG_IF_CODE(
            CPL_ERROR_ILLEGAL_INPUT,
            "", "input noise values, area (if input) and sigma (if input) "
                "must be positive");
        CLIPM_ERROR_SET_MSG_IF_CODE(
            CPL_ERROR_ILLEGAL_OUTPUT,
            "", "memory allocation failed");
        CLIPM_TRY_CHECK_ERROR_STATE();
        
        if (local_error == CPL_ERROR_NONE &&
            used_cov_mat != NULL && *used_cov_mat != NULL)
        {
            if (x0_uncertainty != NULL)
                *x0_uncertainty =   sqrt(   intern_red_chisq *
                                            cpl_matrix_get(
                                                *used_cov_mat, 0, 0));
            if (sigma_uncertainty != NULL)
                *sigma_uncertainty = sqrt(  intern_red_chisq *
                                            cpl_matrix_get(
                                                *used_cov_mat, 1, 1));
            if (fwhm_uncertainty != NULL)
                *fwhm_uncertainty = CPL_MATH_FWHM_SIG *
                                    sqrt(   intern_red_chisq *
                                            cpl_matrix_get(
                                                *used_cov_mat, 1, 1));
        }
    }
    CLIPM_CATCH
    {
        clipm_priv_matrix_null(used_cov_mat);
        
        if (x0_uncertainty != NULL)
            *x0_uncertainty = -1;
        if (sigma_uncertainty != NULL)
            *sigma_uncertainty = -1;
        if (fwhm_uncertainty != NULL)
            *fwhm_uncertainty = -1;
    }
    
    if (used_xvalues != xvalues)
        cpl_vector_delete(used_xvalues);
    if (used_ysigmas != ysigmas)
        cpl_vector_delete(used_ysigmas);
    if (red_chisq != NULL)
        *red_chisq = intern_red_chisq;
    if (used_cov_mat != covariance)
        cpl_matrix_delete(*used_cov_mat);
    
    return CLIPM_ERROR_GET_NEW_SINCE_TRY();
}

/*----------------------------------------------------------------------------*/
/**
 * @brief   Optimize a set of parameters using the Downhill Simplex method by
 *          Nelder and Mead (also known as the "Amoeba" method).
 * @param   fit_params      Vector containing the parameters to be optimized,
 *                          which must be set to an initial guess, and will
 *                          contain the result at the end
 * @param   start_delta     Vector containing a set of parameter step sizes
 *                          for the start of the optimization
 * @param   precision       Optional vector containing respective precision
 *                          limits (> 0) which must be met for convergence,
 *                          can be NULL
 * @param   evalfunc        Pointer to function which evaluates the parameters
 *                          and returns a value which must be minimal in the
 *                          optimum (for example an error value, or the
 *                          inverse or negative of a matching), taking as
 *                          first argument the parameter vector, and as second
 *                          argument an arbitrary pointer (to provide any
 *                          other data)
 * @param   other_func_par  Arbitrary pointer used as the input for the
 *                          second argument of @a evalfunc, for example a
 *                          pointer to a dataset or a reference to a struct
 *                          of pointers and further parameters
 * @param   max_iterations  The maximum number of iterations
 * @param   optimum         Optional output return value of @a evalfunc with
 *                          optimized parameters, can be NULL
 * @param   done_iterations Optional output of number of performed iterations,
 *                          can be NULL
 * @return  CPL error code
 * 
 * @par Indication:
 * - This function optimizes a set of paramaters, which can for example be used
 *   in a fitting process. The feature of the Downhill Simplex method is that
 *   it works without an analytical derivative of the evaluation function. If
 *   a derivative is available, the usage of @a cpl_fit_lvmq()
 *   (Levenberg-Marquard fitting) should be preferred.
 * 
 * @par Usage:
 * Two modes of operation are provided:
 * -# If @a precision is given, then @a fit_params is optimized until the
 *    variations of the parameters stay below their respective @a precision
 *    entry. This is defined as convergence. If @a max_iterations is reached
 *    before convergence, @a fit_params will contain the latest parameter
 *    values, and the error code CPL_ERROR_CONTINUE is set and returned.
 * -# If @a precision is NULL, then @a max_iterations iterations are performed,
 *    and @a fit_params contains the result at the end.
 * 
 * @par Evaluation Function:
 * - If @a evalfunc sets a CPL error, this is interpreted as a real error, and
 *   the optimization is stopped. In this case, this error code is set and
 *   returned. @a fit_params contains the latest parameters.
 * - If there are hard limits for the parameters, these must be implemented in
 *   the way that @a evalfunc returns a very bad (high) value, but does not
 *   set an error.
 * 
 * @note
 * - Although this function works also with only one dimension (i.e. a
 *   parameter vector of size 1), this is not recommended since there are
 *   better optimization methods for the univariate case.
 * - The number of optimization iterations does not equal the number of calls
 *   of @a evalfunc. Some iteration steps require several evaluations, i.e.
 *   the number of @a evalfunc calls is higher than the number of iterations.
 * - @a done_iterations, if not NULL, returns the number of performed
 *   iterations. This is intended for instance for debugging purposes, or for
 *   the optimization of methods using this function.
 * 
 * @a Limitations:
 * - It is possible that a local minimum is found. There are different
 *   approaches to solve this issue, one possibility is for example to restart
 *   the optimization with a reasonable high initial step size
 *   (@a start_delta), or to apply the optimization to a cluster of
 *   parameter guesses. But it should be tried out whether this is necessary
 *   at all, before taking the effort. This depends strongly on the data,
 *   and a general answer to this topic is hardly possible.
 * 
 * @par Error Handling:
 * If an error occurs, but some iteration could be made, then the latest
 * parameter values are stored in @a fit_params.
 * \n\n
 * The following error codes can be set and returned:
 * - CPL_ERROR_NULL_INPUT: @a fit_params or @a start_delta or @a evalfunc
 *   is/are NULL
 * - CPL_ERROR_INCOMPATIBLE_INPUT: @a fit_params and @a start_delta or
 *   @a precision (if not NULL) have different sizes
 * - CPL_ERROR_ILLEGAL_INPUT: an entry of @a precision (if not NULL) is <= 0
 * - CPL_ERROR_CONTINUE: the optimization could not converge to @a precision
 *   (see above)
 * - Any error that is set by @a evalfunc
 * 
 * @see
 * Nelder, J. A. & Mead, R. A: A Simplex Method for Function Minimization.
 * Computer Journal, 1965, 7, 308-313
 */
/*----------------------------------------------------------------------------*/
cpl_error_code  clipm_priv_optimize_downhill_simplex(
                                            cpl_vector          *fit_params,
                                            const cpl_vector    *start_delta,
                                            const cpl_vector    *precision,
                                            double  (*evalfunc)
                                                    (cpl_vector *params,
                                                     void       *other),
                                            void    *other_func_par,
                                            int     max_iterations,
                                            double  *optimum,
                                            int     *done_iterations)
{
    int         nrparams = 0,
                nrpoints = 0,
                bestpos = -1, /* define: bestpos < 0, unless started */
                iteration = 0;
    double      **pointdata = NULL;
    cpl_vector  *values = NULL,
                *newpoint = NULL,
                *newerpoint = NULL,
                *mean_of_good = NULL;
    cpl_vector  **points = NULL;
    
    CLIPM_TRY
    {
        int     ipar,
                ipo,
                below_precision = 0,
                convergence_indicator = 0;
        double  *valuedata = NULL;
        const double
                *precisiondata = NULL;
        
        /* init output */
        if (optimum != NULL)
            *optimum = 0;
        
        /* check input */
        nrparams = cpl_vector_get_size(fit_params);
        {
            int n2 = cpl_vector_get_size(start_delta);
            CLIPM_TRY_CHECK_ERROR_STATE();
            if (n2 != nrparams)
                CLIPM_TRY_EXIT_WITH_ERROR_MSG(
                                            CPL_ERROR_INCOMPATIBLE_INPUT,
                                            "fit_params, start_delta",
                                            CLIPM_MSG_ERR_DIFFSIZES);
        }
        
        if (evalfunc == NULL)
            CLIPM_TRY_EXIT_WITH_ERROR(  CPL_ERROR_NULL_INPUT);
        
        /* precision is optional */
        if (precision != NULL)
        {
            if (nrparams != cpl_vector_get_size(precision))
                CLIPM_TRY_EXIT_WITH_ERROR_MSG(
                                            CPL_ERROR_INCOMPATIBLE_INPUT,
                                            "fit_params, precision",
                                            CLIPM_MSG_ERR_DIFFSIZES);
            precisiondata = cpl_vector_get_data_const(precision);
            for (ipar = 0; ipar < nrparams; ipar++)
                if (precisiondata[ipar] <= 0)
                    CLIPM_TRY_EXIT_WITH_ERROR_MSG(
                                            CPL_ERROR_ILLEGAL_INPUT,
                                            "precision",
                                            "values must be > 0");
        }
        
        /* set up working variables, vectors... */
        nrpoints = nrparams + 1;
        
        values = cpl_vector_new(nrpoints);
        CLIPM_TRY_CHECK_ERROR_STATE();
        valuedata = cpl_vector_get_data(values);
        
        points = cpl_calloc(nrpoints, sizeof(*points));
        pointdata = cpl_calloc(nrpoints, sizeof(*pointdata));
        CLIPM_TRY_ASSERT(points != NULL && pointdata != NULL);

        mean_of_good = cpl_vector_new(nrparams);
        newpoint = cpl_vector_new(nrparams);
        newerpoint = cpl_vector_new(nrparams);
        CLIPM_TRY_ASSERT_ERROR_STATE();

        for (ipo = 0; ipo < nrpoints; ipo++)
        {
            CLIPM_TRY_EXIT_IFN(
                points[ipo] = cpl_vector_duplicate(fit_params));
            CLIPM_TRY_EXIT_IFN(
                pointdata[ipo] = cpl_vector_get_data(points[ipo]));
        }
        
        /* for each new point, vary ONE parameter at beginning */
        for (ipar = 0; ipar < nrparams; ipar++)
        {
            ipo = ipar + 1;
            pointdata[ipo][ipar] += cpl_vector_get(start_delta, ipar);
        }
        
        /* and get the function values at all points */
        for (ipo = 0; ipo < nrpoints; ipo++)
        {
            valuedata[ipo] = (*evalfunc)(points[ipo], other_func_par);
            CLIPM_TRY_CHECK_ERROR_STATE();
        }
        
        iteration = 0;
        while (iteration < max_iterations)
        {
            int         worstpos;
            double      minvalue,
                        maxvalue,
                        newvalue;
            
            /* get best/worst points.
             * this is the only time per iteration when the best is
             * determined, so if something fails later, the best parameter
             * set can still be returned (the algorithm does not modify the
             * best one) */
            minvalue = clipm_priv_vector_get_min(values, &bestpos);
            maxvalue = clipm_priv_vector_get_max(values, &worstpos);
            
            /* check whether we are converging */
            if (precisiondata != NULL)
            {
                below_precision = 1;
                for (ipar = 0; ipar < nrparams; ipar++)
                {
                    if (    fabs(           pointdata[worstpos][ipar] -
                                            pointdata[bestpos][ipar])
                            > precisiondata[ipar])
                    {
                        below_precision = 0;
                    }
                }
                if (below_precision)
                    convergence_indicator++;
                else
                    convergence_indicator = 0;
                
                /* decide on convergence: if we are below the precision for
                 * 2*nrparams times */
                if (convergence_indicator >= 2*nrparams)
                    break;
            }

            /* get the mean of the points which are not the worst one */
            ipar = 0;
            for (ipo = 0; ipo < nrpoints; ipo++)
            {
                if (ipo != worstpos)
                {
                    if (ipar == 0)
                        cpl_vector_copy(mean_of_good, points[ipo]);
                    else
                        cpl_vector_add(mean_of_good, points[ipo]);
                    ipar++;
                }
            }
            cpl_vector_multiply_scalar(mean_of_good, 1/(double)nrparams);
            
            /* mirror worst point at mean of good */
            cpl_vector_copy(newpoint, mean_of_good);
            cpl_vector_multiply_scalar(newpoint, 2);
            cpl_vector_subtract(newpoint, points[worstpos]);
            CLIPM_TRY_CHECK_ERROR_STATE();
            
            newvalue = (*evalfunc)(newpoint, other_func_par);
            CLIPM_TRY_CHECK_ERROR_STATE();
            
            /* Nelder & Mead:
             * - if the new point is the best, then try a bigger step
             *   in the same direction
             * - if the new point is just better, replace the worst by it
             * - if the new point is worse than before,
             *   first try to shrink the old worst point towards the mean of
             *   good, if this doesn't help, then shrink the whole simplex */
            if (newvalue < minvalue)
            {
                double  newervalue;

                cpl_vector_copy(newerpoint, newpoint);
                cpl_vector_multiply_scalar(newerpoint, 2);
                cpl_vector_subtract(newerpoint, mean_of_good);
                newervalue = (*evalfunc)(newerpoint, other_func_par);
                CLIPM_TRY_CHECK_ERROR_STATE();
                
                if (!(newervalue < newvalue))
                {
                    cpl_vector_copy(points[worstpos], newpoint);
                    cpl_vector_set(values, worstpos, newvalue);
                }
                else do
                {
                    cpl_vector_copy(points[worstpos], newerpoint);
                    cpl_vector_set(values, worstpos, newervalue);
                    newvalue = newervalue;
                    
                    /* extend the Nelder-Mead algorithm by this while-loop
                     * to further try values in the good direction */
                    cpl_vector_multiply_scalar(newerpoint, 2);
                    cpl_vector_subtract(newerpoint, mean_of_good);
                    newervalue = (*evalfunc)(newerpoint, other_func_par);
                    CLIPM_TRY_CHECK_ERROR_STATE();
                }
                while (newervalue < newvalue);
            }
            else if (newvalue < maxvalue)
            {
                cpl_vector_copy(points[worstpos], newpoint);
                cpl_vector_set(values, worstpos, newvalue);
            }
            else
            {
                /* get the point between mean of good and the worst */
                cpl_vector_copy(newpoint, points[worstpos]);
                cpl_vector_add(newpoint, mean_of_good);
                cpl_vector_multiply_scalar(newpoint, 0.5);
                newvalue = (*evalfunc)(newpoint, other_func_par);
                CLIPM_TRY_CHECK_ERROR_STATE();
                
                /* if now better, use this */
                if (newvalue < maxvalue)
                {
                    cpl_vector_copy(points[worstpos], newpoint);
                    cpl_vector_set(values, worstpos, newvalue);
                }
                else /* otherwise shrink the simplex around best point */
                {
                    for (ipo = 0; ipo < nrpoints; ipo++)
                    {
                        if (ipo != bestpos)
                        {
                            cpl_vector_add(points[ipo], points[bestpos]);
                            cpl_vector_multiply_scalar(points[ipo], 0.5);
                            newvalue = (*evalfunc)(points[ipo], other_func_par);
                            CLIPM_TRY_CHECK_ERROR_STATE();
                            cpl_vector_set(values, ipo, newvalue);
                        }
                    }
                }
            }
            
            iteration++;
        }
        /* we could get here without error, so get result index */
        clipm_priv_vector_get_min(values, &bestpos);
        
        /* decide on convergence error */
        if (    precisiondata != NULL &&
                iteration >= max_iterations)
            CLIPM_TRY_EXIT_WITH_ERROR_MSG(  CPL_ERROR_CONTINUE,
                                            "",
                                            "could not fall well below given "
                                            "precision in max_iterations");
        
    }
    CLIPM_CATCH
    {
    }
    
    /* save the latest result, and the number of iterations */
    if (points != NULL &&
        bestpos >= 0 &&     /* we did some iterations (convergence unclear) */
        bestpos < nrpoints)
    {
        cpl_vector_copy(fit_params, points[bestpos]);
        
        /* if we did some iterations, we also have a vector "values" */
        if (    optimum != NULL &&
                values != NULL)
            *optimum = cpl_vector_get(values, bestpos);
    }
    if (done_iterations != NULL)
        *done_iterations = iteration;
    
    /* and clean up */
    cpl_vector_delete(mean_of_good);
    cpl_vector_delete(newerpoint);
    cpl_vector_delete(newpoint);
    cpl_vector_delete(values);
    
    if (points != NULL)
    {
        int n;
        for (n = 0; n < nrpoints; n++)
            cpl_vector_delete(points[n]);
        cpl_free(points);
    }
    if (pointdata != NULL)
        cpl_free(pointdata);
    
    return CLIPM_ERROR_GET_NEW_SINCE_TRY();
}

/**@}*/
