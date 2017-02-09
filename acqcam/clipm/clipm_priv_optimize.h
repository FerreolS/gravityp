
/*********************************************************************
 * E.S.O. - VLT project
 *
 * "@(#) $Id: clipm_priv_optimize.h 164853 2008-03-17 16:43:59Z hlorch $"
 *
 * PRIVATE functions for optimization
 *
 * who       when        what
 * --------  ----------  ----------------------------------------------
 * hlorch    2007-02-07  created
 */
 
#ifndef CLIPM_PRIV_OPTIMIZE_H
#define CLIPM_PRIV_OPTIMIZE_H

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
    Prototypes
 -----------------------------------------------------------------------------*/

cpl_error_code  clipm_priv_optimize_gaussian(
                                            const cpl_vector    *xvalues,
                                            const cpl_vector    *xsigmas,
                                            const cpl_vector    *yvalues,
                                            const cpl_vector    *ysigmas,
                                            cpl_fit_mode        fit_pars,
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
                                            int         *upper_ndx);

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
                                            int     *done_iterations);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}   /* extern "C" */
#endif

#endif /* CLIPM_PRIV_OPTIMIZE_H */
