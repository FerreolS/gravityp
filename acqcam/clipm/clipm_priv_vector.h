
/*********************************************************************
 * E.S.O. - VLT project
 *
 * "@(#) $Id: clipm_priv_vector.h 176658 2008-11-26 18:02:24Z hlorch $"
 *
 * PRIVATE functions for vector handling
 *
 * who       when        what
 * --------  ----------  ----------------------------------------------
 * hlorch    2007-04-27  created
 */
 
#ifndef CLIPM_PRIV_VECTOR_H
#define CLIPM_PRIV_VECTOR_H

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

double          clipm_priv_vector_get_min(  const cpl_vector    *v,
                                            int                 *index);

double          clipm_priv_vector_get_max(  const cpl_vector    *v,
                                            int                 *index);

cpl_vector      *clipm_priv_vector_expand(  const cpl_vector    *v,
                                            int                 prepend_nr,
                                            int                 append_nr);

cpl_vector      *clipm_priv_vector_integrate(
                                            const cpl_vector    *v);

cpl_vector      *clipm_priv_vector_differentiate(
                                            const cpl_vector    *v);

cpl_vector      *clipm_priv_vector_get_normal(
                                            const cpl_vector    **dv);

void            clipm_priv_vector_null(     cpl_vector          **v);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}   /* extern "C" */
#endif

#endif /* CLIPM_PRIV_VECTOR_H */
