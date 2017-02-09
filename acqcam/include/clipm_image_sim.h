
/*********************************************************************
 * E.S.O. - VLT project
 *
 * "@(#) $Id: clipm_image_sim.h 177213 2008-12-05 15:33:51Z hlorch $"
 *
 * Functions for image creation and simulation
 *
 * who       when        what
 * --------  ----------  ----------------------------------------------
 * hlorch    2008-12-02  created
 */

#ifndef CLIPM_IMAGE_SIM_H
#define CLIPM_IMAGE_SIM_H

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

cpl_error_code  clipm_image_sim_noise_gaussian(
                                            cpl_image       *img_modified,
                                            const int       window_xxyy[4],
                                            double          sigma);

cpl_error_code  clipm_image_sim_circle(     cpl_image       *img_modified,
                                            double          centre_x,
                                            double          centre_y,
                                            double          radius,
                                            double          brightness);

cpl_error_code  clipm_image_sim_rectangle(  cpl_image       *img_modified,
                                            double          centre_x,
                                            double          centre_y,
                                            double          angle,
                                            double          length,
                                            double          width,
                                            double          brightness);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}   /* extern "C" */
#endif

#endif  /* CLIPM_IMAGE_SIM_H */
