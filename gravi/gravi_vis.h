/* $Id: gravi_vis.h,v 1.12 2014/11/12 06:10:40 nazouaoui Exp $
 *
 * This file is part of the GRAVI Pipeline
 * Copyright (C) 2002,2003 European Southern Observatory
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/*
 * $Author: nazouaoui $
 * $Date: 2011/05/31 06:10:40 $
 * $Revision: 1.12 $
 * $Name:  $
 */

#ifndef GRAVI_VIS_H_
#define GRAVI_VIS_H_

/*-----------------------------------------------------------------------------
                                   Includes
 -----------------------------------------------------------------------------*/

#include <cpl.h>

/*-----------------------------------------------------------------------------
                                       Prototypes
 -----------------------------------------------------------------------------*/

gravi_data * gravi_compute_vis (gravi_data *p2vmred_data,
                                const cpl_parameterlist *,
                                cpl_size * current_frame);

cpl_error_code gravi_compute_vis_qc (gravi_data * vis_data);

cpl_error_code gravi_vis_mjd_to_time (gravi_data * vis_data);

cpl_error_code gravi_flat_flux (gravi_data * vis_calib, gravi_data * p2vm_data);
cpl_error_code gravi_normalize_sc_to_ft (gravi_data * visdata);

cpl_error_code gravi_average_vis (gravi_data * oi_merged);
cpl_error_code gravi_vis_resamp (gravi_data * oi_merged, cpl_size nsamp);
cpl_error_code gravi_vis_smooth (gravi_data * oi_data,
                                 cpl_size nsamp_vis, cpl_size nsamp_flx,
                                 cpl_size maxdeg);

cpl_error_code gravi_vis_flag_threshold (cpl_table * oi_table, const char * data, const char *flag, double value);
cpl_error_code gravi_vis_flag_relative_threshold (cpl_table * oi_table, const char * err, const char * data, const char *flag, double value);
cpl_error_code gravi_vis_flag_lower (cpl_table * oi_table, const char * data, const char *flag, double value);

cpl_error_code gravi_force_uncertainties (gravi_data * oi_data,
                                          const cpl_parameterlist * parlist);

cpl_error_code gravi_vis_copy_fluxdata (gravi_data * oi_data);

cpl_error_code gravi_vis_erase_obs (cpl_table * oi_table, cpl_array *flag_array, cpl_size ntel);
cpl_error_code gravi_vis_force_time (gravi_data * oi_data);

#endif /* GRAVI_VIS_H_ */
