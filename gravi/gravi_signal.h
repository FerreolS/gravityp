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

#ifndef GRAVI_SIGNAL_H_
#define GRAVI_SIGNAL_H_

/*-----------------------------------------------------------------------------
                                   Includes
 -----------------------------------------------------------------------------*/

#include <cpl.h>

/*-----------------------------------------------------------------------------
                              Public prototypes
 -----------------------------------------------------------------------------*/

cpl_error_code gravi_compute_snr (gravi_data * reduce_data,
                                  const cpl_parameterlist * parlist);

cpl_error_code gravi_compute_signals (gravi_data * reduce_data,
									gravi_data * disp_data,
									const cpl_parameterlist * parlist);

cpl_error_code gravi_compute_rejection (gravi_data * reduce_data,
										const cpl_parameterlist * parlist);

cpl_error_code gravi_flux_create_fddllin_sc (cpl_table * flux_SC,
                                             cpl_table * disp_table);

cpl_error_code gravi_vis_create_met_ft (cpl_table * vis_FT,
                                        cpl_table * vis_MET);

cpl_error_code gravi_vis_create_opdsc_ft (cpl_table * vis_FT,
                                          cpl_table * vis_SC,
                                          double dit_sc);

#endif /* GRAVI_SIGNAL_H_ */
