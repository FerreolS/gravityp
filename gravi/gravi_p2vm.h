/* $Id: gravi_utils.h,v 1.12 2011/05/31 06:10:40 nazouaoui Exp $
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
 * History
 * ekw  04/12/2018 add *gravi_data wave_param to gravi_create_p2vm
 */
#ifndef GRAVI_P2VM_H_
#define GRAVI_P2VM_H_

/*-----------------------------------------------------------------------------
                                   Includes
 -----------------------------------------------------------------------------*/

#include <cpl.h>
#include "gravi_data.h"

/*-----------------------------------------------------------------------------
                               Public prototypes
 -----------------------------------------------------------------------------*/

gravi_data * gravi_create_p2vm (gravi_data * wave_cal_data, gravi_data *wave_param);
cpl_table* gravi_create_p2vm_table (cpl_table * detector_table, int n_wave);
cpl_error_code gravi_compute_p2vm (gravi_data * , gravi_data * , int ** , int ** ,
                                   enum gravi_detector_type det_type);
cpl_error_code gravi_p2vm_normalisation (gravi_data * , int ** , int **  );
cpl_error_code gravi_p2vm_phase_correction (gravi_data * p2vm_map, gravi_data * p2vmreduced, int full_phase);
cpl_error_code gravi_p2vm_transmission (gravi_data * p2vm_map, gravi_data * p2vmreduced);

#endif /* GRAVI_P2VM_H_ */
