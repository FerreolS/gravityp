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

#ifndef GRAVI_WAVE_H_
#define GRAVI_WAVE_H_

/*-----------------------------------------------------------------------------
                                   Includes
 -----------------------------------------------------------------------------*/

#include <cpl.h>
#include "gravi_data.h"

/*-----------------------------------------------------------------------------
                                Public prototypes
 -----------------------------------------------------------------------------*/

cpl_table * gravi_compute_argon_wave (cpl_table * spectrum_table);

cpl_error_code gravi_wave_compute_opds (gravi_data * spectrum_data,
                                        cpl_table  * met_table,
                                        const cpl_parameterlist * parlist);

cpl_error_code  gravi_compute_wave (gravi_data * wave_map,
                                    gravi_data * spectrum_data,
                                    int type_data, const cpl_parameterlist * parlist,
                                    gravi_data * wave_param);

cpl_error_code gravi_wave_qc (gravi_data * wave_map,
                              gravi_data * profile_map);

cpl_error_code gravi_wave_correct_color (gravi_data * vis_data);

#endif /* GRAVI_WAVE_H_ */
