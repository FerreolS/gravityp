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
 * ekw 07/12/2018 add wave_param to gravi_compute_argon_pos
 */
#ifndef GRAVI_DISP_H_
#define GRAVI_DISP_H_

/*-----------------------------------------------------------------------------
                                   Includes
 -----------------------------------------------------------------------------*/

#include <cpl.h>
#include "gravi_data.h"

/*-----------------------------------------------------------------------------
                                   Defines
 ----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------
                                 Public prototypes
 -----------------------------------------------------------------------------*/

gravi_data * gravi_compute_disp (gravi_data * vis_data);

cpl_error_code gravi_disp_cleanup (gravi_data * vis_data);

cpl_error_code gravi_compute_argon_pos (gravi_data * preproc_data, gravi_data *wave_param);

#endif /* GRAVI_DISP_H_ */
