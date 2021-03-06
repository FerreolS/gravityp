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
 * History :
 * 12/11/2018  Add one gravi_Data* to  gravi_metrology_reduce
 */

#ifndef METREDLIB_H
#define METREDLIB_H

/*-----------------------------------------------------------------------------
                                   Includes
 -----------------------------------------------------------------------------*/

#include <stdio.h>
#include <cpl.h>
#include "gravi_data.h"

/*-----------------------------------------------------------------------------
                               Public prototypes
 -----------------------------------------------------------------------------*/

cpl_error_code gravi_metrology_reduce (gravi_data *,
                                       gravi_data *,
				       gravi_data *,
                                       gravi_data * met_pos,
                                       const cpl_parameterlist * parlist);

cpl_table * gravi_metrology_create (cpl_table * metrology_table,
                                    cpl_propertylist * header);

cpl_error_code gravi_metrology_drs (cpl_table * metrology_table,
                                    cpl_table * vismet_table,
                                    cpl_propertylist * header,
                                       const cpl_parameterlist * parlist);

cpl_error_code gravi_metrology_tac (cpl_table * metrology_table,
                                    cpl_table * vismet_table,
                                    cpl_propertylist * header);

cpl_table * gravi_metrology_compute_p2vm (cpl_table * metrology_table, double wave_met);

#endif 	/* !METREDLIB_H */
