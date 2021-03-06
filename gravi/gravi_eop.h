/*
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

#ifndef GRAVI_EOP_H
#define GRAVI_EOP_H

/*-----------------------------------------------------------------------------
                                   Includes
 -----------------------------------------------------------------------------*/

#include <cpl.h>
#include "gravi_calib.h"

/*-----------------------------------------------------------------------------
                               Public prototypes
 -----------------------------------------------------------------------------*/

cpl_error_code gravi_eop_pointing_uv (cpl_table * input_table,
                                      cpl_propertylist * header,
                                      cpl_table * eop_table,
                                      cpl_propertylist * eop_header,
							          int save_pointing,
							          cpl_table * array_table);
    
cpl_error_code gravi_compute_pointing_uv (gravi_data * p2vmred_data, 
                                          gravi_data * eop_data);

char * gravity_eop_download_finals2000A (const char * eop_host, 
                                         const char * eop_urlpath, 
                                         int * data_lenght);

cpl_table * gravity_eop_data_totable (const char * eop_data, int data_length);

#endif /* GRAVI_EOP_H */
