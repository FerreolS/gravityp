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

#ifndef GRAVI_CALIB_H_
#define GRAVI_CALIB_H_

/*-----------------------------------------------------------------------------
                                   Includes
 -----------------------------------------------------------------------------*/

#include <cpl.h>
#include "gravi_data.h"

/*-----------------------------------------------------------------------------
                                   Defines
 ----------------------------------------------------------------------------*/

#define BADPIX_DARK		(1 << 0)
#define BADPIX_RMS		(1 << 1)
#define BADPIX_FLAT		(1 << 2)
#define BADPIX_PROFILE_LIMIT		0.05  //level on the profile to flag flat bad pix into the badpix

/*-----------------------------------------------------------------------------
                                 Public prototypes
 -----------------------------------------------------------------------------*/

gravi_data * gravi_compute_dark (gravi_data * );
gravi_data * gravi_average_dark (gravi_data ** data, cpl_size ndata);

gravi_data * gravi_compute_profile (gravi_data ** , gravi_data * , gravi_data * ,
									int , const cpl_parameterlist * );

cpl_propertylist * gravi_compute_gain (gravi_data ** flats_data,
                                       int nrawgain, gravi_data * dark_map);

gravi_data * gravi_compute_badpix (gravi_data * ,
                                   gravi_data ** flats_data,
                                   int nflat,
                                   const cpl_parameterlist * );

gravi_data * gravi_compute_biasmask (gravi_data * dark_map,
                                     gravi_data ** flats_data,
                                     int nflat,
                                     const cpl_parameterlist * params);

gravi_data * gravi_compute_baseline (gravi_data * );



#endif /* GRAVI_PREPROC_H_ */
