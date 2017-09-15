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

#ifndef GRAVI_PREPROC_H_
#define GRAVI_PREPROC_H_

/*-----------------------------------------------------------------------------
                                   Includes
 -----------------------------------------------------------------------------*/

#include <cpl.h>
#include "gravi_data.h"
#include "gravi_pfits.h"

/*-----------------------------------------------------------------------------
                              Public prototypes
 -----------------------------------------------------------------------------*/

int gravi_pixel_is_good (cpl_image * bad_img, int x, int y);

cpl_error_code gravi_remove_badpixel_sc (cpl_imagelist * imglist_sc, cpl_image * bad_img);

gravi_data * gravi_extract_spectrum (gravi_data * raw_data,
									 gravi_data * profile_map,
									 gravi_data * dark_map,
									 gravi_data * bad_map,
									 gravi_data * sky_map,
                                     const cpl_parameterlist * parlist,
                                     enum gravi_detector_type det_type);

cpl_error_code gravi_align_spectrum (gravi_data * spectrum_data,
                                     gravi_data * wave_map,
                                     gravi_data * p2vm_map,
                                     enum gravi_detector_type det_type);

#endif /* GRAVI_CALIB_H_ */
