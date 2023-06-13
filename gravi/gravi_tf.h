/* $Id: gravi_tf.h,v 1.12 2014/11/12 06:10:40 nazouaoui Exp $
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

#ifndef GRAVI_TF_H_
#define GRAVI_TF_H_

/*-----------------------------------------------------------------------------
                                   Includes
 -----------------------------------------------------------------------------*/

#include <cpl.h>

/*-----------------------------------------------------------------------------
                              Public prototypes
 -----------------------------------------------------------------------------*/

gravi_data * gravi_compute_tf (gravi_data *, gravi_data *);

gravi_data * gravi_calibrate_vis (gravi_data *, gravi_data ** , int,  gravi_data *, gravi_data *, const cpl_parameterlist *);

gravi_data * gravi_compute_zp (gravi_data ** , int);

cpl_error_code gravi_compute_tf_qc (gravi_data * oi_vis, gravi_data * diamcat_data);

cpl_error_code gravi_apply_tf_amp (gravi_data * science,
								   gravi_data * science_tf,
								   gravi_data ** used_tf_data,
								   int num_tf_data,
								   const char * extName,
								   const char * insName,
								   const char * ampName,
								   const char * ampErrName,
								   int nbase, double delta_t);

cpl_error_code gravi_apply_tf_phi( gravi_data * science,
 								   gravi_data * science_tf,
								   gravi_data ** used_tf_data,
								   int num_tf_data,
								   const char * extName,
								   const char * insName,
								   const char * phiName,
								   const char * phiErrName,
								   int nbase, double delta_t);

#endif /* GRAVI_TF_H_ */
