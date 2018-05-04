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

#ifndef GRAVI_P2VMRED_H_
#define GRAVI_P2VMRED_H_

/*-----------------------------------------------------------------------------
                                   Includes
 -----------------------------------------------------------------------------*/

#include <cpl.h>

/*-----------------------------------------------------------------------------
                               Public prototypes
 -----------------------------------------------------------------------------*/

gravi_data * gravi_compute_p2vmred (gravi_data *, gravi_data *, const char *, 
                                    const cpl_parameterlist *,
                                    enum gravi_detector_type det_type);

cpl_error_code gravi_compute_opdc_state (gravi_data * p2vmred_data);

cpl_error_code gravi_compute_tau0 (gravi_data * data);

cpl_error_code gravi_compute_qc_injection (gravi_data * data);

cpl_error_code gravi_compute_qc_ft_opd_estimator (gravi_data * data);

#endif /* GRAVI_P2VMRED_H_ */
