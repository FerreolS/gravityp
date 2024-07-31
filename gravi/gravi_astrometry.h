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

#ifndef GRAVI_ASTROMETRY_H
#define GRAVI_ASTROMETRY_H

#include <cpl.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include "gravi_data.h"

typedef struct _astrometry_data_ astro_data;

astro_data *gravi_astrometry_load(gravi_data *data);
cpl_error_code gravi_astrometry_dump(astro_data * self, FILE *handle);
void gravi_astrometry_delete(astro_data *self);

double gravi_astrometry_get_mean_ftflux(astro_data *self);
cpl_error_code gravi_astrometry_filter_ftflux(astro_data *self, double threshold);
cpl_error_code gravi_astrometry_normalise_to_ft(astro_data *self);
cpl_error_code gravi_astrometry_create_phase_reference(astro_data *self, astro_data **phase_refs, cpl_size nphase, astro_data **swaps, cpl_size nswap, cpl_parameterlist *parlist);
cpl_table *gravi_astrometry_get_phase_reference(astro_data *self);

cpl_error_code gravi_astrometry_reduce_swaps(astro_data **swap_data, cpl_size nswap, cpl_parameterlist *parlist);

#endif