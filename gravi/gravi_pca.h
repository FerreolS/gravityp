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

#ifndef GRAVI_PCA_H
#define GRAVI_PCA_H

#include <cpl.h>

typedef struct _gravi_pca_result_ gravi_pca_result;
typedef struct _gravi_pca_model_ gravi_pca_model;

gravi_pca_result *gravi_pca_create_result(const cpl_matrix *data, const cpl_matrix *mask) CPL_ATTR_ALLOC;
cpl_error_code gravi_pca_decomp_matrix_svd(gravi_pca_result *self);
cpl_error_code gravi_pca_set_component_signs(gravi_pca_result *self, const cpl_array *wave, int num_components);

gravi_pca_model *gravi_pca_load_model(const cpl_matrix *components) CPL_ATTR_ALLOC;
gravi_pca_model *gravi_pca_create_model(const gravi_pca_result **results, int num_results, int num_components) CPL_ATTR_ALLOC;

cpl_error_code gravi_pca_fit_components_polynomial(gravi_pca_model *self, const cpl_array *wave, int degree, int num_components);
cpl_error_code gravi_pca_fit_components_bspline(gravi_pca_model *self, const cpl_array *wave, int degree, int num_components);
cpl_error_code gravi_pca_refine_mean(gravi_pca_result *self, const cpl_matrix *residual);
cpl_error_code gravi_pca_fit_model(gravi_pca_result *self, const gravi_pca_model *model, cpl_boolean fit_mean_subtracted, cpl_boolean verbose);

cpl_vector *gravi_pca_get_component(const gravi_pca_result *self, int component) CPL_ATTR_ALLOC;
cpl_vector *gravi_pca_get_component_median(const gravi_pca_model *self, int component) CPL_ATTR_ALLOC;
cpl_vector *gravi_pca_get_component_fit(const gravi_pca_model *self, int component) CPL_ATTR_ALLOC;

cpl_matrix *gravi_pca_get_data_fit(const gravi_pca_result *self) CPL_ATTR_ALLOC;
cpl_matrix *gravi_pca_get_data_residual(const gravi_pca_result *self) CPL_ATTR_ALLOC;

void gravi_pca_result_delete(gravi_pca_result *self);
void gravi_pca_model_delete(gravi_pca_model *self);

#endif // GRAVI_PCA_H
