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

/**
 * @defgroup gravi_pca  PCA decomposition of visibility phase
 *
 * This module contains functions for calibrating the visibility phase data
 * using the principal component analysis method.
 * These functions are used by the gravity_pcacal recipe to calculate a PCA
 * model from calibration data, and subsequently by the gravity_vis recipe to
 * flatten visibility phase data from observations with the model.
 */

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_statistics.h>

#include "gravi_pca.h"
#include "gravi_utils.h"

/* 
 * this function copied from GSL 2.7
 * if pipeline dependency is updated, it will need to be removed from here
 */ 
static double gsl_vector_sum(const gsl_vector *a)
{
    const size_t N = a->size;
    const size_t stride = a->stride;
    double sum = 0.0;
    size_t i;
    
    for (i = 0; i < N; i++)
        sum += a->data[i * stride];
    return sum;
}

/**
 * @brief Type to hold results of a PCA decomposition. This is an opaque type.
 * @note  The PCA functions use GSL types internally so this type should only
 *        be accessed via the public interface which handles the conversions to
 *        the corresponding CPL types.
 */
struct _gravi_pca_result_
{
    int n_obs, n_wave; // matrix dimensions
    int n_valid; // number of nonzero singular values = min(n_obs, n_wave)
    gsl_matrix *data; // original data matrix
    gsl_matrix_int *mask; // mask for invalid data values
    gsl_vector *mean; // columnwise mean of original data matrix
    
    gsl_matrix *Vt; // PCA components
    gsl_vector *S;  // singular values
    gsl_vector *signs; // sign convention for PCA components, size (n_comp,)

    int n_comp; // number of components fits evaluated for
    gsl_matrix *data_fit; // matrix of fit to data, size (n_obs, n_wave)
};

/**
 * @brief Type to hold average (median) components obtained from a set of PCA
 *        decompositions and/or best-fit model to average components.
 *        This is an opaque type.
 * @note  The PCA functions use GSL types internally so this type should only
 *        be accessed via the public interface which handles the conversions to
 *        the corresponding CPL types.
 */
struct _gravi_pca_model_
{
    int n_comp; // number of components fits evaluated for
    int n_wave; // length of wavelength axis
    gsl_matrix *median_components; // median-average of components from all results, size (n_comp, n_wave)
    gsl_matrix *component_fits; // matrix of fits to each PCA component, size (n_comp, n_wave)
};

/**
 * @brief Deallocate a @c gravi_pca_result object.
 * 
 * @param self result object to free.
 */
void gravi_pca_result_delete(gravi_pca_result *self) {
    if (!self)
        return;

    FREE(gsl_matrix_free, self->data);
    FREE(gsl_matrix_int_free, self->mask);
    FREE(gsl_vector_free, self->mean);
    FREE(gsl_matrix_free, self->Vt);
    FREE(gsl_vector_free, self->S);
    FREE(gsl_vector_free, self->signs);
    FREE(gsl_matrix_free, self->data_fit);
    cpl_free(self);
}

/**
 * @brief Deallocate a @c gravi_pca_model object.
 * 
 * @param self model object to free.
 */
void gravi_pca_model_delete(gravi_pca_model *self) {
    if (!self)
        return;

    FREE(gsl_matrix_free, self->median_components);
    FREE(gsl_matrix_free, self->component_fits);
    cpl_free(self);
}

/**
 * @brief Handle computing SVD for MxN matrix, M<N.
 * 
 * @param[inout] A Pointer to the matrix to decompose. Will be replaced with the value U.
 * @param[inout] V Pointer to hold matrix of right singular vectors.
 * @param[out] S Vector of singular values.
 * 
 * @note The arguments V and S are allocated and assigned by this function.
 *       They should be deallocated using gsl_{matrix,vector}_free.
 *       The argument A is modified and may be reassigned by this function.
 */
static int svd_wrapper(gsl_matrix **A, gsl_matrix **V, gsl_vector **S)
{
    int m = (*A)->size1, n = (*A)->size2, ret;
    gsl_vector *w;

    if (m >= n) {
        *V = gsl_matrix_alloc(n, n);
        *S = gsl_vector_alloc(n);
        w = gsl_vector_alloc(n);
        ret = gsl_linalg_SV_decomp(*A /* becomes U */, *V, *S, w);
    } else {
        /* Transpose A */
        gsl_matrix *At = gsl_matrix_alloc(n, m);
        gsl_matrix_transpose_memcpy(At, *A);
        
        *V = gsl_matrix_alloc(m, m);
        *S = gsl_vector_alloc(m);
        w = gsl_vector_alloc(m);
        ret = gsl_linalg_SV_decomp(At /* becomes U */, *V, *S, w);

        /* Rearrange the result matrices: V becomes At, At becomes V */
        gsl_matrix_free(*A);
        *A = *V; *V = At;
    }

    gsl_vector_free(w);
    return ret;
}

/**
 * @brief Adjust U and V so that loadings for most significant components are positive.
 * @param U the U matrix.
 * @param V the V matrix.
 * 
 * @note See: https://github.com/scikit-learn/scikit-learn/blob/main/sklearn/utils/extmath.py#L766
 *          but note that their V is actually Vt, unlike here.
 */
static void svd_flip(gsl_matrix *U, gsl_matrix *V)
{
    /* sign of max value in each column of U */
    gsl_vector *signs = gsl_vector_alloc(U->size2);

    for (int j = 0; j < U->size2; j++) {
        gsl_vector_const_view Ucol = gsl_matrix_const_column(U, j);

        double max_abs_val = -INFINITY;
        int max_abs_ind = -1;
        for(int j = 0; j < U->size2; j++) {
            double abs_val = fabs(gsl_vector_get(&Ucol.vector, j));
            if (abs_val > max_abs_val) max_abs_ind = j;
        }

        double max_val = gsl_vector_get(&Ucol.vector, max_abs_ind);
        gsl_vector_set(signs, j, max_val > 0 ? 1 : -1);
    }
    
    /* adjust U */
    for (int i = 0; i < U->size1; i++) {
        gsl_vector_view Urow = gsl_matrix_row(U, i);
        gsl_vector_mul(&Urow.vector, signs);
    }

    /* adjust V */
    for (int i = 0; i < V->size1; i++) {
        gsl_vector_view Vrow = gsl_matrix_row(V, i);
        gsl_vector_mul(&Vrow.vector, signs);
    }

    gsl_vector_free(signs);
}

/**
 * @brief Construct a new @c gravi_pca_result object from a matrix of visphi data.
 *
 * @param data Source data matrix. The matrix should have size MxN,
 *             with M the number of observations and N the number of wavelengths.
 * @param mask Mask matrix. Should match the size of @c data. Nonzero values denote
 *             invalid elements. If NULL, all data will be used.
 * 
 * @return Newly-allocated @c gravi_pca_result object. This object must be
 *         deallocated with @c gravi_pca_result_delete.
 */
gravi_pca_result *gravi_pca_create_result(const cpl_matrix *data, const cpl_matrix *mask)
{
    cpl_ensure(data != NULL, CPL_ERROR_NULL_INPUT, NULL);

    unsigned int n_obs = cpl_matrix_get_nrow(data);
    unsigned int n_wave = cpl_matrix_get_ncol(data);
    unsigned int nvalid = CPL_MIN(n_obs, n_wave);

    /* Copy from cpl_matrix to gsl_matrix */
    gsl_matrix_const_view data_gsl_view = gsl_matrix_const_view_array(
        cpl_matrix_get_data_const(data), n_obs, n_wave);

    gsl_matrix_int *use_mask = gsl_matrix_int_alloc(n_obs, n_wave);
    if (mask) {
        for(int i = 0; i < n_obs; i++) {
            for (int j = 0; j < n_wave; j++)
                gsl_matrix_int_set(use_mask, i, j, 
                    cpl_matrix_get(mask, i, j) > 0 ? 0 : 1);
        }
    } else {
        gsl_matrix_int_set_all(use_mask, 1);
    }
    
    /* Calculate mean over columns */
    gsl_vector* column_mean = gsl_vector_alloc(n_wave);
    for(int j = 0; j < n_wave; j++) {
        gsl_vector_const_view col = gsl_matrix_const_column(&data_gsl_view.matrix, j);
        double mu = (1.0 / n_obs) * gsl_vector_sum(&col.vector);
        gsl_vector_set(column_mean, j, mu);
    }

    /* Centre data by subtracting mean */
    gsl_matrix* mean_subtracted_data = gsl_matrix_alloc(n_obs, n_wave);
    gsl_matrix_memcpy(mean_subtracted_data, &data_gsl_view.matrix);
    for(int i = 0; i < n_obs; i++) {
        gsl_vector_view mean_subtracted_row_view = gsl_matrix_row(mean_subtracted_data, i);
        gsl_vector_sub(&mean_subtracted_row_view.vector, column_mean);
    }

    gravi_pca_result *result = cpl_malloc(sizeof(gravi_pca_result));
    result->n_obs = n_obs;
    result->n_wave = n_wave;
    result->n_valid = nvalid;
    result->data = mean_subtracted_data;
    result->mask = use_mask;
    result->mean = column_mean;
    
    /* Initialise PCA data to known values */
    result->S = NULL;
    result->Vt = NULL;

    /* Initialise fitting data to known values */
    result->n_comp = -1;
    result->signs = NULL;
    result->data_fit = NULL;

    return result;
}

/**
 * @brief Perform PCA decomposition by calculating singular value decomposition.
 *
 * @param self The @c gravi_pca_result object to perform the decomposition on.
 */
cpl_error_code gravi_pca_decomp_matrix_svd(gravi_pca_result *self)
{
    cpl_ensure_code(self, CPL_ERROR_NULL_INPUT);

    unsigned int n_obs = self->n_obs;
    unsigned int n_wave = self->n_wave;
    unsigned int n_valid = self->n_valid;
    int svd_ret;

    /* Perform SVD */
    gsl_vector *S;
    gsl_matrix *V;
    gsl_matrix *U = gsl_matrix_alloc(n_obs, n_wave);
    gsl_matrix_memcpy(U, self->data);
    svd_ret = svd_wrapper(&U, &V, &S);

    if (svd_ret)
        return CPL_ERROR_ILLEGAL_INPUT;

    svd_flip(U, V);
    gsl_matrix_free(U);

    /* First nvalid singular values */
    self->S = S;

    /* explained variance of each component */
    // const int max_components = 10;
    // gsl_vector *Ssum = gsl_vector_alloc(max_components);
    // gsl_vector_set(Ssum, 0, gsl_vector_get(S, 0));
    // for (int i = 1; i < max_components; i++) {
    //     double accum = gsl_vector_get(Ssum, i - 1) + gsl_vector_get(S, i);
    //     gsl_vector_set(Ssum, i, accum);
    // }
    // gsl_vector_scale(Ssum, 1.0 / gsl_vector_get(Ssum, max_components - 1));
    // printf("Explained Variance:\n");
    // gsl_vector_fprintf(stdout, Ssum, "%f");
    // gsl_vector_free(Ssum);

    /* First nvalid rows of Vt */
    self->Vt = gsl_matrix_alloc(n_valid, n_wave);
    gsl_matrix_transpose_memcpy(self->Vt, V);
    gsl_matrix_free(V);

    return CPL_ERROR_NONE;
}

/**
 * @brief Calculate sign convention on visphi components.
 * 
 * @param visphi Array of visphi component values.
 * @param wave Array of wavelength values.
 * 
 * @return value by which components should be multiplied, either +1 or -1.
 */
static double gravi_pca_get_component_sign(const gsl_vector *visphi, const gsl_vector *wave)
{
    /* wavelength ranges where components should be positive (2 ranges, upper/lower limits) */
    const double ranges[2][2] = {{2.02, 2.04}, {2.34, 2.36}};
    int range_indices[2][2];

    /* Find start and end indices in the wavelength array for each range value */
    int nwave = wave->size, irange = 0;
    for (int w = 0; w < nwave; w++) {
        if (gsl_vector_get(wave, w) > ranges[irange / 2][irange % 2]) {
            range_indices[irange / 2][irange % 2] = w;
            ++irange;
        }
        if (irange == 4)
            break;
    }

    int range_length[2] = {
        range_indices[0][1] - range_indices[0][0],
        range_indices[1][1] - range_indices[1][0]
    };
    int total_length = range_length[0] + range_length[1];

    /* Extract visphi segments for each range */
    int offset = 0;
    gsl_vector *visphi_range = gsl_vector_alloc(total_length);
    for (int r = 0; r < 2; r++) {
        gsl_vector_const_view range = gsl_vector_const_subvector(
            visphi, range_indices[r][0], range_length[r]);
        gsl_vector_view range_dest = gsl_vector_subvector(
            visphi_range, offset, range_length[r]);
        gsl_vector_memcpy(&range_dest.vector, &range.vector);
        offset += range_length[r];
    }

    /* Check sign */
    // double result = gsl_stats_median(visphi_range->data, 1, visphi_range->size) < 0 ? -1.0 : 1.0;
    gsl_sort_vector(visphi_range);
    double result = gsl_stats_median_from_sorted_data(
        visphi_range->data, 1, visphi_range->size
    ) < 0 ? -1.0 : 1.0;

    gsl_vector_free(visphi_range);
    return result;
}

/**
 * @brief Impose sign convention on PCA components.
 * @note  The SVD of a matrix is unique up to the sign of the component values,
 *        so we arbitrarily enforce a convention that the components should be
 *        positive for 2.02 < lambda(um) < 2.04 and 2.34 < lambda(um) < 2.36.
 * 
 * @param result @c gravi_pca_result object.
 * @param wave Array of wavelength values.
 * @param num_components Number of components to compute signs for.
 */
cpl_error_code gravi_pca_set_component_signs(gravi_pca_result *self, const cpl_array *wave, const int ncomp)
{
    cpl_ensure_code(wave, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code(ncomp > 0 && ncomp <= self->n_valid, CPL_ERROR_ILLEGAL_INPUT);

    const int n_wave = self->n_wave;
    gsl_vector_const_view wave_vec = gsl_vector_const_view_array(
        cpl_array_get_data_double_const(wave), n_wave);
    gsl_vector *signs = gsl_vector_alloc(ncomp);

    for (int c = 0; c < ncomp; c++) {
        gsl_vector_const_view components = gsl_matrix_const_row(self->Vt, c);
        double sign = gravi_pca_get_component_sign(&components.vector, &wave_vec.vector);
        gsl_vector_set(signs, c, sign);
    }

    self->signs = signs;
    return CPL_ERROR_NONE;
}

/**
 * @brief Compute median values of PCA components over a set of decomposition results.
 * 
 * @param results Array of @c gravi_pca_result objects.
 * @param num_results Number of PCA results present.
 * @param num_components Number of PCA components to calculate medians for.
 * 
 * @return Newly-allocated @c gravi_pca_model object. This object must be
 *         deallocated with @c gravi_pca_model_delete.
 */
gravi_pca_model *gravi_pca_create_model(const gravi_pca_result **results, int nres, int ncomp)
{
    cpl_ensure(results, CPL_ERROR_NULL_INPUT, NULL);
    cpl_ensure(ncomp > 0 && ncomp <= results[0]->n_valid, CPL_ERROR_ILLEGAL_INPUT, NULL);
    for (int r = 0; r < nres; r++)
        cpl_ensure(results[r]->signs, CPL_ERROR_ILLEGAL_INPUT, NULL);

    const int n_wave = results[0]->n_wave;

    gsl_matrix *all_comp_vals = gsl_matrix_alloc(nres, n_wave);
    gsl_matrix *median_comp_vals = gsl_matrix_alloc(ncomp, n_wave);
    gsl_matrix *flip_vals = gsl_matrix_alloc(nres, n_wave);
    gsl_vector *flip_vec = gsl_vector_alloc(n_wave);

    for (int c = 0; c < ncomp; c++) {
        /* build temp matrix with all values for c'th component */
        for (int r = 0; r < nres; r++) {
            gsl_vector_const_view rth_components = gsl_matrix_const_row(results[r]->Vt, c);
            gsl_matrix_set_row(all_comp_vals, r, &rth_components.vector);

            /* determine flip vals */
            double sign = gsl_vector_get(results[r]->signs, c);
            gsl_vector_set_all(flip_vec, sign);
            gsl_matrix_set_row(flip_vals, r, flip_vec);
        }

        /* multiply by flip values */
        gsl_matrix_mul_elements(all_comp_vals, flip_vals);

        /* take median average down columns */
        for (int j = 0; j < n_wave; j++) {
            gsl_vector_view column_comp_vals = gsl_matrix_column(all_comp_vals, j);
            gsl_sort_vector(&column_comp_vals.vector);
            double median_val = gsl_stats_median_from_sorted_data(
                column_comp_vals.vector.data, n_wave, nres
            );
            gsl_matrix_set(median_comp_vals, c, j, median_val);
        }
    }
    gsl_matrix_free(all_comp_vals);
    gsl_matrix_free(flip_vals);
    gsl_vector_free(flip_vec);

    gravi_pca_model *model = cpl_malloc(sizeof(gravi_pca_model));
    model->n_wave = n_wave;
    model->n_comp = ncomp;
    model->median_components = median_comp_vals;
    model->component_fits = NULL;
    return model;
}

/**
 * @brief Create PCA model from existing set of components.
 * 
 * @param components Matrix of component values, size (n_components, n_wave).
 *                   First row is mean, second is primary component, etc.
 * 
 * @return Newly-allocated @c gravi_pca_model object. This object must be
 *         deallocated with @c gravi_pca_model_delete.
 */
gravi_pca_model *gravi_pca_load_model(const cpl_matrix *components)
{
    int n_comp = cpl_matrix_get_nrow(components);
    int n_wave = cpl_matrix_get_ncol(components);

    /* Copy from cpl_matrix to gsl_matrix */
    gsl_matrix_const_view comps_gsl_view = gsl_matrix_const_view_array(
        cpl_matrix_get_data_const(components), n_comp, n_wave);
    gsl_matrix *my_comps = gsl_matrix_alloc(n_comp, n_wave);
    gsl_matrix_memcpy(my_comps, &comps_gsl_view.matrix);
    
    gravi_pca_model *result = cpl_malloc(sizeof(gravi_pca_model));
    result->n_comp = n_comp;
    result->n_wave = n_wave;
    result->median_components = NULL;
    result->component_fits = my_comps;

    return result;
}

/**
 * @brief Fit B-spline model to each of a set of median-averaged PCA components.
 * 
 * @param self The @c gravi_pca_model object obtained from a previous call to @c gravi_pca_create_model.
 * @param wave Array of wavlength values.
 * @param degree Degree of fit (number of B-spline coefficients).
 * @param ncomp Number of components to fit.
 */
cpl_error_code gravi_pca_fit_components_bspline(gravi_pca_model *self, const cpl_array *wave, int degree, int ncomp)
{
    cpl_ensure_code(self, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code(self->median_components, CPL_ERROR_ILLEGAL_INPUT);

    cpl_ensure_code(wave, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code(degree > 0, CPL_ERROR_ILLEGAL_INPUT);
    cpl_ensure_code(ncomp > 0 && ncomp <= self->n_comp, CPL_ERROR_ILLEGAL_INPUT);
    
    const int n_coeffs = degree;
    const int k = 4; // k=4 => cubic spline
    const int nbreak = n_coeffs - k + 2;
    const int n_wave = self->n_wave;

    gsl_bspline_workspace *bw = gsl_bspline_alloc(k, nbreak);
    const gsl_multifit_robust_type *T = gsl_multifit_robust_default;
    gsl_multifit_robust_workspace *mw = gsl_multifit_robust_alloc(T, n_wave, n_coeffs);
    gsl_multifit_robust_maxiter(1000, mw);
    
    gsl_vector *B = gsl_vector_alloc(n_coeffs);
    gsl_vector *coeffs = gsl_vector_alloc(n_coeffs);
    gsl_matrix *cov = gsl_matrix_alloc(n_coeffs, n_coeffs);
    double wvi, Bj, fit_val, fit_err;

    gsl_bspline_knots_uniform(cpl_array_get(wave, 0, NULL), cpl_array_get(wave, n_wave-1, NULL), bw);

    gsl_matrix *result = gsl_matrix_alloc(ncomp, n_wave);

    /* Loop over components */
    for (int c = 0; c < ncomp; c++) {
        gsl_vector_const_view cv = gsl_matrix_const_row(self->median_components, c);

        /* construct fit matrix X */
        gsl_matrix *X = gsl_matrix_alloc(n_wave, n_coeffs);
        for (int i = 0; i < n_wave; i++) {
            wvi = cpl_array_get(wave, i, NULL);
            gsl_bspline_eval(wvi, B, bw);

            for (int j = 0; j < n_coeffs; j++) {
                Bj = gsl_vector_get(B, j);
                gsl_matrix_set(X, i, j, Bj);
            }
        }

        /* do fit */
        gsl_multifit_robust(X, &cv.vector, coeffs, cov, mw);

        /* store fit results */
        for (int j = 0; j < n_wave; j++) {
            wvi = cpl_array_get(wave, j, NULL);
            gsl_bspline_eval(wvi, B, bw);
            gsl_multifit_robust_est(B, coeffs, cov, &fit_val, &fit_err);
            gsl_matrix_set(result, c, j, fit_val);
        }
        gsl_matrix_free(X);
    }

    self->n_comp = ncomp;
    self->component_fits = result;
    
    gsl_vector_free(B);
    gsl_vector_free(coeffs);
    gsl_matrix_free(cov);
    gsl_bspline_free(bw);
    gsl_multifit_robust_free(mw);
    return CPL_ERROR_NONE;
}

/**
 * @brief Fit polynomial model to each of a set of median-averaged PCA components.
 * 
 * @param self The @c gravi_pca_model object obtained from a previous call to @c gravi_pca_create_model.
 * @param wave Array of wavlength values.
 * @param degree Degree of fit (order of polynomial).
 * @param ncomp Number of components to fit.
 */
cpl_error_code gravi_pca_fit_components_polynomial(gravi_pca_model *self, const cpl_array *wave, int degree, int ncomp)
{
    cpl_ensure_code(self, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code(self->median_components, CPL_ERROR_ILLEGAL_INPUT);

    cpl_ensure_code(wave, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code(degree > 0, CPL_ERROR_ILLEGAL_INPUT);
    cpl_ensure_code(ncomp > 0 && ncomp <= self->n_comp, CPL_ERROR_ILLEGAL_INPUT);
    
    const int n_wave = self->n_wave;
    const int n_coeffs = degree;

    const gsl_multifit_robust_type *T = gsl_multifit_robust_default;
    gsl_multifit_robust_workspace *mw = gsl_multifit_robust_alloc(T, n_wave, n_coeffs);
    gsl_multifit_robust_maxiter(1000, mw);
    
    gsl_vector *coeffs = gsl_vector_alloc(n_coeffs);
    gsl_matrix *cov = gsl_matrix_alloc(n_coeffs, n_coeffs);
    double wvi, Xij, fit_val, fit_err;

    gsl_matrix *result = gsl_matrix_alloc(ncomp, n_wave);

    /* Loop over components */
    for (int c = 0; c < ncomp; c++) {
        gsl_vector_const_view cv = gsl_matrix_const_row(self->median_components, c);

        /* construct fit matrix X */
        gsl_matrix *X = gsl_matrix_alloc(n_wave, n_coeffs);
        for (int i = 0; i < n_wave; i++) {
            double wvi = cpl_array_get(wave, i, NULL);

            for (int j = 0; j < n_coeffs; j++) {
                double Xij = pow(wvi, j);
                gsl_matrix_set(X, i, j, Xij);
            }
        }

        /* do fit */
        gsl_multifit_robust(X, &cv.vector, coeffs, cov, mw);

        /* store fit results */
        for (int j = 0; j < n_wave; j++) {
            gsl_vector_const_view Xi = gsl_matrix_const_row(X, j);
            gsl_multifit_robust_est(&Xi.vector, coeffs, cov, &fit_val, &fit_err);
            gsl_matrix_set(result, c, j, fit_val);
        }
        gsl_matrix_free(X);
    }

    self->n_comp = ncomp;
    self->component_fits = result;
    
    gsl_vector_free(coeffs);
    gsl_matrix_free(cov);
    gsl_multifit_robust_free(mw);
    return CPL_ERROR_NONE;
}

/**
 * @brief Override the mean component calculated from the PCA decomposition. In conjunction with the
 *  @c fit_mean_subtracted option to the @c gravi_pca_fit_model function, this can be used to obtain an improved fit.
 * 
 * @param self The @c gravi_pca_model object obtained from a previous call to @c gravi_pca_create_model.
 * @param residual The matrix of residuals with dimensions (n_obs, n_wave) equal to the original data.
 */
cpl_error_code gravi_pca_refine_mean(gravi_pca_result *self, const cpl_matrix *residual)
{
    cpl_ensure_code(self, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code(residual, CPL_ERROR_NULL_INPUT);

    int n_obs = self->n_obs;
    int n_wave = self->n_wave;

    /* Copy from cpl_matrix to gsl_matrix */
    gsl_matrix_const_view resid_gsl_view = gsl_matrix_const_view_array(
        cpl_matrix_get_data_const(residual), n_obs, n_wave);
    
    /* Calculate new mean over columns */
    gsl_vector* refined_mean = gsl_vector_alloc(n_wave);
    for(int j = 0; j < n_wave; j++) {
        gsl_vector_const_view col = gsl_matrix_const_column(&resid_gsl_view.matrix, j);
        double mu = (1.0 / n_obs) * gsl_vector_sum(&col.vector);
        gsl_vector_set(refined_mean, j, mu);
    }

    /* Re-centre the data using the new mean */
    for (int i = 0; i < n_obs; i++) {
        gsl_vector_view row = gsl_matrix_row(self->data, i);
        gsl_vector_add(&row.vector, self->mean);
        gsl_vector_sub(&row.vector, refined_mean);
    }

    /* Store the new mean */
    gsl_vector_swap(self->mean, refined_mean);
    gsl_vector_free(refined_mean);

    return CPL_ERROR_NONE;
}

/**
 * @brief Get median-averaged component from a set of PCA decompositions.
 * 
 * @param self The @c gravi_pca_model object obtained from a previous call to @c gravi_pca_create_model.
 * @param component The component to return, indexed from 1.
 *
 * @return Vector of components with length n_wave.
 */
cpl_vector *gravi_pca_get_component_median(const gravi_pca_model *self, int component)
{
    cpl_ensure(self, CPL_ERROR_NULL_INPUT, NULL);
    cpl_ensure(self->median_components, CPL_ERROR_ILLEGAL_INPUT, NULL);
    cpl_ensure(component > 0 && component <= self->n_comp, CPL_ERROR_ILLEGAL_INPUT, NULL);

    cpl_vector *result = cpl_vector_new(self->n_wave);
    gsl_vector_const_view row = gsl_matrix_const_row(self->median_components, component-1);
    for(int i = 0; i < self->n_wave; i++)
        cpl_vector_set(result, i, gsl_vector_get(&row.vector, i));
    return result;
}

/**
 * @brief Get fit to components from PCA decomposition.
 * 
 * @param self The @c gravi_pca_model object obtained from a previous call to @c gravi_pca_create_model.
 * One of the functions @c gravi_pca_fit_components_{polynomial,spline} must have previously been called.
 * @param component The component to return the fit for, indexed from 1.
 *
 * @return Vector of fitted values with length n_wave.
 */
cpl_vector *gravi_pca_get_component_fit(const gravi_pca_model *self, int component)
{
    cpl_ensure(self, CPL_ERROR_NULL_INPUT, NULL);
    cpl_ensure(self->component_fits, CPL_ERROR_ILLEGAL_INPUT, NULL);
    cpl_ensure(component > 0 && component <= self->n_comp, CPL_ERROR_ILLEGAL_INPUT, NULL);

    cpl_vector *result = cpl_vector_new(self->n_wave);
    gsl_vector_const_view row = gsl_matrix_const_row(self->component_fits, component-1);
    for(int i = 0; i < self->n_wave; i++)
        cpl_vector_set(result, i, gsl_vector_get(&row.vector, i));
    return result;
}

/**
 * @brief Type to hold parameters for fitting PCA components to data.
 */
typedef struct _gravi_pca_model_params_ {
    int n_data; // number of data values 
    int n_components; // number of PCA components in fit
    const gsl_matrix *components; // component fit values
    gsl_vector_const_view data; // data vector
    gsl_vector_int_const_view mask; // mask vector
} gravi_pca_model_params;

/**
 * @brief Evaluate chi^2 statistic for fit.
 * @param x Vector of coefficients C_i.
 * @param params Fit parameter values, of actual type @c gravi_pca_model_params.
 */
static double gravi_pca_model_chi2(const gsl_vector *x, void *params)
{
    gravi_pca_model_params *p = params;
    gsl_vector *model = gsl_vector_calloc(p->n_data);
    gsl_vector *scratch = gsl_vector_alloc(p->n_data);
    gsl_vector *chi2 = gsl_vector_alloc(p->n_data);
    
    /* model = C1 * c1vals + C2 * c2vals + ... */
    for (int i = 0; i < p->n_components; i++) {
        gsl_vector_const_view cv = gsl_matrix_const_row(p->components, i);
        gsl_vector_memcpy(scratch, &cv.vector);
        gsl_vector_scale(scratch, gsl_vector_get(x, i));
        gsl_vector_add(model, scratch);
    }

    /* form chi2 */
    gsl_vector_memcpy(chi2, &(p->data.vector));
    gsl_vector_sub(chi2, model);
    gsl_vector_mul(chi2, chi2);

    int n_valid = 0;
    double chi2_sum = 0.0;
    for (int i = 0; i < p->n_data; i++) {
        if(gsl_vector_int_get(&(p->mask.vector), i) > 0) {
            chi2_sum += gsl_vector_get(chi2, i);
            n_valid++;
        }
    }

    /* Check for no valid data */
    if (n_valid <= p->n_components)
        return -1;

    chi2_sum /= (n_valid - p->n_components);

    gsl_vector_free(chi2);
    gsl_vector_free(scratch);
    gsl_vector_free(model);
    return chi2_sum;
}

/**
 * @brief Fit model formed from linear combination of PCA components to data.
 * @param self The @c gravi_pca_result object containing data to fit.
 * @param model The @c gravi_pca_model object containing the model PCA component to fit with.
 * @param fit_mean_subtracted If True, fit to the mean-subtracted data, otherwise add the mean back on first.
 */
cpl_error_code gravi_pca_fit_model(gravi_pca_result *self, const gravi_pca_model *model, cpl_boolean fit_mean_subtracted, cpl_boolean verbose)
{    
    cpl_ensure_code(self, CPL_ERROR_NULL_INPUT);
    if (fit_mean_subtracted)
        cpl_ensure_code(self->mean, CPL_ERROR_ILLEGAL_INPUT);
    
    cpl_ensure_code(model, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code(model->component_fits, CPL_ERROR_ILLEGAL_INPUT);

    const int MAX_ITERATIONS = 1000;
    const int n_obs = self->n_obs;
    const int n_wave = model->n_wave;
    const int n_comp = model->n_comp;

    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2rand;
    gsl_multimin_fminimizer *mini = gsl_multimin_fminimizer_alloc(T, n_comp);
    int status;

    gsl_vector *coeffs = gsl_vector_alloc(n_comp);
    gsl_vector *step_size = gsl_vector_alloc(n_comp);

    /* Add mean back on to mean-subtracted data */
    gsl_matrix *data = gsl_matrix_calloc(n_obs, n_wave);
    if (!fit_mean_subtracted) {
        for (int i = 0; i < n_obs; i++)
            gsl_matrix_set_row(data, i, self->mean);
    }
    gsl_matrix_add(data, self->data);

    self->data_fit = gsl_matrix_calloc(n_obs, n_wave);
    self->n_comp = n_comp;

    /* Loop over rows in the original data */
    for (int i = 0; i < n_obs; i++) {
        gsl_vector_set_all(coeffs, 1.0);
        gsl_vector_set_all(step_size, 0.1);

        /* Prepare container of function parameters */
        gravi_pca_model_params params = {
            .n_data = n_wave, .n_components = n_comp,
            .components = model->component_fits,
            .data = gsl_matrix_const_row(data, i),
            .mask = gsl_matrix_int_const_row(self->mask, i)
        };

        /* Prepare minimization function object */
        gsl_multimin_function func = {
            .f = &gravi_pca_model_chi2,
            .n = n_comp,
            .params = &params
        };

        gsl_multimin_fminimizer_set(mini, &func, coeffs, step_size);

        int iterations = 0;
        double size;
        if (verbose)
            printf("iterations c_0 c_1 size chi2\n");
        do {
            iterations++;
            status = gsl_multimin_fminimizer_iterate(mini);

            /* Error taking optimisation step */
            if (status)
                break;
            
            /* All data was invalid */
            if (mini->fval < 0) {
                cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT, "No valid visphi data for fitting\n");
                status = GSL_EINVAL;
                break;
            }

            size = gsl_multimin_fminimizer_size(mini);
            status = gsl_multimin_test_size(size, 1e-3);
            
            if (verbose)
                printf("%d %f %f %f %f\n", iterations, mini->x->data[0], mini->x->data[1], size, mini->fval);
        } while (status == GSL_CONTINUE && iterations < MAX_ITERATIONS);

        /* Store best-fit model */
        if (status == GSL_CONTINUE) {
            cpl_msg_warning(cpl_func, "model-fitting did not converge");
        } else if (status != GSL_SUCCESS) {
            cpl_msg_error(cpl_func, "model-fitting failed with error code %d (%s)", status, gsl_strerror(status));
            break;
        }
        gsl_vector_view fit_vals = gsl_matrix_row(self->data_fit, i);
        for (int c = 0; c < n_comp; c++) {
            gsl_vector *term = gsl_vector_alloc_row_from_matrix(model->component_fits, c);
            gsl_vector_scale(term, gsl_vector_get(mini->x, c));
            gsl_vector_add(&fit_vals.vector, term);
            gsl_vector_free(term);
        }
        
        if (fit_mean_subtracted) {
            gsl_vector_add(&fit_vals.vector, self->mean);
        }
    }

    gsl_matrix_free(data);
    gsl_vector_free(coeffs);
    gsl_vector_free(step_size);
    gsl_multimin_fminimizer_free(mini);
    return status ? CPL_ERROR_ILLEGAL_OUTPUT : CPL_ERROR_NONE;
}

/**
 * @brief Get components from PCA decomposition.
 * 
 * @param self The @c gravi_pca_result object obtained from a previous decomposition.
 * @param component The component to return. 0 will return the mean, 1 the first PCA component, etc.
 *
 * @return Vector of components with length N, the number of observations in the original data.
 */
cpl_vector *gravi_pca_get_component(const gravi_pca_result *self, int component)
{
    cpl_ensure(self, CPL_ERROR_NULL_INPUT, NULL);
    cpl_ensure(component >= 0 && component <= self->n_valid, CPL_ERROR_ILLEGAL_INPUT, NULL);

    cpl_vector *result = cpl_vector_new(self->n_wave);
    if (component == 0) {
        /* Return the mean */
        for(int i = 0; i < self->n_wave; i++)
            cpl_vector_set(result, i, gsl_vector_get(self->mean, i));
    } else {
        /* Return the n'th component */
        gsl_vector_const_view row = gsl_matrix_const_row(self->Vt, component-1);
        for(int i = 0; i < self->n_wave; i++) {
            cpl_vector_set(result, i, gsl_vector_get(&row.vector, i));
        }
        cpl_vector_multiply_scalar(result, gsl_vector_get(self->signs, component-1));
    }
    return result;
}

/**
 * @brief Get noise-free model.
 * 
 * @param self The @c gravi_pca_result object obtained from a previous decomposition.
 *
 * @return Matrix of fitted values with same dimensions as original data.
 */
cpl_matrix *gravi_pca_get_data_fit(const gravi_pca_result *self)
{
    cpl_ensure(self, CPL_ERROR_NULL_INPUT, NULL);
    cpl_ensure(self->data_fit, CPL_ERROR_NULL_INPUT, NULL);

    cpl_matrix *result = cpl_matrix_new(self->n_obs, self->n_wave);
    for(int i = 0; i < self->n_obs; i++) {
        for(int j = 0; j < self->n_wave; j++) {
            cpl_matrix_set(result, i, j, gsl_matrix_get(self->data_fit, i, j));
        }
    }
    return result;
}

/**
 * @brief Get residual (data - model).
 * 
 * @param self The @c gravi_pca_result object obtained from a previous decomposition.
 * 
 * @return Matrix of residual values with same dimensions as original data.
*/
cpl_matrix *gravi_pca_get_data_residual(const gravi_pca_result *self)
{
    cpl_ensure(self, CPL_ERROR_NULL_INPUT, NULL);
    cpl_ensure(self->data, CPL_ERROR_NULL_INPUT, NULL);
    cpl_ensure(self->mean, CPL_ERROR_NULL_INPUT, NULL);
    cpl_ensure(self->data_fit, CPL_ERROR_NULL_INPUT, NULL);

    /* Add mean back onto centred data, then subtract model */
    gsl_matrix *resid = gsl_matrix_alloc(self->n_obs, self->n_wave);
    gsl_matrix_memcpy(resid, self->data);

    for(int i = 0; i < self->n_obs; i++) {
        gsl_vector_view row_view = gsl_matrix_row(resid, i);
        gsl_vector_add(&row_view.vector, self->mean);
    }

    gsl_matrix_sub(resid, self->data_fit);

    cpl_matrix *result = cpl_matrix_new(self->n_obs, self->n_wave);
    for(int i = 0; i < self->n_obs; i++) {
        for(int j = 0; j < self->n_wave; j++) {
            cpl_matrix_set(result, i, j, gsl_matrix_get(resid, i, j));
        }
    }
    return result;
}