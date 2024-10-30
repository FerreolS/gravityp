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
 * @defgroup gravi_astrometry  Astrometric Solutions
 *
 * This module contains functions for computing phase referencing for dual-
 * field astrometric observations. These functions are used by the
 * gravity_astrometry recipe, and are applicable to on-axis, off-axis and swap
 * observing strategies.
 */

#include "gravi_astrometry.h"
#include "gravi_dfs.h"
#include "gravi_pfits.h"
#include "gravi_utils.h"

#include <assert.h>
#include <ctype.h>
#include <math.h>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_multimin.h>

#define PI 3.14159265359
#define TWOPI 6.28318530718
#define MAS_TO_RAD (2.0 * PI / (3600 * 360 * 1000))

/*-----------------------------------------------------------------------------
                              Private prototypes
 -----------------------------------------------------------------------------*/

struct _astrometry_data_ {
    char *filename, *insname;
    int nchannel, ndit, nwave, nwave_ft;
    double sobj_x, sobj_y;
    double dit;
    double mjd;

    double swap_astrometry_guess[2];
    double swap_astrometry_fit[2];
    cpl_boolean swap;

    gsl_vector *wave;
    gsl_vector *u, *v;
    gsl_matrix *ucoord, *vcoord;

    gsl_matrix_complex *visdata;
    gsl_matrix_complex *viserr;
    gsl_matrix_complex *visdata_ft;

    int nflag;
    gsl_matrix_int *flag;
    // gsl_matrix_int **flag_cov;
    
    gsl_matrix *phase_ref;
    gsl_matrix *opd_disp;
    gsl_matrix *phase_met_telfc;
    gsl_matrix_complex *vis_ref;

    gsl_vector **vis_ref_cov;
    gsl_vector_complex **vis_ref_pcov;

    gsl_matrix *amp_ref_astro;
    gsl_matrix_complex *phase_ref_astro;
};

typedef enum _find_closest_mode {
    FIND_MODE_NONE,
    FIND_MODE_BEFORE,
    FIND_MODE_AFTER
} find_closest_mode;

static gsl_vector *load_vector(cpl_table *table, const char *name) CPL_ATTR_ALLOC;
static gsl_matrix *load_matrix(cpl_table *table, const char *name) CPL_ATTR_ALLOC;
static gsl_matrix_complex *load_matrix_complex(cpl_table *table, const char *name) CPL_ATTR_ALLOC;

static cpl_error_code gravi_astrometry_scale_visibilities(astro_data *self, double factor);
static cpl_error_code gravi_astrometry_mul_visibilities(astro_data *self, const gsl_vector *factor);
static cpl_error_code gravi_astrometry_add_phase(astro_data *self, const gsl_matrix *phase);
static cpl_error_code gravi_astrometry_recentre_phase(astro_data *self, double xvalue, double yvalue);
static cpl_size gravi_astrometry_find_closest_mjd(astro_data *self, astro_data **others, cpl_size n_other, find_closest_mode mode);

static gsl_vector *average_vector_over_dits(gsl_vector *v, int ndit, int nchannel) CPL_ATTR_ALLOC;
static gsl_matrix *average_matrix_over_dits(gsl_matrix *m, int ndit, int nchannel, int nwave, gsl_matrix_int *flag) CPL_ATTR_ALLOC;
static gsl_matrix_complex *average_matrix_complex_over_dits(gsl_matrix_complex *m, int ndit, int nchannel, int nwave, gsl_matrix_int *flag) CPL_ATTR_ALLOC;

static astro_data *gravi_astrometry_average_over_dits(astro_data *self) CPL_ATTR_ALLOC;

typedef struct _gravi_astrometry_model_params_ {
    astro_data ** group1, ** group2;
    int n1, n2;
} gravi_astrometry_model_params;

static gsl_matrix_complex *gravi_astrometry_calculate_visref_swap(astro_data **data, cpl_size ndata, double ra, double dec) CPL_ATTR_ALLOC;
static gsl_matrix_complex *gravi_astrometry_create_swap_reference(astro_data **swap_data, cpl_size nswap) CPL_ATTR_ALLOC;
static double gravi_astrometry_calculate_chi2(const gsl_vector *X, void *params);
static cpl_error_code gravi_astrometry_minimise_chi2_grid(gravi_astrometry_model_params *params, const gsl_vector *ra_grid, const gsl_vector *dec_grid, gsl_vector *ra_dec_out, gsl_matrix *chi2_map);
static cpl_error_code gravi_astrometry_minimise_chi2_descent(gravi_astrometry_model_params *params, double ra_guess, double dec_guess, gsl_vector *Xsolve);

/*-----------------------------------------------------------------------------
                              Function code
 -----------------------------------------------------------------------------*/

static int _gsl_vector_int_sum(const gsl_vector_int *a)
{
    const size_t N = a->size;
    const size_t stride = a->stride;
    int sum = 0.0;
    size_t i;
    
    for (i = 0; i < N; i++)
        sum += a->data[i * stride];
    return sum;
}

/**
 * @brief Load data from table into GSL vector.
*/
static gsl_vector *load_vector(cpl_table *table, const char *name)
{
    cpl_ensure(table != NULL, CPL_ERROR_NULL_INPUT, NULL);   
    cpl_ensure(name != NULL, CPL_ERROR_ILLEGAL_INPUT, NULL);

    cpl_size depth = cpl_table_get_column_depth(table, name);
    if (depth == - 1) {
        cpl_error_set_message(cpl_func, CPL_ERROR_DATA_NOT_FOUND, "Field %s not found", name);
        return NULL;
    } else if (depth > 0) {
        cpl_error_set_message(cpl_func, CPL_ERROR_INVALID_TYPE, "Field %s not of scalar type", name);
        return NULL;
    }

    cpl_size nelem = cpl_table_get_nrow(table);
    double * cpl_vec = cpl_table_get_data_double(table, name);

    gsl_vector *gsl_vec = gsl_vector_alloc(nelem);
    memcpy(gsl_vec->data, cpl_vec, nelem * sizeof(double));

    return gsl_vec;
}

/**
 * @brief Load data from table into GSL matrix. Each table row is copied into the corresponding matrix row.
*/
static gsl_matrix *load_matrix(cpl_table *table, const char *name)
{
    cpl_ensure(table != NULL, CPL_ERROR_NULL_INPUT, NULL);   
    cpl_ensure(name != NULL, CPL_ERROR_ACCESS_OUT_OF_RANGE, NULL);

    cpl_size depth = cpl_table_get_column_depth(table, name);
    if (depth == -1) {
        cpl_error_set_message(cpl_func, CPL_ERROR_DATA_NOT_FOUND, "Field %s not found", name);
        return NULL;
    } else if (depth == 0) {
        cpl_error_set_message(cpl_func, CPL_ERROR_INVALID_TYPE, "Field %s not of 2D array type", name);
        return NULL;
    }

    cpl_size nrow = cpl_table_get_nrow(table);
    cpl_size ncol = cpl_table_get_column_dimension(table, name, 0);
    const cpl_array ** cpl_matrix = cpl_table_get_data_array_const(table, name);

    gsl_matrix *gsl_matrix = gsl_matrix_alloc(nrow, ncol);
    for (int i = 0; i < nrow; i++) {
        gsl_vector_const_view row_view = gsl_vector_const_view_array(
            cpl_array_get_data_double_const(cpl_matrix[i]), ncol);
        gsl_matrix_set_row(gsl_matrix, i, &row_view.vector);
    }

    return gsl_matrix;
}

/**
 * @brief Load data from table into GSL matrix. Each table row is copied into the corresponding matrix row.
*/
static gsl_matrix_complex *load_matrix_complex(cpl_table *table, const char *name)
{
    cpl_ensure(table != NULL, CPL_ERROR_NULL_INPUT, NULL);   
    cpl_ensure(name != NULL, CPL_ERROR_ACCESS_OUT_OF_RANGE, NULL);

    cpl_size depth = cpl_table_get_column_depth(table, name);
    if (depth == - 1) {
        cpl_error_set_message(cpl_func, CPL_ERROR_DATA_NOT_FOUND, "Field %s not found", name);
        return NULL;
    } else if (depth == 0) {
        cpl_error_set_message(cpl_func, CPL_ERROR_INVALID_TYPE, "Field %s not of 2D array type", name);
        return NULL;
    }

    cpl_size nrow = cpl_table_get_nrow(table);
    cpl_size ncol = cpl_table_get_column_dimension(table, name, 0);
    const cpl_array ** cpl_matrix = cpl_table_get_data_array_const(table, name);

    gsl_matrix_complex *gsl_matrix = gsl_matrix_complex_alloc(nrow, ncol);
    for (int i = 0; i < nrow; i++) {
        for (int j = 0; j < ncol; j++) {
            double complex cpl_matrix_val = cpl_array_get_complex(cpl_matrix[i], j, NULL);
            gsl_complex gsl_matrix_val = gsl_complex_rect(creal(cpl_matrix_val), cimag(cpl_matrix_val));
            gsl_matrix_complex_set(gsl_matrix, i, j, gsl_matrix_val);
        }
    }

    return gsl_matrix;
}

/**
 * @brief Scale VISDATA, VISREF, errors by scalar factor.
 * 
 * @param self astro_data to modify inplace.
 * @param factor scaling factor.
*/
static cpl_error_code gravi_astrometry_scale_visibilities(astro_data *self, double factor)
{
    cpl_ensure_code(self, CPL_ERROR_NULL_INPUT);

    gsl_complex complex_factor = gsl_complex_rect(factor, 0);
    gsl_matrix_complex_scale(self->visdata, complex_factor);
    gsl_matrix_complex_scale(self->visdata_ft, complex_factor);
    gsl_matrix_complex_scale(self->vis_ref, complex_factor);

    double factor2 = factor * factor;
    for (int i = 0; i < self->ndit * self->nchannel; i++) {
        gsl_vector_scale(self->vis_ref_cov[i], factor2);
        gsl_vector_complex_scale(self->vis_ref_pcov[i], gsl_complex_rect(factor2, 0));
    }

    return CPL_ERROR_NONE;
}

/**
 * @brief Scale VISDATA, VISREF, errors by vector factor, broadcasting over wavelength axis.
 * 
 * @param self astro_data to modify inplace.
 * @param factor vector of scaling factors, of size self->ndit * self->nchannel.
*/
static cpl_error_code gravi_astrometry_mul_visibilities(astro_data *self, const gsl_vector *factor)
{
    cpl_ensure_code(self, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code(factor, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code((size_t) (self->ndit * self->nchannel) == factor->size, CPL_ERROR_INCOMPATIBLE_INPUT);

    gsl_matrix_complex *complex_factor = gsl_matrix_complex_alloc(self->ndit * self->nchannel, self->nwave);
    gsl_matrix_complex *complex_factor_ft = gsl_matrix_complex_alloc(self->ndit * self->nchannel, self->nwave_ft);
    for (int i = 0; i < self->ndit * self->nchannel; i++) {
        gsl_complex cf = gsl_complex_rect(gsl_vector_get(factor, i), 0);
        for (int j = 0; j < self->nwave; j++)
            gsl_matrix_complex_set(complex_factor, i, j, cf);
        for (int j = 0; j < self->nwave_ft; j++)
            gsl_matrix_complex_set(complex_factor_ft, i, j, cf);
    }

    gsl_matrix_complex_mul_elements(self->visdata, complex_factor);
    gsl_matrix_complex_mul_elements(self->vis_ref, complex_factor);
    gsl_matrix_complex_mul_elements(self->visdata_ft, complex_factor_ft);

    for (int i = 0; i < self->ndit * self->nchannel; i++) {
        double factor2 = gsl_vector_get(factor, i);
        factor2 *= factor2;
        gsl_vector_scale(self->vis_ref_cov[i], factor2);
        gsl_vector_complex_scale(self->vis_ref_pcov[i], gsl_complex_rect(factor2, 0));
    }

    FREE(gsl_matrix_complex_free, complex_factor);
    FREE(gsl_matrix_complex_free, complex_factor_ft);

    return CPL_ERROR_NONE;
}

/**
 * @brief Add phase to visRef, updating errors accordingly
*/
static cpl_error_code gravi_astrometry_add_phase(astro_data *self, const gsl_matrix *phase)
{
    gsl_matrix_complex *phi = gsl_matrix_complex_alloc(self->ndit * self->nchannel, self->nwave);
    for (int i = 0; i < self->ndit * self->nchannel; i++) {
        for (int j = 0; j < self->nwave; j++) {
            gsl_complex phi_val = gsl_complex_polar(1.0, gsl_matrix_get(phase, i, j));
            gsl_matrix_complex_set(phi, i, j, phi_val);
        }
    }

    gsl_matrix_complex_mul_elements(self->vis_ref, phi);

    gsl_vector_complex *cov = gsl_vector_complex_alloc(self->nwave);
    gsl_vector_complex *pcov = gsl_vector_complex_alloc(self->nwave);
    gsl_vector_complex *tmp = gsl_vector_complex_alloc(self->nwave);
    
    for (int i = 0; i < self->ndit * self->nchannel; i++) {
        gsl_vector_complex_set_zero(cov);
        gsl_vector_complex_set_zero(pcov);
        gsl_vector_complex_set_zero(tmp);

        // Conjugate vector of phase
        gsl_vector_complex_const_view phi_row = gsl_matrix_complex_const_row(phi, i);
        gsl_vector_complex_memcpy(tmp, &phi_row.vector);
        gsl_vector_view phi_imag = gsl_vector_complex_imag(tmp);
        gsl_vector_scale(&phi_imag.vector, -1);

        // Need complex vector of C
        gsl_vector_view cov_real = gsl_vector_complex_real(cov);
        gsl_vector_memcpy(&cov_real.vector, self->vis_ref_cov[i]);

        // Adjust Cov
        gsl_vector_complex_mul(tmp, cov); // D = conj(phi) * C
        gsl_vector_complex_mul(tmp, &phi_row.vector); // D * phi
        cov_real = gsl_vector_complex_real(tmp);
        gsl_vector_memcpy(self->vis_ref_cov[i], &cov_real.vector);

        // Adjust Pcov
        gsl_vector_complex_memcpy(pcov, self->vis_ref_pcov[i]);
        gsl_vector_complex_mul(pcov, &phi_row.vector);
        gsl_vector_complex_mul(pcov, &phi_row.vector);
        gsl_vector_complex_memcpy(self->vis_ref_pcov[i], pcov);
    }

    FREE(gsl_vector_complex_free, cov);
    FREE(gsl_vector_complex_free, pcov);
    FREE(gsl_vector_complex_free, tmp);
    FREE(gsl_matrix_complex_free, phi);

    return CPL_ERROR_NONE;
}

/**
 * @brief Recentre visRef and visRefError on the given position (in mas) by shifting the phases.
*/
static cpl_error_code gravi_astrometry_recentre_phase(astro_data *self, double xvalue, double yvalue)
{
    gsl_matrix *this_u = gsl_matrix_alloc(self->ndit * self->nchannel, self->nwave);
    gsl_matrix_memcpy(this_u, self->ucoord);
    gsl_matrix_scale(this_u, xvalue);
    
    gsl_matrix *this_v = gsl_matrix_alloc(self->ndit * self->nchannel, self->nwave);
    gsl_matrix_memcpy(this_v, self->vcoord);
    gsl_matrix_scale(this_v, yvalue);

    gsl_matrix_add(this_u, this_v);
    gsl_matrix_scale(this_u, TWOPI * MAS_TO_RAD);

    gravi_astrometry_add_phase(self, this_u);
    
    FREE(gsl_matrix_free, this_u);
    FREE(gsl_matrix_free, this_v);

    return CPL_ERROR_NONE;
}

/**
 * @brief Filter based on FT flux threshold and normalise.
*/
cpl_error_code gravi_astrometry_filter_ftflux(astro_data *self, double threshold)
{
    cpl_ensure_code(self, CPL_ERROR_NULL_INPUT);

    int npoints = 0;
    for (int i = 0; i < self->ndit * self->nchannel; i++) {
        gsl_vector_complex_const_view visft_row = gsl_matrix_complex_const_row(self->visdata_ft, i);
        double row_abs_sum = 0.0;
        for (int j = 0; j < self->nwave_ft; j++)
            row_abs_sum += gsl_complex_abs(gsl_vector_complex_get(&visft_row.vector, j));
        row_abs_sum /= self->nwave_ft;

        if (row_abs_sum < threshold) {
            npoints += 233;
            for (int k = 0; k < self->nwave; k++) {
                gsl_matrix_int_set(self->flag, i, k, 1);
                gsl_matrix_complex_set(self->visdata, i, k, GSL_COMPLEX_ZERO);
                // for (int k = 0; k < self->nwave; k++) {
                //     gsl_matrix_int_set(self->flag_cov[i], j, k, 1);
                //     gsl_matrix_int_set(self->flag_cov[i], k, j, 1);
                // }
            }
        }
        // if (row_abs_sum < threshold) {
        //     npoints += 233;
        //     int dit = i / self->nchannel;
        //     // int cha = i % self->nchannel;
        //     // printf("%d %d %d\n", i, dit, cha);
        //     for (int j = 0; j < self->nchannel; j++) {
        //         int z = j + dit * self->nchannel;
        //         for (int k = 0; k < self->nwave; k++) {
        //             gsl_matrix_int_set(self->flag, z, k, 1);
        //             gsl_matrix_complex_set(self->visdata, z, k, GSL_COMPLEX_ZERO);
        //             // for (int k = 0; k < self->nwave; k++) {
        //             //     gsl_matrix_int_set(self->flag_cov[i], j, k, 1);
        //             //     gsl_matrix_int_set(self->flag_cov[i], k, j, 1);
        //             // }
        //         }
        //     }
        // }
    }
    cpl_msg_debug(cpl_func, "Flagging %d points below FT threshold %.3e", npoints, threshold);
    self->nflag += npoints;

    return CPL_ERROR_NONE;
}

/**
 * @brief Normalise visibilities to average FT flux.
*/
cpl_error_code gravi_astrometry_normalise_to_ft(astro_data *self)
{
    gsl_vector *ftflux_scale = gsl_vector_alloc(self->ndit * self->nchannel);
    for (int i = 0; i < self->ndit * self->nchannel; i++) {
        double scale = 0;
        for (int j = 0; j < self->nwave_ft; j++)
            scale += gsl_complex_abs(gsl_matrix_complex_get(self->visdata_ft, i, j));
        gsl_vector_set(ftflux_scale, i, self->nwave / scale);
    }

    gravi_astrometry_mul_visibilities(self, ftflux_scale);
    FREE(gsl_vector_free, ftflux_scale);
    return CPL_ERROR_NONE;
}

/**
 * @brief Load data for astrometry from a gravi_data.
 */
astro_data *gravi_astrometry_load(gravi_data *data)
{
    cpl_ensure(data != NULL, CPL_ERROR_NULL_INPUT, NULL);

    cpl_propertylist *hdr = gravi_data_get_header(data);

    astro_data *self = cpl_malloc(sizeof(astro_data));

    const char *filename = cpl_propertylist_get_string(hdr, "PIPEFILE");
    self->filename = cpl_strdup(filename);

    const char *insname = cpl_propertylist_get_string(hdr, "TELESCOP");
    self->insname = cpl_strdup(insname);

    self->nchannel = 6;
    if (cpl_propertylist_has(hdr, "ESO TPL NDIT OBJECT"))
        self->ndit = cpl_propertylist_get_int(hdr, "ESO TPL NDIT OBJECT");
    else
        self->ndit = cpl_propertylist_get_int(hdr, "ESO TPL NDIT SKY");
    self->dit = cpl_propertylist_get_double(hdr, "ESO DET2 SEQ1 DIT");

    self->swap = (!strcmp(cpl_propertylist_get_string(hdr, "ESO INS SOBJ SWAP"), "YES"));
    cpl_msg_debug(cpl_func, "SWAP keyword is: %s", cpl_propertylist_get_string(hdr, "ESO INS SOBJ SWAP"));
    // self->swap_astrometry_guess = NULL;
    // self->swap_astrometry_fit = NULL;

    self->sobj_x = cpl_propertylist_get_double(hdr, "ESO INS SOBJ X");
    self->sobj_y = cpl_propertylist_get_double(hdr, "ESO INS SOBJ Y");
    cpl_msg_debug(cpl_func, "Fiber position is: (%f, %f)", self->sobj_x, self->sobj_y);

    self->mjd = cpl_propertylist_get_double(hdr, "MJD-OBS");

    /* Assigned by gravi_astrometry_create_phase_reference */
    self->amp_ref_astro = NULL;
    self->phase_ref_astro = NULL;

    /* Get size of wavelength axis */
    cpl_propertylist *wave_plist = gravi_data_get_oi_wave_plist(data, GRAVI_SC, 0, 1);
    self->nwave = cpl_propertylist_get_int(wave_plist, "NWAVE");

    wave_plist = gravi_data_get_oi_wave_plist(data, GRAVI_FT, 0, 1);
    self->nwave_ft = cpl_propertylist_get_int(wave_plist, "NWAVE");

    /* Get wavelength and cast to double */
    cpl_table *wave_table = gravi_data_get_oi_wave(data, GRAVI_SC, 0, 1);
    cpl_table_cast_column(wave_table, "EFF_WAVE", "EFF_WAVE", CPL_TYPE_DOUBLE);
    self->wave = load_vector(wave_table, "EFF_WAVE");

    /* Load visibility data from VIS_OI */
    cpl_table *oivis_table = gravi_data_get_oi_vis(data, GRAVI_SC, 0, 1);

    self->u = load_vector(oivis_table, "UCOORD");
    self->v = load_vector(oivis_table, "VCOORD");

    /* uv divided by wavelength is more generally useful */
    self->ucoord = gsl_matrix_alloc(self->ndit * self->nchannel, self->nwave);
    self->vcoord = gsl_matrix_alloc(self->ndit * self->nchannel, self->nwave);
    for (int j = 0; j < self->nwave; j++) {
        gsl_vector_view col_view;
        col_view = gsl_matrix_column(self->ucoord, j);
        gsl_vector_memcpy(&col_view.vector, self->u);
        col_view = gsl_matrix_column(self->vcoord, j);
        gsl_vector_memcpy(&col_view.vector, self->v);
    }

    for (int i = 0; i < self->ndit * self->nchannel; i++) {
        gsl_vector_view row_view;
        row_view = gsl_matrix_row(self->ucoord, i);
        gsl_vector_div(&row_view.vector, self->wave);
        row_view = gsl_matrix_row(self->vcoord, i);
        gsl_vector_div(&row_view.vector, self->wave);
    }

    self->visdata = load_matrix_complex(oivis_table, "VISDATA");
    self->viserr = load_matrix_complex(oivis_table, "VISERR");
    self->visdata_ft = load_matrix_complex(oivis_table, "VISDATA_FT");

    self->phase_ref = load_matrix(oivis_table, "PHASE_REF");
    self->opd_disp = load_matrix(oivis_table, "OPD_DISP");
    self->phase_met_telfc = load_matrix(oivis_table, "PHASE_MET_TELFC");

    /* Load FLAG */
    const cpl_array ** cpl_flag = cpl_table_get_data_array_const(oivis_table, "FLAG");
    const int* cpl_rej = NULL;
    if (cpl_table_has_column(oivis_table, "REJECTION_FLAG"))
        cpl_rej = cpl_table_get_data_int_const(oivis_table, "REJECTION_FLAG");
    
    self->nflag = 0;
    self->flag = gsl_matrix_int_alloc(self->ndit * self->nchannel, self->nwave);
    for (int i = 0; i < self->ndit * self->nchannel; i++) {
        // int rej_val = ((cpl_rej ? cpl_rej[i] : 0) & 19) > 0;
        int rej_val = (cpl_rej ? cpl_rej[i] : 0) > 0;
        for (int j = 0; j < self->nwave; j++) {
            int flag_val = cpl_array_get_int(cpl_flag[i], j, NULL);
            flag_val |= rej_val;
            gsl_matrix_int_set(self->flag, i, j, flag_val);
            self->nflag += flag_val;
        }
    }

    /* Clean flagged visData */
    int ncleaned = 0;
    for (int i = 0; i < self->ndit * self->nchannel; i++) {
        for (int j = 0; j < self->nwave; j++) {
            if (gsl_matrix_int_get(self->flag, i, j)) {
                ++ncleaned;
                gsl_matrix_complex_set(self->visdata, i, j, GSL_COMPLEX_ZERO);
            }
        }
    }
    cpl_msg_debug(cpl_func, "Cleaned %d flagged values in VISDATA", ncleaned);

    /* Construct flags for covariance matrices (unused) */
    // self->flag_cov = cpl_malloc(self->ndit * self->nchannel * sizeof(gsl_matrix_int*));
    // for (int i = 0; i < self->ndit * self->nchannel; i++) {
    //     self->flag_cov[i] = gsl_matrix_int_alloc(self->nwave, self->nwave);
    //     for (int j = 0; j < self->nwave; j++) {
    //         if (gsl_matrix_int_get(self->flag, i, j)) {
    //             for  (int k = 0; k < self->nwave; k++) {
    //                 gsl_matrix_int_set(self->flag_cov[i], j, k, 1);
    //                 gsl_matrix_int_set(self->flag_cov[i], k, j, 1);
    //             }
    //         }
    //     }
    // }

    /* Compute VIS_REF (visData_phasedFT per manual §10.26) */
    self->vis_ref = gsl_matrix_complex_alloc(self->ndit * self->nchannel, self->nwave);
    for (int i = 0; i < self->ndit * self->nchannel; i++) {
        for (int j = 0; j < self->nwave; j++) {
            gsl_complex vis_ref_val = gsl_matrix_complex_get(self->visdata, i, j);
            gsl_complex pha = gsl_complex_polar(1.0, gsl_matrix_get(self->phase_ref, i, j));
            vis_ref_val = gsl_complex_mul(vis_ref_val, pha);
            gsl_matrix_complex_set(self->vis_ref, i, j, vis_ref_val);
        }
    }

    /* Construct covariance and psuedo-covariance matrices for VIS_REF */
    /* These are diagonal, so we only store a vector of the diagonal terms */
    self->vis_ref_cov = cpl_malloc(self->ndit * self->nchannel * sizeof(gsl_vector*));
    self->vis_ref_pcov = cpl_malloc(self->ndit * self->nchannel * sizeof(gsl_vector_complex*));
    gsl_vector *vreal = gsl_vector_alloc(self->nwave);
    gsl_vector *vimag = gsl_vector_alloc(self->nwave);
    for (int i = 0; i < self->ndit * self->nchannel; i++) {
        // Diagonal of |err|^2 for this dit and channel
        gsl_vector_complex_const_view err_vals = gsl_matrix_complex_const_row(self->viserr, i);
        gsl_vector_const_view real_vals = gsl_vector_complex_const_real(&err_vals.vector);
        gsl_vector_const_view imag_vals = gsl_vector_complex_const_imag(&err_vals.vector);
        // Real part squared
        gsl_vector_memcpy(vreal, &real_vals.vector);
        gsl_vector_mul(vreal, vreal);
        // Imag part squared
        gsl_vector_memcpy(vimag, &imag_vals.vector);
        gsl_vector_mul(vimag, vimag);

        // Form Cov
        gsl_vector *cov_diag = gsl_vector_alloc(self->nwave);
        gsl_vector_memcpy(cov_diag, vreal);
        gsl_vector_add(cov_diag, vimag);
        self->vis_ref_cov[i] = cov_diag;

        // Form PCov
        gsl_vector_complex *pcov_diag = gsl_vector_complex_calloc(self->nwave);
        gsl_vector_view pcov_diag_real = gsl_vector_complex_real(pcov_diag);
        gsl_vector_memcpy(&pcov_diag_real.vector, vreal);
        gsl_vector_sub(&pcov_diag_real.vector, vimag);

        for (int j = 0; j < self->nwave; j++) {
            // pha = pcov_diag * exp(2i*phase_ref)
            gsl_complex pcov_val = gsl_complex_polar(1.0, 2 * gsl_matrix_get(self->phase_ref, i, j));
            pcov_val = gsl_complex_mul(gsl_vector_complex_get(pcov_diag, j), pcov_val);
            gsl_vector_complex_set(pcov_diag, j, pcov_val);
        }
        self->vis_ref_pcov[i] = pcov_diag;
    }
    FREE(gsl_vector_free, vreal);
    FREE(gsl_vector_free, vimag);

    /* Normalise visibilities by the DIT */
    gravi_astrometry_scale_visibilities(self, 1.0 / self->dit);

    /* Compute metrology and dispersion phase corrections */
    gsl_matrix *wave_opd_disp = gsl_matrix_alloc(self->ndit * self->nchannel, self->nwave);
    gsl_matrix_memcpy(wave_opd_disp, self->opd_disp);
    
    gsl_vector *wavenumber = gsl_vector_alloc(self->nwave);
    gsl_vector_set_all(wavenumber, TWOPI);
    gsl_vector_div(wavenumber, self->wave);

    /* 2*pi/wave * opdDisp */
    for (int i = 0; i < self->ndit * self->nchannel; i++) {
        gsl_vector_view wave_opd_dist_row = gsl_matrix_row(wave_opd_disp, i);
        gsl_vector_mul(&wave_opd_dist_row.vector, wavenumber);
    }

    /* Total correction should be subtracted from phase, so negate first */
    gsl_matrix_add(wave_opd_disp, self->phase_met_telfc);
    gsl_matrix_scale(wave_opd_disp, -1);
    gravi_astrometry_add_phase(self, wave_opd_disp);

    FREE(gsl_vector_free, wavenumber);
    FREE(gsl_matrix_free, wave_opd_disp);

    return self;
}

cpl_error_code gravi_astrometry_dump(astro_data *self, FILE *handle)
{
    cpl_ensure_code(self, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code(handle, CPL_ERROR_NULL_INPUT);

    fprintf(handle, "%s\n", self->filename);
    fprintf(handle, "%s\n", self->insname);

    fprintf(handle, "%d %d %d\n", self->ndit, self->nchannel, self->nwave);
    
    gsl_matrix_fprintf(handle, self->ucoord, "%f");
    gsl_matrix_fprintf(handle, self->vcoord, "%f");
    gsl_matrix_int_fprintf(handle, self->flag, "%d");
    gsl_matrix_complex_fprintf(handle, self->visdata, "%f");
    gsl_matrix_complex_fprintf(handle, self->vis_ref, "%f");

    return CPL_ERROR_NONE;
}

void gravi_astrometry_delete(astro_data *self)
{
    FREE(cpl_free, self->filename);
    FREE(cpl_free, self->insname);
    FREE(gsl_vector_free, self->wave);
    FREE(gsl_vector_free, self->u);
    FREE(gsl_vector_free, self->v);

    FREE(gsl_matrix_complex_free, self->visdata);
    FREE(gsl_matrix_complex_free, self->visdata_ft);
    FREE(gsl_matrix_complex_free, self->viserr);
    FREE(gsl_matrix_free, self->phase_ref);
    FREE(gsl_matrix_free, self->phase_met_telfc);
    FREE(gsl_matrix_free, self->opd_disp);
    FREE(gsl_matrix_complex_free, self->vis_ref);
    FREE(gsl_matrix_int_free, self->flag);

    // FREELOOP(gsl_matrix_int_free, self->flag_cov, self->ndit * self->nchannel);
    FREELOOP(gsl_vector_free, self->vis_ref_cov, self->ndit * self->nchannel);
    FREELOOP(gsl_vector_complex_free, self->vis_ref_pcov, self->ndit * self->nchannel);

    // FREE(gsl_vector_free, self->swap_astrometry_guess);
    // FREE(gsl_vector_free, self->swap_astrometry_fit);
    FREE(gsl_matrix_free, self->amp_ref_astro);
    FREE(gsl_matrix_complex_free, self->phase_ref_astro);

    cpl_free(self);
}

double gravi_astrometry_get_mean_ftflux(astro_data *self)
{
    cpl_ensure(self, CPL_ERROR_NULL_INPUT, -1);

    double abs_val = 0;
    for (int i = 0; i < self->ndit * self->nchannel; i++) {
        for (int j = 0; j < self->nwave_ft; j++) {
            abs_val += gsl_complex_abs(gsl_matrix_complex_get(self->visdata_ft, i, j));
        }
    }
    return abs_val / (self->ndit * self->nchannel * self->nwave_ft);
}

/**
 * @brief Given observation and list of reference observations, find closest in time from the list.
 * 
 * @param self baseline astro_data.
 * @param others list of reference astro_data.
 * @param n_other length of list.
 * @param mode if FIND_MODE_BEFORE, require closest that precedes. If FIND_MODE_AFTER, require closest that succeeds.
 * Else, return closest regardless of sign of time difference.
 * 
 * @return index into passed list that is closest, or -1 if invalid.
*/
cpl_size gravi_astrometry_find_closest_mjd(astro_data *self, astro_data **others, cpl_size n_other, find_closest_mode mode)
{
    cpl_ensure(self, CPL_ERROR_NULL_INPUT, -1);
    cpl_ensure(others, CPL_ERROR_NULL_INPUT, -1);

    cpl_boolean any_before = FALSE, any_after = FALSE;
    for (int i = 0; i < n_other; i++) {
        double delta = self->mjd - others[i]->mjd;
        if (delta >= 0)
            any_before = CPL_TRUE;
        if (delta <= 0)
            any_after = CPL_TRUE;
    }

    cpl_size imin = -1;
    double min_delta = INFINITY;

    for (int i = 0; i < n_other; i++) {
        double delta = self->mjd - others[i]->mjd;
        switch (mode) {
case FIND_MODE_BEFORE:
            if (any_before && delta >= 0 && delta < min_delta) {
                imin = i;
                min_delta = delta;
            }
            break;
case FIND_MODE_AFTER:
            if (any_after && delta <= 0 && -delta < min_delta) {
                imin = i;
                min_delta = -delta;
            }
            break;
default:
            if (delta <= 0 && fabs(delta) < min_delta) {
                imin = i;
                min_delta = delta;
            }
            break;
        }
    }
    return imin;
}

/**
 * @brief Compute the final astrometric phase reference.
 * 
 * @param self astro_data to compute ref for.
 * @param phase_refs list of astro_data to compute ref from.
 * @param nphase length of list.
 * @param swaps list of astro_data to compute swap ref from.
 * @param nswap length of list.
 * @param parlist recipe parameters, esp. "calib_strategy".
 */
cpl_error_code gravi_astrometry_create_phase_reference(astro_data *self, astro_data **phase_refs, cpl_size nphase, astro_data **swaps, cpl_size nswap, cpl_parameterlist *parlist)
{
    cpl_ensure_code(self, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code(phase_refs, CPL_ERROR_NULL_INPUT);
    if (nswap > 0)
        cpl_ensure_code(swaps, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code(parlist, CPL_ERROR_NULL_INPUT);

    gsl_matrix *amp_ref = NULL, *amp_ref_tmp = NULL;
    gsl_matrix_complex *vis_ref = NULL, *vis_ref_tmp = NULL;

    /* Parameters */
    char *calib_strategy = cpl_strdup(cpl_parameter_get_string(
        cpl_parameterlist_find(parlist, "gravity.astrometry.calib-strategy")));
    /* just in case: tolerate mixed/lower case arguments */
    for(int i = 0; calib_strategy[i]; i++)
        calib_strategy[i] = toupper(calib_strategy[i]);
    
    CPLCHECK_INT("Could not get parameters");

    cpl_msg_debug(cpl_func, "Calibrating amplitude reference with '%s' strategy", calib_strategy);

    if (!strcmp(calib_strategy, "NONE")) {
        amp_ref = gsl_matrix_alloc(self->nchannel, self->nwave);
        vis_ref = gsl_matrix_complex_alloc(self->nchannel, self->nwave);
        gsl_matrix_set_all(amp_ref, 1.0);
        gsl_matrix_complex_set_all(vis_ref, GSL_COMPLEX_ONE);
    } else if (!strcmp(calib_strategy, "SELF")) {
        // abs first, then mean over dits
        amp_ref_tmp = gsl_matrix_alloc(self->ndit * self->nchannel, self->nwave);
        for (int i = 0; i < self->ndit * self->nchannel; i++) {
            for (int j = 0; j < self->nwave; j++) {
                gsl_complex vis_ref_val = gsl_matrix_complex_get(self->vis_ref, i, j);
                gsl_matrix_set(amp_ref_tmp, i, j, gsl_complex_abs(vis_ref_val));
            }
        }
        amp_ref = average_matrix_over_dits(amp_ref_tmp, self->ndit, self->nchannel, self->nwave, self->flag);
        FREE(gsl_matrix_free, amp_ref_tmp);

        vis_ref = gsl_matrix_complex_alloc(self->nchannel, self->nwave);
        for (int i = 0; i < self->nchannel; i++) {
            for (int j = 0; j < self->nwave; j++) {
                double amp_ref_val = gsl_matrix_get(amp_ref, i, j);
                gsl_matrix_complex_set(vis_ref, i, j, gsl_complex_rect(amp_ref_val, 0));
            }
        }
    } else {
        int nphase_used;
        gsl_vector_int *phase_indices;

        if (!strcmp(calib_strategy, "ALL")) {
            nphase_used = nphase;
            phase_indices = gsl_vector_int_alloc(nphase_used);
            for (int n = 0; n < nphase_used; n++)
                gsl_vector_int_set(phase_indices, n, n);
        } else if (!strcmp(calib_strategy, "NEAREST")) {
            nphase_used = 2;
            phase_indices = gsl_vector_int_alloc(nphase_used);
            cpl_size idx_before = gravi_astrometry_find_closest_mjd(self, phase_refs, nphase, FIND_MODE_BEFORE);
            cpl_size idx_after = gravi_astrometry_find_closest_mjd(self, phase_refs, nphase, FIND_MODE_AFTER);
            // cpl_msg_debug(cpl_func, "File %s", self->filename);
            // cpl_msg_debug(cpl_func, "Preceding on-star file %s", phase_refs[idx_before]->filename);
            // cpl_msg_debug(cpl_func, "Succeeding on-star file %s", phase_refs[idx_after]->filename);
            if (idx_before == -1 || idx_after == -1) {
                cpl_error_set_message(cpl_func, CPL_ERROR_ACCESS_OUT_OF_RANGE, "Cannot find time-bracketing phase references");
                return CPL_ERROR_ACCESS_OUT_OF_RANGE;
            }
            gsl_vector_int_set(phase_indices, 0, idx_before);
            gsl_vector_int_set(phase_indices, 1, idx_after);
        } else {
            cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT, "Unknown calibration strategy %s", calib_strategy);
            return CPL_ERROR_ILLEGAL_INPUT;
        }

        amp_ref = gsl_matrix_calloc(self->nchannel, self->nwave);
        vis_ref = gsl_matrix_complex_calloc(self->nchannel, self->nwave);
        amp_ref_tmp = gsl_matrix_alloc(self->nchannel, self->nwave);

        for (int n = 0; n < nphase_used; n++) {
            cpl_size idx = gsl_vector_int_get(phase_indices, n);
            astro_data *soi = phase_refs[idx];
            cpl_msg_debug(cpl_func, "Take idx %lld, filename %s", idx, soi->filename);

            // Average over DITs to get object of shape (n_channel, n_wave)
            vis_ref_tmp = average_matrix_complex_over_dits(soi->vis_ref, soi->ndit, soi->nchannel, soi->nwave, soi->flag);

            for (int i = 0; i < soi->nchannel; i++) {
                for (int j = 0; j < soi->nwave; j++) {
                    gsl_complex vis_ref_val = gsl_matrix_complex_get(vis_ref_tmp, i, j);
                    gsl_matrix_set(amp_ref_tmp, i, j, gsl_complex_abs(vis_ref_val));
                }
            }

            gsl_matrix_add(amp_ref, amp_ref_tmp);
            gsl_matrix_complex_add(vis_ref, vis_ref_tmp);
            FREE(gsl_matrix_complex_free, vis_ref_tmp);
        }

        gsl_matrix_scale(amp_ref, 1.0 / nphase_used);
        gsl_matrix_complex_scale(vis_ref, gsl_complex_rect(1.0 / nphase_used, 0.0));

        for (int i = 0; i < self->nchannel; i++) {
            for (int j = 0; j < self->nwave; j++) {
                double theta = gsl_complex_arg(gsl_matrix_complex_get(vis_ref, i, j));
                double r = gsl_matrix_get(amp_ref, i, j);
                gsl_matrix_complex_set(vis_ref, i, j, gsl_complex_polar(r, theta));
            }
        }
        FREE(gsl_matrix_free, amp_ref_tmp);
        FREE(gsl_vector_int_free, phase_indices);
    }
    
    self->amp_ref_astro = amp_ref;
    self->phase_ref_astro = vis_ref;

    /* extract phase reference from swaps */
    if (nswap > 0) {
        cpl_msg_debug(cpl_func, "Extract phase reference from swaps");
        /* Should already exist, but will be modified and/or overwritten */
        amp_ref = self->amp_ref_astro;
        vis_ref = self->phase_ref_astro;

        /* Take the amplitude reference from the swaps */
        if (!strcmp(calib_strategy, "SWAP")) {
            cpl_size n_swapon = 0;

            gsl_matrix_set_zero(amp_ref);
            gsl_matrix_complex_set_zero(vis_ref);

            for (int n = 0; n < nswap; n++) {
                astro_data *soi = swaps[n];

                if (soi->swap) {
                    n_swapon++;
                    // abs first, then mean over dits
                    amp_ref_tmp = gsl_matrix_alloc(soi->ndit * soi->nchannel, soi->nwave);
                    for (int i = 0; i < soi->ndit * soi->nchannel; i++) {
                        for (int j = 0; j < soi->nwave; j++) {
                            gsl_complex vis_ref_val = gsl_matrix_complex_get(soi->vis_ref, i, j);
                            gsl_matrix_set(amp_ref_tmp, i, j, gsl_complex_abs(vis_ref_val));
                        }
                    }
                    gsl_matrix *amp_ref_avg = average_matrix_over_dits(amp_ref_tmp, soi->ndit, soi->nchannel, soi->nwave, soi->flag);
                    gsl_matrix_add(amp_ref, amp_ref_avg);

                    FREE(gsl_matrix_free, amp_ref_tmp);
                    FREE(gsl_matrix_free, amp_ref_avg);
                }
            }

            gsl_matrix_scale(amp_ref, 1.0 / n_swapon);
            for (int i = 0; i < self->nchannel; i++) {
                for (int j = 0; j < self->nwave; j++) {
                    /* Phase ref will be added later. Divide by 2 is because flux is 2x higher in swap */
                    double amp_ref_val = gsl_matrix_get(amp_ref, i, j);
                    gsl_matrix_complex_set(vis_ref, i, j, gsl_complex_rect(0.5 * amp_ref_val, 0));
                }
            }
        }

        /* Finally calculate and add the phase ref */
        gsl_matrix_complex *swap_ref = swaps[0]->phase_ref_astro;
        for (int i = 0; i < self->nchannel; i++) {
            for (int j = 0; j < self->nwave; j++) {
                gsl_complex phi = gsl_matrix_complex_get(swap_ref, i, j);
                double argphi = gsl_complex_arg(phi);

                gsl_complex vis_ref_val = gsl_matrix_complex_get(vis_ref, i, j);
                double vis_ref_abs = gsl_complex_abs(vis_ref_val);
                /* Factor 2 because of the beamsplitter for on-star observations*/
                gsl_matrix_complex_set(vis_ref, i, j, gsl_complex_polar(2 * vis_ref_abs, argphi));
            }
        }
    }

    cpl_free(calib_strategy);
    return CPL_ERROR_NONE;
}

/**
 * Return CPL table with the final astrometric phase reference.
*/
// cpl_table *gravi_astrometry_get_phase_reference(astro_data *self)
// {
//     cpl_table *table = cpl_table_new(self->nchannel);
//     cpl_table_new_column_array(table, "ASTRO_VISREF", CPL_TYPE_DOUBLE_COMPLEX, self->nwave);
//     cpl_table_new_column_array(table, "ASTRO_AMPREF", CPL_TYPE_DOUBLE, self->nwave);
//     cpl_table_new_column_array(table, "ASTRO_COV", CPL_TYPE_DOUBLE, self->ndit * self->nwave);
//     cpl_table_new_column_array(table, "ASTRO_PCOV", CPL_TYPE_DOUBLE_COMPLEX, self->ndit * self->nwave);

//     cpl_array *tmp_arr = cpl_array_new(self->nwave, CPL_TYPE_DOUBLE_COMPLEX);
//     for (int i = 0; i < self->nchannel; i++) {
//         for (int j = 0; j < self->nwave; j++) {
//             gsl_complex val = gsl_matrix_complex_get(self->phase_ref_astro, i, j);
//             cpl_array_set_complex(tmp_arr, j, GSL_REAL(val) + I * GSL_IMAG(val));
//         }
//         cpl_table_set_array(table, "ASTRO_VISREF", i, tmp_arr);
//     }
//     FREE(cpl_array_delete, tmp_arr);

//     tmp_arr = cpl_array_new(self->nwave, CPL_TYPE_DOUBLE);
//     for (int i = 0; i < self->nchannel; i++) {
//         for (int j = 0; j < self->nwave; j++) {
//             double val = gsl_matrix_get(self->amp_ref_astro, i, j);
//             cpl_array_set(tmp_arr, j, val);
//         }
//         cpl_table_set_array(table, "ASTRO_AMPREF", i, tmp_arr);
//     }
//     FREE(cpl_array_delete, tmp_arr);

//     // vis_ref_cov is array of size ndit * nchannel, of vectors of length nwave
//     tmp_arr = cpl_array_new(self->ndit * self->nwave, CPL_TYPE_DOUBLE);
//     for (int i = 0; i < self->nchannel; i++) {
//         for (int j = 0; j < self->ndit; j++) {
//             for (int k = 0; k < self->nwave; k++) {
//                 double val = gsl_vector_get(self->vis_ref_cov[i + j * self->nchannel], k);
//                 cpl_array_set(tmp_arr, k + j * self->nwave, val);
//             }
//         }
//         cpl_table_set_array(table, "ASTRO_COV", i, tmp_arr);
//     }
//     CPLCHECK_NUL("Could not add ASTRO_COV");
//     FREE(cpl_array_delete, tmp_arr);

//     // vis_ref_pcov is array of size ndit * nchannel, of complex vectors of length nwave
//     tmp_arr = cpl_array_new(self->ndit * self->nwave, CPL_TYPE_DOUBLE_COMPLEX);
//     for (int i = 0; i < self->nchannel; i++) {
//         for (int j = 0; j < self->ndit; j++) {
//             for (int k = 0; k < self->nwave; k++) {
//                 gsl_complex val = gsl_vector_complex_get(self->vis_ref_pcov[i + j * self->nchannel], k);
//                 cpl_array_set_complex(tmp_arr, k + j * self->nwave, GSL_REAL(val) + I * GSL_IMAG(val));
//             }
//         }
//         cpl_table_set_array(table, "ASTRO_PCOV", i, tmp_arr);
//     }
//     CPLCHECK_NUL("Could not add ASTRO_PCOV");
//     FREE(cpl_array_delete, tmp_arr);

//     return table;
// }

/**
 * Return CPL table with the final astrometric phase reference.
 * 
 * @param self astro_data to extract table from.
 * 
 * @note Row axis of table corresponds to wavelength.
 */
cpl_table *gravi_astrometry_get_phase_reference(astro_data *self)
{
    cpl_table *table = cpl_table_new(self->nwave);
    cpl_table_new_column_array(table, "ASTRO_VISREF", CPL_TYPE_DOUBLE_COMPLEX, self->nchannel);
    cpl_table_new_column_array(table, "ASTRO_AMPREF", CPL_TYPE_DOUBLE, self->nchannel);
    cpl_table_new_column_array(table, "ASTRO_COV", CPL_TYPE_DOUBLE, self->ndit * self->nchannel);
    cpl_table_new_column_array(table, "ASTRO_PCOV", CPL_TYPE_DOUBLE_COMPLEX, self->ndit * self->nchannel);

    cpl_array *tmp_arr = cpl_array_new(self->nchannel, CPL_TYPE_DOUBLE_COMPLEX);
    for (int j = 0; j < self->nwave; j++) {
        for (int i = 0; i < self->nchannel; i++) {
            gsl_complex val = gsl_matrix_complex_get(self->phase_ref_astro, i, j);
            cpl_array_set_complex(tmp_arr, i, GSL_REAL(val) + I * GSL_IMAG(val));
        }
        cpl_table_set_array(table, "ASTRO_VISREF", j, tmp_arr);
    }
    CPLCHECK_NUL("Could not add ASTRO_VISREF");
    FREE(cpl_array_delete, tmp_arr);

    tmp_arr = cpl_array_new(self->nchannel, CPL_TYPE_DOUBLE);
    for (int j = 0; j < self->nwave; j++) {
        for (int i = 0; i < self->nchannel; i++) {
            double val = gsl_matrix_get(self->amp_ref_astro, i, j);
            cpl_array_set(tmp_arr, i, val);
        }
        cpl_table_set_array(table, "ASTRO_AMPREF", j, tmp_arr);
    }
    CPLCHECK_NUL("Could not add ASTRO_AMPREF");
    FREE(cpl_array_delete, tmp_arr);

    // vis_ref_cov is array of size ndit * nchannel, of vectors of length nwave
    tmp_arr = cpl_array_new(self->ndit * self->nchannel, CPL_TYPE_DOUBLE);
    for (int j = 0; j < self->nwave; j++) {
        for (int i = 0; i < self->ndit * self->nchannel; i++) {
            double val = gsl_vector_get(self->vis_ref_cov[i], j);
            cpl_array_set(tmp_arr, i, val);
        }
        cpl_table_set_array(table, "ASTRO_COV", j, tmp_arr);
    }
    CPLCHECK_NUL("Could not add ASTRO_COV");
    FREE(cpl_array_delete, tmp_arr);

    // vis_ref_pcov is array of size ndit * nchannel, of vectors of length nwave
    tmp_arr = cpl_array_new(self->ndit * self->nchannel, CPL_TYPE_DOUBLE_COMPLEX);
    for (int j = 0; j < self->nwave; j++) {
        for (int i = 0; i < self->ndit * self->nchannel; i++) {
            gsl_complex val = gsl_vector_complex_get(self->vis_ref_pcov[i], j);
            cpl_array_set_complex(tmp_arr, i, GSL_REAL(val) + I * GSL_IMAG(val));
        }
        cpl_table_set_array(table, "ASTRO_PCOV", j, tmp_arr);
    }
    CPLCHECK_NUL("Could not add ASTRO_PCOV");
    FREE(cpl_array_delete, tmp_arr);

    return table;
}

static gsl_vector *average_vector_over_dits(gsl_vector *v, int ndit, int nchannel)
{
    gsl_vector *vavg = gsl_vector_alloc(nchannel);

    for (int i = 0; i < nchannel; i++) {
        double avg = gsl_stats_mean(v->data + i, nchannel, ndit);
        gsl_vector_set(vavg, i, avg);
    }

    return vavg;
}

static gsl_matrix *average_matrix_over_dits(gsl_matrix *m, int ndit, int nchannel, int nwave, gsl_matrix_int *flag)
{
    gsl_matrix *mavg = gsl_matrix_alloc(nchannel, nwave);

    for (int i = 0; i < nchannel; i++) {
        for (int j = 0; j < nwave; j++) {
            cpl_size start = i * nwave + j;
            gsl_vector_const_view dit_view = gsl_vector_const_view_array_with_stride(
                m->data + start, nwave * nchannel, ndit);
            gsl_vector_int_const_view flag_view = flag ?
                gsl_vector_int_const_view_array_with_stride(flag->data + start, nwave * nchannel, ndit) :
                (gsl_vector_int_const_view){};
            
            double avg = 0.0;
            cpl_size Nvalid = 0;
            for (int k = 0; k < ndit; k++) {
                if (!flag || gsl_vector_int_get(&flag_view.vector, k) == 0) {
                    avg += gsl_vector_get(&dit_view.vector, k);
                    Nvalid++;
                }
            }
            if (Nvalid > 0) {
                avg /= Nvalid;
                gsl_matrix_set(mavg, i, j, avg);
            } else {
                gsl_matrix_set(mavg, i, j, 0.0);
            }
        }
    }

    return mavg;
}

static gsl_matrix_complex *average_matrix_complex_over_dits(gsl_matrix_complex *m, int ndit, int nchannel, int nwave, gsl_matrix_int *flag)
{
    gsl_matrix_complex *mavg = gsl_matrix_complex_alloc(nchannel, nwave);

    for (int i = 0; i < nchannel; i++) {
        for (int j = 0; j < nwave; j++) {
            cpl_size start = i * nwave + j;
            gsl_vector_complex_const_view dit_view = gsl_vector_complex_const_view_array_with_stride(
                m->data + 2 * start, nwave * nchannel, ndit); /* factor 2 for complex */
            gsl_vector_int_const_view flag_view = flag ?
                gsl_vector_int_const_view_array_with_stride(flag->data + start, nwave * nchannel, ndit) :
                (gsl_vector_int_const_view){};
            
            gsl_complex avg = GSL_COMPLEX_ZERO;
            cpl_size Nvalid = 0;
            for (int k = 0; k < ndit; k++) {
                if (!flag || gsl_vector_int_get(&flag_view.vector, k) == 0) {
                    avg = gsl_complex_add(avg, gsl_vector_complex_get(&dit_view.vector, k));
                    Nvalid++;
                }
            }
            if (Nvalid > 0) {
                avg = gsl_complex_div_real(avg, Nvalid);
                gsl_matrix_complex_set(mavg, i, j, avg);
            } else {
                gsl_matrix_complex_set(mavg, i, j, GSL_COMPLEX_ZERO);
            }
        }
    }

    return mavg;
}

/**
 * @brief Average input astro_data over DITs and return new flattened astro_data.
*/
static astro_data *gravi_astrometry_average_over_dits(astro_data *self)
{
    cpl_ensure(self != NULL, CPL_ERROR_NULL_INPUT, NULL);
    cpl_msg_debug(cpl_func, "Averaging file %s", self->filename);

    astro_data *flat = cpl_malloc(sizeof(astro_data));

    flat->filename = cpl_strdup(self->filename);
    flat->insname = cpl_strdup(self->insname);

    flat->nchannel = self->nchannel;
    flat->ndit = 1;
    flat->dit = self->dit;
    flat->sobj_x = self->sobj_x;
    flat->sobj_y = self->sobj_y;
    flat->mjd = self->mjd;
    flat->swap = self->swap;
    flat->nwave = self->nwave;
    flat->nwave_ft = self->nwave_ft;

    // flat->swap_astrometry_guess = NULL;
    // flat->swap_astrometry_fit = NULL;

    flat->wave = gsl_vector_alloc(flat->nwave);
    gsl_vector_memcpy(flat->wave, self->wave);

    /* Calculate mean in SC frame */
    gravi_astrometry_recentre_phase(self, self->sobj_x, self->sobj_y);

    flat->u = average_vector_over_dits(self->u, self->ndit, self->nchannel);
    flat->v = average_vector_over_dits(self->u, self->ndit, self->nchannel);
    flat->ucoord = average_matrix_over_dits(self->ucoord, self->ndit, self->nchannel, self->nwave, self->flag);
    flat->vcoord = average_matrix_over_dits(self->vcoord, self->ndit, self->nchannel, self->nwave, self->flag);

    flat->visdata = average_matrix_complex_over_dits(self->visdata, self->ndit, self->nchannel, self->nwave, self->flag);
    flat->vis_ref = average_matrix_complex_over_dits(self->vis_ref, self->ndit, self->nchannel, self->nwave, self->flag);

    /* cov, pcov are array of shape (ndit * nchannel,) of vectors with shape (nwave,)
       sum over dits, divide by ndit**2 */
    flat->vis_ref_cov = cpl_malloc(flat->nchannel * sizeof(gsl_vector *));
    flat->vis_ref_pcov = cpl_malloc(flat->nchannel * sizeof(gsl_vector_complex *));

    for (int i = 0; i < self->nchannel; i++) {
        flat->vis_ref_cov[i] = gsl_vector_calloc(flat->nwave);
        flat->vis_ref_pcov[i] = gsl_vector_complex_calloc(flat->nwave);

        for (int j = 0; j < self->ndit; j++) {
            cpl_size idx = j * self->nchannel + i;
            gsl_vector_add(flat->vis_ref_cov[i], self->vis_ref_cov[idx]);
            gsl_vector_complex_add(flat->vis_ref_pcov[i], self->vis_ref_pcov[idx]);
        }
        gsl_vector_scale(flat->vis_ref_cov[i], 1.0 / (self->ndit * self->ndit));
        gsl_vector_complex_scale(flat->vis_ref_pcov[i], gsl_complex_rect(1.0 / (self->ndit * self->ndit), 0));
    }

    flat->nflag = 0;
    flat->flag = gsl_matrix_int_calloc(flat->nchannel, flat->nwave);
    for (int i = 0; i < self->nchannel; i++) {
        for (int j = 0; j < self->nwave; j++) {
            cpl_size start = i * self->nwave + j;
            gsl_vector_int_const_view flag_view = gsl_vector_int_const_view_array_with_stride(
                self->flag->data + start, self->nwave * self->nchannel, self->ndit);
            
            if (_gsl_vector_int_sum(&flag_view.vector) == self->ndit) {
                gsl_matrix_int_set(flat->flag, i, j, 1);
                flat->nflag++;
            }
        }
    }

    /* Shift flattened back to FT frame */
    gravi_astrometry_recentre_phase(flat, -flat->sobj_x, -flat->sobj_y);

    /* Shift original back to FT frame */
    gravi_astrometry_recentre_phase(self, -self->sobj_x, -self->sobj_y);

    return flat;
}

static gsl_matrix_complex *gravi_astrometry_calculate_visref_swap(astro_data **data, cpl_size ndata, double ra, double dec) {
    gsl_matrix_complex *visref = gsl_matrix_complex_calloc(data[0]->nchannel, data[0]->nwave);
    gsl_matrix_complex *visref_ditsum = gsl_matrix_complex_alloc(data[0]->nchannel, data[0]->nwave);
    
    gsl_matrix_int *ngood = gsl_matrix_int_calloc(data[0]->nchannel, data[0]->nwave);
    gsl_matrix_int *ngood_ditsum = gsl_matrix_int_alloc(data[0]->nchannel, data[0]->nwave);

    for (int n = 0; n < ndata; n++) {
        astro_data *self = data[n];

        for (int j = 0; j < self->nchannel; j++) {
            for (int k = 0; k < self->nwave; k++) {
                int ngood_jk = 0;
                gsl_complex visref_jk = GSL_COMPLEX_ZERO;
                
                for (int i = 0; i < self->ndit; i++) {
                    cpl_size idx = i * self->nchannel + j;
                    double uv = gsl_matrix_get(self->ucoord, idx, k) * ra + gsl_matrix_get(self->vcoord, idx, k) * dec;
                    double phase = uv * TWOPI * MAS_TO_RAD;
        
                    gsl_complex phi = gsl_complex_polar(1.0, phase);
                    if (self->swap)
                        phi = gsl_complex_conjugate(phi);
                    gsl_complex visref_ijk = gsl_complex_mul(phi, gsl_matrix_complex_get(self->vis_ref, idx, k));

                    int flag = gsl_matrix_int_get(self->flag, idx, k);
                    if (flag == 0) {
                        ngood_jk++;
                        visref_jk = gsl_complex_add(visref_jk, visref_ijk);
                    }
                }                
                gsl_matrix_complex_set(visref_ditsum, j, k, visref_jk);
                gsl_matrix_int_set(ngood_ditsum, j, k, ngood_jk);
            }
        }
        gsl_matrix_complex_add(visref, visref_ditsum);
        gsl_matrix_int_add(ngood, ngood_ditsum);
    }

    for (int j = 0; j < data[0]->nchannel; j++) {
        for (int k = 0; k < data[0]->nwave; k++) {
            int ngood_jk = gsl_matrix_int_get(ngood, j, k);
            if (ngood_jk > 0) {
                gsl_complex visref_jk = gsl_matrix_complex_get(visref, j, k);
                visref_jk = gsl_complex_div_real(visref_jk, ngood_jk);
                gsl_matrix_complex_set(visref, j, k, visref_jk);
            }
        }
    }

    FREE(gsl_matrix_complex_free, visref_ditsum);
    return visref;
}

/**
 * @brief Evaluate chi^2 statistic.
 * 
 * @param x Vector of model parameters [ra, dec].
 * @param params Fit parameter values, of actual type @c gravi_astrometry_model_params.
 */
static double gravi_astrometry_calculate_chi2(const gsl_vector *X, void *params)
{
    gravi_astrometry_model_params *p = params;
    astro_data** s1 = p->group1;
    astro_data** s2 = p->group2;

    cpl_size nchannel = s1[0]->nchannel;
    cpl_size nwave = s1[0]->nwave;

    double ra  = gsl_vector_get(X, 0);
    double dec = gsl_vector_get(X, 1);

    gsl_matrix_complex *visref_s1 = gravi_astrometry_calculate_visref_swap(s1, p->n1, ra, dec);
    gsl_matrix_complex *visref_s2 = gravi_astrometry_calculate_visref_swap(s2, p->n2, ra, dec);

    double chi2 = 0.0;
    for (int i = 0; i < nchannel; i++) {
        double chi2_channel = 0;
        for (int k = 0; k < nwave; k++) {
            gsl_complex visref_swap = gsl_complex_mul(
                gsl_matrix_complex_get(visref_s1, i, k),
                gsl_complex_conjugate(gsl_matrix_complex_get(visref_s2, i, k))
            );
            visref_swap = gsl_complex_sqrt(visref_swap);
            chi2_channel += GSL_IMAG(visref_swap) * GSL_IMAG(visref_swap);
        }
        chi2 += chi2_channel;
    }

    FREE(gsl_matrix_complex_free, visref_s1);
    FREE(gsl_matrix_complex_free, visref_s2);
    return chi2;
}

static cpl_error_code gravi_astrometry_minimise_chi2_grid(
    gravi_astrometry_model_params *params, const gsl_vector *ra_grid, const gsl_vector *dec_grid, gsl_vector *ra_dec_out, gsl_matrix *chi2_map)
{
    cpl_ensure_code(ra_grid, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code(dec_grid, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code(ra_dec_out, CPL_ERROR_NULL_INPUT);

    cpl_size n_ra = ra_grid->size;
    cpl_size n_dec = dec_grid->size;
    double ra_dec[2] = {0};
    gsl_vector_const_view ra_dec_view = gsl_vector_const_view_array(ra_dec, 2);
    double best_ra, best_dec, best_chi2 = INFINITY;
    
    for (int i = 0; i < n_ra; i++) {
        ra_dec[0] = gsl_vector_get(ra_grid, i);
        cpl_msg_debug(cpl_func, "calculating chi2 for ra=%f", ra_dec[0]);
        for (int j = 0; j < n_dec; j++) {
            ra_dec[1] = gsl_vector_get(dec_grid, j);

            double chi2 = gravi_astrometry_calculate_chi2(
                &ra_dec_view.vector, params
            );

            if (chi2_map)
                gsl_matrix_set(chi2_map, i, j, chi2);

            if (chi2 < best_chi2) {
                best_ra = ra_dec[0];
                best_dec = ra_dec[1];
                best_chi2 = chi2;
            }
        }
    }

    cpl_msg_debug(cpl_func, "Best chi2 is %f at [%f, %f]",
        best_chi2, best_ra, best_dec);

    gsl_vector_set(ra_dec_out, 0, best_ra);
    gsl_vector_set(ra_dec_out, 1, best_dec);
    return CPL_ERROR_NONE;
}

static cpl_error_code gravi_astrometry_minimise_chi2_descent(gravi_astrometry_model_params *params, double ra_guess, double dec_guess, gsl_vector *Xsolve)
{
    /* The minimisation has 2 parameters: ra, dec */
    const int dim = 2;
    
    const int MAX_ITERATIONS = 1000;
    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2rand;
    gsl_multimin_fminimizer *mini = gsl_multimin_fminimizer_alloc(T, dim);
    int status;
    double size;

    gsl_multimin_function func = {
        .f = &gravi_astrometry_calculate_chi2,
        .n = dim,
        .params = params
    };

    double initial_guess[] = {ra_guess, dec_guess};
    gsl_vector_view guess_view = gsl_vector_view_array(initial_guess, dim);
    double step_size[] = {0.1, 0.1};
    gsl_vector_view step_view = gsl_vector_view_array(step_size, dim);
    gsl_multimin_fminimizer_set(mini, &func, &guess_view.vector, &step_view.vector);

    int iterations = 0;
    do {
        iterations++;
        status = gsl_multimin_fminimizer_iterate(mini);
    
        /* Error taking optimisation step */
        if (status)
            break;
        
        size = gsl_multimin_fminimizer_size(mini);
        status = gsl_multimin_test_size(size, 1e-3);
    } while (status == GSL_CONTINUE && iterations < MAX_ITERATIONS);

    if (status == GSL_SUCCESS)
        gsl_vector_memcpy(Xsolve, mini->x);
    
    if (status == GSL_CONTINUE)
        cpl_msg_warning(cpl_func, "model-fitting did not converge");
    else if (status != GSL_SUCCESS)
        cpl_msg_error(cpl_func, "model-fitting failed with error code %d (%s)", status, gsl_strerror(status));
    else
        cpl_msg_debug(cpl_func, "converged with best-fit:\nra\tdec\tchi2\n%.3f\t%.3f\t%.4f", mini->x->data[0], mini->x->data[1], mini->fval);

    FREE(gsl_multimin_fminimizer_free, mini);
    return status ? CPL_ERROR_ILLEGAL_OUTPUT : CPL_ERROR_NONE;
}

cpl_error_code gravi_astrometry_reduce_swaps(astro_data **swap_data, cpl_size nswap, cpl_parameterlist *parlist)
{
    cpl_ensure_code(swap_data, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code(parlist, CPL_ERROR_NULL_INPUT);

    astro_data **group1 = NULL, **group2 = NULL;
    cpl_size n1 = 0, n2 = 0;

    gsl_vector *ra_grid = NULL, *dec_grid = NULL;

    gsl_vector *ra_dec_fit = NULL;
    gsl_matrix *chi2_map = NULL;

    /* Parameters */
    cpl_boolean use_swap_fiber_pos = cpl_parameter_get_bool(
        cpl_parameterlist_find(parlist, "gravity.astrometry.use-swap-fiber-pos"));

    cpl_boolean go_fast = cpl_parameter_get_bool(
        cpl_parameterlist_find(parlist, "gravity.astrometry.average-over-dits"));

    double ra_lim = cpl_parameter_get_double(
        cpl_parameterlist_find(parlist, "gravity.astrometry.ra-lim-swap"));
    int n_ra = cpl_parameter_get_int(
        cpl_parameterlist_find(parlist, "gravity.astrometry.nra-swap"));

    double dec_lim = cpl_parameter_get_double(
        cpl_parameterlist_find(parlist, "gravity.astrometry.ra-lim-swap"));
    int n_dec = cpl_parameter_get_int(
        cpl_parameterlist_find(parlist, "gravity.astrometry.ndec-swap"));

    double zoom = cpl_parameter_get_double(
        cpl_parameterlist_find(parlist, "gravity.astrometry.zoom-factor"));

    CPLCHECK_CLEAN("Could not get parameters");

    if (use_swap_fiber_pos) {
        cpl_msg_warning(cpl_func, "No astrometric solution is computed for swaps. Default to fiber position");
        for (int i = 0; i < nswap; i++) {
            astro_data *swap = swap_data[i];
            // swap->swap_astrometry_guess = gsl_vector_alloc(2);
            // gsl_vector_set(swap->swap_astrometry_guess, 0, swap->sobj_x);
            // gsl_vector_set(swap->swap_astrometry_guess, 1, swap->sobj_y);
            // swap->swap_astrometry_fit = gsl_vector_alloc(2);
            // gsl_vector_memcpy(swap->swap_astrometry_fit, swap->swap_astrometry_guess);
            swap->swap_astrometry_guess[0] = swap->swap_astrometry_fit[0] = swap->sobj_x;
            swap->swap_astrometry_guess[1] = swap->swap_astrometry_fit[1] = swap->sobj_y;
        }
        return CPL_ERROR_NONE;
    }

    /* Divide the swaps into on/off-axis groups */
    for (int i = 0; i < nswap; i++) {
        if (swap_data[i]->swap)
            n1++;
        else
            n2++;
    }

    if (n1 == 0)
        cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT, "No files with SWAP=YES");
    else if (n2 == 0)
        cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT, "No files with SWAP=NO");
    CPLCHECK_CLEAN("Could not partition swaps");

    group1 = cpl_malloc(n1 * sizeof(astro_data*));
    group2 = cpl_malloc(n2 * sizeof(astro_data*));
    cpl_size i1 = 0, i2 = 0;

    for (int i = 0; i < nswap; i++) {
        if (swap_data[i]->swap) {
            group1[i1++] = go_fast ? gravi_astrometry_average_over_dits(swap_data[i]) : swap_data[i];
        } else {
            group2[i2++] = go_fast ? gravi_astrometry_average_over_dits(swap_data[i]) : swap_data[i];
        }
    }

    gravi_astrometry_model_params params = {
            .group1 = group1, .group2 = group2, .n1 = n1, .n2 = n2
    };

    /* Determine RA/dec grid */
    int field; 
    if (strstr(swap_data[0]->insname, "U1234"))
        field = 60;
    else
        field = 240;

    if (ra_lim < 0)
        ra_lim = field / 2;
    if (dec_lim < 0)
        dec_lim = field / 2;

    double avg_sx = 0.0, avg_sy = 0.0;
    for (int i = 0; i < nswap; i++) {
        avg_sx += swap_data[i]->sobj_x * (swap_data[i]->swap ? -1 : 1);
        avg_sy += swap_data[i]->sobj_y * (swap_data[i]->swap ? -1 : 1);
    }
    avg_sx /= nswap;
    avg_sy /= nswap;
    cpl_msg_debug(cpl_func, "Centre: [%f, %f] mas", avg_sx, avg_sy);

    double ra_fit, dec_fit;
    double ra_min = avg_sx - ra_lim, ra_max = avg_sx + ra_lim;
    double dec_min = avg_sy - dec_lim, dec_max = avg_sy + dec_lim;
    double d_ra = (ra_max - ra_min) / (n_ra - 1);
    double d_dec = (dec_max - dec_min) / (n_dec - 1);

    ra_grid = gsl_vector_alloc(n_ra);
    dec_grid = gsl_vector_alloc(n_dec);
    for (int i = 0; i < n_ra; i++)
        gsl_vector_set(ra_grid, i, ra_min + i * d_ra);
    for (int j = 0; j < n_dec; j++)
        gsl_vector_set(dec_grid, j, dec_min + j * d_dec);

    ra_dec_fit = gsl_vector_alloc(2);
    chi2_map = gsl_matrix_alloc(n_ra, n_dec);

    if (zoom > 1.0) {
        cpl_msg_debug(cpl_func, "RA grid: [%f, %f] with %d points", ra_min, ra_max, n_ra);
        cpl_msg_debug(cpl_func, "dec grid: [%f, %f] with %d points", dec_min, dec_max, n_dec);

        gravi_astrometry_minimise_chi2_grid(&params, ra_grid, dec_grid, ra_dec_fit, chi2_map);
        ra_fit = gsl_vector_get(ra_dec_fit, 0);
        dec_fit = gsl_vector_get(ra_dec_fit, 1);
        cpl_msg_debug(cpl_func, "big chi2 map minimsed for ra=%f, dec=%f", ra_fit, dec_fit);

        // FILE *f = fopen("chi2_map_initial.dat", "w");
        // gsl_matrix_fprintf(f, chi2_map, "%g");
        // fclose(f);

        ra_lim = (ra_max - ra_min) / (4 * zoom);
        dec_lim = (dec_max - dec_min) / (4 * zoom);

        ra_min = ra_fit - ra_lim;
        ra_max = ra_fit + ra_lim;
        dec_min = dec_fit - dec_lim;
        dec_max = dec_fit + dec_lim;

        d_ra = (ra_max - ra_min) / (n_ra - 1);
        d_dec = (dec_max - dec_min) / (n_dec - 1);

        for (int i = 0; i < n_ra; i++)
            gsl_vector_set(ra_grid, i, ra_min + i * d_ra);
        for (int j = 0; j < n_dec; j++)
            gsl_vector_set(dec_grid, j, dec_min + j * d_dec);
    }

    cpl_msg_debug(cpl_func, "RA grid: [%f, %f] with %d points", ra_min, ra_max, n_ra);
    cpl_msg_debug(cpl_func, "dec grid: [%f, %f] with %d points", dec_min, dec_max, n_dec);
    
    gravi_astrometry_minimise_chi2_grid(&params, ra_grid, dec_grid, ra_dec_fit, chi2_map);
    ra_fit = gsl_vector_get(ra_dec_fit, 0);
    dec_fit = gsl_vector_get(ra_dec_fit, 1);
    cpl_msg_debug(cpl_func, "zoomed chi2 map minimised for ra=%f, dec=%f", ra_fit, dec_fit);

    // FILE *f = fopen("chi2_map.dat", "w");
    // gsl_matrix_fprintf(f, chi2_map, "%g");
    // fclose(f);

    for (int i = 0; i < nswap; i++) {
        // swap_data[i]->swap_astrometry_guess = gsl_vector_alloc(2);
        // gsl_vector_memcpy(swap_data[i]->swap_astrometry_guess, ra_dec_fit);
        swap_data[i]->swap_astrometry_guess[0] = ra_fit;
        swap_data[i]->swap_astrometry_guess[1] = dec_fit;
    }

    gravi_astrometry_minimise_chi2_descent(&params, ra_fit, dec_fit, ra_dec_fit);
    CPLCHECK_CLEAN("Could not minimise swap astrometry");

    ra_fit = gsl_vector_get(ra_dec_fit, 0);
    dec_fit = gsl_vector_get(ra_dec_fit, 1);
    cpl_msg_debug(cpl_func, "final astrometric solution after gradient descent is ra=%f, dec=%f", ra_fit, dec_fit);

    for (int i = 0; i < nswap; i++) {
        // swap_data[i]->swap_astrometry_fit = gsl_vector_alloc(2);
        // gsl_vector_memcpy(swap_data[i]->swap_astrometry_fit, ra_dec_fit);
        swap_data[i]->swap_astrometry_fit[0] = ra_fit;
        swap_data[i]->swap_astrometry_fit[1] = dec_fit;
    }

    /* Compute the swap phaseref, which is common to all the swap data */
    cpl_msg_debug(cpl_func, "Create the swap phase reference");
    gsl_matrix_complex *swap_ref = gravi_astrometry_create_swap_reference(swap_data, nswap);
    CPLCHECK_MSG("Could not extract swap phase reference");
    swap_data[0]->phase_ref_astro = swap_ref;
    for (int i = 1; i < nswap; i++)
        swap_data[i]->phase_ref_astro = gsl_matrix_complex_alloc_from_matrix(swap_ref, 0, 0, swap_ref->size1, swap_ref->size2);

cleanup:
    if (go_fast) {
        FREELOOP(gravi_astrometry_delete, group1, n1);
        FREELOOP(gravi_astrometry_delete, group2, n2);
    } else {
        FREE(cpl_free, group1);
        FREE(cpl_free, group2);
    }

    FREE(gsl_matrix_free, chi2_map);
    FREE(gsl_vector_free, ra_grid);
    FREE(gsl_vector_free, dec_grid);
    FREE(gsl_vector_free, ra_dec_fit);

    return CPL_ERROR_NONE;
}

/**
 * @brief Centre swap data on zero OPD and extract average phase.
 * 
 * @param swap_data list of astro_data for swaps.
 * @param nswap length of list.
 */
static gsl_matrix_complex *gravi_astrometry_create_swap_reference(astro_data **swap_data, cpl_size nswap)
{
    cpl_ensure(swap_data, CPL_ERROR_NULL_INPUT, NULL);

    int nchannel = swap_data[0]->nchannel;
    int nwave = swap_data[0]->nwave;

    gsl_matrix_complex *phase_ref_s1  = gsl_matrix_complex_calloc(nchannel, nwave);
    gsl_matrix_complex *phase_ref_s2  = gsl_matrix_complex_calloc(nchannel, nwave);
    gsl_matrix_complex *phase_ref_avg = gsl_matrix_complex_calloc(nchannel, nwave);

    for (int n = 0; n < nswap; n++) {
        astro_data *soi = swap_data[n];
        cpl_ensure(
            soi->nchannel == nchannel &&
            soi->nwave == nwave,
            CPL_ERROR_INCOMPATIBLE_INPUT,
            NULL
        );

        /* First shift swap visibilities to zero OPD using the fitted ra and dec */
        // double swap_ra = gsl_vector_get(soi->swap_astrometry_fit, 0) * (soi->swap ? -1 : 1);
        // double swap_dec = gsl_vector_get(soi->swap_astrometry_fit, 1) * (soi->swap ? -1 : 1);
        double swap_ra = soi->swap_astrometry_fit[0] * (soi->swap ? -1 : 1);
        double swap_dec =soi->swap_astrometry_fit[1] * (soi->swap ? -1 : 1);
        gravi_astrometry_recentre_phase(soi, swap_ra, swap_dec);
        
        /* Now take mean of vis and extract phase, for each position of the swap */
        gsl_matrix_complex *vis_ref_avg = average_matrix_complex_over_dits(soi->vis_ref, soi->ndit, soi->nchannel, soi->nwave, NULL);
        if (soi->swap)
            gsl_matrix_complex_add(phase_ref_s1, vis_ref_avg);
        else
            gsl_matrix_complex_add(phase_ref_s2, vis_ref_avg);
        FREE(gsl_matrix_complex_free, vis_ref_avg);
        
    }

    gsl_matrix_complex_memcpy(phase_ref_avg, phase_ref_s1);
    gsl_matrix_complex_add(phase_ref_avg, phase_ref_s2);
    gsl_matrix_complex_scale(phase_ref_avg, gsl_complex_rect(0.5, 0));

    FREE(gsl_matrix_complex_free, phase_ref_s1);
    FREE(gsl_matrix_complex_free, phase_ref_s2);

    return phase_ref_avg;
}
