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
 * @defgroup gravi_demodulate  Demodulation of metrology data
 *
 * This file contains code to remove the modulation signal from the metrology data.
 * The main function @c gravi_metrology_demodulate is called by the gravity_vis recipe
 * before proceeding with the standard metrology reduction.
 * The demodulation is accomplished by splitting the exposure into chunks of no more
 * than 100s duration (over which the change in modulation parameters due to sky motion
 * can be safely neglected), fitting a model for the modulation to each chunk, and
 * removing the modelled component.
 */

/*-----------------------------------------------------------------------------
                                   Includes
 -----------------------------------------------------------------------------*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cpl.h>
#include <math.h>
#include <time.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_statistics.h>

#include "gravi_demodulate.h"

#include "gravi_pfits.h"
#include "gravi_utils.h"

#define FT 0
#define SC 1

#define STEPS_PER_SECOND 500
#define MAX_SECONDS_PER_CHUNK 100

#define PI 3.14159265359
#define TWOPI 6.28318530718

/*-----------------------------------------------------------------------------
                              Private prototypes
 -----------------------------------------------------------------------------*/

/**
 * @brief Return column index in VOLT table for specific diode.
 * 
 * @param tele Telescope in {0..3}
 * @param diode Diode in {0..3}
 * @param side 0 (FT) or 1 (SC)
 * 
 * @return Column index for xvolts, add one for corresponding yvolts
*/
static cpl_size diode_column_index (int tele, int diode, int side)
{
    cpl_size i, r;
    i = (side == FT) ? 0 : 16;
    r = i + (4 * tele) + diode;
    return 2 * r;
}

const cpl_size ntel = 4;
const cpl_size ndiode = 4;
const cpl_size nside = 2;

/**
 * Diode zero offsets from Stefan Gillesen dated TODO
 * To be updated with new values dated TODO
 * Ordering is such that it matches the VOLT column in the metrology data
 * Used iff no DARK was provided
 **/
const double diode_zeros[] = {
    -0.0039217663072871976, -0.005375368693370768,  -0.0039996565553508424, -0.0025503151499257103,
    -0.004993987370698725,  -0.0025982105556817733, -0.0033492901476214598, -0.0024666122995294915,
    -0.004914056163432626,  -0.004605706229646769,  -0.0024913748841215717, -0.0025406836901383087,
    -0.0035045471990661765, -0.004864493604480354,  -0.0038971489167063875, -0.002645429318536674,
    -0.0031186035580824016, -0.0030236901965512925, -0.00455348514695007,   -0.0027973126198521715,
    -0.00329832479577174,   -0.004717880979213343,  -0.0038105012008050844, -0.0034806513678800706,
    -0.004843769955531035,  -0.0038383890934418244, -0.003413267406206906,  -0.0032447329693567305,
    -0.0019067813627941682, -0.003741090151381262,  -0.0031293194990934074, -0.002917132698761452,
    -0.003846563742841961,  -0.0016110031116351089, -0.004055586210967794,  -0.004665653836460701,
    -0.0035651245279117107, -0.003703624441150295,  -0.0042027507769680774, -0.003685527333212408,
    -0.0018717967368672794, -0.004309454414641543,  -0.004161530517996172,  -0.003603457935318008,
    -0.0031196452749965883, -0.003456638711821881,  -0.005237079636108828,  -0.0029119219882156274,
    -0.002883279151359765,  -0.0044843120806118035, -0.0033190561101935586, -0.004244083297006716,
    -0.005334829455674041,  -0.004158489195313494,  -0.00242824997808483,   -0.0018624941393818074,
    -0.0029737356896519826, -0.004348274099950943,  -0.002654502952395225,  -0.0030576540586404475,
    -0.0032698514683164475, -0.003638017106725761,  -0.0029787962274195018, -0.002629800558439812,
     0.0005056591585058203,  0.0006414248583206548,  0.0005120370824436637,  0.0006245251138490539,
     0.0002781845526483057, -0.00026739218921257777, 0.0005758898707200935,  0.0003888211264503715,
     0.00036289398311494447, 6.419787480760673e-05,  0.0008198190512220566,  0.0007714407655275396,
     0.0008010985218159154,  9.667686713004861e-05,  0.0008396378442216361,  0.0005619108431658763
};

static double model_x (const gsl_vector *X, double offset, int step)
{
    double fstep = (1.0 * step) / STEPS_PER_SECOND;
    double a = gsl_vector_get(X, 0);
    double b = gsl_vector_get(X, 1);
    double pha1 = gsl_vector_get(X, 2);
    double pha2 = gsl_vector_get(X, 3);

    return offset + a * sin(b * sin(fstep * TWOPI + pha1) + pha2);
}

static double model_y (const gsl_vector *X, double offset, int step)
{
    double fstep = (1.0 * step) / STEPS_PER_SECOND;
    double a = gsl_vector_get(X, 0);
    double b = gsl_vector_get(X, 1);
    double pha1 = gsl_vector_get(X, 2);
    double pha2 = gsl_vector_get(X, 3);

    return offset + a * cos(b * sin(fstep * TWOPI + pha1) + pha2);
}

static double modulation_model_chi2 (const gsl_vector *x, void *params);

typedef struct _demodulation_model_params_ {
    const gsl_vector *volts_x, *volts_y;
    double offset_x, offset_y;
} model_params;

static cpl_error_code fit_model_modulation (const gsl_vector *vx, const gsl_vector *vy, double zx, double zy, gsl_vector *Xsolve);

/*-----------------------------------------------------------------------------
                              Function code
 -----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
 * @brief Demodulate the metrology.
 *
 * @param[inout] met_data       input table of metrology data
 * @param[in] zero_subtracted   flag indicating metrology has been zero-subtracted using DARK
 * 
 * @return CPL_ERROR_NONE for successful demodulation.
 */
/*----------------------------------------------------------------------------*/
cpl_error_code gravi_metrology_demodulate (gravi_data *met_data, cpl_boolean zero_subtracted)
{
    gravi_msg_function_start(1);
    cpl_ensure_code(met_data, CPL_ERROR_NULL_INPUT);

    cpl_boolean modulation_active;
    if (cpl_propertylist_has(gravi_data_get_header(met_data), "ESO INS PMC1 MODULATE")) {
	   modulation_active = cpl_propertylist_get_bool(gravi_data_get_header(met_data), "ESO INS PMC1 MODULATE");
    } else {
	    modulation_active = 0;
	    cpl_msg_warning (cpl_func, "Can't find modulation keyword, use false");
    }

    cpl_table *met_data_table = gravi_data_get_table(met_data, GRAVI_METROLOGY_EXT);

    /* Return immediately if demodulation disabled or original metrology data is not modulated */
    if (modulation_active) {
        cpl_msg_info (cpl_func, "Demodulate the metrology");
        if (!zero_subtracted)
            cpl_msg_info(cpl_func, "Demodulating metrology with no DARK. Using hardcoded diode zero offsets");
    } else {
        cpl_msg_info (cpl_func, "Metrology modulation is not enabled");
        gravi_msg_function_exit(1);
        return CPL_ERROR_NONE;
    }

    /* Duplicate the VOLT column to store a copy of original metrology data */
    // cpl_table_duplicate_column(met_data_table, "VOLT_RAW", met_data_table, "VOLT");

    cpl_size ndata = cpl_table_get_nrow(met_data_table);
    cpl_size ncol = cpl_table_get_column_dimension(met_data_table, "VOLT", 0);

    /* Break table into chunks of up to MAX_SECONDS_PER_CHUNK */
    double exptime = cpl_propertylist_get_double(gravi_data_get_header(met_data), "EXPTIME");
    cpl_size nchunks = floor(exptime / MAX_SECONDS_PER_CHUNK) + 1;
    cpl_size seconds_per_chunk = floor(exptime / nchunks);
    cpl_size chunk_size = STEPS_PER_SECOND * seconds_per_chunk;

    /* Loop over chunks */
    for (int ichunk = 0; ichunk < nchunks; ichunk++) {
        /* Start and end indices of chunk, including remainder for last chunk */
        cpl_size start = ichunk * chunk_size;
        cpl_size count = (ichunk < nchunks - 1) ? chunk_size : ndata - start;
        cpl_table * chunk_data = cpl_table_extract(met_data_table, start, count);

        cpl_size nrow = cpl_table_get_nrow(chunk_data);
        cpl_size nsec = nrow / STEPS_PER_SECOND;
        cpl_msg_debug(cpl_func, "chunk %d has %lld rows (%lld-%lld)\n", ichunk, nrow, start, start+count);
        cpl_msg_debug(cpl_func, "seconds in this chunk: %lld\n", nsec);
        cpl_msg_debug(cpl_func, "leftover portion of a second: %lld\n", nrow % seconds_per_chunk);

        /* Extract all volt data into a GSL matrix */
        gsl_matrix * volts = gsl_matrix_alloc(nrow, ncol);
        cpl_table_cast_column(chunk_data, "VOLT", "VOLTf64", CPL_TYPE_DOUBLE);
        const cpl_array ** volt_arrs = cpl_table_get_data_array_const(chunk_data, "VOLTf64");
        for (int irow = 0; irow < nrow; irow++) {
            gsl_vector_const_view row_view = gsl_vector_const_view_array(
                cpl_array_get_data_double_const(volt_arrs[irow]), ncol
            );
            gsl_matrix_set_row(volts, irow, &row_view.vector);
        }
        CPLCHECK_INT("Cannot extract the metrology data");

        /* Calculate average over steps within chunk */
        gsl_matrix * volts_averaged = gsl_matrix_alloc(STEPS_PER_SECOND, ncol);
        for (int step = 0; step < STEPS_PER_SECOND; step++) {
            for (int icol = 0; icol < ncol; icol++) {
                /* Select the data in the appropriate column for the current step within each second */
                double *start = volts->data + (step * ncol + icol);
                gsl_vector_const_view step_col_view = gsl_vector_const_view_array_with_stride(
                    start, STEPS_PER_SECOND * ncol, nsec);
                /* NB: ignores fractional second at the end of last chunk */
                double mean = gsl_stats_mean(
                    step_col_view.vector.data, step_col_view.vector.stride, step_col_view.vector.size);
                gsl_matrix_set(volts_averaged, step, icol, mean);
            }
        }

        /* Fit demodulation model to averaged voltages */
        gsl_matrix *Xsolve = gsl_matrix_alloc(ntel * ndiode * nside, 4);
        int isolve = 0;
        for (int tel = 0; tel < ntel; tel++) {
            for (int diode = 0; diode < ndiode; diode++) {
                for (int side = FT; side <= SC; side++) {
                    int ix = diode_column_index(tel, diode, side);
                    int iy = ix + 1;
                    gsl_vector_const_view vx = gsl_matrix_const_column(volts_averaged, ix);
                    gsl_vector_const_view vy = gsl_matrix_const_column(volts_averaged, iy);

                    double ox = zero_subtracted ? 0.0 : diode_zeros[ix];
                    double oy = zero_subtracted ? 0.0 : diode_zeros[iy];

                    gsl_vector_view Xsolve_row = gsl_matrix_row(Xsolve, isolve);
                    cpl_error_code status = fit_model_modulation(&vx.vector, &vy.vector, ox, oy, &Xsolve_row.vector);
                    if (status)
                        cpl_error_set(cpl_func, status);
                    isolve++;
                }
            }
        }
        CPLCHECK_INT("Cannot fit modulation model");

        /* If DARK provided, metrology has already been zeroed by subtracting values from the DARK */
        /* Otherwise use hardcoded diode_zeros array to subtract offsets */
        gsl_matrix * volts_zeroed = gsl_matrix_alloc(nrow, ncol);
        if (zero_subtracted) {
            gsl_matrix_memcpy(volts_zeroed, volts);
        } else {
            /* Shift voltages by zero-offset */
            gsl_vector_const_view zero_view = gsl_vector_const_view_array(diode_zeros, ncol);
            for (int i = 0; i < nrow; i++) {
                gsl_matrix_set_row(volts_zeroed, i, &zero_view.vector);
            }
            gsl_matrix_scale(volts_zeroed, -1.0);
            gsl_matrix_add(volts_zeroed, volts);
        }
        
        /* Recalculate average over steps of zeroed data */
        gsl_matrix * volts_zeroed_averaged = gsl_matrix_alloc(STEPS_PER_SECOND, ncol);
        for (int step = 0; step < STEPS_PER_SECOND; step++) {
            for (int icol = 0; icol < ncol; icol++) {
                /* Select the data in the appropriate column for the current step within each second */
                double *start = volts_zeroed->data + (step * ncol + icol);
                gsl_vector_const_view step_col_view = gsl_vector_const_view_array_with_stride(
                    start, STEPS_PER_SECOND * ncol, nsec);
                /* NB: ignores fractional second at the end of last chunk */
                double mean = gsl_stats_mean(
                    step_col_view.vector.data, step_col_view.vector.stride, step_col_view.vector.size);
                gsl_matrix_set(volts_zeroed_averaged, step, icol, mean);
            }
        }

        /* Calculate phase from average */
        gsl_matrix * volts_phase = gsl_matrix_alloc(STEPS_PER_SECOND, ntel * ndiode * nside);
        int icol = 0;
        for (int tel = 0; tel < ntel; tel++) {
            for (int diode = 0; diode < ndiode; diode++) {
                for (int side = FT; side <= SC; side++) {
                    double pha2 = gsl_matrix_get(Xsolve, icol, 3);
                    int ix = diode_column_index(tel, diode, side);
                    int iy = ix + 1;
                    gsl_vector_const_view vx = gsl_matrix_const_column(volts_zeroed_averaged, ix);
                    gsl_vector_const_view vy = gsl_matrix_const_column(volts_zeroed_averaged, iy);
                    for (int step = 0; step < STEPS_PER_SECOND; step++) {
                        double phase = atan2(
                            gsl_vector_get(&vy.vector, step), gsl_vector_get(&vx.vector, step));
                        gsl_matrix_set(volts_phase, step, icol, -phase - pha2);
                    }
                    icol++;
                }
            }
        }

        /* Remove modulation component from phase */
        icol = 0;
        for (int tel = 0; tel < ntel; tel++) {
            for (int diode = 0; diode < ndiode; diode++) {
                for (int side = FT; side <= SC; side++) {
                    int ix = diode_column_index(tel, diode, side);
                    int iy = ix + 1;
                    gsl_vector_view vx = gsl_matrix_column(volts_zeroed, ix);
                    gsl_vector_view vy = gsl_matrix_column(volts_zeroed, iy);
                    for (int step = 0; step < STEPS_PER_SECOND; step++) {
                        double phase = gsl_matrix_get(volts_phase, step, icol);
                        double s = sin(phase), c = cos(phase);
                        int idx = step;
                        /* weird loop structure to correctly handle fractional second at end of chunk */
                        while (idx < nrow) {
                            double vx_orig = gsl_vector_get(&vx.vector, idx);
                            double vy_orig = gsl_vector_get(&vy.vector, idx);
                            double vx_demod = c * vx_orig - s * vy_orig;
                            double vy_demod = s * vx_orig + c * vy_orig;
                            gsl_vector_set(&vx.vector, idx, vx_demod);
                            gsl_vector_set(&vy.vector, idx, vy_demod);
                            idx += STEPS_PER_SECOND;
                        }
                    }
                    icol++;
                }
            }
        }

        /* Write back the demodulated data */
        cpl_array * cpl_row = cpl_array_new(ncol, CPL_TYPE_FLOAT);
        for (int irow = 0; irow < nrow; irow++) {
            gsl_vector_const_view row_view = gsl_matrix_const_row(volts_zeroed, irow);
            for (cpl_size i = 0; i < ncol; i++) {
                cpl_array_set_float(cpl_row, i, (float) gsl_vector_get(&row_view.vector, i));
            }
            cpl_table_set_array(met_data_table, "VOLT", start + irow, cpl_row);
        }
        cpl_array_delete(cpl_row);
        CPLCHECK_INT("Cannot store demodulated metrology data");
        
        /* Clean up intermediate values */
        gsl_matrix_free(volts_phase);
        gsl_matrix_free(volts_zeroed_averaged);
        gsl_matrix_free(volts_zeroed);
        gsl_matrix_free(Xsolve);
        gsl_matrix_free(volts_averaged);
        gsl_matrix_free(volts);
        cpl_table_delete(chunk_data);
    } /* End loop over chunks */

    gravi_msg_function_exit(1);
    return CPL_ERROR_NONE;
}

/**
 * @brief Calculate chi^2 statistic for modulation model.
 * 
 * @param X         vector of model parameters [a, b, pha1, pha2]
 * @param params    pointer to @c model_params struct
 * 
 * @return chi^2 statistic evaluated for the supplied model parameters.
*/
double modulation_model_chi2 (const gsl_vector *X, void *params)
{
    model_params *p = params;

    /* model */
    gsl_vector *model_vx = gsl_vector_alloc(STEPS_PER_SECOND);
    gsl_vector *model_vy = gsl_vector_alloc(STEPS_PER_SECOND);
    for (int i = 0; i < STEPS_PER_SECOND; i++) {
        gsl_vector_set(model_vx, i, model_x(X, p->offset_x, i));
        gsl_vector_set(model_vy, i, model_y(X, p->offset_y, i));
    }

    /* form chi2 */
    gsl_vector *chi2_x = gsl_vector_alloc(STEPS_PER_SECOND);
    gsl_vector_memcpy(chi2_x, p->volts_x);
    gsl_vector_sub(chi2_x, model_vx);
    gsl_vector_mul(chi2_x, chi2_x);

    gsl_vector *chi2_y = gsl_vector_alloc(STEPS_PER_SECOND);
    gsl_vector_memcpy(chi2_y, p->volts_y);
    gsl_vector_sub(chi2_y, model_vy);
    gsl_vector_mul(chi2_y, chi2_y);

    double chi2 = 0.0;
    for (int i = 0; i < STEPS_PER_SECOND; i++)
        chi2 += gsl_vector_get(chi2_x, i) + gsl_vector_get(chi2_y, i);

    gsl_vector_free(model_vx);
    gsl_vector_free(model_vy);
    gsl_vector_free(chi2_x);
    gsl_vector_free(chi2_y);
    
    return chi2;
}

/**
 * @brief Fit the modulation model to voltage data for a single diode.
 * 
 * @param vx x voltages
 * @param vy y voltages
 * @param zx x zero offset
 * @param zy y zero offset
 * @param[out] Xsolve returned best-fit parameters, which must already be allocated.
 * 
 * @return CPL_ERROR_NONE if fit successful
*/
cpl_error_code fit_model_modulation (const gsl_vector *vx, const gsl_vector *vy, double zx, double zy, gsl_vector *Xsolve)
{
    cpl_ensure_code(vx, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code(vy, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code(Xsolve, CPL_ERROR_NULL_INPUT);

    /* The minimisation has 4 free parameters: a, b, pha1, pha2 */
    const int dim = 4;
    
    const int MAX_ITERATIONS = 1000;
    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2rand;
    gsl_multimin_fminimizer *mini = gsl_multimin_fminimizer_alloc(T, dim);
    int status;
    double size;

    model_params params = {
        .volts_x = vx,
        .volts_y = vy,
        .offset_x = zx,
        .offset_y = zy,
    };

    gsl_multimin_function func = {
        .f = &modulation_model_chi2,
        .n = dim,
        .params = &params
    };

    double best_chi2 = INFINITY; 
    double a_guess = gsl_stats_sd(vx->data, vx->stride, vx->size);
    double b_guess = 0.25;

    double phase_range = PI;
    int nphase = 1;
    double phase_step = phase_range / (2 * nphase);

    /* Try multiple starting points for minimisation to ensure global minimum */
    for (int iph1 = -nphase; iph1 < nphase ; iph1++) {
        for (int iph2 = -nphase; iph2 < nphase; iph2++) {
            double pha1_guess = iph1 * phase_step;
            double pha2_guess = iph2 * phase_step;

            double initial_guess[] = {a_guess, b_guess, pha1_guess, pha2_guess};
            gsl_vector_view guess_view = gsl_vector_view_array(initial_guess, dim);
            double step_size[] = {0.1, 0.1, 1, 1};
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

            if (status == GSL_SUCCESS && mini->fval < best_chi2) {
                best_chi2 = mini->fval;
                gsl_vector_memcpy(Xsolve, mini->x);
            }
        }
    }

    /* Normalise fit parameters: a positive, adjust pha2 as required */
    double a, pha1, pha2;
    if ((a = gsl_vector_get(Xsolve, 0)) < 0) {
        gsl_vector_set(Xsolve, 0, -a);
        pha2 = gsl_vector_get(Xsolve, 3);
        gsl_vector_set(Xsolve, 3, pha2 + PI);
    }

    /* wrap phases */
    if (fabs(pha1 = gsl_vector_get(Xsolve, 2)) > PI)
        gsl_vector_set(Xsolve, 2, pha1 - SIGN(TWOPI, pha1));
    if (fabs(pha2 = gsl_vector_get(Xsolve, 3)) > PI)
        gsl_vector_set(Xsolve, 3, pha2 - SIGN(TWOPI, pha2));
    
    if (status == GSL_CONTINUE)
        cpl_msg_warning(cpl_func, "model-fitting did not converge");
    else if (status != GSL_SUCCESS)
        cpl_msg_error(cpl_func, "model-fitting failed with error code %d (%s)", status, gsl_strerror(status));
    else
        cpl_msg_debug(cpl_func, "converged with best-fit:\na\tb\tpha1\tpha2\n%.4f\t%.4f\t%.4f\t%.4f\n", mini->x->data[0], mini->x->data[1], mini->x->data[2], mini->x->data[3]);

    return status ? CPL_ERROR_ILLEGAL_OUTPUT : CPL_ERROR_NONE;
}
