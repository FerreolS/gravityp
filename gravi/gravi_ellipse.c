/* $Id: gravi_ellipse.c,v 1.10 2011/05/31 06:10:40 nazouaoui Exp $
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

/**
 * @defgroup gravi_ellipse  Phase and OPD estimation with ellipse
 *
 * This module implements the computation of the phase with the ellipse methode.
 * See @c gravi_ellipse_phase_create(). It also implement the function
 * @c gravi_ellipse_meanopd_create() to compute the OPD estimated overs the
 * bandpass.
 *
 */
/**@{*/

/*-----------------------------------------------------------------------------
                                   Includes
 -----------------------------------------------------------------------------*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cpl.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "gravi_data.h"

#include "gravi_dfs.h"
#include "gravi_pfits.h"
#include "gravi_cpl.h"

#include "gravi_utils.h"
#include "gravi_ellipse.h"

/*-----------------------------------------------------------------------------
                                   Function code
 -----------------------------------------------------------------------------*/

static int ellipse(const double x_in[], const double v[], double *result) {

	double x = x_in[0], y = x_in[1], a = v[0], b = v[1], c = v[2], d = v[3],
			e = v[4];
	*result = sqrt(pow(a * x + b * y + c, 2) + pow(d * y + e, 2));
    
	return (0);
}

static int dfda_ellipse(const double x_in[], const double v[], double result[]) {

	double x = x_in[0], y = x_in[1], a = v[0], b = v[1], c = v[2], d = v[3],
		   e = v[4];
	float inv_sqrtf;
	inv_sqrtf = 1 / sqrt(pow(a * x + b * y + c, 2) + pow(d * y + e, 2));

	result[0] = inv_sqrtf * (2 * a * x * x + 2 * x * b * y + 2 * x * c);
	result[1] = inv_sqrtf * (2 * b * y * y + 2 * a * x * y + 2 * y * c);
	result[2] = inv_sqrtf * (2 * a * x + 2 * b * y + 2 * c);
	result[3] = inv_sqrtf * (2 * y * e + 2 * y * y * d);
	result[4] = inv_sqrtf * (2 * e + 2 * d * y);

	return (0);
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Compute the phase atan{X',Y'}, unwraped from first sample
 * 
 * @param vectCA     Input vector
 * @param vectCA     Input vector
 * @param envelope   Input vector
 *
 * @return A newly allocated vector with the computed unwrapped phase [rad]
 *
 * The quantities vectCA and vectDB are fitted with an ellipse of amplitude
 * envelope, and converted into X' and Y' (see pipeline design document).
 * Then the phase is computed and unwrapped.
 *
 * If the variable USE_LINEAR_ELLIPSE is 1, the functions uses a linear
 * fit. Otherwise it uses a longuer and possibly less stable minimisation.
 * 
 * vectCA, vectDB and envelope shall have the same size and ideally shall
 * modulate over several fringes with good phase diversity in order for
 * the ellipse fit to be well defined.
 */
/*---------------------------------------------------------------------------*/

cpl_vector * gravi_ellipse_phase_create (cpl_vector * vectCA,
                                         cpl_vector * vectDB,
                                         cpl_vector * envelope)
{
    gravi_msg_function_start(0);
	cpl_ensure (vectCA, CPL_ERROR_NULL_INPUT, NULL);
	cpl_ensure (vectDB, CPL_ERROR_NULL_INPUT, NULL);

	cpl_size size = cpl_vector_get_size (vectCA);

    /* Linear ellipse fitting */
    if (USE_LINEAR_ELLIPSE) {
        cpl_errorstate prestate = cpl_errorstate_get();
        
        /* Recenter the vector in-place */
        cpl_vector_subtract_scalar (vectCA, cpl_vector_get_mean (vectCA));
        cpl_vector_subtract_scalar (vectDB, cpl_vector_get_mean (vectDB));
        
        /* Normalise the vector in-place */
        cpl_vector_divide_scalar (vectCA, cpl_vector_get_stdev (vectCA));
        cpl_vector_divide_scalar (vectDB, cpl_vector_get_stdev (vectDB));
        
        /* Fill matrix */
        cpl_matrix * coeff = cpl_matrix_new (size, 5);
        for (cpl_size v=0; v<size; v++) {
            double x = cpl_vector_get (vectCA, v);
            double y = cpl_vector_get (vectDB, v);
            cpl_matrix_set (coeff, v, 0, x*x);
            cpl_matrix_set (coeff, v, 1, x*y);
            cpl_matrix_set (coeff, v, 2, x);
            cpl_matrix_set (coeff, v, 3, y);
            cpl_matrix_set (coeff, v, 4, y*y);
        }

        /* Fill right-hand-side */
        cpl_matrix * rhs   = cpl_matrix_new (size, 1);
        for (cpl_size v=0; v<size; v++) {
            double value = envelope ? cpl_vector_get (envelope, v) : 1.0;
            cpl_matrix_set (rhs, v, 0, value);
        }
        
        /* Solve */
        cpl_matrix * res = cpl_matrix_solve_normal (coeff, rhs);
        
        /* Delete */
        FREE (cpl_matrix_delete, coeff);
        FREE (cpl_matrix_delete, rhs);
        
        /* Recover on error but return NULL */
        if (cpl_error_get_code()) {
            cpl_errorstate_set (prestate);
            cpl_msg_warning (cpl_func, "Error during fit of ellipse");
            return NULL;
        }
        
        /* Convert to Lacour coefficients */
        double b  = cpl_matrix_get (res,4,0);
        if (b <= 0) {cpl_msg_warning (cpl_func, "Error during fit of ellipse");
            FREE (cpl_matrix_delete, res); return NULL;}
        
        double d  = cpl_matrix_get (res,1,0) / (2*b);
        double a  = cpl_matrix_get (res,0,0) - b*d*d;
        if (a <= 0) {cpl_msg_warning (cpl_func, "Error during fit of ellipse");
            FREE (cpl_matrix_delete, res); return NULL;}
        
        double c2 = cpl_matrix_get (res,3,0) / (2*b);
        double c1 = (cpl_matrix_get (res,2,0) - 2*b*d*c2) / (2*a);
        
        FREE (cpl_matrix_delete, res);
        
        /* Replace vectCA and vectDB by their corrected version */
        a = sqrt(a);
        b = sqrt(b);
        for (cpl_size v=0; v<size; v++) {
            double x = cpl_vector_get (vectCA, v);
            double y = cpl_vector_get (vectDB, v);
            cpl_vector_set (vectCA, v, a * (x+c1));
            cpl_vector_set (vectDB, v, b * (y+d*x+c2));
        }
        
    }
    /* Non-linear ellipse fitting */
    else {
        /* Initialization of th init_val */
        cpl_vector * init_val = cpl_vector_new(5);
        cpl_vector_set(init_val, 0, 1.);
        cpl_vector_set(init_val, 1, 0.);
        cpl_vector_set(init_val, 2, (-1) * cpl_vector_get_mean(vectCA));
        cpl_vector_set(init_val, 3, 1);
        cpl_vector_set(init_val, 4, (-1) * cpl_vector_get_mean(vectDB));
        
        /* Construct the matrix_XY and vector_1 */
        cpl_matrix * matrix_XY = get_matrix_from_vector(vectCA, vectDB);
        cpl_vector * vector_1 = cpl_vector_new (size);
        for (cpl_size v=0; v<size; v++) {
            double value = envelope ? cpl_vector_get (envelope, v) : 1.0;
            cpl_vector_set (vector_1, v, value);
        }
                
        /* Fit ellipse */
        cpl_errorstate prestate = cpl_errorstate_get();
        
        double mse = 0;
        int val_to_fit[] = {1,1,1,1,1};
        cpl_fit_lvmq(matrix_XY, NULL, vector_1, NULL, init_val, val_to_fit,
                     &ellipse, &dfda_ellipse, CPL_FIT_LVMQ_TOLERANCE,
                     CPL_FIT_LVMQ_COUNT, CPL_FIT_LVMQ_MAXITER, &mse, NULL, NULL);
        
        /* Recover on error but return NULL */
        if (cpl_error_get_code()){
            cpl_vector_delete(vector_1);
            cpl_vector_delete(init_val);
            cpl_matrix_delete(matrix_XY);
            cpl_errorstate_set(prestate);
            cpl_msg_warning(cpl_func, "Error during fit of ellipse");
            return NULL;
        }
        
        /* Replace the new vector_X values after the fit where
         * vector_X = init_val(0) * vector_X + init_val(1) *  vector_Y + init_val(2)
         * vector_Y = init_val(3) * vector_Y + init_val(4) */
        
        /* Return the corrected vector */
        double vector_Yp_i, vector_X_i, vector_Y_i;
        for (cpl_size i_data = 0; i_data < size; i_data++) {
            
            vector_Yp_i = cpl_vector_get(vectDB, i_data) *
                cpl_vector_get(init_val, 1);
            vector_X_i = cpl_vector_get(vectCA, i_data) *
                cpl_vector_get(init_val, 0) + cpl_vector_get(init_val, 2);
            
            cpl_vector_set (vectCA, i_data, vector_X_i + vector_Yp_i);
            
            vector_Y_i = cpl_vector_get (vectDB, i_data) *
                cpl_vector_get(init_val, 3) + cpl_vector_get (init_val, 4);
            
            cpl_vector_set(vectDB, i_data, vector_Y_i);
        }
        
        cpl_matrix_delete (matrix_XY);
        cpl_vector_delete (init_val);
        cpl_vector_delete (vector_1);
        
    } /* End non-linear ellipse fitting */
        
    

	/* Initialization of the OPD vector and compute the phase*/
	cpl_vector * OPD_rad = cpl_vector_new(size);
	cpl_vector_set (OPD_rad, 0,  atan2(cpl_vector_get(vectDB, 0),
                                       cpl_vector_get(vectCA, 0)));

	/* Dewrap opd */
    double d_opd_i, d_opd_ii, d_opd;
	for (cpl_size i_data = 1; i_data < size; i_data++){

		/* Evaluation of the OPD(i_data) and
		 * OPD(i_data-1) */
		d_opd_i = atan2 (cpl_vector_get(vectDB, i_data),
                         cpl_vector_get(vectCA, i_data));

		d_opd_ii = atan2 (cpl_vector_get(vectDB, i_data - 1),
                          cpl_vector_get(vectCA, i_data - 1));
        
		d_opd = d_opd_i - d_opd_ii;

		if (d_opd > M_PI) d_opd = d_opd - 2 * M_PI;
		if (d_opd < - M_PI) d_opd = d_opd + 2 * M_PI;

		cpl_vector_set (OPD_rad, i_data, cpl_vector_get (OPD_rad, i_data-1) +
                        d_opd);
	}

	/* End fit opd */
    gravi_msg_function_exit(0);
	return OPD_rad;
}

/*----------------------------------------------------------------------------*/

cpl_vector * gravi_ellipse_phase_create_fast (cpl_vector * vectCA, cpl_vector * vectDB)
{
    gravi_msg_function_start(0);
	cpl_ensure (vectCA, CPL_ERROR_NULL_INPUT, NULL);
	cpl_ensure (vectDB, CPL_ERROR_NULL_INPUT, NULL);

	cpl_size size = cpl_vector_get_size (vectCA);

    /* Recenter the vector in-place */
    cpl_vector_subtract_scalar (vectCA, cpl_vector_get_mean (vectCA));
    cpl_vector_subtract_scalar (vectDB, cpl_vector_get_mean (vectDB));

	/* Initialization of the OPD vector and compute the phase*/
	cpl_vector * OPD_rad = cpl_vector_new(size);
	cpl_vector_set (OPD_rad, 0,  atan2(cpl_vector_get(vectDB, 0),
                                       cpl_vector_get(vectCA, 0)));

	/* Dewrap opd */
    double d_opd_i, d_opd_ii, d_opd;
	for (cpl_size i_data = 1; i_data < size; i_data++){

		/* Evaluation of the OPD(i_data) and
		 * OPD(i_data-1) */
		d_opd_i = atan2 (cpl_vector_get(vectDB, i_data),
                         cpl_vector_get(vectCA, i_data));

		d_opd_ii = atan2 (cpl_vector_get(vectDB, i_data - 1),
                          cpl_vector_get(vectCA, i_data - 1));
        
		d_opd = d_opd_i - d_opd_ii;

		if (d_opd > M_PI) d_opd = d_opd - 2 * M_PI;
		if (d_opd < - M_PI) d_opd = d_opd + 2 * M_PI;

		cpl_vector_set (OPD_rad, i_data, cpl_vector_get (OPD_rad, i_data-1) +
                        d_opd);
	}

	/* End fit opd */
    gravi_msg_function_exit(0);
	return OPD_rad;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Compute the OPD modulation of a baseline from spectrum
 * 
 * @param spectrum_table     The input spectrum
 * @param detector_table     The corresponding detector_table
 * @param oiwave_tables      The wavelength tables(s) 
 * @param guess_vector       An optional vector with guess [m]
 *
 * The routine compute the OPD of each spectral channel with the
 * ellipse. These are averaged over channel and polarisations to
 * provide a single value per baseline and per DIT.
 * 
 * The group-delay is also computed internally, to then ensure that
 * the returned OPD = 0 [m] corresponds to group-delay = 0 [m].
 *
 * If a guess is provided (shall be in [m]), the phase of each channel 
 * is first unwrap following this guess before averaging.
 */
/*----------------------------------------------------------------------------*/

cpl_vector * gravi_ellipse_meanopd_create (cpl_table * spectrum_table,
                                           cpl_table * detector_table,
                                           cpl_table ** oiwave_tables,
                                           cpl_vector * guess_vector,
                                           int base)
{
	gravi_msg_function_start(0);
    cpl_ensure (spectrum_table, CPL_ERROR_NULL_INPUT, NULL);
    cpl_ensure (detector_table, CPL_ERROR_NULL_INPUT, NULL);
    cpl_ensure (oiwave_tables,  CPL_ERROR_NULL_INPUT, NULL);

    /* Get the size of the vectors */
    cpl_size nwave = gravi_spectrum_get_nwave (spectrum_table);
    cpl_size nrow  = cpl_table_get_nrow (spectrum_table);
    int npol = (cpl_table_get_nrow (detector_table) > 24 ? 2 : 1);
    CPLCHECK_NUL ("Cannot get size");

    /* Check input */
    for (int pol = 0; pol < npol; pol++)
        cpl_ensure (oiwave_tables[pol], CPL_ERROR_NULL_INPUT, NULL);

    /* Spectral band to integrate */
    int wave_start = nwave > GRAVI_LBD_FTSC ? nwave/4 : 1;
    int wave_end   = nwave > GRAVI_LBD_FTSC ? (nwave*3)/4 : nwave-2;
    
    /* Init mean_opd to average over channel and polar */
    cpl_vector * mean_opd = cpl_vector_new (nrow);
    cpl_vector_fill (mean_opd, 0.0);
    
    /*
    if (base == 0)
    {
                                FILE *fAst=fopen("opd_text.txt","a");
    fprintf(fAst, "%i\n" , wave_start);
        fprintf(fAst, "%i\n" , wave_end);
        fprintf(fAst, "%i\n" , 555);
                                fclose(fAst);
    }*/

    /* Do the computation twice, first without enveloppe and 
     * then with enveloppe correction */
    for (int loop = 0; loop < 2; loop++) {
        
        cpl_msg_debug (cpl_func, "Now run %s envelope calibration", loop ? "WITH" : "WITHOUT");
        
        cpl_vector * opd_vector = cpl_vector_duplicate (mean_opd);
        cpl_vector_fill (mean_opd, 0.0);
        
        /* Save the OPD to compute the enveloppe
         * accurately in the second loop */
        
        /* Loop on polarisations */
        for (int pol = 0; pol < npol; pol ++) {
            
            /* Sign of this baseline */
            double phi_sign = gravi_region_get_base_sign (detector_table, base);
            
            /* Get the region index of this base and pol */
            int iA = gravi_get_region (detector_table, base, 'A', pol);
            int iB = gravi_get_region (detector_table, base, 'B', pol);
            int iC = gravi_get_region (detector_table, base, 'C', pol);
            int iD = gravi_get_region (detector_table, base, 'D', pol);
            if (iA<0 || iB<0 || iC<0 || iD<0){
                cpl_msg_error (cpl_func, "Don't found the A, B, C or D !!!");
            }

            /* Init the mean phase over channel for this pol */
            cpl_vector * pol_opd = cpl_vector_new (nrow);
            cpl_vector_fill (pol_opd, 0.0);
		
            /* Init the mean inter-spectra over channel for this pol */
            cpl_vector * phase_previous = NULL;
            cpl_array * is_array = cpl_array_new (nrow, CPL_TYPE_DOUBLE_COMPLEX);
            cpl_array_fill_window_double_complex (is_array, 0, nrow, 0 + 0.0*I);

            /* Loop on wave */
            for (cpl_size wave = wave_start; wave <= wave_end; wave++){
                cpl_vector * vector_T = NULL;
                                
                /* Define the vector_X = C - A */
                cpl_vector * vector_X;
                vector_X = gravi_table_get_vector (spectrum_table, wave, GRAVI_DATA[iC]);
                vector_T = gravi_table_get_vector (spectrum_table, wave, GRAVI_DATA[iA]);
                cpl_vector_subtract (vector_X, vector_T);
                FREE (cpl_vector_delete, vector_T);
                
                /* Define the vector_Y = D - B */
                cpl_vector * vector_Y;
                vector_Y = gravi_table_get_vector (spectrum_table, wave, GRAVI_DATA[iD]);
                vector_T = gravi_table_get_vector (spectrum_table, wave, GRAVI_DATA[iB]);
                cpl_vector_subtract (vector_Y, vector_T);
                FREE (cpl_vector_delete, vector_T);
                
                if (wave_end < 10)
                {
                    cpl_vector* vector_X1= cpl_vector_duplicate (vector_X);
                    cpl_vector* vector_Y1= cpl_vector_duplicate (vector_Y);
                    double new_scalar;
                    for (cpl_size n=1; n < cpl_vector_get_size(vector_X)-1;n++)
                    {
                        new_scalar=cpl_vector_get(vector_X1,n-1);
                        new_scalar+=cpl_vector_get(vector_X1,n);
                        new_scalar+=cpl_vector_get(vector_X1,n+1);
                        cpl_vector_set(vector_X,n,new_scalar/3);
                    }
                    for (cpl_size n=1; n < cpl_vector_get_size(vector_Y)-1;n++)
                    {
                        new_scalar=cpl_vector_get(vector_Y1,n-1);
                        new_scalar+=cpl_vector_get(vector_Y1,n);
                        new_scalar+=cpl_vector_get(vector_Y1,n+1);
                        cpl_vector_set(vector_Y,n,new_scalar/3);
                    }
                    FREE (cpl_vector_delete, vector_X1);
                    FREE (cpl_vector_delete, vector_Y1);
                }
                /*
                if ((base == 0)&&(wave_end == 4))
                {
                    FILE *fAst=fopen("opd_text.txt","a");
                    fprintf(fAst, "%i\n" , 500);
                    for (cpl_size n=0; n < cpl_vector_get_size(vector_X);n++)
                        fprintf(fAst, "%g\n" , cpl_vector_get(vector_X,n) );
                    fprintf(fAst, "%i\n" , 600);
                    for (cpl_size n=0; n < cpl_vector_get_size(vector_Y);n++)
                        fprintf(fAst, "%g\n" , cpl_vector_get(vector_Y,n) );
                    fclose(fAst);
                }
                 */
                
                CPLCHECK_NUL ("Cannot extract the X and Y vectors");

                /* Compute envelope from OPD for this channel */
                cpl_vector * envelope_vector = gravi_compute_envelope (opd_vector, wave, nwave);
                
                /* Compute the phase of this channel from ellipse */
                cpl_vector * phase = gravi_ellipse_phase_create (vector_X, vector_Y, envelope_vector);
                cpl_vector_delete (vector_X);
                cpl_vector_delete (vector_Y);
                cpl_vector_delete (envelope_vector);
                CPLCHECK_NUL ("Error during the fit of the ellipse");
                
                /* Multiply by the sign */    
                cpl_vector_multiply_scalar (phase, phi_sign);

                /* Wavelength of this channel */
                double lbd_channel = cpl_table_get (oiwave_tables[pol], "EFF_WAVE", wave, NULL);
                
                /* Unwrap if opd_guess is provided */
                if (guess_vector) {
                    gravi_vector_unwrap_with_guess (phase, guess_vector, CPL_MATH_2PI / lbd_channel);
                    CPLCHECK_NUL ("Error during the unwrap");
                }

                /* Accumulate the inter-spectra over the spectral channels */
                if (wave != wave_start) {
                    cpl_vector_subtract (phase_previous, phase);
                    for (cpl_size row=0; row<nrow; row++) {
                        cpl_array_set_complex (is_array, row,
                                               cpl_array_get_complex (is_array, row, NULL) +
                                               cexp ( 1.*I*cpl_vector_get (phase_previous, row)));
                        CPLCHECK_NUL ("Cannot accumulate is_array");
                    }
                    FREE (cpl_vector_delete, phase_previous);
                }
                if (wave != wave_end)
                    phase_previous = cpl_vector_duplicate (phase);
                
                /*
                if ((base == 0)&&(wave_end == 4))
                {
                    FILE *fAst=fopen("opd_text.txt","a");
                    fprintf(fAst, "%i\n" , 400);
                    for (cpl_size n=0; n < cpl_vector_get_size(phase);n++)
                        fprintf(fAst, "%g\n" , cpl_vector_get(phase,n) );
                    fclose(fAst);
                }*/
                
                /* Convert phase to opd and accumulate the mean
                 * phase opd over the spectral channels */
                cpl_vector_multiply_scalar (phase, lbd_channel / CPL_MATH_2PI);
                cpl_vector_add (pol_opd, phase);
                CPLCHECK_NUL ("Cannot accumulate phase");
                
                
                FREE (cpl_vector_delete, phase);
            } /* end loop on wave */
		

            /* Compute mean phase over the channels [rad] */
            cpl_vector_divide_scalar (pol_opd, wave_end - wave_start + 1);
            
            /* Compute the group-delay [rad/cannal] */
            cpl_array_arg (is_array);

            /* Fit group-delay [rad/cannal] versus opd [um] 
             * with a linear relation */
            cpl_polynomial * fit_lin =	cpl_polynomial_new(1);
            cpl_matrix * matrix = cpl_matrix_wrap (1, nrow, cpl_array_get_data_double(is_array));
            const cpl_size  mindeg = 0, maxdeg = 1;
            
            cpl_polynomial_fit (fit_lin , matrix, NULL, pol_opd, NULL,
                                CPL_FALSE, &mindeg, &maxdeg);
            
            /* Get the opd where the group-delay == 0 */
            double opd0 = cpl_polynomial_get_coeff (fit_lin, &mindeg);

            cpl_matrix_unwrap (matrix);
            cpl_array_delete (is_array);
            cpl_polynomial_delete (fit_lin);
            
            CPLCHECK_NUL ("Cannot fit");
            
            /* Remove opd0 from the mean phase to ensure
             * pol_opd = 0   <->    gdelay = 0 */
            cpl_vector_subtract_scalar (pol_opd, opd0);
            
            cpl_msg_info (cpl_func, "opd0 = %g [um] (base %d pol %d %s envelope)",
                          opd0*1e6, base + 1, pol + 1, loop ? "with" : "without");

            /* Integrathe the two polarisation */
            cpl_vector_add (mean_opd, pol_opd);
            cpl_vector_delete (pol_opd);
            
        } /* End loop on pol */
        
        cpl_vector_divide_scalar (mean_opd, npol);

        FREE (cpl_vector_delete, opd_vector);
    } /* End loop on with/without correction */
    
	gravi_msg_function_exit(0);
    return mean_opd;
}

/**@}*/
