/* $Id: gravi_calib.c,v 1.10 2012/03/23 15:10:40 nazouaoui Exp $
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
 * @defgroup gravi_wave     Spectral calibration
 *
 * This module implements functions involved in the spectral calibration. The
 * function @c gravi_compute_wave() do the calibration itself :
 *
 * The wavelength of each spectral element is computed by comparing the measured
 * phases of this spectral element with the realized OPD: OPDFT or OPDSC.
 * The measured phases are computed from the A, B, C and D measurements with ellipse
 * methode (see module Phase and OPD estimation with ellipse).
 *
 * For each computed phase we know the expected OPD, OPDFT or OPD_SC from the metrology (@c gravi_compute_opds()).
 * The slope of the phase versus OPD gives us the wavelength of the spectral element.
 *
 * When all spectral element wavelengths are computed we have two sets of calibrated points,
 * one for each polarization. On each of these two sets, a model of lambda versus
 * position on the detector is fitted. And from this the wavelength of each spectral
 * element of each spectrum is computed and put in the wavelength map.
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
#include "cpl_wlcalib.h"
#include "irplib_wavecal.h"
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <complex.h>

#include "gravi_data.h"
#include "gravi_dfs.h"
#include "gravi_pfits.h"
#include "gravi_cpl.h"

#include "gravi_utils.h"
#include "gravi_ellipse.h"
#include "gravi_signal.h"

#include "gravi_preproc.h"
#include "gravi_calib.h"
#include "gravi_wave.h"
#include "gravi_metrology.h"

/*----------------------------------------------------------------------------
                                    DEBUG
 -----------------------------------------------------------------------------*/

#define PLOT_WAVE_PHASE_VS_OPD 0
#define WAVE_TO_PLOT 5

/*-----------------------------------------------------------------------------
                                 Private prototypes
 -----------------------------------------------------------------------------*/

cpl_table * gravi_opds_compute_guess (cpl_table * spectrumsc_table,
                                      cpl_table * opdft_table,
                                      cpl_table * met_table,
                                      double dit_sc,
                                      double lbd_met);

cpl_table * gravi_opds_calibration (cpl_table * spectrum_table,
                                    cpl_table * detector_table,
                                    cpl_table * opdguess_table);

cpl_error_code gravi_opds_correct_closures (cpl_table * opd_sc,
                                            const char *name);

cpl_vector * gravi_opds_fit_opdmet (cpl_table *opd_ft, double lbd_met);

cpl_table * gravi_wave_fibre (cpl_table * spectrum_table,
                              cpl_table * detector_table,
                              cpl_table * opd_table);

cpl_table * gravi_wave_fit_2d (cpl_table * wavefibre_table,
                               cpl_table * detector_table,
                               gravi_data * wave_param,
                               cpl_size fullstartx,
                               int spatial_order,
                               int spectral_order,
                               double * rms_residuals);

cpl_table * gravi_wave_fit_individual (cpl_table * wave_individual_table,
                                          cpl_table * weight_individual_table,
                                          cpl_table * wave_fitted_table,
                                          cpl_table * opd_table,
                                          cpl_table * spectrum_table,
                                          cpl_table * detector_table,
                                          cpl_size fullstartx,
                                       double n0, double n1, double n2,
                                       double * rms_residuals);

cpl_error_code gravi_wave_correct_dispersion (cpl_table * wave_fibre,
                                              double n0, double n1, double n2);

cpl_imagelist * gravi_wave_test_image (cpl_table * wavedata_table,
                                       cpl_table * wavefibre_table,
                                       cpl_table * profile_table,
                                       cpl_table * detector_table);




/*-----------------------------------------------------------------------------
                                Functions code
 -----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
 * @brief Compute a WAVE calibration from the ARGON data (SC only)
 *
 * @param spectrum_table   The input spectrum of each region
 *
 * @return wave_table      The output WAVE table of each region
 * 
 * The spectra of each region are adjusted independently with an argon
 * lines and a slit model.
 */
/*----------------------------------------------------------------------------*/

cpl_table * gravi_compute_argon_wave (cpl_table * spectrum_table)
{
    gravi_msg_function_start(1);
    cpl_ensure (spectrum_table, CPL_ERROR_NULL_INPUT, NULL);
	
	/* Get the SC spectrums */
	cpl_size nwave = cpl_table_get_column_depth (spectrum_table,"DATA1");
	cpl_size nregion = gravi_spectrum_get_nregion (spectrum_table);
	CPLCHECK_NUL ("Cannot get data");

	cpl_msg_info (cpl_func, "spectrum_table = %lld",cpl_table_get_column_depth (spectrum_table, "DATA1"));

	/* Prepare the outputs */
	cpl_table * wave_data_sc = cpl_table_new (1);
	CPLCHECK_NUL ("Cannot create output");

	/* Init the list of lines*/
	cpl_size nlines = 14;
	double plines_str[] = {60,
				   160,
			       90,
			       1800,
			       50,
			       30,
			       30,
			       1500,
			       130,
			       510,
			       240,
			       169,
			       64,
			       73
	};
	
	double plines_lbd[] = {1.982291,
				   1.997118,
			       2.032256,
			       2.062186,
			       2.065277,
			       2.073922,
			       2.081672,
			       2.099184,
			       2.133871,
			       2.154009,
			       2.208321,
			       2.313952,
			       2.385154,
			       2.397306
	};
	
	cpl_vector * lines_lbd = cpl_vector_wrap (nlines,plines_lbd);
	cpl_vector * lines_str = cpl_vector_wrap (nlines,plines_str);
	cpl_bivector * lines = cpl_bivector_wrap_vectors (lines_lbd, lines_str);
	
	/* Init the slit and lamp model */
	cpl_wlcalib_slitmodel * model = cpl_wlcalib_slitmodel_new ();
	cpl_wlcalib_slitmodel_set_wslit (model, 0.1);
	cpl_wlcalib_slitmodel_set_wfwhm (model, 1.25);
	cpl_wlcalib_slitmodel_set_threshold (model, 5.0);
	cpl_wlcalib_slitmodel_set_catalog (model, lines);

	/* Starting point of dispersion */
	double values_hr[] = {1.97048,0.2228e-3,2.7728e-08,-4.72e-12};
	double values_mr[] = {1.9697,0.2179e-2,2.5e-07,-5e-10};
	double * values = nwave > 1000 ? values_hr : values_mr;
	
	/* Fill the dispersion polynomial */
	cpl_polynomial * dispersion0 = cpl_polynomial_new (1);
	for (cpl_size order = 0 ; order < 4 ; order ++) {
	  cpl_polynomial_set_coeff (dispersion0, &order, values[order]);
	}

	/* QC parameters */
	double minwave = -1e10;
	double maxwave = +1e10;

	cpl_vector * spectrum = cpl_vector_new (nwave);
	cpl_vector * input_spectrum = cpl_vector_new (nwave);
	cpl_vector * model_spectrum = cpl_vector_new (nwave);

	/* Prepare output table */
	cpl_table * fit_table = cpl_table_new (nregion);
    gravi_table_new_column_array (fit_table, "DATA", "adu", CPL_TYPE_DOUBLE, nwave);
    gravi_table_new_column_array (fit_table, "DATA_MODEL", "adu", CPL_TYPE_DOUBLE, nwave);
    gravi_table_new_column_array (fit_table, "DATA_WAVE", "adu", CPL_TYPE_DOUBLE, nwave);

	/* The starting point will be the fit of previous region */
	cpl_polynomial * dispersion = cpl_polynomial_duplicate (dispersion0);

	/* Loop on regions */
	for (cpl_size reg = 0; reg < nregion ; reg ++ ) {
	  cpl_msg_info (cpl_func,"Fit region %lld over %lld", reg+1, nregion);

	  /* Copy the spectrum of this region into a vector */
	  const cpl_array * data = cpl_table_get_array (spectrum_table, GRAVI_DATA[reg], 0);
	  for (cpl_size wave=0; wave < nwave; wave ++)  {
	    cpl_vector_set (spectrum, wave, log (CPL_MAX (cpl_array_get (data, wave, NULL), 1.0)));
	    cpl_vector_set (input_spectrum, wave, cpl_array_get (data, wave, NULL));
	  }
	  
	  CPLCHECK_NUL ("Cannot prepare data");

	  /* Parameters: Explore a large hsize [pixel]
	   * only for the first region, since we keep
	   * the starting point */
	  int maxdeg = 3, nmaxima = 0, linelim = 99;
	  int maxite = 100 * maxdeg, maxfail = 10, maxcont = 10;
	  int hsize = (reg == 0 ? 40 : 5);
	  double pixtol = 1e-6, pixstep = 0.5, pxc = 0.0;

	  /* Use the log of spectrum to put
	   * less weight on the intensity */
	  irplib_polynomial_find_1d_from_correlation_all (dispersion, maxdeg, spectrum,
													  nmaxima, linelim,
													  (irplib_base_spectrum_model*)model,
													  &irplib_vector_fill_logline_spectrum,
													  pixtol, pixstep,
													  hsize, maxite, maxfail, maxcont,
													  CPL_FALSE,&pxc);
	  CPLCHECK_NUL ("Cannot fit polynomial");

	  /* Fill in the final fit */
	  irplib_vector_fill_line_spectrum (model_spectrum, dispersion, (void *)model);
	  cpl_polynomial_dump (dispersion, stdout);

	  /* Evaluate the polynomial for the output */
	  cpl_array * wave_data = cpl_array_new (nwave, CPL_TYPE_DOUBLE);
	  for (cpl_size pix = 0; pix < nwave ; pix++) 
	    cpl_array_set (wave_data, pix, cpl_polynomial_eval_1d (dispersion, (double)pix, NULL) * 1e-6);

	  /* Check the limits for the QC */
	  minwave = CPL_MAX (cpl_array_get_min (wave_data), minwave);
	  maxwave = CPL_MIN (cpl_array_get_max (wave_data), maxwave);
	  
	  /* Fill the WAVE_DATA_SC */
	  cpl_table_new_column_array (wave_data_sc, GRAVI_DATA[reg], CPL_TYPE_DOUBLE, nwave);
	  cpl_table_set_array (wave_data_sc, GRAVI_DATA[reg], 0, wave_data);
	  FREE (cpl_array_delete, wave_data);

	  /* Save the test tables */
      cpl_array * tmp_array;
      tmp_array = cpl_array_wrap_double (cpl_vector_get_data (input_spectrum), nwave);
      cpl_table_set_array (fit_table, "DATA", reg, tmp_array);
      cpl_array_unwrap (tmp_array);
      tmp_array = cpl_array_wrap_double (cpl_vector_get_data (model_spectrum), nwave);
      cpl_table_set_array (fit_table, "DATA_MODEL", reg, tmp_array);
      cpl_array_unwrap (tmp_array);
      tmp_array = cpl_array_new (nwave, CPL_TYPE_DOUBLE);
	  for (cpl_size wave = 0; wave < nwave; wave ++) 
          cpl_array_set (tmp_array, wave, cpl_polynomial_eval_1d (dispersion, (double)wave, NULL) * 1e-6);
      cpl_table_set_array (fit_table, "DATA_WAVE", reg, tmp_array);
      cpl_array_delete (tmp_array);

	  CPLCHECK_NUL ("Cannot set data");
	} /* End loop on regions */

	// cpl_table_save (fit_table, NULL, NULL, "fit_wave.fits", CPL_IO_CREATE);
    FREE (cpl_table_delete, fit_table);

	/* Delete */
	FREE (cpl_vector_delete, spectrum);
	FREE (cpl_vector_delete, input_spectrum);
	FREE (cpl_vector_delete, model_spectrum);
	FREE (cpl_polynomial_delete, dispersion0);
	FREE (cpl_polynomial_delete, dispersion);

	/* This memory desalocation does not work: */
	// FREE (cpl_wlcalib_slitmodel_delete, model);
	
	gravi_msg_function_exit(1);
	return wave_data_sc;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Correct the input table to ensure constant closure phase
 *
 * @param phase_table    The table to be corrected in-place
 * @param name           The column to correct (PHASE or OPD)
 *
 * The raw broad-band phase from gravi_opds_calibration doesn't necessarely
 * form a closing system (constant closure phase versus modulation), if the
 * the effective broad-band wavelength of the baselines are differents. For
 * instance because of detector distortion. This routine determines the scaling
 * coefficients and reconstructs a closing system.
 *
 * The routine find the 4 k and 5 coefficients solving the following
 * linear system of closure phase:
 *
 *     B0 + B3 - B1  =  f0 B0 + f3 B3 - f1 B1 + k0
 *     B0 + B4 - B2  =  f0 B0 + f4 B3 - f2 B2 + k1
 *     B1 + B5 - B2  =  f1 B1             - f2 B2 + k2
 *     B3 + B5 - B4  =  f3 B3             - f4 B4 + k3
 *
 * Then it re-scales the first 5 baselines:
 *     B0' = (1-f0).B0
 *     B1' = (1-f1).B1
 *     B2' = (1-f2).B2
 *     B3' = (1-f3).B3
 *     B4' = (1-f4).B4
 * 
 * The 'name' column shall be a scalar colum, and have NDIT * NBASE rows.
 */
/*----------------------------------------------------------------------------*/
    
cpl_error_code gravi_opds_correct_closures (cpl_table * phase_table,
                                            const char *name)
{
    gravi_msg_function_start(1);
	cpl_ensure_code (phase_table, CPL_ERROR_NULL_INPUT);
	cpl_ensure_code (name,        CPL_ERROR_NULL_INPUT);
    
    int nbase = 6, nclo = 4;
    cpl_size nrow = cpl_table_get_nrow (phase_table) / nbase;

    /* Get the data */
    double * phase = cpl_table_get_data_double (phase_table, name);
    CPLCHECK_MSG ("Cannot get data");

    /* Model and right-hand-side for the lineary system
     * (unfilled matrix are 0.0) */
    cpl_matrix * rhs_matrix   = cpl_matrix_new (nrow * nclo, 1);
    cpl_matrix * model_matrix = cpl_matrix_new (nrow * nclo, nbase-1 + nclo);

    for (cpl_size row = 0; row < nrow; row++) {
        for (int clo = 0; clo < nclo; clo++) {
            int b0 = GRAVI_CLO_BASE[clo][0];
            int b1 = GRAVI_CLO_BASE[clo][1];
            int b2 = GRAVI_CLO_BASE[clo][2];

            /* Fill the rhs with measured closure */
            double closure = phase[row*nbase+b0] +
                             phase[row*nbase+b1] -
                             phase[row*nbase+b2];
            cpl_matrix_set (rhs_matrix, row*nclo+clo, 0, closure);
            CPLCHECK_MSG ("Cannot fill rhs_matrix");

            /* Fill the k (unfilled are 0.0) */
            cpl_matrix_set (model_matrix, row*nclo+clo, clo, 1.0);
            CPLCHECK_MSG ("Cannot fill k in model_matrix");

            /* Fill the f in of the model_matrix */
            if (b0 < nbase-1) cpl_matrix_set (model_matrix, row*nclo+clo, nclo+b0,
                                              +phase[row*nbase+b0]);
            if (b1 < nbase-1) cpl_matrix_set (model_matrix, row*nclo+clo, nclo+b1,
                                              +phase[row*nbase+b1]);
            if (b2 < nbase-1) cpl_matrix_set (model_matrix, row*nclo+clo, nclo+b2,
                                              -phase[row*nbase+b2]);
            CPLCHECK_MSG ("Cannot fill phase in model_matrix");
        }
    } /* End loop on clo and rows */

    /* Solve the system */
    cpl_matrix * res_matrix = cpl_matrix_solve_normal (model_matrix, rhs_matrix);
    CPLCHECK_MSG ("Cannot solve system");

    /* Compute residuals */
    cpl_matrix * residual_matrix = cpl_matrix_product_create (model_matrix, res_matrix);
    cpl_matrix_subtract (residual_matrix, rhs_matrix);
	double rms_fit = cpl_matrix_get_stdev (residual_matrix);

    /* Correct baseline with (1-f). Last baseline is kept unchanged */
	for (int base = 0; base < nbase - 1; base ++) {
        double f = cpl_matrix_get (res_matrix, nclo + base, 0);
        cpl_msg_info (cpl_func,"correction factor f = 1 %+.20f", -1*f);
        for (cpl_size row = 0; row < nrow; row++) phase[row*nbase+base] *= 1 - f;
	}

    /* Dump the residual with their units */
    const char * unit = cpl_table_get_column_unit (phase_table, name);
    cpl_msg_info (cpl_func, "residual stdev = %.20g [%s]", rms_fit, unit);
    
    /* Delete matrix */
    FREE (cpl_matrix_delete, residual_matrix);
    FREE (cpl_matrix_delete, res_matrix);
    FREE (cpl_matrix_delete, model_matrix);
    FREE (cpl_matrix_delete, rhs_matrix);
    
    gravi_msg_function_exit(1);
	return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Compute the absolute scaling coefficients of SC and FT OPDs
 * 
 * @param ft_table    Input table at the sampling rate of FT
 * 				    
 * @return coef       A new cpl_vector with the coeficients and the 
 *                    chi2 of the fit {a,b,chi2}  in [m]
 *
 * The routine determine the absolute scaling of SC (a) and FT (b) OPDs by
 * solving the linear system:
 *
 *    2pi./LBD_MET * PHASE_MET_ijt = a.OPD_SC_ijt - b.OPD_FT_ijt + c_ij
 *
 * The routine discards samples that are outside the SC DITs (defined with
 * OPD_SC==0). Then it solves the systems and compute residuals.
 * 
 * Since the routine uses a single coefficient (a) for all 6 baselines, the
 * 6 corresponding OPDs shall form a system with constant closure phase
 * versus modulation.
 *
 * The ft_table shall contain a column OPD with the FT OPDs, a 
 * column OPD_SC with the corresponding SC OPD (correctly resampled),
 * and a column PHASE_MET_FC with the correscponding FC phase (also resampled).
 * It shall have NDIT_SC * NBASE rows.
 */
/*----------------------------------------------------------------------------*/

cpl_vector * gravi_opds_fit_opdmet (cpl_table * ft_table, double lbd_met)
{
    gravi_msg_function_start(1);
	cpl_ensure (ft_table, CPL_ERROR_NULL_INPUT, NULL);
    
	int nbase = 6;

	/* Get the number of acquisitions */
	cpl_size nrow = cpl_table_get_nrow (ft_table) / nbase;

    /* Get the pointer to data */
    double * opd_sc    = cpl_table_get_data_double (ft_table, "OPD_SC");
    double * opd_ft    = cpl_table_get_data_double (ft_table, "OPD");
    double * phase_met = cpl_table_get_data_double (ft_table, "PHASE_MET_FC");
    CPLCHECK_NUL ("Cannot get data");
    
    /* Number of valid rows */
    cpl_size nrow_valid = 0;
    for (cpl_size row = 0; row < nrow; row++) if (opd_sc[row*nbase] != 0) nrow_valid++;
    cpl_msg_info (cpl_func,"nrow_valid = %lld", nrow_valid);

    cpl_ensure (nrow_valid > 100, CPL_ERROR_ILLEGAL_INPUT, NULL);

    /* Model and right-hand-side for the lineary system
     * (unfilled matrix are 0.0) */
    cpl_matrix * rhs_matrix   = cpl_matrix_new (nrow_valid*nbase, 1);
    cpl_matrix * model_matrix = cpl_matrix_new (nrow_valid*nbase, nbase + 2);

    for (int base = 0; base < nbase; base++) {
        cpl_size row_valid = 0;
        
        for (cpl_size row=0; row<nrow; row++) {
            if (opd_sc[row*nbase] == 0) continue;
            
            int idv = row_valid * nbase + base;
            int id  = row * nbase + base;

            /* Fill the OPD metrology */
            cpl_matrix_set (rhs_matrix, idv, 0, phase_met[id] * lbd_met / CPL_MATH_2PI);
            CPLCHECK_NUL ("Cannot set OPD");

            /* Fill the model b and c */
            cpl_matrix_set (model_matrix, idv, 0,  1*opd_sc[id]);
            cpl_matrix_set (model_matrix, idv, 1, -1*opd_ft[id]);
            CPLCHECK_NUL ("Cannot set SC or FT");

            /* Fill the model Aij (unfilled matrix are 0.0) */
            cpl_matrix_set (model_matrix, idv, 2 + base, 1.0);
            CPLCHECK_NUL ("Cannot set the zero points");

            row_valid++;
        } /* End loop rows */
    } /* End loop on base*/
    
    /* Solve the system */
    cpl_matrix * res_matrix = cpl_matrix_solve_normal (model_matrix, rhs_matrix);

    /* Compute residuals */
    cpl_matrix * residual_matrix = cpl_matrix_product_create (model_matrix, res_matrix);
    cpl_matrix_subtract (residual_matrix, rhs_matrix);
	double rms_fit = cpl_matrix_get_stdev (residual_matrix);

    /* Verbose */
    cpl_msg_info (cpl_func, "coeff SC = %.20g ", cpl_matrix_get (res_matrix, 0, 0));
    cpl_msg_info (cpl_func, "coeff FT = %.20g ", cpl_matrix_get (res_matrix, 1, 0));
    cpl_msg_info (cpl_func, "residual stdev = %.20g [m]", rms_fit);
    
    /* Set the Results */
    cpl_vector * opd_coeff = cpl_vector_new(3);
    cpl_vector_set (opd_coeff, GRAVI_SC, cpl_matrix_get (res_matrix, 0, 0));
    cpl_vector_set (opd_coeff, GRAVI_FT, cpl_matrix_get (res_matrix, 1, 0));
    cpl_vector_set (opd_coeff, 2, rms_fit);
    
    /* Delete matrix */
    FREE (cpl_matrix_delete, residual_matrix);
    FREE (cpl_matrix_delete, res_matrix);
    FREE (cpl_matrix_delete, model_matrix);
    FREE (cpl_matrix_delete, rhs_matrix);

    gravi_msg_function_exit(1);
	return opd_coeff;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Compute a guess of the OPD modulation of SC from FT and MET
 * 
 * @param spectrumsc_table      Input SC spectrum table (for DITs)
 * @param ft_table              The opd from FT
 * @param vismet_table          The retuced MET signal
 *
 * For each baseline and each row of spectrum, the routine compute 
 *
 *    guess = OPD_FT  -  LBD_MET * PHASE_MET / 2pi [m]
 *
 * The result is saved in the column OPD of the output table, which has
 * NDIT_SC * NBASE rows. Note that the FT and FT are averaged over each SC DIT.
 *
 * All tables shall contain a TIME columns in integer [us]. The ft_table
 * shall have a OPD column [m] and be of size NDIT_FT * NBASE. The vismet_table
 * shall have a PHASE_FC column [rad] and be of size NDIT_MET * NTEL.
 *
 * FIXME: this routine can be replaced by manipulation of columns
 *        created with the gravi_signal package.
 */
/*----------------------------------------------------------------------------*/

cpl_table * gravi_opds_compute_guess (cpl_table * spectrumsc_table,
                                      cpl_table * ft_table,
                                      cpl_table * vismet_table,
                                      double dit_sc,
                                      double lbd_met)
{
    gravi_msg_function_start(1);
    cpl_ensure (spectrumsc_table, CPL_ERROR_NULL_INPUT, NULL);
    cpl_ensure (ft_table,         CPL_ERROR_NULL_INPUT, NULL);
    cpl_ensure (vismet_table,        CPL_ERROR_NULL_INPUT, NULL);
    cpl_ensure (dit_sc>0,         CPL_ERROR_ILLEGAL_INPUT, NULL);

    int nbase = 6, ntel = 4;

    cpl_size nrow = cpl_table_get_nrow (spectrumsc_table);
    cpl_size nrow_met = cpl_table_get_nrow (vismet_table) / ntel;
    cpl_size nrow_ft  = cpl_table_get_nrow (ft_table) / nbase;
    int * time_SC  = cpl_table_get_data_int (spectrumsc_table, "TIME");
    CPLCHECK_NUL ("Cannot get data");

    /* Pointer on MET data */
    int * time_MET     = cpl_table_get_data_int (vismet_table, "TIME");
    double * phase_MET = cpl_table_get_data_double (vismet_table, "PHASE_FC");
    CPLCHECK_NUL ("Cannot get MET data");

    /* Pointer on OPD data */
    int * time_FT   = cpl_table_get_data_int (ft_table, "TIME");
    double * opd_FT = cpl_table_get_data_double (ft_table, "OPD");
    CPLCHECK_NUL ("Cannot get FT data");
    
    /* Create table */
    cpl_table * guess_table = cpl_table_new (nrow * nbase);    
    gravi_table_new_column (guess_table, "OPD", "m", CPL_TYPE_DOUBLE);

    /* Loop on base */
    for (int base = 0; base < nbase; base++) {
        int tel0 = GRAVI_BASE_TEL[base][0];
        int tel1 = GRAVI_BASE_TEL[base][1];
        
        /* Loop on the SC frames */
        cpl_size row_met = 0, row_ft = 0;
        for (cpl_size row = 0 ; row < nrow ; row++) {
                    
            /*
             * Average the MET 
             */
            int counter_met = 0;
            double opd_met = 0.0;
            while   (time_MET[row_met*ntel] < (time_SC[row] + dit_sc/2.)) {
                if ((time_MET[row_met*ntel] > (time_SC[row] - dit_sc/2.)) && (row_met < nrow_met)) {
                    opd_met += phase_MET[row_met*ntel+tel1] - phase_MET[row_met*ntel+tel0];
                    counter_met ++;
                }
                
                /* If not enough data to synchronize */
                if (row_met > nrow_met - 2) {
                    cpl_msg_warning (cpl_func,"Not enough data to synchronise the MET with SC");
                    break;
                }
                row_met ++;
            }
            /* End sum the MET over the current SC DIT */
                    
            opd_met = opd_met * lbd_met / CPL_MATH_2PI / counter_met; // [m]
                    
            /*
             * Average the FT
             */
            int counter_ft = 0;
            double opd_ft = 0.0;
            while   (time_FT[row_ft*nbase+base] < (time_SC[row] + dit_sc/2.)) {
                if ((time_FT[row_ft*nbase+base] > (time_SC[row] - dit_sc/2.)) && (row_ft < nrow_ft)) {
                    opd_ft += opd_FT[row_ft*nbase+base];
                    counter_ft ++;
                }
                
                /* If not enough data to synchronize */
                if (row_ft > nrow_ft - 2) {
                    cpl_msg_warning (cpl_func,"Not enough data to synchronise the FT with SC");
                    break;
                }
                row_ft ++;
            }
            /* End sum the FT over the current SC DIT */
            
            opd_ft = opd_ft / counter_ft;
                    
            /* Set the total guess OPD as OPD_FT - OPD_MET */
            cpl_table_set (guess_table, "OPD", row*nbase+base, opd_ft - opd_met);
            CPLCHECK_NUL ("Cannot compute opd guess");
        } /* End loop on row */
    } /* End loop on base */

    /* 
     * Remove the mean of each base
     */
    cpl_msg_info (cpl_func, "Remove the mean from opdguess table");
    
    for (int base = 0; base < 6; base++) {
        double mean = gravi_table_get_column_mean (guess_table, "OPD", base, nbase);
        gravi_table_add_scalar (guess_table, "OPD", base, nbase, -1*mean);
    }
    
    gravi_msg_function_exit(1);
	return guess_table;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Compute the mean opd of each baseline from spectrum
 * 
 * @param spectrum_table   Input spectrum table
 * @param detector_table   Corresponding IMAGING_DETECTOR
 * @param guess_table      Table with a guess of OPD [m], or NULL
 *
 * For each baseline, the routine compute the mean opd over the spectral
 * channels and polarisation using the routine gravi_ellipse_meanopd_create.
 * 
 * The results are saved in the column OPD [m] of the output table, which
 * has NDIT * NBASE rows, where NDIT in the size of the input spectrum table.
 * The output table also contains a TIME columns in integer [us].
 *
 * The guess_table shall be of size NDIT * NBASE and contains a column
 * OPD [m], filled with expected opd modulation in the spectrum table.
 * If provided, this tables is used to unwrap the computed phase before
 * averaging the spectral channels.
 *
 * Note that the routine uses a static wavelength calibration to compute the
 * OPDs, hence the resulting value are to-a-scaling-factor.
 */
/*----------------------------------------------------------------------------*/

cpl_table * gravi_opds_calibration (cpl_table * spectrum_table,
                                    cpl_table * detector_table,
                                    cpl_table * guess_table)
{
	int nbase = 6;
    
	/* Verbose */
	gravi_msg_function_start(1);
	cpl_ensure (spectrum_table, CPL_ERROR_NULL_INPUT, NULL);
	cpl_ensure (detector_table, CPL_ERROR_NULL_INPUT, NULL);

    /* Verbose the phase algorithm selected */
    if (USE_LINEAR_ELLIPSE) {
        cpl_msg_info (cpl_func, "Solve linearly the ellipse X,Y");
    }
    else {
        cpl_msg_info (cpl_func, "Solve non-linearly the ellipse X,Y");
    }
    
    /* Get the size of the vectors */
    cpl_size nrow  = cpl_table_get_nrow (spectrum_table);
    cpl_size nwave = gravi_spectrum_get_nwave (spectrum_table);
    cpl_size npol  = gravi_spectrum_get_npol (spectrum_table);
    CPLCHECK_NUL ("Cannot get size");

    /* Create a fake OI_WAVE table, to average the different
     * spectral channels [m] */
    
    cpl_table ** oiwave_tables = cpl_calloc (npol, sizeof (cpl_table *));
    for (int pol = 0; pol < npol; pol++) {
        oiwave_tables[pol] = cpl_table_new (nwave);
        cpl_table_new_column (oiwave_tables[pol], "EFF_WAVE", CPL_TYPE_DOUBLE);
        for (cpl_size wave = 0; wave < nwave; wave++) {
            double value = 1.97e-6 + (2.48e-6 - 1.97e-6) / nwave * wave;
            cpl_table_set (oiwave_tables[pol], "EFF_WAVE", wave, value);
        }
    }

    /* Create column in phase table */
    cpl_table * output_table = cpl_table_new (nbase * nrow);
    gravi_table_new_column (output_table,"OPD", "m", CPL_TYPE_DOUBLE);
    gravi_table_new_column (output_table,"TIME", "us", CPL_TYPE_INT);

    /* Set the time */
    for (cpl_size row = 0; row < nrow; row++) {
        double value = cpl_table_get (spectrum_table, "TIME", row, NULL);
        for (int base = 0; base < nbase; base++)
            cpl_table_set (output_table, "TIME", row*nbase+base, value);
    }

	/* Loop on base */
	for (int base = 0; base < 6; base ++) {
        
        /* Check if a guess exists */
        cpl_vector * opd_guess = NULL;
        if (guess_table) {
            opd_guess = cpl_vector_new (nrow);
            for (cpl_size row = 0; row < nrow; row++) {
                double value = cpl_table_get (guess_table, "OPD", row*nbase+base, NULL);
                cpl_vector_set (opd_guess, row, value);
            }
        }

        /* Recover the mean OPD modulation (averaved over pol
         * and channels) using the ellipses. In [m]. However since oiwave
         * is not known... this is to a scaling factor correction. */
        cpl_vector * mean_opd;
        mean_opd = gravi_ellipse_meanopd_create (spectrum_table, detector_table,
                                                 oiwave_tables, opd_guess, base);
        CPLCHECK_NUL ("Cannot compute opd");
        
        /* Save the mean opd [m] */
        for (cpl_size row = 0; row < nrow; row++) {
            double value = cpl_vector_get (mean_opd, row);
            cpl_table_set (output_table, "OPD", row*nbase+base, value);
        }
        CPLCHECK_NUL ("Cannot set phase");

		FREE (cpl_vector_delete, mean_opd);
        FREE (cpl_vector_delete, opd_guess);

	} /* End loop on base */

    FREELOOP (cpl_table_delete, oiwave_tables, npol);

	/* Verbose */
	gravi_msg_function_exit(1);
	return output_table;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Recover the OPD modulation from a spectrum_data and vismet_table
 *
 * @param spectrum_table    The input spectrum, with region DATA#
 * @param met_table         The corresponding METROLOGY table.
 * @param plist             The plist where to dump QC parameters
 *
 * @return Two allocated table containing the OPDs of SC and FT [m]
 *
 * The routine first reduces the METROLOGY table into metrology OPDs.
 * Then it computes the mean opd of each baseline for the FT, averaged
 * over spectral channels and polarisation, by fitting the ellipses.
 * 
 * It uses these OPDs from MET and FT to guess the opd modulation
 * on the SC, in order to unwrap it properly. Then it computes the mean
 * OPD of each baseline of the SC (averaged over spectral channels),
 * using this guess, and fitting the ellipses.
 * 
 * The routine then solves the linear system:
 *
 *     \f$\lambda_{MET}.OPD_{MET}  / 2\pi = a.OPD_{SC} - b.OPD_{FT} + c\f$
 * 
 * in order to determine the true scaling of the SC (a)
 * and of the SC (b). A special care is taken to only consider the
 * samples inside the SC DIT, to have a proper phase relation.
 *
 * The return tables are of size NDIT_SC * NBASE, and NDIT_SC * NBASE,
 * where NDIT in the size of the input spectrum tables. They contain a TIME
 * column [us] and a OPD table [m]. These tables are added to the input
 * spectrum_data with EXTNAME OPD_FT and OPD_SC.
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_wave_compute_opds (gravi_data * spectrum_data,
                                        cpl_table  * met_table,
                                        enum gravi_detector_type det_type)
{
	gravi_msg_function_start(1);
	cpl_ensure_code (spectrum_data, CPL_ERROR_NULL_INPUT);
	cpl_ensure_code (met_table,     CPL_ERROR_NULL_INPUT);

    /* Get DIT */
	cpl_propertylist * spectrum_header = gravi_data_get_header (spectrum_data);
    double dit_sc = gravi_pfits_get_dit_sc (spectrum_header)*1e6; // [us]
    CPLCHECK_MSG ("Cannot get DIT");
    
    /* Get the input */
    cpl_table * spectrumft_table = gravi_data_get_spectrum_data (spectrum_data, GRAVI_FT);
    cpl_table * detectorft_table = gravi_data_get_imaging_detector (spectrum_data, GRAVI_FT);
    cpl_table * spectrumsc_table = gravi_data_get_spectrum_data (spectrum_data, GRAVI_SC);
    cpl_table * detectorsc_table = gravi_data_get_imaging_detector (spectrum_data, GRAVI_SC);
    CPLCHECK_MSG ("Cannot get data");

    
	/*
	 * Reduce the raw METROLOGY into VIS_MET
	 */
	cpl_msg_info (cpl_func, "Compute the phase of MET_FC");
    
    cpl_table * vismet_table = gravi_metrology_create (met_table, spectrum_header);
    gravi_metrology_drs (met_table, vismet_table, spectrum_header);
    
    /* Compute the mean LBD_MET for this file */
    double lbd_met = gravi_pfits_get_met_wavelength_mean (spectrum_header, met_table);

	/*
	 * Compute the phase of SC and FT. For the SC, we use a guess 
     * of the opd modulation, computed from FT and MET, to unwrap.
	 */
	cpl_msg_info (cpl_func, "Compute OPD of FT from ellipse");
    
	cpl_table * ft_table;
	ft_table = gravi_opds_calibration (spectrumft_table, detectorft_table, NULL);
    gravi_opds_correct_closures (ft_table, "OPD");

	cpl_msg_info (cpl_func, "Compute OPD of SC from ellipse");
    
	cpl_table * guess_table;
    guess_table = gravi_opds_compute_guess (spectrumsc_table, ft_table, vismet_table, dit_sc, lbd_met);
    
	cpl_table * sc_table;
	sc_table  = gravi_opds_calibration (spectrumsc_table, detectorsc_table, guess_table);
    gravi_opds_correct_closures (sc_table, "OPD");
    FREE (cpl_table_delete, guess_table);

	CPLCHECK_MSG ("Cannot calibrate phase");

    /*
     * Interpolate MET and SC at the sampling frequency of FT
     * Results are saved in the ft_table table
     */
	cpl_msg_info (cpl_func, "Fit MET = a.SC - b.FT + c to get absolute modulation");
    
    gravi_vis_create_met_ft (ft_table, vismet_table);
    gravi_vis_create_opdsc_ft (ft_table, sc_table, dit_sc);
    
	CPLCHECK_MSG ("Cannot resample SC or MET at the FT frequency");

    /* 
     * Compute the scaling coefficients of OPDs by fitting:
     * OPD_MET_ijt = a.OPD_SC_ijt - b.OPD_FT_ijt + c_ij
     */
    cpl_vector * coeff_vector = gravi_opds_fit_opdmet (ft_table, lbd_met);
    
	CPLCHECK_MSG ("Cannot fit opdmet");
    
	/* Set the CHI2 of the fit in the QC parameters */
	cpl_propertylist_update_double (spectrum_header, QC_PHASECHI2,
                                    cpl_vector_get (coeff_vector, 2));
	cpl_propertylist_set_comment (spectrum_header, QC_PHASECHI2,
                                  "chi2 of a.SC-b.FT+c=MET");

    /* Add the OPD COEFF in QC parameters */
    for (int type_data = 0; type_data < 2; type_data++) {
        double tmp = cpl_vector_get (coeff_vector, type_data);
        cpl_propertylist_update_float (spectrum_header, OPD_COEFF_SIGN(type_data), tmp);
        cpl_propertylist_set_comment (spectrum_header, OPD_COEFF_SIGN(type_data), "wavelength correction");
    }

    /* Set a warning */
    if (cpl_vector_get (coeff_vector, 2) > 1e-7) {
        gravi_pfits_add_check (spectrum_header,"Residu of fit MET=a.SC-b.FT+c are high");
        
		cpl_msg_info (cpl_func,    "*************************************************");
		cpl_msg_warning (cpl_func, "****  !!! residuals of the fit too high !!!  ****");
		cpl_msg_warning (cpl_func, "****     SC and RMN may be desynchronized    ****");
		cpl_msg_warning (cpl_func, "****     (or out of the envelope in LOW)     ****");
		cpl_msg_info (cpl_func,    "*************************************************");
    }
    
    CPLCHECK_MSG ("Cannot set QC parameter");
    
    /*
     * Correct opd from overall scaling by fitting MET
     */
    double coeff_sc = cpl_vector_get (coeff_vector, GRAVI_SC);
    cpl_table_multiply_scalar (sc_table, "OPD", coeff_sc);
    
    double coeff_ft = cpl_vector_get (coeff_vector, GRAVI_FT);
    cpl_table_multiply_scalar (ft_table, "OPD", coeff_ft);

    FREE (cpl_vector_delete, coeff_vector);
	CPLCHECK_MSG ("Cannot correct OPDs from scaling coefficients");
    
    /* 
     * Fill the output 
     */
    if ((det_type == GRAVI_DET_SC || det_type == GRAVI_DET_ALL))
    {
        gravi_data_add_table (spectrum_data, NULL, "OPD_SC", sc_table);
        CPLCHECK_MSG ("Cannot set OPD_SC table");
    }
    else
        FREE (cpl_table_delete, sc_table);
        

    if ((det_type == GRAVI_DET_FT || det_type == GRAVI_DET_ALL))
    {
        gravi_data_add_table (spectrum_data, NULL, "OPD_FT", ft_table);
        CPLCHECK_MSG ("Cannot set OPD_FT table");
    }
    else
        FREE (cpl_table_delete, ft_table);

    gravi_data_add_table (spectrum_data, NULL, GRAVI_OI_VIS_MET_EXT, vismet_table);
	CPLCHECK_MSG ("Cannot set OI_VIS_MET table");
    
	gravi_msg_function_exit(1);
    return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Compute the wavelength of each channel for each 6 baselines
 *
 * @param spectrum_table    The input spectrum, with region DATA#
 * @param detector_table    The corresponding detector table
 * @param opd_table         The corresponding OPD modulation [m]
 *
 * @return A allocated table WAVE_FIBRE  [m]
 * 
 * For each baseline/polar, the phase of each channel is computed by
 * fitting the ellipse, using the opd_table modulation to unwrap the results.
 * Then the effective wavelength of this channel is computed by fitting
 * the relation:
 *
 *     OPD = lbd . PHASE / 2pi
 * 
 * The results are stored into 6 (or 12) columns of the form DATA_12_S
 * in the returned table.
 *
 * The opd_table shall be of size nrow * nbase, where nrow is the size
 * of the spectrum_table (number of DITs). It shall contains a column
 * OPD [m], filled with the opd modulation in the spectrum table.
 */
/*----------------------------------------------------------------------------*/

cpl_table * gravi_wave_fibre (cpl_table * spectrum_table,
                              cpl_table * detector_table,
                              cpl_table * opd_table)
{
	gravi_msg_function_start(1);
    cpl_ensure (spectrum_table, CPL_ERROR_NULL_INPUT, NULL);
    cpl_ensure (detector_table, CPL_ERROR_NULL_INPUT, NULL);
    cpl_ensure (opd_table,      CPL_ERROR_NULL_INPUT, NULL);

    int nbase = 6;
    char name[100];
    
    /* Create the output table */
    cpl_table * wave_fibre = cpl_table_new (1);

    /* Get the number of wavelength, region, polarisation... */
    cpl_size nwave = cpl_table_get_column_depth (spectrum_table, "DATA1");
    cpl_size n_region = cpl_table_get_nrow (detector_table);
    cpl_size nrow = cpl_table_get_nrow (spectrum_table);

    int npol = (n_region > 24 ? 2 : 1);

    /*
     * Calibration of each polarization and base
     */

    for (int pol = 0; pol < npol; pol++) {
        for (int base = 0; base < nbase; base ++) {
            cpl_msg_info (cpl_func, "Compute wave fibre for pol %i over %i, base %i over %i",
                          pol+1, npol, base+1, nbase);

            /* Get the index of the ABCD. */
            int iA = gravi_get_region (detector_table, base, 'A', pol);
            int iB = gravi_get_region (detector_table, base, 'B', pol);
            int iC = gravi_get_region (detector_table, base, 'C', pol);
            int iD = gravi_get_region (detector_table, base, 'D', pol);
            if (iA<0 || iB<0 || iC<0 || iD<0){
                cpl_msg_warning (cpl_func, "Don't found the A, B, C or D !!!");
            }
            
            /* Sign of this baseline */
            double phi_sign = gravi_region_get_base_sign (detector_table, base);
            
            /* Get the OPD into various flarous */
            cpl_matrix * opd_matrix = cpl_matrix_new (1, nrow);
            cpl_vector * opd_vector = cpl_vector_new (nrow);
            for (cpl_size row = 0; row < nrow; row ++ ) {
                double value = cpl_table_get (opd_table, "OPD", row*nbase+base, NULL);
                cpl_matrix_set (opd_matrix, 0, row, value);
                cpl_vector_set (opd_vector, row, value);
            }
            CPLCHECK_NUL ("Cannot extract the OPD");

            /* Create array to fill the wavelenghts */
            cpl_array * wavelength = cpl_array_new (nwave, CPL_TYPE_DOUBLE);
            cpl_array * wavechi2   = cpl_array_new (nwave, CPL_TYPE_DOUBLE);
            cpl_array_fill_window (wavelength, 0, nwave, 0.0);
            cpl_array_fill_window (wavechi2,   0, nwave, 1e10);
				
            /* 
             * Loop on spectral channel 
             */
            for (cpl_size wave = 0; wave < nwave; wave++) {
                    
                cpl_vector * vector_T = NULL, * vector_X, * vector_Y;

                /* Define the vector_X = C - A */
                vector_X = gravi_table_get_vector (spectrum_table, wave, GRAVI_DATA[iC]);
                vector_T = gravi_table_get_vector (spectrum_table, wave, GRAVI_DATA[iA]);
                cpl_vector_subtract (vector_X, vector_T);
                FREE (cpl_vector_delete, vector_T);
                
                /* Define the vector_Y = D - B */
                vector_Y = gravi_table_get_vector (spectrum_table, wave, GRAVI_DATA[iD]);
                vector_T = gravi_table_get_vector (spectrum_table, wave, GRAVI_DATA[iB]);
                cpl_vector_subtract (vector_Y, vector_T);
                FREE (cpl_vector_delete, vector_T);
                
                CPLCHECK_NUL ("Cannot compute vector_X or Y");

                /* Compute envelope from OPD for this channel */
                cpl_vector * envelope_vector = gravi_compute_envelope (opd_vector, wave, nwave);
                                
                /* Compute the phase from the ellipse */
                cpl_vector * phase;
                phase = gravi_ellipse_phase_create (vector_X, vector_Y,
                                                    envelope_vector);
                FREE (cpl_vector_delete, vector_X);
                FREE (cpl_vector_delete, vector_Y);
                FREE (cpl_vector_delete, envelope_vector);
                
                /* If the computation of the ellipse fails, we continue with next wave */
                if (phase == NULL) {
                    cpl_msg_warning (cpl_func, "Cannot compute wave for %lld and base %d", wave, base);
                    continue;
                }
                
                /* Multiply by the sign */    
                cpl_vector_multiply_scalar (phase, phi_sign);

                /* Unwrap phase from the OPD */
                double lbd_channel = 1.95e-6 + (2.46e-6 - 1.95e-6) / nwave * wave;
                gravi_vector_unwrap_with_guess (phase, opd_vector, CPL_MATH_2PI / lbd_channel);

                /* Fit the slope of the phase versus OPD gives the
                 * wavelength of the spectral element */
                const cpl_size mindeg = 0, maxdeg = 1;
                cpl_polynomial * fit_slope = cpl_polynomial_new (1);
                cpl_polynomial_fit (fit_slope, opd_matrix, NULL, phase, NULL, CPL_FALSE, &mindeg, &maxdeg);
                
                /* Compute residuals */
                cpl_vector * residuals = cpl_vector_new (nwave);
                double rechisq;
                cpl_vector_fill_polynomial_fit_residual	(residuals, phase, NULL, fit_slope, opd_matrix, &rechisq);
                
                if(PLOT_WAVE_PHASE_VS_OPD && (wave == WAVE_TO_PLOT))
                {
                    char gnuplot_str[200];
                    sprintf (gnuplot_str, "set title 'Wavelength (base %d)'; set xlabel 'Phase [rad]'; set ylabel 'OPD (m)';", base);
                    cpl_plot_vector (gnuplot_str, NULL, NULL, opd_vector);
                    sprintf (gnuplot_str, "set title 'Wavelength residuals (base %d)'; set xlabel 'Phase [rad]'; set ylabel 'OPD (m)';", base);
                    cpl_plot_vector (gnuplot_str, NULL, NULL, residuals);
                    CPLCHECK_NUL ("Cannot plot OPD versus phase");
                }

                CPLCHECK_NUL ("Cannot fit the OPD and phase");

                /* Get the slope */
                const cpl_size pow_slope = 1;
                double slope = cpl_polynomial_get_coeff (fit_slope, &pow_slope);

                /* Check slope sign. Should be positive */
                if (slope < 0.0 && wave == 0) {
                    cpl_msg_warning (cpl_func, "Negative wavelength!! "
                                     "Report to DRS team");
                }
					
                /* Set the wavelength */
                cpl_array_set (wavechi2, wave, sqrt(rechisq));
                cpl_array_set (wavelength, wave, CPL_MATH_2PI / fabs(slope));

                /* Free memory */
                cpl_vector_delete (phase);
                cpl_vector_delete (residuals);
                cpl_polynomial_delete (fit_slope);
                
            } /* End loop on wave */
				
            cpl_matrix_delete (opd_matrix);
            cpl_vector_delete (opd_vector);

            /* Set the wavelength array in the output table */
            sprintf (name, "BASE_%s_%s", GRAVI_BASE_NAME[base], GRAVI_POLAR(pol,npol));
            cpl_table_new_column_array (wave_fibre, name, CPL_TYPE_DOUBLE, nwave);
            cpl_table_set_column_unit (wave_fibre, name, "m");
            cpl_table_set_array (wave_fibre, name, 0, wavelength);
            cpl_array_delete (wavelength);

            /* Set the chi2 array in the output table */
            sprintf (name, "BASE_%s_%s_CHI2", GRAVI_BASE_NAME[base], GRAVI_POLAR(pol,npol));
            cpl_table_new_column_array (wave_fibre, name, CPL_TYPE_DOUBLE, nwave);
            cpl_table_set_array (wave_fibre, name, 0, wavechi2);
            cpl_array_delete (wavechi2);
				
        } /* End loop on base*/ 
    } /* End loop on polar */

	gravi_msg_function_exit(1);
	return wave_fibre;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Compute the WAVE_DATA table (1 wavelength per region) from
 *        the WAVE_FIBRE table (1 wavelength per baseline)
 * 
 * @param wavefibre_table     Input WAVE_FIBRE table
 * @param detector_table      Input IMAGING_DETECTOR table
 * @param fullstartx          FULLSTARTX of the spectrum
 *
 * @return The WAVE_DATA table, with interpolated DATA# column.
 *
 * The input wavefibre_table shall have columns of form BASE_12_S.
 * The rountine interpolate the wavelegength of these baseline as a
 * 2d map on pixel space, using the pixel position of each baseline
 * from detector_table. Then this polynomial fit it evaluated from
 * each region individually in order to build the DATA# column of
 * the wave_data output table.
 */
/*----------------------------------------------------------------------------*/

cpl_table * gravi_wave_fit_2d (cpl_table * wavefibre_table,
                               cpl_table * detector_table,
                               gravi_data * wave_param,
                               cpl_size fullstartx,
                               int spatial_order,
                               int spectral_order,
                               double * rms_residuals)
{
	gravi_msg_function_start(1);
	cpl_ensure (wavefibre_table, CPL_ERROR_NULL_INPUT, NULL);
	cpl_ensure (detector_table,  CPL_ERROR_NULL_INPUT, NULL);

	int nbase = 6;
	char name[100];
	*rms_residuals = 0;

    /* Get numbers */
    cpl_size n_region = cpl_table_get_nrow (detector_table);
    int npol = n_region > 24 ? 2 : 1;
    cpl_size nwave = cpl_table_get_column_depth (wavefibre_table, npol > 1 ? "BASE_12_S" : "BASE_12_C");
    CPLCHECK_NUL ("Cannot get data");
    
    /* get the wave param */
    cpl_propertylist * wave_param_plist = gravi_data_get_plist(wave_param,GRAVI_PRIMARY_HDR_EXT);

    /* Odd index, for SC only */
    cpl_vector * odd_index = cpl_vector_new (nwave);
    for (int i = fullstartx; i < fullstartx + nwave; i++)  {
        if (nwave > GRAVI_LBD_FTSC) cpl_vector_set (odd_index, i - fullstartx, ((i/64)%2 == 0) ? 0 : 1);
        else cpl_vector_set (odd_index, i - fullstartx, 0);
    }
    
    CPLCHECK_NUL ("Cannot buid odd_index");

    /* Save the 2D coeficient of each polar */
    cpl_polynomial ** coef_poly = cpl_calloc (npol, sizeof (cpl_polynomial*));
    
    /* Loop on polarisation */
    for (int pol = 0; pol < npol; pol++) {

		/* Prepare the 2D coordinates 
		 * and values to fit */
		cpl_vector * coord_X = cpl_vector_new (nbase * nwave);
		cpl_vector * coord_Y = cpl_vector_new (nbase * nwave);
        
		cpl_vector * all_wavelength = cpl_vector_new (nbase * nwave);
		cpl_vector * all_wavechi2 = cpl_vector_new (nbase * nwave);
		cpl_vector * all_valid = cpl_vector_new (nbase * nwave);
        
		/* 
		 * Loop on base and wave to get all
		 * wavelenght and coordinates 
		 */
		for (int base = 0; base < nbase; base ++) {
            
            /* Mean position of this baseline */
            int iA = gravi_get_region (detector_table, base, 'A', pol);
            int iB = gravi_get_region (detector_table, base, 'B', pol);
            int iC = gravi_get_region (detector_table, base, 'C', pol);
            int iD = gravi_get_region (detector_table, base, 'D', pol);
            
            /* WAVE_FIBRE data */
            sprintf (name, "BASE_%s_%s", GRAVI_BASE_NAME[base], GRAVI_POLAR(pol,npol));
            cpl_array * wavelength = cpl_table_get_data_array (wavefibre_table, name)[0];
            
            sprintf (name, "BASE_%s_%s_CHI2", GRAVI_BASE_NAME[base], GRAVI_POLAR(pol,npol));
            cpl_array * wavechi2 = cpl_table_get_data_array (wavefibre_table, name)[0];
            
            /* Loop on wave */
            for (cpl_size wave = 0; wave < nwave; wave++) {
                cpl_size pos = base * nwave + wave;
                
                /* Get the values */
                int nv = 0;
                double wave_value = cpl_array_get (wavelength, wave, &nv);
                double chi2_value = cpl_array_get (wavechi2, wave, &nv);

                /* The FT accept all channel for the 2D fit */
                cpl_vector_set (all_valid, pos, 1);

                if (nwave > 1000) {
                    /* SC HIGH */
                    if ((chi2_value > M_PI_4 ||
                         wave_value < gravi_pfits_get_double_default(wave_param_plist, "ESO OIWAVE HIGH LBD MIN", 2.01e-6) ||
                         wave_value > gravi_pfits_get_double_default(wave_param_plist, "ESO OIWAVE HIGH LBD MAX", 2.43e-6)))
                        cpl_vector_set (all_valid, pos, 0);
                }
                else if (nwave > 100) {
                    /* SC MEDIUM */
                    if ((chi2_value > M_PI_4 ||
                         wave_value < gravi_pfits_get_double_default(wave_param_plist, "ESO OIWAVE MED LBD MIN", 2.01e-6) ||
                         wave_value > gravi_pfits_get_double_default(wave_param_plist, "ESO OIWAVE MED LBD MAX", 2.43e-6)))
                        cpl_vector_set (all_valid, pos, 0);
                }
               else if (nwave > GRAVI_LBD_FTSC) {
                    /* SC LOW */
                    if ((chi2_value > M_PI_4 ||
                         //wave_value < 2.01e-6 ||
                         wave_value < gravi_pfits_get_double_default(wave_param_plist, "ESO OIWAVE LOW LBD MIN", 1.99e-6) ||
                         wave_value > gravi_pfits_get_double_default(wave_param_plist, "ESO OIWAVE LOW LBD MAX", 2.5e-6) ||
                         wave == 0 || wave == nwave-1))
                        cpl_vector_set (all_valid, pos, 0);
                }
                else if (nwave == GRAVI_LBD_FTSC) {
                    /* FT */
                    if ((chi2_value > M_PI_4 ||
                         //wave_value < 2.01e-6 ||
                         wave_value < 1.99e-6 ||
                         wave_value > 2.5e-6))
                        cpl_vector_set (all_valid, pos, 0);
                }
                
                /* Set the wavelength */
                cpl_vector_set (all_wavelength, pos, wave_value);
                cpl_vector_set (all_wavechi2,   pos, chi2_value);
                
                /* Set the X position as the mean of the 4 regions */
                cpl_vector_set (coord_X, pos, (double)(iA + iB + iC + iD) / 4.);
                
                /* Set the Y position. Add a shift of 0.15 pixels
                 * on odd output for SC detector */
                cpl_vector_set (coord_Y, pos, wave + cpl_vector_get (odd_index, wave)*0.15);
                
            } /* End loop on wave */
		} /* End loop on base */
        
		CPLCHECK_NUL ("Error get all wavelength");
		
		/* 
		 * Reformat the valid point in vector and matrix
		 */
		cpl_size nvalid = cpl_vector_get_sum (all_valid);
		
		cpl_msg_info (cpl_func, "Remove %lld/%lld bad wave",
					  nbase*nwave - nvalid, nbase * nwave);
		
		cpl_vector * vector = cpl_vector_new (nvalid);
		cpl_matrix * matrix = cpl_matrix_new (2, nvalid);
		
		for (cpl_size c = 0, i = 0 ; i < nwave * nbase; i ++) {
            if (!cpl_vector_get (all_valid, i)) continue;
            cpl_vector_set (vector, c, cpl_vector_get (all_wavelength, i));
            cpl_matrix_set (matrix, 0, c, cpl_vector_get (coord_X, i));
            cpl_matrix_set (matrix, 1, c, cpl_vector_get (coord_Y, i));
            c++;
		}
        
		CPLCHECK_NUL ("Error cropping all wavelength");
		
		/* 
		 * Perform a 2D fit with a polynomial model
		 * between the position and the wavelength = F(y, x)
		 */  
		cpl_size  deg2d[2] = {2, 3};
        if ( (nwave < 20) && (nwave > 8) ) {deg2d[0] = 2; deg2d[1] = 7;} // FIXME Useless line ?
        deg2d[0] = spatial_order;
        deg2d[1] = spectral_order;
        
		cpl_msg_info (cpl_func, "Fit a 2d polynomial {%lli..%lli} to the wavelengths map", deg2d[0], deg2d[1]);
        
		cpl_polynomial * fit2d = cpl_polynomial_new (2);
		cpl_polynomial_fit (fit2d, matrix, NULL, vector, NULL, CPL_TRUE, NULL, deg2d);
		coef_poly[pol] = fit2d;
        
		CPLCHECK_NUL ("Cannot fit 2D");
		
		/*
		 * Compute residuals
		 */
		double rechisq = 0.0;
		cpl_vector * residuals = cpl_vector_new (nvalid);
		cpl_vector_fill_polynomial_fit_residual	(residuals, vector, NULL, fit2d, matrix, &rechisq);
		*rms_residuals += cpl_vector_get_stdev(residuals)/npol;
        FREE (cpl_vector_delete, residuals);
		CPLCHECK_NUL ("Cannot compute residuals");

        
		FREE (cpl_matrix_delete, matrix);
		FREE (cpl_vector_delete, vector);
		FREE (cpl_vector_delete, all_wavelength);
		FREE (cpl_vector_delete, all_wavechi2);
		FREE (cpl_vector_delete, all_valid);
        FREE (cpl_vector_delete, coord_X);
        FREE (cpl_vector_delete, coord_Y);
        
    }	/* End loop on polarisation */
    
	
    /*
     * Create and fill the interpolated WAVE_DATA table
     */
    
    cpl_table * wavedata_table = cpl_table_new (1);
    cpl_vector * pos = cpl_vector_new (2);
    cpl_array * value = cpl_array_new (nwave, CPL_TYPE_DOUBLE);
	
    for (cpl_size region = 0 ; region < n_region; region ++) {
        
		int pol = gravi_region_get_pol (detector_table, region);
		/* Loop on wave to evaluate the 2D polynome */
		for (cpl_size wave = 0; wave < nwave; wave ++) {
            
            /* Evaluate */
            cpl_vector_set (pos, 0, region);
            cpl_vector_set (pos, 1, wave + cpl_vector_get (odd_index, wave)*0.15);
            
            double result = cpl_polynomial_eval (coef_poly[pol], pos);
            cpl_array_set (value, wave, result);

            }
		/* ensure cresent wavelength */
		double previous_wave = cpl_array_get(value, nwave/2, NULL);
		for (cpl_size wave_loop = nwave/2 ; wave_loop >= 0 ; wave_loop --){
            if (previous_wave < cpl_array_get(value, wave_loop, NULL))
                cpl_array_set(value, wave_loop, previous_wave);
            else previous_wave = cpl_array_get(value, wave_loop, NULL);
		}

        previous_wave = cpl_array_get(value, nwave/2, NULL);
        for (cpl_size wave_loop = nwave/2 ; wave_loop < nwave ; wave_loop ++){
            if (previous_wave > cpl_array_get(value, wave_loop, NULL))
                cpl_array_set(value, wave_loop, previous_wave);
            else previous_wave = cpl_array_get(value, wave_loop, NULL);
        }

		/* Add column */
		char * data_x = GRAVI_DATA[region];
		cpl_table_new_column_array (wavedata_table, data_x, CPL_TYPE_DOUBLE, nwave);
		cpl_table_set_array (wavedata_table, data_x, 0, value);
		
    } /* End loop on regions */

    /* Delete allocations */
    FREE (cpl_vector_delete, pos);
    FREE (cpl_array_delete, value);
    FREELOOP (cpl_polynomial_delete, coef_poly, npol);
    FREE (cpl_vector_delete, odd_index);
	
	gravi_msg_function_exit(1);
	return wavedata_table;
}



/*----------------------------------------------------------------------------*/
/*
 * @brief TBD
 */
/*----------------------------------------------------------------------------*/

cpl_table * gravi_wave_fit_individual (cpl_table * wave_individual_table,
                                          cpl_table * weight_individual_table,
                                          cpl_table * wave_fitted_table,
                                              cpl_table * opd_table,
                                              cpl_table * spectrum_table,
                                              cpl_table * detector_table,
                                              cpl_size fullstartx,
                                       double n0, double n1, double n2,
                                       double * rms_residuals)
{
    
    gravi_msg_function_start(1);
    
    cpl_ensure (wave_individual_table, CPL_ERROR_NULL_INPUT, NULL);
    cpl_ensure (weight_individual_table,  CPL_ERROR_NULL_INPUT, NULL);
    cpl_ensure (wave_fitted_table, CPL_ERROR_NULL_INPUT, NULL);
    cpl_ensure (opd_table,  CPL_ERROR_NULL_INPUT, NULL);
    cpl_ensure (spectrum_table,  CPL_ERROR_NULL_INPUT, NULL);
    cpl_ensure (detector_table,  CPL_ERROR_NULL_INPUT, NULL);
    
    /* Get the number of wavelength, region, polarisation... */
    cpl_size nwave = cpl_table_get_column_depth (spectrum_table, "DATA1");
    cpl_size n_region = cpl_table_get_nrow (detector_table);
    cpl_size nrow = cpl_table_get_nrow (spectrum_table);
    int npol = (n_region > 24 ? 2 : 1);
    int nbase = 6;
    cpl_size nwave_ref=3000;
    if (nwave<10) nwave_ref=600;
    
    CPLCHECK_NUL ("Cannot buid odd_index");
    
    /* Get OPD Table */
    
    
    
    /* create temporary variables */
    
    cpl_array * wave_individual_array = cpl_array_new (nwave, CPL_TYPE_DOUBLE);
    cpl_array * weight_individual_array = cpl_array_new (nwave, CPL_TYPE_DOUBLE);
    cpl_matrix * data_flux_matrix = cpl_matrix_new (nrow, nwave);
    cpl_matrix * vis_to_flux_matrix = cpl_matrix_new (nrow, 3);
    cpl_matrix * signal_matrix = cpl_matrix_new (nwave, nwave_ref);
    cpl_matrix * residual_matrix = cpl_matrix_new (nwave, nwave_ref);
    cpl_array * wave_reference_array = cpl_array_new (nwave_ref, CPL_TYPE_DOUBLE);
    
    /*initialise arrays */
    
    cpl_matrix_fill_column(vis_to_flux_matrix,1,2);
    for (cpl_size wave_ref = 0; wave_ref < nwave_ref; wave_ref+=1)
    {
        // make matrix todo
        double wave_value=1.95e-6+wave_ref*0.6e-6/((double) nwave_ref);
        cpl_array_set_double(wave_reference_array,wave_ref,wave_value);
    }
    
    CPLCHECK_NUL ("Cannot initialize arrays for wavelength fit");
    
    
    for (cpl_size region = 0 ; region < n_region; region ++) {
        
        int base=gravi_region_get_base (detector_table, region);
        char * data_x = GRAVI_DATA[region];
    
        cpl_msg_info_overwritable (cpl_func, "Least square fitting of wavelength for region %s", data_x);
        // get data_flux_matrix
        
        for (cpl_size row = 0; row < nrow; row ++ ) {
            cpl_array * flux_array= cpl_table_get_data_array(spectrum_table,data_x)[row];
            for (cpl_size wave = 0; wave < nwave; wave ++) {
                cpl_matrix_set (data_flux_matrix, row, wave, cpl_array_get(flux_array,wave, NULL));
            }
        }
        
    
        for (cpl_size wave_ref = 0; wave_ref < nwave_ref; wave_ref+=1)
        {
            // make matrix
            double wave_value=cpl_array_get(wave_reference_array,wave_ref,NULL);
            
            for (cpl_size row = 0; row < nrow; row ++ ) {
                double opd = cpl_table_get (opd_table, "OPD", row*nbase+base, NULL);
                double coherence_loss=1;
                if (fabs(opd) > 1e-9)
                {
                    if (nwave <30)
                    {
                        coherence_loss=sin(opd*19500)/(opd*19500);
                    }
                }
                cpl_matrix_set(vis_to_flux_matrix,row,0,cos(opd*6.28318530718/wave_value)*coherence_loss);
                cpl_matrix_set(vis_to_flux_matrix,row,1,sin(opd*6.28318530718/wave_value)*coherence_loss);
            }
            
            cpl_matrix * coef_vis = cpl_matrix_solve_normal(vis_to_flux_matrix,data_flux_matrix); // coef_vis is 3xnwave
            cpl_matrix * data_flux_fit = cpl_matrix_product_create(vis_to_flux_matrix,coef_vis);
            cpl_matrix * residuals_fit = cpl_matrix_duplicate(data_flux_fit);
            cpl_matrix_subtract (residuals_fit,data_flux_matrix); //residuals_fit is nrowxnwave
            
            for (cpl_size wave = 0; wave < nwave; wave ++ ) {
                cpl_matrix * temp_matrix = cpl_matrix_extract_column (data_flux_fit, wave);
                cpl_matrix_set(signal_matrix,wave,wave_ref,cpl_matrix_get_stdev(temp_matrix));
                FREE (cpl_matrix_delete, temp_matrix);
                cpl_matrix * temp_matrix2 = cpl_matrix_extract_column (residuals_fit, wave);
                cpl_matrix_set(residual_matrix,wave,wave_ref,cpl_matrix_get_stdev(temp_matrix2));
                FREE (cpl_matrix_delete, temp_matrix2);
            }
            
            FREE (cpl_matrix_delete, coef_vis);
            FREE (cpl_matrix_delete, data_flux_fit);
            FREE (cpl_matrix_delete, residuals_fit);
            
            CPLCHECK_NUL ("Cannot do Matrix inversion to calculate optimum wavelength");
            
        }
        
        // get minimum chi2 and amplitude signal for each wave
        cpl_size wave_ref=1;
        cpl_size discarded=1;
        for (cpl_size wave = 0; wave < nwave; wave ++ ) {
            
            cpl_matrix * chi2_extract=cpl_matrix_extract_row(residual_matrix,wave);
            
            cpl_matrix_get_minpos(chi2_extract,&discarded, &wave_ref );
            
            double wave_value = cpl_array_get(wave_reference_array, wave_ref, NULL );
            double weight_value = cpl_matrix_get(signal_matrix, wave , wave_ref )/(0.1+cpl_matrix_get(residual_matrix, wave, wave_ref ));
            
            cpl_array_set_double (wave_individual_array, wave, wave_value);
            cpl_array_set_double (weight_individual_array, wave, weight_value);

            FREE (cpl_matrix_delete, chi2_extract);
            
        }
        
        /* Add column */
        cpl_table_new_column_array (wave_individual_table, data_x, CPL_TYPE_DOUBLE, nwave);
        cpl_table_set_array (wave_individual_table, data_x, 0, wave_individual_array);
        cpl_table_new_column_array (weight_individual_table, data_x, CPL_TYPE_DOUBLE, nwave);
        cpl_table_set_array (weight_individual_table, data_x, 0, weight_individual_array);
        cpl_table_new_column_array (wave_fitted_table, data_x, CPL_TYPE_DOUBLE, nwave);
        cpl_table_set_array (wave_fitted_table, data_x, 0, wave_individual_array);
    }
    
    CPLCHECK_NUL ("Cannot get individual wavelength for each pixel");
    
    FREE (cpl_array_delete   ,wave_individual_array);
    FREE (cpl_array_delete   ,weight_individual_array);
    FREE (cpl_array_delete   ,wave_reference_array);
    FREE (cpl_matrix_delete  ,data_flux_matrix);
    FREE (cpl_matrix_delete  ,vis_to_flux_matrix);
    FREE (cpl_matrix_delete  ,signal_matrix);
    FREE (cpl_matrix_delete  ,residual_matrix);
    
    cpl_msg_info (cpl_func, "Now fitting polynomials on wavelength channels");
    
    cpl_matrix * coef_to_wave = cpl_matrix_new (n_region / npol,5);
    cpl_matrix * coef_to_wave_weight = cpl_matrix_new (n_region / npol,n_region / npol);
    cpl_matrix * wavelength = cpl_matrix_new(n_region / npol,1);
    
    // set coordinates
    for (cpl_size region = 0 ; region < n_region/ npol; region ++)
    {
        double mean_region = region - (n_region/npol-1)*0.5;
        cpl_matrix_set (coef_to_wave, region, 0, 1);
        cpl_matrix_set (coef_to_wave, region, 1, mean_region);
        cpl_matrix_set (coef_to_wave, region, 2, mean_region*mean_region);
        cpl_matrix_set (coef_to_wave, region, 3, mean_region*mean_region*mean_region);
        cpl_matrix_set (coef_to_wave, region, 4, mean_region*mean_region*mean_region*mean_region);
    }
    
    for (int pol = 0; pol < npol; pol++) {
        
        cpl_msg_info (cpl_func, "Looping for polyfit now, with pol: %i",(int) pol);
        
        for (cpl_size wave = 0; wave < nwave; wave ++) {
            
            // get the data for a common wave row
            for (cpl_size region = 0 ; region < n_region/ npol; region ++) {
     
                // Get the pointers to the table arrays
                char * data_x = GRAVI_DATA[region*npol+pol];
                
                
                const cpl_array * wave_array = cpl_table_get_array (wave_individual_table, data_x, 0);
                const cpl_array * weight_array = cpl_table_get_array (weight_individual_table, data_x, 0);
                
                
                cpl_matrix_set (wavelength, region, 0, cpl_array_get(wave_array,wave,NULL));
                double weight_value=cpl_array_get(weight_array,wave,NULL);
                cpl_matrix_set (coef_to_wave_weight, region, region, weight_value*weight_value);
            
            }
            
            
            cpl_matrix * coef_to_wave2 = cpl_matrix_product_create(coef_to_wave_weight,coef_to_wave);
            cpl_matrix * wavelength2 = cpl_matrix_product_create(coef_to_wave_weight,wavelength);
            
            // Fit a second order polynomial
            cpl_matrix * coeff = cpl_matrix_solve_normal(coef_to_wave2, wavelength2); // 5 x 1
            cpl_matrix * wavelength_fitted = cpl_matrix_product_create(coef_to_wave, coeff); // n_region / npol x 1
            
            CPLCHECK_NUL ("Cannot do Matrix inversion to calculate optimum polynomial for wavelength");
            
            // Write result to new wave table
            for (cpl_size region = 0 ; region < n_region/ npol; region ++) {
                
                char * data_x = GRAVI_DATA[region*npol+pol];
                //cpl_msg_info (cpl_func, "Writing: region  %i",region);
                cpl_array * wave_array = cpl_table_get_data_array (wave_fitted_table, data_x)[0];
                
                cpl_array_set_double(wave_array,wave,cpl_matrix_get(wavelength_fitted,region,0));

            }
            
            FREE (cpl_matrix_delete  ,coef_to_wave2);
            FREE (cpl_matrix_delete  ,wavelength2);
            FREE (cpl_matrix_delete  ,coeff);
            FREE (cpl_matrix_delete  ,wavelength_fitted);
        }
 
    }
    FREE (cpl_matrix_delete  ,coef_to_wave);
    FREE (cpl_matrix_delete  ,coef_to_wave_weight);
    FREE (cpl_matrix_delete  ,wavelength);
    
    CPLCHECK_NUL ("Cannot fit individual wavelength with 3rd order polynomial");
    
    cpl_msg_info (cpl_func, "Correcting for wavelength error");
    
    
    /* Correct the computed wavelength from the dispersion */
    for (cpl_size region = 0 ; region < n_region; region ++)
    {
            char * data_x = GRAVI_DATA[region];
            cpl_array * wavelength = cpl_table_get_data_array (wave_fitted_table, data_x)[0];
            cpl_size nwave = cpl_array_get_size (wavelength);
            for (cpl_size wave = 0 ; wave < nwave ; wave ++ ) {
                
                double result = cpl_array_get (wavelength, wave, NULL);
                double d_met  = (result - LAMBDA_MET) / LAMBDA_MET;
                cpl_array_set (wavelength, wave, result * (n0 + n1*d_met + n2*d_met*d_met));
                
            }
            for (cpl_size wave = nwave/2 ; wave < nwave-1 ; wave ++ ) {
                double result = cpl_array_get (wavelength, wave, NULL);
                double result2 = cpl_array_get (wavelength, wave+1, NULL);
                if (result2<result+2e-10) {
                    result2=result+2e-10;
                    cpl_array_set (wavelength, wave+1, result2);
                }
            }
            for (cpl_size wave = nwave/2 ; wave > 0 ; wave -- ) {
                double result = cpl_array_get (wavelength, wave, NULL);
                double result2 = cpl_array_get (wavelength, wave-1, NULL);
                if (result2>result-2e-10) {
                        result2=result-2e-10;
                    cpl_array_set (wavelength, wave-1, result2);
                }
            }
        }
    
    CPLCHECK_NUL ("Error in correcting for dispersion");
    
    gravi_msg_function_exit(1);
    return wave_fitted_table;
}




/*----------------------------------------------------------------------------*/
/**
 * @brief Correct the WAVE_FIBRE table from a harcoded dispersion model.
 *
 * @param wave_fibre      Input table, modified inplace
 * @param n0, n1, n2      Dispersion parameter
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_wave_correct_dispersion (cpl_table * wave_fibre,
                                              double n0, double n1, double n2)
{
	gravi_msg_function_start(1);
	cpl_ensure_code (wave_fibre, CPL_ERROR_NULL_INPUT);
    
	int nbase = 6;
	char name[100];
	  
    /* Get the number of polarisation */
    int npol = cpl_table_has_column (wave_fibre, "BASE_12_P") ? 2 : 1;

    /* Loop on columns in the table */
    for (int pol = 0; pol < npol; pol ++) {
		for (int base = 0; base < nbase; base ++) {
            
            /* Get data of this region */
            sprintf (name, "BASE_%s_%s", GRAVI_BASE_NAME[base], GRAVI_POLAR(pol, npol));
            cpl_array * wavelength = cpl_table_get_data_array (wave_fibre, name)[0];
            CPLCHECK_MSG ("Cannot get data");
            
            /* Loop on wave */
            cpl_size nwave = cpl_array_get_size (wavelength);
            for (cpl_size wave = 0 ; wave < nwave ; wave ++ ) {
                
                /* Correct the computed wavelength from the dispersion */
                double result = cpl_array_get (wavelength, wave, NULL);
                double d_met  = (result - LAMBDA_MET) / LAMBDA_MET;
                cpl_array_set (wavelength, wave, result * (n0 + n1*d_met + n2*d_met*d_met));
                
            }
		} /* End loop on base */
    } /* End loop on pol */
    
	gravi_msg_function_exit(1);
	return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Create a OI_WAVELENGTH_CORR table with color corrected wavelength.
 *
 * @param vis_data      Input gravi_data, modified inplace
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_wave_correct_color (gravi_data * vis_data)
{
    gravi_msg_function_start(1);
    cpl_ensure_code (vis_data, CPL_ERROR_NULL_INPUT);

    cpl_propertylist * primary_header = gravi_data_get_header (vis_data);
    cpl_propertylist * oiwave_header = NULL;
    cpl_table * oiwave_table = NULL;

    /* For each type of data SC / FT */
    int ntype_data = 2;
    for (int type_data = 0; type_data < ntype_data ; type_data ++) {

        /* Loop on polarisation */
        int npol = gravi_pfits_get_pola_num (primary_header, type_data);
        for (int pol = 0 ; pol<npol ; pol++) {
            oiwave_table = cpl_table_duplicate ( gravi_data_get_oi_wave ( vis_data, type_data, pol, npol ) );
            oiwave_header = cpl_propertylist_duplicate ( gravi_data_get_oi_wave_plist ( vis_data, type_data, pol, npol ) );

            /* here you can do what you want on this duplicated oi_wave_table
             * to get OI_FLUX table :
             *     cpl_table * oiflux_table = gravi_data_get_oi_flux(vis_data, type_data, pol, npol)
             * */


            /* save the new extension with name OI_WAVELENGTH_CORR */
            gravi_data_add_table(vis_data, oiwave_header, "OI_WAVELENGTH_CORR", oiwave_table);

            CPLCHECK_MSG("Cannot apply color wave correction");
        }
        /* End loop on polarisation */
    }
    /* End loop on data_type */

    gravi_msg_function_exit(1);
    return CPL_ERROR_NONE;
}


/*----------------------------------------------------------------------------*/
/**
 * @brief Compute the (useless) TEST_WAVE table from the WAVE_FIBRE
 *        and the PROFILE maps. 
 * 
 * @param wavedata_table      Input WAVE_DATA table (with DATA# columns)
 * @param wavefibre_table     Input WAVE_FIBRE table (with BASE_# columns)
 * @param profile_table       Input profile table
 * @param detector_table      Input IMAGING_DETECTOR table
 *
 * @return A allocated imagelist with the test images.
 */
/*----------------------------------------------------------------------------*/

cpl_imagelist * gravi_wave_test_image (cpl_table * wavedata_table,
                                       cpl_table * wavefibre_table,
                                       cpl_table * profile_table,
                                       cpl_table * detector_table)
{
	gravi_msg_function_start(1);
	cpl_ensure (wavedata_table,   CPL_ERROR_NULL_INPUT, NULL);
	cpl_ensure (wavefibre_table,  CPL_ERROR_NULL_INPUT, NULL);
	cpl_ensure (detector_table,   CPL_ERROR_NULL_INPUT, NULL);
	cpl_ensure (profile_table,    CPL_ERROR_NULL_INPUT, NULL);

	int nv = 0;
	char name[100];

	cpl_size n_region = cpl_table_get_nrow (detector_table);
	int npol = (n_region > 24 ? 2 : 1);

	CPLCHECK_NUL ("Cannot get data");

	/* Init the output images */
	cpl_size sizex = cpl_table_get_column_dimension (profile_table, "DATA1", 0);
	cpl_size sizey = cpl_table_get_column_dimension (profile_table, "DATA1", 1);
	
	cpl_image * profilesum_image = cpl_image_new (sizex, sizey, CPL_TYPE_DOUBLE);
	cpl_image_fill_window (profilesum_image, 1, 1, sizex, sizey, 0.0);
	
	cpl_image * wave_image = cpl_image_new (sizex, sizey, CPL_TYPE_DOUBLE);
	cpl_image_fill_window (wave_image, 1, 1, sizex, sizey, 0.0);

	cpl_image * realwave_image = cpl_image_new (sizex, sizey, CPL_TYPE_DOUBLE);
	cpl_image_fill_window (realwave_image, 1, 1, sizex, sizey, 0.0);
	
	CPLCHECK_NUL ("Cannot prepare output image");

	/* Loop on regions */
	for (cpl_size region = 0 ; region < n_region; region ++) {

	  /* Load the profile of this region as an image */
      cpl_imagelist * profile_imglist = gravi_imagelist_wrap_column (profile_table, GRAVI_DATA[region]);
	  cpl_image * profile_image = cpl_imagelist_get (profile_imglist, 0);

	  CPLCHECK_NUL ("Cannot get data");

	  /* Sum all profils of all regions */
	  cpl_image_add (profilesum_image, profile_image);

	  /*
	   * Define an image of each region with its WAVE_DATA 
	   */
	  const cpl_array * wavelength;
	  wavelength = cpl_table_get_array (wavedata_table, GRAVI_DATA[region], 0);
	  CPLCHECK_NUL ("Cannot get data");
	  
	  for (cpl_size x = 0; x < sizex; x ++){
          for (cpl_size y = 0; y < sizey; y ++){
              if (cpl_image_get (profile_image, x+1, y+1, &nv) > 0.01)
                  cpl_image_set (wave_image, x+1, y+1,
                                 cpl_array_get (wavelength, x, NULL));
          }
	  }
	  CPLCHECK_NUL ("Cannot make image of wave_map");

	  /* 
	   * Define an image of each region with its WAVE_FIBER 
	   */
	  int base = gravi_region_get_base (detector_table, region);
	  int pol  = gravi_region_get_pol (detector_table, region);
	  sprintf (name, "BASE_%s_%s", GRAVI_BASE_NAME[base], GRAVI_POLAR(pol, npol));
	  wavelength = cpl_table_get_array (wavefibre_table, name, 0);
	  CPLCHECK_NUL ("Cannot get data");
	  
	  for (cpl_size x = 0; x < sizex; x ++){
          for (cpl_size y = 0; y < sizey; y ++){
              if (cpl_image_get (profile_image, x+1, y+1, &nv) > 0.01)
                  cpl_image_set (realwave_image, x+1, y+1,
                                 cpl_array_get (wavelength, x, NULL));
          }
	  }
	  
      gravi_imagelist_unwrap_images (profile_imglist);
	  CPLCHECK_NUL ("Cannot make image of wave_fibre");
	} /* End loop on regions */

    /* Create an imagelist */
	cpl_imagelist *  testwave_imglist = cpl_imagelist_new ();
	cpl_imagelist_set (testwave_imglist, wave_image, 0);
	cpl_imagelist_set (testwave_imglist, profilesum_image, 1);
	cpl_imagelist_set (testwave_imglist, realwave_image, 2);
   	  
	gravi_msg_function_exit(1);
	return  testwave_imglist;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Compute the QC parameters of the WAVE product
 *
 * @param wave_map     Input/output gravi_data update
 * @param profile_map  Input profile calibration
 *
 * Compute the QC parameters from the WAVE_DATA of each combiner (SC/FT). Also
 * computes a uselss image of the detector with the wavelength... for visual.
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_wave_qc (gravi_data * wave_map, gravi_data * profile_map)
{
	gravi_msg_function_start(1);
	cpl_ensure_code (wave_map, CPL_ERROR_NULL_INPUT);

	char name[100];

	cpl_propertylist * wave_header = gravi_data_get_header (wave_map);
	
	/* Loop on extensions (thus SC/FT) */
	for (int type_data = 0; type_data < 2; type_data ++ ) {

	  /* Get WAVE_FIBRE and IMAGING_DETECTOR tables */
	  cpl_table * detector_table = gravi_data_get_imaging_detector (wave_map, type_data);
	  cpl_table * wave_data = gravi_data_get_wave_data (wave_map, type_data);
	  cpl_size n_region = cpl_table_get_nrow (detector_table);

	  /* Get the full STARTX */
	  cpl_propertylist * plist = gravi_data_get_wave_data_plist (wave_map, type_data);
	  int fullstartx = gravi_pfits_get_fullstartx (plist);

	  CPLCHECK_MSG ("Cannot get data");

	  /* 
	   * Compute the min and max wave over all regions
	   */
	  double minwave = -1e10;
	  double maxwave =  1e10;
	  
	  for (int region = 0; region < n_region; region++) {
          const cpl_array * wavelength = cpl_table_get_array (wave_data, GRAVI_DATA[region], 0);
          
		minwave = CPL_MAX (minwave, cpl_array_get_min (wavelength));
		maxwave = CPL_MIN (maxwave, cpl_array_get_max (wavelength));
	  } /* End loop on regions */

      cpl_msg_info (cpl_func,"%s = %g [m]", QC_MINWAVE(type_data), minwave);
      cpl_msg_info (cpl_func,"%s = %g [m]", QC_MAXWAVE(type_data), maxwave);
	  cpl_propertylist_update_double (wave_header, QC_MINWAVE(type_data), minwave);
	  cpl_propertylist_update_double (wave_header, QC_MAXWAVE(type_data), maxwave);
	  
	  CPLCHECK_MSG ("Cannot compute minwave or maxwave");

	  
	  /* 
	   * Compute the pixel position of QC specified argon wavelength.
	   * Having two lines allow to check the position and dispersion 
	   */
	  int qc_reg[3] = {0,23,47};
	  double qc_wave[2] = {2.099184e-6, 2.313952e-6};

	  /* Loop on regions */
	  for (int reg = 0; reg < (n_region > 24 ? 3 : 2); reg++) {
		cpl_size region = qc_reg[reg];

		const cpl_array * wavelength = cpl_table_get_array (wave_data, GRAVI_DATA[region], 0);
		cpl_size nwave = cpl_array_get_size (wavelength);
		const double * wave_tab = cpl_array_get_data_double_const (wavelength);
		CPLCHECK_MSG ("Cannot get wavelength data");

		/* Loop on the two argon lines */
		for (int iqc = 0 ; iqc < 2 ; iqc++) {
		  cpl_size l2 = 0;
		  if ( wave_tab[0] < wave_tab[nwave-1]) {while (wave_tab[l2] < qc_wave[iqc]) l2 ++;}
		  else                          {while (wave_tab[l2] > qc_wave[iqc]) l2 ++;}

		  if (l2-1 < 0 || l2 > nwave-1) {
			cpl_msg_error (cpl_func, "Cannot find the QC position for lbd=%g", qc_wave[iqc]);
			continue;
		  }
		  
		  /* Position on full detector */
		  double qc_pos = 0.0;
		  qc_pos = fullstartx + (l2-1) + (qc_wave[iqc] - wave_tab[l2-1]) / (wave_tab[l2] - wave_tab[l2-1]);
		  
		  sprintf (name, "ESO QC REFWAVE%i", iqc+1);
		  cpl_propertylist_update_double (wave_header, name, qc_wave[iqc]);
		  cpl_propertylist_set_comment (wave_header, name, "[m] value of ref wave");
		  
		  sprintf (name, "ESO QC REFPOS%i %s%lli", iqc+1, GRAVI_TYPE(type_data),region+1);
		  cpl_propertylist_update_double (wave_header, name, qc_pos);
		  cpl_propertylist_set_comment (wave_header, name, "[pix] position of ref wave");
		  
		  cpl_msg_info (cpl_func, "%s = %f [pix] for %e [m]", name, qc_pos, qc_wave[iqc]);
		  
		  CPLCHECK_MSG ("Cannot set QC");
		}
		/* End loop on the 2 argon lines */
	  } /* End loop on 2 or 3 regions */
	  
	} /* End loop on SC / FT */

    
    /* 
     * Create the test image for SC (only used for debug)
     */
	cpl_table * wavedata_table = gravi_data_get_wave_data (wave_map, GRAVI_SC);
	cpl_table * wavefibre_table = gravi_data_get_wave_fibre (wave_map, GRAVI_SC);
	cpl_table * profile_table = gravi_data_get_table (profile_map, GRAVI_PROFILE_DATA_EXT);
	cpl_table * detector_table = gravi_data_get_imaging_detector (wave_map, GRAVI_SC);

    cpl_imagelist * testwave_imglist;
    testwave_imglist = gravi_wave_test_image (wavedata_table, wavefibre_table,
                                              profile_table, detector_table);
	CPLCHECK_MSG ("Cannot compute TEST_WAVE");

	/* Set TEST_WAVE in output */
	cpl_propertylist * plist = gravi_data_get_plist (profile_map, GRAVI_PROFILE_DATA_EXT);
	gravi_data_add_cube (wave_map, cpl_propertylist_duplicate (plist),
                         "TEST_WAVE", testwave_imglist);
    
	CPLCHECK_MSG ("Cannot set data");
	
	gravi_msg_function_exit(1);
	return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Create the WAVE calibration map. 
 *
 * @param wave_map          Output wave_map, already allocated.
 * @param spectrum_data     Input spectrum_data
 * @param type_data         GRAVI_SC or GRAVI_FT
 * @param parlist           Input parameter list with :
 *                          - force-wave-ft-equal : Force the spatial order of
 *                          the wavelength 2D fit for FT to zero (so all region
 *                          share the same calibration). This is used to build
 *                          the P2VM calibration of the TAC real-time code running
 *                          on the instrument ifself.
 *
 *
 * The output WAVE map is filled with WAVE_FIBRE and WAVE_DATA
 * tables, as well as QC parameters in the main header.
 * 
 * The input spectrum_data shall contain the IMAGING_DATA_FT, IMAGING_DETECTOR_FT
 * and OPD_FT tables if type_data is GRAVI_FT (respectively SC). The latter 
 * table corresponds to the OPD modulation computed with gravi_wave_compute_opds.
 *
 * The function first compute the WAVE_FIBRE table by fitting the modulation
 * for each baseline and spectral channel.
 *
 * Then the 6 baselines of WAVE_FIBRE are converted into a wavelength map for
 * each 24 output, by the mean of a spatial fit over the detector.
 */
/*----------------------------------------------------------------------------*/

cpl_error_code  gravi_compute_wave (gravi_data * wave_map,
                                    gravi_data * spectrum_data,
                                    int type_data, const cpl_parameterlist * parlist,
                                    gravi_data * wave_param)
{
	gravi_msg_function_start(1);
	cpl_ensure_code (wave_map,      CPL_ERROR_NULL_INPUT);
	cpl_ensure_code (spectrum_data, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code (type_data == GRAVI_SC || type_data == GRAVI_FT,
                     CPL_ERROR_ILLEGAL_INPUT);

    /* Get headers */
	cpl_propertylist * raw_header = gravi_data_get_header (spectrum_data);
    cpl_propertylist * wave_header = gravi_data_get_header (wave_map);

    /* Dump full header of wave data */
    if (type_data == GRAVI_FT) {
        cpl_propertylist_append (wave_header, raw_header);
    }

    
    /* Copy IMAGING_DETECTOR in output WAVE */
    const char * extname = (type_data == GRAVI_SC) ? GRAVI_IMAGING_DETECTOR_SC_EXT : GRAVI_IMAGING_DETECTOR_FT_EXT;
    gravi_data_copy_ext (wave_map, spectrum_data, extname);

    /* 
     * Compute WAVE_FIBRE map 
     */
    cpl_msg_info (cpl_func, "Compute WAVE_FIBRE for %s", GRAVI_TYPE(type_data));
        
    cpl_table * spectrum_table = gravi_data_get_spectrum_data (spectrum_data, type_data);
    cpl_table * detector_table = gravi_data_get_imaging_detector (spectrum_data, type_data);
    cpl_table * opd_table      = gravi_data_get_table (spectrum_data, type_data==GRAVI_SC ? "OPD_SC" : "OPD_FT");

    cpl_table * wavefibre_table;
    wavefibre_table = gravi_wave_fibre (spectrum_table, detector_table, opd_table);
    CPLCHECK_MSG ("Cannot compute the WAVE_FIBRE");

    /* 
     * Correct WAVE_FIBRE from the dispersion model 
     */
    cpl_msg_info (cpl_func, "Correct dispersion in WAVE_FIBRE of %s", GRAVI_TYPE(type_data));

	/* Hardcoded values to correct for the fiber dispersion */
	double n0 = 1.0, n1 = -0.0165448, n2 = 0.00256002;
	cpl_msg_info (cpl_func,"Rescale wavelengths with dispersion (%g,%g,%g)",n0,n1,n2);
	
	cpl_propertylist_update_string (wave_header, "ESO QC WAVE_CORR", "lbd*(N0+N1*(lbd-lbd0)/lbd0+N2*(lbd-lbd0)^2/lbd0^2)");
	cpl_propertylist_update_double (wave_header, "ESO QC WAVE_CORR N0", n0);
	cpl_propertylist_update_double (wave_header, "ESO QC WAVE_CORR N1", n1);
	cpl_propertylist_update_double (wave_header, "ESO QC WAVE_CORR N2", n2);
	CPLCHECK_MSG ("Cannot set Keywords");
        
    gravi_wave_correct_dispersion (wavefibre_table, n0, n1, n2);
    CPLCHECK_MSG ("Cannot correct dispersion in wave");
    
    /* 
     * Fit the 2d dispersion WAVE_FIBRE into WAVE_DATA
     */
    cpl_msg_info (cpl_func,"Fit WAVE_DATA map for %s", GRAVI_TYPE(type_data));
    
    /* Get the fullstartx */
    cpl_propertylist * spectrum_plist = gravi_data_get_spectrum_data_plist (spectrum_data, type_data);
    int fullstartx = gravi_pfits_get_fullstartx (spectrum_plist);

    /* Interpolate table 2D */
    cpl_table * wavedata_table;
    int spatial_order=2; // default spatial order
    int spectral_order=3; // default spectral order
    double rms_residuals;
    if (type_data == GRAVI_FT && gravi_param_get_bool(parlist, "gravity.calib.force-wave-ft-equal")) {
    	spatial_order = 0;
    	cpl_msg_info (cpl_func, "Option force-waveFT-equal applied");
    }
// Keep default value
//    if (type_data == GRAVI_SC) {
//        spectral_order = gravi_param_get_int(parlist, "gravity.calib.wave-spectral-order");
//        cpl_msg_info (cpl_func, "Option set_spectral order to %d", spectral_order);
//    }

    wavedata_table = gravi_wave_fit_2d (wavefibre_table,
                                        detector_table,
                                        wave_param,
                                        fullstartx, spatial_order, spectral_order,  &rms_residuals);
	cpl_propertylist_update_double (wave_header, QC_RMS_RESIDUALS(type_data), rms_residuals);
    
    
    double rms_fit=cpl_propertylist_get_double (raw_header, QC_PHASECHI2);
    cpl_propertylist_update_double (wave_header, QC_CHI2WAVE(type_data),rms_fit*1e9);
    cpl_propertylist_set_comment (wave_header, QC_CHI2WAVE(type_data),"[nm]rms a.SC-b.FT+c=MET");

    CPLCHECK_MSG ("Cannot fit 2d data");
    
    /*
     * New Wavelength interpolation made by sylvestre on January 30 2018 
     */
    
    if (type_data == GRAVI_SC && !strcmp (gravi_pfits_get_spec_res(raw_header), "LOW"))
    {
        cpl_msg_info (cpl_func, "Additional Wavelength Fit");
        cpl_table * wave_individual_table = cpl_table_new (1);
        cpl_table * weight_individual_table   = cpl_table_new (1);
        cpl_table * wave_fitted_table = cpl_table_new (1);
        
        gravi_wave_fit_individual (wave_individual_table,
                                   weight_individual_table,
                                   wave_fitted_table,
                                   opd_table,
                                   spectrum_table,
                                   detector_table,
                                   fullstartx,
                                   n0,n1,n2,
                                   &rms_residuals);
        
        cpl_msg_info (cpl_func,"Add tables in wave_map");

        gravi_data_add_table (wave_map, cpl_propertylist_duplicate (spectrum_plist),
                              "WAVE_INDIV_SC", wave_individual_table);
        gravi_data_add_table (wave_map, cpl_propertylist_duplicate (spectrum_plist),
                              "WAVE_WEIGHT_SC", weight_individual_table);
        gravi_data_add_table (wave_map, cpl_propertylist_duplicate (spectrum_plist),
                              "WAVE_FITTED_SC", wavedata_table);
        
        gravi_data_add_table (wave_map, cpl_propertylist_duplicate (spectrum_plist),
                              GRAVI_WAVE_FIBRE_EXT(type_data), wavefibre_table);
        gravi_data_add_table (wave_map, cpl_propertylist_duplicate (spectrum_plist),
                              GRAVI_WAVE_DATA_EXT(type_data), wave_fitted_table);
    } else {
        /*
         * Add the WAVE_FIBRE and WAVE_DATA table in the wave_map
         */
        cpl_msg_info (cpl_func,"Add WAVE_FIBRE and WAVE_DATA in wave_map");

        gravi_data_add_table (wave_map, cpl_propertylist_duplicate (spectrum_plist),
                              GRAVI_WAVE_FIBRE_EXT(type_data), wavefibre_table);

        gravi_data_add_table (wave_map, cpl_propertylist_duplicate (spectrum_plist),
                              GRAVI_WAVE_DATA_EXT(type_data), wavedata_table);

        }



    
    CPLCHECK_MSG ("Cannot set data");
    
	gravi_msg_function_exit(1);
	return CPL_ERROR_NONE;
}


/**@}*/
