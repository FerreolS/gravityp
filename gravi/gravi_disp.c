/* $Id: gravi_disp.c,v 1.10 2012/03/23 15:10:40 nazouaoui Exp $
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
 * @defgroup gravi_disp  Dispersion and Argon calibration
 *
 * This module contains two main functions : @c gravi_compute_argon_pos() and
 * @c gravi_compute_disp() which are called respectively by the recipe
 * @c gravity_wavelamp and @c gravity_disp.
 *
 * The @c gravi_compute_argon_pos computes the position of the known argon
 * line wavelengths on which the @c gravi_compute_disp function relies as absolute
 * wavelength calibration. Based on this, a precise dispersion calibration can be
 * done.
 *
 */
/**@{*/

/*
 * History
 * 07/12/208 add wave_param to  gravi_compute_argon_pos
 */
/*----------------------------------------------------------------------------
                                    DEBUG
 -----------------------------------------------------------------------------*/

#define INFO_DEBUG 0
#define GRAVI_ACOEFF_RANGE 0.02

/*-----------------------------------------------------------------------------
                                   Includes
 -----------------------------------------------------------------------------*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cpl.h>
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
#include "gravi_signal.h"

#include "gravi_vis.h"
#include "gravi_disp.h"

/*-----------------------------------------------------------------------------
                                 Private prototypes
 -----------------------------------------------------------------------------*/

cpl_table * gravi_fit_fddl_lin (cpl_table * oiflux_table);

cpl_table * gravi_fit_dispersion (cpl_table * oiflux_table,
                                  cpl_table * oivis_table,
                                  cpl_table * oiwave_table,
                                  double * GDrms,
                                  double * Amin,
                                  double * Amax);

/*-----------------------------------------------------------------------------
                             Functions code
 -----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
 * @brief Compute the DISP_MODEL calibration map
 *
 * @param vis_data  The input vis_data with several observations
 * 
 * @return a gravi_data with the DISP_MODEL table extention.
 * 
 * \exception CPL_ERROR_NULL_INPUT input data is missing
 *
 *
 * The recipe first compute the linearisation coefficients of the FDDL
 * and store them into column LIN_FDDL_SC and LIN_FDDL_FT. Then
 * it computes the dispersion index of the FDDL and store them into
 * columns BETA and GAMMA. This table is then stored as extention
 * DISP_MODEL into the returned, newly allocated, gravi_data.
 *
 * The input vis_data shall have at least 10 observations (60 rows
 * in OI_VIS tables).
 *
 */
/*----------------------------------------------------------------------------*/

gravi_data * gravi_compute_disp (gravi_data * vis_data)
{
    gravi_msg_function_start(1);
	cpl_ensure (vis_data, CPL_ERROR_NULL_INPUT, NULL);

    /* Get data */
    cpl_size ntel = 4;
    cpl_propertylist * vis_header = gravi_data_get_header (vis_data);
    cpl_size npol = gravi_pfits_get_pola_num (vis_header, GRAVI_SC);
	cpl_table * oiflux_table = gravi_data_get_oi_flux (vis_data, GRAVI_SC, 0, npol);
    CPLCHECK_NUL ("Cannot get data");

    /* Get the number of observations */
    cpl_size nrow = cpl_table_get_nrow (oiflux_table) / ntel;
    cpl_msg_info (cpl_func,"Input vis_data has %lld observation",nrow);

    /* Check the number of observation */
    if (nrow < 10) {
        cpl_error_set_message (cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                               "Not enough observations to compute"
                               "a dispersion model. Check input SOF");
        return NULL;
    }

    /* 
     * Create output data
     */
    gravi_data * disp_map = gravi_data_new (0);

    /* Set the input header */
    cpl_propertylist * disp_header = gravi_data_get_header (disp_map);
    cpl_propertylist_append (disp_header, vis_header);

    /* Set a QC parameter with the number of observations */
    const char * qc_name = "ESO QC DISP NEXP";
    cpl_propertylist_update_int (disp_header, qc_name, nrow);
    cpl_propertylist_set_comment (disp_header, qc_name, "Number of exposures used");
    

    /* 
     * Compute the coefficient of FDDL linearity
     */
    
    cpl_table * linearity_table;
    linearity_table = gravi_fit_fddl_lin (oiflux_table);
    CPLCHECK_NUL ("Cannot compute the FDDL linearity");

    
    /* 
     * Compute the coefficients for FDDL index dispersion
     */
    
    cpl_table * dispwave_table;
    double GDrms = 0.0, Amin = 1e4, Amax = -1e4;

    /* Loop on polarisations */
    for (int pol = 0; pol < npol; pol++) {

        /* Get table for this polarisation */
        cpl_table * oiflux_table = gravi_data_get_oi_flux (vis_data, GRAVI_SC, pol, npol);
        cpl_table * oivis_table  = gravi_data_get_oi_vis (vis_data, GRAVI_SC, pol, npol);
        cpl_table * oiwave_table = gravi_data_get_oi_wave (vis_data, GRAVI_SC, pol, npol);
        
        /* (Re) create the column   FDDLi = (FDDL_FTi + FDDL_SCi)/2 */
        gravi_flux_create_fddllin_sc (oiflux_table, linearity_table);

        /* Compute the BETA and GAMMA for each wavelength */
        cpl_table * dispwave_table0;
        dispwave_table0 = gravi_fit_dispersion (oiflux_table, oivis_table,
                                                oiwave_table, &GDrms,
                                                &Amin, &Amax);
        CPLCHECK_NUL ("Cannot compute dispersion");

        /* Co-add the two polarisation. So here we assume the wavelength
         * table are the same for the two polarisation */
        if (pol == 0) {
            dispwave_table = dispwave_table0;
        } else {
            gravi_msg_warning ("FIXME", "Assumes same OI_WAVE for both polar of SC");
            gravi_table_add_columns (dispwave_table, "BETA", dispwave_table0, "BETA");
            gravi_table_multiply_scalar (dispwave_table, "BETA", 0, 1, 0.5);
            gravi_table_add_columns (dispwave_table, "GAMMA", dispwave_table0, "GAMMA");
            gravi_table_multiply_scalar (dispwave_table, "GAMMA", 0, 1, 0.5);
            cpl_table_delete (dispwave_table0);
        }
    } /* End loop on polarisations */


    /* Set a QC parameters */
    qc_name = "ESO QC DISP GDELAY_RMS";
    cpl_propertylist_update_double (disp_header, qc_name, GDrms);
    cpl_propertylist_set_comment (disp_header, qc_name, "[m] GDELAY rms over files");

    qc_name = "ESO QC DISP BETA_CORRECTION MIN";
    cpl_propertylist_update_double (disp_header, qc_name, Amin);
    cpl_propertylist_set_comment (disp_header, qc_name, "Fine correction");

    qc_name = "ESO QC DISP BETA_CORRECTION MAX";
    cpl_propertylist_update_double (disp_header, qc_name, Amax);
    cpl_propertylist_set_comment (disp_header, qc_name, "Fine correction");

    qc_name = "ESO QC DISP BETA_CORRECTION RANGE";
    cpl_propertylist_update_double (disp_header, qc_name, GRAVI_ACOEFF_RANGE);
    cpl_propertylist_set_comment (disp_header, qc_name, "Fine correction");
    
    
    /* 
     * Interpolate BETA and GAMMA at known Argon wavelength
     * WAVE is the position of the argon line on the current OI_WAVE 
     * WAVE_TH is the true, vaccum line wavelength 
     */
    
    cpl_table * pos_table = gravi_data_get_table (vis_data, "POS_ARGON");
    cpl_table * dispth_table = cpl_table_duplicate (pos_table);
    cpl_size nline = cpl_table_get_nrow (dispth_table);

    gravi_table_interpolate_column (dispth_table, "WAVE", "BETA",
                                    dispwave_table, "EFF_WAVE", "BETA");

    gravi_table_interpolate_column (dispth_table, "WAVE", "GAMMA",
                                    dispwave_table, "EFF_WAVE", "GAMMA");

    CPLCHECK_NUL ("Cannot interpolate into argon lines");

    /*
     * Compute the optical index N_MEAN and N_DIFF from BETA and GAMMA 
     * N_MEAN = BETA  * WAVE_TH / LAMBDA_MET
     * N_DIFF = GAMMA * WAVE_TH / LAMBDA_MET
     */
    
    cpl_table_duplicate_column (dispth_table, "N_MEAN", dispth_table, "BETA");
    cpl_table_duplicate_column (dispth_table, "N_DIFF", dispth_table, "GAMMA");
    
    for (cpl_size line = 0; line < nline; line++) {
        double value = cpl_table_get (dispth_table, "WAVE_TH", line, NULL) / LAMBDA_MET;
        cpl_array_multiply_scalar (cpl_table_get_data_array (dispth_table, "N_MEAN")[line], value);
        cpl_array_multiply_scalar (cpl_table_get_data_array (dispth_table, "N_DIFF")[line], value);
    }
    
    
    /* 
     * Create the output table from the linearity table 
     */
    cpl_table * disp_table = linearity_table;
    
    
    /* 
     * Fit dispersion by a polynomial of order 3 and fill
     * the output table
     * SG 2019-08-02: increased order from 2 to 3
     */
    
    cpl_size mindeg = 0, maxdeg = 3;
    gravi_table_new_column_array (disp_table, "N_MEAN", NULL, CPL_TYPE_DOUBLE, maxdeg+1);
    gravi_table_new_column_array (disp_table, "N_DIFF", NULL, CPL_TYPE_DOUBLE, maxdeg+1);
    gravi_table_new_column (disp_table, "WAVE0", "m", CPL_TYPE_DOUBLE);
    gravi_table_new_column_array (disp_table, "BETA",  NULL, CPL_TYPE_DOUBLE, maxdeg+1);
    gravi_table_new_column_array (disp_table, "GAMMA", NULL, CPL_TYPE_DOUBLE, maxdeg+1);

    /* Allocation of the fit */
    cpl_matrix * matrix = cpl_matrix_new (1, nline);
    cpl_vector * vector = cpl_vector_new (nline);
    cpl_polynomial * poly = cpl_polynomial_new (1);
    cpl_array * coeff = cpl_array_new (maxdeg+1, CPL_TYPE_DOUBLE);
    
    /* The axis for the 2.2e-6/wave_th - 1 */
    double wave0 = 2.2e-6;
    for (cpl_size line = 0; line < nline; line++) {
        double wave_th = cpl_table_get (dispth_table, "WAVE_TH", line, NULL);
        cpl_matrix_set (matrix, 0, line, wave0/wave_th - 1.);
    }

    for (cpl_size tel = 0; tel < ntel; tel++) {
        cpl_table_set (disp_table, "WAVE0", tel, wave0);
        
        /* Fit the BETA */
        for (cpl_size line = 0; line < nline; line++)
            cpl_vector_set (vector, line, gravi_table_get_value (dispth_table, "BETA", line, tel));
        cpl_polynomial_fit (poly, matrix, NULL, vector, NULL, CPL_FALSE, &mindeg, &maxdeg);
        for (cpl_size order = 0; order <= maxdeg; order ++)
            cpl_array_set (coeff, order, cpl_polynomial_get_coeff (poly, &order));
        cpl_table_set_array (disp_table, "BETA", tel, coeff);
        
        /* Fit the GAMMA */
        for (cpl_size line = 0; line < nline; line++)
            cpl_vector_set (vector, line, gravi_table_get_value (dispth_table, "GAMMA", line, tel));
        cpl_polynomial_fit (poly, matrix, NULL, vector, NULL, CPL_FALSE, &mindeg, &maxdeg);
        for (cpl_size order = 0; order <= maxdeg; order ++)
            cpl_array_set (coeff, order, cpl_polynomial_get_coeff (poly, &order));
        cpl_table_set_array (disp_table, "GAMMA", tel, coeff);

        /* Fit the N_MEAN */
        for (cpl_size line = 0; line < nline; line++)
            cpl_vector_set (vector, line, gravi_table_get_value (dispth_table, "N_MEAN", line, tel));
        cpl_polynomial_fit (poly, matrix, NULL, vector, NULL, CPL_FALSE, &mindeg, &maxdeg);
        for (cpl_size order = 0; order <= maxdeg; order ++)
            cpl_array_set (coeff, order, cpl_polynomial_get_coeff (poly, &order));
        cpl_table_set_array (disp_table, "N_MEAN", tel, coeff);
        
        /* Fit the N_DIFF */
        for (cpl_size line = 0; line < nline; line++)
            cpl_vector_set (vector, line, gravi_table_get_value (dispth_table, "N_DIFF", line, tel));
        cpl_polynomial_fit (poly, matrix, NULL, vector, NULL, CPL_FALSE, &mindeg, &maxdeg);
        for (cpl_size order = 0; order <= maxdeg; order ++)
            cpl_array_set (coeff, order, cpl_polynomial_get_coeff (poly, &order));
        cpl_table_set_array (disp_table, "N_DIFF", tel, coeff);
    }
    CPLCHECK_NUL ("Cannot fit the dispersion coefficients");
    
    FREE (cpl_vector_delete, vector);
    FREE (cpl_matrix_delete, matrix);
    FREE (cpl_polynomial_delete, poly);
    FREE (cpl_array_delete, coeff);
    
    
    /* 
     * Output data 
     */

    /* Add the DISP_MODEL in the output gravi_data */
    gravi_data_add_table (disp_map, NULL, "DISP_MODEL", disp_table);

    /* Add the DISP_WAVE in the output gravi_data */
    gravi_data_add_table (disp_map, NULL, "DISP_WAVE", dispwave_table);

    /* Add the DISP_WAVETH in the output gravi_data */
    gravi_data_add_table (disp_map, NULL, "DISP_WAVETH", dispth_table);

    
    gravi_msg_function_exit(1);
    return disp_map;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Cleanup a VIS gravi_data before calibrating the dispersion
 *
 * @param vis_data     The VIS data, modified in-place
 * 
 * \exception CPL_ERROR_NULL_INPUT input data is missing
 *
 * The function erase all observation with low FT visibilities (could be
 * tracking on second lobes), and then keep only those with the longuest
 * LKDT sequence.
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_disp_cleanup (gravi_data * vis_data)
{
    gravi_msg_function_start(1);
    cpl_ensure_code (vis_data, CPL_ERROR_NULL_INPUT);

    cpl_size nbase = 6, ntel = 4, nclo = 4;

    cpl_propertylist * header = gravi_data_get_header (vis_data);
    cpl_size npol_sc = gravi_pfits_get_pola_num (header, GRAVI_SC);
    cpl_size npol_ft = gravi_pfits_get_pola_num (header, GRAVI_FT);

    /* Get one FT OI_VIS table */
    cpl_table * oivis_table  = gravi_data_get_oi_vis (vis_data, GRAVI_FT, 0, npol_ft);
    cpl_table * oiflux_table = gravi_data_get_oi_flux (vis_data, GRAVI_SC, 0, npol_sc);
    cpl_size nrow = cpl_table_get_nrow (oivis_table) / nbase;

    /* Rejection flag (1 = rejected) */
    cpl_array * flag_array = cpl_array_new (nrow, CPL_TYPE_INT);
    cpl_array_fill_window (flag_array, 0, nrow, 0);

    /* 
     * Verify the visibility amplitude 
     */
    
    for (cpl_size row = 0; row < nrow; row++) {
        for (cpl_size base = 0; base < nbase; base++) {
            cpl_size id = row * nbase + base;
            double vis = cpl_array_get_median (cpl_table_get_array (oivis_table, "VISAMP", id));
            cpl_msg_debug ("TEST", "vis = %g", vis);
            if ( vis < 0.35) cpl_array_set (flag_array, row, 1);
        }
    }
    CPLCHECK_MSG ("Cannot compute flag_array");

    /* 
     * Flag on lockdate LKDT
     */

    /* Get longuest sequence */
    cpl_size first = 0, nobs = 0;
    gravi_lkdt_get_sequence (oiflux_table, 4, &first, &nobs);

    /* Flag all observations outside this sequence */
    for (cpl_size row = 0; row < nrow; row++) {
        if (row < first || row >= first+nobs)
            cpl_array_set (flag_array, row, 1);
    }

    if (nobs != nrow) {
        cpl_msg_warning (cpl_func, "LKDT not stable over all files "
                         "(keep %lld over %lld)", nobs, nrow);
    } else {
        cpl_msg_info (cpl_func, "LKDT stable over all files");
    }

    /* 
     * Cleanup all tables with this flag_array (rejection flag)
     */
    
    for (int type_data = 0; type_data < 2; type_data ++) {
        cpl_size npol = gravi_pfits_get_pola_num (header, type_data);
        for (int pol = 0; pol < npol; pol ++) {
            gravi_vis_erase_obs (gravi_data_get_oi_flux (vis_data, type_data, pol, npol), flag_array, ntel);
            gravi_vis_erase_obs (gravi_data_get_oi_vis (vis_data, type_data, pol, npol), flag_array, nbase);
            gravi_vis_erase_obs (gravi_data_get_oi_vis2 (vis_data, type_data, pol, npol), flag_array, nbase);
            gravi_vis_erase_obs (gravi_data_get_oi_t3 (vis_data, type_data, pol, npol), flag_array, nclo);
            CPLCHECK_MSG ("Cannot erase flagged observations");
        }
    }
    FREE (cpl_array_delete, flag_array);

    /* Verbose */
    cpl_size nrow_new = cpl_table_get_nrow (oivis_table) / nbase;
    cpl_msg_info (cpl_func, "Initial data had %lld obs, now %lld", nrow, nrow_new);

    gravi_msg_function_exit(1);
    return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Compute the linearity coefficient of FDDLs
 * 
 * @param oiflux_table   The input OI_FLUX table
 * @return a table with the coefficients [um/V^i]
 * \exception CPL_ERROR_NULL_INPUT input data is missing
 *
 * Found the 22 parameters 6 Aij + 4 B1i + 4 B2i + 4 C1i + 4 C2i
 * by solving the linear system:
 * METit - METjt = Aij + B1i FT_POSit    - B1j FT_POSjt
 *                     + B2i FT_POSit**2 - B2j FT_POSjt**2  
 *                     - C1i SC_POSit    + C1j SC_POSjt
 *                     - C2i SC_POSit**2 + C2j SC_POSjt**2
 *
 * The input OI_FLUX table shall contain several observation at various
 * FDDL position, so that the system is invertible. It shall contains the
 * columns OPD_MET_FC, FT_POS, SC_POS.
 *
 * The output table contains 4 rows, one per beam i, and columns
 * LIN_FDDL_FT with the B0 (=0), B1, B2 coefficients
 * LIN_FDDL_SC with the C0 (=0), C1, C2 coefficients
 */
/*----------------------------------------------------------------------------*/

cpl_table * gravi_fit_fddl_lin (cpl_table * oiflux_table)
{
    gravi_msg_function_start(1);
	cpl_ensure (oiflux_table, CPL_ERROR_NULL_INPUT, NULL);
    
    cpl_size ntel = 4, nbase = 6;
	cpl_size nrow  = cpl_table_get_nrow (oiflux_table) / ntel;

    /* Model and right-hand-side for the lineary system */
    cpl_matrix * rhs_matrix   = cpl_matrix_new (nrow * nbase, 1);
    cpl_matrix * model_matrix = cpl_matrix_new (nrow * nbase, nbase + 4 * ntel);

    for (int base = 0; base < nbase; base++) {
        int i = GRAVI_BASE_TEL[base][0];
        int j = GRAVI_BASE_TEL[base][1];
    
        for (cpl_size row=0; row<nrow; row++) {
            int id  = row * nbase + base;
            int idi = row * ntel + i;
            int idj = row * ntel + j;

            /* Fill the MET [um] */
            double meti = cpl_table_get (oiflux_table, "OPD_MET_FC", idi, NULL);
            double metj = cpl_table_get (oiflux_table, "OPD_MET_FC", idj, NULL);
            cpl_matrix_set (rhs_matrix, id, 0, (meti - metj)*1e6);
            
            /* Fill the model Aij (unfilled matrix are 0.0) */
            cpl_matrix_set (model_matrix, id, base, 1.0);
            
            /* Fill the model Bi, Bj, Ci, Cj */
            double ft_posi = cpl_table_get (oiflux_table, "FT_POS", idi, NULL);
            double ft_posj = cpl_table_get (oiflux_table, "FT_POS", idj, NULL);
            cpl_matrix_set (model_matrix, id, 6 +i,    ft_posi);
            cpl_matrix_set (model_matrix, id, 6 +j, -1*ft_posj);
            cpl_matrix_set (model_matrix, id, 10+i,    ft_posi*ft_posi);
            cpl_matrix_set (model_matrix, id, 10+j, -1*ft_posj*ft_posj);
            
            /* Fill the model Di, Dj, Ei, Ej */
            double sc_posi = cpl_table_get (oiflux_table, "SC_POS", idi, NULL);
            double sc_posj = cpl_table_get (oiflux_table, "SC_POS", idj, NULL);
            cpl_matrix_set (model_matrix, id, 14+i, -1*sc_posi);
            cpl_matrix_set (model_matrix, id, 14+j,    sc_posj);
            cpl_matrix_set (model_matrix, id, 18+i, -1*sc_posi*sc_posi);
            cpl_matrix_set (model_matrix, id, 18+j,    sc_posj*sc_posj);
        } /* End loop on rows */
    } /* End loop on bases */
    
    /* Solve the system */
    cpl_matrix * res_matrix = cpl_matrix_solve_normal (model_matrix, rhs_matrix);
    FREE (cpl_matrix_delete, model_matrix);
    FREE (cpl_matrix_delete, rhs_matrix);


    /* 
     * Fill the linearity coefficients in the output table
     */
    cpl_table * lin_table = cpl_table_new (ntel);
    gravi_table_new_column_array (lin_table, "LIN_FDDL_SC", "um/V^i", CPL_TYPE_DOUBLE, 3);
    gravi_table_new_column_array (lin_table, "LIN_FDDL_FT", "um/V^i", CPL_TYPE_DOUBLE, 3);

    cpl_array * coeff = cpl_array_new (3, CPL_TYPE_DOUBLE);
    for (cpl_size tel = 0; tel < ntel; tel++) {
        cpl_array_set (coeff, 0, 0);
        cpl_array_set (coeff, 1, cpl_matrix_get (res_matrix, 6 +tel, 0));
        cpl_array_set (coeff, 2, cpl_matrix_get (res_matrix, 10+tel, 0));
        cpl_table_set_array (lin_table, "LIN_FDDL_FT", tel, coeff);
        cpl_array_set (coeff, 0, 0);
        cpl_array_set (coeff, 1, cpl_matrix_get (res_matrix, 14+tel, 0));
        cpl_array_set (coeff, 2, cpl_matrix_get (res_matrix, 18+tel, 0));
        cpl_table_set_array (lin_table, "LIN_FDDL_SC", tel, coeff);
        CPLCHECK_NUL ("Cannot set dispersion coeff");
    }
    
    /* Free results */
    FREE (cpl_matrix_delete, res_matrix);
    FREE (cpl_array_delete, coeff);

    gravi_msg_function_exit(1);
    return lin_table;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Compute the dispersion coefficient of FDDLs
 * 
 * @param oivis_table    The input OI_VIS table
 * @param oiflux_table   The input OI_FLUX table
 * @param oiwave_table   The input OI_WAVELENGTH table
 * @param GDrms          GDELAY RMS over obs (max value over base) [m]
 * @param Amin           Minimum fine correction in wavenumber [m^-1]
 * @param Amin           Maximum fine correction in wavenumber [m^-1]
 * 
 * @return a table with the coefficients
 *
 * \exception CPL_ERROR_NULL_INPUT input data is missing
 *
 * Found the  14 parameters 6Aij + 4Bi + CGi in 
 * by solving the linear system:
 *
 * PHASEijt * LAMBDA_MET / 2pi = Aij + Ci (FDDL_FTit + FDDL_SCit)/2 -
 *  Cj (FDDL_FTjt + FDDL_SCjt)/2 + Bi METit - Bj METjt
 *
 * for all wavelength independently
 *
 * The routine takes a great care to unwrap the phase 
 * before actually performing this fit.
 *
 * The input tables shall contain several observation at various
 * FDDL position, so that the system is invertible. The OI_VIS table
 * shall contain the columns VISDATA, OPD_MET_FC. The OI_FLUX table
 * shall contain the columns OPD_MET_FC, FDDL.
 * 
 * Note that the VISDATA are modified in-place by the routine.
 *
 * The output table contains nwave rows (size of OI_WAVELENGTH table)
 * with the columns BETA (the 4 Bi coefficients) and GAMMA (the 4
 * Ci coefficients). It also have a column EFF_WAVE duplicated from
 * the input OI_WAVELENGTH table.
 */
/*----------------------------------------------------------------------------*/

cpl_table * gravi_fit_dispersion (cpl_table * oiflux_table,
                                  cpl_table * oivis_table,
                                  cpl_table * oiwave_table,
                                  double * GDrms,
                                  double * Amin,
                                  double * Amax)
{
    gravi_msg_function_start(1);
	cpl_ensure (oiflux_table, CPL_ERROR_NULL_INPUT, NULL);
	cpl_ensure (oivis_table,  CPL_ERROR_NULL_INPUT, NULL);
	cpl_ensure (oiwave_table, CPL_ERROR_NULL_INPUT, NULL);
    
    cpl_size nbase = 6, ntel = 4;
	cpl_size nrow  = cpl_table_get_nrow (oiflux_table) / 4;
    cpl_size nwave = cpl_table_get_column_depth (oiflux_table, "FLUX");
    CPLCHECK_NUL ("Cannot get data");

    /* 
     * Compute a guess of the BETA dispersion coefficient
     */
    double beta0 = 0.8651, beta1 = 0.8814;
    
    cpl_table_new_column (oiwave_table, "BETA", CPL_TYPE_DOUBLE);
    for (cpl_size wave = 0; wave < nwave; wave ++) {
        double lbd  = cpl_table_get (oiwave_table, "EFF_WAVE", wave, NULL);
        double beta = beta0 + beta1 * (2.2e-6/lbd - 1.0);
        cpl_table_set (oiwave_table, "BETA", wave,  beta);
    }
    CPLCHECK_NUL ("Cannot create BETA column");
    
    /* Get direct pointer to data */
    double * metdata   = cpl_table_get_data_double (oivis_table, "OPD_MET_FC");
    double complex ** visdata = gravi_table_get_data_array_double_complex (oivis_table, "VISDATA");
    double * beta    = cpl_table_get_data_double (oiwave_table, "BETA");
    float  * effwave = cpl_table_get_data_float (oiwave_table, "EFF_WAVE");
    CPLCHECK_NUL ("Cannot get data");

    
    /* 
     * Correction par la metrologie  (Correction # 1)
     * VIS_ijlt *= exp (-2ipi * BETA_l * METC_ijt / LAMBDA_MET)
     */
    cpl_msg_info (cpl_func, "Correction #1");
    
    /* Loop on base, rows and wave */
    for (cpl_size base = 0; base < nbase ; base ++) {
        for (cpl_size row = 0; row < nrow ; row ++) {
            int id  = row * nbase + base;
            for (cpl_size wave = 0; wave < nwave; wave ++) {
                visdata[id][wave] *= cexp (- 2*I*CPL_MATH_PI * beta[wave] * metdata[id] / LAMBDA_MET);
            }
        }
    }
        
    
    /* 
     * Correction par le groupe delay  (Correction # 2)
     * VIS_ijlt *= exp (-2ipi * GD_ij / LBD_l)
     * with GD_ij = <GD_ijt>
     */
    cpl_msg_info (cpl_func, "Correction #2");

    /* Compute the GD of all base and rows (with wavelength in glass) */
	gravi_table_compute_group_delay (oivis_table, "VISDATA", "GDELAY", oiwave_table);

    /* Allocate memory for result */
    cpl_vector * GDb = cpl_vector_new (nbase);
    
    for (cpl_size base = 0; base < nbase ; base ++) {
        double mean = gravi_table_get_column_mean (oivis_table, "GDELAY", base, nbase);
        double std  = gravi_table_get_column_std (oivis_table, "GDELAY", base, nbase);
        cpl_vector_set (GDb, base, mean);

        /* Correct from the mean group-delay */
        for (cpl_size row = 0; row < nrow ; row ++) {
            int id  = row * nbase + base;
            for (cpl_size wave = 0; wave < nwave; wave ++) {
                visdata[id][wave] *= cexp (-2*I*CPL_MATH_PI * mean / effwave[wave]);
            }
        }
        
        cpl_msg_info (cpl_func, "GD mean = %g [um]", mean*1e6);
        cpl_msg_info (cpl_func, "GD std  = %g [um]", std*1e6);

        /* Save the overall worst value of GD rms */
        *GDrms = CPL_MAX (std, *GDrms);
    }

    
    /*
     * Correction from the residual slope versus met  (Correction # 3)
     * VIS_ijlt *= exp (-2ipi * A_bl * METC_ijt / LAMBDA_MET)
     * 
     * Where A_bl is computed such that the following is maximum:
     *  | Sum_t[ VIS_ijlt * exp (-2ipi * A_bl * METC_ijt / LAMBDA_MET) ] |
     *
     * Hence the A value is a BETA coefficient fine correction.
     */
    cpl_msg_info (cpl_func, "Correction #3");

    /* Allocate memory for force-brut exploration of A values */
    cpl_size nA = 1000;
    double complex * phasor = cpl_calloc (nrow * nA, sizeof (double complex));
    // cpl_vector * plot_vector = cpl_vector_new (nA);

    /* Allocate memory for results */
    cpl_matrix * Abl = cpl_matrix_new (nbase, nwave);

    /* Loop on base and wave */
    for (cpl_size base = 0; base < nbase ; base ++) {
        for (cpl_size wave = 0; wave < nwave ; wave ++) {

            /* Test various possible A value */
            cpl_size iAmax;
            double maxV = 0.0;
            for (cpl_size iA = 0; iA < nA; iA++) {
                double A = GRAVI_ACOEFF_RANGE * (2.* iA / nA - 1.0);

                /* Accumulate the re-phased complex Note that we 
                 * compute the exp(i.a.METbm) only if needed */
                double complex currentV = 0.0;
                for (cpl_size row = 0; row < nrow; row++) {
                    if (wave==0) phasor[row*nA+iA] = cexp (-2.* CPL_MATH_PI * I * A *
                                                           metdata[row*nbase+base] /
                                                           LAMBDA_MET);
                    currentV += phasor[row*nA+iA] * visdata[row*nbase+base][wave];
                }

                // if (base == 0 && wave == 1700) cpl_vector_set (plot_vector, iA, cabs (currentV));
                
                /* Check if better fit */
                if (cabs (currentV) > maxV) {
                    cpl_matrix_set (Abl, base, wave, A);
                    iAmax = iA;
                    maxV = cabs (currentV);
                }
            }/* End exploration in A values */

            /* Correct the visdata of this base and wave 
             * from the best-fit A value */
            for (cpl_size row = 0; row < nrow; row++) {
                visdata[row*nbase+base][wave] *= phasor[row*nA+iAmax];
            }

        }
    } /* End loop on base and wave */
    FREE (cpl_free, phasor);

    /* Some verbose */
    cpl_msg_info (cpl_func, "Abl range = %g (beta correction)",
                  GRAVI_ACOEFF_RANGE);
    cpl_msg_info (cpl_func, "Abl mean  = %g (beta correction)",
                  cpl_matrix_get_mean (Abl));
    cpl_msg_info (cpl_func, "Abl std   = %g (beta correction)",
                  cpl_matrix_get_stdev (Abl));

    *Amax = CPL_MAX (*Amax, cpl_matrix_get_max (Abl));
    *Amin = CPL_MIN (*Amin, cpl_matrix_get_min (Abl));


    // cpl_plot_vector (NULL, NULL, NULL, plot_vector);
    // FREE (cpl_vector_delete, plot_vector);


    /* 
     * Remove the mean phase over the rows  (Correction # 4)
     * VIS_ijlt *= exp (-i*O_ijl)
     * with O_ijl = arg (<VIS_ijlt>)
     */
    cpl_msg_info (cpl_func, "Correction #4");

    /* Allocate memory for results */
    cpl_matrix * Obl = cpl_matrix_new (nbase, nwave);
    
    /* Loop on base and wave */
    for (cpl_size base = 0; base < nbase ; base ++) {
        for (cpl_size wave = 0; wave < nwave ; wave ++) {

            /* Compute the mean phase */
            double complex currentV = 0.0;
            for (cpl_size row = 0; row < nrow; row++) {
                currentV += visdata[row*nbase+base][wave];
            }
            cpl_matrix_set (Obl, base, wave, carg (currentV));
            
            /* Correct the visdata of this base and wave 
             * from this mean phase */
            for (cpl_size row = 0; row < nrow; row++) {
                visdata[row*nbase+base][wave] *= conj (currentV);
            }
        }
    } /* End loop on base and wave */

    
    
    /* 
     * Search for 6Aij + 4Bi + 4Ci solving the linear system:
     * PHASEijt * LAMBDA_MET / 2pi = Aij + Ci (FDDL_FTit + FDDL_SCit)/2
     *                                   - Cj (FDDL_FTjt + FDDL_SCjt)/2
     *                                   + Bi METit - Bj METjt
     * for all wavelength and both polaristions
     */
    cpl_msg_info (cpl_func, "Fit dispersion model");

    /* Output of all wavelenths */
    cpl_matrix * disp_fits = cpl_matrix_new (nbase + ntel * 2, nwave);
    
    /* Model and right-hand-side for the lineary system (unfilled matrix are 0.0) */
    cpl_matrix * rhs_matrix   = cpl_matrix_new (nrow * nbase, 1);
    cpl_matrix * model_matrix = cpl_matrix_new (nrow * nbase, nbase + ntel * 2);

    for (cpl_size wave = 0; wave < nwave; wave ++) {

        for (int base = 0; base < nbase; base++) {
            int i = GRAVI_BASE_TEL[base][0];
            int j = GRAVI_BASE_TEL[base][1];
            
            for (cpl_size row=0; row<nrow; row++) {
                int id  = row * nbase + base;
                int idi = row * ntel + i;
                int idj = row * ntel + j;

                /* Fill with unwrap phases from all the corrections
                 * PHIblt = angle(visdata) + Obl + 2pi*Abl*METbt*Abl/LBD_l
                 *          2pi*GDb/LBD_l + 2pi*beta_l*METbt/LAMBDA_MET */
                double phi = carg (visdata[id][wave]);
                phi += cpl_matrix_get (Obl, base, wave);
                phi += CPL_MATH_2PI * cpl_matrix_get (Abl, base, wave) * metdata[id] / LAMBDA_MET;
                phi += CPL_MATH_2PI * cpl_vector_get (GDb, base) / effwave[wave];
                phi += CPL_MATH_2PI * beta[wave] * metdata[id] / LAMBDA_MET;
                cpl_matrix_set (rhs_matrix, id, 0, phi * LAMBDA_MET / CPL_MATH_2PI);
                
                /* Fill the model Aij (unfilled matrix are 0.0) */
                cpl_matrix_set (model_matrix, id, base, 1);
                
                /* Fill the model GAMMAi, GAMMAj */
                double fddli = cpl_table_get (oiflux_table, "FDDL", idi, NULL);
                double fddlj = cpl_table_get (oiflux_table, "FDDL", idj, NULL);
                cpl_matrix_set (model_matrix, id, 6+i,    fddli);
                cpl_matrix_set (model_matrix, id, 6+j, -1*fddlj);

                /* Fill the model BETAi, BETAj */
                double meti = cpl_table_get (oiflux_table, "OPD_MET_FC", idi, NULL);
                double metj = cpl_table_get (oiflux_table, "OPD_MET_FC", idj, NULL);
                cpl_matrix_set (model_matrix, id, 10+i,    meti);
                cpl_matrix_set (model_matrix, id, 10+j, -1*metj);
            } /* End loop on rows */
        } /* End loop on bases */

        /* Solve the system */
        cpl_matrix * res_matrix = cpl_matrix_solve_normal (model_matrix, rhs_matrix);

        /* Save result in disp_fits */
        for (cpl_size param = 0; param < nbase + ntel * 2; param++)
            cpl_matrix_set (disp_fits, param, wave, cpl_matrix_get (res_matrix,param,0));
        FREE (cpl_matrix_delete, res_matrix);
        
    } /* End loop on waves */
    FREE (cpl_matrix_delete, model_matrix);
    FREE (cpl_matrix_delete, rhs_matrix);

    /* Delete pointer to data */
    FREE (cpl_free, visdata);
    // FREE (cpl_table_delete, oiwavefb_table);

    /* Delete corrections */
    FREE (cpl_vector_delete, GDb);
    FREE (cpl_matrix_delete, Abl);
    FREE (cpl_matrix_delete, Obl);

    
    /* Convert the result into DISP_WAVE table,
     * inspired from OI_WAVE */
    cpl_table * dispwave_table = cpl_table_duplicate (oiwave_table);

    /* Add the BETA and GAMMA columns */
    gravi_table_init_column_array (dispwave_table, "BETA",  NULL, CPL_TYPE_DOUBLE, ntel);
    gravi_table_init_column_array (dispwave_table, "GAMMA", NULL, CPL_TYPE_DOUBLE, ntel);
    
    /* Fill the BETA and GAMMA columns */
    for (cpl_size tel = 0; tel < ntel; tel++) {
        for (cpl_size wave = 0; wave < nwave; wave ++) {
            gravi_table_set_value (dispwave_table,"BETA",wave,tel,
                                   cpl_matrix_get (disp_fits, 10+tel, wave));
            gravi_table_set_value (dispwave_table,"GAMMA",wave,tel,
                                   cpl_matrix_get (disp_fits, 6+tel, wave));
            CPLCHECK_NUL ("Cannot fill the dispwave_table");
        }
    }

    /* Delete the matrix */
    FREE (cpl_matrix_delete, disp_fits);
    
	gravi_msg_function_exit (1);
	return dispwave_table;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Compute position of argon lines in SC spectrum
 * 
 * @param preproc_data:   ARGON frame already preproc
 * 
 * \exception CPL_ERROR_NULL_INPUT input data is missing
 * \exception CPL_ERROR_ILLEGAL_INPUT mission table (spectrum or OI_WAVE) in
 * the input data
 *
 * The function computes the position of various argon lines in the SPECTRUM
 * of the SC. Then creates a new table POS_ARGON with the position in
 * [pixel], in wavelength [m] and the expected theoretical wavelength [m] of
 * each line.
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_compute_argon_pos (gravi_data * preproc_data, gravi_data *wave_param)
{
	gravi_msg_function_start(1);
	cpl_ensure_code (preproc_data,  CPL_ERROR_NULL_INPUT);

    /* Get data */
	cpl_table * spectrum_table = gravi_data_get_spectrum_data (preproc_data, GRAVI_SC);
	cpl_size n_region = gravi_spectrum_get_nregion (spectrum_table);
    cpl_size npol = gravi_spectrum_get_npol (spectrum_table);
	cpl_size nwave = gravi_spectrum_get_nwave (spectrum_table);
    CPLCHECK_MSG ("Cannot get data");

    /* Get the OI_WAVE */
    gravi_msg_warning ("FIXME", "Assumes same OI_WAVE for both polar of SC");
	cpl_table * oi_wave = gravi_data_get_oi_wave (preproc_data, GRAVI_SC, 0, npol);

    /* Ensure */
	cpl_ensure_code (spectrum_table,   CPL_ERROR_ILLEGAL_INPUT);
	cpl_ensure_code (oi_wave,          CPL_ERROR_ILLEGAL_INPUT);

    /*
	 * Compute the mean of all the spectra of the argon
	 */

	const cpl_array * array_data;
	cpl_array * argon = cpl_array_new (nwave, CPL_TYPE_DOUBLE);
	cpl_array_fill_window (argon, 0, nwave, 0.0);

	/* Loop region */
	for (cpl_size region = 0; region < n_region; region ++) {
		array_data = cpl_table_get_array (spectrum_table, GRAVI_DATA[region], 0);
		cpl_array_add (argon, array_data);
	}

	cpl_array_divide_scalar (argon, n_region);

	/* 
     * Wavelengths of the argon emission lines [m]
     */
#if(0)
	int nlines = 10;
	double  line_wave[] = {/*1.982291e-6,*/
                   1.997118e-6,
                   2.032256e-6,
                   2.062186e-6,
                   /*2.065277e-6,*/
                   /*2.073922e-6,*/
                   /*2.081672e-6,*/
                   2.099184e-6,
                   2.133871e-6,
                   2.154009e-6,
                   2.208321e-6,
                   2.313952e-6,
                   2.385154e-6,
                   2.397306e-6};
#endif

    double *line_wave;
    int nlines;
    cpl_msg_info (cpl_func, "Vor line_wave2");

	cpl_table * line_wave_table = gravi_data_get_table (wave_param, "ARGON_TAB");
    if ( cpl_table_has_column (line_wave_table , "ARGON_LINES") ) {
        
        line_wave = cpl_table_get_data_double (line_wave_table, "ARGON_LINES");
        nlines = cpl_table_get_nrow (line_wave_table);

        cpl_msg_info (cpl_func,"line_wave [0] : %e", line_wave[0]);
        cpl_msg_info (cpl_func,"line_wave [1] : %e", line_wave[1]);
        cpl_msg_info (cpl_func,"line_wave [2] : %e", line_wave[2]);
        cpl_msg_info (cpl_func,"line_wave [3] : %e", line_wave[3]);
        cpl_msg_info (cpl_func,"line_wave [4] : %e", line_wave[4]);
        cpl_msg_info (cpl_func,"line_wave [5] : %e", line_wave[5]);
        cpl_msg_info (cpl_func,"line_wave [6] : %e", line_wave[6]);
        cpl_msg_info (cpl_func,"line_wave [7] : %e", line_wave[7]);
        cpl_msg_info (cpl_func,"line_wave [8] : %e", line_wave[8]);
        cpl_msg_info (cpl_func,"line_wave [9] : %e", line_wave[9]);
        cpl_msg_info (cpl_func,"nlines        : %d", nlines);
    }
    else {
        /* We cannot continue without line information */
        cpl_msg_error(cpl_func,"Cannot get the default values for Argon line_wave");
        return CPL_ERROR_ILLEGAL_INPUT;
    }

	/*
	 * Fit the position of each emission line
     */

    /* Number of pixels to fit around the line */
	int fitwidth  = nwave > 500 ? 10 : 3;
    int nfitwidth = fitwidth * 2;
    
    /* Create output tables */
	cpl_table * outTable = cpl_table_new (nlines);
	gravi_table_new_column (outTable, "WAVE_TH", "m", CPL_TYPE_DOUBLE);
	gravi_table_new_column (outTable, "WAVE", "m", CPL_TYPE_DOUBLE);
	gravi_table_new_column (outTable, "SIGMA", "m", CPL_TYPE_DOUBLE);
	gravi_table_new_column (outTable, "DIFF", "m", CPL_TYPE_DOUBLE);
	gravi_table_new_column (outTable, "DIFF_PIX", "pix", CPL_TYPE_DOUBLE);
	gravi_table_new_column_array (outTable, "DATA_MEAN", "adu", CPL_TYPE_DOUBLE, nfitwidth);

    /* Allocate vector to extract only sub-part of spectrum */
    cpl_vector * vector_x = cpl_vector_new (nfitwidth);
    cpl_vector * vector_y = cpl_vector_new (nfitwidth);

	for (cpl_size list = 0; list < nlines; list ++) {

        /* Expected position */
        double waveI = line_wave[list];

        /* Expected position in integer [pix] */
        cpl_size pixI = 0;
        while ( cpl_table_get (oi_wave, "EFF_WAVE", pixI, NULL) < waveI) {
            CPLCHECK_MSG ("Cannot get the expected position");
            pixI++;
        }
        
        /* Fill the extracted sub-vector */
		for (cpl_size i = 0; i < nfitwidth; i++) {
            cpl_size w =  pixI - fitwidth + i;
			cpl_vector_set (vector_x, i, cpl_table_get (oi_wave, "EFF_WAVE", w, NULL));
			cpl_vector_set (vector_y, i, cpl_array_get (argon, w, NULL));
		}

		/* Fit Gaussian */
		cpl_errorstate prestate = cpl_errorstate_get();
		double w0 = waveI, sigma, area, offset, mse;
		cpl_vector_fit_gaussian (vector_x, NULL, vector_y, NULL,
                                 CPL_FIT_ALL, &w0, &sigma, &area,
								 &offset, &mse, NULL, NULL);

		if (cpl_error_get_code() == CPL_ERROR_CONTINUE){
			cpl_errorstate_set (prestate);
			cpl_msg_warning(cpl_func, "The gaussian fit did not converge");
            cpl_vector_multiply (vector_x, vector_y);
            w0 = cpl_vector_get_mean (vector_x) /
                 cpl_vector_get_mean (vector_y);
			sigma = 100.0;
		}

        /* compute difference in [m] and [pix] */
        double diff = w0 - waveI;
        double scale = (cpl_vector_get_max (vector_x) - cpl_vector_get_min (vector_x)) / (nfitwidth - 1);
        double diff_pix = diff / scale;
        
        /* Print results */
        cpl_msg_info (cpl_func,"Argon line %lld: %.3g [nm] %.3g [pix] (over %i)",
                      list, 1e9*diff, diff_pix, fitwidth);

        /* Set the result */
		cpl_table_set (outTable, "WAVE_TH", list, waveI);
		cpl_table_set (outTable, "WAVE", list, w0);
		cpl_table_set (outTable, "SIGMA", list, sigma);
		cpl_table_set (outTable, "DIFF", list, diff);
		cpl_table_set (outTable, "DIFF_PIX", list, diff_pix);

        /* Set the extracted part of spectrum for this line */
        cpl_array * tmp_array = cpl_array_wrap_double (cpl_vector_get_data (vector_y), nfitwidth);
        cpl_table_set_array (outTable, "DATA_MEAN", list, tmp_array);
        FREE (cpl_array_unwrap, tmp_array);

		CPLCHECK_MSG ("Error during the computation");
	} /* End loop on list of lines */

    /* Delete vector extraction */
    FREE (cpl_vector_delete, vector_y);
    FREE (cpl_vector_delete, vector_x);
	FREE (cpl_array_delete, argon);
    
	/*
	 * Compute RMS of difference
	 */
	cpl_msg_info (cpl_func, "MIN=%e MAX=%e RMS=%e [nm]",
                  cpl_table_get_column_min (outTable, "DIFF") * 1e9,
                  cpl_table_get_column_max (outTable, "DIFF") * 1e9,
                  cpl_table_get_column_stdev (outTable, "DIFF") * 1e9);

    /* Set the table in gravi_data */
    gravi_data_add_table (preproc_data, NULL, "POS_ARGON", outTable);
    
	gravi_msg_function_exit(1);
	return CPL_ERROR_NONE;
}


/*----------------------------------------------------------------------------*/


/**@}*/
