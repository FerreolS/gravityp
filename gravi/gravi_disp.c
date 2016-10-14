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
 * @defgroup gravi_disp  TBD
 */
/**@{*/

/*----------------------------------------------------------------------------
                                    DEBUG
 -----------------------------------------------------------------------------*/

#define INFO_DEBUG 0
#define GRAVI_ACOEFF_RANGE 2e3

/*-----------------------------------------------------------------------------
                                   Includes
 -----------------------------------------------------------------------------*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#define _XOPEN_SOURCE
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

cpl_vector * gravi_fit_fddl_lin (cpl_table * oiflux_table);

cpl_matrix * gravi_fit_dispersion (cpl_table * oiflux_table,
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
 * The recipe first compute the linearisation coefficients of the FDDL
 * and store them into column LIN_FDDL_SC and LIN_FDDL_FT. Then
 * it computes the dispersion index of the FDDL and store them into
 * columns BETA and GAMMA. This table is then stored as extention
 * DISP_MODEL into the returned, newly allocated, gravi_data.
 *
 * The input vis_data shall have at least 10 observation (60 rows
 * in OI_VIS tables).
 */
/*----------------------------------------------------------------------------*/

gravi_data * gravi_compute_disp (gravi_data * vis_data)
{
    gravi_msg_function_start(1);
	cpl_ensure (vis_data, CPL_ERROR_NULL_INPUT, NULL);

    gravi_msg_warning ("FIXME", "Assumes same OI_WAVE for both polar of SC");

    /* Get data */
    cpl_size ntel = 4, npol = 2;
    cpl_propertylist * vis_header = gravi_data_get_header (vis_data);
	cpl_table * oiflux_table = gravi_data_get_oi_flux (vis_data, GRAVI_SC, 0, npol);
	cpl_table * oiwave_table = gravi_data_get_oi_wave (vis_data, GRAVI_SC, 0, npol);
    CPLCHECK_NUL ("Cannot get data -- need VIS file with 2 polarisation");

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
    
    /* Create the output table */
    cpl_table * disp_table = cpl_table_new (ntel);

    /* 
     * Compute the coefficient of FDDL linearity
     */
    
    cpl_vector * lin_vector; // (22)  in [um/V**i]
    lin_vector = gravi_fit_fddl_lin (oiflux_table);
    CPLCHECK_NUL ("Cannot compute lin_vector");

    /* 
     * Fill the linearity coefficients 
     */
    cpl_size lin_maxdeg = 2;
    gravi_table_new_column_array (disp_table, "LIN_FDDL_SC", "um/V^i", CPL_TYPE_DOUBLE, lin_maxdeg+1);
    gravi_table_new_column_array (disp_table, "LIN_FDDL_FT", "um/V^i", CPL_TYPE_DOUBLE, lin_maxdeg+1);

    cpl_array * coeff = cpl_array_new (lin_maxdeg+1, CPL_TYPE_DOUBLE);
    for (cpl_size tel = 0; tel < ntel; tel++) {
        cpl_array_set (coeff, 0, 0);
        cpl_array_set (coeff, 1, cpl_vector_get (lin_vector, 6 +tel));
        cpl_array_set (coeff, 2, cpl_vector_get (lin_vector, 10+tel));
        cpl_table_set_array (disp_table, "LIN_FDDL_FT", tel, coeff);
        cpl_array_set (coeff, 0, 0);
        cpl_array_set (coeff, 1, cpl_vector_get (lin_vector, 14+tel));
        cpl_array_set (coeff, 2, cpl_vector_get (lin_vector, 18+tel));
        cpl_table_set_array (disp_table, "LIN_FDDL_SC", tel, coeff);
        CPLCHECK_NUL ("Cannot set dispersion coeff");
    }
    FREE (cpl_array_delete, coeff);
    
    /* 
     * Compute the coefficients for FDDL index dispersion
     */
    
    cpl_matrix * disp_matrix; // (14, nwave)  in [refractive index]
    double GDrms = 0.0, Amin = 1e4, Amax = -1e4;

    /* Loop on polarisations */
    for (int pol = 0; pol < npol; pol++) {

        /* Get table for this polarisation */
        cpl_table * oiflux_table = gravi_data_get_oi_flux (vis_data, GRAVI_SC, pol, npol);
        cpl_table * oivis_table  = gravi_data_get_oi_vis (vis_data, GRAVI_SC, pol, npol);
        
        /* (Re) create the column   FDDLi = (FDDL_FTi + FDDL_SCi)/2 */
        gravi_flux_create_fddllin_sc (oiflux_table, disp_table);
        
        cpl_matrix * disp_matrix0;
        disp_matrix0 = gravi_fit_dispersion (oiflux_table, oivis_table,
                                             oiwave_table, &GDrms,
                                             &Amin, &Amax);
        CPLCHECK_NUL ("Cannot compute disp_matrix");

        /* Co-add the two polarisation. So here we assume the wavelength
         * table are the same for the two polarisation */
        if (pol == 0) {
            disp_matrix = disp_matrix0;
        } else {
            cpl_matrix_add (disp_matrix, disp_matrix0);
            cpl_matrix_divide_scalar (disp_matrix, 2.0);
            cpl_matrix_delete (disp_matrix0);
        }
    } /* End loop on polarisations */


    /* Set a QC parameters */
    qc_name = "ESO QC DISP GDELAY_RMS";
    cpl_propertylist_update_double (disp_header, qc_name, GDrms);
    cpl_propertylist_set_comment (disp_header, qc_name, "[m] GDELAY rms over files");

    qc_name = "ESO QC DISP ACOEFF MIN";
    cpl_propertylist_update_double (disp_header, qc_name, Amin);
    cpl_propertylist_set_comment (disp_header, qc_name, "[m^-1] Wavenumber fine correction");

    qc_name = "ESO QC DISP ACOEFF MAX";
    cpl_propertylist_update_double (disp_header, qc_name, Amax);
    cpl_propertylist_set_comment (disp_header, qc_name, "[m^-1] Wavenumber fine correction");

    qc_name = "ESO QC DISP ACOEFF RANGE";
    cpl_propertylist_update_double (disp_header, qc_name, GRAVI_ACOEFF_RANGE);
    cpl_propertylist_set_comment (disp_header, qc_name, "[m^-1] Wavenumber fine correction");
    
    
    /* 
     * Interpolate dispersion at known Argon wavelength
     */

    /* WAVE is the the position of the argon line on the current OI_WAVE 
     * WAVE_TH is the true, vaccum line wavelength */
    cpl_table * pos_table = gravi_data_get_table (vis_data, "POS_ARGON");
    cpl_size nline = cpl_table_get_nrow (pos_table);
    cpl_size nwave = cpl_table_get_nrow (oiwave_table);

    /* Input vector of OI_WAVE */
    cpl_vector * oiwave_vector = cpl_vector_new (nwave);
    for (cpl_size wave = 0; wave < nwave; wave++)
        cpl_vector_set (oiwave_vector, wave, cpl_table_get (oiwave_table, "EFF_WAVE", wave, NULL));

    /* Output vector of line position */
    cpl_vector * wave_vector = cpl_vector_new (nline);
    for (cpl_size line = 0; line < nline; line++)
        cpl_vector_set (wave_vector, line, cpl_table_get (pos_table, "WAVE", line, NULL));
    
    /* Interpolate all coefficients */
    cpl_matrix * displine_matrix; // (14, nline)
    displine_matrix = gravi_matrix_interpolate_col (disp_matrix, oiwave_vector, wave_vector);
    FREE (cpl_vector_delete, oiwave_vector);
    FREE (cpl_vector_delete, wave_vector);


    /* 
     * Fit by Polynomial of order 2 and fill
     * the output table
     */
    
    cpl_size mindeg = 0, maxdeg = 2;
    gravi_table_new_column_array (disp_table, "BETA",  "um^i", CPL_TYPE_DOUBLE, maxdeg+1);
    gravi_table_new_column_array (disp_table, "GAMMA", "um^i", CPL_TYPE_DOUBLE, maxdeg+1);
    gravi_table_new_column (disp_table, "WAVE0", "um", CPL_TYPE_DOUBLE);

    /* Allocation of the fit */
    cpl_matrix * matrix = cpl_matrix_new (1, nline);
    cpl_vector * vector = cpl_vector_new (nline);
    cpl_polynomial * poly = cpl_polynomial_new (1);
    coeff = cpl_array_new (maxdeg+1, CPL_TYPE_DOUBLE);
    
    /* The axis for the 1/wave_th - 1/2.2 [um^-1] */
    double wave0 = 2.2;
    for (cpl_size line = 0; line < nline; line++) {
        double wave_th = cpl_table_get (pos_table, "WAVE_TH", line, NULL);
        cpl_matrix_set (matrix, 0, line, 1.e-6/wave_th - 1./wave0);
    }

    for (cpl_size tel = 0; tel < ntel; tel++) {
        cpl_table_set (disp_table, "WAVE0", tel, wave0);
        
        /* Fit the Bi (FDDL index) */
        for (cpl_size line = 0; line < nline; line++)
            cpl_vector_set (vector, line, cpl_matrix_get (displine_matrix, 10+tel, line));
        cpl_polynomial_fit (poly, matrix, NULL, vector, NULL, CPL_FALSE, &mindeg, &maxdeg);
        for (cpl_size order = 0; order < 10; order ++)
            
        for (cpl_size order = 0; order <= maxdeg; order ++)
            cpl_array_set (coeff, order, cpl_polynomial_get_coeff (poly, &order));
        cpl_table_set_array (disp_table, "BETA", tel, coeff);
        
        /* Fit the Ci (MET index) */
        for (cpl_size line = 0; line < nline; line++)
            cpl_vector_set (vector, line, cpl_matrix_get (displine_matrix, 6+tel, line));
        cpl_polynomial_fit (poly, matrix, NULL, vector, NULL, CPL_FALSE, &mindeg, &maxdeg);
        for (cpl_size order = 0; order <= maxdeg; order ++)
            cpl_array_set (coeff, order, cpl_polynomial_get_coeff (poly, &order));
        cpl_table_set_array (disp_table, "GAMMA", tel, coeff);
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

    /* Duplicate the OI_WAVE into DISP_WAVE */
    cpl_table * dispwave_table = cpl_table_duplicate (oiwave_table);
    gravi_data_add_table (disp_map, NULL, "DISP_WAVE", dispwave_table);

    /* Add the BETA and GAMMA columns */
    gravi_table_init_column_array (dispwave_table, "BETA",  NULL, CPL_TYPE_DOUBLE, ntel);
    gravi_table_init_column_array (dispwave_table, "GAMMA", NULL, CPL_TYPE_DOUBLE, ntel);
    
    /* Fill the BETA and GAMMA columns */
    for (cpl_size tel = 0; tel < ntel; tel++) {
        for (cpl_size wave = 0; wave < nwave; wave ++) {
            gravi_table_set_value (dispwave_table,"BETA",wave,tel,
                                   cpl_matrix_get (disp_matrix, 10+tel, wave));
            gravi_table_set_value (dispwave_table,"GAMMA",wave,tel,
                                   cpl_matrix_get (disp_matrix, 6+tel, wave));
            CPLCHECK_NUL ("Cannot fill the DISP_WAVE");
        }
    }

    
    /* Duplicate the POS_ARGON into DISP_WAVETH */
    cpl_table * dispth_table = cpl_table_duplicate (pos_table);
    gravi_data_add_table (disp_map, NULL, "DISP_WAVETH", dispth_table);

    /* Add the BETA and GAMMA columns */
    gravi_table_init_column_array (dispth_table, "BETA",  NULL, CPL_TYPE_DOUBLE, ntel);
    gravi_table_init_column_array (dispth_table, "GAMMA", NULL, CPL_TYPE_DOUBLE, ntel);
    
    /* Fill the BETA and GAMMA columns */
    for (cpl_size tel = 0; tel < ntel; tel++) {
        for (cpl_size wave = 0; wave < nline; wave ++) {
            gravi_table_set_value (dispth_table,"BETA",wave,tel,
                                   cpl_matrix_get (displine_matrix, 10+tel, wave));
            gravi_table_set_value (dispth_table,"GAMMA",wave,tel,
                                   cpl_matrix_get (displine_matrix, 6+tel, wave));
            CPLCHECK_NUL ("Cannot fill the DISP_WAVETH");
        }
    }

    
    /* Free results */
    FREE (cpl_vector_delete, lin_vector);
    FREE (cpl_matrix_delete, disp_matrix);
    FREE (cpl_matrix_delete, displine_matrix);
    
    gravi_msg_function_exit(1);
    return disp_map;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Cleanup a VIS gravi_data before calibrating the dispersion
 *
 * @param vis_data     The VIS data, modified in-place
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

    cpl_size npol = 2;
    cpl_size nbase = 6, ntel = 4, nclo = 4;

    /* Get one FT OI_VIS table */
    cpl_table * oivis_table  = gravi_data_get_oi_vis (vis_data, GRAVI_FT, 0, npol);
    cpl_table * oiflux_table = gravi_data_get_oi_flux (vis_data, GRAVI_SC, 0, npol);
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
 * @return a vector with 22 coefficients [um/V^i]
 *
 * Found the 22 parameters 6Aij  4Bi  4Ci  4Di  4Ei
 * by solving the linear system:
 * METit - METjt = Aij + Bi FT_POSit    - Bj FT_POSjt
 *                     + Ci FT_POSit**2 - Cj FT_POSjt**2  
 *                     - Di SC_POSit    + Dj SC_POSjt
 *                     - Ei SC_POSit**2 + Ej SC_POSjt**2
 *
 * The input OI_FLUX table shall contain several observation at various
 * FDDL position, so that the system is invertible. It shall contains the
 * columns OPD_MET_FC, FT_POS, SC_POS.
 */
/*----------------------------------------------------------------------------*/

cpl_vector * gravi_fit_fddl_lin (cpl_table * oiflux_table)
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

    /* Convert to vector */
    cpl_vector * fddl_fit = cpl_vector_wrap (22, cpl_matrix_get_data (res_matrix));
    FREE (cpl_matrix_unwrap, res_matrix);
 
    gravi_msg_function_exit(1);
    return fddl_fit;
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
 * @return a matrix with nwave x 14 coefficients
 *
 * Found the  14 parameters 6Aij + 4Bi + 4Ci in 
 * by solving the linear system:
 * PHASEijt * LAMBDA_MET / 2pi = Aij + Bi (FDDL_FTit + FDDL_SCit)/2
 *                                   - Bj (FDDL_FTjt + FDDL_SCjt)/2
 *                                   + Ci METit - Cj METjt
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
 */
/*----------------------------------------------------------------------------*/

cpl_matrix * gravi_fit_dispersion (cpl_table * oiflux_table,
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
     * Calcul of wavelength [m] and wavenumber in glass [m-1]
     * FIXME: verify the scaling !!!
     */
    cpl_table * oiwavefb_table = cpl_table_duplicate (oiwave_table);
    cpl_table_new_column (oiwavefb_table, "EFF_SIGMA", CPL_TYPE_DOUBLE);
    
    for (cpl_size wave = 0; wave < nwave; wave ++) {
        double lbd = cpl_table_get (oiwave_table, "EFF_WAVE", wave, NULL) * 1e6;
        lbd = LAMBDA_MET / (0.8651+1.9396*(1/lbd-1/2.2));
        cpl_table_set (oiwavefb_table, "EFF_WAVE", wave, lbd);
        cpl_table_set (oiwavefb_table, "EFF_SIGMA", wave, 1./lbd);
    }

    /* Get direct pointer to data */
    double * metdata   = cpl_table_get_data_double (oivis_table, "OPD_MET_FC");
    double * sigmadata = cpl_table_get_data_double (oiwavefb_table, "EFF_SIGMA");
    double complex ** visdata = gravi_table_get_data_array_double_complex (oivis_table, "VISDATA");
    CPLCHECK_NUL ("Cannot get data");

    
    /* 
     * Correction par la metrologie  (Correction # 1)
     * VIS_ijlt *= exp (-2ipi * METC_ijt * SIGMA_l)
     */
    cpl_msg_info (cpl_func, "Correction #1");
    
    /* Loop on base, rows and wave */
    for (cpl_size base = 0; base < nbase ; base ++) {
        for (cpl_size row = 0; row < nrow ; row ++) {
            int id  = row * nbase + base;
            for (cpl_size wave = 0; wave < nwave; wave ++) {
                visdata[id][wave] *= cexp (- 2*I*CPL_MATH_PI * metdata[id] * sigmadata[wave]);
            }
        }
    }
        
    
    /* 
     * Correction par le groupe delay  (Correction # 2)
     * VIS_ijlt *= exp (-2ipi * GD_ij * SIGMA_l)
     * with GD_ij = <GD_ijt>
     */
    cpl_msg_info (cpl_func, "Correction #2");

    /* Compute the GD of all base and rows (with wavelength in glass) */
	gravi_table_compute_group_delay (oivis_table, "VISDATA", "GDELAY", oiwavefb_table);

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
                visdata[id][wave] *= cexp (- 2*I*CPL_MATH_PI * mean * sigmadata[wave]);
            }
        }
        
        cpl_msg_info (cpl_func, "GD mean = %g [um]", mean*1e6);
        cpl_msg_info (cpl_func, "GD std  = %g [um]", std*1e6);

        /* Save the overall worst value of GD rms */
        *GDrms = CPL_MAX (std, *GDrms);
    }

    
    /*
     * Correction from the slope versus met  (Correction # 3)
     * VIS_ijlt *= exp (-2ipi * METC_ijt * A_bl)
     * 
     * Where A_bl is computed such that the following is maximum:
     *  | Sum_t[ VIS_ijlt * exp (-2ipi * METC_ijt * A_bl) ] |
     *
     * Hence the A value is a wavenumber fine correction [m^-1].
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
                    if (wave==0) phasor[row*nA+iA] = cexp (-1.* CPL_MATH_2PI * I * A *
                                                           metdata[row*nbase+base]);
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
    cpl_msg_info (cpl_func, "Abl mean = %g [m^-1]",
                  cpl_matrix_get_mean (Abl));
    cpl_msg_info (cpl_func, "Abl std  = %g [m^-1]",
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
     * PHASEijt * LAMBDA_MET / 2pi = Aij + Bi (FDDL_FTit + FDDL_SCit)/2
     *                                   - Bj (FDDL_FTjt + FDDL_SCjt)/2
     *                                   + Ci METit - Cj METjt
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
                 * PHIblt = angle(visdata) + Obl + 2pi/lbdmet*METbt*Abl + 
                 *          2pi*(GDb + METbt) * sigmafb */
                double phi = carg (visdata[id][wave]);
                phi += cpl_matrix_get (Obl, base, wave);
                phi += CPL_MATH_2PI * cpl_matrix_get (Abl, base, wave) * metdata[id];
                phi += CPL_MATH_2PI * (cpl_vector_get (GDb, base) + metdata[id]) * sigmadata[wave];
                cpl_matrix_set (rhs_matrix, id, 0, phi * LAMBDA_MET / CPL_MATH_2PI);
                
                /* Fill the model Aij (unfilled matrix are 0.0) */
                cpl_matrix_set (model_matrix, id, base, 1);
                
                /* Fill the model Bi, Bj */
                double fddli = cpl_table_get (oiflux_table, "FDDL", idi, NULL);
                double fddlj = cpl_table_get (oiflux_table, "FDDL", idj, NULL);
                cpl_matrix_set (model_matrix, id, 6+i,    fddli);
                cpl_matrix_set (model_matrix, id, 6+j, -1*fddlj);

                /* Fill the model Ci, Cj */
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
    FREE (cpl_table_delete, oiwavefb_table);

    /* Delete corrections */
    FREE (cpl_vector_delete, GDb);
    FREE (cpl_matrix_delete, Abl);
    FREE (cpl_matrix_delete, Obl);

	gravi_msg_function_exit (1);
	return disp_fits;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Compute position of argon lines in SC spectrum
 * 
 * @param preproc_data:   ARGON frame already preproc
 * 
 * The function computes the position of various argon lines in the SPECTRUM
 * of the SC. Then creates a new table POS_ARGON with the position in
 * [pixel], in wavelength [m] and the expected theoretical wavelength [m] of
 * each line.
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_compute_argon_pos (gravi_data * preproc_data)
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
