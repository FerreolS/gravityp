/* $Id: gravi_vis.c,v 1.10 2014/11/12 15:10:40 nazouaoui Exp $
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
 * @defgroup gravi_signal     Manipulated signals amoung real-time tables.
 *
 * This module computes all the quantities from the individual coherent flux
 * measurements that are useful for the final processing of the visibilities.
 * It includes the signal synchronization, resampling and averaged over SC DIT, as well
 * as qualitative parameters.
 * The main function called by @c gravity_vis are :
 * - @c gravi_compute_snr() : see Algorithms/Computation of SNR
 * - @c gravi_compute_signals() : see Algorithms/Computing the vFactor,
 * Algorithms/Computing the pFactor, Algorithms/Phase referencing, Algorithms/Geometric flux
 * - @c gravi_compute_rejection() : see Algorithms/Frame rejection
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
#include <complex.h>

#include "gravi_data.h"
#include "gravi_dfs.h"
#include "gravi_pfits.h"
#include "gravi_cpl.h"

#include "gravi_utils.h"
#include "gravi_signal.h"

/*-----------------------------------------------------------------------------
                              Private prototypes
 -----------------------------------------------------------------------------*/

cpl_error_code gravi_vis_bootstrap_snr_and_delay(cpl_table * oi_vis,
                                                 const char * name_snr,
                                                 const char * name_gdl);
cpl_error_code gravi_vis_correct_phasediff(cpl_table * oi_vis1, const char *name1,
                                           cpl_table * oi_vis2, const char *name2,
                                           double * phasediff);
cpl_error_code gravi_vis_compute_mean_phasor(cpl_table * oi_vis,
                                             const char * name_vis,
                                             const char * name_err,
                                             const char * name_pha,
                                             const char * name_var,
					     const char * name_flag);
cpl_error_code gravi_vis_compute_interspectre (cpl_table * oi_vis,
                                               const char * name_vis,
                                               const char * name_is,
                                               const char * name_flag);
cpl_error_code gravi_vis_compute_snr(cpl_table * oi_vis,
                                     const char * name_pha,
                                     const char * name_var,
                                     const char * name_snr);
cpl_error_code gravi_vis_compute_isdelay(cpl_table * oi_vis,
                                         const char * name_isp,
                                         const char * name_gdl,
                                         cpl_table * oi_wavelength);
cpl_error_code gravi_vis_create_pfactor_sc (cpl_table * vis_SC, cpl_table * flux_FT);
cpl_error_code gravi_vis_create_f1f2_sc (cpl_table * vis_SC, cpl_table * flux_SC);
cpl_error_code gravi_vis_create_f1f2_ft (cpl_table * vis_FT, cpl_table * flux_FT);
cpl_error_code gravi_vis_create_phaseref_ft (cpl_table * vis_FT);
cpl_error_code gravi_vis_create_met_sc (cpl_table * vis_SC, cpl_table * vis_MET, cpl_table * wave_table);
cpl_error_code gravi_create_outlier_flag_sc (cpl_table * flux_SC, cpl_table * vis_SC,
					     double chi2r_threshold, double chi2r_sigma);
cpl_error_code gravi_create_outlier_flag_ft (cpl_table * flux_SC, cpl_table * vis_SC);
cpl_error_code gravi_flux_create_fddlpos_sc (cpl_table * flux_SC, cpl_table * fddl_table);
cpl_error_code gravi_flux_create_totalflux_sc (cpl_table * flux_SC, cpl_table * flux_FT);
cpl_error_code gravi_flux_create_met_sc (cpl_table * flux_SC, cpl_table * vis_MET);
cpl_error_code gravi_flux_create_acq_sc (cpl_table * vis_SC,
                                         cpl_table * vis_ACQ);

cpl_error_code gravi_vis_create_acq_sc (cpl_table * vis_SC,
                                        cpl_table * vis_ACQ);

cpl_error_code gravi_vis_create_vfactor_sc (cpl_table * vis_SC,
                                            cpl_table * wave_table_sc,
                                            cpl_table * vis_FT,
                                            cpl_table * wave_table_ft);
cpl_error_code gravi_vis_create_lockratio_sc (cpl_table * vis_SC,
                                              cpl_table * vis_FT);

cpl_error_code gravi_vis_create_phaseref_sc (cpl_table * vis_SC,
                                             cpl_table * wave_table_sc,
                                             cpl_table * wave_table_ft,
					     cpl_propertylist * header,
					     const cpl_parameterlist * parlist);

cpl_error_code gravi_vis_create_opddisp_sc (cpl_table * vis_SC,
                                            cpl_table * flux_SC,
                                            cpl_table * wave_table,
                                            cpl_table * disp_table,
                                            cpl_propertylist * header,
                                            const cpl_parameterlist * parlist);

cpl_error_code gravi_vis_create_imagingref_sc (cpl_table * vis_SC,
                                               cpl_table * wave_table,
                                               cpl_propertylist * header,
                                               const cpl_parameterlist * parlist);

/*-----------------------------------------------------------------------------
                                 Function code
 -----------------------------------------------------------------------------*/

/* -------------------------------------------------------------------------- */
/**
 * @brief Boostrap real-time SNR and GDELAY
 * 
 * @param oi_vis:    input OI_VIS table
 * @param name_snr:  name of input/output SNR column
 * @param name_gdl:  name of input/output group-delay column
 *
 * Compute a 'bootstraped' version of the SNR and GDELAY for each baseline,
 * by looking at the information of other baselines forming closing triangles.
 * Then update the SNR and Group-Delay of this baseline accordingly to the
 * best combination found.
 */
/* -------------------------------------------------------------------------- */

cpl_error_code gravi_vis_bootstrap_snr_and_delay(cpl_table * oi_vis,
                                                 const char * name_snr,
                                                 const char * name_gdl)
{
  gravi_msg_function_start(0);
  cpl_ensure_code (oi_vis,   CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (name_snr, CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (name_gdl, CPL_ERROR_NULL_INPUT);

  /* Get poointer to data */
  double * snr = cpl_table_get_data_double (oi_vis, name_snr);
  double * gdl = cpl_table_get_data_double (oi_vis, name_gdl);
  
  CPLCHECK_MSG("Cannot load data to compute SNR column");

  /* Loop on base */
  cpl_size nrow = cpl_table_get_nrow (oi_vis) / 6;
  for (int base=0; base<6; base++ ) {
	  
	/* Loop on the two closing triangles of this baseline */
	for (int tri=0; tri<2; tri++ ) {
	  
	  /* Search a closing triangle */
	  int b1 = GRAVI_TRI_BASE[base][tri][0];
	  int b2 = GRAVI_TRI_BASE[base][tri][1];
	  int sign1 = GRAVI_TRI_SIGN[base][tri][0];
	  int sign2 = GRAVI_TRI_SIGN[base][tri][1];
	  
	  cpl_msg_debug(cpl_func, "Found triangle tels %i%i -> bases (%i,%i, %i,%i)",
					GRAVI_BASE_TEL[base][0],GRAVI_BASE_TEL[base][1],b1,b2,sign1,sign2);
	  
	  /* Loop on rows */
	  for (cpl_size row = 0; row < nrow; row++) {
		/* Get the bootstraped SNR as the min over the closing baseline */
		double snrN = CPL_MIN( snr[row*6+b1], snr[row*6+b2] );
		if ( snrN > snr[row*6+base] ) {
		  snr[row*6+base] = snrN;
		  gdl[row*6+base] = sign1 * gdl[row*6+b1] + sign2 * gdl[row*6+b2];
		}
	  }
	  /* End loop on rows */
	}
	/* End loop on closing triangle */

	cpl_msg_debug(cpl_func,"TODO: loop on 3-baseline bootstrap");
  }
  /* End loop on base */

  gravi_msg_function_exit(0);
  return CPL_ERROR_NONE;
}

/* -------------------------------------------------------------------------- */
/**
 * @brief Correct for mean phase-difference between coherent fluxes
 * 
 * @param oi_vis1:   first input OI_VIS table (modified inplace)
 * @param name1:     name of coherent flux in oi_vis1
 * @param oi_vis2:   second input OI_VIS table (untouched)
 * @param name2:     name of coherent flux in oi_vis2
 * @param phasedif:  mean phase-difference measured
 *
 * Compute the mean phasor difference between the column name1 and name2.
 * Multiply column name1 by the conjugate of this phasor. Note that name1
 * and name2 should be DOUBLE COMPLEX.
 */
/* -------------------------------------------------------------------------- */

cpl_error_code gravi_vis_correct_phasediff(cpl_table * oi_vis1, const char *name1,
                                           cpl_table * oi_vis2, const char *name2,
                                           double * phasediff)
{
  gravi_msg_function_start(0);
  cpl_ensure_code (oi_vis1, CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (oi_vis2, CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (name1,   CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (name2,   CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (phasediff, CPL_ERROR_ILLEGAL_OUTPUT);

  int nbase = 6;
  
  cpl_type type1 = cpl_table_get_column_type (oi_vis1, name1 );
  cpl_type type2 = cpl_table_get_column_type (oi_vis2, name2 );
  cpl_size nrow1 = cpl_table_get_nrow (oi_vis1) / nbase;
  cpl_size nrow2 = cpl_table_get_nrow (oi_vis2) / nbase;

  if ( type1 != CPL_TYPE_DOUBLE_COMPLEX ||
	   type2 != CPL_TYPE_DOUBLE_COMPLEX ||
	   nrow1 != nrow2)
	return cpl_error_set_message (cpl_func,CPL_ERROR_ILLEGAL_INPUT,"input columns not conformables or not DOUBLE COMPLEX");
  
  double complex * data1 = cpl_table_get_data_double_complex (oi_vis1, name1);
  double complex * data2 = cpl_table_get_data_double_complex (oi_vis2, name2);
  CPLCHECK_MSG("Cannot load data");

  /* Loop on base */
  for (int base = 0; base < nbase; base++) {
	double complex phasor = 0.0;
	
	/* Multiply and sum */
	for (cpl_size row = 0; row < nrow1; row++) {
	  phasor += data1[row*nbase+base] * conj( data2[row*6+base] );
	}

	phasediff[base] = carg( phasor );

	/* Correct phase from oi_vis1 to be the one of oi_vis2 */
	for (cpl_size row = 0; row < nrow1; row++) {
	  data1[row*6+base] *= cexp( -I*phasediff[base] );
	}
  }
  
  gravi_msg_function_exit(0);
  return CPL_ERROR_NONE;
}

/* -------------------------------------------------------------------------- */
/**
 * @brief Compute real-time mean phasor of a VISDATA by
 *        averaging all spectral elements
 * 
 * @param oi_vis:    input OI_VIS table
 * @param name_vis:  name of coherent flux column to be used (VISDATA)
 * @param name_err:  corresponding error column
 * @param name_pha:  name of the column to be created with mean phasor
 * @param name_var:  corresponding variance column to be created
 * @param name_flag: name of flag column
 */
/* -------------------------------------------------------------------------- */

cpl_error_code gravi_vis_compute_mean_phasor(cpl_table * oi_vis,
                                             const char * name_vis,
                                             const char * name_err,
                                             const char * name_pha,
                                             const char * name_var,
					     const char * name_flag)
{
  gravi_msg_function_start(0);
  cpl_ensure_code (oi_vis,   CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (name_vis, CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (name_err, CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (name_pha, CPL_ERROR_ILLEGAL_OUTPUT);
  cpl_ensure_code (name_var, CPL_ERROR_ILLEGAL_OUTPUT);
  cpl_ensure_code (name_flag,CPL_ERROR_ILLEGAL_OUTPUT);

  cpl_size nrow = cpl_table_get_nrow (oi_vis);
  
  /* Init memory for the output array. */
  gravi_table_new_column (oi_vis, name_pha, "e", CPL_TYPE_DOUBLE_COMPLEX);
  double complex * phasorRaw = cpl_table_get_data_double_complex (oi_vis, name_pha);

  gravi_table_new_column (oi_vis, name_var, "e^2", CPL_TYPE_DOUBLE);
  double * varRaw = cpl_table_get_data_double (oi_vis, name_var);
  
  /* Get input pointer to speed up */  
  cpl_array** tVis  = cpl_table_get_data_array (oi_vis, name_vis);
  cpl_array** tErr  = cpl_table_get_data_array (oi_vis, name_err);
  cpl_array** tFlag = cpl_table_get_data_array (oi_vis, name_flag);
  cpl_ensure_code (tFlag, CPL_ERROR_ILLEGAL_INPUT);

  /* Get nwave */
  int nwave = cpl_array_get_size (tVis[0]);

  CPLCHECK_MSG ("Cannot get data");
  
  /* Loop on rows and bases */
  for (cpl_size row = 0; row < nrow; row++ ) {
	
	double complex * cpx = cpl_array_get_data_double_complex (tVis[row]);
	double complex * err = cpl_array_get_data_double_complex (tErr[row]);
	cpl_ensure_code (cpx, CPL_ERROR_ILLEGAL_INPUT);
	cpl_ensure_code (err, CPL_ERROR_ILLEGAL_INPUT);
	phasorRaw[row] = 0.0;
	varRaw[row] = 0.0;

	int * flag = cpl_array_get_data_int (tFlag[row]);

	/* Compute integrated phasor and variance over the spectral channels,
	   ignore flagged values if a flag is provided */
	for (cpl_size wave = 0; wave < nwave; wave++)
	{
	  if (!flag[wave]) {
	    phasorRaw[row] += cpx[wave];
	    varRaw[row] += cabs(err[wave]) * cabs(err[wave]);
	  }
	}
  }
  /* End loop on rows */

  CPLCHECK_MSG("Cannot fill IS and PHASOR column");

  gravi_msg_function_exit(0);
  return CPL_ERROR_NONE;
}

/* -------------------------------------------------------------------------- */
/**
 * @brief Compute the real-time interspectra.
 *   
 * @param oi_vis:    input OI_VIS table
 * @param name_vis:  name of VISDATA column to be used
 * @param name_is:   name of interspectra column to be created
 * @param name_flag: name of flag column
 *
 * The real-time interspectral is computed by averaging the interspectra
 * of all consecutive spectral channels pairs (for each DIT).
 */
/* -------------------------------------------------------------------------- */

cpl_error_code gravi_vis_compute_interspectre (cpl_table * oi_vis,
                                               const char * name_vis,
                                               const char * name_is,
					       const char * name_flag)
{
  gravi_msg_function_start(0);
  cpl_ensure_code (oi_vis,    CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (name_vis,  CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (name_is,   CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (name_flag, CPL_ERROR_NULL_INPUT);
  
  cpl_size nrow = cpl_table_get_nrow (oi_vis);
  
  /* Get input pointer to speed up */  
  cpl_array** tVis = cpl_table_get_data_array (oi_vis, name_vis);
  cpl_ensure_code (tVis, CPL_ERROR_ILLEGAL_INPUT);

  cpl_array** tFlag = cpl_table_get_data_array (oi_vis, name_flag);
  cpl_ensure_code (tFlag, CPL_ERROR_ILLEGAL_INPUT);

  /* Get nwave */
  int nwave = cpl_array_get_size (tVis[0]);

  CPLCHECK_MSG ("Cannot get data");

  /* Init memory for the output array. */
  gravi_table_new_column (oi_vis, name_is, "e^2", CPL_TYPE_DOUBLE_COMPLEX);
  double complex * interSpectreRaw = cpl_table_get_data_double_complex (oi_vis, name_is);
  
  CPLCHECK_MSG ("Cannot create columns");
  
  /* Loop on rows and bases */
  for (cpl_size row = 0; row < nrow; row++) {
	
	double complex * cpx = cpl_array_get_data_double_complex (tVis[row]);
	cpl_ensure_code (cpx, CPL_ERROR_ILLEGAL_INPUT);

	int * flag = cpl_array_get_data_int (tFlag[row]);
	cpl_ensure_code (flag, CPL_ERROR_ILLEGAL_INPUT);

	interSpectreRaw[row] = 0.0;
	
	/* Compute integrated interspectre over the spectral channels,
	   ignore the flagged values if a flag is provided */
	for (cpl_size wave = 0; wave < nwave-1; wave++)
	{
	  if (!flag[wave] && !flag[wave+1]) {
	    interSpectreRaw[row] += cpx[wave] * conj( cpx[wave+1] );
	  }
	}
  }
  /* End loop on rows */

  CPLCHECK_MSG ("Cannot fill IS column");

  gravi_msg_function_exit(0);
  return CPL_ERROR_NONE;
}


/* -------------------------------------------------------------------------- */
/**
 * @brief Compute the real-time SNR
 * 
 * @param oi_vis:   input OI_VIS table
 * @param name_pha: name of phasor column to be used
 * @param name_var: name of variance column to be used
 * @param name_snr: name of output column for snr to be created
 *
 * The real-time SNR is computed as |pha| / sqrt{var}
 */
/* -------------------------------------------------------------------------- */

cpl_error_code gravi_vis_compute_snr (cpl_table * oi_vis,
                                      const char * name_pha,
                                      const char * name_var,
                                      const char * name_snr)
{
  gravi_msg_function_start(0);
  cpl_ensure_code (oi_vis,   CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (name_pha, CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (name_var, CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (name_snr, CPL_ERROR_ILLEGAL_OUTPUT);

  cpl_size nrow = cpl_table_get_nrow (oi_vis);

  /* Get input pointer to data */
  double complex * cpx = cpl_table_get_data_double_complex (oi_vis, name_pha);
  double * dbl = cpl_table_get_data_double (oi_vis, name_var);
  
  CPLCHECK_MSG ("Cannot load data to compute SNR column");

  /* Init memory for the output array. */
  gravi_table_new_column (oi_vis, name_snr, NULL, CPL_TYPE_DOUBLE);
  double * snr = cpl_table_get_data_double (oi_vis, name_snr);

  CPLCHECK_MSG ("Cannot createl SNR column");
  
  /* Loop on rows and bases */
  for (cpl_size row = 0; row < nrow; row++) {
	snr[row] = cabs( cpx[row] ) / sqrt (fabs(dbl[row]));
  }

  CPLCHECK_MSG ("Cannot fill SNR column");

  gravi_msg_function_exit(0);
  return CPL_ERROR_NONE;
}

/* -------------------------------------------------------------------------- */
/**
 * @brief Compute the group-delay from interspectra
 * 
 * @param oi_vis:        input/output OI_VIS table
 * @param name_isp:      name of interspectra column to be used
 * @param name_gdl:      name of output group-delay column to be created
 * @param oi_wavelength: wavelength table corresponding to OI_VIS
 *
 * The group-delay is computed as arg{isp} / 2pi / delta_sigma
 * where delta_sigma = 1/lbd[n/2] - 1/lbd[n/2+1]
 */
/* -------------------------------------------------------------------------- */

cpl_error_code gravi_vis_compute_isdelay (cpl_table * oi_vis,
                                          const char * name_isp,
                                          const char * name_gdl,
                                          cpl_table * oi_wavelength)
{
  gravi_msg_function_start(0);
  cpl_ensure_code (oi_vis,   CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (name_isp, CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (name_gdl, CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (oi_wavelength, CPL_ERROR_NULL_INPUT);
  
  cpl_size nrow = cpl_table_get_nrow (oi_vis);

  /* Get pointer to data */
  double complex * cpx1 = cpl_table_get_data_double_complex (oi_vis, name_isp);

  CPLCHECK_MSG("Cannot load data to compute GDELAY column");

  /* Init memory for the output array. */
  gravi_table_new_column (oi_vis, name_gdl, "m", CPL_TYPE_DOUBLE);
  double * gdl = cpl_table_get_data_double (oi_vis, name_gdl);
  CPLCHECK_MSG("Cannot set GDELAY column");

  /* Compute the coherence */
  cpl_size wave = cpl_table_get_ncol (oi_wavelength)/2;
  double factor = 1./fabs(1./cpl_table_get (oi_wavelength, "EFF_WAVE", wave, NULL) - 
			  1./cpl_table_get (oi_wavelength, "EFF_WAVE", wave+1, NULL));

  cpl_msg_debug (cpl_func, "Compute the coherence length (%e m)", factor);

  CPLCHECK_MSG("Cannot compute coherence length");

  /* Loop on rows and bases */
  for (cpl_size row = 0; row < nrow; row++ ) {
	/* Compute the GD */
	gdl[row] = (double)carg( cpx1[row] ) / (2.0*M_PI) * factor;
  }

  CPLCHECK_MSG("Cannot fill GDELAY column");

  gravi_msg_function_exit(0);
  return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Compute the outliers flags
 * 
 * @param p2vmred_data is the P2VMREDUCED data (modified in-place)
 * 
 * Fill the FLAG columns in the P2VMREDUCED data
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_compute_outliers (gravi_data * p2vmred_data,
				       const cpl_parameterlist * parlist)
{
  gravi_msg_function_start(1);
  cpl_ensure_code (p2vmred_data, CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (parlist,      CPL_ERROR_NULL_INPUT);

  char qc_name[100];
  int ntel = 4;

  cpl_propertylist * p2vmred_header = gravi_data_get_header (p2vmred_data);
  cpl_ensure_code (p2vmred_header, CPL_ERROR_ILLEGAL_INPUT);

  /* Loop on polarisations and type_data */
  for (int type_data = 0; type_data < 2; type_data ++) {
    
    int npol = gravi_pfits_get_pola_num (p2vmred_header, type_data);
    for (int pol = 0; pol < npol; pol++) {

        cpl_table * oi_vis  = gravi_data_get_oi_vis (p2vmred_data, type_data, pol, npol);
        cpl_table * oi_flux = gravi_data_get_oi_flux (p2vmred_data, type_data, pol, npol);
	cpl_size nrow = cpl_table_get_nrow (oi_flux) / ntel;
	cpl_size nwave = cpl_table_get_column_depth (oi_flux, "FLUX");
        CPLCHECK_MSG ("Cannot get data");

	if (type_data == GRAVI_FT)
	{
	  gravi_create_outlier_flag_ft (oi_flux, oi_vis);
	}
	else
	{
	  double chi2r_threshold = gravi_param_get_double (parlist, "gravity.signal.chi2r-threshold");
	  double chi2r_sigma = gravi_param_get_double (parlist, "gravity.signal.chi2r-sigma");
	  
	  gravi_create_outlier_flag_sc (oi_flux, oi_vis, chi2r_threshold, chi2r_sigma);
	  
	  /* Fraction of chi2 outliers per channel */
	  cpl_array * array = gravi_table_get_column_sum_array (oi_flux, "FLAG", 0, ntel);
	  
	  int stat[2] = {0,0};
	  for (cpl_size wave = 0; wave < nwave ; wave ++) {
	    double value = (double)cpl_array_get (array, wave, NULL) / nrow;
	    if (value > 0)   stat[0]++;
	    if (value > 0.5) stat[1]++;
	  }
	  
	  FREE (cpl_array_delete, array);  
	  CPLCHECK_MSG ("Cannot compute stat");

	  /* QC */
  	  sprintf (qc_name, "ESO QC OUTLIER_SC_P%i", pol+1);
  	  cpl_propertylist_update_int (p2vmred_header, qc_name, stat[0]);
  	  cpl_propertylist_set_comment (p2vmred_header, qc_name, "Channels with at least one outlier");
  	  cpl_msg_info (cpl_func,"%s = %i ", qc_name, stat[0]);
	  
  	  sprintf (qc_name, "ESO QC OUTLIER50_SC_P%i",pol+1);
  	  cpl_propertylist_update_int (p2vmred_header, qc_name, stat[1]);
  	  cpl_propertylist_set_comment (p2vmred_header, qc_name, "Channels with more than 50% outlier");
  	  cpl_msg_info (cpl_func,"%s = %i ", qc_name, stat[1]);
	}
	
	CPLCHECK_MSG ("Cannot create outlier flag");
    }
  }
  /* End loop on polarisation and type_data */
  
  gravi_msg_function_exit(1);
  return CPL_ERROR_NONE;
}
/*----------------------------------------------------------------------------*/
/**
 * @brief Compute real-time SNR and Group-Delay of the observation.
 * 
 * @param p2vmred_data is the P2VMREDUCED data (modified in-place)
 * 
 * Create the SNR_SMT and GDELAY_SMT columns in the P2VMREDUCED data
 * These are the estimate considering both polarisation if any, all
 * closing baseline and a running sum of the complex coherent flux.
 * These quantities are computed for SC and FT.
 *
 * Create column: SNR, SNR_BOOT, GDELAY_BOOT
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_compute_snr (gravi_data * p2vmred_data,
                                  const cpl_parameterlist * parlist)
{
  /* Message and timer */
  gravi_msg_function_start(1);
  cpl_ensure_code (p2vmred_data, CPL_ERROR_NULL_INPUT);

  int nbase = 6;
  char qc_name[100];
  double qc_value;
  
  cpl_propertylist * p2vmred_header = gravi_data_get_header (p2vmred_data);
  double periode_sc = gravi_pfits_get_period_sc (p2vmred_header); // [s]
  CPLCHECK_MSG ("Cannot get header");

  /*
   * Compute the real-time group-delay and SNR
   * as well as a running average.
   */
  
  /* Loop on polarisations and type_data */
  for (int type_data = 0; type_data < 2; type_data ++) {

    if (gravi_data_has_type (p2vmred_data, GRAVI_TYPE(type_data)) < 2 ) {
        cpl_msg_info (cpl_func, "Cannot compute real-time snr for %s "
                      "(no data)", GRAVI_TYPE(type_data));
        continue;
    }

    int nsmooth;
    if (type_data == GRAVI_FT) {
        /* For FT, the DIT is selected to match somehow the atmospheric
         * coherence time, hence we shall smooth over 'a few' samples */
        nsmooth = gravi_param_get_int (parlist, "gravity.signal.nsmooth-snr-ft");
    }
    else {
        /* For SC, we select a smoothing of 1s. Note that the SNR of SC is
         * not used for selecting frame, thus this choice is not critical */
        nsmooth = 1./ periode_sc;
        nsmooth = CPL_MIN (CPL_MAX (nsmooth, 0), 20);
    }

    /* Loop on polarisation */
    int npol = gravi_pfits_get_pola_num (p2vmred_header, type_data);
    for (int pol = 0; pol < npol; pol++) {
  
        cpl_msg_info (cpl_func, "Compute SNR for pol %i of %s (smoothed over %i frames)",
		      pol+1,GRAVI_TYPE(type_data),2*nsmooth+1);
	
        cpl_table * oi_vis  = gravi_data_get_oi_vis (p2vmred_data, type_data, pol, npol);
        CPLCHECK_MSG ("Cannot get data");
	
        /* Compute interspectre and phasor */
        gravi_vis_compute_interspectre (oi_vis, "VISDATA", "IPHASOR", "FLAG");
        gravi_vis_compute_mean_phasor (oi_vis, "VISDATA", "VISERR", "PHASOR", "PHASOR_VAR", "FLAG");

        /* Compute real-time SNR */
        gravi_vis_compute_snr (oi_vis, "PHASOR", "PHASOR_VAR", "SNR");
        CPLCHECK_MSG ("Cannot compute PHASOR and IPHASOR");
        
        /* Compute a smoothed version of quantities. Note that this is a running SUM,
         * not a running MEAN. Hence the SNR is enlarged by the smoothing */
        gravi_table_runint_column (oi_vis, "IPHASOR", "IPHASOR_SMT", nsmooth, nbase);
        gravi_table_runint_column (oi_vis, "PHASOR", "PHASOR_SMT", nsmooth, nbase);
        gravi_table_runint_column (oi_vis, "PHASOR_VAR", "PHASOR_VAR_SMT", nsmooth, nbase);
        CPLCHECK_MSG ("Cannot running integration");
        
    } /* End loop on pol */
  } /* End loop on type_data */
  
  /*
   * Compute bootstraped GD and SNR for the FT.
   * Bootstraped meant with all possible information merging
   * (polar, baseline, time)
   */
  
  /* Loop on type_data */
  for (int type_data = 0; type_data < 2; type_data ++) {

    if (gravi_data_has_type (p2vmred_data, GRAVI_TYPE(type_data)) < 2 ) {
        cpl_msg_info (cpl_func, "Cannot compute bootstraped snr %s "
                      "(no data)", GRAVI_TYPE(type_data));
        continue;
    }

    int npol = gravi_pfits_get_pola_num (p2vmred_header, type_data);
	
    cpl_msg_info(cpl_func, "Compute bootstraped GDELAY_BOOT and SNR for %s...",(type_data==GRAVI_FT?"FT":"SC"));
  
    /* Get the FT data */
    cpl_table * oi_wave = gravi_data_get_oi_wave (p2vmred_data, type_data, 0, npol);
    cpl_table * oi_vis_p2, * oi_vis_p1;
    oi_vis_p1 = gravi_data_get_oi_vis (p2vmred_data, type_data, 0, npol);
  
    CPLCHECK_MSG("Cannot get OI_VIS_P1 table");
  
    /* Duplicate interspectre to get the bootstraped column */
    cpl_table_duplicate_column (oi_vis_p1, "IPHASOR_BOOT", oi_vis_p1, "IPHASOR_SMT");
    CPLCHECK_MSG("Cannot duplicate columns");
  
    /* Duplicate phasor and variance to get the bootstraped column. */
    cpl_table_duplicate_column (oi_vis_p1, "PHASOR_BOOT", oi_vis_p1, "PHASOR_SMT");
    CPLCHECK_MSG("Cannot duplicate columns");
    
    cpl_table_duplicate_column (oi_vis_p1, "PHASOR_VAR_BOOT", oi_vis_p1, "PHASOR_VAR_SMT");
    CPLCHECK_MSG("Cannot duplicate columns");
  
    /* If two polarisations */
    if ( npol > 1) {
  
  	cpl_msg_info(cpl_func, "Add the signal of both polarisation to enhance SNR and GDELAY accuracy");
  	oi_vis_p2 = gravi_data_get_oi_vis (p2vmred_data, type_data, 1, npol);
  	
  	CPLCHECK_MSG("Cannot get OI_VIS_P2 table");
  
  	/* Sum the interspectre of the second polarisation. No need to deal
  	 * with polarisation shifts as they probably have the same GD */
  	gravi_table_add_columns (oi_vis_p1, "IPHASOR_BOOT", oi_vis_p2, "IPHASOR_SMT");
  	CPLCHECK_MSG("Cannot add columns");
  	
  	/* Removing the mean phase difference from PHASOR of first polarisation 
  	 * to be able to sum with maximum SNR */
  	double phasediff[6];
  	gravi_vis_correct_phasediff (oi_vis_p1, "PHASOR_BOOT", oi_vis_p2, "PHASOR_SMT", phasediff);
  	
  	/* Sum the phasor and the variance of the second polarisation. */
  	gravi_table_add_columns (oi_vis_p1, "PHASOR_BOOT", oi_vis_p2, "PHASOR_SMT");
  	gravi_table_add_columns (oi_vis_p1, "PHASOR_VAR_BOOT", oi_vis_p2, "PHASOR_VAR_SMT");
  	
  	/* Dump these phase difference as QC parameters */
  	cpl_msg_info(cpl_func,"Add the phase difference as QC parameters");
  	for (int base = 0; base < nbase; base++) {
  	  sprintf (qc_name, "ESO QC PHASE_POLDIFF_%s%d%d",(type_data==GRAVI_FT?"FT":"SC"),
  			   GRAVI_BASE_TEL[base][0]+1, GRAVI_BASE_TEL[base][1]+1);
  	  cpl_propertylist_update_double (p2vmred_header, qc_name, phasediff[base] * CPL_MATH_DEG_RAD);
  	  cpl_propertylist_set_comment (p2vmred_header, qc_name, "[deg] differential phase" );
  	  cpl_msg_info (cpl_func,"%s=%f [deg]",qc_name,phasediff[base] * CPL_MATH_DEG_RAD);
  	  CPLCHECK_MSG("QC PHASE_POLDIFF");
  	}
  	
  	CPLCHECK_MSG("Cannot add columns of second polarisation");
    }
    /* End case two polarisations */
  
    /* Compute the highest possible GDDELAY and SNR per baseline. */
    gravi_vis_compute_isdelay (oi_vis_p1, "IPHASOR_BOOT", "GDELAY_BOOT",oi_wave);
    gravi_vis_compute_snr (oi_vis_p1, "PHASOR_BOOT", "PHASOR_VAR_BOOT", "SNR_BOOT_TMP");
    
    /* Bootstrap over the baselines. Note the aggressive smooth of SNR to avoid
     * bootstrapoing on noisy SNR (stabilization). */
    int nstabilize = (type_data == GRAVI_FT ? 50 : 1);
    gravi_table_smooth_column (oi_vis_p1, "SNR_BOOT_TMP", "SNR_BOOT", nstabilize, nbase);
    gravi_vis_bootstrap_snr_and_delay (oi_vis_p1, "SNR_BOOT", "GDELAY_BOOT");
    
    CPLCHECK_MSG("Cannot compute and fill SNR_BOOT or GDELAY_BOOT column");
    
    /* If two polarisations. Copy this SNR and GDELAY to the second polarisation */
    if (npol > 1) {
        cpl_msg_info (cpl_func, "Duplicate columns in polarisation 2");
        if (cpl_table_has_column (oi_vis_p2, "SNR_BOOT"))
            cpl_table_erase_column(oi_vis_p2, "SNR_BOOT");
        cpl_table_duplicate_column (oi_vis_p2, "SNR_BOOT", oi_vis_p1, "SNR_BOOT");
        
        if (cpl_table_has_column (oi_vis_p2, "GDELAY_BOOT"))
            cpl_table_erase_column(oi_vis_p2, "GDELAY_BOOT");
        cpl_table_duplicate_column (oi_vis_p2, "GDELAY_BOOT", oi_vis_p1, "GDELAY_BOOT");
        CPLCHECK_MSG("Cannot duplicate column for polarisation 2");
    }
  
    /* Remove useless column (IPHASOR, PHASOR, PHASOR_VAR) when the
     * smoothing and bootstraping are done */
    cpl_msg_info (cpl_func, "Erase useless columns in P2VMREDUCED");
  
    cpl_table_erase_column (oi_vis_p1, "IPHASOR");
    cpl_table_erase_column (oi_vis_p1, "IPHASOR_SMT");
    cpl_table_erase_column (oi_vis_p1, "IPHASOR_BOOT");
    cpl_table_erase_column (oi_vis_p1, "PHASOR");
    cpl_table_erase_column (oi_vis_p1, "PHASOR_SMT");
    cpl_table_erase_column (oi_vis_p1, "PHASOR_BOOT");
    cpl_table_erase_column (oi_vis_p1, "PHASOR_VAR");
    cpl_table_erase_column (oi_vis_p1, "PHASOR_VAR_SMT");
    cpl_table_erase_column (oi_vis_p1, "PHASOR_VAR_BOOT");
    cpl_table_erase_column (oi_vis_p1, "SNR_BOOT_TMP");
    CPLCHECK_MSG("Cannot erase column for polarisation 1");
  
    if ( npol > 1) {
      cpl_table_erase_column (oi_vis_p2, "IPHASOR");
      cpl_table_erase_column (oi_vis_p2, "IPHASOR_SMT");
      cpl_table_erase_column (oi_vis_p2, "PHASOR");
      cpl_table_erase_column (oi_vis_p2, "PHASOR_SMT");
      cpl_table_erase_column (oi_vis_p2, "PHASOR_VAR");
      cpl_table_erase_column (oi_vis_p2, "PHASOR_VAR_SMT");
      CPLCHECK_MSG("Cannot erase column for polarisation 2");
    }

    /* Add the QC on SNR */
    for (cpl_size base = 0; base < nbase; base++) {
        snprintf (qc_name, 100, "ESO QC SNRB_%s%s AVG", GRAVI_TYPE(type_data), GRAVI_BASE_NAME[base]);
        qc_value = gravi_table_get_column_mean (oi_vis_p1, "SNR_BOOT", base, nbase);
        cpl_propertylist_update_double (p2vmred_header, qc_name, qc_value);
        cpl_propertylist_set_comment (p2vmred_header, qc_name, "mean bootstrapped SNR");
    }
  
  } /* End loop on type_data */

  gravi_msg_function_exit(1);
  return CPL_ERROR_NONE;
}


/*----------------------------------------------------------------------------*/
/**
 * @brief Compute synchronisation indices between OIFITS tables
 * 
 * @param vis_SC:   table for which to compute the indices (SC)
 * @param nbase_sc: integer, number of base (6 for SC)
 * @param dit_sc:   double, search window for synch
 * 
 * @param vis_FT:   table where to search for synch (FT, MET)
 * @param nbase_ft: integer, number of base (6 for FT, 1 for MET)
 * @param name:     string, to define the LAST_## and FIRST_## column names
 *
 * Search for the first and last indices in VIS_FT for each frame in
 * vis_SC, based on the TIME columns and the +-dit_sc/2 window, and creates
 * the LAST_## and FIRST_## new columns with those index.
 * The accepted frames are {first to last-1}
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_signal_create_sync (cpl_table * vis_SC, int nbase_sc, double dit_sc,
										 cpl_table * vis_FT, int nbase_ft,
										 const char * name)
{
  gravi_msg_function_start(1);
  cpl_ensure_code (vis_SC, CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (vis_FT, CPL_ERROR_NULL_INPUT);

  /* Get the number of rows */
  cpl_size nrow_sc = cpl_table_get_nrow (vis_SC) / nbase_sc;
  cpl_size nrow_ft = cpl_table_get_nrow (vis_FT) / nbase_ft;
  
  CPLCHECK_MSG ("Cannot get data");

  /* Create columns in vis_SC */
  char full_name[90];
  sprintf (full_name, "FIRST_%s", name);
  gravi_table_new_column (vis_SC, full_name, NULL, CPL_TYPE_INT);
  int * first_ft = cpl_table_get_data_int (vis_SC, full_name);
		
  sprintf (full_name, "LAST_%s", name);
  gravi_table_new_column (vis_SC, full_name, NULL, CPL_TYPE_INT);
  int * last_ft  = cpl_table_get_data_int (vis_SC, full_name);

  sprintf (full_name, "NFRAME_%s", name);
  gravi_table_new_column (vis_SC, full_name, NULL, CPL_TYPE_INT);
  int * nframe_ft = cpl_table_get_data_int (vis_SC, full_name);

  CPLCHECK_MSG ("Cannot create columns");
  
  /* Loop on base (here assume base may have different timing) */
  for (cpl_size base_sc = 0; base_sc < nbase_sc; base_sc++) {

	/* Start info on the second table 
	 * base_ft is always 0 if nbase_ft is 1 */
	cpl_size row_ft = 0;
	cpl_size base_ft = base_sc % nbase_ft;

	/* Loop on SC rows */
	for (cpl_size row_sc = 0; row_sc < nrow_sc; row_sc ++) {

	  /* Get the first FT sample */
	  while ( cpl_table_get (vis_FT, "TIME", row_ft * nbase_ft + base_ft, NULL) <
              cpl_table_get (vis_SC, "TIME", row_sc * nbase_sc + base_sc, NULL) - dit_sc/2.) {
		if (row_ft >= nrow_ft-1) break; // avoid reading outside of the table
		row_ft ++;
	  }
	  first_ft[row_sc * nbase_sc + base_sc] = row_ft;

	  /* Get the last sample */
	  while ( cpl_table_get (vis_FT, "TIME", row_ft * nbase_ft + base_ft, NULL) <
			  cpl_table_get (vis_SC, "TIME", row_sc * nbase_sc + base_sc, NULL) + dit_sc/2.) {
		if (row_ft >= nrow_ft-1) break; // avoid reading outside of the table
		row_ft ++;
	  }
	  last_ft[row_sc * nbase_sc + base_sc] = row_ft;

	  /* Check if enough data */
	  if ( first_ft[row_sc * nbase_sc + base_sc] < 2 ||
		   last_ft[row_sc * nbase_sc + base_sc] > nrow_ft - 2 ) {
        if (!strcmp (name, "ACQ"))
            {
	        cpl_msg_info (cpl_func,"Not enough %s data to synchronise with DIT %lli over %lli", name, row_sc+1, nrow_sc);
            } else {
            cpl_msg_warning (cpl_func,"Not enough %s data to synchronise with DIT %lli over %lli", name, row_sc+1, nrow_sc);
            }
		first_ft[row_sc * nbase_sc + base_sc] = 0;
		last_ft[row_sc * nbase_sc + base_sc]  = 0;
	  }

	  /* Number of frames */
	  nframe_ft[row_sc * nbase_sc + base_sc] = last_ft[row_sc * nbase_sc + base_sc] - first_ft[row_sc * nbase_sc + base_sc];
	}
  }
  /* End loop on SC row and base */
		
  gravi_msg_function_exit(1);
  return CPL_ERROR_NONE;
}

/* -------------------------------------------------------------------------- */
/**
 * @brief Compute the PFACTOR for the SC
 * 
 * @param vis_SC:   output OI_VIS table of the SC
 * @param flux_FT:  input OI_FLUX table of the FT
 *
 * The PFACTOR is computed for each SC DIT and saved in a newly created column
 * P_FACTOR in the vis_SC table. The synchronisation info shall be available.
 */
/* -------------------------------------------------------------------------- */

cpl_error_code gravi_vis_create_pfactor_sc (cpl_table * vis_SC, cpl_table * flux_FT)
{
  gravi_msg_function_start(1);
  cpl_ensure_code (vis_SC,  CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (flux_FT, CPL_ERROR_NULL_INPUT);

  /* Get the number of rows */
  cpl_size nbase = 6, ntel = 4;
  cpl_size nrow_sc = cpl_table_get_nrow (vis_SC) / nbase;

  /* Create the column */
  gravi_table_new_column (vis_SC, "P_FACTOR", NULL, CPL_TYPE_DOUBLE);
  double * pFactor = cpl_table_get_data_double (vis_SC, "P_FACTOR");

  /* Get SC data */
  int * first_ft = cpl_table_get_data_int (vis_SC, "FIRST_FT");
  int * last_ft  = cpl_table_get_data_int (vis_SC, "LAST_FT");

  /* Get FT data (used already computed and smoothed total-flux) */
  double * flux = cpl_table_get_data_double (flux_FT, "TOTALFLUX");

  CPLCHECK_MSG ("Cannot get pointer to data");

  /* Loop on base and SC rows */
  for (cpl_size base = 0; base < nbase; base++) {
	for (cpl_size row_sc = 0; row_sc < nrow_sc; row_sc ++) {
	  cpl_size nsc = row_sc*nbase+base;

	  /* Loop on the sync FT frames */
	  double sf0f1 = 0.0, f0 = 0.0, f1 = 0.0;
	  for (cpl_size row_ft = first_ft[nsc] ; row_ft < last_ft[nsc]; row_ft++) {
		int t0 = GRAVI_BASE_TEL[base][0] + row_ft * ntel;
		int t1 = GRAVI_BASE_TEL[base][1] + row_ft * ntel;
		sf0f1 += sqrt (CPL_MAX(flux[t0] * flux[t1], 0));
		f0   += flux[t0];
		f1   += flux[t1];
	  }

	  /* Discard unused */
	  if (f0==0 || f1==0) continue;
	  
	  /* Compute the pFactor */
	  pFactor[nsc] = sf0f1 * sf0f1 / (f0 * f1);
	}
  }
  
  gravi_msg_function_exit(1);
  return CPL_ERROR_NONE;
}

/* -------------------------------------------------------------------------- */
/**
 * @brief Compute the photometric normalisation for the SC
 * 
 * @param vis_SC:   output OI_VIS table of the SC
 * @param flux_SC:  input OI_FLUX table of the SC
 *
 * The normalisation is computed for each SC DIT and saved in a newly
 * created column F1F2 in the vis_SC table. It is simply the product of
 * the fluxes of the 2 telescopes corresponding to each base.
 */
/* -------------------------------------------------------------------------- */

cpl_error_code gravi_vis_create_f1f2_sc (cpl_table * vis_SC, cpl_table * flux_SC)
{
  gravi_msg_function_start(1);
  cpl_ensure_code (vis_SC,  CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (flux_SC, CPL_ERROR_NULL_INPUT);

  cpl_size nbase = 6, ntel = 4;
  cpl_size nrow_sc = cpl_table_get_nrow (flux_SC) / ntel;
  cpl_size nwave_sc = cpl_table_get_column_depth (flux_SC, "FLUX");

  /* Get pointer to data */
  cpl_array ** flux_sc = cpl_table_get_data_array (flux_SC, "FLUX");

  CPLCHECK_MSG ("Cannot get data");
  
  /* New column */
  gravi_table_new_column_array (vis_SC, "F1F2", "e^2", CPL_TYPE_DOUBLE, nwave_sc);
  cpl_array ** f1f2_sc = cpl_table_get_data_array (vis_SC, "F1F2");

  for (cpl_size base = 0; base < nbase; base++) {
	for (cpl_size row = 0; row < nrow_sc; row ++) {
	  
	  /* Compute the photometric normalisation of SC
	   * This is simply FLUX1 * FLUX2 */
	  f1f2_sc[row*nbase+base] = cpl_array_cast (flux_sc[row*ntel+GRAVI_BASE_TEL[base][0]], CPL_TYPE_DOUBLE);
	  cpl_array_multiply (f1f2_sc[row*nbase+base], flux_sc[row*ntel+GRAVI_BASE_TEL[base][1]]);
	}
  }
		
  gravi_msg_function_exit(1);
  return CPL_ERROR_NONE;
}

/* -------------------------------------------------------------------------- */
/**
 * @brief Compute the photometric normalisation for the FT
 * 
 * @param vis_FT:   output OI_VIS table of the FT
 * @param flux_FT:  input OI_FLUX table of the FT
 *
 * The normalisation is computed for each FT DIT and saved in a newly
 * created column F1F2 in the vis_FT table. To enhance SNR and avoid division
 * by zero, it uses a combination of time-averaged and channel-averaged
 * signals to reconstruct the real-time per-channel photometric flux. The time
 * averaged signal is from TOTALFLUX which shall thus be already computed.
 */
/* -------------------------------------------------------------------------- */

cpl_error_code gravi_vis_create_f1f2_ft (cpl_table * vis_FT, cpl_table * flux_FT)
{
  gravi_msg_function_start(1);
  cpl_ensure_code (vis_FT,  CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (flux_FT, CPL_ERROR_NULL_INPUT);

  cpl_size ntel = 4, nbase = 6;
  cpl_size nrow_ft = cpl_table_get_nrow (flux_FT) / ntel;
  cpl_size nwave_ft = cpl_table_get_column_depth (flux_FT, "FLUX");

  /* Get pointer to data */
  double * total_flux_ft = cpl_table_get_data_double (flux_FT, "TOTALFLUX");
  cpl_array ** flux_ft = cpl_table_get_data_array (flux_FT, "FLUX");

  CPLCHECK_MSG ("Cannot get data");
  
  /* Compute the four mean_spectra */
  cpl_array ** mean_spectra = cpl_malloc (4 * sizeof (cpl_array*));
  
  for (cpl_size tel = 0; tel < ntel; tel ++) {
	mean_spectra[tel] = cpl_array_duplicate (flux_ft[tel]);
	for (cpl_size n = 1; n < nrow_ft; n ++) cpl_array_add (mean_spectra[tel], flux_ft[n*ntel+tel]);
	cpl_array_divide_scalar (mean_spectra[tel], cpl_array_get_mean (mean_spectra[tel]) * nwave_ft);
	CPLCHECK_MSG ("Cannot compute mean spectra");
  }

  /* Compute the photometric normalisation for the FT:
   * F1F2 = Fsmooth1(t) * Fsmooth2(t) * mean_spectra1(lbd) * mean_spectra2(lbd) */
  gravi_table_new_column_array (vis_FT, "F1F2", "e^2", CPL_TYPE_DOUBLE, nwave_ft);
  cpl_array ** f1f2_ft = cpl_table_get_data_array (vis_FT, "F1F2");
		
  CPLCHECK_MSG ("Cannot create columns");

  /* Loop on base */
  for (cpl_size base = 0; base < nbase; base++) {
	int t0 = GRAVI_BASE_TEL[base][0];
	int t1 = GRAVI_BASE_TEL[base][1];
	for (cpl_size n = 0; n < nrow_ft; n ++) {
	  f1f2_ft[n*nbase+base] = cpl_array_duplicate (mean_spectra[t0]);
	  cpl_array_multiply (f1f2_ft[n*nbase+base], mean_spectra[t1]);
	  cpl_array_multiply_scalar (f1f2_ft[n*nbase+base], total_flux_ft[n*ntel+t0] * total_flux_ft[n*ntel+t1]);
	  // f1f2_ft[n*nbase+base] = cpl_array_duplicate (flux_ft[n*ntel+t0]);
	  // cpl_array_multiply (f1f2_ft[n*nbase+base], flux_ft[n*ntel+t1]);
	}
  }
  
  FREELOOP (cpl_array_delete, mean_spectra, 4);
  
  gravi_msg_function_exit(1);
  return CPL_ERROR_NONE;
}

/* -------------------------------------------------------------------------- */
/**
 * @brief Compute the self-reference phase for each FT DIT.
 * 
 * @param vis_FT:   input/output OI_VIS table of the FT
 *
 * The reference phase is computed for each FT DIT and saved in a newly
 * created column PHASE_REF in the vis_FT table. It is a running mean of 
 * few DITs, performed independently for each spectral channel.
 */
/* -------------------------------------------------------------------------- */

cpl_error_code gravi_vis_create_phaseref_ft (cpl_table * vis_FT)
{
  gravi_msg_function_start(1);
  cpl_ensure_code (vis_FT,  CPL_ERROR_NULL_INPUT);
  
  /* Get the FT data */
  cpl_size nbase = 6;
  cpl_size nrow_ft  = cpl_table_get_nrow (vis_FT) / nbase;
  cpl_size nwave_ft = cpl_table_get_column_depth (vis_FT, "VISDATA");
  cpl_array **  visData_ft = cpl_table_get_data_array (vis_FT, "VISDATA");
  
  CPLCHECK_MSG ("Cannot get data");
  
  /* Create the column */
  gravi_table_new_column_array (vis_FT, "SELF_REF", "rad", CPL_TYPE_DOUBLE, nwave_ft);
  cpl_array ** phaseref_ft = cpl_table_get_data_array (vis_FT, "SELF_REF");
  
  CPLCHECK_MSG ("Cannot create columns");
  
  /* Compute the PHASE_REF for the FT as a weighted
   * mean of -3..+3 frames */
  for (cpl_size base = 0; base < nbase; base ++) {

	/* Loop on running frames */
	for (cpl_size n = 3; n < nrow_ft - 3; n ++) {
	  phaseref_ft[n*nbase+base] = cpl_array_cast (visData_ft[(n-1)*nbase+base], CPL_TYPE_DOUBLE_COMPLEX);
	  cpl_array_add (phaseref_ft[n*nbase+base], visData_ft[(n+1)*nbase+base]);
	  cpl_array_multiply_scalar (phaseref_ft[n*nbase+base], 2.0);
	  cpl_array_add (phaseref_ft[n*nbase+base], visData_ft[(n-2)*nbase+base]);
	  cpl_array_add (phaseref_ft[n*nbase+base], visData_ft[(n+2)*nbase+base]);
	  cpl_array_multiply_scalar (phaseref_ft[n*nbase+base], 2.0);
	  cpl_array_add (phaseref_ft[n*nbase+base], visData_ft[(n-3)*nbase+base]);
	  cpl_array_add (phaseref_ft[n*nbase+base], visData_ft[(n+3)*nbase+base]);
	  // cpl_array_fill_window_complex (phaseref_ft[n*nbase+base], 0, nwave_ft, cpl_array_get_mean_complex (phaseref_ft[n*nbase+base]));
	  cpl_array_arg (phaseref_ft[n*nbase+base]);
	  cpl_array_multiply_scalar (phaseref_ft[n*nbase+base], -1.0);
	  
	  CPLCHECK_MSG("Cannot compute the PHASE_REF for the FT");
	}
	phaseref_ft[0*nbase+base] = cpl_array_duplicate (phaseref_ft[3*nbase+base]);
	phaseref_ft[1*nbase+base] = cpl_array_duplicate (phaseref_ft[3*nbase+base]);
	phaseref_ft[2*nbase+base] = cpl_array_duplicate (phaseref_ft[3*nbase+base]);
	phaseref_ft[(nrow_ft-3)*nbase+base] = cpl_array_duplicate (phaseref_ft[(nrow_ft-4)*nbase+base]);
	phaseref_ft[(nrow_ft-2)*nbase+base] = cpl_array_duplicate (phaseref_ft[(nrow_ft-4)*nbase+base]);
	phaseref_ft[(nrow_ft-1)*nbase+base] = cpl_array_duplicate (phaseref_ft[(nrow_ft-4)*nbase+base]);
  }
  /* End loop on FT rows and base */

  cpl_msg_warning (cpl_func,"Change of PHASE_REF_FT common channel new... to be decided !!");

  gravi_msg_function_exit(1);
  return CPL_ERROR_NONE;
}

/* -------------------------------------------------------------------------- */
/**
 * @brief Compute the averaged ACQ signal for each SC DIT per base
 * 
 * @param vis_SC:   input/output OI_VIS table of the SC
 * @param vis_ACQ:  input OI_VIS_ACQ table
 *
 * The averaged quantities are stored in new columns in the vis_SC table. The
 * ACQ_CAM signals are averaged with flat weighting inside each SC DIT. The
 * synchronisation info shall already be computed.
 */
/* -------------------------------------------------------------------------- */

cpl_error_code gravi_vis_create_acq_sc (cpl_table * vis_SC,
                                        cpl_table * vis_ACQ)
{
  gravi_msg_function_start(1);
  cpl_ensure_code (vis_SC,  CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (vis_ACQ, CPL_ERROR_NULL_INPUT);
  
  cpl_size nbase = 6, ntel = 4;
  cpl_size nrow_sc  = cpl_table_get_nrow (vis_SC) / nbase;

  /* Get SC data */
  int * first = cpl_table_get_data_int (vis_SC, "FIRST_ACQ");
  int * last  = cpl_table_get_data_int (vis_SC, "LAST_ACQ");
  CPLCHECK_MSG("Cannot get data");

  /* Get ACQ data */
  int * pup_n = cpl_table_get_data_int (vis_ACQ, "PUPIL_NSPOT");
  double * pup_x = cpl_table_get_data_double (vis_ACQ, "PUPIL_X");
  double * pup_y = cpl_table_get_data_double (vis_ACQ, "PUPIL_Y");
  double * pup_z = cpl_table_get_data_double (vis_ACQ, "PUPIL_Z");
  CPLCHECK_MSG("Cannot get direct pointer to data");

  /* New columns */
  gravi_table_new_column (vis_SC, "PUPIL_X", NULL, CPL_TYPE_DOUBLE);
  double * pup_x_sc = cpl_table_get_data_double (vis_SC, "PUPIL_X");
  gravi_table_new_column (vis_SC, "PUPIL_Y", NULL, CPL_TYPE_DOUBLE);
  double * pup_y_sc = cpl_table_get_data_double (vis_SC, "PUPIL_Y");
  gravi_table_new_column (vis_SC, "PUPIL_Z", NULL, CPL_TYPE_DOUBLE);
  double * pup_z_sc = cpl_table_get_data_double (vis_SC, "PUPIL_Z");

  /* Loop on base and rows */
  for (cpl_size base = 0; base < nbase; base++) {
	for (cpl_size row_sc = 0; row_sc < nrow_sc; row_sc ++) {
	  cpl_size nsc = row_sc * nbase + base;
  
	  /* Sum over synch ACQ frames, only valid frames */
      cpl_size nframe = 0;
	  for (cpl_size row = first[nsc] ; row < last[nsc]; row++) {
          cpl_size row0 = row * ntel + GRAVI_BASE_TEL[base][0];
          cpl_size row1 = row * ntel + GRAVI_BASE_TEL[base][1];

          if (pup_n[row0] != 0 && pup_n[row1] !=0 ) {
              pup_x_sc[nsc] += pup_x[row0] - pup_x[row1];
              pup_y_sc[nsc] += pup_y[row0] - pup_y[row1];
              pup_z_sc[nsc] += pup_z[row0] - pup_z[row1];
              nframe ++;
          }
          
          CPLCHECK_MSG ("Fail to integrate the ACQ frames");
      }
      
	  /* Normalize the means  (if nframe == 0, values are zero) */
	  if (nframe != 0 ){
          pup_x_sc[nsc] /= (double)nframe;
          pup_y_sc[nsc] /= (double)nframe;
          pup_z_sc[nsc] /= (double)nframe;
      }
      
	} /* End loop on SC frames */
  }/* End loop on bases */

  gravi_msg_function_exit(1);
  return CPL_ERROR_NONE;
}

/* -------------------------------------------------------------------------- */
/**
 * @brief Compute the averaged ACQ signal for each SC DIT per beam
 * 
 * @param flux_SC:  input/output OI_VIS table of the SC
 * @param vis_ACQ:  input OI_VIS_ACQ table
 *
 * The averaged quantities are stored in new columns in the vis_SC table. The
 * ACQ_CAM signals are averaged with flat weighting inside each SC DIT. The
 * synchronisation info shall already be computed.
 */
/* -------------------------------------------------------------------------- */

cpl_error_code gravi_flux_create_acq_sc (cpl_table * flux_SC,
                                         cpl_table * vis_ACQ)
{
  gravi_msg_function_start(1);
  cpl_ensure_code (flux_SC,  CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (vis_ACQ, CPL_ERROR_NULL_INPUT);
  
  cpl_size ntel = 4;
  cpl_size nrow_sc  = cpl_table_get_nrow (flux_SC) / ntel;

  /* Get SC data */
  int * first = cpl_table_get_data_int (flux_SC, "FIRST_ACQ");
  int * last  = cpl_table_get_data_int (flux_SC, "LAST_ACQ");
  CPLCHECK_MSG ("Cannot get data");

  /* Get ACQ data */
  int * pup_n = cpl_table_get_data_int (vis_ACQ, "PUPIL_NSPOT");
  double * pup_x = cpl_table_get_data_double (vis_ACQ, "PUPIL_X");
  double * pup_y = cpl_table_get_data_double (vis_ACQ, "PUPIL_Y");
  double * pup_z = cpl_table_get_data_double (vis_ACQ, "PUPIL_Z");
  CPLCHECK_MSG ("Cannot get direct pointer to data");

  /* New columns -- filled with zero */
  gravi_table_new_column (flux_SC, "PUPIL_X", NULL, CPL_TYPE_DOUBLE);
  double * pup_x_sc = cpl_table_get_data_double (flux_SC, "PUPIL_X");
  gravi_table_new_column (flux_SC, "PUPIL_Y", NULL, CPL_TYPE_DOUBLE);
  double * pup_y_sc = cpl_table_get_data_double (flux_SC, "PUPIL_Y");
  gravi_table_new_column (flux_SC, "PUPIL_Z", NULL, CPL_TYPE_DOUBLE);
  double * pup_z_sc = cpl_table_get_data_double (flux_SC, "PUPIL_Z");

  /* Loop on base and rows */
  for (cpl_size tel = 0; tel < ntel; tel++) {
	for (cpl_size row_sc = 0; row_sc < nrow_sc; row_sc ++) {
	  cpl_size nsc = row_sc * ntel + tel;
  
	  /* Sum over synch ACQ frames, only valid frames */
      cpl_size nframe = 0;
	  for (cpl_size row = first[nsc] ; row < last[nsc]; row++) {
          cpl_size row0 = row * ntel + tel;

          if (pup_n[row0] != 0) {
              pup_x_sc[nsc] += pup_x[row0];
              pup_y_sc[nsc] += pup_y[row0];
              pup_z_sc[nsc] += pup_z[row0];
              nframe ++;
          }
          
          CPLCHECK_MSG ("Fail to integrate the ACQ frames");
      }
      
	  /* Normalize the means  (if nframe == 0, values are zero) */
	  if (nframe != 0 ){
          pup_x_sc[nsc] /= (double)nframe;
          pup_y_sc[nsc] /= (double)nframe;
          pup_z_sc[nsc] /= (double)nframe;
      }
      
	} /* End loop on SC frames */
  }/* End loop on bases */

  gravi_msg_function_exit(1);
  return CPL_ERROR_NONE;
}

/* -------------------------------------------------------------------------- */
/**
 * @brief Compute the averaged MET signal for each SC DIT per base
 * 
 * @param vis_SC:   input/output OI_VIS table of the SC
 * @param vis_MET:  input OI_VIS_MET table
 *
 * The averaged quantities are stored in new columns in the vis_SC table. The
 * metrology signals are averaged with flat weighting inside each SC DIT. The
 * synchronisation info shall already be computed.
 */
/* -------------------------------------------------------------------------- */

cpl_error_code gravi_vis_create_met_sc (cpl_table * vis_SC, cpl_table * vis_MET,
                                        cpl_table * wave_table)
{
  gravi_msg_function_start(1);
  cpl_ensure_code (vis_SC,  CPL_ERROR_NULL_INPUT);
    cpl_ensure_code (vis_MET, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code (wave_table, CPL_ERROR_NULL_INPUT);

  cpl_size nbase = 6, ndiode = 4, ntel = 4;
  cpl_size nrow_sc  = cpl_table_get_nrow (vis_SC) / nbase;

   /* get wavenummer from wave table */
    
    cpl_size nwave_sc = cpl_table_get_column_depth (vis_SC, "VISDATA");
    double * twopi_wavenumber_sc  = cpl_malloc (nwave_sc * sizeof(double));
    for (cpl_size wave = 0; wave < nwave_sc ; wave ++) {
        twopi_wavenumber_sc[wave] = CPL_MATH_2PI / cpl_table_get (wave_table, "EFF_WAVE", wave, NULL);
    }
    
  /* Get SC data */
  int * first_met = cpl_table_get_data_int (vis_SC, "FIRST_MET");
  int * last_met  = cpl_table_get_data_int (vis_SC, "LAST_MET");
  
  CPLCHECK_MSG("Cannot get data");

  /* Get MET data */
  double * opd_met_fc           = cpl_table_get_data_double (vis_MET, "OPD_FC");
  cpl_array ** opd_met_tel      = cpl_table_get_data_array (vis_MET, "OPD_TEL");

  double * opd_met_fc_corr        = cpl_table_get_data_double (vis_MET, "OPD_FC_CORR");
  double * opd_met_telfc_mcorr    = cpl_table_get_data_double (vis_MET, "OPD_TELFC_MCORR");
  cpl_array ** opd_met_telfc_corr = cpl_table_get_data_array (vis_MET, "OPD_TELFC_CORR");

  CPLCHECK_MSG("Cannot get direct pointer to data");

  /* New columns */
  gravi_table_new_column_array (vis_SC, "PHASE_MET_TELFC", "rad", CPL_TYPE_DOUBLE, nwave_sc);
  cpl_array ** phase_metdit_telfc = cpl_table_get_data_array (vis_SC, "PHASE_MET_TELFC");

  gravi_table_new_column (vis_SC, "OPD_MET_FC", "m", CPL_TYPE_DOUBLE);
  double * opd_metdit_fc = cpl_table_get_data_double (vis_SC, "OPD_MET_FC");
    
  gravi_table_new_column_array (vis_SC, "OPD_MET_TEL", "m", CPL_TYPE_DOUBLE, ndiode);
  cpl_array ** opd_metdit_tel = cpl_table_get_data_array (vis_SC, "OPD_MET_TEL");

  gravi_table_new_column (vis_SC, "OPD_MET_FC_CORR", "m", CPL_TYPE_DOUBLE);
  double * opd_metdit_fc_corr = cpl_table_get_data_double (vis_SC, "OPD_MET_FC_CORR");

  gravi_table_new_column (vis_SC, "OPD_MET_TELFC_MCORR", "m", CPL_TYPE_DOUBLE);
  double * opd_metdit_telfc_mcorr = cpl_table_get_data_double (vis_SC, "OPD_MET_TELFC_MCORR");

  gravi_table_new_column_array (vis_SC, "OPD_MET_TELFC_CORR", "m", CPL_TYPE_DOUBLE, ndiode);
  cpl_array ** opd_metdit_telfc_corr = cpl_table_get_data_array (vis_SC, "OPD_MET_TELFC_CORR");
    
  CPLCHECK_MSG("Cannot create columns");

  /* Loop on base and rows */
  for (cpl_size base = 0; base < nbase; base++) {
	for (cpl_size row_sc = 0; row_sc < nrow_sc; row_sc ++) {
	  cpl_size nsc = row_sc * nbase + base;
        
      /* init cpl arrays */
	  opd_metdit_tel[nsc]      = gravi_array_init_double (ndiode, 0.0);
      opd_metdit_telfc_corr[nsc] = gravi_array_init_double (ndiode, 0.0);
      cpl_array * phasor_metdit_telfc = gravi_array_init_double_complex (nwave_sc, 0.0+I*0.0);
        
	  /* Sum over synch MET frames */
	  for (cpl_size row_met = first_met[nsc] ; row_met < last_met[nsc]; row_met++) {
	    cpl_size nmet0 = row_met * ntel + GRAVI_BASE_TEL[base][0];
	    cpl_size nmet1 = row_met * ntel + GRAVI_BASE_TEL[base][1];
          
        /* compute phasor of astrometric quantity */
          double opd_astro = opd_met_fc_corr[nmet0] - opd_met_fc_corr[nmet1] + opd_met_telfc_mcorr[nmet0] - opd_met_telfc_mcorr[nmet1];
          for (cpl_size wave = 0; wave < nwave_sc ; wave ++) {
              cpl_array_set_complex(phasor_metdit_telfc,wave,
                                    cpl_array_get_complex(phasor_metdit_telfc,wave,NULL)+
                                      cexp(opd_astro*I*twopi_wavenumber_sc[wave]));
                                                          };
          
	    /* Mean OPD_FC_CORR and OPD_TELFC_MCORR for each BASELINE */
	    opd_metdit_fc_corr[nsc] += opd_met_fc_corr[nmet0] - opd_met_fc_corr[nmet1];
	    opd_metdit_telfc_mcorr[nsc] += opd_met_telfc_mcorr[nmet0] - opd_met_telfc_mcorr[nmet1];

	    /* Mean OPD_TELFC_CORR for each BASELINE and diode */
	    cpl_array_add (opd_metdit_telfc_corr[nsc], opd_met_telfc_corr[nmet0]);
	    cpl_array_subtract (opd_metdit_telfc_corr[nsc], opd_met_telfc_corr[nmet1]);

	    /* Mean OPD_MET at Telescope (each diode) */
	    cpl_array_add (opd_metdit_tel[nsc], opd_met_tel[nmet0]);
	    cpl_array_subtract (opd_metdit_tel[nsc], opd_met_tel[nmet1]);
	    
	    /* Mean OPD_MET_FC at Beam Combiner */
	    opd_metdit_fc[nsc] += opd_met_fc[nmet0] - opd_met_fc[nmet1];

	    CPLCHECK_MSG ("Fail to integrate the metrology");
	  }
	  
	  /* Normalize the means  (if nframe == 0, values are zero) */
	  cpl_size nframe = last_met[nsc] - first_met[nsc];
	  if (nframe != 0 ){
                
	    opd_metdit_fc_corr[nsc] /= nframe;
	    opd_metdit_telfc_mcorr[nsc] /= nframe;
	    cpl_array_divide_scalar (opd_metdit_telfc_corr[nsc], (double)nframe);
	    cpl_array_divide_scalar (opd_metdit_tel[nsc], (double)nframe);
	    opd_metdit_fc[nsc]  /= nframe;
          
        /* get the astro phase by taking the argument of astro phasor */
        cpl_array_arg (phasor_metdit_telfc);
        phase_metdit_telfc[nsc]=cpl_array_cast(phasor_metdit_telfc, CPL_TYPE_DOUBLE);
	  }
      cpl_array_delete(phasor_metdit_telfc);
	  CPLCHECK_MSG ("Fail to compute metrology per base from metrology per tel");
	  
	} /* End loop on SC frames */
  }/* End loop on bases */

  /* Compute the information comming from VIS_ACQ
   * camera... through the VIS_MET */
  if (cpl_table_has_column (vis_MET,"FIELD_FIBER_DX")) {
      
      double * fdx_met = cpl_table_get_data_double (vis_MET, "FIELD_FIBER_DX");
      double * fdy_met = cpl_table_get_data_double (vis_MET, "FIELD_FIBER_DY");

      gravi_table_new_column (vis_SC, "FIELD_FIBER_DX", "pix", CPL_TYPE_DOUBLE);
      double * fdx_metdit = cpl_table_get_data_double (vis_SC, "FIELD_FIBER_DX");

      gravi_table_new_column (vis_SC, "FIELD_FIBER_DY", "pix", CPL_TYPE_DOUBLE);
      double * fdy_metdit = cpl_table_get_data_double (vis_SC, "FIELD_FIBER_DY");

      /* Loop on base and rows */
      for (cpl_size base = 0; base < nbase; base++) {
          for (cpl_size row_sc = 0; row_sc < nrow_sc; row_sc ++) {
              cpl_size nsc = row_sc * nbase + base;
              
              /* Sum over synch MET frames */
              for (cpl_size row_met = first_met[nsc] ; row_met < last_met[nsc]; row_met++) {
                  cpl_size nmet0 = row_met * ntel + GRAVI_BASE_TEL[base][0];
                  cpl_size nmet1 = row_met * ntel + GRAVI_BASE_TEL[base][1];
                  
                  /* Mean FIELD_FIBER */
                  fdx_metdit[nsc] += fdx_met[nmet0] - fdx_met[nmet1];
                  fdy_metdit[nsc] += fdy_met[nmet0] - fdy_met[nmet1];
              }
              CPLCHECK_MSG ("Fail to compute metrology per base from metrology per tel");

              /* Normalize the means  (if nframe == 0, values are zero) */
              cpl_size nframe = last_met[nsc] - first_met[nsc];
              if (nframe != 0 ){
                  fdx_metdit[nsc]  /= nframe;
                  fdy_metdit[nsc]  /= nframe;
              }
              
          } /* End loop on SC frames */
      }/* End loop on bases */
  }

  /* Compute mean OPD_MET_PUPIL, OPD_MET_PUPIL_STDDEV and OPD_MET_TTPUP 
  only if OPD_PUPIL exists (ACQ option)*/
  if (cpl_table_has_column (vis_MET,"OPD_PUPIL") && cpl_table_has_column (vis_MET,"OPD_TTPUP")) {

	  double * opd_met_pupil = cpl_table_get_data_double (vis_MET, "OPD_PUPIL");
	  double * opd_met_ttpup          = cpl_table_get_data_double (vis_MET, "OPD_TTPUP");

	  gravi_table_new_column (vis_SC, "OPD_MET_PUPIL", "m", CPL_TYPE_DOUBLE);
	  double * opd_metdit_pupil = cpl_table_get_data_double (vis_SC, "OPD_MET_PUPIL");

	  gravi_table_new_column (vis_SC, "OPD_MET_PUPIL_STDDEV", "m", CPL_TYPE_DOUBLE);
	  double * opd_metdit_pupil_stddev = cpl_table_get_data_double (vis_SC, "OPD_MET_PUPIL_STDDEV");

	  gravi_table_new_column (vis_SC, "OPD_MET_TTPUP", "m", CPL_TYPE_DOUBLE);
	  double * opd_metdit_ttpup = cpl_table_get_data_double (vis_SC, "OPD_MET_TTPUP");

	  CPLCHECK_MSG("Cannot create columns");

	  /* Loop on base and rows */
	  for (cpl_size base = 0; base < nbase; base++) {
		for (cpl_size row_sc = 0; row_sc < nrow_sc; row_sc ++) {
		  cpl_size nsc = row_sc * nbase + base;

		  /* Sum over synch MET frames */
		  for (cpl_size row_met = first_met[nsc] ; row_met < last_met[nsc]; row_met++) {
		    cpl_size nmet0 = row_met * ntel + GRAVI_BASE_TEL[base][0];
		    cpl_size nmet1 = row_met * ntel + GRAVI_BASE_TEL[base][1];

		    /* Mean OPD_PUPIL */
		    opd_metdit_pupil[nsc] += opd_met_pupil[nmet0] - opd_met_pupil[nmet1];

		    /* Mean OPD_TTPUP */
		    opd_metdit_ttpup[nsc] += opd_met_ttpup[nmet0] - opd_met_ttpup[nmet1];

		    CPLCHECK_MSG ("Fail to integrate the metrology");
		  }

		  /* Normalize the means  (if nframe == 0, values are zero) */
		  cpl_size nframe = last_met[nsc] - first_met[nsc];
		  if (nframe != 0 ){
		    opd_metdit_pupil[nsc]  /= nframe;
		    opd_metdit_ttpup[nsc]  /= nframe;
		  }

		  /* Sum over synch MET frames again, to compute STDDEV */
		  for (cpl_size row_met = first_met[nsc] ; row_met < last_met[nsc]; row_met++) {
		    cpl_size nmet0 = row_met * ntel + GRAVI_BASE_TEL[base][0];
		    cpl_size nmet1 = row_met * ntel + GRAVI_BASE_TEL[base][1];

		    /* Mean OPD_PUPIL_STDDEV (still sum variance here, sqrt below) */
		    double tmp = opd_met_pupil[nmet0] - opd_met_pupil[nmet1] - opd_metdit_pupil[nsc];
		    opd_metdit_pupil_stddev[nsc] += tmp * tmp;

		    CPLCHECK_MSG ("Fail to integrate the metrology");
		  }

		  /* Normalize the STDDEV  (if nframe <= 1, values are zero) */
		  if (nframe > 1 ){
		    opd_metdit_pupil_stddev[nsc]  = sqrt (opd_metdit_pupil_stddev[nsc]  / (nframe-1));
		  }

		  

          CPLCHECK_MSG ("Fail to compute metrology per base from metrology per tel");

		} /* End loop on SC frames */
	  } /* End loop on bases */

  }


  FREE (cpl_free, twopi_wavenumber_sc);
  gravi_msg_function_exit(1);
  return CPL_ERROR_NONE;
}


/* -------------------------------------------------------------------------- */
/**
 * @brief Compute the resampled MET signal for each FT DIT per base
 * 
 * @param vis_FT:   input/output OI_VIS table of the FT
 * @param vis_MET:  input OI_VIS_MET table
 * @param dit_ft:   the FT DIT  (exposure being TIME -> TIME+PERIOD
 *
 * The resampled quantities are stored in new columns in the vis_FT table.
 * The routine doesn't need the synchronisation columns. The PERIOD of FT
 * is computed as TIME[sample1] - TIME[sample0]. The TIME of the MET table
 * is supposed to be the same for all beam.
 *
 * If one or several MET sample are inside the FT DIT, they are averaged.
 * If no MET samples are inside the FT DIT, the metrology is interpolated
 * linearly at the time of FT.
 *
 * Create table is PHASE_MET_FC (fiber coupler) in [rad]
 */
/* -------------------------------------------------------------------------- */

cpl_error_code gravi_vis_create_met_ft (cpl_table * vis_FT, cpl_table * vis_MET)
{
  gravi_msg_function_start(1);
  cpl_ensure_code (vis_FT,  CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (vis_MET, CPL_ERROR_NULL_INPUT);

  cpl_size nbase = 6, ntel = 4;
  cpl_size nrow_ft  = cpl_table_get_nrow (vis_FT) / nbase;
  cpl_size nrow_met = cpl_table_get_nrow (vis_MET) / ntel;

  CPLCHECK_MSG("Cannot get data");

  /* Get the FT period in table [us] */
  double periode_ft = cpl_table_get (vis_FT, "TIME", nbase, NULL) -
                      cpl_table_get (vis_FT, "TIME", 0, NULL);
  cpl_msg_info (cpl_func, "PERIOD FT = %g [us]", periode_ft);

  /* Get MET data */
  double * phase_met_fc = cpl_table_get_data_double (vis_MET, "PHASE_FC_DRS");

  CPLCHECK_MSG("Cannot get direct pointer to data");

  /* New columns */
  gravi_table_new_column (vis_FT, "PHASE_MET_FC", "rad", CPL_TYPE_DOUBLE);
  double * phase_metdit_fc = cpl_table_get_data_double (vis_FT, "PHASE_MET_FC");
		
  CPLCHECK_MSG("Cannot create columns");

  /* Loop on base and rows */
  for (cpl_size base = 0; base < nbase; base++) {
      int tel0 = GRAVI_BASE_TEL[base][0];
      int tel1 = GRAVI_BASE_TEL[base][1];

      for (cpl_size last_met = 0, first_met = 0, row_ft = 0; row_ft < nrow_ft; row_ft ++) {
          cpl_size nft = row_ft * nbase + base;
          double time_ft = cpl_table_get (vis_FT, "TIME", row_ft*nbase+base, NULL);

          /* 
           * FIXME: to respect previous implementation in gravi_wave
           * while ((time_met  < (time_ft + exptime_ft))){
           * we use this shiffted window time_ft -> time_ft + periode_ft
           */

          /* First sample of MET (assume same DIT for all beam in vis_MET) */
          first_met = CPL_MAX (CPL_MIN (last_met - 5, nrow_met - 1), 0);
          while ((cpl_table_get (vis_MET, "TIME", first_met*ntel, NULL) <= (time_ft + 0.0))) {
              first_met++;
              if (first_met == nrow_met) break;
          }
          CPLCHECK_MSG ("Cannot get first");

          /* Last sample of MET (assume same DIT for all beam in vis_MET) */
          last_met = CPL_MAX (CPL_MIN (first_met - 1, nrow_met - 1), 0);
          while ((cpl_table_get (vis_MET, "TIME", last_met*ntel, NULL) < (time_ft + periode_ft))) {
              last_met++;
              if (last_met == nrow_met) break;
          }
          CPLCHECK_MSG ("Cannot get last");

          /* For first few FT samples, we use the first MET if none found */
          if (row_ft < 5 && last_met == 0) last_met  = 1;
          
          /* For last few FT samples, we use the last MET if none found 
           * also avoid going outside the table */
          if (row_ft > nrow_ft-5 && first_met == nrow_met) first_met = nrow_met - 1;
          if (row_ft > nrow_ft-5 && last_met  == nrow_met) last_met  = nrow_met;

          /* Check if enough data */
          if ( last_met == 0 || last_met > nrow_met ) {
              return cpl_error_set_message (cpl_func,CPL_ERROR_ILLEGAL_INPUT,
                                            "Not enough MET data to synchronise "
                                            "with FT DIT %lli over %lli", row_ft+1, nrow_ft);
          }
          
          if ( last_met - first_met > 0 ) {
              /* If at least one MET samples inside, we average 
               * FIXME: maybe we better always interpolate ?? */
              
              for (cpl_size rin_met = first_met; rin_met < last_met; rin_met ++)
                  phase_metdit_fc[nft] += phase_met_fc[rin_met*ntel+tel0] - phase_met_fc[rin_met*ntel+tel1];
              
              phase_metdit_fc[nft] /= (last_met - first_met);
              
          } else {
              /* If no MET inside sample, we interpolate linear */
              cpl_size rowa_met = first_met-1;
              cpl_size rowb_met = last_met;
              
              double phia_met = phase_met_fc[rowa_met*ntel+tel0] - phase_met_fc[rowa_met*ntel+tel1];
              double phib_met = phase_met_fc[rowb_met*ntel+tel0] - phase_met_fc[rowb_met*ntel+tel1];
              double timea_met = cpl_table_get (vis_MET, "TIME",rowa_met*ntel, NULL);
              double timeb_met = cpl_table_get (vis_MET, "TIME",rowb_met*ntel, NULL);
              
              phase_metdit_fc[nft] = phia_met +
                  (phib_met - phia_met) * (time_ft - timea_met) / (timeb_met - timea_met);
              CPLCHECK_MSG ("Cannot interpolate");
          }
          
          CPLCHECK_MSG ("Fail to compute metrology per FT base from metrology per tel");
      } /* End loop on FT frames */
      
  }/* End loop on bases */

  gravi_msg_function_exit(1);
  return CPL_ERROR_NONE;
}


/* -------------------------------------------------------------------------- */
/**
 * @brief Compute the resampled SC signal for each FT DIT per base
 * 
 * @param vis_FT:   input/output OI_VIS table of the FT
 * @param vis_SC:   input OI_VIS table of SC
 * @param dit_sc:   the SC DIT  (exposure being TIME-DIT/2 -> TIME+DIT/2)
 *
 * The resampled quantities are stored in new columns in the vis_FT table.
 * The routine doesn't need the synchronisations columns.
 * For all FT samples inside an SC DIT, the routine filled the PHASE_SC column
 * with the value of PHASE column from vis_SC. Samples outside an SC DIT are
 * filled with 0.0
 *
 * Create table PHASE_SC [rad]
 */
/* -------------------------------------------------------------------------- */

cpl_error_code gravi_vis_create_opdsc_ft (cpl_table * vis_FT,
                                          cpl_table * vis_SC,
                                          double dit_sc)
{
  gravi_msg_function_start(1);
  cpl_ensure_code (vis_FT, CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (vis_SC, CPL_ERROR_NULL_INPUT);

  cpl_size nbase = 6;
  cpl_size nrow_ft = cpl_table_get_nrow (vis_FT) / nbase;
  cpl_size nrow_sc = cpl_table_get_nrow (vis_SC) / nbase;

  CPLCHECK_MSG("Cannot get data");

  /* Get OPD data */
  double * phase_sc = cpl_table_get_data_double (vis_SC, "OPD");

  CPLCHECK_MSG("Cannot get direct pointer to data");

  /* New columns */
  gravi_table_new_column (vis_FT, "OPD_SC", "rad", CPL_TYPE_DOUBLE);
  double * phase_scdit = cpl_table_get_data_double (vis_FT, "OPD_SC");
		
  CPLCHECK_MSG("Cannot create columns");

  /* Loop on base and rows */
  for (cpl_size base = 0; base < nbase; base++) {
      cpl_size row_ft = 0;
      
      for (cpl_size row_sc = 0; row_sc < nrow_sc; row_sc ++) {
          cpl_size nsc = row_sc * nbase + base;
          double time_sc = cpl_table_get (vis_SC, "TIME", nsc, NULL);

          /* First sample of FT in the SC DIT */
          while ((cpl_table_get (vis_FT, "TIME", row_ft*nbase+base, NULL) < (time_sc - dit_sc/2))) {
              row_ft++;
              if (row_ft >= nrow_ft) break;
          }

          /*  While inside the SC DIT -- we fill the OPD_FT column */
          while ((cpl_table_get (vis_FT, "TIME", row_ft*nbase+base, NULL) < (time_sc + dit_sc/2))) {
              
              phase_scdit[row_ft*nbase+base] = phase_sc[nsc];
              
              row_ft++;
              if (row_ft >= nrow_ft) break;
          }

          if (row_ft >= nrow_ft) {
              return cpl_error_set_message (cpl_func,CPL_ERROR_ILLEGAL_INPUT,
                                            "Not enough FT data to synchronise "
                                            "with SC DIT %lli over %lli", row_sc+1, nrow_sc);
          }

      } /* End loop on SC frames */
      
  }/* End loop on bases */

  gravi_msg_function_exit(1);
  return CPL_ERROR_NONE;
}


/* -------------------------------------------------------------------------- */
/**
 * @brief Compute the averaged MET signal for each SC DIT per beam
 * 
 * @param flux_SC:   input/output OI_FLUX table of the SC
 * @param vis_MET:   input OI_VIS_MET table
 *
 * The averaged quantities are stored in new columns in the vis_SC table. The
 * metrology signals are averaged with flat weighting inside each SC DIT. The
 * synchronisation info shall already be computed.
 */
/* -------------------------------------------------------------------------- */

cpl_error_code gravi_flux_create_met_sc (cpl_table * flux_SC, cpl_table * vis_MET)
{
  gravi_msg_function_start(1);
  cpl_ensure_code (flux_SC,  CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (vis_MET, CPL_ERROR_NULL_INPUT);

  cpl_size ntel = 4, ndiode = 4;
  cpl_size nrow_sc  = cpl_table_get_nrow (flux_SC) / ntel;

  /* Get SC data */
  int * first_met = cpl_table_get_data_int (flux_SC, "FIRST_MET");
  int * last_met  = cpl_table_get_data_int (flux_SC, "LAST_MET");
  
  CPLCHECK_MSG("Cannot get data");

  /* Get MET data */
  double * opd_met_fc      = cpl_table_get_data_double (vis_MET, "OPD_FC");
  cpl_array ** opd_met_tel = cpl_table_get_data_array (vis_MET, "OPD_TEL");

  double * opd_met_fc_corr        = cpl_table_get_data_double (vis_MET, "OPD_FC_CORR");
  double * opd_met_telfc_mcorr    = cpl_table_get_data_double (vis_MET, "OPD_TELFC_MCORR");
  cpl_array ** opd_met_telfc_corr = cpl_table_get_data_array (vis_MET, "OPD_TELFC_CORR");
    double * phase_met_fc           = cpl_table_get_data_double (vis_MET, "PHASE_FC_DRS");
    double * phase_met_fcft         = cpl_table_get_data_double (vis_MET, "PHASE_FCFT_DRS");
    double * phase_met_fcsc         = cpl_table_get_data_double (vis_MET, "PHASE_FCSC_DRS");


  CPLCHECK_MSG("Cannot get direct pointer to data");

  /* New columns */
    gravi_table_new_column (flux_SC, "PHASE_MET_FC", "rad", CPL_TYPE_DOUBLE);
    double * phase_metdit_fc = cpl_table_get_data_double (flux_SC, "PHASE_MET_FC");
      
      gravi_table_new_column (flux_SC, "PHASE_MET_FCFT", "rad", CPL_TYPE_DOUBLE);
      double * phase_metdit_fcft = cpl_table_get_data_double (flux_SC, "PHASE_MET_FCFT");
      
      gravi_table_new_column (flux_SC, "PHASE_MET_FCSC", "rad", CPL_TYPE_DOUBLE);
      double * phase_metdit_fcsc = cpl_table_get_data_double (flux_SC, "PHASE_MET_FCSC");
    
  gravi_table_new_column (flux_SC, "OPD_MET_FC", "m", CPL_TYPE_DOUBLE);
  double * opd_metdit_fc = cpl_table_get_data_double (flux_SC, "OPD_MET_FC");

  gravi_table_new_column_array (flux_SC, "OPD_MET_TEL", "m", CPL_TYPE_DOUBLE, ndiode);
  cpl_array ** opd_metdit_tel = cpl_table_get_data_array (flux_SC, "OPD_MET_TEL");

  gravi_table_new_column (flux_SC, "OPD_MET_FC_CORR", "m", CPL_TYPE_DOUBLE);
  double * opd_metdit_fc_corr = cpl_table_get_data_double (flux_SC, "OPD_MET_FC_CORR");

  gravi_table_new_column (flux_SC, "OPD_MET_TELFC_MCORR", "m", CPL_TYPE_DOUBLE);
  double * opd_metdit_telfc_mcorr = cpl_table_get_data_double (flux_SC, "OPD_MET_TELFC_MCORR");

  gravi_table_new_column_array (flux_SC, "OPD_MET_TELFC_CORR", "m", CPL_TYPE_DOUBLE, ndiode);
  cpl_array ** opd_metdit_telfc_corr = cpl_table_get_data_array (flux_SC, "OPD_MET_TELFC_CORR");
  
  CPLCHECK_MSG("Cannot create columns");

  /* Loop on base and rows */
  for (cpl_size tel = 0; tel < ntel; tel++) {
	for (cpl_size row_sc = 0; row_sc < nrow_sc; row_sc ++) {
	  cpl_size nsc = row_sc * ntel + tel;

	  opd_metdit_tel[nsc] = gravi_array_init_double (ndiode, 0.0);
	  opd_metdit_telfc_corr[nsc] = gravi_array_init_double (ndiode, 0.0);

	  /* Sum over synch MET frames */
	  for (cpl_size row_met = first_met[nsc] ; row_met < last_met[nsc]; row_met++) {
        cpl_size nmet = row_met * ntel + tel;
        
		/* Mean OPD_MET_TEL at Telescope (each diode) */
        cpl_array_add (opd_metdit_tel[nsc], opd_met_tel[nmet]);
        
		/* Mean OPD_MET_FC at Beam Combiner */
		opd_metdit_fc[nsc] += opd_met_fc[nmet];

	    /* Mean OPD_FC_CORR and OPD_TELFC_MCORR for each BASELINE */
	    opd_metdit_fc_corr[nsc] += opd_met_fc_corr[nmet];
	    opd_metdit_telfc_mcorr[nsc] += opd_met_telfc_mcorr[nmet];
        
	    /* Mean OPD_TELFC_CORR for each BASELINE and diode */
	    cpl_array_add (opd_metdit_telfc_corr[nsc], opd_met_telfc_corr[nmet]);
          
        /* Mean PHASE_MET_FC at Beam Combiner */
        phase_metdit_fc[nsc] += phase_met_fc[nmet];
          
        /* Mean PHASE_MET_FCFT at Beam Combiner */
        phase_metdit_fcft[nsc] += phase_met_fcft[nmet];
          
        /* Mean PHASE_MET_FCSC at Beam Combiner */
        phase_metdit_fcsc[nsc] += phase_met_fcsc[nmet];
		
		CPLCHECK_MSG ("Fail to integrate the metrology");
	  }

	  /* Normalize the means  (if nframe == 0, values are zero) */
	  cpl_size nframe = last_met[nsc] - first_met[nsc];
	  if (nframe != 0 ){
          opd_metdit_fc_corr[nsc] /= nframe;
          opd_metdit_telfc_mcorr[nsc] /= nframe;
          phase_metdit_fc[nsc]  /= nframe;
          phase_metdit_fcft[nsc]  /= nframe;
          phase_metdit_fcsc[nsc]  /= nframe;
          cpl_array_divide_scalar (opd_metdit_telfc_corr[nsc], (double)nframe);
          
		  cpl_array_divide_scalar (opd_metdit_tel[nsc], (double)nframe);
		  opd_metdit_fc[nsc] /= nframe;
	  }
	  CPLCHECK_MSG ("Fail to integrate the metrology");

	} /* End loop on SC frames */
  }/* End loop on bases */

  
  /* Compute the information comming from VIS_ACQ
   * camera... through the VIS_MET */
  if (cpl_table_has_column (vis_MET,"FIELD_FIBER_DX")) {
      
      double * fdx_met = cpl_table_get_data_double (vis_MET, "FIELD_FIBER_DX");
      double * fdy_met = cpl_table_get_data_double (vis_MET, "FIELD_FIBER_DY");

      gravi_table_new_column (flux_SC, "FIELD_FIBER_DX", "pix", CPL_TYPE_DOUBLE);
      double * fdx_metdit = cpl_table_get_data_double (flux_SC, "FIELD_FIBER_DX");

      gravi_table_new_column (flux_SC, "FIELD_FIBER_DY", "pix", CPL_TYPE_DOUBLE);
      double * fdy_metdit = cpl_table_get_data_double (flux_SC, "FIELD_FIBER_DY");

      /* Loop on tel and rows */
      for (cpl_size tel = 0; tel < ntel; tel++) {
          for (cpl_size row_sc = 0; row_sc < nrow_sc; row_sc ++) {
              cpl_size nsc = row_sc * ntel + tel;
              
              /* Sum over synch MET frames */
              for (cpl_size row_met = first_met[nsc] ; row_met < last_met[nsc]; row_met++) {
                  cpl_size nmet0 = row_met * ntel + tel;
                  
                  /* Mean FIELD_FIBER */
                  fdx_metdit[nsc] += fdx_met[nmet0];
                  fdy_metdit[nsc] += fdy_met[nmet0];
              }
              CPLCHECK_MSG ("Fail to compute metrology per base from metrology per tel");

              /* Normalize the means  (if nframe == 0, values are zero) */
              cpl_size nframe = last_met[nsc] - first_met[nsc];
              if (nframe != 0 ){
                  fdx_metdit[nsc]  /= nframe;
                  fdy_metdit[nsc]  /= nframe;
              }
              
          } /* End loop on SC frames */
      }/* End loop on tels */
  }

  gravi_msg_function_exit(1);
  return CPL_ERROR_NONE;
}

/* -------------------------------------------------------------------------- */
/**
 * @brief Compute the averaged FDDL signal for each SC DIT per beam.
 * 
 * @param flux_SC:    input/output OI_FLUX table of the SC
 * @param vis_MET:    input OI_VIS_MET table
 * @param fddl_table: input FDDL table
 *
 * The averaged quantities are stored in new columns in the flux_SC table. The
 * MET and FDDL signals are averaged with flat weighting inside each SC DIT. The
 * synchronisation info shall already be computed.
 */
/* -------------------------------------------------------------------------- */

cpl_error_code gravi_flux_create_fddlpos_sc (cpl_table * flux_SC, cpl_table * fddl_table)
{
  gravi_msg_function_start(1);
  cpl_ensure_code (flux_SC,    CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (fddl_table, CPL_ERROR_NULL_INPUT);

  cpl_size ntel = 4;
  cpl_size nrow_sc  = cpl_table_get_nrow (flux_SC) / ntel;

  /* Get SC data */
  int * first_fddl = cpl_table_get_data_int (flux_SC, "FIRST_FDDL");
  int * last_fddl  = cpl_table_get_data_int (flux_SC, "LAST_FDDL");
  
  cpl_array ** ftpos  = cpl_table_get_data_array (fddl_table, "FT_POS");
  cpl_array ** scpos  = cpl_table_get_data_array (fddl_table, "SC_POS");
  cpl_array ** oplair = cpl_table_get_data_array (fddl_table, "OPL_AIR");
  
  CPLCHECK_MSG("Cannot get data");

  /* New columns */
  gravi_table_new_column (flux_SC, "FT_POS", "V", CPL_TYPE_DOUBLE);
  double * ftpos_flux = cpl_table_get_data_double (flux_SC, "FT_POS");

  gravi_table_new_column (flux_SC, "SC_POS", "V", CPL_TYPE_DOUBLE);
  double * scpos_flux = cpl_table_get_data_double (flux_SC, "SC_POS");
		
  gravi_table_new_column (flux_SC, "OPL_AIR", "m", CPL_TYPE_DOUBLE);
  double * oplair_flux = cpl_table_get_data_double (flux_SC, "OPL_AIR");

  CPLCHECK_MSG("Cannot create columns");
		
  /* Loop on tel and frames */
  for (cpl_size tel = 0; tel < ntel; tel++) {
	for (cpl_size row_sc = 0; row_sc < nrow_sc; row_sc ++) {
	  cpl_size nsc = row_sc * ntel + tel;

	  /* Mean FDDL and OPL_AIR during this frame for this tel */
      for (cpl_size row_fddl = first_fddl[nsc] ; row_fddl < last_fddl[nsc]; row_fddl++) {
			ftpos_flux[nsc]  += cpl_array_get (ftpos[row_fddl],  tel, NULL);
			scpos_flux[nsc]  += cpl_array_get (scpos[row_fddl],  tel, NULL);
			oplair_flux[nsc] += cpl_array_get (oplair[row_fddl], tel, NULL);
      }

      /* Normalise the means (if nframe == 0, values are zero) */
      cpl_size nframe = last_fddl[nsc] - first_fddl[nsc];
	  if (nframe != 0 ) {
          ftpos_flux[nsc]  /= nframe;
          scpos_flux[nsc]  /= nframe;
          oplair_flux[nsc] /= nframe;
      }

	} /* End loop on SC frames */
  }/* End loop on tels */

  gravi_msg_function_exit(1);
  return CPL_ERROR_NONE;
}

/* -------------------------------------------------------------------------- */
/**
 * @brief Compute total flux of each DIT for the SC and of the FT
 * 
 * @param flux_SC:    input/output OI_FLUX table of the SC
 * @param flux_FT:    input/output OI_FLUX table of the FT
 *
 * Create new columns for the total flux of the SC for each SC DIT, the total
 * flux of the FT for each FT DIT, and the total flux of the FT inside
 * each SC DIT. The synchronisation info shall already be computed.
 */
/* -------------------------------------------------------------------------- */

cpl_error_code gravi_flux_create_totalflux_sc (cpl_table * flux_SC, cpl_table * flux_FT)
{
  gravi_msg_function_start(1);
  cpl_ensure_code (flux_SC, CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (flux_FT, CPL_ERROR_NULL_INPUT);

  cpl_size ntel = 4;
  cpl_size nrow_sc  = cpl_table_get_nrow (flux_SC) / ntel;
  cpl_size nwave_sc = cpl_table_get_column_depth (flux_SC, "FLUX");
  cpl_size nwave_ft = cpl_table_get_column_depth (flux_FT, "FLUX");

  /* Get SC and FT data */
  int * first_ft = cpl_table_get_data_int (flux_SC, "FIRST_FT");
  int * last_ft  = cpl_table_get_data_int (flux_SC, "LAST_FT");
  
  cpl_array ** flag_sc = cpl_table_get_data_array (flux_SC, "FLAG");
  cpl_array ** flux_sc = cpl_table_get_data_array (flux_SC, "FLUX");
  cpl_array ** flux_ft = cpl_table_get_data_array (flux_FT, "FLUX");
  
  CPLCHECK_MSG("Cannot get data");

  /* New columns */
  gravi_table_new_column (flux_SC, "TOTALFLUX_SC", "e", CPL_TYPE_DOUBLE);
  double * total_flux_scdit = cpl_table_get_data_double (flux_SC, "TOTALFLUX_SC");

  gravi_table_new_column (flux_SC, "TOTALFLUX_FT", "e", CPL_TYPE_DOUBLE);
  double * total_flux_ftdit = cpl_table_get_data_double (flux_SC, "TOTALFLUX_FT");

  CPLCHECK_MSG("Cannot create columns");

  for (cpl_size tel = 0; tel < ntel; tel++) {
	for (cpl_size row_sc = 0; row_sc < nrow_sc; row_sc ++) {
	  cpl_size nsc = row_sc * ntel + tel;

	  /* Store the total flux of SC over the SC frame [e] */
	  total_flux_scdit[nsc] = 0.0;
	  int nvalid = 0;
	  for (int wave = 0; wave < nwave_sc; wave++) {
	    if (!cpl_array_get (flag_sc[nsc], wave, NULL)) {
	      total_flux_scdit[nsc] += cpl_array_get (flux_sc[nsc], wave, NULL);
	      nvalid++;
	    }
	  }

	  /* Normalise to replace the rejected values by the mean */
	  if (nvalid > 0) {
	    total_flux_scdit[nsc] *= (double)nwave_sc / (double)nvalid;
	  }
	  
	  /* Store the total flux of FT over the SC frame [e] */
	  for (cpl_size row_ft = first_ft[nsc] ; row_ft < last_ft[nsc]; row_ft++) {
		total_flux_ftdit[nsc] += cpl_array_get_mean (flux_ft[row_ft * ntel + tel]) * nwave_ft;
	  }
	  
	  CPLCHECK_MSG("Issue in the loop to average FT and SC flux per frame");

	} /* End loop on SC frames */
  }/* End loop on tels */

  gravi_msg_function_exit(1);
  return CPL_ERROR_NONE;
}

/* -------------------------------------------------------------------------- */
/**
 * @brief Compute the VFACTOR for each SC DIT based on real-time FT
 * 
 * @param vis_SC:        input/output OI_VIS table of the SC
 * @param wave_table_sc: wavelength table corresponding to vis_SC
 * @param vis_FT:        input/output OI_VIS table of the FT
 * @param wave_table_ft: wavelength table corresponding to vis_FT
 *
 * Create new columns in vis_SC with the VFACTOR of each SC DIT computed from
 * the ratio between coherence and incoherent integration of the FT fringes
 * inside each DIT. The synchronisation info shall already be computed.
 * The VFACTOR is measured for a pseudo-broad band light (averaging all FT
 * channels) and then extrapolated to all SC wavelength with a
 * theoretical model. The synchronisation info shall already be computed.
 */
/* -------------------------------------------------------------------------- */

cpl_error_code gravi_vis_create_vfactor_sc (cpl_table * vis_SC,
											cpl_table * wave_table_sc,
											cpl_table * vis_FT,
											cpl_table * wave_table_ft)
{
  gravi_msg_function_start(1);
  cpl_ensure_code (vis_SC,        CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (vis_FT,        CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (wave_table_sc, CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (wave_table_ft, CPL_ERROR_NULL_INPUT);

  /* Get variables */
  cpl_size nbase = 6;
  cpl_size nrow_sc  = cpl_table_get_nrow (vis_SC) / nbase;
  cpl_size nwave_sc = cpl_table_get_column_depth (vis_SC, "VISDATA");
  cpl_size nwave_ft = cpl_table_get_column_depth (vis_FT, "VISDATA");

  /* Get SC and FT data */
  int * first_ft = cpl_table_get_data_int (vis_SC, "FIRST_FT");
  int * last_ft  = cpl_table_get_data_int (vis_SC, "LAST_FT");

  cpl_array ** visErr_ft   = cpl_table_get_data_array (vis_FT, "VISERR");
  cpl_array ** visData_ft  = cpl_table_get_data_array (vis_FT, "VISDATA");

  CPLCHECK_MSG ("Cannot get data");

  /* Create wavenumber to be faster */
  double meanwave_ft = cpl_table_get_column_mean (wave_table_ft, "EFF_WAVE");
  double * wavenumber_sc  = cpl_malloc (nwave_sc * sizeof(double));
  for (cpl_size wave = 0; wave < nwave_sc ; wave ++) {
	wavenumber_sc[wave] = 1. / cpl_table_get (wave_table_sc, "EFF_WAVE", wave, NULL);
  }

  CPLCHECK_MSG ("Cannot get data");
  
  /* New columns */
  gravi_table_new_column_array (vis_SC, "VISDATA_FT", "e", CPL_TYPE_DOUBLE_COMPLEX, nwave_ft);
  cpl_array ** visData_ftdit = cpl_table_get_data_array (vis_SC, "VISDATA_FT");

  gravi_table_new_column_array (vis_SC, "VISVAR_FT", "e^2", CPL_TYPE_DOUBLE, nwave_ft);
  cpl_array ** visVar_ftdit = cpl_table_get_data_array (vis_SC, "VISVAR_FT");

  gravi_table_new_column_array (vis_SC, "VISPOWER_FT", "e^2", CPL_TYPE_DOUBLE, nwave_ft);
  cpl_array ** visPower_ftdit = cpl_table_get_data_array (vis_SC, "VISPOWER_FT");

  gravi_table_new_column_array (vis_SC, "V_FACTOR", NULL, CPL_TYPE_DOUBLE, nwave_sc);
  cpl_array ** vFactor = cpl_table_get_data_array (vis_SC, "V_FACTOR");

  gravi_table_new_column_array (vis_SC, "V_FACTOR_FT", NULL, CPL_TYPE_DOUBLE, nwave_ft);
  cpl_array ** vFactor_ftdit = cpl_table_get_data_array (vis_SC, "V_FACTOR_FT");

  gravi_table_new_column (vis_SC, "V_FACTOR_WL", NULL, CPL_TYPE_DOUBLE);
  double * vFactor_wl = cpl_table_get_data_double (vis_SC, "V_FACTOR_WL");
  
  CPLCHECK_MSG ("Cannot create columns");

  /* Loop on base and row SC */
  for (cpl_size base = 0; base < nbase; base++) {
	for (cpl_size row_sc = 0; row_sc < nrow_sc; row_sc++) {
	  int nsc = row_sc * nbase + base;

	  /* Init integration of FT quantities during the SC frame */
	  visData_ftdit[nsc]  = gravi_array_init_double_complex (nwave_ft, 0.0 + I*0.0);
	  visVar_ftdit[nsc]   = gravi_array_init_double (nwave_ft, 0.0);
	  visPower_ftdit[nsc] = gravi_array_init_double (nwave_ft, 0.0);
	  
	  double vFactor_incoh = 0.0;
	  double complex vFactor_coh = 0.0 + I * 0.0;
	  double vFactor_var = 0.0;

	  /* Integrate quantities during the SC frame */
	  cpl_size nframe = last_ft[nsc] - first_ft[nsc];
	  for (cpl_size row_ft = first_ft[nsc] ; row_ft < last_ft[nsc]; row_ft++) {
		cpl_array * tmp_data;

		/* Integrate visVar_ftdit = < |visDataErr|^2 > over the current SC frame [e^2] */
		tmp_data = visErr_ft[row_ft * nbase + base];
		for (cpl_size wave = 0; wave < nwave_ft; wave ++){
		  cpl_array_set (visVar_ftdit[nsc], wave,
						 cpl_array_get (visVar_ftdit[nsc], wave, NULL) +
						 pow (cabs( cpl_array_get_complex (tmp_data, wave, NULL) ), 2));
		}

		/* Integrate visData_ftdit[nsc] = < visData > over the current SC frame [e] */
		tmp_data = visData_ft[row_ft * nbase + base];
		for (cpl_size wave = 0; wave < nwave_ft; wave ++){
			cpl_array_set_complex (visData_ftdit[nsc], wave,
								   cpl_array_get_complex (visData_ftdit[nsc], wave, NULL) +
								   cpl_array_get_complex (tmp_data, wave, NULL));
		}

		/* Integrate visPower_ftdit = < |visData|^2 > over the current SC frame [e^2] */
		for (cpl_size wave = 0; wave < nwave_ft; wave ++){
			cpl_array_set (visPower_ftdit[nsc], wave,
						   cpl_array_get (visPower_ftdit[nsc], wave, NULL) +
						   pow (cabs( cpl_array_get_complex (tmp_data, wave, NULL)), 2));
		}

		/* Integrate the same quantities for the white light vFactor 
		 * This is first doing a coherent integration over the wavelengths */
		vFactor_incoh += pow (cabs (cpl_array_get_mean_complex (visData_ft[row_ft * nbase + base])) * nwave_ft, 2);
		vFactor_coh += cpl_array_get_mean_complex (visData_ft[row_ft * 6 + base]) * nwave_ft;
		for (cpl_size wave = 0; wave < nwave_ft; wave ++) {
			vFactor_var += pow (cabs (cpl_array_get_complex (visErr_ft[row_ft * nbase + base], wave, NULL)), 2);
		}
		
		CPLCHECK_MSG("Issue in the loop to build average FT quantities");
	  } /* End loop on FT frame within this SC frame */

	  /* Compute the vFactor as the contrast attenuation within the SC DIT
	   *  (|<visData>|^2 - <|visErr|^2>) / (<|visData|^2> - <|visErr|^2>) / nframe */
	  vFactor_ftdit[nsc] = cpl_array_new (nwave_ft, CPL_TYPE_DOUBLE);
	  if (nframe != 0) {
		  for (cpl_size wave = 0; wave < nwave_ft; wave ++) {
			cpl_array_set (vFactor_ftdit[nsc], wave,
						   (pow (cabs (cpl_array_get_complex (visData_ftdit[nsc], wave, NULL)),2) -
							cpl_array_get (visVar_ftdit[nsc], wave, NULL)) /
						   (cpl_array_get (visPower_ftdit[nsc], wave, NULL) -
							cpl_array_get (visVar_ftdit[nsc], wave, NULL)) / (double)nframe);
		  }

		  /* Compute the white light vFactor, that is fist performing a coherent
		   * integration of the spectral channels. To lower the bias */
		  vFactor_wl[nsc] = (cabs(vFactor_coh)*cabs(vFactor_coh) - vFactor_var) /
							(vFactor_incoh - vFactor_var) / (double)nframe;
	  }
	  /* if nframe == 0  set the vFactor to 0 */
	  else {
		  for (cpl_size wave = 0; wave < nwave_ft; wave ++) {
			cpl_array_set (vFactor_ftdit[nsc], wave,0);
		  }
		  vFactor_wl[nsc] = 0;
	  }
	  
	  /* Compute the mean vFactor of the FT, and fit with a function exp(-a2/lbd2) to 
	   * project on the SC wavelength. This has only one free parameter */
	  vFactor[nsc] = cpl_array_new (nwave_sc, CPL_TYPE_DOUBLE);
	  double a2 = log (CPL_MAX(CPL_MIN (vFactor_wl[nsc], 1.0),1e-10)) * meanwave_ft*meanwave_ft;
	  for (cpl_size wave = 0; wave < nwave_sc; wave ++) {
		cpl_array_set (vFactor[nsc], wave, exp ( a2 * pow (wavenumber_sc[wave], 2)));
	  }
	  
	}
  }
  /* End loop on SC frames and base */

  FREE (cpl_free, wavenumber_sc);
  
  gravi_msg_function_exit(1);
  return CPL_ERROR_NONE;
}

/* -------------------------------------------------------------------------- */
/**
 * @brief Create the clockratio for each SC DIT and baseline.
 * 
 * @param vis_SC:          input/output OI_VIS table of the SC
 * @param vis_FT:          input OI_VIS table of the FT
 *
 * Compute the fraction of FT DIT not-rejected inside each SC DIT,
 * looking at the REJECTION_FLAG column in vis_FT.
 * The result is saved in a newly created column column FRINGEDET_RATIO
 * in the vis_SC column.
 * 
 * The synchronisation info shall already be computed.
 */
/* -------------------------------------------------------------------------- */

cpl_error_code gravi_vis_create_lockratio_sc (cpl_table * vis_SC,
											  cpl_table * vis_FT)
{
    gravi_msg_function_start(1);
    cpl_ensure_code (vis_SC, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code (vis_FT, CPL_ERROR_NULL_INPUT);
    
    /* Get variables */
    cpl_size nbase = 6;
    cpl_size nrow_sc  = cpl_table_get_nrow (vis_SC) / nbase;
    
    /* Get SC and FT data */
    int * reject_flag_ft = cpl_table_get_data_int (vis_FT, "REJECTION_FLAG");
    int * first_ft = cpl_table_get_data_int (vis_SC, "FIRST_FT");
    int * last_ft  = cpl_table_get_data_int (vis_SC, "LAST_FT");
    
    CPLCHECK_MSG ("Cannot get data");
    
    /* Create columns */
    gravi_table_new_column (vis_SC, "FRINGEDET_RATIO", NULL, CPL_TYPE_DOUBLE);
    double * fringedet_ftdit = cpl_table_get_data_double (vis_SC, "FRINGEDET_RATIO");
    
    CPLCHECK_MSG ("Cannot create columns");
    
    /* Loop on base and row SC */
    for (cpl_size base = 0; base < nbase; base++) {
        for (cpl_size row_sc = 0; row_sc < nrow_sc; row_sc++) {
            int nsc = row_sc * nbase + base;
            
            /* Integrate fraction of the SC frame with detected FT fringes */
            for (cpl_size row_ft = first_ft[nsc] ; row_ft < last_ft[nsc]; row_ft++) {
                fringedet_ftdit[nsc] += (reject_flag_ft[row_ft * nbase + base] == 0 ? 1 : 0);
            }
            
            /* Normalize the mean (if nframe == 0, value is zero) */
            cpl_size nframe = last_ft[nsc] - first_ft[nsc];
            if (nframe != 0) fringedet_ftdit[nsc] /= (double)nframe;
        }
    } /* End loop on base and SC frames */
    
    gravi_msg_function_exit(1);
    return CPL_ERROR_NONE;
}

/* -------------------------------------------------------------------------- */
/**
 * @brief Compute the reference phase for each SC DIT
 * 
 * @param vis_SC:   input/output OI_VIS table of the SC
 * @param wave_table_sc: wavelength table corresponding to OI_VIS
 * @param wave_table_ft: wavelength table corresponding to OI_VIS FT
 *
 * The reference phase is computed for each SC DIT and saved in a newly
 * created column PHASE_REF in the vis_SC table. It is constructed from
 * the coherent-flux of the FT, already averaged inside each SC DIT
 * (VISDATA_FT). It is then interpolatated to the SC wavelength table.
 *
 * The polynomial coefficients used to extrapolate the PHASE_REF from
 * the FT wavelengths to the SC wavelengths are saved in column
 * PHASE_REF_COEFF.
 */
/* -------------------------------------------------------------------------- */

cpl_error_code gravi_vis_create_phaseref_sc (cpl_table * vis_SC,
					     cpl_table * wavesc_table,
					     cpl_table * waveft_table,
					     cpl_propertylist * header,
					     const cpl_parameterlist * parlist)
{
  gravi_msg_function_start(1);
  cpl_ensure_code (vis_SC,        CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (wavesc_table, CPL_ERROR_NULL_INPUT);

  cpl_array ** visData_ftdit;
  double * gdelay_ftdit;
  const char * output_name = NULL, * coeff_name = NULL;
  if (waveft_table != NULL) {
      cpl_msg_info (cpl_func, "Compute reference phase of SC from the FT data");
      visData_ftdit = cpl_table_get_data_array (vis_SC, "VISDATA_FT");
      gdelay_ftdit  = cpl_table_get_data_double (vis_SC, "GDELAY_FT");
      output_name = "PHASE_REF";
      coeff_name  = "PHASE_REF_COEFF";
  } else {
      cpl_msg_info (cpl_func, "Compute reference phase of SC from the SC");
      visData_ftdit = cpl_table_get_data_array (vis_SC, "VISDATA");
      gdelay_ftdit  = cpl_table_get_data_double (vis_SC, "GDELAY");
      waveft_table  = wavesc_table;
      output_name = "SELF_REF";
      coeff_name  = "SELF_REF_COEFF";
  }
  
  /* Get general data */
  cpl_size nbase = 6;
  cpl_size nrow_sc  = cpl_table_get_nrow (vis_SC) / nbase;
  cpl_size nwave_sc = cpl_table_get_nrow (wavesc_table);
  cpl_size nwave_ft = cpl_table_get_nrow (waveft_table);
  
  CPLCHECK_MSG ("Cannot get data");

  /* Variable for fit */

  /* Create maxdeg as an option for user */
  cpl_size mindeg = 0;
  cpl_size maxdeg = gravi_param_get_int (parlist, "gravity.signal.phase-ref-sc-maxdeg");

  /* FE 2019-08-01: proper imaging phase requires higher order phase reference */
  // cpl_size mindeg = 0, maxdeg = 3;  
  /* SG 2019-08-07: trying order 2 and phase-calibration=FULL */
  /* cpl_size mindeg = 0, maxdeg = 2;  */

  cpl_msg_info (cpl_func, "phaseref with polynomial mindeg=%lli to maxdeg=%lli", mindeg, maxdeg);
  cpl_propertylist_update_int (header, "ESO QC PHASEREF_SC MINDEG", mindeg);
  cpl_propertylist_set_comment (header, "ESO QC PHASEREF_SC MINDEG", "fit of FT phase");
  cpl_propertylist_update_int (header, "ESO QC PHASEREF_SC MAXDEG", maxdeg);
  cpl_propertylist_set_comment (header, "ESO QC PHASEREF_SC MAXDEG", "fit of FT phase");

  
  cpl_polynomial * fit = cpl_polynomial_new (1);

  /* Create the vectors and matrix only once to be faster */
  cpl_matrix * sigma_ft = cpl_matrix_new (1,nwave_ft);
  cpl_vector * wave_sc = cpl_vector_new (nwave_sc);
  cpl_array * wavenumber_ft = cpl_array_new (nwave_ft, CPL_TYPE_DOUBLE);

  double lbd0   = cpl_table_get_column_mean (wavesc_table, "EFF_WAVE");
  double delta0 = cpl_table_get_column_max (wavesc_table, "EFF_WAVE") -
                  cpl_table_get_column_min (wavesc_table, "EFF_WAVE");
  for (cpl_size wave = 0; wave < nwave_ft; wave ++) {
      double lbd = cpl_table_get (waveft_table, "EFF_WAVE", wave, NULL);
      cpl_matrix_set (sigma_ft, 0, wave, (lbd0/lbd - 1.) * lbd0/delta0 );
      cpl_array_set (wavenumber_ft, wave, 1./lbd);
  }
  for (cpl_size wave = 0; wave < nwave_sc; wave ++) {
	cpl_vector_set (wave_sc, wave, cpl_table_get (wavesc_table, "EFF_WAVE", wave, NULL));
  }

  CPLCHECK_MSG ("Cannot create wave arrays");

  /* Create columns */
  gravi_table_new_column_array (vis_SC, output_name, "rad", CPL_TYPE_DOUBLE, nwave_sc);
  cpl_array ** phaseref = cpl_table_get_data_array (vis_SC, output_name);

  /* Create columns */
  gravi_table_new_column_array (vis_SC, coeff_name, NULL, CPL_TYPE_DOUBLE, maxdeg+1);
  cpl_array ** phase_coeff = cpl_table_get_data_array (vis_SC, coeff_name);
  
  CPLCHECK_MSG ("Cannot create column");
		
  /* Loop on base and SC frames */
  for (cpl_size base = 0; base < nbase; base++) {
	for (cpl_size row_sc = 0; row_sc < nrow_sc; row_sc ++) {
	  int nsc = row_sc * nbase + base;
			  
	  /* phaseref_ftdit is the FT phase arg{visData_ftdit}
	   * Need to be unwrapped before interpolation. */
	  cpl_array * phaseref_ftdit = cpl_array_cast (visData_ftdit[nsc], CPL_TYPE_DOUBLE_COMPLEX);
			  
	  /* Remove mean group-delay and phase-delay to unwrap */
	  gravi_array_multiply_phasor (phaseref_ftdit, - 2*I*CPL_MATH_PI * gdelay_ftdit[nsc], wavenumber_ft);
	  double mean_phase = carg (cpl_array_get_mean_complex (phaseref_ftdit));
	  cpl_array_multiply_scalar_complex (phaseref_ftdit, cexp(- I * mean_phase));
			  
	  /* Compute argument and add back the delay and the phase [rad] */
	  cpl_array_arg (phaseref_ftdit);
	  gravi_array_add_phase (phaseref_ftdit, 2.*CPL_MATH_PI*gdelay_ftdit[nsc], wavenumber_ft);
	  cpl_array_add_scalar (phaseref_ftdit, mean_phase);

	  /* Interpolate the FT phase at the SC wavelengths with a polynomial of order 2
	   * rewrap the phase, and make it phase_ref = - phase_ft */

      /* Polynomial fit */
      cpl_vector * input = cpl_vector_wrap (nwave_ft, cpl_array_get_data_double (phaseref_ftdit));
      cpl_polynomial_fit (fit, sigma_ft, NULL, input, NULL, CPL_FALSE, &mindeg, &maxdeg);
      cpl_vector_unwrap (input);
	  cpl_array_delete (phaseref_ftdit);

      /* Save fit coefficients */
      phase_coeff[nsc] = cpl_array_new (maxdeg+1, CPL_TYPE_DOUBLE);
      for (cpl_size d = 0; d < maxdeg+1; d++) 
          cpl_array_set (phase_coeff[nsc], d, cpl_polynomial_get_coeff (fit, &d));
      
      /* Evaluate polynomial at the output sampling */
      phaseref[nsc] = cpl_array_new (nwave_sc, CPL_TYPE_DOUBLE);
      for (cpl_size w = 0; w < nwave_sc; w++) {
          double delta = (lbd0/cpl_vector_get(wave_sc, w) - 1.) * lbd0/delta0 ;
          cpl_array_set (phaseref[nsc], w, cpl_polynomial_eval_1d (fit, delta, NULL));
      }

	  gravi_array_phase_wrap (phaseref[nsc]);
	  cpl_array_multiply_scalar (phaseref[nsc], -1.0);
			  
	  CPLCHECK_MSG ("Cannot compute the PHASE_REF for SC");
	} /* End loop on SC frames */
  } /* End loop on base */

  FREE (cpl_vector_delete, wave_sc);
  FREE (cpl_matrix_delete, sigma_ft);
  FREE (cpl_array_delete, wavenumber_ft);
  FREE (cpl_polynomial_delete, fit);

  gravi_msg_function_exit(1);
  return CPL_ERROR_NONE;
}

/* -------------------------------------------------------------------------- */
/**
 * @brief Create the list of outlier based on the values of the chi2.
 * 
 * @param flux_SC:          input/output OI_FLUX table of the SC
 * @param vis_SC:           input/output OI_VIS table of the SC
 * @param chi2r_threshold:  threshold for absolute value of chi2r
 * @param chi2r_sigma:      threshold in number of sigma
 * @param stat[2]:          stat[0] is number of channels with at least 
 *                          one outlier. stat[1] is number of channels
 *                          with more then 50% outliers.
 *
 */
/* -------------------------------------------------------------------------- */

cpl_error_code gravi_create_outlier_flag_sc (cpl_table * flux_SC,
					     cpl_table * vis_SC,
					     double chi2r_threshold,
					     double chi2r_sigma)
{
  gravi_msg_function_start(1);
  cpl_ensure_code (flux_SC, CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (vis_SC,  CPL_ERROR_NULL_INPUT);

  cpl_msg_info (cpl_func, "chi2_threshold = %f", chi2r_threshold);
  cpl_msg_info (cpl_func, "chi2_sigma = %f", chi2r_sigma);

  cpl_size ntel = 4;
  cpl_size nrow = cpl_table_get_nrow (flux_SC) / ntel;
  cpl_size nwave = cpl_table_get_column_depth (flux_SC, "FLUX");

  cpl_array ** chi2 = cpl_table_get_data_array (flux_SC, "CHI2");

  /* Create column */
  cpl_array ** flux_outliers = cpl_table_get_data_array (flux_SC, "FLAG");

  CPLCHECK_MSG ("Cannot create column");
  
  /* vector to store data */
  cpl_vector * vector = cpl_vector_new (nwave);
  cpl_vector * med = NULL;
  
  /* Loop on rows */
  for (cpl_size row = 0; row < nrow; row ++) {

    /* Chi2 data as vector so that we can apply filtering */
    for (cpl_size wave = 0; wave < nwave ; wave ++) {
      double value = cpl_array_get (chi2[row*ntel+0], wave, NULL);
      cpl_vector_set (vector, wave, value);
    }

    /* Normalise the chi2 by median in spectral direction */
    med = gravi_vector_median (vector, 50);
    /* remove risk of divide by zero */
    cpl_vector_add_scalar(med, cpl_vector_get_mean(med)*1e-9);
    cpl_vector_divide (vector, med);
    FREE (cpl_vector_delete, med);

    /* Threshold on normalised value */
    for (cpl_size wave = 0; wave < nwave ; wave ++) {
      double value = cpl_vector_get (vector, wave);
      if (value > chi2r_threshold) cpl_array_set (flux_outliers[row*ntel+0], wave, 1);
    }

    /* Compute the distance to median,
       in unit of variance */
    cpl_vector_subtract_scalar (vector, 1.0);
    gravi_vector_abs (vector);
    med = gravi_vector_median (vector, 50);
    /* remove risk of divide by zero */
    cpl_vector_add_scalar(med, cpl_vector_get_mean(med)*1e-9);
    cpl_vector_divide (vector, med);
    FREE (cpl_vector_delete, med);
    
    /* Flag on distance in unit of local variance */
    for (cpl_size wave = 0; wave < nwave ; wave ++) {
      double value = cpl_vector_get (vector, wave);
      if (value > chi2r_sigma) cpl_array_set (flux_outliers[row*ntel+0], wave, 1);
    }
  }
              
  /* Free memory */
  FREE (cpl_vector_delete, vector);

  CPLCHECK_MSG ("Cannot fill outliers");

  /* Duplicate the detection of outliers in flux_SC,
     based on the one of the first beam */
  
  /* Loop on row */
  for (cpl_size row = 0; row < nrow; row ++) {
    for (cpl_size tel = 1; tel < ntel; tel++) {
      cpl_array_add (flux_outliers[row*ntel+tel], flux_outliers[row*ntel+0]);
    }
  }    

  CPLCHECK_MSG ("Cannot duplicate in OI_FLUX");

  /* Duplicate the detection of outliers in vis_SC,
     based on the one of the first beam */
  
  /* Create column */
  cpl_array ** vis_outliers = cpl_table_get_data_array (vis_SC, "FLAG");

  /* Loop on row */
  int nbase = 6;
  for (cpl_size row = 0; row < nrow; row ++) {
    for (cpl_size base = 0; base < nbase; base++) {
      cpl_array_add (vis_outliers[row*nbase+base], flux_outliers[row*ntel+0]);
    }
  }

  CPLCHECK_MSG ("Cannot duplicate in OI_VIS");

  gravi_msg_function_exit(1);
  return CPL_ERROR_NONE;
}

/* -------------------------------------------------------------------------- */
/**
 * @brief Create the list of outlier. For the FT, this is filled with 0
 * 
 * @param flux_FT:          input/output OI_FLUX table of the FT
 * @param vis_FT:           input/output OI_VIS table of the FT
 *
 */
/* -------------------------------------------------------------------------- */

cpl_error_code gravi_create_outlier_flag_ft (cpl_table * flux_FT,
					     cpl_table * vis_FT)
{
  gravi_msg_function_start(1);
  cpl_ensure_code (flux_FT, CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (vis_FT,  CPL_ERROR_NULL_INPUT);

  /* Nothing is checked for FT */

  gravi_msg_function_exit(1);
  return CPL_ERROR_NONE;
}

/* -------------------------------------------------------------------------- */
/**
 * @brief Compute the (FDDL_SC + FDDL_FT)/2 position in [m]
 * 
 * @param flux_SC:     input/output OI_FLUX table of the SC
 * @param disp_table:  FDDL dispersion model
 *
 * Create new columns in flux_SC (FDDL)
 */
/* -------------------------------------------------------------------------- */

cpl_error_code gravi_flux_create_fddllin_sc (cpl_table * flux_SC,
                                             cpl_table * disp_table)
{
  gravi_msg_function_start(1);
  cpl_ensure_code (flux_SC,    CPL_ERROR_NULL_INPUT);
  
  cpl_size ntel = 4;
  cpl_size nrow = cpl_table_get_nrow (flux_SC) / ntel;

  /* Create the columns */
  gravi_table_new_column (flux_SC, "FDDL", "m", CPL_TYPE_DOUBLE);
  CPLCHECK_MSG ("Cannot create columns");
    gravi_table_new_column (flux_SC, "FDDL_FT", "m", CPL_TYPE_DOUBLE);
    CPLCHECK_MSG ("Cannot create columns");
    gravi_table_new_column (flux_SC, "FDDL_SC", "m", CPL_TYPE_DOUBLE);
    CPLCHECK_MSG ("Cannot create columns");
  
  /* If not DISP_DATA, we just create the columns */
  if (disp_table == NULL) {
	gravi_msg_function_exit(1);
	return CPL_ERROR_NONE;
  }

  /* Get the list of coeficients for FDDL linearity */
  cpl_msg_info (cpl_func, "Load the linearity model coeficients");
		  
  double ** lin_fddl_sc = gravi_table_get_data_array_double (disp_table, "LIN_FDDL_SC");
  double ** lin_fddl_ft = gravi_table_get_data_array_double (disp_table, "LIN_FDDL_FT");
  cpl_size disp_order   = cpl_table_get_column_depth (disp_table, "LIN_FDDL_SC");
  CPLCHECK_MSG ("Cannot get linearity model data");
  
  /* Get data */
  double * ftpos = cpl_table_get_data_double (flux_SC, "FT_POS");
  double * scpos = cpl_table_get_data_double (flux_SC, "SC_POS");
    double * fddl  = cpl_table_get_data_double (flux_SC, "FDDL");
    double * fddl_ft  = cpl_table_get_data_double (flux_SC, "FDDL_FT");
    double * fddl_sc  = cpl_table_get_data_double (flux_SC, "FDDL_SC");
  CPLCHECK_MSG ("Cannot get POS data");
  
  /* Loop on tel and frames */
  for (cpl_size tel = 0; tel < ntel; tel++) {
      for (cpl_size row = 0; row < nrow; row ++) {
          cpl_size nsc = row * ntel + tel;

          /* Apply the non-linearity to the FDDL while 
           * computing the mean of SC and FT: [V] -> [m] */
          fddl_ft[nsc] = 0.0;
          fddl_sc[nsc] = 0.0;
          for (int o = 0; o < disp_order; o++) {
              fddl_sc[nsc] += lin_fddl_sc[tel][o] * pow (scpos[nsc], (double)o) * 1.0e-6;
              fddl_ft[nsc] += lin_fddl_ft[tel][o] * pow (ftpos[nsc], (double)o) * 1.0e-6;
          }
          fddl[nsc]=0.5*(fddl_ft[nsc]+fddl_sc[nsc]);
      } /* End loop on rows */
  } /* End loop on base */

  /* Free the temporary allocations */
  FREE (cpl_free, lin_fddl_sc);
  FREE (cpl_free, lin_fddl_ft);
		  
  gravi_msg_function_exit(1);
  return CPL_ERROR_NONE;
}

/* -------------------------------------------------------------------------- */
/**
 * @brief Compute the MET opd including the dispersion from the FDDL
 * 
 * @param vis_SC:      input/output OI_VIS table of the SC
 * @param flux_SC:     input/output OI_FLUX table of the SC
 * @param wave_table:  wavelength table corresponding to OI_VIS
 * @param disp_table:  FDDL dispersion model
 * @param parlist:     parameter list of the recipe
 *
 * Create new columns in vis_SC (OPD per base) and flux_SC (OPD per beam)
 * by combining the already computed MET and FDDL signals averaged inside
 * each SC DIT and the dispersion coefficient from disp_table.
 */
/* -------------------------------------------------------------------------- */

cpl_error_code gravi_vis_create_opddisp_sc (cpl_table * vis_SC,
					    cpl_table * flux_SC,
					    cpl_table * wave_table,
					    cpl_table * disp_table,
					    cpl_propertylist * header,
					    const cpl_parameterlist * parlist)
{
  gravi_msg_function_start(1);
  cpl_ensure_code (vis_SC,     CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (wave_table, CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (header,     CPL_ERROR_NULL_INPUT);  
  
  cpl_size nwave_sc = cpl_table_get_column_depth (vis_SC, "VISDATA");

  /* If not DISP_DATA, we cannot create these columns */
  if (disp_table == NULL) {
      cpl_msg_info (cpl_func,"Cannot create OPD_DISP, not DISP_MODEL table");
      gravi_msg_function_exit(1);
      return CPL_ERROR_NONE;
  }

		  
  /* Compute the N_MEAN(lbd) and N_DIFF(lbd) for each beam from their
   * polynomial model: n = Sum[ coef_i x ((lbd-lbd_met)/lbd_met)^i ]  */
  cpl_msg_info (cpl_func, "Load the dispersion model coeficients");

  /* Allocate memory */
  double ** n_mean = cpl_malloc (4 * sizeof(double*));
  double ** n_diff = cpl_malloc (4 * sizeof(double*));
  for (int t = 0; t < 4; t++) {
      n_mean[t] = cpl_calloc (nwave_sc, sizeof(double));
      n_diff[t] = cpl_calloc (nwave_sc, sizeof(double));
  }

  if ( cpl_table_has_column (disp_table,"N_MEAN") &&
       cpl_table_has_column (disp_table,"WAVE0")) {
      /* The N_MEAN and N_DIFF are described as n(lbd) / n(lbdmet)
       * with a polynomial versus (wave0/wave-1) */
      cpl_msg_info (cpl_func, "Dispersion model as N/Nmet  [wave0/wave-1]");
      cpl_size disp_order = cpl_table_get_column_depth (disp_table, "N_MEAN");
      for (int t = 0; t < 4; t++) {
          double lbd0 = cpl_table_get (disp_table, "WAVE0", t, NULL);
          for (int w = 0; w < nwave_sc; w++ ) {
              double lbd = cpl_table_get (wave_table, "EFF_WAVE", w, NULL);
              double xfit = lbd0/lbd - 1.0;
              for (int o = 0; o < disp_order; o++) {
                  n_mean[t][w] += gravi_table_get_value (disp_table, "N_MEAN", t, o) * pow (xfit, o);
                  n_diff[t][w] += gravi_table_get_value (disp_table, "N_DIFF", t, o) * pow (xfit, o);
              }
          }
      }
  } else {
      FREELOOP (cpl_free, n_mean, 4);
      FREELOOP (cpl_free, n_diff, 4);
      cpl_msg_error (cpl_func,"The DISP_MODEL is not recognized... contact the DRS team");
      return cpl_error_set_message (cpl_func,CPL_ERROR_ILLEGAL_INPUT,
                                    "The DISP_MODEL is not recognized... contact the DRS team");
  }
  
  CPLCHECK_MSG ("Cannot compute N_MEAN or N_DIFF");

  cpl_size nbase = 6, ntel = 4;
  cpl_size nrow_sc = cpl_table_get_nrow (vis_SC) / nbase;

  /* Create the columns */
  gravi_table_new_column_array (vis_SC, "OPD_DISP", "m", CPL_TYPE_DOUBLE, nwave_sc);
  cpl_array ** opd_disp = cpl_table_get_data_array (vis_SC, "OPD_DISP");

  gravi_table_new_column (vis_SC, "GDELAY_DISP", "m", CPL_TYPE_DOUBLE);
  double * gd_disp = cpl_table_get_data_double (vis_SC, "GDELAY_DISP");

  gravi_table_new_column_array (vis_SC, "PHASE_DISP", "rad", CPL_TYPE_DOUBLE, nwave_sc);
  cpl_array ** phase_disp = cpl_table_get_data_array (vis_SC, "PHASE_DISP");
  
  CPLCHECK_MSG ("Cannot create columns");

  /* Get data */
  double * opd_met = cpl_table_get_data_double (flux_SC, "OPD_MET_FC");
  double * fddl    = cpl_table_get_data_double (flux_SC, "FDDL");

  CPLCHECK_MSG ("Cannot get data");

  /* Create wavenumber to be faster */
  double * wavenumber_sc  = cpl_malloc (nwave_sc * sizeof(double));
  for (cpl_size wave = 0; wave < nwave_sc ; wave ++) {
	wavenumber_sc[wave] = 1. / cpl_table_get (wave_table, "EFF_WAVE", wave, NULL);
  }

  /* Determine OPD_MET_ZERO_FC, based on the content of the header */
  int t;
  char name[100];
  double * opl_zero_fc = cpl_malloc (4 * sizeof(double));
  // double gd_zero_fc;
  
  /* Initialise opl_zero_fc to zero */
  for (t = 0; t < ntel; t++) {
    opl_zero_fc[t] = 0.0;
  }
  if ( gravi_param_get_bool (parlist, "gravity.signal.use-met-zero") ) {
    cpl_msg_info (cpl_func, "Metrology zero calculation is enabled!");
    
    /* Replace by OCS MET OPL_ZERO_FC, if available */
    for (t = 0; t < ntel; t++) {
      sprintf (name, "ESO OCS MET OPL_ZERO_FC%i", t+1); 
      if (cpl_propertylist_has (header, name)) {
        opl_zero_fc[t] = cpl_propertylist_get_double (header, name)*1e-3;
        sprintf (name, "ESO FDDL MET OFFSET%i", t+1);
        opl_zero_fc[t] += cpl_propertylist_get_double (header, name)*1e-3;
        cpl_msg_info (cpl_func, "Updating metrology zero with OCS MET OPL_ZERO_FC%i and FDDL MET OFFSET%i: %f [mm]", t+1, t+1, opl_zero_fc[t]);
      }
    }
    
    /* Replace by PRO MET GD_ZERO_FC, if available */
    for (t = 0; t < ntel; t++) {
      if ( gravi_pfits_has_gdzero (header, t+1) ) {
          opl_zero_fc[t] =
                  gravi_pfits_get_gdzero (header, t+1)*1e-3
                  * (wavenumber_sc[nwave_sc/2-1] - wavenumber_sc[nwave_sc/2+1])
                  / (n_mean[t][nwave_sc/2-1] * wavenumber_sc[nwave_sc/2-1] - n_mean[t][nwave_sc/2+1] * wavenumber_sc[nwave_sc/2+1]);
          cpl_msg_info (cpl_func, "Updating metrology zero with QC/PRO MET GD_ZERO_FC%i: %f [mm]", t+1, opl_zero_fc[t]);
      }
    }
    // Replace by PRO MET OPL_ZERO_FC, if available
    for (t = 0; t < ntel; t++) {
      if ( gravi_pfits_get_oplzero (header, t+1) ){
          opl_zero_fc[t] = gravi_pfits_get_oplzero (header, t+1)*1e-3;
          cpl_msg_info (cpl_func, "Updating metrology zero with QC/PRO MET OPL_ZERO_FC%i: %f [mm]", t+1, opl_zero_fc[t]);
      }
    }
  } else {
    cpl_msg_info (cpl_func, "Metrology zero calculation is disabled!");
  }
    
    
  /* Loop on tel and frames */
  for (cpl_size base = 0; base < nbase; base++) {
	for (cpl_size row_sc = 0; row_sc < nrow_sc; row_sc ++) {
        
	  cpl_size nsc = row_sc * nbase + base;
	  cpl_size t0 = GRAVI_BASE_TEL[base][0], t0f = row_sc * ntel + t0;
	  cpl_size t1 = GRAVI_BASE_TEL[base][1], t1f = row_sc * ntel + t1;

	  /* Compute the group-delay introduced by FDDL in the middle of the band
	   * This is obviously to a constant */
	  cpl_size wave = nwave_sc / 2;
	  double s1 = wavenumber_sc[wave-1];
	  double s2 = wavenumber_sc[wave+1];
	  double o1 =
	    (n_mean[t0][wave-1] * (opd_met[t0f]-opl_zero_fc[t0]) + n_diff[t0][wave-1] * fddl[t0f]) -
	    (n_mean[t1][wave-1] * (opd_met[t1f]-opl_zero_fc[t1]) + n_diff[t1][wave-1] * fddl[t1f]);
	  double o2 =
	    (n_mean[t0][wave+1] * (opd_met[t0f]-opl_zero_fc[t0]) + n_diff[t0][wave+1] * fddl[t0f]) -
	    (n_mean[t1][wave+1] * (opd_met[t1f]-opl_zero_fc[t1]) + n_diff[t1][wave+1] * fddl[t1f]);
				  
	  gd_disp[nsc] = (o1*s1 - o2*s2) / (s1-s2);

	  /* Compute the phase introduced by FDDL. This is also to a dispersive constant 
	   * since the amount of fiber and vaccum is not known */
	  phase_disp[nsc] = cpl_array_new (nwave_sc, CPL_TYPE_DOUBLE);
				
	  for (int w = 0; w < nwave_sc ; w++) {
		cpl_array_set (phase_disp[nsc], w,
			       ((n_mean[t0][w] * (opd_met[t0f]-opl_zero_fc[t0]) + n_diff[t0][w] * fddl[t0f]) -
				(n_mean[t1][w] * (opd_met[t1f]-opl_zero_fc[t1]) + n_diff[t1][w] * fddl[t1f]) -
						gd_disp[nsc]) * wavenumber_sc[w] * CPL_MATH_2PI);
	  }

	  gravi_array_phase_wrap (phase_disp[nsc]);

	  /* Compute the OPD_DISP for each wavelength :
	   * T0 = N_DIFF0 * (FDDL_SC0+FDDL_FT0)/2 + N_MEAN0 * (MET_SC0-MET_FT0)
	   * T1 = N_DIFF1 * (FDDL_SC1+FDDL_FT1)/2 + N_MEAN1 * (MET_SC1-MET_FT1) */
	  opd_disp[nsc] = cpl_array_new (nwave_sc, CPL_TYPE_DOUBLE);
				
	  for (int w = 0; w < nwave_sc ; w++) {
		cpl_array_set (opd_disp[nsc], w,
			       (n_mean[t0][w] * (opd_met[t0f]-opl_zero_fc[t0]) + n_diff[t0][w] * fddl[t0f]) -
			       (n_mean[t1][w] * (opd_met[t1f]-opl_zero_fc[t1]) + n_diff[t1][w] * fddl[t1f]) );
	  }

	} /* End loop on rows */
  } /* End loop on base */

  /* Free the temporary allocations */
  FREELOOP (cpl_free, n_mean, 4);
  FREELOOP (cpl_free, n_diff, 4);
  FREE (cpl_free, wavenumber_sc);
  FREE (cpl_free, opl_zero_fc);
		  
  gravi_msg_function_exit(1);
  return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Create phase-referenced imaging data in the P2VMREDUCED file
 *
 * @param vis_SC:      input/output OI_VIS table of the SC
 * @param wave_table:  wavelength table corresponding to the OI_VIS above
 * @param header:      main header
 *
 * Create IMAGING_REF column in vis_SC based on the following:
 * IMAGING_REF = PHASE_REF - OPD_DISP * (2pi/EFF_WAVE)
 *         + (UCOORD * SOBJ_X + VCOORD * SOBJ_Y) * (2pi/EFF_WAVE)
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_vis_create_imagingref_sc (cpl_table * vis_SC,
                                                cpl_table * wave_table,
                                                cpl_propertylist * header,
                                                const cpl_parameterlist * parlist)
{
    gravi_msg_function_start(1);
    cpl_ensure_code (vis_SC,     CPL_ERROR_NULL_INPUT);
    cpl_ensure_code (wave_table, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code (header,     CPL_ERROR_NULL_INPUT);

    /* Get from OI_VIS the number of wavelength channels and rows */
    cpl_size nwave = cpl_table_get_column_depth (vis_SC, "VISDATA");
    cpl_size nrow  = cpl_table_get_nrow (vis_SC);

    /* Get from OI_VIS the data needed for the calculation */
    double * ucoord        = cpl_table_get_data_double (vis_SC, "UCOORD");
    double * vcoord        = cpl_table_get_data_double (vis_SC, "VCOORD");
    cpl_array ** phase_ref = cpl_table_get_data_array  (vis_SC, "PHASE_REF");
    double * opd_met_telfc_mcorr = NULL;
    double * opd_met_fc_corr     = NULL;

    CPLCHECK_MSG ("Cannot get input data");

    /* If OPD_DISP is not computed, we cannot compute IMAGING_REF neither */
    if ( !cpl_table_has_column (vis_SC, "OPD_DISP") ) {
        cpl_msg_info (cpl_func,"Cannot compute IMAGING_REF, not column OPD_DISP");
        gravi_msg_function_exit(1);
        return CPL_ERROR_NONE;
    }

    cpl_array ** opd_disp  = cpl_table_get_data_array  (vis_SC, "OPD_DISP");
    CPLCHECK_MSG ("Cannot get OPD_DISP data");
    
    cpl_array ** phase_met_telfc  = cpl_table_get_data_array  (vis_SC, "PHASE_MET_TELFC");
    CPLCHECK_MSG ("Cannot get PHASE_MET_TELFC data");

    /* Get from header the separation, converted from mas to radian */
    double sep_U = gravi_pfits_get_sobj_x (header)*1e-3/3600.0/CPL_MATH_DEG_RAD;
    double sep_V = gravi_pfits_get_sobj_y (header)*1e-3/3600.0/CPL_MATH_DEG_RAD;

    /* Create new PHASE_REF_IMG column array */
    gravi_table_new_column_array (vis_SC, "IMAGING_REF", "rad", CPL_TYPE_DOUBLE, nwave);
    cpl_array ** imaging_ref = cpl_table_get_data_array (vis_SC, "IMAGING_REF");
    CPLCHECK_MSG ("Cannot create column");

    /* Which megtrology should be used for IMAGING_REF calculation? */
    const char * imaging_ref_met = gravi_param_get_string (parlist, "gravity.signal.imaging-ref-met");
    cpl_msg_info (cpl_func, "imaging-ref-met = %s", imaging_ref_met);
    if (!strcmp (imaging_ref_met,"TEL")) {
        cpl_msg_info (cpl_func,"Use telescope metrology for IMAGING_REF computation");
        opd_met_telfc_mcorr = cpl_table_get_data_double (vis_SC, "OPD_MET_TELFC_MCORR");
        opd_met_fc_corr = cpl_table_get_data_double (vis_SC, "OPD_MET_FC_CORR");
    } else if (!strcmp (imaging_ref_met,"FC_CORR")) {
        cpl_msg_info (cpl_func,"Use corrected fiber coupler metrology for IMAGING_REF computation");
        opd_met_fc_corr = cpl_table_get_data_double (vis_SC, "OPD_MET_FC_CORR");
    } else if (!strcmp (imaging_ref_met,"FC")) {
        cpl_msg_info (cpl_func,"Use fiber coupler metrology for IMAGING_REF computation");
    } else {
        cpl_msg_error (cpl_func,"Unknown metrology source for IMAGING_REF calculation!");
        cpl_ensure_code (imaging_ref_met, CPL_ERROR_ILLEGAL_INPUT);
    }

    /* Compute the reference phase for each row */
    for (cpl_size row = 0; row < nrow; row ++) {

        /* New array */
        imaging_ref[row] = cpl_array_new (nwave, CPL_TYPE_DOUBLE);

        /* If requested, use fiber coupler correction or telescope metrology */
	// TODO : FE 20190509 added opd_met_ttpup (commented)?
        double opd_met_corr = (opd_met_telfc_mcorr ? opd_met_telfc_mcorr[row] : 0.0)
	                          + (opd_met_fc_corr ? opd_met_fc_corr[row] : 0.0);
    //                            + (opd_met_ttpup ? opd_met_ttpup[row] : 0.0);

        /* Compute VISPHI for each wavelength */
        for (int w = 0; w < nwave; w++) {
            
            double wavelength = cpl_table_get (wave_table, "EFF_WAVE", w, NULL);
            
            /* IMAGING_REF = PHASE_REF - OPD_DISP * (2PI/EFF_WAVE) + 
               (UCOORD*SOBJ_X + VCOORD*SOBJ_Y) * (2PI/EFF_WAVE) */
            cpl_array_set (imaging_ref[row], w,
                           cpl_array_get (phase_ref[row], w, NULL)
                           - cpl_array_get (opd_disp[row],  w, NULL)  * CPL_MATH_2PI / wavelength
                           + (ucoord[row] * sep_U + vcoord[row] * sep_V) * CPL_MATH_2PI / wavelength
                           - opd_met_corr * CPL_MATH_2PI / wavelength);
            
            /* Here, we are using the new computation of phase_met_telfc (overide previous calculation) */
            if (!strcmp (imaging_ref_met,"TEL"))
                cpl_array_set (imaging_ref[row], w,
                               cpl_array_get (phase_ref[row], w, NULL)
                               - cpl_array_get (opd_disp[row],  w, NULL)  * CPL_MATH_2PI / wavelength
                               + (ucoord[row] * sep_U + vcoord[row] * sep_V) * CPL_MATH_2PI / wavelength
                               - cpl_array_get (phase_met_telfc[row],  w, NULL) );
                
            CPLCHECK_MSG ("Cannot compute the imaging phase");
        }

        /* Wrap this phase in [rad] */
        gravi_array_phase_wrap (imaging_ref[row]);
    }

    gravi_msg_function_exit(1);
    return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Create intermediate signal in the P2VMREDUCED file
 * 
 * @param p2vmred_data:  the P2VMREDUCED data (modified in-place)
 * @param disp_data:     DISP_PARAM file, to compute the OPD_DISP
 * @param parlist:       parameter list of the recipe (unused so far)
 *
 * Create many intermediate signal in the P2VMREDUCED file, necesary
 * for the further averaging of frames. Especially phase references,
 * vFactor...
 * These computations are mandatory for the further averaging of the
 * frames into an OIFITS file with gravi_compute_vis
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_compute_signals (gravi_data * p2vmred_data,
                                      gravi_data * disp_data,
                                      const cpl_parameterlist * parlist)
{
  gravi_msg_function_start(1);
  cpl_ensure_code (p2vmred_data, CPL_ERROR_NULL_INPUT);

  int nbase = 6, ntel = 4;

  /* Get header data */
  cpl_propertylist * p2vmred_header = gravi_data_get_header (p2vmred_data);
  int npol_ft = gravi_pfits_get_pola_num (p2vmred_header, GRAVI_FT);
  int npol_sc = gravi_pfits_get_pola_num (p2vmred_header, GRAVI_SC);
  
  CPLCHECK_MSG ("Cannot get header");

  /* FIXME: probably doesn't work with npol_sc != npol_ft */
  if ( npol_sc != npol_ft ) {
      gravi_msg_warning ("FIXME", "Not sure this function works with npol_sc != npol_ft");
  }

  
  /* Get the SC DIT */
  double dit_sc = gravi_pfits_get_dit_sc (p2vmred_header) * 1e6;
  cpl_msg_info (cpl_func,"dit_sc = %g [us]", dit_sc);
  

  /*
   * Loop on polarisations 
   */
  for (int pol = 0; pol < CPL_MAX(npol_sc,npol_ft); pol++) {
	
	/* verbose */
	cpl_msg_info (cpl_func, "Start polarisation %d over %d",pol+1, CPL_MAX(npol_sc,npol_ft));
	cpl_msg_info(cpl_func, "Insname FT : %s, pol %d npol %d",
                 GRAVI_INSNAME(GRAVI_FT,pol,npol_ft), pol+1, npol_ft);
	cpl_msg_info(cpl_func, "Insname SC : %s, pol %d npol %d",
                 GRAVI_INSNAME(GRAVI_SC,pol,npol_sc), pol+1, npol_sc);

	/* Get the table of reduced data from FT */
	cpl_table * vis_FT = gravi_data_get_oi_vis (p2vmred_data, GRAVI_FT, pol, npol_ft);
	cpl_table * flux_FT = gravi_data_get_oi_flux (p2vmred_data, GRAVI_FT, pol, npol_ft);
	CPLCHECK_MSG ("Cannot get the FT tables");
	
	/* Get the table of reduced data from SC */
	cpl_table * vis_SC = gravi_data_get_oi_vis (p2vmred_data, GRAVI_SC, pol, npol_sc);
	cpl_table * flux_SC = gravi_data_get_oi_flux (p2vmred_data, GRAVI_SC, pol, npol_sc);
	CPLCHECK_MSG ("Cannot get the SC tables");
	
	/* Get the metrology and FDDL */
	cpl_table * vis_met = gravi_data_get_table (p2vmred_data, GRAVI_OI_VIS_MET_EXT);
	cpl_table * fddl_table = gravi_data_get_table (p2vmred_data, GRAVI_FDDL_EXT);
	CPLCHECK_MSG ("Cannot get the VIS_MET and FDDL tables");
	
	/* Get the OIFITS tables that are alredy in the data */
	cpl_table * oi_wavelengthft = gravi_data_get_oi_wave (p2vmred_data, GRAVI_FT, pol, npol_ft);
	cpl_table * oi_wavelengthsc = gravi_data_get_oi_wave (p2vmred_data, GRAVI_SC, pol, npol_sc);
	CPLCHECK_MSG ("Cannot get the OI_WAVELENGTH tables");


	/* 
	 * (1) Create synchronisation information
	 */
	
	/* Create the FIRST_FT and LAST_FT for each SC frame */
	gravi_signal_create_sync (vis_SC, 6, dit_sc, vis_FT, 6, "FT");
	gravi_signal_create_sync (vis_SC, 6, dit_sc, vis_met, 4, "MET");
	
	CPLCHECK_MSG ("Cannot sync vis_SC");
	
	/* Create the FIRST_MET and LAST_MET for each SC frame */
	gravi_signal_create_sync (flux_SC, 4, dit_sc, flux_FT, 4, "FT");
	gravi_signal_create_sync (flux_SC, 4, dit_sc, vis_met, 4, "MET");
	gravi_signal_create_sync (flux_SC, 4, dit_sc, fddl_table, 1, "FDDL");
	
	CPLCHECK_MSG ("Cannot sync flux_SC");
	
	/* 
	 * (2) Create the signals for FLUX_FT
	 */
	
	cpl_array ** flux_ft = cpl_table_get_data_array (flux_FT, "FLUX");
	cpl_size nwave_ft = cpl_table_get_column_depth (vis_FT, "VISDATA");
	cpl_size nrow_ft  = cpl_table_get_nrow (vis_FT) / nbase;

	gravi_table_new_column (flux_FT, "TOTALFLUX", "e", CPL_TYPE_DOUBLE);
	double * total_flux_ft = cpl_table_get_data_double (flux_FT, "TOTALFLUX");
	
	CPLCHECK_MSG ("Cannot create columns");

	/* Total flux in [e] */
	for (cpl_size row = 0; row < nrow_ft * ntel; row ++) {
	  total_flux_ft[row] = cpl_array_get_mean (flux_ft[row]) * nwave_ft;
	}

	/* Smooth TOTALFLUX of the FT over few samples. Maybe make sure
	 * this is a constant frequency, not a constant nb of samples */
	gravi_table_smooth_column (flux_FT, "TOTALFLUX", "TOTALFLUX", 10, ntel);

	CPLCHECK_MSG("Cannot compute TOTALFLUX for FT");

	/* 
	 * (3) Create the signals for VIS_FT
	 */
	
	/* Create F1F2 and PHASE_REF column */
	gravi_vis_create_f1f2_ft (vis_FT, flux_FT);
	gravi_vis_create_phaseref_ft (vis_FT);

	CPLCHECK_MSG ("Cannot create signals for VIS_FT");
	
	/* 
	 *  (4) Create the signal for FLUX_SC
	 */
    gravi_flux_create_met_sc (flux_SC, vis_met);
	gravi_flux_create_fddlpos_sc (flux_SC, fddl_table);
	gravi_flux_create_totalflux_sc (flux_SC, flux_FT);

	CPLCHECK_MSG ("Cannot create the signal for FLUX_SC");

	/* 
	 *  (5) Create the signals for VIS_SC
	 */
	gravi_vis_create_pfactor_sc (vis_SC, flux_FT);
	gravi_vis_create_f1f2_sc (vis_SC, flux_SC);
	gravi_vis_create_met_sc (vis_SC, vis_met, oi_wavelengthsc);
	
	gravi_vis_create_vfactor_sc (vis_SC, oi_wavelengthsc,
								 vis_FT, oi_wavelengthft);

	CPLCHECK_MSG ("Cannot create signals for VIS_SC");

    /* 
     * Create QC for PFACTOR and VFACTOR
     */

    for (int base = 0; base < nbase; base++) {        
        char qc_name[100];
        
        sprintf (qc_name, "ESO QC VFACTOR%s_P%d AVG", GRAVI_BASE_NAME[base], pol+1);
        double vmean = gravi_table_get_column_mean (vis_SC, "V_FACTOR_WL", base, nbase);
        cpl_propertylist_update_double (p2vmred_header, qc_name, vmean);
        cpl_propertylist_set_comment (p2vmred_header, qc_name, "mean v-factor");
        
        sprintf (qc_name, "ESO QC PFACTOR%s_P%d AVG", GRAVI_BASE_NAME[base], pol+1);
        double pmean = gravi_table_get_column_mean (vis_SC, "P_FACTOR", base, nbase);
        cpl_propertylist_update_double (p2vmred_header, qc_name, pmean);
        cpl_propertylist_set_comment (p2vmred_header, qc_name, "mean p-factor");
    }

    /* 
     * If available, create the signal from ACQ camera
     */
    if (gravi_data_has_extension (p2vmred_data, GRAVI_OI_VIS_ACQ_EXT)) {
        
        cpl_table * vis_ACQ = gravi_data_get_table (p2vmred_data, GRAVI_OI_VIS_ACQ_EXT);
        gravi_signal_create_sync (vis_SC, 6, dit_sc, vis_ACQ, 4, "ACQ");
        gravi_signal_create_sync (flux_SC, 4, dit_sc, vis_ACQ, 4, "ACQ");

        gravi_vis_create_acq_sc (vis_SC, vis_ACQ);
        gravi_flux_create_acq_sc (flux_SC, vis_ACQ);
    }

	/* 
     * Compute FDDL = (FDDL_SC + FDDL_FT)/2 in [m]
     * from the POS_SC, POS_FT and the linearity coeficients 
     *
     * Compute the OPD_DISP = n(lbd) * ODD_MET 
     * from the OPD_MET and the dispersion coeficients 
     */
     cpl_table * disp_table =  disp_data ? gravi_data_get_table (disp_data, "DISP_MODEL") : NULL;
     
     gravi_flux_create_fddllin_sc (flux_SC, disp_table);
     gravi_vis_create_opddisp_sc (vis_SC, flux_SC, oi_wavelengthsc, disp_table, p2vmred_header, parlist);

     CPLCHECK_MSG ("Cannot compute the OPD_DISP");

	
     /* 
      * Compute the GDELAY of the FT and SC with the proper algorithm.
      * Critical since these quantities are used for debuging astrometry 
      */
    
     /* Recompute the GDELAY of the FT (probably useless since one use GDELAY_FT) */
     gravi_table_compute_group_delay (vis_FT, "VISDATA", "FLAG",
				      "GDELAY", oi_wavelengthft);
	
     /* Recompute the GDELAY of the SC (probably usefull) */
     gravi_table_compute_group_delay (vis_SC, "VISDATA", "FLAG",
				      "GDELAY", oi_wavelengthsc);

     /* Compute the GDELAY_FT of VISDATA_FT in SC table (critical) */
     gravi_table_compute_group_delay (vis_SC, "VISDATA_FT", "FLAG",
				      "GDELAY_FT", oi_wavelengthft);
	
     CPLCHECK_MSG ("Cannot compute the GDELAYs");

     /* Compute the SELF_REF phase reference for SC (first) */
     gravi_vis_create_phaseref_sc (vis_SC, oi_wavelengthsc, NULL, p2vmred_header, parlist);
	
     /* Compute the PHASE_REF reference for SC */
     gravi_vis_create_phaseref_sc (vis_SC, oi_wavelengthsc, oi_wavelengthft, p2vmred_header, parlist);

     /* Create the IMAGING_REF phase ref, need PHASE_REF */
     gravi_vis_create_imagingref_sc (vis_SC, oi_wavelengthsc, p2vmred_header, parlist);

     CPLCHECK_MSG ("Cannot compute the PHASE_REF");
  } /* End loop on pol */
  

  gravi_msg_function_exit(1);
  return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Create rejection flags P2VMREDUCED file
 * 
 * @param p2vmred_data:  the P2VMREDUCED data (modified in-place)
 * @param parlist:       parameter list of the recipe
 *
 * Create the rejection flags in the P2VMREDUCED file...
 * These computations are mandatory for the further averaging of the
 * frames into an OIFITS file with gravi_compute_vis
 *
 * Create column in vis_FT with the REJECTION_FLAT computed by comparing
 * the SNR and STATE to thresholds provided in parlist.
 *
 * Create column in vis_SC with the FRINGEDE_RATIO.
 *
 * Create column in vis_SC with the REJECTION_FLAT computed by comparing
 * the VFACTOR and the FRINGEDET_RATIO of each SC DIT to the specified
 * thresholds in parlist.
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_compute_rejection (gravi_data * p2vmred_data,
										const cpl_parameterlist * parlist)
{
  gravi_msg_function_start(1);
  cpl_ensure_code (p2vmred_data, CPL_ERROR_NULL_INPUT);
    
  char qc_name[100];
  int nbase = 6;

  /* Get header data */
  cpl_propertylist * p2vmred_header = gravi_data_get_header (p2vmred_data);
  CPLCHECK_MSG ("Cannot get header");

  /* create array to store rejection rate */
  cpl_array * rejected_array = cpl_array_new(nbase,CPL_TYPE_DOUBLE);
  cpl_array_fill_window_double (rejected_array, 0, nbase, 0.0);
    
  /* 
   * (1) Do the FT rejection flags 
   */
  if (gravi_data_has_type (p2vmred_data, "_FT") < 2 ) {
      cpl_msg_info (cpl_func, "Cannot compute rejection flags for FT (no FT data)");
  }
  else {
      cpl_msg_info (cpl_func, "Compute rejection flags for FT");
      
      /* Get the SNR thresholds from parameter */
      double threshold_SNR_ft = gravi_param_get_double (parlist, "gravity.signal.snr-min-ft");
      double threshold_STATE_ft = gravi_param_get_double (parlist, "gravity.signal.state-min-ft");
      double min_GSTATE_ft = gravi_param_get_double (parlist, "gravity.signal.global-state-min-ft");
      double max_GSTATE_ft = gravi_param_get_double (parlist, "gravity.signal.global-state-max-ft");

      cpl_msg_info (cpl_func,"SNR threshold to define fringe-detection in FT: %g", threshold_SNR_ft);
      cpl_msg_info (cpl_func,"STATE threshold for FT: %g", threshold_STATE_ft);
      cpl_msg_info (cpl_func,"GLOBAL_STATE threshold for FT: %g - %g", min_GSTATE_ft, max_GSTATE_ft);

      /* Loop on polarisations */
      int npol_ft = gravi_pfits_get_pola_num (p2vmred_header, GRAVI_FT);
      cpl_size nrow_ft = 0;
      for (int pol = 0; pol < npol_ft; pol++) {
          
          /* Get the pointer to data */
          cpl_table * vis_FT = gravi_data_get_oi_vis (p2vmred_data, GRAVI_FT, pol, npol_ft);
          double * snr = cpl_table_get_data_double (vis_FT, "SNR_BOOT");
          int  * state  = cpl_table_get_data_int (vis_FT, "STATE");
          int  * gstate = cpl_table_get_data_int (vis_FT, "OPDC_STATE");
          
          gravi_table_new_column (vis_FT, "REJECTION_FLAG", NULL, CPL_TYPE_INT);
          int * reject_flag_ft = cpl_table_get_data_int (vis_FT, "REJECTION_FLAG");

          CPLCHECK_MSG ("Cannot create columns");

          /* Loop on base and rows */
          nrow_ft  = cpl_table_get_nrow (vis_FT) / nbase;
          for (cpl_size row = 0; row < nrow_ft * nbase; row ++) {

              /* Rejection based on SNR (first bit) */
              int snr_bit = 0;
              int reject = 0;
              if ( snr[row] < threshold_SNR_ft )
              {
                  gravi_bit_set (reject_flag_ft[row], snr_bit);
                  reject = 1;
              }
              else
              {
                  gravi_bit_clear (reject_flag_ft[row], snr_bit);
              }
              
              /* Rejection based on STATE (2sd bit) */
              int state_bit = 1;
              if ( state[row] < threshold_STATE_ft ||
                   gstate[row] < min_GSTATE_ft ||
                   gstate[row] > max_GSTATE_ft )
              {
                  gravi_bit_set (reject_flag_ft[row], state_bit);
                  reject = 1;
              }
              else
              {
                  gravi_bit_clear (reject_flag_ft[row], state_bit);
              }
              
              /* add 1 if rejected to array */
              cpl_array_set_double(rejected_array,row%nbase,
                                 cpl_array_get_double(rejected_array,row%nbase,NULL)+reject );
          }
          
      } /* End loop on polarisation */
    
    /* normalize the rejection ratio as percent */
    cpl_array_multiply_scalar (rejected_array,100./(npol_ft*nrow_ft));
      
      /* store rejection ratio in QC */
    for (int base = 0; base < nbase; base++) {
        
        sprintf (qc_name, "ESO QC REJECTED RATIO FT%s", GRAVI_BASE_NAME[base]);
        double ratio = cpl_array_get_double (rejected_array, base, NULL);
        cpl_msg_info (cpl_func, "%s = %f", qc_name, ratio);
        cpl_propertylist_update_int (p2vmred_header, qc_name, round(ratio));
        cpl_propertylist_set_comment (p2vmred_header, qc_name, "[%] ratio of FT data flagged");
    }
      
      sprintf (qc_name, "ESO QC REJECTED RATIO FT");
      double ratio = cpl_array_get_mean (rejected_array);
      cpl_msg_info (cpl_func, "%s = %f", qc_name, ratio);
      cpl_propertylist_update_int (p2vmred_header, qc_name, round(ratio));
      cpl_propertylist_set_comment (p2vmred_header, qc_name, "[%] ratio of FT data flagged");
      cpl_array_fill_window_double (rejected_array, 0, nbase, 0.0);
      
  } /* End FT */
  
  /* 
   * (2) Ratio of FT valid inside each SC DIT
   */
  if (gravi_data_has_type (p2vmred_data, "_FT") < 1 ||
      gravi_data_has_type (p2vmred_data, "_SC") < 1) {
      cpl_msg_info (cpl_func, "Cannot compute tracking ratio for SC (no FT or SC data)");
  }
  else {
      cpl_msg_info (cpl_func, "Compute tracking ratio for each SC DIT");

      /* Loop on polarisations */
      int npol_sc = gravi_pfits_get_pola_num (p2vmred_header, GRAVI_SC);
      int npol_ft = gravi_pfits_get_pola_num (p2vmred_header, GRAVI_FT);
      for (int pol = 0; pol < npol_sc; pol++) {

          cpl_table * vis_SC = gravi_data_get_oi_vis (p2vmred_data, GRAVI_SC, pol, npol_sc);
          cpl_table * vis_FT = gravi_data_get_oi_vis (p2vmred_data, GRAVI_FT, CPL_MIN(pol,npol_ft), npol_ft);
          gravi_vis_create_lockratio_sc (vis_SC, vis_FT);
          CPLCHECK_MSG ("Cannot compute lockratio_sc");
      }
  }

  /* 
   * (3) Do the SC rejection flags 
   */
  if (gravi_data_has_type (p2vmred_data, "_SC") < 2 ) {
      cpl_msg_info (cpl_func, "Cannot compute rejection flags for SC (no SC data)");
  }
  else {
      cpl_msg_info (cpl_func, "Compute rejection flags for SC");
    
      /* Get the SC rejection parameters */
      double minlockratio_sc = gravi_param_get_double (parlist, "gravity.signal.tracking-min-sc");
      double minvfactor_sc = gravi_param_get_double (parlist, "gravity.signal.vfactor-min-sc");
      
      cpl_msg_info (cpl_func,"Fringe-detection ratio to discard frame on SC: %g", minlockratio_sc);
      cpl_msg_info (cpl_func,"vFactor threshold to discard frame on SC: %g", minvfactor_sc);

      double opd_pupil_max_sc = gravi_param_get_double (parlist, "gravity.signal.opd-pupil-max-sc");
      double opd_pupil_stddev_max_sc = gravi_param_get_double (parlist, "gravity.signal.opd-pupil-stddev-max-sc");
 
      cpl_msg_info (cpl_func,"OPD_PUPIL threshold (abs) to discard frame on SC: %g", opd_pupil_max_sc);
      cpl_msg_info (cpl_func,"OPD_PUPIL_STDDEV threshold to discard frame on SC: %g", opd_pupil_stddev_max_sc);

      /* Loop on polarisations */
      int npol_sc = gravi_pfits_get_pola_num (p2vmred_header, GRAVI_SC);
      cpl_size nrow_sc = 0;
      for (int pol = 0; pol < npol_sc; pol++) {

          /* Get the table of reduced data */
          cpl_table * vis_SC = gravi_data_get_oi_vis (p2vmred_data, GRAVI_SC, pol, npol_sc);
          double * vFactor_wl = cpl_table_get_data_double (vis_SC, "V_FACTOR_WL");
          double * fringedet_ftdit = cpl_table_get_data_double (vis_SC, "FRINGEDET_RATIO");

          CPLCHECK_MSG("Cannot load data");

          double * opd_metdit_pupil = NULL;
          if ( cpl_table_has_column (vis_SC, "OPD_MET_PUPIL") )
        	  opd_metdit_pupil = cpl_table_get_data_double (vis_SC, "OPD_MET_PUPIL");

          double * opd_metdit_pupil_stddev = NULL;
          if ( cpl_table_has_column (vis_SC, "OPD_MET_PUPIL_STDDEV") )
        	  opd_metdit_pupil_stddev = cpl_table_get_data_double (vis_SC, "OPD_MET_PUPIL_STDDEV");

          CPLCHECK_MSG("Cannot load data");

          gravi_table_new_column (vis_SC, "REJECTION_FLAG", NULL, CPL_TYPE_INT);
          int * reject_flag_sc = cpl_table_get_data_int (vis_SC, "REJECTION_FLAG");

          CPLCHECK_MSG ("Cannot create columns");

          /* Loop on base and row SC */
          nrow_sc  = cpl_table_get_nrow (vis_SC) / nbase;
          for (cpl_size row_sc = 0; row_sc < nrow_sc * nbase; row_sc++) {
              
              /* Rejection based in lockratio (first bit) */
              int lock_bit = 0;
              int reject = 0;
              if ( fringedet_ftdit[row_sc] < minlockratio_sc )
              {
                  gravi_bit_set (reject_flag_sc[row_sc], lock_bit);
                  reject = 1;
              }
              else
                  gravi_bit_clear (reject_flag_sc[row_sc], lock_bit);
              
              /* Rejection based in the white-light vFactor (second bit) */
              int vfactor_bit = 1;
              if ( vFactor_wl[row_sc] < minvfactor_sc )
              {
                  gravi_bit_set (reject_flag_sc[row_sc], vfactor_bit);
                  reject = 1;
              }
              else
              {
                  gravi_bit_clear (reject_flag_sc[row_sc], vfactor_bit);
              }

              /* Rejection based on OPD_PUPIL (third bit) */
              int opd_pupil_bit = 2;
    	      if ( opd_metdit_pupil ) {
    	    	  if ( fabs(opd_metdit_pupil[row_sc]) > opd_pupil_max_sc )
                  {
                      gravi_bit_set (reject_flag_sc[row_sc], opd_pupil_bit);
                      reject = 1;
                  }
                  else
                  {
                      gravi_bit_clear (reject_flag_sc[row_sc], opd_pupil_bit);
                  }
    	      }

              /* Rejection based on OPD_PUPIL_STDDEV (fourth bit) */
              int opd_pupil_stddev_bit = 3;
    	      if ( opd_metdit_pupil_stddev ) {
    	    	  if ( opd_metdit_pupil_stddev[row_sc] > opd_pupil_stddev_max_sc )
                  {
                      gravi_bit_set (reject_flag_sc[row_sc], opd_pupil_stddev_bit);
                      reject = 1;
                  }
                  else
                  {
                      gravi_bit_clear (reject_flag_sc[row_sc], opd_pupil_stddev_bit);
                  }
    	      }
              
              /* add 1 if rejected to array */
              cpl_array_set_double(rejected_array,row_sc%nbase,
                                 cpl_array_get_double(rejected_array,row_sc%nbase,NULL)+reject );
          }
      } /* End loop on pol */
      
      /* normalize the rejection ratio as percent */
      cpl_array_multiply_scalar (rejected_array,100./(npol_sc*nrow_sc));
        
        /* store rejection ratio in QC */
      for (int base = 0; base < nbase; base++) {
          
          sprintf (qc_name, "ESO QC REJECTED RATIO SC%s", GRAVI_BASE_NAME[base]);
          double ratio = cpl_array_get_double (rejected_array, base, NULL);
          cpl_msg_info (cpl_func, "%s = %f", qc_name, ratio);
          cpl_propertylist_update_int (p2vmred_header, qc_name, ratio);
          cpl_propertylist_set_comment (p2vmred_header, qc_name, "[%] ratio of SC data flagged");
      }
        
        sprintf (qc_name, "ESO QC REJECTED RATIO SC");
        double ratio = cpl_array_get_mean (rejected_array);
        cpl_msg_info (cpl_func, "%s = %f", qc_name, ratio);
        cpl_propertylist_update_int (p2vmred_header, qc_name, ratio);
        cpl_propertylist_set_comment (p2vmred_header, qc_name, "[%] ratio of SC data flagged");
  } /* End SC */

    cpl_array_delete(rejected_array);
  gravi_msg_function_exit(1);
  return CPL_ERROR_NONE;
}

/**@}*/
