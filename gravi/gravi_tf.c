/* $Id: gravi_tf.c,v 1.10 2014/11/12 15:10:40 nazouaoui Exp $
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
 * @defgroup gravi_tf
 */
/**@{*/

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
#include <math.h>
#include <time.h>
#include <complex.h>

#include "gravi_data.h"
#include "gravi_dfs.h"
#include "gravi_pfits.h"
#include "gravi_cpl.h"

#include "gravi_utils.h"

#include "gravi_vis.h"
#include "gravi_tf.h"

/*-----------------------------------------------------------------------------
                              Private prototypes
 -----------------------------------------------------------------------------*/

int gravi_array_set_invalid_negative (cpl_array * array);
cpl_error_code gravi_vis_flag_negative (cpl_table * oi_table,
										const char * data, const char *flag);
cpl_error_code gravi_vis_flag_invalid (cpl_table * oi_table,
									   const char * data, const char *flag);
char * gravi_calib_setupstring (gravi_data * data);
double gravi_visibility_UD (double uv, double diam, double lbd);
cpl_size gravi_get_row_in_cat (cpl_table * diam_table, double ra, double dec, double *separation);

/*-----------------------------------------------------------------------------
                                 Function code
 -----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------*/
/**
 * @brief Set invalid to all negative elements of an array.
 *
 * @param array: the cpl_array to modify in-place
 */
/*-----------------------------------------------------------------------------*/

int gravi_array_set_invalid_negative (cpl_array * array)
{
  cpl_ensure (array, CPL_ERROR_NULL_INPUT, -1);
  
  /* Check the TF is positive */
  if ( cpl_array_get_min ( array ) > 0.0 ) return 0;
  
  cpl_size indx, size = cpl_array_get_size (array);
  int nv = 0, num = 0;

  /* Loop on array */
  for ( indx = 0 ; indx < size ; indx ++ ) {
	if ( cpl_array_get (array, indx, &nv) < 0.0 ) {
	  cpl_array_set_invalid (array, indx);
	  num++;
	}
  }

  return num;
}

/*-----------------------------------------------------------------------------*/
/**
 * @brief Flag negative element of an OIFITS table
 * 
 * @param oi_table:  the OIFITS table to modify in-place
 * @param data:      the column name (ex: VIS2DATA)
 * @param flag:      the corresponding flag column (ex: FLAG)
 */
/*-----------------------------------------------------------------------------*/

cpl_error_code gravi_vis_flag_negative (cpl_table * oi_table,
										const char * data, const char *flag)
{
  gravi_msg_function_start(0);
  cpl_ensure_code (oi_table, CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (data, CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (flag, CPL_ERROR_ILLEGAL_OUTPUT);
  
  /* Get pointer to speed up */
  int nv = 0;
  cpl_size row, nrow = cpl_table_get_nrow (oi_table);
  cpl_array ** pdata = cpl_table_get_data_array (oi_table, data);
  cpl_array ** pflag = cpl_table_get_data_array (oi_table, flag);
  cpl_size indx, size = cpl_array_get_size (pdata[0]);

  CPLCHECK_MSG ("Cannot get data");

  /* Loop on row and index. Add to FLAG if data is NULL */
  for ( row = 0 ; row < nrow ; row ++ ) {
	for ( indx = 0 ; indx < size ; indx ++ ) {
	  if ( cpl_array_get (pdata[row], indx, &nv) < 0.0) {
		cpl_array_set (pflag[row], indx, cpl_array_get (pflag[row], indx, &nv) + 1 );
	  }
	}
  }

  gravi_msg_function_exit(0);
  return CPL_ERROR_NONE;
}

/*-----------------------------------------------------------------------------*/
/**
 * @brief Flag invalid element of an OIFITS table
 * 
 * @param oi_table:  the OIFITS table to modify in-place
 * @param data:      the column name (ex: VIS2DATA)
 * @param flag:      the corresponding flag column (ex: FLAG)
 */
/*-----------------------------------------------------------------------------*/

cpl_error_code gravi_vis_flag_invalid (cpl_table * oi_table,
									   const char * data, const char *flag)
{
  gravi_msg_function_start(0);
  cpl_ensure_code (oi_table, CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (data, CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (flag, CPL_ERROR_ILLEGAL_OUTPUT);
  
  /* Get pointer to speed up */
  int nv = 0;
  cpl_size row, nrow = cpl_table_get_nrow (oi_table);
  cpl_array ** pdata = cpl_table_get_data_array(oi_table, data);
  cpl_array ** pflag = cpl_table_get_data_array(oi_table, flag);
  cpl_size indx, size = cpl_array_get_size (pdata[0]);
  
  CPLCHECK_MSG ("Cannot get data");
  
  /* Loop on row and index. Add to FLAG if data is NULL */
  for ( row = 0 ; row < nrow ; row ++ ) {
        if (pdata[row]==NULL) continue;
	for ( indx = 0 ; indx < size ; indx ++ ) {
	  if ( !cpl_array_is_valid (pdata[row], indx) ) {
		cpl_array_set (pflag[row], indx, cpl_array_get (pflag[row], indx, &nv) + 1 );
	  }
	}
  }

  gravi_msg_function_exit(0);
  return CPL_ERROR_NONE;
}

/*-----------------------------------------------------------------------------*/
/**
 * @brief Build a unique setup string per calibratable setup.
 * 
 * @param data:      the input gravi_data 
 * 
 * @return setupstring: the allocated setupstring coding the setup
 * 
 * The setup string contains: INS.SPEC.RES, INS.POLA.MODE, FT.POLA.MODE
 * FT.DIT and SC.DIT. The setupstring shall be desallocated properly.
 */
/*-----------------------------------------------------------------------------*/

char * gravi_calib_setupstring (gravi_data * data)
{
  gravi_msg_function_start(0);
  cpl_ensure (data, CPL_ERROR_NULL_INPUT, NULL);

  /* Main HEADER */
  cpl_propertylist * hdr_data = gravi_data_get_header (data);
  
  /* Build the string */
  char* setupstring = cpl_sprintf( "%10s %s %s %.2fms %.2fs",
								   gravi_pfits_get_spec_res (hdr_data),
								   gravi_pfits_get_pola_mode (hdr_data, GRAVI_SC),
								   gravi_pfits_get_pola_mode (hdr_data, GRAVI_FT),
								   gravi_pfits_get_ft_period (hdr_data),
								   gravi_pfits_get_sc_dit (hdr_data) );

  cpl_msg_debug (cpl_func, "Get setup string: %s", setupstring);
  
  CPLCHECK_NUL ("Cannot compute the setupstring");
  
  gravi_msg_function_exit(0);
  return setupstring;
}

/*-----------------------------------------------------------------------------*/
/**
 * @brief Interpolate the TF at the time of the science observation for
 *        an amplitude quantity.
 * 
 * @param science:      the science data to be calibrated, inplace
 * @param science_tf:   an already allocated data (duplication of science)
 *                      in which the TF interpolated point are stored.
 * @param used_tf_data: list of the TF data to be interpolated
 * @param num_tf_data:  number of TF data
 * @param extName:      EXTNAME of the OIFITS extension to be calibrated
 * @param insName:      INSNAME of the OIFITS extension to be calibrated
 * @param ampName:      column name of the amplitude to be calibrated
 * @param ampErrName:   corresponding error column
 * 
 * All dataset shall be conformable in SETUP, but may contain
 * a different number of observation (rnow/nbase)
 * amp is calibrated as a real number: < AMP . W >
 * The weight W is computed as a mixture between the SNR of each
 * TF measurement and its time distance to the science observation.
 */ 
/*-----------------------------------------------------------------------------*/

cpl_error_code gravi_apply_tf_amp( gravi_data * science,
								   gravi_data * science_tf,
								   gravi_data ** used_tf_data,
								   int num_tf_data,
								   const char * extName,
								   const char * insName,
								   const char * ampName,
								   const char * ampErrName,
								   int nbase, double delta_t)
{
  gravi_msg_function_start(0);
  cpl_ensure_code (science,      CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (insName,      CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (extName,      CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (ampName,      CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (ampErrName,   CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (used_tf_data, CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (num_tf_data>0, CPL_ERROR_ILLEGAL_INPUT);
  
  cpl_msg_debug (cpl_func, "%s %s amp=%s ampErr=%s nbase=%i",extName,insName,ampName,ampErrName,nbase);

  int i, row_sc, row_cal, nv=0;

  /* Get correct table */
  cpl_table * sci_table = gravi_data_get_oi_table (science, extName, insName);
  cpl_table * sci_tf_table = (science_tf?gravi_data_get_oi_table (science_tf, extName, insName):NULL);

  /* Some generic info */
  int nrow_sc = cpl_table_get_nrow (sci_table);
  int nwave = cpl_table_get_column_depth (sci_table, ampName);

  CPLCHECK_MSG( "Cannot get data on SCIENCE" );
  
  /* Loop on row of SCIENCE */
  for (row_sc = 0; row_sc < nrow_sc; row_sc ++){
	double time_sci = cpl_table_get_double (sci_table, "TIME", row_sc, &nv);
	int base = row_sc % nbase;

	/* Init the TF */
	cpl_array * tf_mean = cpl_array_new (nwave, CPL_TYPE_DOUBLE);
	cpl_array_fill_window (tf_mean, 0, nwave, 0.0);
	
	/* Loop on all possible TF = loop on TF files and loop
	 * on rows in each of this TF file */
	double norm = 0.0;
	for (i = 0; i < num_tf_data; i++){
	  
	  cpl_table * tf_table = gravi_data_get_oi_table (used_tf_data[i], extName, insName);
	  int nrow_tf = cpl_table_get_nrow (tf_table);
	  CPLCHECK_MSG("Cannot get the table on TF");
	  
	  for (row_cal = base; row_cal < nrow_tf; row_cal += nbase){
		double time_tf = cpl_table_get_double (tf_table, "TIME", row_cal, &nv);

		/* Get the data of this TF */
		cpl_array * tf_data = cpl_array_duplicate ( cpl_table_get_array (tf_table, ampName, row_cal) );
		cpl_array * tf_err  = cpl_array_duplicate ( cpl_table_get_array (tf_table, ampErrName, row_cal) );
		CPLCHECK_MSG("Cannot get data on TF");

		/* Compute the mean error. Don't give added advantage
		   for 2% error, Idea from John Monnier */
		double sigma = cpl_array_get_median (tf_err);
		sigma = CPL_MAX (sigma, 0.02*fabs(cpl_array_get_median (tf_data)));
		sigma = CPL_MAX (sigma, 1e-10);
					
		/* Compute the weighted mean of TF for this baseline */
		double coeff = exp (-2 * fabs (time_sci - time_tf) / delta_t) / (sigma*sigma);
		cpl_array_multiply_scalar (tf_data, coeff);
		cpl_array_add (tf_mean, tf_data);
		norm += coeff;
		
		cpl_array_delete (tf_data);
		cpl_array_delete (tf_err);
		CPLCHECK_MSG( "Error while integrating the TF" );
	  }
	  /* End loop on row in this TF file */
	}
	/* End loop on TF files */
			
	/* Divide by the sum of the weights of used CAL */
	cpl_array_divide_scalar (tf_mean, norm);

	/* Check the TF is positive */
	// ninv = gravi_array_set_invalid_negative (tf_mean);
	// if ( ninv ) cpl_msg_info (cpl_func, "Invalidate %i negative channels over %i", ninv, nwave);

	/* Dump the TF mean in the output table of TF_SCI*/
	if (sci_tf_table) cpl_table_set_array (sci_tf_table, ampName, row_sc, tf_mean);
	
	/* Apply the TF to the science_calibrated data */
	cpl_array_divide (cpl_table_get_data_array (sci_table, ampName)[row_sc], tf_mean);
	
	/* Apply the TF to the error on the science_calibrated data -- FIXME: error on TF not propagated */
	cpl_array_divide (cpl_table_get_data_array (sci_table, ampErrName)[row_sc], tf_mean);

	/* Free tf_mean */
	cpl_array_delete (tf_mean);
  }
  /* End loop on row of SCIENCE */

  /* Flag invalid data (NULL data, NULL error, negative errors) */
  if (sci_tf_table) {
    gravi_vis_flag_invalid (sci_tf_table, ampName, "FLAG");
    gravi_vis_flag_negative (sci_tf_table, ampErrName, "FLAG");
    gravi_vis_flag_invalid (sci_tf_table, ampErrName, "FLAG");
  }

  gravi_vis_flag_invalid (sci_table, ampName, "FLAG");
  gravi_vis_flag_invalid (sci_table, ampErrName, "FLAG");
  gravi_vis_flag_negative (sci_table, ampErrName, "FLAG");
  
  gravi_msg_function_exit(0);
  return CPL_ERROR_NONE;
}

/*-----------------------------------------------------------------------------*/
/**
 * @brief Interpolate the TF at the time of the science observation for
 *        a phase quantity (in deg).
 * 
 * @param science:      the science data to be calibrated, inplace
 * @param science_tf:   an already allocated data (duplication of science)
 *                      in which the TF interpolated point are stored.
 * @param used_tf_data: list of the TF data to be interpolated
 * @param num_tf_data:  number of TF data
 * @param extName:      EXTNAME of the OIFITS extension to be calibrated
 * @param insName:      INSNAME of the OIFITS extension to be calibrated
 * @param phiName:      column name of the phase to be calibrated
 * @param phiErrName:   corresponding error column
 * 
 * All dataset shall be conformable in SETUP, but may contain
 * a different number of observation (rnow/nbase)
 * phi is calibrated as a phasor: arg{< exp(i PHI ) . W >}
 * The weight W is computed as a mixture between the SNR of each
 * TF measurement and its time distance to the science observation.
 */ 
/*-----------------------------------------------------------------------------*/

cpl_error_code gravi_apply_tf_phi( gravi_data * science,
 								   gravi_data * science_tf,
								   gravi_data ** used_tf_data,
								   int num_tf_data,
								   const char * extName,
								   const char * insName,
								   const char * phiName,
								   const char * phiErrName,
								   int nbase, double delta_t)
{
  gravi_msg_function_start(0);
  cpl_ensure_code (science,      CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (insName,      CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (extName,      CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (phiName,      CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (phiErrName,   CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (used_tf_data, CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (num_tf_data>0, CPL_ERROR_ILLEGAL_INPUT);
  
  cpl_msg_debug(cpl_func, "%s %s phi=%s phiErr=%s nbase=%i",extName,insName,phiName,phiErrName,nbase);

  int i, row_sc, row_cal, nv=0;

  /* Get correct table */
  cpl_table * sci_table = gravi_data_get_oi_table (science, extName, insName);
  cpl_table * sci_tf_table = (science_tf?gravi_data_get_oi_table (science_tf, extName, insName):NULL);

  /* Some generic info */
  int nrow_sc = cpl_table_get_nrow (sci_table);
  int nwave = cpl_table_get_column_depth (sci_table, phiName);

  CPLCHECK_MSG("Cannot get data on SCIENCE");
  
  /* Loop on row of SCIENCE */
  for (row_sc = 0; row_sc < nrow_sc; row_sc ++){
	double time_sci = cpl_table_get_double (sci_table, "TIME", row_sc, &nv);
	int base = row_sc % nbase;

	/* Init the TF */
	cpl_array * tf_mean = cpl_array_new (nwave, CPL_TYPE_DOUBLE_COMPLEX);
	cpl_array_fill_window_complex (tf_mean, 0, nwave, 0.0 + 0.0*I);
	
	/* Loop on all possible TF = loop on TF files and loop
	 * on rows in each of this TF file */
	for (i = 0; i < num_tf_data; i++){
	  
	  cpl_table * tf_table = gravi_data_get_oi_table (used_tf_data[i], extName, insName );
	  int nrow_tf = cpl_table_get_nrow (tf_table);
	  CPLCHECK_MSG("Cannot get the table on TF");
	  
	  for (row_cal = base; row_cal < nrow_tf; row_cal += nbase){
		double time_tf = cpl_table_get_double (tf_table, "TIME", row_cal, &nv);

		/* Get the data of this TF in the form exp( i phase ) */
		cpl_array * tf_data = gravi_array_cexp ( I*CPL_MATH_PI/180.0, cpl_table_get_array (tf_table, phiName, row_cal) );
		cpl_array * tf_err  = cpl_array_duplicate (cpl_table_get_array (tf_table, phiErrName, row_cal) );
		CPLCHECK_MSG("Cannot get data on TF");

		/* Compute the mean error. Don't give added advantage
		   for 0.2deg for phase, Idea from John Monnier */
		double sigma = CPL_MAX (cpl_array_get_median (tf_err), 0.2);
					
		/* Compute the weighted mean of TF for this baseline */
		double coeff = exp (-2 * fabs (time_sci - time_tf) / delta_t) / (sigma*sigma);
		cpl_array_multiply_scalar (tf_data, coeff);
		cpl_array_add (tf_mean, tf_data);
		
		cpl_array_delete (tf_data);
		cpl_array_delete (tf_err);
		CPLCHECK_MSG("Error while integrating the TF");
	  }
	  /* End loop on row in this TF file */
	}
	/* End loop on TF files */

	/* Apply the TF to the science_calibrated data 
	 * VISPHI = arg{  exp(i VISPHI/180*pi) / tf_mean } * 180/pi */
	cpl_array * sci_data = gravi_array_cexp ( I*CPL_MATH_PI/180.0, cpl_table_get_array (sci_table, phiName, row_sc) );
	cpl_array_divide (sci_data, tf_mean);
	cpl_array_arg (sci_data);
	cpl_array_multiply_scalar (sci_data, 180./CPL_MATH_PI);
	cpl_table_set_array (sci_table, phiName, row_sc, sci_data);
	
			
	/* Dump the TF mean in the output table of TF_SCI
	 * VIPHI = arg{ tf_mean } * 180 / pi */
	cpl_array_arg (tf_mean);
	cpl_array_multiply_scalar (tf_mean, 180./CPL_MATH_PI);
	if (sci_tf_table) cpl_table_set_array (sci_tf_table, phiName, row_sc, tf_mean );
	
	/* No need to change error -- FIXME: error on TF not propagated */

	/* Free data */
	cpl_array_delete (tf_mean);
	cpl_array_delete (sci_data);
  }
  /* End loop on row of SCIENCE */
  
  /* Flag invalid data (NULL data, NULL error, negative errors) */
  if (sci_tf_table) {
    gravi_vis_flag_invalid (sci_tf_table, phiName, "FLAG");
    gravi_vis_flag_invalid (sci_tf_table, phiErrName, "FLAG");
    gravi_vis_flag_negative (sci_tf_table, phiErrName, "FLAG");
  }
  gravi_vis_flag_invalid (sci_table, phiName, "FLAG");
  gravi_vis_flag_invalid (sci_table, phiErrName, "FLAG");
  gravi_vis_flag_negative (sci_table, phiErrName, "FLAG");
  
  gravi_msg_function_exit(0);
  return CPL_ERROR_NONE;
}


/*----------------------------------------------------------------------------*/
/**
 * @brief Computes the calibrated visibility from science a single
 * 	  	  data and several previously evaluated instrumental visibility
 * 
 * @param vis_data:   The science OIFITS data to be calibrated
 * @param tf_data:    The list of OIFITS transfer function data
 * @param num_tf:     The number of transfer function data
 * @param zero:       The astrometric zero-point (unused so far)
 * @param tf_science: The already allocated gravi_data in which
 *                    the interpolated TF at the time of the science 
 *                    is returned.
 * @param parlist:   The parameter list
 * 
 * @return The calibrated OIFITS data.
 */
/*----------------------------------------------------------------------------*/

gravi_data * gravi_calibrate_vis(gravi_data * vis_data, gravi_data ** tf_data, int num_tf,
							 gravi_data * zero, gravi_data * tf_science,
							 const cpl_parameterlist * parlist)
{
    gravi_msg_function_start(1);
    cpl_ensure (vis_data, CPL_ERROR_NULL_INPUT, NULL);
    cpl_ensure (tf_data,  CPL_ERROR_NULL_INPUT, NULL);
    cpl_ensure (parlist,  CPL_ERROR_NULL_INPUT, NULL);
    cpl_ensure (num_tf>0, CPL_ERROR_NULL_INPUT, NULL);

	int i;
	gravi_data * vis_calib;
	char * setup_science, * setup_tf;

	/* Check the inputs */
	cpl_ensure( (vis_data != NULL) && (tf_data != NULL), CPL_ERROR_NULL_INPUT, NULL );

	/* 
	 * Find out the TF files who have the same setup keywords
	 * with the input SCIENCE file 
	 */
	
	/* This will store the pointer to usefull TF data */
	gravi_data ** used_tf_data = cpl_malloc (sizeof( gravi_data *) * num_tf);
	cpl_msg_info (cpl_func,"Number of possible TF: %i", num_tf);

	/* Verbose the options */
	double delta_t = gravi_param_get_double (parlist, "gravity.viscal.delta-time-calib");
	cpl_msg_info (cpl_func, "Delta time to interpolate the TF: %f s (%f h)", delta_t, delta_t/3600.0);

	int force_calib = gravi_param_get_bool (parlist, "gravity.viscal.force-calib");
	cpl_msg_info (cpl_func,"Force calibration of the SCI by CALs: %s", force_calib?"T":"F");
	
	/* Get setup keywords of the input SCIENCE file */
	setup_science = gravi_calib_setupstring (vis_data);
	cpl_msg_info (cpl_func,"Setup of SCIENCE: %s", setup_science);
	
	CPLCHECK_NUL("Cannot build the setup string for SCIENCE");

	/* Loop on the TF files */
	int num_used_tf = 0;
	for (i = 0; i < num_tf; i++){

	  /* Get the setup string of this TF */
	  setup_tf = gravi_calib_setupstring (tf_data[i]);

	  /* Check if compatible with SCIENCE */
	  if (!(strcmp (setup_tf, setup_science )) ) {
		/* case same setup */
		used_tf_data[num_used_tf] = tf_data[i];
		num_used_tf ++;
		cpl_msg_info (cpl_func,"Setup of TF file: %s -> keep", setup_tf);
	  } else if (force_calib) {
		/* case different setups but forced */
		used_tf_data[num_used_tf] = tf_data[i];
		num_used_tf ++;
		cpl_msg_info (cpl_func,"Setup of TF file: %s -> keep  (force_calib)", setup_tf);
	  } else 
		/* case different setups */
		cpl_msg_info (cpl_func,"Setup of TF file: %s -> discard", setup_tf);
		
	  cpl_free (setup_tf);
	}
	/* End loop on TF files */

	if (num_used_tf == 0) {
		cpl_free (used_tf_data);
		cpl_error_set_message (cpl_func, CPL_ERROR_NULL_INPUT, "No calib file with the same keywords");
		return NULL;
	}

	/* Duplicate the data to create a calibrated dataset */
	vis_calib = gravi_data_duplicate (vis_data);

	/* Get the header */
	cpl_propertylist * hdr_data = gravi_data_get_header (vis_calib);
	
	/* For each type of data SC / FT */
	int type_data, ntype_data = 2;
	for (type_data = 0; type_data < ntype_data ; type_data ++) {

	  /* Loop on polarisation */
	  int pol, npol = gravi_pfits_get_pola_num( hdr_data, type_data );
	  for ( pol= 0 ; pol < npol ; pol++ ) {

		/* Calibrate the VIS2 as a real quantity */
		gravi_apply_tf_amp (vis_calib, tf_science, used_tf_data, num_used_tf,
							GRAVI_OI_VIS2_EXT,
							GRAVI_INSNAME(type_data, pol, npol),
							"VIS2DATA", "VIS2ERR", 6, delta_t);
		  
		CPLCHECK_NUL("Cannot apply tf to VIS2DATA");

		/* Calibrate the VISAMP as a real quantity --> to be discussed */
		gravi_apply_tf_amp (vis_calib, tf_science, used_tf_data, num_used_tf,
							GRAVI_OI_VIS_EXT,
							GRAVI_INSNAME(type_data, pol, npol),
							"VISAMP", "VISAMPERR", 6, delta_t);

		CPLCHECK_NUL("Cannot apply tf to VISAMP");

		/* Calibrate the VISPHI --> to be discussed for the astrometry */
		gravi_apply_tf_phi (vis_calib, tf_science, used_tf_data, num_used_tf,
							GRAVI_OI_VIS_EXT,
							GRAVI_INSNAME(type_data, pol, npol),
							"VISPHI", "VISPHIERR", 6, delta_t);
		
		CPLCHECK_NUL("Cannot apply tf to VISPHI");

		/* Calibrate the T3AMP as a scalar quantity --> to be discussed */
		gravi_apply_tf_amp (vis_calib, tf_science, used_tf_data, num_used_tf,
							GRAVI_OI_T3_EXT,
							GRAVI_INSNAME(type_data, pol, npol),
							"T3AMP", "T3AMPERR", 4, delta_t);

		CPLCHECK_NUL("Cannot apply tf to VISAMP");
		
		/* Calibrate the T3PHI as a phasor */
		gravi_apply_tf_phi (vis_calib, tf_science, used_tf_data, num_used_tf,
							GRAVI_OI_T3_EXT,
							GRAVI_INSNAME(type_data, pol, npol),
							"T3PHI", "T3PHIERR", 4, delta_t);
		
		CPLCHECK_NUL("Cannot apply tf to T3PHI");
		
		/* Calibrate the FLUX as a real quantity  --> not calibrated */
		// gravi_apply_tf_amp (vis_calib, tf_science, used_tf_data, num_used_tf,
		// 					GRAVI_OI_FLUX_EXT,
		// 					GRAVI_INSNAME(type_data, pol, npol),
		// 					"FLUX", "FLUXERR", 4, delta_t);
		
		// CPLCHECK_NUL("Cannot apply tf to FLUX");

	  }
	  /* End loop on polarisation */
	}
	/* End loop on data_type */

	/* Free */
	cpl_free (used_tf_data);
	cpl_free (setup_science);
	
	gravi_msg_function_exit(1);
	return vis_calib;
}


/*----------------------------------------------------------------------------*/
/**
 * @brief Compute the expected visibility from a UD model
 * 
 * @param uv           baseline lenght [m]
 * @param diam         diameter [mas]
 * @param lbd          wavelength [m]
 * 
 * @return   2 * J1 (pi.x) / (pi.x)  with x = uv * diam / lbd
 */
/*----------------------------------------------------------------------------*/

double gravi_visibility_UD (double uv, double diam, double lbd)
{
  if (lbd <=0) return 0.0;
  double x = CPL_MATH_PI * uv / lbd * (diam * 1e-3 / 3600 / 180 * CPL_MATH_PI);
  return ( (x<=0) ? 1.0 : 2.0 * j1 (x) / x );
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Get the row in a cpl_table matching closest RAJ2000 and DEJ2000
 * 
 * @param diam_table   cpl_table with RAJ2000 and DECJ2000 in 
 *                     '+DD MM SS.SS' and 'HH MM SS.SS' format
 * @param ra           ra to search in [rad]
 * @param dec          dec to search in [rad]
 * 
 * @return closest match if <5", or -2 if no match
 */
/*----------------------------------------------------------------------------*/

cpl_size gravi_get_row_in_cat (cpl_table * diam_table, double ra, double dec, double *separation)
{
  gravi_msg_function_start(0);
  cpl_ensure (diam_table, CPL_ERROR_NULL_INPUT, -1);
  
  /* Assume the data are on the first extension */
  const char ** p_ra  = cpl_table_get_data_string_const (diam_table, "RAJ2000");
  const char ** p_dec = cpl_table_get_data_string_const (diam_table, "DEJ2000");
  cpl_size nrow = cpl_table_get_nrow (diam_table);

  cpl_ensure (p_ra,   CPL_ERROR_ILLEGAL_INPUT, -1);
  cpl_ensure (p_dec,  CPL_ERROR_ILLEGAL_INPUT, -1);
  cpl_ensure (nrow>0, CPL_ERROR_ILLEGAL_INPUT, -1);

  /* Init search */
  cpl_size row0 = -2;
  double dis0 = 1e20;

  /* Loop on rows */
  for (cpl_size row = 0; row < nrow ; row ++) {

	/* RA and DEC of this row in [rad] */
	double c_ra  = gravi_ra_to_rad (p_ra[row]);
	double c_dec = gravi_dec_to_rad (p_dec[row]);

	/* Compute distance in [rad] 
	 * FIXME: we assume J2000 everywere */
	double dist = acos ( sin (c_dec) * sin (dec) + cos (c_dec) * cos (dec) * cos (ra - c_ra) );	

	/* Closest so far */
	if ( dist < dis0 ) {
	  row0 = row;
	  dis0 = dist;
	}
  } /* End loop on rows */

  /* Best separation in arcsec */
  if (separation != NULL) *separation = dis0 / CPL_MATH_PI * 180 * 3600;
  
  gravi_msg_function_exit(0);
  return row0;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief This function evaluates the transfer function from the observation
 *        of a reference star whose diameter can be determined.
 * 
 * @param vis_data:      The input OIFITS observation
 * @param diamcat_data:  The input catalog of calibrator diameters (optional)
 * 
 * @return The OIFITS transfer function.
 */
/*----------------------------------------------------------------------------*/

gravi_data * gravi_compute_tf (gravi_data * vis_data, gravi_data * diamcat_data)
{
    gravi_msg_function_start(1);
	cpl_ensure (vis_data, CPL_ERROR_NULL_INPUT, NULL);
  
	gravi_data * tf_data = NULL;
	cpl_propertylist * plist;
	cpl_table * tf_vistable;
	cpl_table * tf_vis2table;
	cpl_table * oi_wavelength;
	int nv, type_data;
	double diameter, diameter_err, lambda, r_vis, r_vis2, tf_v, tf_v2;
	char qc_name[90];

	/* Get the OIFITS table */
	tf_data = gravi_data_duplicate (vis_data);
	plist = gravi_data_get_header (tf_data);

	/* Loop on SC / FT */
	for (type_data = 0; type_data < 2; type_data ++) {
	    cpl_table * diam_table;
	    double sep = 99.0;
	    cpl_size row_cat = -1;

	    /* Search in catalogue. Actually this should only
	     * be for SINGLE since catalogue won't be accurate for dual */
	    if (diamcat_data) {
	      /* Get the matching row */
	      diam_table = gravi_data_get_table_x (diamcat_data, 0);
	      row_cat = gravi_get_row_in_cat (diam_table,
	    							  gravi_pfits_get_type_raep (plist, type_data),
	    							  gravi_pfits_get_type_decep (plist, type_data),
									  &sep);

		  /* Use it only if closest than 5" */
	      if ( (row_cat > -1) && (sep < 5.0) ) {
			cpl_msg_info (cpl_func, "Find match in DIAMETER_CAT (best sep=%.2f arcsec)", sep);
			diameter = cpl_table_get (diam_table, "UDDK", row_cat, NULL);
			diameter_err = cpl_table_get (diam_table, "e_LDD", row_cat, NULL);
	      } else {
			cpl_msg_warning (cpl_func, "No match in DIAMETER_CAT (best sep=%.2f arcsec), use the HEADER value instead", sep);
			diameter = gravi_pfits_get_diameter (plist, type_data);
			diameter_err = 0.15;
		  }
		  
	    } else {
	      cpl_msg_info (cpl_func, "No DIAMETER_CAT, use the HEADER value if any");
		  diameter = gravi_pfits_get_diameter (plist, type_data);
		  diameter_err = 0.15;
	    }

		/* Verbose about diameter value with a warning if weird value 
		 * Or a return NULL is value is zero */
		if ( diameter>=4.0 ) {
		  cpl_msg_warning (cpl_func,"Diameter for %s target: %.3f (+-%.3f) mas. Expected ?", GRAVI_TYPE(type_data), diameter, diameter_err);
		}
		else if (diameter <= 0.0 && !gravi_data_is_internal(vis_data)) {
		  cpl_msg_warning (cpl_func, "Diameter for %s target: 0.0 mas. Probably wrong, check parameter value", GRAVI_TYPE(type_data));
		  /* FREE (gravi_data_delete, tf_data);
		  return NULL;*/
		}
		else {
		  cpl_msg_info(cpl_func,"Diameter for %s target: %.3f (+-%.3f) mas", GRAVI_TYPE(type_data), diameter, diameter_err);
		}
		
		/* Loop on polarisations */
		int pol,npol = gravi_pfits_get_pola_num (plist, type_data);
		for (pol = 0; pol < npol; pol ++){

			oi_wavelength = gravi_data_get_oi_wave (vis_data, type_data, pol,npol);
			tf_vistable = gravi_data_get_oi_vis (tf_data, type_data, pol,npol);
			tf_vis2table = gravi_data_get_oi_vis2 (tf_data, type_data, pol,npol);

			if ((tf_vistable == NULL) || (tf_vis2table == NULL) || (oi_wavelength == NULL)){
			  FREE (gravi_data_delete, tf_data);
			  cpl_error_set_message (cpl_func, CPL_ERROR_NULL_INPUT,
									 "Missing OI_VIS or OI_VIS2 or OI_WAVELENGTH");
			  return NULL;
			}

			/* Construction of the data */

			/* Loop on row in the table -- warning, we here suppose the OI_VIS and
			 * OI_VIS2 have the same number of rows */
			int nrow = cpl_table_get_nrow (tf_vistable);
			int nwave = cpl_table_get_column_depth (tf_vis2table, "VIS2DATA" );
			for (cpl_size row = 0; row < nrow; row ++){

				/* Compute norm for the baseline in meters */
				r_vis  = sqrt (pow (cpl_table_get_double (tf_vistable, "UCOORD", row, &nv), 2) +
							   pow (cpl_table_get_double (tf_vistable, "VCOORD", row, &nv), 2));
				
				r_vis2 = sqrt (pow (cpl_table_get_double (tf_vis2table, "UCOORD", row, &nv), 2) +
							   pow (cpl_table_get_double (tf_vis2table, "VCOORD", row, &nv), 2));

				CPLCHECK_NUL ("Cannot extract UVCOORD");

				/* Compute the model visibility from uniform disk */
				cpl_array * model_vis = cpl_array_new (nwave, CPL_TYPE_DOUBLE);
				cpl_array * model_vis2 = cpl_array_new (nwave, CPL_TYPE_DOUBLE);

				/* Loop on wave */
				for (cpl_size wave = 0; wave < nwave ; wave ++){
					/* This computation is validated */
					lambda = cpl_table_get (oi_wavelength, "EFF_WAVE", wave, &nv);
					tf_v  = gravi_visibility_UD (r_vis, diameter, lambda);
					tf_v2 = pow (gravi_visibility_UD (r_vis2, diameter, lambda), 2.0);

					/* Set into the temporary array */
					cpl_array_set_double (model_vis, wave, tf_v);
					cpl_array_set_double (model_vis2, wave, tf_v2);
				}
				/* End loop on wave */

				/* Divide the observed visibilities by the model for VISAMP */
				cpl_array * tf_vis  = cpl_table_get_data_array (tf_vistable, "VISAMP")[row];
				cpl_array * tf_visErr  = cpl_table_get_data_array (tf_vistable, "VISAMPERR")[row];
				cpl_array_divide (tf_vis, model_vis);
				cpl_array_divide (tf_visErr, model_vis);
				
				/* Divide the observed visibilities by the model for VIS2 */
				cpl_array * tf_vis2 = cpl_table_get_data_array (tf_vis2table, "VIS2DATA")[row];
				cpl_array * tf_vis2Err = cpl_table_get_data_array (tf_vis2table, "VIS2ERR")[row];
				cpl_array_divide (tf_vis2, model_vis2);
				cpl_array_divide (tf_vis2Err, model_vis2);

				CPLCHECK_NUL("Cannot set the tf array");

				/* Compute the relative error on the TF at 2.2 microns 
				 * due to the diameter uncertainty */
				double errRelTF = fabs (( gravi_visibility_UD (r_vis, (diameter+diameter_err), 2.2e-6) -
										  gravi_visibility_UD (r_vis, (diameter-diameter_err), 2.2e-6) ) /
										gravi_visibility_UD (r_vis, diameter, 2.2e-6) );
				
				/* Note that VISPHI, T3PHI, T3AMP, VISDATA are already in
				 * output table since they are duplicated */

				/* Add QC params. FIXME: deal with QC parameter if multiple files
				 * - *mean* over file QC TF VISAMP_SC12_P1 AVG
				 * - *p2p* over file 
				 * FIXME: put both polarisation in a single QC */

				/* Get baseline name and id */
				int base=row%6;

				sprintf (qc_name, "ESO QC TF VISMOD_%s%s RELERR", GRAVI_TYPE(type_data), GRAVI_BASE_NAME[base]);
				gravi_pfits_update_double (plist, qc_name, errRelTF);
				cpl_propertylist_set_comment (plist, qc_name, "TF rel.err from diameter at 2.2um");
				CPLCHECK_NUL ("QC TF VISMOD RELERR");

				sprintf (qc_name, "ESO QC TF VISAMP_%s%s_P%d MED", GRAVI_TYPE(type_data), GRAVI_BASE_NAME[base], pol+1);
				gravi_pfits_update_double (plist, qc_name, cpl_array_get_median (tf_vis));
				cpl_propertylist_set_comment (plist, qc_name, "TF. VIS median over lbd.");
				CPLCHECK_NUL ("QC TF VIS AVG");

				sprintf (qc_name, "ESO QC TF VIS2_%s%s_P%d MED", GRAVI_TYPE(type_data), GRAVI_BASE_NAME[base], pol+1);
				gravi_pfits_update_double (plist, qc_name, cpl_array_get_median (tf_vis2));
				cpl_propertylist_set_comment (plist, qc_name, "TF. VIS2 median over lbd.");
				CPLCHECK_NUL ("QC TF VIS2 AVG");

				/* Free Memory */
				cpl_array_delete (model_vis2);
				cpl_array_delete (model_vis);
			}
			/* End loop on rows */
		}
		/* End loop on pol */

		
		/* Compute the transmission QC only if a
		 * valid match with catalog is found */
		if ( (row_cat > -1) && (sep < 5.0) ) {

		  /* Check if consistent setup (all ATs or all UTs) */
		  cpl_table * oi_array = gravi_data_get_table (vis_data, GRAVI_OI_ARRAY_EXT);
		  const char sta = *cpl_table_get_string (oi_array, "TEL_NAME", 0);
		  if (*cpl_table_get_string (oi_array, "TEL_NAME", 1) != sta ||
			  *cpl_table_get_string (oi_array, "TEL_NAME", 2) != sta ||
			  *cpl_table_get_string (oi_array, "TEL_NAME", 3) != sta)
		  {
			cpl_error_set_message (cpl_func, CPL_ERROR_ILLEGAL_INPUT, "AT/UT mode is not supported");
			return NULL;
		  }

		  /* Compute in [photons/s/m2] */
		  double diam = 1.0;
		  if      (sta=='U') { diam = 8.0; cpl_msg_info (cpl_func, "Use UTs"); }
		  else if (sta=='A') { diam = 1.8; cpl_msg_info (cpl_func, "Use ATs"); }
		  else               { diam = 1.0; cpl_msg_warning (cpl_func, "Cannot find the diameter of telescope"); }

		  /* Get the expected flux in [photons/s] */
		  double Kmag = cpl_table_get (diam_table, "Kmag", row_cat, NULL);
		  double flux0 = 1.71173e+09 * pow (10.0, -Kmag/2.5) * (CPL_MATH_PI * pow (diam/2, 2));
		  cpl_msg_info (cpl_func, "Use Kmag=%.2f for QC.TRANS", Kmag);

		  /* Loop on beam */
		  for (int tel = 0; tel < 4; tel++) {
			double flux = 0.0;

			/* Get the total observed flux in [e/s] */
			sprintf (qc_name, "ESO QC FLUXRATE_%s%i_P1 SUM",GRAVI_TYPE(type_data),tel+1);
			flux += gravi_pfits_get_double_silentdefault (plist, qc_name, 0.);
			sprintf (qc_name, "ESO QC FLUXRATE_%s%i_P2 SUM",GRAVI_TYPE(type_data),tel+1);
			flux += gravi_pfits_get_double_silentdefault (plist, qc_name, 0.);
			CPLCHECK_NUL ("Cannot get fluxrate");

			/* Add the QC parameter */
			sprintf (qc_name, "ESO QC TF TRANS_%s%i",GRAVI_TYPE(type_data),tel+1);
			gravi_pfits_update_double (plist, qc_name, flux / flux0 * 100.0);
			cpl_propertylist_set_comment (plist, qc_name, "[%] Total transmission");
			cpl_msg_info (cpl_func, "%s = %.2f%%", qc_name, flux / flux0 * 100.0);
			CPLCHECK_NUL ("QC TRANS");
		  }
		} /* End computation of TF */
		
	}
	/* End loop on type_data */

	gravi_msg_function_exit(1);
	return tf_data;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Fill QC parameters related to transfer function.
 * 
 * @param oi_vis:        The input OIFITS observation
 * @param diamcat_data:  The input catalog of calibrator diameters (optional)
 * 
 * This function evaluates the transfer function and create
 * QC parameters in the header of the input OIFITS observation.
 * These QC parameters are created only if the TF could be computed.
 * The TF OIFITS data itself is actually erased.
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_compute_tf_qc (gravi_data * oi_vis, gravi_data * diamcat_data)
{
  gravi_msg_function_start(0);
  cpl_ensure_code (oi_vis, CPL_ERROR_NULL_INPUT);
  
  cpl_errorstate prestate = cpl_errorstate_get();

  cpl_msg_info (cpl_func, "Compute the QC TF parameters");
  
  /* Compute TF */
  gravi_data * oi_tf = gravi_compute_tf (oi_vis, diamcat_data);

  /* If an error is catch when computing the QC parameters, dump this error but continue */
  if ( cpl_error_get_code() || oi_tf == NULL) {
	cpl_msg_warning (cpl_func, "Cannot compute the QC TF parameters for this observation... continue.");
	cpl_errorstate_set (prestate);
	return CPL_ERROR_NONE;
  }

  /* Copy to header */
  cpl_propertylist_copy_property_regexp (gravi_data_get_header (oi_vis),
										 gravi_data_get_header (oi_tf),
										 ".*QC TF.*", 0);
  FREE (gravi_data_delete,oi_tf);
  
  gravi_msg_function_exit(0);
  return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
  @brief   Compute the ZP data
 */
/*----------------------------------------------------------------------------*/

gravi_data * gravi_compute_zp (gravi_data ** vis_calib, int num_calib)
{
    gravi_msg_function_start(1);
    cpl_ensure (vis_calib, CPL_ERROR_NULL_INPUT, NULL);
    cpl_ensure (num_calib>0, CPL_ERROR_ILLEGAL_INPUT, NULL);
    
	gravi_data * zero_met = gravi_data_duplicate (vis_calib[0]);

    gravi_msg_function_exit(1);
	return zero_met;
}


/**@}*/
