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
 * @defgroup gravi_p2vmred  Compute the visibilities
 *
 * This module implement the functions to compute the coherent flux of all the
 * frames of the beam combiners (FT and SC). The main function is @c gravi_compute_p2vmred().
 * The algorithms involved in these fonction are description in the sections :
 * - algorithms/Extraction of the coherent fluxes and telescope fluxes via P2VM
 * - algorithms/Computation of SNR
 *
 */
/**@{*/

/*
 * History :
 *      11/01/2019 Fix Warnings , unused parameter :  npol_ft, npol_sc, p2vmreducedFlag ,init_type_data, end_type_data
 */

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

#include "gravi_p2vmred.h"
#include "gravi_vis.h"

/*-----------------------------------------------------------------------------
                              Private prototypes
 -----------------------------------------------------------------------------*/

cpl_table * gravi_create_oitarget_table (const cpl_propertylist *header,
                                         const char * mode);

cpl_table * gravi_create_oiarray_table (const cpl_table * array_geometry,
                                        int is_cal);

/*-----------------------------------------------------------------------------
                              Function code
 -----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
 * @brief Create the OI_TARGET table from the main header
 *
 * @param header      Main header of an observation
 * @return The OI_TARGET table.
 *
 * The SC refers to TARGET_ID 1 is mode 'single' and to
 * TARGET_ID 2 in mode 'dual'. The FT refers to TARGET_ID 1
 * is all cases.
 */
/*----------------------------------------------------------------------------*/

cpl_table * gravi_create_oitarget_table (const cpl_propertylist *header,
                                         const char * mode)
{
	/* Message and timer */
	gravi_msg_function_start(1);
	cpl_ensure (header, CPL_ERROR_NULL_INPUT, NULL);
	cpl_ensure (mode,   CPL_ERROR_NULL_INPUT, NULL);

	/* Check if this observation has a reference object in FT and a science object in SC */
	int n_target = 1;
	if (!(strcmp (mode, "gravi_dual"))) n_target = 2;
	else if (!(strcmp (mode, "gravi_single"))) n_target = 1;
	else {cpl_error_set_message(cpl_func, CPL_ERROR_NULL_INPUT, "Invalid mode (bug)");return NULL;}

    /* Create the table */
	cpl_table * oi_target = cpl_table_new (n_target);

	/* TARGET_ID */
	cpl_table_new_column (oi_target, "TARGET_ID", CPL_TYPE_INT);
	cpl_table_set_column_savetype (oi_target, "TARGET_ID", CPL_TYPE_SHORT);
	/* TARGET */
	cpl_table_new_column (oi_target, "TARGET", CPL_TYPE_STRING);
	/* EQUINOX */
	cpl_table_new_column (oi_target, "EQUINOX", CPL_TYPE_FLOAT);
	cpl_table_set_column_unit (oi_target, "EQUINOX", "yr");
	/* RAEP0 and DECEP0 */
	cpl_table_new_column (oi_target, "RAEP0", CPL_TYPE_DOUBLE);
	cpl_table_set_column_unit (oi_target, "RAEP0", "deg");
	cpl_table_new_column (oi_target, "DECEP0", CPL_TYPE_DOUBLE);
	cpl_table_set_column_unit (oi_target, "DECEP0", "deg");
	/* RA_ERR and DEC_ERR */
	cpl_table_new_column (oi_target, "RA_ERR", CPL_TYPE_DOUBLE);
	cpl_table_set_column_unit (oi_target, "RA_ERR", "deg");
	cpl_table_new_column (oi_target, "DEC_ERR", CPL_TYPE_DOUBLE);
	cpl_table_set_column_unit (oi_target, "DEC_ERR", "deg");
	/* SYSVEL, VELTYP, VELDEF */
	cpl_table_new_column (oi_target, "SYSVEL", CPL_TYPE_DOUBLE);
	cpl_table_set_column_unit (oi_target, "SYSVEL", "m/s");
	cpl_table_new_column (oi_target, "VELTYP", CPL_TYPE_STRING);
	cpl_table_new_column (oi_target, "VELDEF", CPL_TYPE_STRING);
	/* PMRA, PMDEC and associated errors */
	cpl_table_new_column (oi_target, "PMRA", CPL_TYPE_DOUBLE);
	cpl_table_set_column_unit (oi_target, "PMRA", "deg/yr");
	cpl_table_new_column (oi_target, "PMDEC", CPL_TYPE_DOUBLE);
	cpl_table_set_column_unit (oi_target, "PMDEC", "deg/yr");
	cpl_table_new_column (oi_target, "PMRA_ERR", CPL_TYPE_DOUBLE);
	cpl_table_set_column_unit (oi_target, "PMRA_ERR", "deg/yr");
	cpl_table_new_column (oi_target, "PMDEC_ERR", CPL_TYPE_DOUBLE);
	cpl_table_set_column_unit (oi_target, "PMDEC_ERR", "deg/yr");
	/* PARALLAX and PARA_ERR */
	cpl_table_new_column (oi_target, "PARALLAX", CPL_TYPE_FLOAT);
	cpl_table_set_column_unit (oi_target, "PARALLAX", "deg");
	cpl_table_new_column (oi_target, "PARA_ERR", CPL_TYPE_FLOAT);
	cpl_table_set_column_unit (oi_target, "PARA_ERR", "deg");
	/* SPECTYP */
	cpl_table_new_column (oi_target, "SPECTYP", CPL_TYPE_STRING);
	CPLCHECK_NUL("Cannot create columns in OI_TARGET");

	/* loop on target to fill the columns*/
	for (int ti=0; ti<n_target; ti++) {
		
		/* TARGET_ID=1 is for FT and TARGET_ID=2 is for SC (if any) */
		cpl_table_set_int (oi_target, "TARGET_ID", ti, ti+1);

		/* TARGET name, forced to match 16 chars */
		const char *name = (ti?gravi_pfits_get_sobj(header):gravi_pfits_get_robj(header));
		gravi_table_set_string_fixlen (oi_target, "TARGET", ti, name, 16);

		CPLCHECK_NUL("Cannot set Target");
		
		/* At the moment, we can't fill the following columns. We need to put something for
		 * standard compliance. Some third party software fail if those columns are missing.
		 * Should be checked whether better values are available. */

		/* EQUINOX, RAEP0 and DECEP0 at mean equinox, for FT or SC  */
		float equinox = (float)gravi_pfits_get_double_silentdefault (header, "EQUINOX", 0.);
		double raep = (ti?gravi_pfits_get_sobj_raep (header):gravi_pfits_get_robj_raep (header));
		double decp = (ti?gravi_pfits_get_sobj_decep (header):gravi_pfits_get_robj_decep (header));
        if (ti==0) gravi_dump_the_boss (raep, decp);
        
		double ra_err =  0.0;
		double dec_err = 0.0;
		raep *= CPL_MATH_DEG_RAD; // [deg]
		decp *= CPL_MATH_DEG_RAD; // [deg]
		CPLCHECK_NUL("Cannot get EQUINOX, RA and DEC");

		/* Some verbose */
		cpl_msg_info (cpl_func,"Found target '%s' with ra=%fd and dec=%fd at equinox %.2f", name, raep, decp, equinox);

		/* SYSVEL, VELTYP, VELDEF */
		double sysvel = gravi_pfits_get_double_silentdefault (header, "ESO FT ROBJ RADVEL", 0.);
		const char * veltyp = "UNKNOWN ";
		const char * veldef = "OPTICAL ";
		CPLCHECK_NUL("Cannot get SYSVEL, VELTYPE...");

		/* PMRA, PMDEC and associated errors in deg/years 
		 * These quantities are in as/years in the HEADER */
		double pmra      = gravi_pfits_get_pmra (header) / 3600.0; // [deg/year]
		double pmdec     = gravi_pfits_get_pmdec (header) / 3600.0; // [deg/year]
		double pmra_err  = 0.0;
		double pmdec_err = 0.0;
		CPLCHECK_NUL("Cannot get PMRA...");
		
		/* PARALLAX and PARA_ERR */
		float parallax = gravi_pfits_get_plx (header); // [as]
		float para_err = 0.0;
		CPLCHECK_NUL("Cannot get PARRALAX...");

		/* SPECTYP */
		const char * spectyp = "UNKNOWN         ";

		/* Fill common columns */
		cpl_table_set_float  (oi_target, "EQUINOX",  ti, equinox);
		cpl_table_set_double (oi_target, "RAEP0",    ti, raep);
		cpl_table_set_double (oi_target, "DECEP0",   ti, decp);
		cpl_table_set_double (oi_target, "RA_ERR",   ti, ra_err);
		cpl_table_set_double (oi_target, "DEC_ERR",  ti, dec_err);
		cpl_table_set_double (oi_target, "SYSVEL",   ti, sysvel);
		cpl_table_set_string (oi_target, "VELTYP",   ti, veltyp);
		cpl_table_set_string (oi_target, "VELDEF",   ti, veldef);
		cpl_table_set_double (oi_target, "PMRA",     ti, pmra);
		cpl_table_set_double (oi_target, "PMDEC",    ti, pmdec);
		cpl_table_set_double (oi_target, "PMRA_ERR", ti, pmra_err);
		cpl_table_set_double (oi_target, "PMDEC_ERR",ti, pmdec_err);
		cpl_table_set_float  (oi_target, "PARALLAX", ti, parallax);
		cpl_table_set_float  (oi_target, "PARA_ERR", ti, para_err);
		cpl_table_set_string (oi_target, "SPECTYP",  ti, spectyp);
	}
	/* End loop on target_id */
    
	gravi_msg_function_exit(1);
    return oi_target;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Create the OI_ARRAY table from the ARRAY_GEOMETRY
 * 
 * @param array_geometry    The input ARRAY_GEOMETRY table
 * @param is_cal            0/1, to select the verbosity level
 *
 * This function duplicate the ARRAY_GEOMETRY table into a OI_ARRAY table
 * and makes several checks of OIFITS standart. It also ensures that the
 * constructed OI_ARRAY has four rows, one per beam. Missing beam are
 * built as T1, T2...
 */
/*----------------------------------------------------------------------------*/

cpl_table * gravi_create_oiarray_table (const cpl_table * array_geometry,
                                        int is_cal)
{
    gravi_msg_function_start(1);
    cpl_ensure (array_geometry, CPL_ERROR_NULL_INPUT, NULL);

    cpl_table * oi_array = cpl_table_duplicate (array_geometry);
    
    /* Rename column STA to STAXYZ */
    if (cpl_table_has_column (oi_array, "STA") ) {
    	cpl_table_name_column (oi_array, "STA", "STAXYZ");
    }
    CPLCHECK_NUL ("Cannot rename STA into STAXYZ");
    
    /* Make sure there are 4 rows (one per beam) */
    cpl_size n_row = cpl_table_get_nrow (oi_array);
    if (n_row < 4) {
        
        if (is_cal) cpl_msg_info (cpl_func, "ARRAY_GEOMETRY has %lld rows instead of 4", n_row);
        else        cpl_msg_warning (cpl_func, "ARRAY_GEOMETRY has %lld rows instead of 4", n_row);
        cpl_table_set_size (oi_array, 4);
        
        cpl_array * zero_array = gravi_array_init_double (3,0.0);
    
        /* Fill these new rows with default values. */
        for ( int row=n_row; row<4; row++ ) {
            char name16[17];
            cpl_table_set_int (oi_array, "STA_INDEX", row, 100+row);
            sprintf (name16, "S%i", row);
            cpl_table_set_string (oi_array, "STA_NAME", row, name16);
            sprintf (name16, "T%i", row);
            cpl_table_set_string (oi_array, "TEL_NAME", row, name16);
            cpl_table_set_array (oi_array, "STAXYZ", row, zero_array);
            cpl_table_set_float (oi_array, "DIAMETER", row, 0.0);
            cpl_table_set_int (oi_array, "MNTSTA", row, 0);
            if (is_cal) cpl_msg_info (cpl_func,"Add STA_INDEX %i with TEL_NAME=%s in OI_ARRAY",100+row,name16);
            else        cpl_msg_warning (cpl_func,"Add STA_INDEX %i with TEL_NAME=%s in OI_ARRAY",100+row,name16);
            CPLCHECK_NUL ("Cannot insert rows");
        }

        FREE (cpl_array_delete, zero_array);
      
    }/* End case there are missing rows */
    
    /* Update TEL_NAME and STA_NAME to use 16 char for OIFITS compliance */
    for ( int row=0; row<cpl_table_get_nrow (oi_array); row++ ) {
        const char * tel_name  = cpl_table_get_string (oi_array, "TEL_NAME", row);
        gravi_table_set_string_fixlen (oi_array, "TEL_NAME", row, tel_name, 16);
        const char * sta_name  = cpl_table_get_string (oi_array, "STA_NAME", row);
        gravi_table_set_string_fixlen (oi_array, "STA_NAME", row, sta_name, 16);
    }
    CPLCHECK_NUL ("Cannot force 16 chars in OI_ARRAY");
    
    gravi_msg_function_exit(1);
    return oi_array;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Converts preprocessed data into coherent fluxes using the P2VM.
 * 
 * @param preproc_data   The spectrum data input
 * @param p2vm_map		  The p2vm table with with the values of the transmission,
 * 	  	  	  	  	  	  phase and coherence
 * @param mode            gravi_single / gravi_dual
 * @param parlist 	  	  The parameter list containing the variables defining
 * 	  	  	  	  	      the size of the profile
 * @param det_type        The detector to process. If GRAVI_DET_SC, then only
 *                        the science detector extensions will be processed.
 *                        GRAVI_DET_FT will do the same for FT detector
 *                        and GRAVI_DET_ALL will do it for both. 
 * 
 * @return  The data containing the OI_FLUX and OI_VIS tables who compute
 * 	        the visibilies of each acquisition. 
 *
 * \exception CPL_ERROR_NULL_INPUT input data is missing
 * \exception CPL_ERROR_ILLEGAL_INPUT input data are not compatible
 * 
 * It inverts the v2pm with a singular
 * value decomposition, and compute the product spectrum_values*p2vm.
 * The P2VM algorithm is applied for each SPECTRUM DIT and each wavelength.
 * The coherence fluxes are saved in OI_VIS tables, and the photometric
 * fluxes are saved in OI_FLUX tables.
 * This function also computes the u,v vectors, TIME...
 * This function also creates the OI_TARGET and OI_ARRAY tables.
 *
 * IMPORTANT NOTE:
 * To save memory space, the data in preproc_data are deleted 
 * by the routine, while looping on the rows.
 */
/*----------------------------------------------------------------------------*/

gravi_data * gravi_compute_p2vmred (gravi_data * preproc_data, gravi_data * p2vm_map,
                                    const char * mode, const cpl_parameterlist * parlist)
{
	/* Message and timer */
	gravi_msg_function_start(1);
	cpl_ensure (preproc_data, CPL_ERROR_NULL_INPUT, NULL);
	cpl_ensure (p2vm_map,      CPL_ERROR_NULL_INPUT, NULL);
	cpl_ensure (mode,          CPL_ERROR_NULL_INPUT, NULL);
	cpl_ensure (parlist,       CPL_ERROR_NULL_INPUT, NULL);

	char qc_name[90];
    int ntel = 4, nclo = 4, nbase = 6;
	int nv; /* npol_ft, npol_sc; */

	/* Timers */
	clock_t timer_ft = 0, timer_sc = 0, timer_other = 0, timer_start = 0;
	timer_start = clock();
	
	/* Global properties */
	cpl_propertylist * preproc_header = gravi_data_get_header (preproc_data);

	/* Flag if the p2vmreduced-file is requested */
	/* int p2vmreducedFlag = 0; */
    if ( (cpl_parameterlist_find_const (parlist, "gravity.dfs.p2vmred-file") != NULL &&
          gravi_param_get_bool (parlist,"gravity.dfs.p2vmred-file")) ||
         (cpl_parameterlist_find_const (parlist, "gravity.dfs.astro-file") != NULL &&
          gravi_param_get_bool (parlist,"gravity.dfs.astro-file")) ) {
        cpl_msg_info (cpl_func, "p2vmreduced-file is requested, will add computation time.");
        /* p2vmreducedFlag = 1; */
    } else {
        cpl_msg_info (cpl_func, "p2vmreduced-file is not requested.");
	}

	CPLCHECK_NUL("Cannot check the parameter");

	/* Get this flag if internal calib or on-sky */
	int is_cal = gravi_pfits_is_calib (preproc_header);

	/* Build the output gravi_data. */
	gravi_data * p2vmred_data = gravi_data_new (0);
    cpl_propertylist * p2vmred_header = gravi_data_get_header (p2vmred_data);
    
    /* Dump all header from RAW product */
	cpl_propertylist_append (p2vmred_header, preproc_header);

    /* Duplicate necessary tables in product */
    gravi_data_copy_ext (p2vmred_data, preproc_data, GRAVI_OI_WAVELENGTH_EXT);

	/* Create OI_TARGET from header and add it to product */
    cpl_table * oi_target = gravi_create_oitarget_table (preproc_header, mode);
	gravi_data_add_table (p2vmred_data, NULL, GRAVI_OI_TARGET_EXT, oi_target);
	CPLCHECK_NUL("Cannot create OI_TARGET");

    /* Duplicate ARRAY_GEOMETRY as OI_ARRAY */
    cpl_table * array_geo = gravi_data_get_table (preproc_data, GRAVI_ARRAY_GEOMETRY_EXT);
    CPLCHECK_NUL("Cannot get ARRAY_GEOMETRY");
    cpl_table * oi_array  = gravi_create_oiarray_table (array_geo, is_cal);
    CPLCHECK_NUL("Cannot convert ARRAY_GEOMETRY into OI_ARRAY");

    /* Build the header of OI_ARRAY from the header of OI_ARRAY */
    cpl_propertylist * oiarray_plist = gravi_data_get_plist (preproc_data, GRAVI_ARRAY_GEOMETRY_EXT);
    oiarray_plist = cpl_propertylist_duplicate (oiarray_plist);

    /* Set OI_ARRAY in output data */
    gravi_data_add_table (p2vmred_data, oiarray_plist,
                          GRAVI_OI_ARRAY_EXT, oi_array);
    CPLCHECK_NUL ("Cannot add OI_ARRAY");
    
    /* Update the keyword ARRAYX/Y/Z to ensure double */
    cpl_msg_info (cpl_func, "Update the ARRAYX/Y/Z keywords as double precision");
    gravi_pfits_ensure_double (oiarray_plist, "ARRAYX");
    gravi_pfits_ensure_double (oiarray_plist, "ARRAYY");
    gravi_pfits_ensure_double (oiarray_plist, "ARRAYZ");

	/* Read OPTICAL_TRAIN */
	cpl_table * optical_train_table = gravi_data_get_table (preproc_data, GRAVI_OPTICAL_TRAIN_EXT);

	/* Loop on lab input corresponding to GRAVITY to make sure all beams exists
     * If they don't, they are filled with fake value to macth the OI_ARRAY */
	int labInput[4] = {GRAVI_LABINPUT_1,GRAVI_LABINPUT_2,GRAVI_LABINPUT_3,GRAVI_LABINPUT_4};
	for (int lab = 0 ; lab < 4 ; lab++) {

	  /* Search if this LABINIT already exist */
	  cpl_table_unselect_all (optical_train_table);
	  if ( cpl_table_or_selected_int (optical_train_table, "VALUE2", CPL_EQUAL_TO, labInput[lab]) == 1 ) {
		continue;
	  }
	  
	  /* If not, add a rows and fill it with default values */
	  char name16[17];
	  int nrow = cpl_table_get_nrow (optical_train_table);
	  cpl_table_set_size (optical_train_table, nrow + 1);
	  cpl_table_set_int (optical_train_table, "VALUE2", nrow, labInput[lab]);
	  sprintf (name16, "T%i", nrow);
	  cpl_table_set_string (optical_train_table, "TEL_NAME", nrow, name16);
	  if (is_cal) cpl_msg_info (cpl_func,"Add LABINPUT %i with TEL_NAME=%s in OPTICAL_TRAIN",labInput[lab],name16);
	  else        cpl_msg_warning (cpl_func,"Add LABINPUT %i with TEL_NAME=%s in OPTICAL_TRAIN",labInput[lab],name16);
	  CPLCHECK_NUL("Cannot insert rows");
	} /* End loop in labInput */

	/*
	 * For each type of data SC / FT 
	 */
    /* int init_type_data = 1;
    int end_type_data = 1;

    if(det_type == GRAVI_DET_SC || det_type == GRAVI_DET_ALL)
        init_type_data = 0;
    if(det_type == GRAVI_DET_FT || det_type == GRAVI_DET_ALL)
        end_type_data = 1;
    */
	for (int type_data = 0; type_data < 2; type_data ++){

        /* Check if this data is present */
        if (!gravi_data_has_extension (preproc_data, GRAVI_SPECTRUM_DATA_EXT(type_data))) {
            cpl_msg_info (cpl_func,"No data for %s, skip it", GRAVI_TYPE(type_data));
            continue;
        }
	  
		/* Get the table and property list utilities */
		cpl_table * p2vm_table = gravi_data_get_p2vm_data (p2vm_map, type_data);
		cpl_table * spectrum_table = gravi_data_get_spectrum_data (preproc_data, type_data);
		cpl_table * detector_table = gravi_data_get_imaging_detector (preproc_data, type_data);
		cpl_table * detector_p2vm = gravi_data_get_imaging_detector (p2vm_map, type_data);
        CPLCHECK_NUL ("Cannot get data");

        /* Verify the data and P2VM are conformable */
        cpl_ensure (cpl_table_get_nrow (detector_table) ==
                    cpl_table_get_nrow (detector_p2vm),
                    CPL_ERROR_ILLEGAL_INPUT, NULL);
        
		/* Read nrow and nregion */
		cpl_size nrow = cpl_table_get_nrow (spectrum_table);
        int n_detregion = cpl_table_get_nrow (detector_table);

		/* Guess the number of polarisation by looking at the number of region 
         * Save the number of polarisation for further use */
		int npol = (n_detregion>24 ? 2 : 1);
		/*if (type_data == GRAVI_FT) npol_ft = npol;
		if (type_data == GRAVI_SC) npol_sc = npol; */

		/* FIXME: Here we assume the two polarisation (if any)
		 * have the same number of channels ? */
		cpl_size nwave = gravi_spectrum_get_nwave (spectrum_table);
		
		/* 
		 * Loop on polarisations
		 */
		for (int p = 0; p < npol; p ++) {
            
			cpl_msg_info(cpl_func, "Visibility extraction: polarisation %d over %d for %s",
                         p+1, npol, type_data==GRAVI_FT?"FT":"SC");
			
			timer_other += (clock() - timer_start);
			timer_start = clock();
			
			/* 
             *
             * Create output tables and set them in output data
             * 
             */
		    cpl_msg_info(cpl_func, "Create and set output tables...");

            /* Create table */
			cpl_table * oi_vis = gravi_table_oi_create (nwave, nrow, GRAVI_OI_VIS_EXT);
			cpl_table * oi_t3 = gravi_table_oi_create (nwave, nrow, GRAVI_OI_T3_EXT);
 			cpl_table * oi_flux = gravi_table_oi_create (nwave, nrow, GRAVI_OI_FLUX_EXT);

			/* Delete useless columns in this product */
			cpl_table_erase_column (oi_vis, "VISAMP");
			cpl_table_erase_column (oi_vis, "VISPHI");
			cpl_table_erase_column (oi_vis, "VISAMPERR");
			cpl_table_erase_column (oi_vis, "VISPHIERR");
			cpl_table_erase_column (oi_vis, "RVIS");
			cpl_table_erase_column (oi_vis, "RVISERR");
			cpl_table_erase_column (oi_vis, "IVIS");
			cpl_table_erase_column (oi_vis, "IVISERR");
		
			cpl_table_erase_column (oi_t3, "T3AMP");
			cpl_table_erase_column (oi_t3, "T3PHI");
			cpl_table_erase_column (oi_t3, "T3AMPERR");
			cpl_table_erase_column (oi_t3, "T3PHIERR");
			cpl_table_erase_column (oi_t3, "U1COORD");
			cpl_table_erase_column (oi_t3, "V1COORD");
			cpl_table_erase_column (oi_t3, "U2COORD");
			cpl_table_erase_column (oi_t3, "V2COORD");

			/* Set table in p2vmred */
			cpl_propertylist * oi_plist = cpl_propertylist_new ();
            cpl_propertylist_copy_property (oi_plist, oiarray_plist, "ARRNAME");
            cpl_propertylist_copy_property (oi_plist, preproc_header, "DATE-OBS");
			cpl_propertylist_update_int (oi_plist, "OI_REVN", 1);
			cpl_propertylist_update_int (oi_plist, "NWAVE", nwave);
			cpl_propertylist_update_string (oi_plist, "INSNAME", GRAVI_INSNAME(type_data,p,npol));
			cpl_propertylist_update_int (oi_plist, "EXTVER", GRAVI_EXTVER(type_data,p,npol));
			
			gravi_data_add_table (p2vmred_data, oi_plist,
                                  GRAVI_OI_VIS_EXT, oi_vis);

            oi_plist = cpl_propertylist_duplicate (oi_plist);

			gravi_data_add_table (p2vmred_data, oi_plist,
                                  GRAVI_OI_T3_EXT, oi_t3);

            oi_plist = cpl_propertylist_duplicate (oi_plist);

			gravi_data_add_table (p2vmred_data, oi_plist,
                                  GRAVI_OI_FLUX_EXT, oi_flux);

			cpl_msg_info(cpl_func, "Total time to create tables: %.4f s", (double)(clock()-timer_start)/(double)CLOCKS_PER_SEC );

			/*
			 *
			 * Compute the P2VM and the V2PM for all wavelengths
			 *
			 */
			timer_other += (clock() - timer_start);
			timer_start = clock();
			
		    cpl_msg_info (cpl_func, "Compute the invers V2P2M -> P2VM matrix...");
			
			/* Construction of output tables */
			cpl_matrix ** p2vm = cpl_calloc (nwave,sizeof (cpl_matrix*));
			cpl_matrix ** v2pm = cpl_calloc (nwave,sizeof (cpl_matrix*));

            /* Get the list of regions having this pol */
            int all_region[48];
            int n_region = 0;
            for (int detreg = 0; detreg < n_detregion; detreg++) {
                if (gravi_region_get_pol (detector_p2vm, detreg) != p) continue;
                all_region[n_region] = detreg;
                n_region++;
            }

            /* Get direct pointer to arrays */
            const cpl_array ** trans, ** coh, ** phase;
            trans = cpl_table_get_data_array_const (p2vm_table, TRANSMISSION);
            coh = cpl_table_get_data_array_const (p2vm_table, COHERENCE);
            phase = cpl_table_get_data_array_const (p2vm_table, PHASE);
            CPLCHECK_NUL ("Cannot load the p2vm data");

            /* Allocate memory for direct access to P2VM */
            double ** pP2VM = cpl_malloc (sizeof(double*) * nwave);
            double ** pV2PM = cpl_malloc (sizeof(double*) * nwave);
                
			/* Loop on wave */
			for (cpl_size wave = 0; wave < nwave; wave ++) {

				/* Construction of the v2pm matrix */
				v2pm[wave] = cpl_matrix_new (n_region, 16);

                /* Loop on region */
				for (cpl_size region = 0; region < n_region; region++) {
                    int detregion = all_region[region];
                    
					/* Replace the four first columns of the v2pm by
					 * the transmission */
					for (int i = 0; i < 4; i++){
						cpl_matrix_set (v2pm[wave], region, i,
								cpl_array_get (trans[detregion], wave + i * nwave, &nv));
					}

					/* Replace the twelve other columns of the v2pm by
					 * the real part and the imaginary part of the coherence
					 * and the phase */
					for (int i = 0; i < 6; i++) {
						cpl_matrix_set (v2pm[wave], region, i + 4,
							cpl_array_get (coh[detregion], wave + i * nwave, &nv) *
							  cos(cpl_array_get (phase[detregion], wave + i * nwave, &nv)));
						cpl_matrix_set (v2pm[wave], region, i + 10,
							cpl_array_get (coh[detregion], wave + i * nwave, &nv) *
							  sin(cpl_array_get (phase[detregion], wave + i * nwave, &nv)));
					}
				}
				CPLCHECK_NUL("Cannot fill the V2PM");
				
				/* Ensure the V2PM is flux conservative */
				cpl_matrix_multiply_scalar (v2pm[wave], 1./n_region);

				/* Compute the matrix inversion of the v2pm using the
				 * singular value decomposition method */
				p2vm[wave] = gravi_matrix_invertSV_create (v2pm[wave]);

                /* Keep a pointer to the data */
                pP2VM[wave] = cpl_matrix_get_data (p2vm[wave]);
                pV2PM[wave] = cpl_matrix_get_data (v2pm[wave]);
				
				CPLCHECK_NUL ("Cannot invers V2PM");
			}
			/* End loop on wavelength */

			cpl_msg_info (cpl_func, "Total time to invers matrix: %.4f s", (double)(clock()-timer_start)/(double)CLOCKS_PER_SEC );

			/*
			 * Fill the STA_INDEX and TIME.
			 */
		    cpl_msg_info (cpl_func, "Fill the STA_INDEX and TIME.");
			timer_other += (clock() - timer_start);
			timer_start = clock();
			
			/* Wrap the TIME column into an array. TIME is in [us] from
             * the PCR.ACQ.START (start of RMN recording) */
			cpl_array * times = cpl_array_wrap_int (cpl_table_get_data_int (spectrum_table, "TIME"),
                                                    cpl_table_get_nrow (spectrum_table));
            cpl_ensure (times, CPL_ERROR_ILLEGAL_INPUT, NULL);
            
			/* Compute the MJDs for each row, by adding PCR.ACQ.START
             * (in MJD) to the TIME column */
			cpl_array * mjds = cpl_array_cast (times, CPL_TYPE_DOUBLE);
			cpl_array_divide_scalar (mjds, 86400.E6);
            double mjd0 = gravi_convert_to_mjd (gravi_pfits_get_start_prcacq (preproc_header));
			cpl_array_add_scalar (mjds, mjd0);

			/* We keep the TIME column of the P2VMRED in [us] from
             * the PCR.ACQ.START, in order to correlate with RMN tables.
             * Thus this TIME column is *not* at the OIFITS standart */
			cpl_table_set_column_unit (oi_vis,  "TIME", "us");
			cpl_table_set_column_unit (oi_t3,  "TIME", "us");
			cpl_table_set_column_unit (oi_flux, "TIME", "us");

			/* Fill STA_INDEX and times for OI_VIS */
			cpl_array * sta_index = cpl_array_new (2, CPL_TYPE_INT);
            
			for (int base=0; base < nbase; ++base) {
			    /* Build sta_index */
                int sta0 = gravi_sta_index(GRAVI_BASE_TEL[base][0]+1, optical_train_table, oi_array);
                int sta1 = gravi_sta_index(GRAVI_BASE_TEL[base][1]+1, optical_train_table, oi_array);
				cpl_array_set_int (sta_index, 0, sta0);
				cpl_array_set_int (sta_index, 1, sta1);
                CPLCHECK_NUL ("Cannot find the sta_index");

				/* loop on rows */
				for (int row=0; row < nrow; ++row) {
					int idx = row * nbase + base;
					cpl_table_set_array  (oi_vis, "STA_INDEX", idx, sta_index);
					cpl_table_set (oi_vis, "TIME", idx, cpl_array_get (times, row, &nv));
					cpl_table_set (oi_vis, "MJD", idx, cpl_array_get (mjds, row, &nv));
				}
			} /* End loop on base */

			cpl_array_delete (sta_index);
			CPLCHECK_NUL ("Cannot fill sta_index or time in OI_VIS");

			/* Fill STA_INDEX and times for OI_T3 */
			sta_index = cpl_array_new (3, CPL_TYPE_INT);
            
			for (int closure=0; closure < nclo; ++closure) {
			    /* Build sta_index */
                int sta0 = gravi_sta_index(GRAVI_CLO_TEL[closure][0]+1, optical_train_table, oi_array);
                int sta1 = gravi_sta_index(GRAVI_CLO_TEL[closure][1]+1, optical_train_table, oi_array);
                int sta2 = gravi_sta_index(GRAVI_CLO_TEL[closure][2]+1, optical_train_table, oi_array);
				cpl_array_set_int (sta_index, 0, sta0);
				cpl_array_set_int (sta_index, 1, sta1);
				cpl_array_set_int (sta_index, 2, sta2);
                CPLCHECK_NUL ("Cannot find the sta_index");

				/* loop on rows */
				for (int row=0; row < nrow; ++row) {
					int idx = row * nclo + closure;
					cpl_table_set_array  (oi_t3, "STA_INDEX", idx, sta_index);
					cpl_table_set (oi_t3, "TIME", idx, cpl_array_get (times, row, &nv));
					cpl_table_set (oi_t3, "MJD", idx, cpl_array_get (mjds, row, &nv));
				}
			} /* End loop on closures */

			cpl_array_delete (sta_index);
			CPLCHECK_NUL ("Cannot fill sta_index or time in OI_T3");

			/* Fill STA_INDEX and times for OI_FLUX */
			for (int tel = 0; tel < ntel; tel++){
				int sta0 = gravi_sta_index(tel+1, optical_train_table, oi_array);
				for (int row=0; row < nrow; ++row) {
					int idx = row * ntel + tel;
					cpl_table_set_int (oi_flux, "STA_INDEX", idx, sta0);
					cpl_table_set (oi_flux, "TIME", idx, cpl_array_get (times, row, &nv));
					cpl_table_set (oi_flux, "MJD", idx, cpl_array_get (mjds, row, &nv));
				}
			} /* End loop on tel */
            
			cpl_array_delete (mjds);
			cpl_array_unwrap (times);
			CPLCHECK_NUL ("Cannot fill sta_index or time in OI_FLUX");

			/* Fill the exposure time */
			double exptime = gravi_pfits_get_dit (gravi_data_get_header (preproc_data), type_data);
			cpl_table_fill_column_window  (oi_vis,  "INT_TIME", 0, nrow * nbase, exptime);
			cpl_table_fill_column_window  (oi_t3,  "INT_TIME", 0, nrow * nclo, exptime);
			cpl_table_fill_column_window  (oi_flux, "INT_TIME", 0, nrow * ntel, exptime);
			CPLCHECK_NUL ("Cannot fill exptime");

            /* target_id is 1 unless type_data==SC and mode==dual_field 
             * Same definition applies when building the OI_TARGET */
            int target_id = (!strcmp(mode, "gravi_dual") && (type_data==GRAVI_SC) )?2:1;
			cpl_table_fill_column_window_int (oi_vis, "TARGET_ID", 0, nrow * nbase, target_id);
			cpl_table_fill_column_window_int (oi_t3, "TARGET_ID", 0, nrow * nclo, target_id);
			cpl_table_fill_column_window_int (oi_flux, "TARGET_ID", 0, nrow * ntel, target_id);
			CPLCHECK_NUL ("Cannot fill target_id");
			
			cpl_msg_info (cpl_func, "Total time to fill STA_INDEX and TIME: %.4f s",
                          (double)(clock()-timer_start)/(double)CLOCKS_PER_SEC );

			/* 
			 * 
			 * Extract visibility with p2vm
			 * 
			 */
			
			cpl_msg_info (cpl_func,"Apply p2vm to the data...");
			timer_other += (clock() - timer_start);
			timer_start = clock();

			/* Get the pointers to the input data */
			double ** pReg = cpl_malloc (n_region * sizeof(double*));
			double ** pErr = cpl_malloc (n_region * sizeof(double*));
            cpl_array *** pRegArr = cpl_malloc (n_region * sizeof(cpl_array**));
            cpl_array *** pErrArr = cpl_malloc (n_region * sizeof(cpl_array**));
			for (int reg = 0; reg < n_region; reg++) {
			  pRegArr[reg] = cpl_table_get_data_array (spectrum_table, GRAVI_DATA[all_region[reg]]);
			  pErrArr[reg] = cpl_table_get_data_array (spectrum_table, GRAVI_DATAERR[all_region[reg]]);
			}
			CPLCHECK_NUL ("Cannot get data");

			/* Get the pointers to the output data */
			cpl_array** tFlux    = cpl_table_get_data_array (oi_flux, "FLUX");
			cpl_array** tFluxErr = cpl_table_get_data_array (oi_flux, "FLUXERR");
			cpl_array** tVis    = cpl_table_get_data_array (oi_vis, "VISDATA");
			cpl_array** tVisErr = cpl_table_get_data_array (oi_vis, "VISERR");
			CPLCHECK_NUL ("Cannot get data");

			/* Create column for chi2 */
			gravi_table_new_column_array (oi_flux, "CHI2", NULL, CPL_TYPE_DOUBLE, nwave);
			cpl_array** tChi2 = cpl_table_get_data_array (oi_flux, "CHI2");
			int ndof = n_region - 16;
			
			/* Temporary matrix output memory */
			double* pOut    = cpl_malloc (16 * nwave * sizeof(double));
			double* pOutVar = cpl_malloc (16 * nwave * sizeof(double));
			double* pChi2   = cpl_malloc (nwave * sizeof(double));

            /* Quantities to test flux conservation */
			double full_flux_reg = 0.0, full_flux_tel = 0.0;

			cpl_msg_debug (cpl_func, "Matrix multiplication");

			/* Loop on the frames */
            for (cpl_size row = 0; row < nrow; row ++) {

                /* Get pointers to the input data of this row */
                for (int reg = 0; reg < n_region; reg++) {
                    pReg[reg] = cpl_array_get_data_double (pRegArr[reg][row]);
                    pErr[reg] = cpl_array_get_data_double (pErrArr[reg][row]);
                }

                /* Loop on  wavelength */
                for (cpl_size wave = 0 ; wave < nwave ; wave++ ) {

                    /* Integration the input flux of the row */
                    for (int reg = 0; reg < n_region; reg++) full_flux_reg += pReg[reg][wave];
                    
                    /* Matrix multiplication, loop on outputs and regions 
                     * We neglect the input correlation and don't compute
                     * the output correlations. FIXME: we need to compute
                     * the output covariance of {R,I}  */
                    for (int out=0; out<16; out++) {
                        pOut[out*nwave+wave] = 0.0;
                        pOutVar[out*nwave+wave] = 0.0;
                        for (cpl_size reg = 0; reg < n_region; reg++) {
                            pOut[out*nwave+wave]    += pReg[reg][wave] * pP2VM[wave][out*n_region+reg];
                            pOutVar[out*nwave+wave] += gravi_pow2 (pErr[reg][wave] * pP2VM[wave][out*n_region+reg]);
                        }
                    } /* End outputs and regions */

                    /* We compute the reduced chi2. We loop on regions
                     * to recompute the expected value, and accumulate the chi2 */
                    pChi2[wave] = 0.0;
                    for (cpl_size reg = 0; reg < n_region; reg++) {
                        double value = 0.0;
                        for (int out=0; out<16; out++)
                            value += pV2PM[wave][reg*16+out] * pOut[out*nwave+wave];
                        pChi2[wave] += gravi_pow2 ((value-pReg[reg][wave]) / pErr[reg][wave]) / ndof;
                    }

                    /* Integration the output flux of the row */
                    for (int tel = 0; tel < 4; tel++) full_flux_tel += pOut[tel*nwave+wave];
                    
                } /* End loop on wavelengths */

                /* Set CHI2  (wrap is the fastest to create an array).
                   Same for all telescope since all-together */
                for (int tel = 0; tel < ntel; tel++){
                    double * data = cpl_malloc (nwave * sizeof(double));
                    for (cpl_size wave = 0 ; wave < nwave ; wave++ ) 
                        data[wave] = pChi2[wave];
                    tChi2[row*ntel+tel] = cpl_array_wrap_double (data, nwave);
                }

                /* Set FLUX  (wrap is the fastest to create an array) */
                for (int tel = 0; tel < ntel; tel++){
                    double * data = cpl_malloc (nwave * sizeof(double));
                    for (cpl_size wave = 0 ; wave < nwave ; wave++ ) 
                        data[wave] = pOut[tel*nwave+wave];
                    tFlux[row*ntel+tel] = cpl_array_wrap_double (data, nwave);
                }

                /* Set FLUXERR */
                for (int tel = 0; tel < ntel; tel++){
                    double * data = cpl_malloc (nwave * sizeof(double));
                    for (cpl_size wave = 0 ; wave < nwave ; wave++ ) 
                        data[wave] = sqrt (pOutVar[tel*nwave+wave]);
                    tFluxErr[row*ntel+tel] = cpl_array_wrap_double (data, nwave);
                }

                /* Set VISDATA */
                for (int base = 0; base < nbase; base++){
                    double complex * data = cpl_malloc (nwave * sizeof(double complex));
                    for (cpl_size wave = 0 ; wave < nwave ; wave++ ) 
                        data[wave] = (double complex)( pOut[(base+4)*nwave+wave] + 1.*I * pOut[(base+10)*nwave+wave] );
                    tVis[row*nbase+base] = cpl_array_wrap_double_complex (data, nwave);
                }

                /* Set VISDATAERR */
                for (int base = 0; base < nbase; base++){
                    double complex * data = cpl_malloc (nwave * sizeof(double complex));
                    for (cpl_size wave = 0 ; wave < nwave ; wave++ ) 
                        data[wave] = (double complex)( sqrt(pOutVar[(base+4)*nwave+wave]) + 1.*I * sqrt(pOutVar[(base+10)*nwave+wave]) );
                    tVisErr[row*nbase+base] = cpl_array_wrap_double_complex (data, nwave);
                }

                /* Free the PREPROC data to save memory
                 * (ex: max pointer 34065038 versus 48947400) */
                for (cpl_size reg = 0; reg < n_region; reg++) {
                    FREE (cpl_array_delete, pRegArr[reg][row]);
                    FREE (cpl_array_delete, pErrArr[reg][row]);
                }
                
 				
            }/* End loop on frames */

			/* Deallocation of variables */
			cpl_msg_debug (cpl_func, "Free pointers to data");
			FREE (cpl_free, pReg);
			FREE (cpl_free, pErr);
			FREE (cpl_free, pRegArr);
			FREE (cpl_free, pErrArr);
			FREE (cpl_free, pOut);
			FREE (cpl_free, pOutVar);
			FREE (cpl_free, pChi2);
			FREE (cpl_free, pP2VM);
			FREE (cpl_free, pV2PM);
			FREELOOP (cpl_matrix_delete, p2vm, nwave);
			FREELOOP (cpl_matrix_delete, v2pm, nwave);

			/* Check how "flux conservative" is the P2VM */
			cpl_msg_info (cpl_func, "Total flux in TELs: %.2f [e], in REGIONs:%.2f [e]  (ratio=%.5f)",
						  full_flux_tel,full_flux_reg,full_flux_tel/full_flux_reg);
			
			sprintf (qc_name, "ESO QC TRANS P2VM %s",GRAVI_TYPE(type_data));
			cpl_propertylist_update_double (p2vmred_header, qc_name, full_flux_tel/full_flux_reg);
			cpl_propertylist_set_comment (p2vmred_header, qc_name, "[e/e] at P2VM extraction");
            
			/* Count time */
			cpl_msg_info (cpl_func, "Total time to apply matrix: %.4f s", (double)(clock()-timer_start)/(double)CLOCKS_PER_SEC );
            
			if( type_data==GRAVI_FT )
			  timer_ft += (clock() - timer_start);
			else 
			  timer_sc += (clock() - timer_start);
			timer_start = clock();
		}
		/* End loop on polarisation */
	}
	/* Loop on data type SC / FT */

	/*
	 * Add QC for lambda met wavelength
	 */
	if (gravi_data_has_extension(preproc_data, GRAVI_METROLOGY_EXT)){
		double lambda_met = gravi_pfits_get_met_wavelength_mean(preproc_header, gravi_data_get_table(preproc_data, GRAVI_METROLOGY_EXT));
		sprintf (qc_name, "ESO QC MET LAMBDA MEAN");
		cpl_propertylist_update_double (p2vmred_header, qc_name, lambda_met);
		cpl_propertylist_set_comment (p2vmred_header, qc_name, "[m] Effective mean metrology wavelength");
	}


	/*
	 * Print timing 
	 */
	timer_other += (clock() - timer_start);
	cpl_msg_info (cpl_func, "Total time for FT: %10.4f s", (double)(timer_ft)/(double)CLOCKS_PER_SEC);
	cpl_msg_info (cpl_func, "Total time for SC: %10.4f s", (double)(timer_sc)/(double)CLOCKS_PER_SEC);
	cpl_msg_info (cpl_func, "Total time for OTHER: %7.4f s", (double)(timer_other)/(double)CLOCKS_PER_SEC);

	/* Message and timer */
	gravi_msg_function_exit(1);
	return p2vmred_data;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Compute the real-time tracking state from OPDC
 *
 * @param p2vmred_data     The input/output gravi_data
 *
 * Compute the FT tracking state per baseline (saved in the existing
 * OI_VIS extension) and telescope (saved in the existing OI_FLUX
 * table). It also computes the the target phase of the OPDC per
 * baeline (saved in the existing OI_VIS table).
 * Compute various QC parameters related to tracking quality.
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_compute_opdc_state (gravi_data * p2vmred_data)
{
	/* Message and timer */
	gravi_msg_function_start(1);
    cpl_ensure_code (p2vmred_data, CPL_ERROR_NULL_INPUT);
    
    int nbase = 6, ntel = 4;
    char qc_name[90];

    /* Get necessary data */
    cpl_propertylist * header = gravi_data_get_header (p2vmred_data);
    int npol_ft = gravi_pfits_get_pola_num (header, GRAVI_FT);
	cpl_table * oi_vis  = gravi_data_get_oi_vis (p2vmred_data, GRAVI_FT, 0, npol_ft);
	cpl_table * oi_flux = gravi_data_get_oi_flux (p2vmred_data, GRAVI_FT, 0, npol_ft);
	cpl_table * opdc   = gravi_data_get_table (p2vmred_data, GRAVI_OPDC_EXT);
	cpl_size nrow_ft   = cpl_table_get_nrow (oi_vis) / nbase;
	cpl_size nrow_opdc = cpl_table_get_nrow (opdc);
    CPLCHECK_MSG ("Cannot get data");

	/* Create the target phase for each baseline in the OI_VIS table */
	gravi_table_new_column (oi_vis, "TARGET_PHASE", "rad", CPL_TYPE_DOUBLE);
	
	/* Create the OPDC state for each telescope in the OI_FLUX table */
	cpl_table_new_column (oi_flux, "STATE", CPL_TYPE_INT);
	cpl_table_fill_column_window (oi_flux, "STATE", 0, nrow_ft * ntel, -1);
	
	/* Create the OPDC state for each baseline in the OI_VIS table */
	cpl_table_new_column (oi_vis,  "STATE", CPL_TYPE_INT);
	cpl_table_fill_column_window (oi_vis,  "STATE", 0, nrow_ft * nbase, -1);

    /* Create the Global OPDC state in the OI_VIS table */
	cpl_table_new_column (oi_vis, "OPDC_STATE", CPL_TYPE_INT);
	cpl_table_fill_column_window (oi_vis, "OPDC_STATE", 0, nrow_ft * nbase, -1);
    
	if (nrow_opdc < nrow_ft) 
	  cpl_msg_warning (cpl_func,"Missing FT or OPDC data:  nrow_ft - nrow_opdc = %lli", nrow_ft-nrow_opdc);


    /* Set the global OPDC state */
    int * time_opdc  = cpl_table_get_data_int (opdc, "TIME");
    double * time_ft = cpl_table_get_data_double (oi_vis, "TIME");
    int * global_state_opdc = cpl_table_get_data_int (opdc, "STATE");
    int * global_state = cpl_table_get_data_int (oi_vis, "OPDC_STATE");
    CPLCHECK_MSG ("Cannot get data");
    
    /* Loop on FT rows */
    for (cpl_size row_opdc=0, row_ft=0 ; row_ft<nrow_ft ; row_ft++) {

		/* Check bounds or find the OPDC sample just following the current FT 
		 * FIXME: We should use the closesd OPDC sample in the past, not future */
		if ( (time_ft[row_ft*nbase] < time_opdc[0]) || (time_ft[row_ft*nbase] > time_opdc[nrow_opdc-1]) ) continue;
		while ( time_ft[row_ft*nbase] > time_opdc[row_opdc] ) row_opdc ++;
        
        /* Set the global OPDC state */
        for (int base = 0; base < nbase; base++)
            global_state[row_ft*nbase+base] = global_state_opdc[row_opdc];
    }

    /* BASELINE_STATE was not in the original data of the instrument */
	if ( cpl_table_has_column (opdc,"BASELINE_STATE") ) {
	  
	  int * steps_opdc = cpl_table_get_data_int (opdc, "STEPS");
	  int * state_opdc = cpl_table_get_data_int (opdc, "BASELINE_STATE");
	  int * base_state = cpl_table_get_data_int (oi_vis, "STATE");
	  int * tel_state  = cpl_table_get_data_int (oi_flux, "STATE");
	  double * base_steps = cpl_table_get_data_double (oi_vis, "TARGET_PHASE");
	  CPLCHECK_MSG ("Cannot get data");
	  
	  /* Loop on FT rows */
	  for (cpl_size row_opdc=0, row_ft=0 ; row_ft<nrow_ft ; row_ft++) {

		/* Check bounds or find the OPDC sample just following the current FT 
		 * FIXME: We should use the closesd OPDC sample in the past, not future */
		if ( (time_ft[row_ft*nbase] < time_opdc[0]) || (time_ft[row_ft*nbase] > time_opdc[nrow_opdc-1]) ) continue;
		while ( time_ft[row_ft*nbase] > time_opdc[row_opdc] ) row_opdc ++;

        /* Get the flag */
        int state_flag = state_opdc[row_opdc];
		
		/* Disentangle the state of each telescope:
         * The 4 beams are the first 4 bits */
        for (int tel = 0; tel < ntel; tel++) 
            tel_state[row_ft*ntel+tel] = gravi_bit_get (state_flag, tel);

		/* Disentangle the state of each baseline:
         * The 6 baselines are the bits from 5 to 10 */
        for (int base = 0; base < nbase; base++)
            base_state[row_ft*nbase+base] = gravi_bit_get (state_flag, ntel+base);

        /* Bit 11 (>>10) is the Kalman state */
        // kalman = gravi_bit_get (state_opdc[row_opdc], ntel+nbase);

        /* Remaing bits are the off-load system */
        // offload = state_opdc[row_opdc] >> 11;
        
		/* Use the closing triangles to recover the tracking on some baselines */
		for (cpl_size base = 0; base < nbase; base++) {
		  base_state[row_ft*nbase+base] = CPL_MAX( base_state[row_ft*nbase+base],
                                                   base_state[row_ft*nbase+GRAVI_TRI_BASE[base][0][0]] *
                                                   base_state[row_ft*nbase+GRAVI_TRI_BASE[base][0][1]] );
		  base_state[row_ft*nbase+base] = CPL_MAX( base_state[row_ft*nbase+base],
                                                   base_state[row_ft*nbase+GRAVI_TRI_BASE[base][1][0]] *
                                                   base_state[row_ft*nbase+GRAVI_TRI_BASE[base][1][1]] );
		}

		/* Disentangle the piezo steps, in units of [pi/8 rad] 
         * Each beam is coded over 4 bits (16 values over the circle)
         * *BUT* the last telescope is coded in wrong memory place !
         * beam0 = 0-3, beam1 = 4-7, beam2 = 8-11, beam3 = 16-19 */
        int tmp_arr[4];
        for (int tel = 0; tel < ntel; tel++) {
            int pos = (tel<3) ? 4*tel : 4*(tel+1);
            tmp_arr[tel] = 15 & ((steps_opdc[row_opdc]) >> pos);
        }

		/* Compute the FT target phase of each baseline, in [rad] */
		for (int base=0; base<nbase; base++) {
		  base_steps[row_ft*nbase+base] = (tmp_arr[GRAVI_BASE_TEL[base][1]] - tmp_arr[GRAVI_BASE_TEL[base][0]]) * CPL_MATH_PI / 8.0;
		}
        
	  } /* End loop on FT rows */
      
	} else {
	  cpl_msg_warning (cpl_func,"No column BASELINE_STATE in OPDC... old data ?");
	  cpl_msg_warning (cpl_func,"The STATE flags are set to 1 (valid) although the information is unavailable");
      cpl_table_fill_column_window (oi_flux, "STATE", 0, nrow_ft * ntel, 1);
      cpl_table_fill_column_window (oi_vis,  "STATE", 0, nrow_ft * nbase, 1);
	}

	/* Duplicate in the second polarisation if any */
	if ( npol_ft>1 ) {
      cpl_table * tmp;
	  cpl_msg_debug (cpl_func, "Duplicate the FT tracking state for 2sd polarisation");
      
      tmp = gravi_data_get_oi_vis (p2vmred_data, GRAVI_FT, 1, npol_ft);
	  cpl_table_duplicate_column (tmp, "TARGET_PHASE", oi_vis, "TARGET_PHASE");
	  cpl_table_duplicate_column (tmp, "STATE", oi_vis, "STATE");
	  cpl_table_duplicate_column (tmp, "OPDC_STATE", oi_vis, "OPDC_STATE");
      
      tmp = gravi_data_get_oi_flux (p2vmred_data, GRAVI_FT, 1, npol_ft);
	  cpl_table_duplicate_column (tmp, "STATE", oi_flux, "STATE");
	}
	

	/* 
	 * QC: Compute the phase RMS per base when in tracking state 
	 */
	
	cpl_array **visdata   = cpl_table_get_data_array (oi_vis,"VISDATA");
	int * state           = cpl_table_get_data_int (oi_vis,"STATE");
	double * target_phase = cpl_table_get_data_double (oi_vis,"TARGET_PHASE");
	CPLCHECK_MSG ("Cannot get data");
	
	double complex * tmp_cpx = cpl_malloc (sizeof(double complex)*nrow_ft);

	/* Loop on base */
	for (int base = 0; base < nbase; base ++) {

	  /* Compute real-time phase and its mean (as phasors) */
	  double complex mean_cpx = 0.0 + I * 0.0;
	  for (cpl_size row_ft=0 ; row_ft<nrow_ft ; row_ft++) {
		if (state[row_ft*nbase+base] < 1) continue;
		tmp_cpx[row_ft] = cpl_array_get_mean_complex (visdata[row_ft*nbase+base]) * cexp (-I*target_phase[row_ft*nbase+base]);
		mean_cpx += tmp_cpx[row_ft];
	  }

	  /* Compute <phi**2> and the number of valid point */
	  double sum = 0.0000001, sum2 = 0.0;
	  for (cpl_size row_ft=0 ; row_ft<nrow_ft ; row_ft++) {
		if (state[row_ft*nbase+base] < 1) continue;
		sum2 += pow (carg (tmp_cpx[row_ft] * conj(mean_cpx)), 2);
		sum += 1.0;
	  }

	  /* Write the QC parameter for this base */
	  sprintf (qc_name, "ESO QC PHASE_FT%d%d RMS",GRAVI_BASE_TEL[base][0]+1, GRAVI_BASE_TEL[base][1]+1);
	  cpl_propertylist_update_double (header, qc_name, sqrt (sum2 / sum));
	  cpl_propertylist_set_comment (header, qc_name, "[rad] residuals when tracking");

	  /* Write the QC parameter for this base */
	  sprintf (qc_name, "ESO QC TRACKING_RATIO_FT%d%d",GRAVI_BASE_TEL[base][0]+1, GRAVI_BASE_TEL[base][1]+1);
	  cpl_propertylist_update_int (header, qc_name, sum*100.0/nrow_ft);
	  cpl_propertylist_set_comment (header, qc_name, "[%] ratio of time with tracking");
	}
	/* End loop on base */

	FREE (cpl_free,tmp_cpx);

	/* 
	 * QC: Compute fraction of time traking 
	 */
	
	int* tab_state = cpl_table_get_data_int (gravi_data_get_table (p2vmred_data, GRAVI_OPDC_EXT), "STATE");
	CPLCHECK_MSG ("Cannot get data");
	
	long tracking = 0;
	for (cpl_size row_opdc=0; row_opdc < nrow_opdc; row_opdc++) {
		if (tab_state[row_opdc] == 3 || tab_state[row_opdc] == 2) tracking ++;
	}

	sprintf (qc_name, "ESO QC TRACKING_RATIO");
	cpl_propertylist_update_int (header, qc_name, tracking*100/nrow_ft);
	cpl_propertylist_set_comment (header, qc_name, "[%] ratio of time with full FT tracking");

	CPLCHECK_MSG ("Cannot put QC parameters");

	/* Message and timer */
	gravi_msg_function_exit(1);
	return CPL_ERROR_NONE;
}


/*----------------------------------------------------------------------------*/
/**
 * @brief Compute the QC TAU0 parameter
 * 
 * @param data:   The input gravi_data, shall contain an OPDC table
 *
 * This function computes the QC TAU0 OPDC## parameters from the 
 * PIEZO signal stored in the OPDC table. The parameter are 
 * added to the main primary header of the data.
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_compute_tau0 (gravi_data * data)
{
  /* verbose */
  gravi_msg_function_start(1);
  cpl_ensure_code (data, CPL_ERROR_NULL_INPUT);

  int nbase = 6;
  double gain = 16.8; // [rad/V]
  char qc_name[90];
  
  /* Get the OPDC table */
  cpl_propertylist * header = gravi_data_get_header (data);
  cpl_table * opdc = gravi_data_get_table (data, GRAVI_OPDC_EXT);
  cpl_size nrow = cpl_table_get_nrow (opdc);

  CPLCHECK_MSG("Cannot load the OPDC table");

  /* Compute the delays to explore */
  cpl_size max_delay = 200;

  /* Time in microseconds [us] and offset in [V] */
  int *time    = cpl_table_get_data_int (opdc, "TIME");
  CPLCHECK_MSG("Cannot load the TIME column");
  
  float ** piezo = gravi_table_get_data_array_float (opdc, "PIEZO_DL_OFFSET");  
  CPLCHECK_MSG("Cannot load the PIEZO_OFFSET columns");

  /* Loop on base */
  for (int base=0; base<nbase; base++) {
	int t1 = GRAVI_BASE_TEL[base][0];
	int t2 = GRAVI_BASE_TEL[base][1];

	/* Init the QC parameters */
	sprintf (qc_name, "ESO QC TAU0 OPDC%d%d", t1+1,t2+1);
	cpl_propertylist_update_double (header, qc_name, GRAVI_NAN_DOUBLE);
	cpl_propertylist_set_comment (header, qc_name, "[s] tau0 for variance of 1 rad2");
	
	
	/* Loop on delays */
	for (cpl_size delay=0; delay < max_delay; delay++) {

	  /* srtf(d) = < (opd(t) - opd(t+d))^2> in [rad^2] */
	  double strf = 0.0;
	  for (cpl_size row=0; row < nrow - max_delay - 1 ; row++) {
		float diff = (piezo[row][t1]-piezo[row][t2]) - (piezo[row+delay][t1]-piezo[row+delay][t2]);
		strf += (double)(diff * diff);
	  }
	  strf = strf * gain * gain / (nrow - max_delay - 1);

	  /* If the variance is larger than 1rad2,
	   * we found the tau0 */
	  if ( strf > 1.0 ) {
		double tau0 = (time[delay+1] - time[1]) * 1e-6;
		cpl_msg_info (cpl_func,"Compute %s = %.2f [ms]", qc_name, 1e3 * tau0);
		cpl_propertylist_update_double (header, qc_name, tau0);
		break;
	  }
	  
	} /* End loop on delay */

  } /* End loop on base */

  /* Free the list of arrays */
  FREE (cpl_free, piezo);

  gravi_msg_function_exit(1);
  return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Compute the QC for the injection stability
 * 
 * @param data:   The input data, shall contain an OI_FLUX table
 *
 * For each telescope, the injection stability is calculated as the ratio
 * between the polarisation-averaged 9 percentile and 95 percentile
 * of the inejcted flux for each telescope.
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_compute_qc_injection (gravi_data * data)
{
  gravi_msg_function_start(1);
  cpl_ensure_code (data, CPL_ERROR_NULL_INPUT);
  
  char qc_name[100];
  double p05 = 0.0, p95 = 0.0;
  
  cpl_propertylist * header = gravi_data_get_header (data);
  int npol = gravi_pfits_get_pola_num (header, GRAVI_FT);
  int ntel = 4, hw = 2;

  /* Create the kernel to smooth the flux series */
  cpl_msg_info (cpl_func, "Smooth flux sery over %i samples", hw);

  /* How many rows do we have per telescope */
  cpl_table * table = gravi_data_get_oi_flux (data, GRAVI_FT, 0, npol);
  cpl_size nrow = cpl_table_get_nrow (table) / ntel;


  /* For each telescope */
  for (int tel = 0; tel<ntel; tel++) {

    /* Build a vector of total flux */
    cpl_vector * raw_flux = cpl_vector_new (nrow);
    cpl_vector_fill (raw_flux, 0.0);
    
    for (int pol = 0; pol<npol; pol++) {

        /* Get table */
        cpl_table * table = gravi_data_get_oi_flux (data, GRAVI_FT, pol, npol);
        const cpl_array ** flux_array = cpl_table_get_data_array_const (table, "FLUX");
        
        for (int row = 0; row<nrow; row++)
            cpl_vector_set (raw_flux, row, cpl_vector_get (raw_flux, row) +
                            cpl_array_get_mean (flux_array[tel+row*ntel]));
    }
    
    /* Convolve to filter noise, delete raw flux */
    cpl_vector * flux;
    flux = cpl_vector_filter_lowpass_create (raw_flux, CPL_LOWPASS_GAUSSIAN, hw);
    cpl_vector_delete (raw_flux);

    /* Define positive */
    for (int row = 0; row<nrow; row++) 
        cpl_vector_set (flux, row, CPL_MAX (0.0, cpl_vector_get (flux, row)));

    /* Sort the flux array and normalise */
    cpl_vector_sort (flux, CPL_SORT_ASCENDING);
    /* Do not divide if flux equal to 0*/
    if (cpl_vector_get_max (flux) > 0)
        cpl_vector_divide_scalar (flux, cpl_vector_get_max (flux) / 100);
        
    /* Histogram as a string, values are 0-100 */
    char qc_value[100];
    sprintf (qc_value,"%.0f,%.0f,%.0f,%.0f,%.0f,%.0f,%.0f,%.0f,%.0f",
             cpl_vector_get (flux, (cpl_size)(0.1*(nrow-1))),
             cpl_vector_get (flux, (cpl_size)(0.2*(nrow-1))),
             cpl_vector_get (flux, (cpl_size)(0.3*(nrow-1))),
             cpl_vector_get (flux, (cpl_size)(0.4*(nrow-1))),
             cpl_vector_get (flux, (cpl_size)(0.5*(nrow-1))),
             cpl_vector_get (flux, (cpl_size)(0.6*(nrow-1))),
             cpl_vector_get (flux, (cpl_size)(0.7*(nrow-1))),
             cpl_vector_get (flux, (cpl_size)(0.8*(nrow-1))),
             cpl_vector_get (flux, (cpl_size)(0.9*(nrow-1))));

    /* Create the QC entry in the FITS header */
    cpl_msg_info (cpl_func, "Flux Histo: %s", qc_value);
    sprintf (qc_name, "ESO QC FLUX_FT%d HISTO", tel+1);
    cpl_propertylist_update_string (header, qc_name, qc_value);
    cpl_propertylist_set_comment (header, qc_name, "decile in percent");
      
    /* Compute the 5 percentile and 95 percentile */
    p05 = cpl_vector_get (flux, (cpl_size)(0.05*(nrow-1)));
    p95 = cpl_vector_get (flux, (cpl_size)(0.95*(nrow-1)));
    /* remove error if p95 equal to 0*/
    if (p95<1) p95=1;

    /* Create the QC entry in the FITS header */
    sprintf (qc_name, "ESO QC FLUX_FT%d P05P95", tel+1);
    cpl_propertylist_update_double (header, qc_name, p05/p95);
    cpl_propertylist_set_comment (header, qc_name, "injected flux 5 percentile over 95 percentile");

    /* Deleve vector */
    cpl_vector_delete (flux);
  }

  gravi_msg_function_exit(1);
  return CPL_ERROR_NONE;
}


/*----------------------------------------------------------------------------*/
/**
 * @brief Compute the QC for the FT linearity
 *
 * @param data:   The input data, shall contain an OI_VIS and OPDC table
 *
 * For each baseline, the linearity is calculated between the OPD measured by the P2VM
 * (from the pipeline), and the OPD measured by the FT (in the OPDC table).
 * Biases as a function of OPD phase is measured and stored as QC parameters
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_compute_qc_ft_opd_estimator (gravi_data * p2vmred_data)
{
    gravi_msg_function_start(1);
    cpl_ensure_code (p2vmred_data, CPL_ERROR_NULL_INPUT);
    
    int nbase = 6;
    char qc_name[100];
    int nv;
    
    /* Get necessary data */
    cpl_propertylist * header = gravi_data_get_header (p2vmred_data);
    int npol = gravi_pfits_get_pola_num (header, GRAVI_FT);
    cpl_table * opdc   = gravi_data_get_table (p2vmred_data, GRAVI_OPDC_EXT);
    int * time_opdc  = cpl_table_get_data_int (opdc, "TIME");
    if (cpl_table_has_column(opdc, "OPD") == 1){
        cpl_array **opd_data_opdc   = cpl_table_get_data_array (opdc,"OPD");
        cpl_size nrow = cpl_table_get_nrow (opdc);
        CPLCHECK_MSG ("Cannot get data");
        
        /* initialize arrays */
        cpl_array * opd_res_array = cpl_array_new(nrow, CPL_TYPE_DOUBLE_COMPLEX);
        cpl_array * opd_array = cpl_array_new(nrow,  CPL_TYPE_DOUBLE);
        cpl_array * snr_array = cpl_array_new(nrow,  CPL_TYPE_DOUBLE);
        
        /* do it for both polarisations (if needed) */
        for (int pol = 0; pol<npol; pol++) {
            cpl_table * oi_vis  = gravi_data_get_oi_vis (p2vmred_data, GRAVI_FT, pol, npol);
            double * time_oivis  = cpl_table_get_data_double (oi_vis, "TIME");
            cpl_array **visdata  = cpl_table_get_data_array (oi_vis,"VISDATA");
            cpl_array **visdata_err   = cpl_table_get_data_array (oi_vis,"VISERR");

            cpl_size nrow_oivis  = cpl_table_get_nrow (oi_vis);
            cpl_size nwave= cpl_array_get_size (visdata[0]);


            /* For each baseline */
            for (int base = 0; base<nbase; base++) {
                cpl_size row_oivis = base;
                /* For each row */
                for (cpl_size row = 0; row<nrow; row++)
                {

                    /* FIXME: get rid of the 0.0011/2 (half step) */
                    while (( fabs (time_oivis[row_oivis+nbase]-time_opdc[row]+0.0011/2*1e6 )< fabs (time_oivis[row_oivis]-time_opdc[row]+0.0011/2*1e6 ) ) && (row_oivis<nrow_oivis-2*nbase))
                    {
                        row_oivis+=nbase;
                    }

                    double complex visdata_mean=cpl_array_get_mean_complex (visdata[row_oivis]);
                    cpl_array_set_double_complex(opd_res_array,row, visdata_mean*cexp ( I * cpl_array_get_float (opd_data_opdc[row], base, &nv) ));
                    cpl_array_set_double(opd_array, row, carg(visdata_mean));


                    /* get SNR */
                    double signal =0.0;
                    double noise  =0.0;
                    for (cpl_size wave = 0 ; wave < nwave ; wave ++)
                    {
                        signal+=cabs(cpl_array_get_double_complex (visdata[row_oivis],wave, &nv));
                        double err=cabs(cpl_array_get_double_complex (visdata_err[row_oivis],wave, &nv));
                        noise+=err*err;
                    }
                    noise=sqrt(noise+1);
                    if (signal > noise)
                        cpl_array_set_double(snr_array,row, signal /noise -1 );
                    else
                        cpl_array_set_double(snr_array,row, 0 );

                }
                CPLCHECK_MSG ("Cannot compute new array for cross-referencing OPD and FT OPD");
                
                double complex opd_mean = cpl_array_get_mean_complex (opd_res_array);
                cpl_array_multiply_scalar_complex(opd_res_array,conj(opd_mean));
                
                double c_o=0,s_o=0;
                double a_num=0,b_num=0,c_num=0,d_num=0,a_den=1,b_den=1,c_den=1,d_den=1;
                double opd_std_num=0,opd_std_den=1;
                double opd, opd_ref, snr_ref;
                
                for (cpl_size row = 0; row<nrow; row++)
                {
                    opd=cpl_array_get_double(opd_array,row,&nv);
                    opd_ref=carg(cpl_array_get_double_complex(opd_res_array,row,&nv));
                    snr_ref=cpl_array_get_double(snr_array,row,&nv);
                    c_o=cos(opd);
                    a_num+=snr_ref*opd_ref*c_o;
                    a_den+=snr_ref*c_o*c_o;
                    s_o=sin(opd);
                    b_num+=snr_ref*opd_ref*s_o;
                    b_den+=snr_ref*s_o*s_o;
                    c_o=c_o*c_o;
                    c_num+=snr_ref*opd_ref*c_o;
                    c_den+=snr_ref*c_o*c_o;
                    s_o=s_o*s_o;
                    d_num+=snr_ref*opd_ref*s_o;
                    d_den+=snr_ref*s_o*s_o;

                    /* calculate the standard deviation: numerator and denominator*/
                    opd_std_num+=snr_ref*opd_ref*snr_ref*opd_ref;
                    opd_std_den+=snr_ref;

                }
                
                CPLCHECK_MSG ("Cannot calculate variance of linearity");

                double a=a_num/a_den;
                double b=b_num/b_den;
                double c=c_num/c_den;
                double d=d_num/d_den;
                /* calculate standard deviation */
                double opd_std=sqrt(0.1+opd_std_num)/opd_std_den;

                /* calculate the average standard error of cos, sin, sin2 and cos2 */
                double res_car=a*a+b*b+c*c+d*d-16*opd_std*opd_std;
                if (res_car > 0) res_car=sqrt(res_car);
                else res_car=0;
            
                /* Create the QC entry in the FITS header*/
                if (npol == 2) sprintf (qc_name, "ESO QC LIN_FT P%d_B%d", pol,base+1);
                if (npol == 1) sprintf (qc_name, "ESO QC LIN_FT P%d_B%d", 3,base+1);
                cpl_propertylist_update_double (header, qc_name, res_car);
                cpl_propertylist_set_comment (header, qc_name, "FT nonlinearity biases [rad]");
                
                CPLCHECK_MSG ("Cannot store QCs");
                
            }
        }
        
        cpl_array_delete (opd_array);
        cpl_array_delete (opd_res_array);
        cpl_array_delete (snr_array);
        CPLCHECK_MSG ("Cannot close arrays");
    }
    else{
        cpl_msg_info(cpl_func, "Cannot found the OPD column (probably old data). Skip computation of QC LIN_FT");
    }
    
    gravi_msg_function_exit(1);
    return CPL_ERROR_NONE;
}


/**@}*/
