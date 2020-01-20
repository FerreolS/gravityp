/* $Id: gravi_preproc.c,v 1.12 2011/04/31 06:10:40 nazouaoui Exp $
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
 * @defgroup gravi_p2vm  P2VM calibration
 *
 * This modules implements the calibration of the P2VM. The main functions are called
 * by the @c gravity_p2vm recipe. After allocating the P2VM table with @c gravi_create_p2vm()
 * each file of the P2VM data set with one or two shutters opened (FLAT and P2VM) are processed with
 * @c gravi_compute_p2vm() to fill the TRANSMISSION, COHERENCE and PHASE columns of the P2VM
 * table. Then the TRANSMISSION are normalized with @c gravi_p2vm_normalisation() and
 * the file with 4 shutters opened (WAVE) is processed to calibrated the closure with
 * @c gravi_p2vm_phase_correction().
 *
 * The algorithms involved in this reduction are detailed in section Algorithms/Computation
 *  of the P2VM.
 */
/**@{*/

/*
 *  History
 *    04/12/2018 use GRAVITY_WAVE.fits calibration file instead of hardcoded values
 *    10/01/2019 fix a few warnings : roof_pos, qc_min, qc_max parameter usage
 *    14/03/2019 fix the selection for the medium wavelength range; remove unused DEFINE
 */
/*-----------------------------------------------------------------------------
                                   Includes
 -----------------------------------------------------------------------------*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cpl.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <complex.h>

#include "gravi_data.h"
#include "gravi_dfs.h"
#include "gravi_pfits.h"
#include "gravi_cpl.h"

#include "gravi_utils.h"
#include "gravi_ellipse.h"
#include "gravi_p2vm.h"


/*-----------------------------------------------------------------------------
                              Private prototypes
 -----------------------------------------------------------------------------*/

cpl_table* gravi_create_p2vm_table (cpl_table * detector_table,
                                    int nwave);

cpl_table * gravi_create_oiwave_table_sc (cpl_table * wave_table,
                                          cpl_propertylist * header,
                                          gravi_data *wave_param);

cpl_table * gravi_create_oiwave_table_ft (cpl_table * wave_table,
                                          cpl_table * detector_table,
                                          int pol);

/*-----------------------------------------------------------------------------
                              Function code
 -----------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/**
 * @brief  Create a new p2vm table
 * 
 * @param  n_region      The number of region
 * @param  nwave         The number of wavelength channels
 * @param  ntel          The number of telescopes
 * 
 * @return The p2vm table who must contain all of values of the p2vm
 * 
 * The function returns a table with regname, transmission,
 * coherence and phase columns. All columns are initialized 
 * with valid zero values.
 */
/*---------------------------------------------------------------------------*/

cpl_table* gravi_create_p2vm_table (cpl_table * detector_table,
                                    int nwave)
{
    gravi_msg_function_start(0);
    cpl_ensure (detector_table, CPL_ERROR_NULL_INPUT, NULL);
	cpl_ensure (nwave>0,   CPL_ERROR_ILLEGAL_INPUT, NULL);

    /* Get the number of region */
    int ntel = 4;
    int n_region = cpl_table_get_nrow (detector_table);

	/* Creating the table with nRegion rows */
	cpl_table * p2vm_table = cpl_table_new (n_region);

    /** REGNAME **/
    cpl_table_duplicate_column (p2vm_table, "REGNAME",
                                detector_table, "REGNAME");
	/** Transmission **/
	cpl_table_new_column_array (p2vm_table,"TRANSMISSION", CPL_TYPE_FLOAT,
                                nwave * ntel);
	/** Coherence **/
	cpl_table_new_column_array (p2vm_table,"COHERENCE", CPL_TYPE_FLOAT,
                                nwave * (ntel * (ntel - 1) / 2));
	/** Phase **/
	cpl_table_new_column_array (p2vm_table,"PHASE", CPL_TYPE_FLOAT,
                                nwave * (ntel * (ntel - 1) / 2));
	/** C_matrix **/
	cpl_table_new_column_array (p2vm_table,"C_MATRIX", CPL_TYPE_FLOAT,
                                nwave * (ntel * (ntel - 1) / 2));

	/* Initialization of all these columns  */
	cpl_array * zero_array_transmission, * zero_array;
	zero_array_transmission = cpl_array_new (nwave*ntel, CPL_TYPE_FLOAT);
	cpl_array_fill_window_float (zero_array_transmission, 0, nwave*ntel, 0);

	zero_array = cpl_array_new (nwave*(ntel*(ntel-1)/2),
                                CPL_TYPE_FLOAT );
	cpl_array_fill_window_float (zero_array, 0,
                                 nwave * (ntel * (ntel - 1) / 2), 0);

	/* Fill in the arrays on the columns */
	for (int i = 0; i < n_region; i++){
		cpl_table_set_array (p2vm_table, "TRANSMISSION", i, zero_array_transmission);
		cpl_table_set_array (p2vm_table, "COHERENCE", i, zero_array);
		cpl_table_set_array (p2vm_table, "PHASE", i, zero_array);
		cpl_table_set_array (p2vm_table, "C_MATRIX", i, zero_array);
	}

    cpl_array * dimensions;
	dimensions = cpl_array_new (2, CPL_TYPE_INT);

	/* Define the dimension of the transmission column */
	cpl_array_set_int (dimensions, 0, nwave);
	cpl_array_set_int (dimensions, 1, ntel);
	cpl_table_set_column_dimensions	(p2vm_table, "TRANSMISSION", dimensions);

	/* Define the dimension of the coherence column */
	cpl_array_set_int (dimensions, 1, (ntel*(ntel-1)/2));
	cpl_table_set_column_dimensions	(p2vm_table, "COHERENCE", dimensions);

	/* Define the dimension of the phase column */
	cpl_array_set_int (dimensions, 1, (ntel*(ntel-1)/2));
	cpl_table_set_column_dimensions	(p2vm_table, "PHASE", dimensions);

	cpl_table_set_column_dimensions	(p2vm_table, "C_MATRIX", dimensions);

	FREE (cpl_array_delete, zero_array_transmission);
	FREE (cpl_array_delete, zero_array);
	FREE (cpl_array_delete, dimensions);

    gravi_msg_function_exit(0);
	return p2vm_table;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief  Create a new oiwave table for SC
 * 
 * @param  wave_table    The input wavelength of each region in WAVE map
 * @param  header        The mean header of WAVE map
 * @param  params        The recipe parameters
 * 
 * @return The oiwave table
 * 
 * The function returns a table with EFF_WAVE, EFF_BAND.
 */
/*---------------------------------------------------------------------------*/

cpl_table * gravi_create_oiwave_table_sc (cpl_table * wave_table,
                                          cpl_propertylist * header,
                                          gravi_data *wave_param)
{
    gravi_msg_function_start(1);
    cpl_ensure (wave_table, CPL_ERROR_NULL_INPUT, NULL);
    cpl_ensure (header,     CPL_ERROR_NULL_INPUT, NULL);


    /* EKW START 04/12/2018 */
    /* set the calibrated eff_wave for LOW res*/
    /* OP 2018-02-12: new bandwidths measured for 3pix interpolation 
                      average of the six P1 baselines
		      also valid for P2
    */
/*
    double calib_eff_wave[14] = {7.046E-08,
                                 1.243E-07,
                                 1.342E-07,
                                 1.278E-07,
                                 1.255E-07,
                                 1.345E-07,
                                 1.423E-07,
                                 1.332E-07,
                                 1.235E-07,
                                 1.165E-07,
                                 1.158E-07,
                                 1.157E-07,
                                 1.120E-07,
                                 9.597E-08};

*/
    /* double *roof_pos; */
    double * calib_eff_wave;
    cpl_table * calib_eff_table = gravi_data_get_table (wave_param, "WAVE_TAB");
   // CPLCHECK_MSG ("STATIC_PARAM not available in the SOF. It is mandatory for acqcam reduction.");

     if ( cpl_table_has_column(calib_eff_table , "FBAND_WAVE") ) {
         calib_eff_wave = cpl_table_get_data_double (calib_eff_table, "FBAND_WAVE");
         cpl_msg_info(cpl_func,"calib_eff_wave [0] : %e", calib_eff_wave[0] );
         cpl_msg_info(cpl_func,"calib_eff_wave [1] : %e", calib_eff_wave[1] );
         cpl_msg_info(cpl_func,"calib_eff_wave [2] : %e", calib_eff_wave[2] );
         cpl_msg_info(cpl_func,"calib_eff_wave [3] : %e", calib_eff_wave[3] );
     }
     else {
       cpl_msg_error(cpl_func,"Cannot get the default values for calib_eff_wave");
     }

     /* Read the minimum and maximum wavelength values from the wave_param fits file
      * */
     double gravi_high_lbd_min = cpl_propertylist_get_double (gravi_data_get_plist(wave_param,GRAVI_PRIMARY_HDR_EXT), "ESO OIWAVE HIGH LBD MIN");
     cpl_msg_info (cpl_func,"gravi_high_lbd_min   : %e", gravi_high_lbd_min);

     double gravi_high_lbd_max = cpl_propertylist_get_double (gravi_data_get_plist(wave_param,GRAVI_PRIMARY_HDR_EXT), "ESO OIWAVE HIGH LBD MAX");
     cpl_msg_info (cpl_func,"gravi_high_lbd_max   : %e", gravi_high_lbd_max);

     double gravi_med_lbd_min = cpl_propertylist_get_double (gravi_data_get_plist(wave_param,GRAVI_PRIMARY_HDR_EXT), "ESO OIWAVE MED LBD MIN");
     cpl_msg_info (cpl_func,"gravi_med_lbd_min   : %e", gravi_med_lbd_min);

     double gravi_med_lbd_max = cpl_propertylist_get_double (gravi_data_get_plist(wave_param,GRAVI_PRIMARY_HDR_EXT), "ESO OIWAVE MED LBD MAX");
     cpl_msg_info (cpl_func,"gravi_med_lbd_max   : %e", gravi_med_lbd_max);

     double gravi_low_lbd_min = cpl_propertylist_get_double (gravi_data_get_plist(wave_param,GRAVI_PRIMARY_HDR_EXT), "ESO OIWAVE LOW LBD MIN");
     cpl_msg_info (cpl_func,"gravi_low_lbd_min   : %e", gravi_low_lbd_min);

     double gravi_low_lbd_max = cpl_propertylist_get_double (gravi_data_get_plist(wave_param,GRAVI_PRIMARY_HDR_EXT), "ESO OIWAVE LOW LBD MAX");
     cpl_msg_info (cpl_func,"gravi_low_lbd_max   : %e", gravi_low_lbd_max);

    /* EKW END 04/12/2018 */
    
    /* Get the QC */
    /*double qc_min, qc_max;
    qc_min = cpl_propertylist_get_double (header, QC_MINWAVE(GRAVI_SC));
    qc_max = cpl_propertylist_get_double (header, QC_MAXWAVE(GRAVI_SC));
    */
    CPLCHECK_NUL ("Cannot read the QC MINWAVE MAXWAVE");
    
    /* Get the max_wave and min_wave*/
    int n_element = cpl_table_get_column_depth (wave_table, "DATA1");
    CPLCHECK_NUL ("Cannot get the number of elements");
    double max_wave ,min_wave ;
    if (n_element<20)
    {
       /* EKW  max_wave = ESO OIWAVE LOW LBD MAX
        min_wave = ESO OIWAVE LOW LBD MIN; */
    	max_wave = gravi_low_lbd_max;
    	min_wave = gravi_low_lbd_min;
        cpl_msg_info (cpl_func,"Using Low resolution wavelength table");
        
    } else if (n_element<500)
    {
        /* EKW  max_wave = ESO OIWAVE MED LBD MAX;
         min_wave = ESO OIWAVE MED LBD MAX; */
        max_wave = gravi_med_lbd_max;
        min_wave = gravi_med_lbd_min;
        cpl_msg_info (cpl_func,"Using Med resolution wavelength table");
        
    } else
    {
        /*max_wave = ESO OIWAVE HIGH LBD MAX;
         min_wave = ESO OIWAVE HIGH LBD MAX; */
        max_wave = gravi_high_lbd_max;
        min_wave = gravi_high_lbd_min;
        cpl_msg_info (cpl_func,"Using High resolution wavelength table");
    }
    CPLCHECK_NUL ("Cannot get the max_wave and min_wave");
    
    /* Get the nwave */
    cpl_size nwave = gravi_wave_get_nlambda (wave_table, min_wave, max_wave);
    CPLCHECK_NUL ("Cannot define the limits");
    
    /* Create the OI_WAVELENGTH table */
    cpl_table * oiwave_table = cpl_table_new (nwave);
    
    cpl_table_new_column (oiwave_table, "EFF_WAVE", CPL_TYPE_FLOAT);
    cpl_table_set_column_unit (oiwave_table, "EFF_WAVE", "m");
    
    cpl_table_new_column (oiwave_table, "EFF_BAND", CPL_TYPE_FLOAT);
    cpl_table_set_column_unit (oiwave_table, "EFF_BAND", "m");
    
    /* Construction of the new wavelength table
     * To be improved if needed (non linear computation) */
    for (cpl_size wave = 0; wave < nwave; wave++){
        double lambda = min_wave + wave * (max_wave-min_wave)/(nwave-1);
        cpl_table_set (oiwave_table, "EFF_WAVE", wave, lambda);
        // introduce *2 to be closer to the real Band pass
        cpl_table_set (oiwave_table, "EFF_BAND", wave, (max_wave-min_wave)/(nwave-1)*2);
        //for LOW mode (nwave=11) set to measured bandpass
        if (nwave == 14) cpl_table_set (oiwave_table, "EFF_BAND", wave, calib_eff_wave[wave]);
    }

	//FREE(cpl_table_delete, calib_eff_table); /* EKW 04/12/2018 */

    //free(calib_eff_wave2); /* EKW 04/12/2018 */
    gravi_msg_function_exit(1);
    return oiwave_table;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief  Create a new oiwave table for FT
 * 
 * @param  wave_table        The input wavelength of each region in WAVE map
 * @param  detector_table    The corresponding detector table
 * @param  pol               The selected polarisation
 * 
 * @return The oiwave table
 * 
 * The function returns a table with EFF_WAVE, EFF_BAND.
 */
/*---------------------------------------------------------------------------*/

cpl_table * gravi_create_oiwave_table_ft (cpl_table * wave_table,
                                          cpl_table * detector_table,
                                          int pol)
{
    gravi_msg_function_start(1);
    cpl_ensure (wave_table,     CPL_ERROR_NULL_INPUT, NULL);
    cpl_ensure (detector_table, CPL_ERROR_NULL_INPUT, NULL);
    cpl_ensure (pol == 0 || pol == 1, CPL_ERROR_ILLEGAL_INPUT, NULL);
    
    /* Get the nwave */
    cpl_size nwave = gravi_spectrum_get_nwave (wave_table);
    cpl_size nreg = cpl_table_get_nrow (detector_table);

    /* Create the OI_WAVELENGTH table */
    cpl_table * oiwave_table = cpl_table_new (nwave);
    
    cpl_table_new_column (oiwave_table, "EFF_WAVE", CPL_TYPE_FLOAT);
    cpl_table_set_column_unit (oiwave_table, "EFF_WAVE", "m");
    
    cpl_table_new_column (oiwave_table, "EFF_BAND", CPL_TYPE_FLOAT);
    cpl_table_set_column_unit (oiwave_table, "EFF_BAND", "m");

    /* Fill each channel */
    for (cpl_size wave = 0; wave < nwave; wave++) {
        double lambda = (wave == nwave-1) ? 1.0 : 0.0;
        
        /* Loop on regions  */
        for (int reg = 0; reg < nreg; reg++) {
            if (gravi_region_get_pol (detector_table, reg) != pol) continue;
            double value  = gravi_table_get_value (wave_table, GRAVI_DATA[reg], 0, wave);
            
            if (wave == 0)             lambda = CPL_MAX (lambda, value + 0.00001e-6);
            else if (wave == nwave-1)  lambda = CPL_MIN (lambda, value - 0.00001e-6);
            else                       lambda += value / 24.0;
        }
        cpl_msg_info (cpl_func, "lbd[%lld] = %g [um]", wave, lambda*1e6);
        
        cpl_table_set (oiwave_table, "EFF_WAVE", wave, lambda);
        cpl_table_set (oiwave_table, "EFF_BAND", wave, 8.5e-8);
    }

    //cpl_msg_info (cpl_func, "Force linear wave");
    //double min_wave = cpl_table_get (oiwave_table, "EFF_WAVE", 0, NULL);
    //double max_wave = cpl_table_get (oiwave_table, "EFF_WAVE", nwave-1, NULL);
    //for (cpl_size wave = 0; wave < nwave; wave++){
    //    double lambda = min_wave + wave * (max_wave-min_wave)/(nwave-1);
    //    cpl_table_set (oiwave_table, "EFF_WAVE", wave, lambda);
    //}
    
    gravi_msg_function_exit(1);
    return oiwave_table;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Create a new P2VM map
 * 
 * @param preproc_data   one of the loaded data of this P2VM calibration
 * @param wave_map       the WAVE calibration map of this P2VM
 * 
 * @return a p2vm_map
 *
 * Create the P2VM map, including the P2VM_SC table, the P2VM_FT
 * table and a copy of the P2VM_MET table from wave_map. It also creates
 * the necessary OI_WAVELENGTHs table in the P2VM map.
 */
/*----------------------------------------------------------------------------*/

gravi_data * gravi_create_p2vm (gravi_data * wave_map, gravi_data *wave_param)
{
  cpl_table * oiwave_table;
  
  /* Verbose */
  gravi_msg_function_start(1);
  cpl_ensure (wave_map, CPL_ERROR_NULL_INPUT, NULL);
  
  cpl_propertylist * wave_header = gravi_data_get_header (wave_map);
  CPLCHECK_NUL ("Cannot load data");

  /* Create output map */
  gravi_data * p2vm_map = gravi_data_new(0);
  
  /* Duplicate the P2VM_MET in product */
  gravi_data_copy_ext (p2vm_map, wave_map, GRAVI_P2VM_MET_EXT);
  CPLCHECK_NUL ("Cannot copy tables");

  /* Loop on SC/FT */
  for (int type_data = 0; type_data < 2 ; type_data++) {

      /* Check if this data is present */
      if (!gravi_data_has_wave (wave_map, type_data)) {
          cpl_msg_info (cpl_func,"No data for %s, skip it", GRAVI_TYPE(type_data));
          continue;
      }

      /* Copy the IMAGING_DETECTOR table into product */
      gravi_data_copy_ext (p2vm_map, wave_map, GRAVI_IMAGING_DETECTOR_EXT(type_data));
      CPLCHECK_NUL ("Cannot copy tables");

      /* Get the WAVE_DATA and IMAGIGN_DETECTOR */
      cpl_table * wave_table = gravi_data_get_wave_data (wave_map, type_data);
      cpl_table * detector_table = gravi_data_get_imaging_detector (wave_map, type_data);
      CPLCHECK_NUL ("Cannot get data");
      
      /* 
       * Build the OI_WAVELENGTHs
       */
      
      /* Loop on polarisation */
      int n_pol  = gravi_spectrum_get_npol (wave_table);
      for (int pol = 0 ; pol<n_pol ; pol++) {

          if (type_data == GRAVI_SC)
              /* Create the SC OI_WAVELENGTH, will be the same for
               * both polarisations, and static from user requirement */
              oiwave_table = gravi_create_oiwave_table_sc (wave_table,
                                                           wave_header,
                                                           wave_param);
          else 
              /* Create the FT OI_WAVELENGTH, as a meach of each channel */
              oiwave_table = gravi_create_oiwave_table_ft (wave_table,
                                                           detector_table,
                                                           pol);
          
          /* Create plist */
          cpl_propertylist * oiwave_plist = cpl_propertylist_new ();
          cpl_propertylist_append_int (oiwave_plist, "NWAVE", cpl_table_get_nrow (oiwave_table));
          cpl_propertylist_update_int (oiwave_plist, "OI_REVN", 1);
          cpl_propertylist_append_string (oiwave_plist, "INSNAME", GRAVI_INSNAME(type_data,pol,n_pol));
          cpl_propertylist_append_int (oiwave_plist, "EXTVER", GRAVI_EXTVER(type_data,pol,n_pol));
          
          /* Add this OI_WAVELENGTH to the spectrum_data */
          gravi_data_add_table (p2vm_map, oiwave_plist,
                                GRAVI_OI_WAVELENGTH_EXT, oiwave_table);
          CPLCHECK_NUL ("Cannot add the table");
      }

      /* 
       * Build the P2VM
       */
      cpl_size n_wave = cpl_table_get_nrow (oiwave_table);

      /* Init the table */
      cpl_table * p2vm_table;
      p2vm_table = gravi_create_p2vm_table (detector_table, n_wave);
      
      /* Init the header of this table */
      cpl_propertylist * p2vm_plist = cpl_propertylist_new();
      CPLCHECK_NUL ("Cannot create propertylist");

      /* Add the p2vm_table */
      const char * extname = (type_data == GRAVI_SC ? GRAVI_P2VM_DATA_SC_EXT : GRAVI_P2VM_DATA_FT_EXT);
      gravi_data_add_table (p2vm_map, p2vm_plist, extname, p2vm_table);
      
      CPLCHECK_NUL ("Cannot create P2VM");
  } /* End loop SC/FT */
  
  gravi_msg_function_exit(1);
  return p2vm_map; 
}

/*----------------------------------------------------------------------------*/
/**
 * @brief The given output FITS file contain a p2vm table with the values
 *        of the transmission, phase and coherence extract using the p2vm matrix
 * 
 * @param p2vm_map		The p2vm table with the values of the transmission,
 * 	  	  	  	  		phase and coherence. The table must be allocated
 * 	  	  	  	  		before using this
 * @param preproc_data	The the raw data after preprocessing and will be
 * 	  	  	  	  	    identified to compute the p2vm
 * @param valid_trans	Integer array of 2 dimensions the first dimension
 * 	  	  	  	  	  	represents the type of the data (SC = 0 or FT = 1) ,
 * 	  	  	  	  	  	it means that it has 2 elements.
 * 	  	  	  	  	  	the second represents the information about the files
 * 	  	  	  	  	  	witch compute the transmission of witch telescope,
 * 	  	  	  	  	  	it means that it has 4 elements. (if
 * 	  	  	  	  	  	the file has a shutter 3 open of the data SC
 * 	  	  	  	  	  	valid_trans[0][2] = 1)
 * @param valid_pair	Integer array of 2 dimensions the first dimension
 * 	  	  	  	  	  	represents the type of the data (SC = 0 or FT = 1),
 * 	  	  	  	  	  	it means that it has 2 elements.
 * 	  	  	  	  	  	the second represents the information about the files
 * 	  	  	  	  	  	witch compute the phase and coherence of witch
 * 	  	  	  	  	  	base, it means that it has 4 elements.
 * 	  	  	  	  	  	(if the file has a shutter 3 and 4 open of the
 * 	  	  	  	  	  	data SC valid_pair[0][5] = 1)
 * @param det_type      The detector to align. If GRAVI_DET_SC, then only
 *                      the science detector extensions will be processed.
 *                      GRAVI_DET_FT will do the same for FT detector
 *                      and GRAVI_DET_ALL will do it for both. 
 * 
 * \exception CPL_ERROR_NULL_INPUT input data is missing
 * \exception CPL_ERROR_ILLEGAL_INPUT A table is missing in the input data,
 *  or is does not correspond with a file with one or 2 shutters opened.
 *
 * The function will compute the transmission, phase and coherence and save
 * them in the p2vm_map giving in the inputs (shall be initialized).
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_compute_p2vm (gravi_data * p2vm_map, gravi_data * preproc_data,
								   int ** valid_trans, int ** valid_pair,
                                   enum gravi_detector_type det_type)
{
    gravi_msg_function_start(1);
	cpl_ensure_code (p2vm_map,     CPL_ERROR_NULL_INPUT);
	cpl_ensure_code (preproc_data, CPL_ERROR_NULL_INPUT);

    int nv = 0;
	
	/* 
	 * Loop on type data  (SC/FT)
	 */
    int init_type_data = 1;
    int end_type_data = 1;
    if(det_type == GRAVI_DET_SC || det_type == GRAVI_DET_ALL)
        init_type_data = 0;
    if(det_type == GRAVI_DET_FT || det_type == GRAVI_DET_ALL)
        end_type_data = 1;
    for (int type_data = init_type_data; type_data <= end_type_data; type_data++ ) {

		/* Check if SPECTRUM data exists */
        if (!gravi_data_has_spectrum (preproc_data, type_data)) {
            cpl_msg_info (cpl_func,"No data for %s, skip it", GRAVI_TYPE(type_data));
            continue;
		}

		/* Primary header of input data */
		cpl_propertylist * preproc_header = gravi_data_get_header (preproc_data);
        
        /* Get the tables extension */
        cpl_table * detector_table, * spectrum_table, * p2vm_table;
        detector_table = gravi_data_get_imaging_detector (preproc_data, type_data);
        spectrum_table = gravi_data_get_spectrum_data (preproc_data, type_data);
        p2vm_table = gravi_data_get_p2vm_data (p2vm_map, type_data);
        
        cpl_ensure_code (detector_table, CPL_ERROR_ILLEGAL_INPUT);
        cpl_ensure_code (spectrum_table, CPL_ERROR_ILLEGAL_INPUT);
        cpl_ensure_code (p2vm_table,     CPL_ERROR_ILLEGAL_INPUT);
		
		/* Get the number of regions */
		int n_region = cpl_table_get_nrow (detector_table);
		int npol = n_region > 24 ? npol = 2 : 1;

        /* Get the oiwave_tables */
        cpl_table ** oiwave_tables = gravi_data_get_oiwave_tables (p2vm_map, type_data, npol);
        
        /* Get the number of waves */
		int nwave = cpl_table_get_nrow (oiwave_tables[0]);

		/* Get the P2VM column (where to store the data) */
        cpl_array ** transmission, ** coherence, ** phase, ** norma_m;
        coherence = cpl_table_get_data_array (p2vm_table, "COHERENCE");
        phase = cpl_table_get_data_array (p2vm_table, "PHASE");
        norma_m = cpl_table_get_data_array (p2vm_table, "C_MATRIX");
        transmission = cpl_table_get_data_array (p2vm_table, "TRANSMISSION");

        /* Get if a single telescope or baseline is open */
        int tel = -1, base = -1;
        if (gravi_check_shutter (preproc_header,1,0,0,0)) tel = 0;
        if (gravi_check_shutter (preproc_header,0,1,0,0)) tel = 1;
        if (gravi_check_shutter (preproc_header,0,0,1,0)) tel = 2;
        if (gravi_check_shutter (preproc_header,0,0,0,1)) tel = 3;
        if (gravi_check_shutter (preproc_header,1,1,0,0)) base = 0;
        if (gravi_check_shutter (preproc_header,1,0,1,0)) base = 1;
        if (gravi_check_shutter (preproc_header,1,0,0,1)) base = 2;
        if (gravi_check_shutter (preproc_header,0,1,1,0)) base = 3;
        if (gravi_check_shutter (preproc_header,0,1,0,1)) base = 4;
        if (gravi_check_shutter (preproc_header,0,0,1,1)) base = 5;
        cpl_ensure_code (tel>=0 || base>=0, CPL_ERROR_ILLEGAL_INPUT);

        /* 
         * A single telescope is open 
         */
        if (tel >= 0) {
            valid_trans[type_data][tel] = 1;

            cpl_msg_info (cpl_func, "Compute the transmission of "
                          "tel %d for %s", tel+1, GRAVI_TYPE(type_data));

            /* Init for integration */
            cpl_image * imgflux = NULL;
            cpl_image ** spectrum_img;
            spectrum_img = cpl_calloc (n_region,sizeof(cpl_image*));
            
            for (cpl_size region = 0; region < n_region; region++){
                
                /* Get the data of this region */
                cpl_imagelist * spectrum_imglist;
                spectrum_imglist = gravi_imagelist_wrap_column (spectrum_table, GRAVI_DATA[region]);
                
                /* Compute the median over the rows = un-normalized transmission */
                spectrum_img[region] = cpl_imagelist_collapse_median_create (spectrum_imglist);
                gravi_imagelist_unwrap_images (spectrum_imglist);
                
                /* Compute the total flux for this telescope (sum over all regions)
                 * for further normalization */
                if (imgflux == NULL)
                    imgflux = cpl_image_duplicate (spectrum_img[region]);
                else
                    cpl_image_add (imgflux, spectrum_img[region]);
            } /* loop on regions */
            
            /* Compute the total flux over the spectrum */
            double max_flux = cpl_image_get_flux (imgflux);
            
            /* Set the normalized transmission */
            for (cpl_size region = 0; region < n_region; region++){
                
                /* Normalisation of the transmission*/
                cpl_image_divide_scalar (spectrum_img[region], max_flux);

                /* Set the transmission in the 2D map of transmission */
                for (cpl_size wave = 0; wave < nwave; wave++) {
                    double value = cpl_image_get (spectrum_img[region], wave+1, 1, &nv);
                    cpl_array_set (transmission[region],
                                   wave + nwave * tel, value);
                }
                
                CPLCHECK_MSG ("Cannot compute transmission");
            }/* End loop on region */
            
            FREELOOP (cpl_image_delete, spectrum_img, n_region);
            FREE (cpl_image_delete, imgflux);
        } /* valid test for transmission (single tel) */
			

        /*
         * One pair is open
         */
        if (base >= 0) {
            int tel1 = GRAVI_BASE_TEL[base][0];
            int tel2 = GRAVI_BASE_TEL[base][1];
            
            valid_pair[type_data][base] = 1;
            
            cpl_msg_info(cpl_func, "Compute the coherence and "
                         "phase of pair (%d,%d) for %s", tel1+1,
                         tel2+1, GRAVI_TYPE(type_data));

            /* Get the size of the vectors */
            int nrow = cpl_table_get_nrow (spectrum_table);

            /* 
             * Recover the mean OPD modulation (averaved over pol
             * and channels) using the ellipses. In [m]
             */

            cpl_vector * mean_opd;
            mean_opd = gravi_ellipse_meanopd_create (spectrum_table, detector_table,
                                                     oiwave_tables, NULL, base);
            CPLCHECK_MSG ("Cannot compute opd");

            /* 
             * Compute the P2VM complex coherence for each channel and region
             */
            double residuals_avg = 0;
					
            for (cpl_size wave = 0; wave < nwave; wave++){

                /* Compute envelope from OPD for this channel 
                 * Wavelength is hard-coded in the enveloppe, thus
                 * doesn' depend on polarisation */
                cpl_vector * envelope_vector = gravi_compute_envelope (mean_opd, wave, nwave);

                /* Fill a matrix with nrows and 3 columns:
                 * [ 1.0, sin(2.pi.opd/lbd)*env(opd), cos(2.pi.opd/lbd)*env(opd) ] 
                 * different for each polarisation */
                cpl_matrix ** inv_matrixes = cpl_calloc (npol, sizeof(cpl_matrix *));
                cpl_matrix ** model_matrices = cpl_calloc (npol, sizeof(cpl_matrix *));
                
                for (int pol = 0; pol < npol; pol ++) {
                    model_matrices[pol] = cpl_matrix_new (nrow, 3);
                    
                    for (cpl_size row = 0; row < nrow; row++) {
                        cpl_matrix_set (model_matrices[pol], row, 0, 1.0);
                        double lambda = cpl_table_get (oiwave_tables[pol], "EFF_WAVE", wave, NULL);
                        double phi = cpl_vector_get (mean_opd, row) / lambda * CPL_MATH_2PI;
                        double enveloppe = cpl_vector_get (envelope_vector, row);
                        cpl_matrix_set (model_matrices[pol], row, 1, sin(phi)*enveloppe);
                        cpl_matrix_set (model_matrices[pol], row, 2, cos(phi)*enveloppe);
                    }

                    /* Invers the matrix of the carrying-wave */
                    inv_matrixes[pol] = gravi_matrix_invertSV_create (model_matrices[pol]);
                }

                /* Loop on region to apply this fit */
                for (cpl_size region = 0; region < n_region; region ++){

                    int pol = gravi_region_get_pol (detector_table, region);

                    cpl_vector * y_window;
                    y_window = gravi_table_get_vector (spectrum_table, wave,
                                                                           GRAVI_DATA[region]);
                    cpl_vector * yerr_window = gravi_table_get_vector (spectrum_table, wave,
                                                                           GRAVI_DATAERR[region]);

                    /* Vector init_val contains the best fit coeficient of the fit,
                     * that is the mean flux c and the complex coherence flux,
                     * computed as: init_val = inv_matrix * data 
                     * FIXME: can be done in CPL ? */
                    cpl_vector * init_val = cpl_vector_new(3);

                    for (cpl_size j = 0; j < 3; j++){
                        double comp = 0;
                        for (cpl_size i = 0; i < nrow; i++){
                            comp += cpl_matrix_get (inv_matrixes[pol], j, i) * cpl_vector_get (y_window, i);
                        }
                        cpl_vector_set (init_val, j, comp);
                    }
                    
                    /* compute the residuals */
                    cpl_vector * residuals = cpl_vector_new(nrow);
//                    cpl_vector * fit = cpl_vector_new(nrow);
                    for (cpl_size i = 0; i < nrow; i++){
                        double comp = 0;
                        for (cpl_size j = 0; j < 3; j++){
                            comp += cpl_matrix_get (model_matrices[pol], i, j) * cpl_vector_get (init_val, j);
                        }
                        cpl_vector_set (residuals, i, pow((cpl_vector_get (y_window, i)-comp)/cpl_vector_get (yerr_window, i) , 2));
//                        cpl_vector_set (fit, i, comp);
                    }

//                    if (region == 1)
//                    	if (wave == 7){
//                    		const cpl_vector ** vectors = malloc(4 * sizeof(cpl_vector*));
//                    		cpl_vector * error = cpl_vector_duplicate(y_window);
//                    		cpl_vector_subtract(error, fit);
//                    		vectors[0]=NULL;
//                    		vectors[1]=fit;
//                    		vectors[2] = error;
//                    		vectors[3] = y_window;
//                    		cpl_plot_vectors(NULL, NULL, NULL, vectors, 4);
//                    		cpl_plot_vector(NULL, NULL, NULL, error);
//                    		printf("Chi2 : %g \n", cpl_vector_get_mean(residuals));
//                    	}
                    residuals_avg += cpl_vector_get_mean(residuals);

                    /* Compute the P2VM coherence [e]. */
                    double coherence_fit =
                        sqrt( pow (cpl_vector_get(init_val, 2), 2) +
                              pow (cpl_vector_get(init_val, 1), 2));

                    cpl_array_set (coherence[region], wave + nwave *
                                   base, coherence_fit);

                    /* Compute the P2VM phase [rad] */
                    double phase_fit;
                    phase_fit = atan2( cpl_vector_get(init_val, 2),
                                       cpl_vector_get(init_val, 1));
                                        
                    cpl_array_set (phase[region], wave + nwave * base,
                                  phase_fit);

                    /* Save the c parameter for normalisation purpose */
                    cpl_array_set (norma_m[region], wave + nwave * base,
                                   cpl_vector_get (init_val,0));
                    
                    FREE (cpl_vector_delete, init_val);
                    FREE (cpl_vector_delete, y_window);
                    FREE (cpl_vector_delete, yerr_window);
                    FREE (cpl_vector_delete, residuals);
                } /* loop on region */

                FREELOOP (cpl_matrix_delete, inv_matrixes, npol);
                FREELOOP (cpl_matrix_delete, model_matrices, npol);
                FREE (cpl_vector_delete, envelope_vector);
            } /* loop on wave */

            /* Write QC in header */
        	cpl_propertylist * p2vm_header = gravi_data_get_header (p2vm_map);
    		char qc_name[90];
    		cpl_msg_info (cpl_func, "Averaged of CHI2 of p2vm fit for %s base %i = %f", GRAVI_TYPE(type_data), base+1, residuals_avg / (nwave*n_region));
    		sprintf (qc_name, "ESO QC P2VM%s%i CHI2",  GRAVI_TYPE(type_data), base+1);
    		cpl_propertylist_append_double (p2vm_header, qc_name, residuals_avg / (nwave*n_region));
    		cpl_propertylist_set_comment (p2vm_header, qc_name, "Chi2 avg. of p2vm fit");

            FREE (cpl_vector_delete, mean_opd);
        } /* End case valid_pair (baseline) */

        FREE (cpl_free, oiwave_tables);

		/* If an error is catched */
		CPLCHECK_MSG ("Cannot compute the P2VM");		

	} /* Loop on type_data SC/FT */

    gravi_msg_function_exit(1);
	return CPL_ERROR_NONE;
}


/*----------------------------------------------------------------------------*/
/**
 * @brief The given output FITS file contain a p2vm table with the values
 *        of the transmission, phase and coherence extract using the p2vm matrix so
 *        this function will normalise the p2vm map
 * 
 * @param p2vm_map      The given P2VM calibration already computed using
 *                      the function gravi_compute_p2vm
 * @param valid_trans	Integer array of 2 dimensions the first dimension
 * 	  	  	  	  	  	representes the type of the data (SC = 0 or FT = 1) ,
 * 	  	  	  	  	  	it means that it has 2 elements.
 * 	  	  	  	  	  	the second representes the information about the files
 * 	  	  	  	  	  	witch compute the transmission of witch telescope,
 * 	  	  	  	  	  	it means that it has 4 elements. (if
 * 	  	  	  	  	  	the file has a shutter 3 open of the data SC
 * 	  	  	  	  	  	valid_trans[0][2] = 1). This array will be useful
 * 	  	  	  	  	  	to check witch par of the p2vm map in the
 * 	  	  	  	  	  	transmission must be normalised
 * @param valid_pair	Integer array of 2 dimensions the first dimension
 * 	  	  	  	  	  	representes the type of the data (SC = 0 or FT = 1),
 * 	  	  	  	  	  	it means that it has 2 elements.
 * 	  	  	  	  	  	the second representes the information about the files
 * 	  	  	  	  	  	witch compute the phase and coherence of witch
 * 	  	  	  	  	  	base, it means that it has 4 elements.
 * 	  	  	  	  	  	(if the file has a shutter 3 and 4 open of the
 * 	  	  	  	  	  	data SC valid_pair[0][5] = 1)This array will be useful
 * 	  	  	  	  	  	to check witch par of the p2vm map in the phase
 * 	  	  	  	  	  	and coherence must be normalised
 * 
 * \exception CPL_ERROR_NULL_INPUT input data is missing
 * \exception CPL_ERROR_ILLEGAL_INPUT Cannot retrieve the number of wavelength
 *
 *
 * The function normalises the given p2vm map
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_p2vm_normalisation (gravi_data * p2vm_map,
                                         int ** valid_trans,
                                         int ** valid_pair )
{
	int ntel = 4, n_base = 6;
	int nv = 0;

	/* Message and timer */
	gravi_msg_function_start(1);
	cpl_ensure_code (p2vm_map,    CPL_ERROR_NULL_INPUT);
	cpl_ensure_code (valid_trans, CPL_ERROR_NULL_INPUT);
	cpl_ensure_code (valid_pair,  CPL_ERROR_NULL_INPUT);

	cpl_propertylist * p2vm_header = gravi_data_get_header (p2vm_map);

	/*
	 * Loop on FT and SC
	 */
	for (int type_data = 0; type_data < 2; type_data ++){

        /* Check if this data is present */
        if (!gravi_data_has_p2vm (p2vm_map, type_data)) {
            cpl_msg_info (cpl_func,"No data for %s, skip it", GRAVI_TYPE(type_data));
            continue;
        }

        cpl_msg_info (cpl_func, "Normalisation of P2VM for %s",
                      GRAVI_TYPE(type_data));
        
        /* Load P2VM and detector tables */
        cpl_table * p2vm_table, * detector_table;
        detector_table =  gravi_data_get_imaging_detector (p2vm_map, type_data);
        p2vm_table     =  gravi_data_get_p2vm_data (p2vm_map, type_data);

        /* Get info */
		cpl_size n_region = cpl_table_get_nrow (detector_table);
        int n_pol = n_region > 24 ? 2 : 1;

        /* Get the number of wave */
        cpl_table * oiwave_table;
        oiwave_table   =  gravi_data_get_oi_wave (p2vm_map, type_data, 0, n_pol);
		cpl_size nwave = cpl_table_get_nrow (oiwave_table);
        cpl_ensure_code (nwave, CPL_ERROR_ILLEGAL_INPUT);
        
		/* Get the P2VM column (where to store the data) */
        cpl_array ** transmission, ** coherence, ** phase;
        coherence = cpl_table_get_data_array (p2vm_table, "COHERENCE");
        phase = cpl_table_get_data_array (p2vm_table, "PHASE");
        transmission = cpl_table_get_data_array (p2vm_table, "TRANSMISSION");

		/*
		 * Set P2VM A-phase to zero
		 */
        
		/* Loop on polarisation  and baseline */
		for (int pol = 0; pol < n_pol; pol++) {
			for (int base = 0; base < 6; base++){
                
				/* Check if this baseline is computed */
				if ( !valid_pair[type_data][base]) continue;

                /* Get the index of the k band ABCD */
				int iA = gravi_get_region (detector_table, base, 'A', pol);
				int iB = gravi_get_region (detector_table, base, 'B', pol);
				int iC = gravi_get_region (detector_table, base, 'C', pol);
				int iD = gravi_get_region (detector_table, base, 'D', pol);
                
                for (cpl_size i = 0; i < cpl_array_get_size(phase[0]); i++){
                    
                    double valA = cpl_array_get (phase[iA], i, NULL);
                    cpl_array_set (phase[iA], i, 0.0);
                    
                    double valB = cpl_array_get (phase[iB], i, NULL) - valA;
                    if (valB < 0) valB += CPL_MATH_2PI;
                    cpl_array_set (phase[iB], i, valB);

                    double valC = cpl_array_get (phase[iC], i, NULL) - valA;
                    if (valC < 0) valC += CPL_MATH_2PI;
                    cpl_array_set (phase[iC], i, valC);
                    
                    double valD = cpl_array_get (phase[iD], i, NULL) - valA;
                    if (valD < 0) valD += CPL_MATH_2PI;
                    cpl_array_set (phase[iD], i, valD);
                    
                    CPLCHECK_MSG ("Cannot unwrap AC and BD of P2VM");		
                } /* End loop on wave */
			} /* end loop on baseline */
		} /* end loop on polarisation */

		CPLCHECK_MSG ("Cannot compute P2VM");		

		/*
		 * Normalisation of the P2VM transmission
		 */
        
		/* Compute the mean transmission */
		int ntrans = cpl_array_get_size (transmission[0]);
		cpl_array * mean_transmission = cpl_array_new (ntrans, CPL_TYPE_FLOAT);
		cpl_array_fill_window (mean_transmission, 0, ntrans, 0.);
        
		for (cpl_size region = 0; region < n_region; region++) {
			cpl_array_add (mean_transmission,transmission[region]);
		}
		cpl_array_divide_scalar (mean_transmission, n_region);

		/* Normalized the transmissions */
		for (cpl_size region = 0; region < n_region; region++){
			cpl_array_divide (transmission[region], mean_transmission);
		}
		FREE (cpl_array_delete, mean_transmission);

		/*
		 * Normalisation of the P2VM coherence
		 */

		/* Init quantities to compute mean flux per tel */
		double mean_flux[4] = {0.0,0.0,0.0,0.0};
		int n_flux[4] = {0,0,0,0};

        /* Loop on baseline */
        for (int base = 0; base < 6; base++) {
            int tel0 = GRAVI_BASE_TEL[base][0];
            int tel1 = GRAVI_BASE_TEL[base][1];
            
            /* Compute the input flux solving c=TI for each wave */
            for (cpl_size wave = 0; wave < nwave; wave++){
                
                /* Construction of the transmission matrix */
                cpl_matrix * matrix_T = cpl_matrix_new (n_region, 4);
                
                /* Get index of each tranmission on the matrix */
                cpl_vector * vector_c  = gravi_table_get_vector (p2vm_table, wave+nwave*base, "C_MATRIX");
                cpl_matrix * matrix_c = cpl_matrix_new (n_region, 1);
                for(cpl_size region = 0; region < n_region; region++) {
                    cpl_matrix_set (matrix_c, region, 0, cpl_vector_get (vector_c, region));
                    for (int tel = 0; tel < ntel; tel++){
                        cpl_matrix_set (matrix_T, region, tel,
                                        cpl_array_get (transmission[region],
                                                       tel*nwave+wave ,&nv));
                    }
                } /* End loop on region */
                FREE (cpl_vector_delete, vector_c);
                
                /* Compute the matrix_I */
                cpl_errorstate prestate = cpl_errorstate_get();
                
                cpl_matrix * matrix_I = cpl_matrix_solve_normal (matrix_T, matrix_c);
                
                /* Check if singular value */
                if (! strcmp("Singular matrix", cpl_error_get_message())){
                    cpl_msg_warning(cpl_func, "matrix_c or matrix_T "
                                    "are singular for tel1 = %d tel2 = %d  and wave = %lld",
                                    tel0, tel1, wave);
                    cpl_errorstate_set (prestate);
                    cpl_matrix_delete (matrix_T);
                    cpl_matrix_delete (matrix_I);
                    cpl_matrix_delete (matrix_c);
                } else {
                
                /* Store the total flux for each tel (will be used to compute
                 * the mean over wave, output and files) in [e/dit/output] */
                mean_flux[tel0] += cpl_matrix_get (matrix_I, tel0, 0);
                mean_flux[tel1] += cpl_matrix_get (matrix_I, tel1, 0);
                n_flux[tel0] ++;
                n_flux[tel1] ++;
				
                /* Get the coefficient of normalisation */
                double F0 = cpl_matrix_get (matrix_I, tel0, 0);
                double F1 = cpl_matrix_get (matrix_I, tel1, 0);
                double coeff = sqrt (fabs(F0*F1));
                
                /* Normalisation */
                if (coeff !=0 ) {
                    for (cpl_size region = 0; region < n_region; region++){
                        double value = cpl_array_get (coherence[region],
                                                      wave+nwave*base, &nv);
                        cpl_array_set (coherence[region],
                                       wave+nwave*base,
                                       value / coeff);
                    }
                }
                
                FREE (cpl_matrix_delete, matrix_c);
                FREE (cpl_matrix_delete, matrix_I);
                FREE (cpl_matrix_delete, matrix_T);
                }
            } /* end loop on WAVE */
        } /* end loop on baseline */



		/* 
		 * Compute the QC parameters for the mean flux 
		 */
		char qc_name[90];
		for (int tel = 0; tel<4; tel++ ) {
		  cpl_msg_info (cpl_func, "Mean flux tel%i = %f  (n=%i)", tel, mean_flux[tel] / n_flux[tel], n_flux[tel]);
		  sprintf (qc_name, "ESO QC FLUX_%s%i AVG",  GRAVI_TYPE(type_data), tel);
		  cpl_propertylist_append_double (p2vm_header, qc_name, mean_flux[tel] / n_flux[tel]);
		  cpl_propertylist_set_comment (p2vm_header, qc_name, "[e/DIT/chanel/output/file] mean flux");
		}


		/* Compute the QC parameters:
		 *  - Mean Coherence
		 *  - Mean of rms of Coherence
		 *  - Mean of rms of Phase */

		/* Allocate stat tables */
		cpl_array * sig_phi_arr = cpl_array_new(n_region, CPL_TYPE_DOUBLE);
		cpl_array * sig_coh_arr = cpl_array_new(n_region, CPL_TYPE_DOUBLE);
		cpl_array * mean_coh_arr = cpl_array_new(n_region, CPL_TYPE_DOUBLE);
		cpl_array * mean_coh_base_arr = cpl_array_new(n_base, CPL_TYPE_DOUBLE);
		cpl_array * min_coh_base_arr = cpl_array_new(n_base, CPL_TYPE_DOUBLE);
		cpl_array * mean_trans_arr = cpl_array_new(n_region, CPL_TYPE_DOUBLE);
		cpl_array * sig_transdiff_arr = cpl_array_new(n_region, CPL_TYPE_DOUBLE);

		cpl_array_fill_window (mean_coh_base_arr, 0, n_base, 0);
		cpl_array_fill_window (min_coh_base_arr, 0, n_base, 10);
        
		/* Compute the values per region */
		for (cpl_size region = 0; region < n_region; region ++ ) {
            
			/* Get the base of region */
			int base = gravi_region_get_base (detector_table, region);
            int tel0 = GRAVI_BASE_TEL[base][0];
            int tel1 = GRAVI_BASE_TEL[base][1];
            
            /* Get quantities for this region and base */
            cpl_array * coh_region = cpl_array_extract (coherence[region],
                                                        nwave*base, nwave);
            cpl_array * trans_tel1 = cpl_array_extract (transmission[region],
                                                        nwave*tel0, nwave);
            cpl_array * trans_tel2 = cpl_array_extract (transmission[region],
                                                        nwave*tel1, nwave);
            cpl_array * ph_region = cpl_array_extract (phase[region],
                                                       nwave*base, nwave);
            
            /* Compute spectral averaged quantities */
            cpl_array_set (mean_coh_arr, region,
                           cpl_array_get_mean (coh_region));
            
            cpl_array_set (sig_coh_arr, region,
                           cpl_array_get_stdev (coh_region));
            
            cpl_array_set (mean_trans_arr, region,
                           (cpl_array_get_mean (trans_tel1) +
                            cpl_array_get_mean (trans_tel2)) /2);
            
            cpl_array_set (sig_phi_arr, region,
                          cpl_array_get_stdev (ph_region));

            cpl_array * diff_trans;
            diff_trans = cpl_array_duplicate (trans_tel1);
            cpl_array_subtract (diff_trans, trans_tel2);
            cpl_array_set (sig_transdiff_arr, region,
                           cpl_array_get_stdev (diff_trans));
            cpl_array_delete (diff_trans);
            
            /* Normalize coherence with trans */
            cpl_array_power (trans_tel1, 0.5);
            cpl_array_power (trans_tel2, 0.5);
            cpl_array_divide (coh_region, trans_tel1);
            cpl_array_divide (coh_region, trans_tel2);
            cpl_array_divide_scalar (coh_region, 2.);
            
            /* Compute averaged quantities per baseline */
            cpl_array_set (mean_coh_base_arr, base,
                           cpl_array_get (mean_coh_base_arr, base, NULL) +
                           cpl_array_get_mean (coh_region)/(n_region/6));
            
            /* Compute percentil */
            double min_percentile = gravi_array_get_quantile (coh_region, 0.05);
            
            if (cpl_array_get(min_coh_base_arr, base, NULL) > min_percentile)
                cpl_array_set(min_coh_base_arr, base, min_percentile);
            
            /* Delete variables */
            cpl_array_delete (coh_region);
            cpl_array_delete (trans_tel2);
            cpl_array_delete (trans_tel1);
            cpl_array_delete (ph_region);
            
		} /* End loop on regions */
        
		CPLCHECK_MSG ("Cannot compute_the averages values per region");

		/* Compute the full-averaged QC */
		double mean_coh, sig_coh, sig_phi;
		mean_coh = cpl_array_get_mean (mean_coh_arr) /
                   (2*cpl_array_get_mean (mean_trans_arr));
		sig_coh  = cpl_array_get_mean (sig_coh_arr) /
                   (2*cpl_array_get_mean (mean_trans_arr));
		sig_phi  = cpl_array_get_mean (sig_phi_arr);
		CPLCHECK_MSG ("Cannot compute averaged values");

		/* Set these full-averaged QC parameters */
		cpl_propertylist_update_double (p2vm_header, (type_data)?QC_MEANCOH_FT:QC_MEANCOH_SC,
                                       mean_coh);
		cpl_propertylist_update_double (p2vm_header, (type_data)?QC_RMSCOH_FT:QC_RMSCOH_SC,
                                       sig_coh);
		cpl_propertylist_update_double (p2vm_header, (type_data)?QC_RMSPHASE_FT:QC_RMSPHASE_SC,
                                        sig_phi);
		cpl_msg_info (cpl_func, "QC %s COH_AVG = %e COH_RMS %e PHASE_RMS = %e",
                      GRAVI_TYPE(type_data), mean_coh, (sig_coh), sig_phi);
		
		/* QC per baseline */
		for (int base=0; base< n_base; base++){
			sprintf (qc_name, "ESO QC P2VM_COHERENCE_%s%s", GRAVI_TYPE(type_data), GRAVI_BASE_NAME[base]);
			cpl_propertylist_update_double (p2vm_header, qc_name,
                                            cpl_array_get (mean_coh_base_arr, base, NULL));
			cpl_propertylist_set_comment (p2vm_header, qc_name,
                                          "Avg coh. over lbd per baseline");
            
			sprintf (qc_name, "ESO QC P2VM_MINCOHERENCE_%s%s", GRAVI_TYPE(type_data), GRAVI_BASE_NAME[base]);
			cpl_propertylist_update_double (p2vm_header, qc_name,
                                            cpl_array_get (min_coh_base_arr, base, NULL));
			cpl_propertylist_set_comment (p2vm_header, qc_name,
                                          "Min coh. (5 perc) over lbd per baseline");
		}
		CPLCHECK_MSG ("Cannot append QC params");
		
		FREE (cpl_array_delete, sig_phi_arr);
		FREE (cpl_array_delete, sig_coh_arr);
		FREE (cpl_array_delete, mean_coh_arr);
		FREE (cpl_array_delete, mean_trans_arr);
		FREE (cpl_array_delete, sig_transdiff_arr);
		FREE (cpl_array_delete, mean_coh_base_arr);
		FREE (cpl_array_delete, min_coh_base_arr);

	} /* end loop on FT and SC */

	gravi_msg_function_exit(1);
	return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Correct the phase of the P2VM from internal closure-phases
 * 
 * @param p2vm_map     The given P2VM already computed using the function
 * 	  	  	  	  	   gravi_compute_p2vm
 * 
 * @param p2vmred_data The P2VMRED of the 4-shutter open file
 * 
 * @param full_phase 
 *  0: Force phiA(lbd) to be zero for baselines (01,02,03) = keep only closures
 *  1: Force phiA(lbd) to have zero mean and minimum GD for baselines (01,02,03)
 *  2: Force phiA(lbd) to have zero-GD for baselines (01,02,03)
 *
 * \exception CPL_ERROR_NULL_INPUT input data is missing
 * \exception CPL_ERROR_ILLEGAL_INPUT full_phase parameter out of range
 *
 * The p2vmreduced of a 4-shutter open observation
 * (generally the WAVE) can be used to compute the p2vm closure phase
 * and chromatic phase, which are then removed from the P2VM phases.
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_p2vm_phase_correction (gravi_data * p2vm_map,
                                            gravi_data * p2vmred_data,
                                            int full_phase)
{
	gravi_msg_function_start(1);
	cpl_ensure_code (p2vm_map,   CPL_ERROR_NULL_INPUT);
	cpl_ensure_code (p2vmred_data, CPL_ERROR_NULL_INPUT);
	
	int nv;

	/* Get the header */
	cpl_propertylist * p2vmred_header = gravi_data_get_header (p2vmred_data);
	
	/*
	 * Loop on FT and SC
	 */
	for (int type_data = 0; type_data < 2; type_data ++) {

		/* Check if P2VM has this type */
        if (!gravi_data_has_p2vm (p2vm_map, type_data)) {
            cpl_msg_info (cpl_func,"No data for %s in P2VM, skip it", GRAVI_TYPE(type_data));
            continue;
		}
	  
	  cpl_msg_info (cpl_func, "Calibrate the internal phase of %s (%s)",
			GRAVI_TYPE(type_data),full_phase?"full phases":"closure phases");

	  /* Where the instrumental phase will be stored */
	  cpl_array ** visphase = cpl_calloc (12, sizeof(cpl_array*));
	  
	  /* Loop on polarisation */
	  int npol = gravi_pfits_get_pola_num (p2vmred_header, type_data);
	  for ( int pol= 0 ; pol < npol ; pol++ ) {

		cpl_msg_info (cpl_func, "Compute correction for pol %i over %i", pol+1, npol);
		
		/* Get the OI_WAVE and the OI_VIS tables */
		cpl_table * oi_wave = gravi_data_get_oi_wave (p2vmred_data, type_data, pol, npol);
		cpl_table * oi_vis =  gravi_data_get_oi_vis (p2vmred_data, type_data, pol, npol);
		cpl_array * sigma = gravi_table_create_sigma_array (oi_wave);
		cpl_size nrow  = cpl_table_get_nrow (oi_vis) / 6;

		/* Get a pointer of the VISDATA column */
		cpl_array ** visdata = cpl_table_get_data_array (oi_vis, "VISDATA");
		float * wavedata = cpl_table_get_data_float (oi_wave, "EFF_WAVE");
		cpl_size nwave = cpl_array_get_size (visdata[0]);
		cpl_size wave0 = nwave/2;

	    /* 
		 * Compute OPLs of the 4 beam as the phase of the mean channel 
		 * for the baselines 01, 02, 03, in [rad]
		 */
		cpl_msg_debug (cpl_func,"Compute OPLs of beam");

		double ** opl  = cpl_malloc (sizeof(double*) * nrow);
		for (cpl_size row=0 ; row<nrow ; row++) {
		  opl[row] = cpl_malloc (sizeof(double) * 4);
		  opl[row][0] = 0.0;
		  for ( int base= 1 ; base < 4 ; base++ ) {
			opl[row][base] = carg (cpl_array_get_complex (visdata[row * 6 + base-1], wave0, &nv));
		  }
		} /* End compute OPLs */

		/* Unwrap OPLs */
		for ( int base= 1 ; base < 4 ; base++ ) {
		  double wrap = 0.0, ref = opl[0][base];
		  for (cpl_size row=1 ; row<nrow ; row++) {
			if ( (opl[row][base] - ref) < -CPL_MATH_PI ) wrap += 2.*CPL_MATH_PI;
			if ( (opl[row][base] - ref) >  CPL_MATH_PI ) wrap -= 2.*CPL_MATH_PI;
			ref = opl[row][base];
			opl[row][base] += wrap;
		  }
		} /* End unwrap OPLs */

		/* Remove mean OPLs */
		for ( int base= 1 ; base < 4 ; base++ ) {
		  double mean = 0.0;
		  for (cpl_size row=0 ; row<nrow ; row++) mean += opl[row][base];
		  for (cpl_size row=0 ; row<nrow ; row++) opl[row][base] -= mean / nrow;
		}

		CPLCHECK_MSG("Cannot compute OPL");
		
		/*
		 * Coherent integration of visdata of each baseline with the current P2VM
		 */
		cpl_msg_debug (cpl_func,"Compute coherent integration of VISDATA");

		for ( int base = 0 ; base < 6 ; base++ ) {
		  /* Allocate memory for visphase */
		  visphase[base + pol*6] = cpl_array_new (nwave, CPL_TYPE_FLOAT_COMPLEX);
		  cpl_array_fill_window_complex (visphase[base + pol*6], 0, nwave, 0.0 * I + 0.0);

		  /* Coherent integration of the visdata */
		  for (cpl_size row=0 ; row<nrow ; row++) {
			double x = opl[row][GRAVI_BASE_TEL[base][1]] - opl[row][GRAVI_BASE_TEL[base][0]];
			for (cpl_size wave=0; wave<nwave; wave++)
			  cpl_array_set_complex (visphase[base + pol*6], wave,
									 cpl_array_get_complex (visphase[base + pol*6], wave, &nv) +
									 cpl_array_get_complex (visdata[row * 6 + base], wave, &nv) *
									 cexp (-1.*I * x * wavedata[wave0] / wavedata[wave]) );
		  } /* End loop on rows and waves */

		  /* Normalize to pure phasor, to give equal weight to all channels */
		  // gravi_array_normalize_complex (visphase[base + pol*6]);
		  
		} /* End loop on base */


		/*
		 * Force the phase correction of the 3 first baselines to be zero (01,02,03).
		 * Update other baselines to keep the measured closure phases. Thus
		 * we only correct chromatic closure phases, not chromatic phase.
		 * Here we *assume* the 3 first baseline are {01, 02, 03}
		 */

		cpl_array ** ref = cpl_malloc (4 * sizeof(cpl_array*));
		ref[0] = cpl_array_new (nwave, CPL_TYPE_DOUBLE_COMPLEX);
		cpl_array_fill_window_complex (ref[0], 0, nwave, 0.0 * I + 1.0);

		if (full_phase == 0) {
		  
		  cpl_msg_info (cpl_func,"Force phiA(lbd) to be zero for baselines (01,02,03) = keep only closures");
		  
		  for (int base = 0; base<3 ; base++)
			ref[base+1] = cpl_array_duplicate (visphase[base + pol*6]);
		  
		} else if (full_phase == 1) {
		  
		  cpl_msg_info (cpl_func,"Force phiA(lbd) to have zero mean and minimum GD for baselines (01,02,03)");
		  
		  for (int base = 0; base<3 ; base++) {
              double gd = 0.0;
              gravi_array_get_group_delay_loop (&visphase[base + pol*6], NULL, sigma,
                                                &gd, 1, 2e-3, CPL_FALSE);
			cpl_array * tmp = cpl_array_duplicate (visphase[base + pol*6]);
			gravi_array_multiply_phasor (tmp, -CPL_MATH_2PI*I*gd, sigma);
			gd += carg (cpl_array_get_mean_complex (tmp)) / CPL_MATH_2PI / cpl_array_get_mean (sigma);
			cpl_array_delete (tmp);
			ref[base+1] = gravi_array_cexp (CPL_MATH_2PI * I * gd, sigma);
		  }
		  
		} else if (full_phase == 2) {
			
		  cpl_msg_info (cpl_func,"Force phiA(lbd) to have zero-GD for baselines (01,02,03)");
		  
		  for (int base = 0; base<3 ; base++) {
              double gd = 0.0;
              gravi_array_get_group_delay_loop (&visphase[base + pol*6], NULL, sigma,
                                                &gd, 1, 2e-3, CPL_FALSE);
			ref[base+1] = gravi_array_cexp (CPL_MATH_2PI * I * gd, sigma);
		  }
		} else {
		  return cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT, "Option for phase calibration out of range");
		}

		/* Remove the corresponding REF in each base 
		 * visphase_ij *= ref_i * conj(ref_j) */
		for ( int base = 0 ; base < 6 ; base++ ) {
		  cpl_array_multiply (visphase[base + pol*6], ref[GRAVI_BASE_TEL[base][0]]);
		  gravi_array_multiply_conj (visphase[base + pol*6], ref[GRAVI_BASE_TEL[base][1]]);
		}
		
          
          /* SL: 7 Sept 2019
           * Filter visphase to remove "wiggles", ie, high frequency oscilations
           */
          if (nwave>100) /* Only for MED and HIGH with SCIENCE SPECTROMETER */
              for ( int base = 0 ; base < 6 ; base++ ) {
                  
                  cpl_msg_info (cpl_func,"Performing wiggle removal on baseline %i (pol %i over %i)",
                                base+1, pol+1, npol);
                  
                  /* initialise polynomial of order 6 */
                  cpl_size mindeg = 0, maxdeg = 6;
                  cpl_polynomial * fit = cpl_polynomial_new (1);
                  
                  /* visphase will be added to the phase of the P2VM
                   * Need to be unwrapped before interpolation, and removal of the wiggles */
                  cpl_array * phase_unwraped = cpl_array_cast (visphase[base + pol*6], CPL_TYPE_DOUBLE_COMPLEX);
                  
                  /* Compute group-delay */
                  double gd = 0.0;
                  gravi_array_get_group_delay_loop (&phase_unwraped, NULL, sigma,
                                                    &gd, 1, 2e-3, CPL_FALSE);
                  
                  /* Remove mean group-delay and phase-delay to unwrap */
                  gravi_array_multiply_phasor (phase_unwraped, - 2*I*CPL_MATH_PI * gd, sigma);
                  double mean_phase = carg (cpl_array_get_mean_complex (phase_unwraped));
                  cpl_array_multiply_scalar_complex (phase_unwraped, cexp(- I * mean_phase));
                  
                  /* Compute argument and add back the delay and the phase [rad] */
                  cpl_array_arg (phase_unwraped);
                  gravi_array_add_phase (phase_unwraped, 2.*CPL_MATH_PI * gd, sigma);
                  cpl_array_add_scalar (phase_unwraped, mean_phase);
                  
                  /* Create the wavenumber matrix */
                  
                  double sigma0 = cpl_array_get_mean (sigma);
                  double delta0 = cpl_array_get_max (sigma)-cpl_array_get_min (sigma);
                  
                  cpl_matrix * sigma_matrix = cpl_matrix_new (1,nwave);
                  for (cpl_size wave = 0; wave < nwave; wave ++) {
                      cpl_matrix_set (sigma_matrix, 0, wave, (cpl_array_get (sigma,wave,&nv) - sigma0)/delta0 );
                  }
                  
                  /* Polynomial fit */
                  cpl_vector * input = cpl_vector_wrap (nwave, cpl_array_get_data_double (phase_unwraped));
                  cpl_polynomial_fit (fit, sigma_matrix, NULL, input, NULL, CPL_FALSE, &mindeg, &maxdeg);
                  cpl_vector_unwrap (input);
                  cpl_array_delete (phase_unwraped);
                  
                  /* Evaluate phase of visPhase */
                  for (cpl_size wave=0; wave<nwave;wave++) {
                      cpl_array_set_complex (visphase[base + pol*6], wave, cexp( I * cpl_polynomial_eval_1d (fit, cpl_matrix_get (sigma_matrix,0,wave), NULL)));
                  }

                  FREE (cpl_polynomial_delete, fit);
                  FREE (cpl_matrix_delete, sigma_matrix);
              }
          
		FREELOOP (cpl_array_delete, ref, 4);
		FREELOOP (cpl_free, opl, nrow);
		FREE (cpl_array_delete, sigma);
	  } /* End loop on pol */

	  CPLCHECK_MSG("Cannot perform coherent integration");
	  
        
            
        
        
	  /*
	   * Apply these phase corrections to the P2VM phases
	   */
	  cpl_msg_info (cpl_func,"Apply correction to P2VM phases");
	  
	  cpl_table * p2vm = gravi_data_get_p2vm_data (p2vm_map, type_data);
	  cpl_array ** phase = cpl_table_get_data_array (p2vm, "PHASE");
	  cpl_size nreg = cpl_table_get_nrow (p2vm);
	  cpl_size nwave = cpl_array_get_size (visphase[0]);

	  /* Loop on region */
	  for (cpl_size reg=0; reg<nreg; reg++) {
          
		/* Get the polarisation and base of this region */
		int base = gravi_region_get_base (p2vm, reg);
		int pol  = gravi_region_get_pol (p2vm, reg);

		/* Subtract the corresponding phase */
		for (cpl_size wave=0; wave<nwave;wave++) {
		  double phi = carg ( cexp (1.*I*cpl_array_get (phase[reg], wave + nwave * base, &nv)) *
							  conj(cpl_array_get_complex (visphase[base + pol*6], wave, &nv) ));
		  cpl_array_set (phase[reg], wave + nwave * base, (phi>=0.0?phi:phi+CPL_MATH_2PI));
		}

	  } /* End loop on region */

	  CPLCHECK_MSG("Cannot correct P2VM phases");
	  FREELOOP (cpl_array_delete, visphase, 12);
	  
	} /* End loop on FT and SC */
  
	gravi_msg_function_exit(1);
	return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Compute the flux normalisation in the P2VM
 * 
 * @param p2vm_map      The given p2vm map already computed using the
 *                      function gravi_compute_p2vm
 * 
 * @param p2vmred_data  The P2VMREDUCED of the 4-shutter open file
 *
 * \exception CPL_ERROR_NULL_INPUT input data is missing
 *
 * Compute the mean OI_FLUX from the p2vmreduced data.
 * Normalize to mean=1.0 and store it into the p2vm_map.
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_p2vm_transmission (gravi_data * p2vm_map, gravi_data * p2vmred_data)
{
	gravi_msg_function_start(1);
	cpl_ensure_code (p2vm_map,     CPL_ERROR_NULL_INPUT);
	cpl_ensure_code (p2vmred_data, CPL_ERROR_NULL_INPUT);
	
	cpl_propertylist * p2vmred_header = gravi_data_get_header (p2vmred_data);

	/* Loop on SC / FT */
	for (int type_data = 0; type_data < 2; type_data ++) {

	  /* Loop on polarisation */
	  int npol = gravi_pfits_get_pola_num (p2vmred_header, type_data);
	  for (int pol = 0; pol < npol; pol++) {

		cpl_msg_info (cpl_func, "Compute the internal transmission of %s (pol %i over %i)",
					  GRAVI_TYPE(type_data), pol+1, npol);
		
		/* Get FLUX data */
		cpl_table * flux_tbl = gravi_data_get_oi_flux (p2vmred_data, type_data, pol, npol);
		cpl_array ** flux = cpl_table_get_data_array (flux_tbl, "FLUX");
		cpl_size nwave = cpl_array_get_size (flux[0]);
		cpl_size nrow  = cpl_table_get_nrow (flux_tbl) / 4;

		/* Create OI_FLUX table */
		cpl_table * oi_flux = gravi_table_oi_create (nwave, 1, GRAVI_OI_FLUX_EXT);
		cpl_array * flux_mean = cpl_array_new (nwave, CPL_TYPE_DOUBLE);

		CPLCHECK_MSG("Cannot get data");
		
		for (int tel = 0; tel < 4; tel++) {

		  /* Compute total flux for this tel */
		  cpl_array_fill_window (flux_mean, 0, nwave, 0.0);
		  for (cpl_size row = 0; row < nrow; row++) {
			cpl_array_add (flux_mean, flux[tel + row*4]);
		  }

		  /* Normalize to 1 */
		  cpl_array_divide_scalar (flux_mean, cpl_array_get_mean (flux_mean));

		  /* Set the transmission */
		  cpl_table_set_array (oi_flux, "FLUX", tel, flux_mean);

		  /* Set its uncertainty */
		  cpl_array_fill_window (flux_mean, 0, nwave, 0.0);
		  cpl_table_set_array (oi_flux, "FLUXERR", tel, flux_mean);
		  
		} /* End loop on tel*/

		FREE (cpl_array_delete, flux_mean);

		/* Set this transmission to the output data
         * (need to duplicate header but strickly equal) */
		cpl_propertylist * plist = gravi_data_get_oi_flux_plist (p2vmred_data, type_data, pol, npol);
        plist = cpl_propertylist_duplicate (plist);
		gravi_data_add_table (p2vm_map, plist, NULL, oi_flux);
		
		CPLCHECK_MSG("Cannot set data");
	  } /* End loop on pols */	 
	} /* End loop on SC/FT */
	
	gravi_msg_function_exit(1);
	return CPL_ERROR_NONE;
}

/**@}*/
