/* $Id: gravity_vis.c,v 1.29 2011/12/3 09:16:12 nazouaoui Exp $
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

/*
 * $Author: nazouaoui $
 * $Date: 2011/12/3 09:16:12 $
 * $Revision: 1.29 $
 * $Name:  $
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/*-----------------------------------------------------------------------------
                                Includes
 -----------------------------------------------------------------------------*/

#include <cpl.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "gravi_data.h"
#include "gravi_pfits.h"
#include "gravi_dfs.h"

#include "gravi_utils.h"

#include "gravi_calib.h"
#include "gravi_p2vmred.h"
#include "gravi_acqcam.h"
#include "gravi_eop.h"
#include "gravi_metrology.h"

#include "gravi_signal.h"
#include "gravi_vis.h"
#include "gravi_tf.h"

#include "gravi_preproc.h"

/*-----------------------------------------------------------------------------
                            Private function prototypes
 -----------------------------------------------------------------------------*/

static int gravity_vis_create(cpl_plugin *);
static int gravity_vis_exec(cpl_plugin *);
static int gravity_vis_destroy(cpl_plugin *);
static int gravity_vis(cpl_frameset *, const cpl_parameterlist *);

/*-----------------------------------------------------------------------------
                            Static variables
 -----------------------------------------------------------------------------*/
static char gravity_vis_short[] = "Compute the visibilities from raw observation of OBJECT.";
static char gravity_vis_description[] =
    "This recipe is associated to the observations template. Its reduces the raw data acquired on calibrator or science targets and computes the uncalibrated visibilities, saved in an OIFITS file. If several OBJECT are provided, the recipe will reduce all of them and merge the resulting data into a single OIFITS. If several SKY_RAW are provided, the recipe reduces the first OBJECT with the first SKY file. Then each new OBJECT with the next SKY. When the number of SKYs is reached, the recipe loops back to first SKY file (so if the number of SKYs is larger than the number of OBJECTs, the last SKY won't be used). The recipe will reduce the data even if no SKY or no DARK is provided. However this will lead to wrong estimate of the visibility and squared visibility of the object. If the file DIAMETER_CAT is not provided, the recipe will use the diameter provided in the header to compute the transfer function QC parameters."
    "\n"
    "The tag in the DO category can be SINGLE/DUAL and CAL/SCI. They should reflect the instrument mode (SINGLE or DUAL) and the DPR.CATG of the observation (SCIENCE or CALIB). The tag in the PRO.CATG category will be SINGLE/DUAL and CAL/SCI depending on the input tag.\n"
    GRAVI_RECIPE_FLOW"\n"
    "* Load the input file (loop on input OBJECT files)\n"
    "* Extract the spectra (use BAD, DARK, SKY, FLAT files)\n"
    "* Interpolate the spectra into a common wavelength table (use WAVE file)\n"
    "* Compute the real-time visibilities (use P2VM file)\n"
    "* Compute additional real-time signals (SNR, GDELAY...)\n"
    "* Compute selection flags (= flag frames with SNR lower than threshold, vFactor lower than threshold...)\n"
    "* Average the real-time visibilities, considering the selection flag\n"
    "* Write the product\n"
    GRAVI_RECIPE_INPUT"\n"    
    GRAVI_FLAT_MAP"               : flat calibration (PRO.CATG="GRAVI_FLAT_MAP")\n"
    GRAVI_BAD_MAP"                : badpixel calibration (PRO.CATG="GRAVI_BAD_MAP") \n"
    GRAVI_WAVE_MAP"               : wave calibration (PRO.CATG="GRAVI_WAVE_MAP")\n"
    GRAVI_P2VM_MAP"               : p2vm calibration (PRO.CATG="GRAVI_P2VM_MAP")\n"
    GRAVI_DARK_MAP"               : dark calibration  (PRO.CATG="GRAVI_DARK_MAP")\n"
    GRAVI_SINGLE_SCIENCE_RAW"     : raw object (DPR.TYPE=OBJECT,SINGLE)\n"
    GRAVI_SINGLE_SKY_RAW"     : raw sky (DPR.TYPE=SKY,SINGLE)\n"
    GRAVI_DISP_MODEL" (opt)   : fiber dispersion model (PRO.CATG="GRAVI_DISP_MODEL")\n"
    GRAVI_DIODE_POSITION" (opt)  : met receiver position (PRO.CATG="GRAVI_DIODE_POSITION")\n"
    GRAVI_DIAMETER_CAT" (opt) : catalog of diameter (PRO.CATG="GRAVI_DIAMETER_CAT")\n"
    GRAVI_RECIPE_OUTPUT"\n"
    GRAVI_VIS_SINGLE_SCIENCE"           : OIFITS file with uncalibrated visibilities\n"
    GRAVI_SINGLE_SKY_MAP" (opt)         : sky map\n"
    GRAVI_P2VMRED_SINGLE_SCIENCE" (opt) : intermediate product (see detailled description of data)\n"
    GRAVI_SPECTRUM" (opt)           : intermediate product (see detailled description of data)\n"
    GRAVI_PREPROC" (opt)            : intermediate product (see detailled description of data)\n"
    "";

/*-----------------------------------------------------------------------------
                                Function code
 -----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
  @brief    Build the list of available plugins, for this module. 
  @param    list    the plugin list
  @return   0 if everything is ok, 1 otherwise
  @note     Only this function is exported

  Create the recipe instance and make it available to the application using the 
  interface. 
 */
/*----------------------------------------------------------------------------*/
int cpl_plugin_get_info(cpl_pluginlist * list)
{
    cpl_recipe  *   recipe = cpl_calloc(1, sizeof *recipe );
    cpl_plugin  *   plugin = &recipe->interface;

    if (cpl_plugin_init(plugin,
                    CPL_PLUGIN_API,
                    GRAVI_BINARY_VERSION,
                    CPL_PLUGIN_TYPE_RECIPE,
                    "gravity_vis",
                    gravity_vis_short,
                    gravity_vis_description,
                    "Nabih Azouaoui, Vincent Lapeyrere, JB. Le Bouquin",
                    PACKAGE_BUGREPORT,
                    gravi_get_license(),
                    gravity_vis_create,
                    gravity_vis_exec,
                    gravity_vis_destroy)) {
        cpl_msg_error(cpl_func, "Plugin initialization failed");
        (void)cpl_error_set_where(cpl_func);                          
        return 1;                                               
    }                                                    

    if (cpl_pluginlist_append(list, plugin)) {                 
        cpl_msg_error(cpl_func, "Error adding plugin to list");
        (void)cpl_error_set_where(cpl_func);                         
        return 1;                                              
    }                                                          
    
    return 0;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Setup the recipe options    
  @param    plugin  the plugin
  @return   0 if everything is ok

  Defining the command-line/configuration parameters for the recipe.
 */
/*----------------------------------------------------------------------------*/
static int gravity_vis_create(cpl_plugin * plugin)
{
    cpl_recipe    * recipe;                                               
    cpl_parameter * p;
                                                                       
    /* Do not create the recipe if an error code is already set */     
    if (cpl_error_get_code() != CPL_ERROR_NONE) {                      
        cpl_msg_error(cpl_func, "%s():%d: An error is already set: %s",
                      cpl_func, __LINE__, cpl_error_get_where());      
        return (int)cpl_error_get_code();                              
    }                                                                  
                                                                       
    if (plugin == NULL) {                                              
        cpl_msg_error(cpl_func, "Null plugin");                        
        cpl_ensure_code(0, (int)CPL_ERROR_NULL_INPUT);                 
    }                                                                  
                                                                       
    /* Verify plugin type */                                           
    if (cpl_plugin_get_type(plugin) != CPL_PLUGIN_TYPE_RECIPE) {       
        cpl_msg_error(cpl_func, "Plugin is not a recipe");             
        cpl_ensure_code(0, (int)CPL_ERROR_TYPE_MISMATCH);              
    }                                                                  
                                                                       
    /* Get the recipe */                                               
    recipe = (cpl_recipe *)plugin;                                     
                                                                       
    /* Create the parameters list in the cpl_recipe object */          
    recipe->parameters = cpl_parameterlist_new();                      
    if (recipe->parameters == NULL) {                                  
        cpl_msg_error(cpl_func, "Parameter list allocation failed");   
        cpl_ensure_code(0, (int)CPL_ERROR_ILLEGAL_OUTPUT);             
    }                                                                  

    /* Fill the parameters list */
    int isCalib = 0;

    /* Use static names (output_procatg.fits) */
    gravi_parameter_add_static_name (recipe->parameters);

    /* Intermediate files */
    gravi_parameter_add_biassub_file (recipe->parameters);
    gravi_parameter_add_spectrum_file (recipe->parameters);
    gravi_parameter_add_preproc_file (recipe->parameters);
    gravi_parameter_add_p2vmred_file (recipe->parameters);
    gravi_parameter_add_astro_file (recipe->parameters);

    /* Averaging */
    gravi_parameter_add_average_vis (recipe->parameters);

    /* Bias-method */
    gravi_parameter_add_biasmethod (recipe->parameters);

    /* Extraction */
    gravi_parameter_add_extract (recipe->parameters);
    gravi_parameter_add_metrology (recipe->parameters);
    
    /* Snr, signal, rejectio flags, vis */
    gravi_parameter_add_compute_snr (recipe->parameters, isCalib);
    gravi_parameter_add_compute_signal (recipe->parameters, isCalib);
    gravi_parameter_add_rejection (recipe->parameters, isCalib);
    gravi_parameter_add_compute_vis (recipe->parameters, isCalib);

    /* Correct from internal transmission */
	p = cpl_parameter_new_value ("gravity.vis.flat-flux", CPL_TYPE_BOOL,
                                 "Normalize the flux (stored in OI_FLUX binary extension) with "
                                 "instrument transmission recorded in the \n"
                                 "input P2VM calibration map. Consequently, the flux quantity is either "
                                 "the intensity level recorded \n"
                                 "in the detector, thus including the instrument transmission (FALSE); "
                                 "or the intensity level at the instrument entrance (TRUE).",
                                 "gravity.vis", FALSE);
	cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "flat-flux");
	cpl_parameter_disable (p, CPL_PARAMETER_MODE_ENV);
	cpl_parameterlist_append (recipe->parameters, p);

    /* Average sky */
	p = cpl_parameter_new_value ("gravity.preproc.average-sky", CPL_TYPE_BOOL,
                                 "Average the SKYs into a master SKY. If FALSE, the recipe loops\n "
                                 "over the SKY to reduce each OBJECT with a different SKY",
                                 "gravity.preproc", FALSE);
	cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "average-sky");
	cpl_parameter_disable (p, CPL_PARAMETER_MODE_ENV);
	cpl_parameterlist_append (recipe->parameters, p);

    /* Reduce ACQ_CAM */
	p = cpl_parameter_new_value ("gravity.test.reduce-acq-cam", CPL_TYPE_BOOL,
                                 "If TRUE, reduced ACQ_CAM images",
                                 "gravity.test", FALSE);
	cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "reduce-acq-cam");
	cpl_parameter_disable (p, CPL_PARAMETER_MODE_ENV);
	cpl_parameterlist_append (recipe->parameters, p);
    
	return 0;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Execute the plugin instance given by the interface
  @param    plugin  the plugin
  @return   0 if everything is ok
 */
/*----------------------------------------------------------------------------*/
static int gravity_vis_exec(cpl_plugin * plugin)
{

    cpl_recipe * recipe;                                                   
    int recipe_status;                                                     
    cpl_errorstate initial_errorstate = cpl_errorstate_get();

                                                                           
    /* Return immediately if an error code is already set */               
    if (cpl_error_get_code() != CPL_ERROR_NONE) {                          
        cpl_msg_error(cpl_func, "%s():%d: An error is already set: %s",    
                      cpl_func, __LINE__, cpl_error_get_where());          
        return (int)cpl_error_get_code();                                  
    }                                                                      
                                                                           
    if (plugin == NULL) {                                                  
        cpl_msg_error(cpl_func, "Null plugin");                            
        cpl_ensure_code(0, (int)CPL_ERROR_NULL_INPUT);                     
    }                                                                      
                                                                           
    /* Verify plugin type */                                               
    if (cpl_plugin_get_type(plugin) != CPL_PLUGIN_TYPE_RECIPE) {           
        cpl_msg_error(cpl_func, "Plugin is not a recipe");                 
        cpl_ensure_code(0, (int)CPL_ERROR_TYPE_MISMATCH);                  
    }                                                                      
                                                                           
    /* Get the recipe */                                                   
    recipe = (cpl_recipe *)plugin;                                         
                                                                           
    /* Verify parameter and frame lists */                                 
    if (recipe->parameters == NULL) {
        cpl_msg_error(cpl_func, "Recipe invoked with NULL parameter list");
        cpl_ensure_code(0, (int)CPL_ERROR_NULL_INPUT);
    }
    if (recipe->frames == NULL) {                                          
        cpl_msg_error(cpl_func, "Recipe invoked with NULL frame set");     
        cpl_ensure_code(0, (int)CPL_ERROR_NULL_INPUT);                     
    }                                                                      
                                                                           
    /* Invoke the recipe */                                                
    recipe_status = gravity_vis(recipe->frames, recipe->parameters);

    /* Ensure DFS-compliance of the products */
    if (cpl_dfs_update_product_header(recipe->frames)) {
        if (!recipe_status){
        	recipe_status = (int)cpl_error_get_code();
        }
    }                                                                      

    if (!cpl_errorstate_is_equal(initial_errorstate)) {                    
        /* Dump the error history since recipe execution start.            
           At this point the recipe cannot recover from the error */       
        cpl_errorstate_dump(initial_errorstate, CPL_FALSE, NULL);          
    }                                                                      

    return recipe_status;                                                  
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Destroy what has been created by the 'create' function
  @param    plugin  the plugin
  @return   0 if everything is ok
 */
/*----------------------------------------------------------------------------*/
static int gravity_vis_destroy(cpl_plugin * plugin)
{
    cpl_recipe * recipe;                                          
                                                                  
    if (plugin == NULL) {                                         
        cpl_msg_error(cpl_func, "Null plugin");                   
        cpl_ensure_code(0, (int)CPL_ERROR_NULL_INPUT);            
    }                                                             
                                                                  
    /* Verify plugin type */
    if (cpl_plugin_get_type(plugin) != CPL_PLUGIN_TYPE_RECIPE) {  
        cpl_msg_error(cpl_func, "Plugin is not a recipe");        
        cpl_ensure_code(0, (int)CPL_ERROR_TYPE_MISMATCH);         
    }                                                             
                                                                  
    /* Get the recipe */
    recipe = (cpl_recipe *)plugin;
                                                                  
    cpl_parameterlist_delete(recipe->parameters);             
                                                                  
    return 0;                                                    
}


/*----------------------------------------------------------------------------*/
/**
  @brief    Compute the visibilities, and closure phase and create the io
  	  	  	fits file
  @param    frameset   the frames list
  @param    parlist    the parameters list
  @return   0 if everything is ok
 */
/*----------------------------------------------------------------------------*/
static int gravity_vis(cpl_frameset * frameset,
					  const cpl_parameterlist * parlist)
{
    cpl_frameset * recipe_frameset=NULL, * wavecalib_frameset=NULL, * dark_frameset=NULL,
	  * darkcalib_frameset=NULL, * sky_frameset=NULL, * flatcalib_frameset=NULL, * p2vmcalib_frameset=NULL,
	  * badcalib_frameset=NULL, *used_frameset=NULL, * current_frameset=NULL, * dispcalib_frameset=NULL,
	  * metpos_frameset=NULL, * diamcat_frameset = NULL, *eop_frameset = NULL;
	
	cpl_frame * frame=NULL;
	
	const char * frame_tag=NULL;
	char * proCatg = NULL, * mode=NULL, * redCatg = NULL, * skyCatg = NULL;
	
	gravi_data * p2vm_map=NULL, * data=NULL, * wave_map=NULL, * dark_map=NULL,
        * profile_map=NULL, * badpix_map=NULL, * preproc_data=NULL, * p2vmred_data=NULL, * tmpvis_data=NULL,
        * vis_data=NULL, * disp_map=NULL, * diodepos_data=NULL, * diamcat_data=NULL, *eop_map=NULL;
	gravi_data ** sky_maps = NULL;
	
	int nb_frame, nb_sky, isky;

	/* Message */
	gravity_print_banner (); 
	cpl_msg_set_time_on();
	cpl_msg_set_component_on();
	gravi_msg_function_start(1);

    cpl_ensure_code(gravi_dfs_set_groups(frameset) == CPL_ERROR_NONE,
					cpl_error_get_code()) ;

    /* Dispatch the frameset */
    p2vmcalib_frameset = gravi_frameset_extract_p2vm_map (frameset);
    darkcalib_frameset = gravi_frameset_extract_dark_map (frameset);
    wavecalib_frameset = gravi_frameset_extract_wave_map (frameset);
    dark_frameset = gravi_frameset_extract_dark_data (frameset);
    flatcalib_frameset = gravi_frameset_extract_flat_map (frameset);
    badcalib_frameset = gravi_frameset_extract_bad_map (frameset);
	dispcalib_frameset = gravi_frameset_extract_disp_map (frameset);
	metpos_frameset = gravi_frameset_extract_met_pos (frameset);
	diamcat_frameset = gravi_frameset_extract_diamcat_map (frameset);
	eop_frameset = gravi_frameset_extract_eop_map (frameset);
	
    recipe_frameset = gravi_frameset_extract_fringe_data (frameset);
    sky_frameset = gravi_frameset_extract_sky_data (frameset);
    
	/* To use this recipe the frameset must contain the p2vm, wave and
	 * gain calibration file. */
    if ( cpl_frameset_get_size (p2vmcalib_frameset) !=1 ||
		 cpl_frameset_get_size (wavecalib_frameset) !=1 ||
		 cpl_frameset_get_size (flatcalib_frameset) !=1 ||
		 cpl_frameset_get_size (badcalib_frameset) != 1 ||
		 cpl_frameset_get_size (recipe_frameset) < 1 ||
		 (cpl_frameset_is_empty (dark_frameset) &&
		  cpl_frameset_is_empty (darkcalib_frameset) &&
		  cpl_frameset_is_empty (sky_frameset)) ) {
	  cpl_error_set_message (cpl_func, CPL_ERROR_ILLEGAL_INPUT,
							 "Illegal number of P2VM, FLAT, WAVE, BAD, DARK or SKY, OBJECT file on the frameset");
	  cpl_error_set_message (cpl_func, CPL_ERROR_ILLEGAL_INPUT,
							 "See online help:  esorex --man gravity_vis");
	  goto cleanup;
    }

	/* Insert calibration frame into the used frameset */
	used_frameset = cpl_frameset_new();

    /* 
	 * Identify the DARK in the input frameset 
	 */
	
    if (!cpl_frameset_is_empty (dark_frameset)) {
	  
    	frame = cpl_frameset_get_position (dark_frameset, 0);
		data = gravi_data_load_rawframe (frame, used_frameset);
        gravi_data_detector_cleanup (data, parlist);

		/* Compute dark */
		dark_map = gravi_compute_dark (data);
		FREE (gravi_data_delete, data);
		
		CPLCHECK_CLEAN ("Could not compute the DARK map");

		/* Save the dark map */
		gravi_data_save_new (dark_map, frameset, NULL, parlist,
							 NULL, frame, "gravity_vis",
                             NULL, GRAVI_DARK_MAP);
		
		CPLCHECK_CLEAN ("Could not save the DARK map");
    }
    else if (!cpl_frameset_is_empty (darkcalib_frameset)) {
	  
    	frame = cpl_frameset_get_position (darkcalib_frameset, 0);
		dark_map = gravi_data_load_frame (frame, used_frameset);

		CPLCHECK_CLEAN ("Could not load the DARK map");
    }
    else
		cpl_msg_info (cpl_func, "There is no DARK in the frame set");

	/* Identify the BAD in the input frameset */
	frame = cpl_frameset_get_position (badcalib_frameset, 0);
	badpix_map = gravi_data_load_frame (frame, used_frameset);
		
    /* Identify the FLAT in the input frameset */
	frame = cpl_frameset_get_position (flatcalib_frameset, 0);
	profile_map = gravi_data_load_frame (frame, used_frameset);

    /* Identify the WAVE in the input frameset */
	frame = cpl_frameset_get_position (wavecalib_frameset, 0);
	wave_map = gravi_data_load_frame (frame, used_frameset);

    /* Identify the P2VM in the input frameset */
	frame = cpl_frameset_get_position (p2vmcalib_frameset, 0);
	p2vm_map = gravi_data_load_frame (frame, used_frameset);

    /* Load the DISP_MODEL in the input frameset */
    if (!cpl_frameset_is_empty (dispcalib_frameset)) {
	  frame = cpl_frameset_get_position (dispcalib_frameset, 0);
	  disp_map = gravi_data_load_frame (frame, used_frameset);
	}
	else
	  cpl_msg_info (cpl_func, "There is no DISP_MODEL in the frameset");

    /* Load the DIODE_POSITION in the input frameset */
    if (!cpl_frameset_is_empty (metpos_frameset)) {
	  frame = cpl_frameset_get_position (metpos_frameset, 0);
	  diodepos_data = gravi_data_load_frame (frame, used_frameset);
	}
	else
	  cpl_msg_info (cpl_func, "There is no DIODE_POSITION in the frameset");

	/* Load the EOP_PARAM */
	if ( !cpl_frameset_is_empty (eop_frameset) ) {
		frame = cpl_frameset_get_position (eop_frameset, 0);
		eop_map = gravi_data_load_frame (frame, used_frameset);
	}
	else
	  cpl_msg_info (cpl_func, "There is no EOP_PARAM in the frameset");
	
	/* Load the DIAMETER_CAT */
	if ( !cpl_frameset_is_empty (diamcat_frameset) ) {
		frame = cpl_frameset_get_position (diamcat_frameset, 0);
		diamcat_data = gravi_data_load_frame (frame, used_frameset);
	}
	else
	  cpl_msg_info (cpl_func, "There is no DIAMETER_CAT in the frameset");

	
	CPLCHECK_CLEAN ("Error while loading the calibration maps");

	/* 
	 * Select the PRO CATG (based on first frame) 
	 */
	
	frame_tag = cpl_frame_get_tag (cpl_frameset_get_position (recipe_frameset, 0));

	if ((strcmp(frame_tag, GRAVI_DUAL_CALIB_RAW) == 0)) {
	  redCatg = cpl_sprintf (GRAVI_P2VMRED_DUAL_CALIB);
	  proCatg = cpl_sprintf (GRAVI_VIS_DUAL_CALIB);
      skyCatg = cpl_sprintf (GRAVI_DUAL_SKY_MAP);
	  mode = cpl_sprintf ("gravi_dual");
	}
	else if ((strcmp(frame_tag, GRAVI_DUAL_SCIENCE_RAW) == 0)) {
	  redCatg = cpl_sprintf (GRAVI_P2VMRED_DUAL_SCIENCE);
	  proCatg = cpl_sprintf (GRAVI_VIS_DUAL_SCIENCE);
      skyCatg = cpl_sprintf (GRAVI_DUAL_SKY_MAP);
	  mode = cpl_sprintf ("gravi_dual"); 
	}
	else if ((strcmp(frame_tag, GRAVI_SINGLE_CALIB_RAW) == 0)) {
	  redCatg = cpl_sprintf (GRAVI_P2VMRED_SINGLE_CALIB);
	  proCatg = cpl_sprintf (GRAVI_VIS_SINGLE_CALIB);
      skyCatg = cpl_sprintf (GRAVI_SINGLE_SKY_MAP);
	  mode = cpl_sprintf ("gravi_single");
	}
	else if ((strcmp(frame_tag, GRAVI_SINGLE_SCIENCE_RAW) == 0)) {
	  redCatg = cpl_sprintf (GRAVI_P2VMRED_SINGLE_SCIENCE);
	  proCatg = cpl_sprintf (GRAVI_VIS_SINGLE_SCIENCE);
      skyCatg = cpl_sprintf (GRAVI_SINGLE_SKY_MAP);
	  mode = cpl_sprintf ("gravi_single");
	}
	else {
        cpl_error_set_message (cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                               "Cannot recognize the input DO.CATG");
        goto cleanup;
	}

	cpl_msg_info (cpl_func,"Mode of the first frame is: %s (will be used for all frames)", mode);

    /*
     * Mode for the SKY
     */
    int averageSky = gravi_param_get_bool (parlist,"gravity.preproc.average-sky");

	/*
	 * Loop on input SKY frames to be reduced 
	 */
	nb_sky   = cpl_frameset_get_size (sky_frameset);
	sky_maps = cpl_calloc (CPL_MAX(nb_sky,1), sizeof(gravi_data*));
    
    for (int isky = 0; isky < nb_sky; isky++){
        
		  /* Load the raw SKY  */
		  cpl_msg_info (cpl_func, " ***** SKY %d over %d ***** ", isky+1, nb_sky);
		  frame = cpl_frameset_get_position (sky_frameset, isky);
		  data = gravi_data_load_rawframe (frame, used_frameset);
          gravi_data_detector_cleanup (data, parlist);

		  /* Compute the SKY map */
		  sky_maps[isky] = gravi_compute_dark (data);
		  FREE (gravi_data_delete, data);
		  
		  CPLCHECK_CLEAN ("Error while computing the sky_map");
		  
		  /* Save the SKY map */
          if (averageSky == 0) {
              gravi_data_save_new (sky_maps[isky], frameset, NULL, parlist,
                                   NULL, frame, "gravity_vis",
                                   NULL, skyCatg);
              CPLCHECK_CLEAN ("Could not save the sky");
          }
    }

    /* 
     * Average the sky if requested
     */

    if (averageSky == 1) {
        cpl_msg_info (cpl_func, "Do a MASTER SKY from the %d skys", nb_sky);
        
        gravi_data * msky_map;
        msky_map = gravi_average_dark (sky_maps, nb_sky);
        CPLCHECK_CLEAN ("Cannot do master sky");

        gravi_data_save_new (msky_map, frameset, NULL, parlist, sky_frameset,
                             cpl_frameset_get_position (sky_frameset, 0),
                             "gravity_vis", NULL, skyCatg);
        CPLCHECK_CLEAN ("Cannot save master sky");

        /* Add all sky to used_frameset, and move pointers */
        cpl_frameset_join (used_frameset, sky_frameset);
        for (int isky = 0; isky < nb_sky; isky++)
            FREE (gravi_data_delete, sky_maps[isky]);
        sky_maps[0] = msky_map;
        nb_sky = 1;
    }

	/*
	 * Loop on input RAW frames to be reduced
	 */
	
    nb_frame = cpl_frameset_get_size (recipe_frameset);
	
    for (int ivis = 0; ivis < nb_frame; ivis++){
		current_frameset = cpl_frameset_duplicate (used_frameset);

		cpl_msg_info (cpl_func, " ***** OBJECT %d over %d ***** ", ivis+1, nb_frame);

		/* 
		 * Identify the SKY for this OBJECT
		 */
		isky = nb_sky>0 ? ivis % nb_sky : 0;
		
		if (nb_sky == 0) {
		  /* No SKY */
		  cpl_msg_info (cpl_func, "There is no SKY in the frameset");
		}
		else if (averageSky) {
		  /* Use master SKY already computed, already in frameset */
		  cpl_msg_info (cpl_func, "Use MASTER SKY (already reduced)");
        }
        else {
		  /* SKY already computed, add in the used_frameset */
		  cpl_msg_info (cpl_func, "Use SKY %i over %i (already reduced)", isky+1, nb_sky);
		  frame = cpl_frameset_get_position (sky_frameset, isky);

		  /* Add this frame to the current_frameset as well */
		  cpl_frameset_insert (current_frameset, cpl_frame_duplicate (frame));
		}

		/* 
		 * Reduce the OBJECT
		 */
		
		frame = cpl_frameset_get_position (recipe_frameset, ivis);
		data = gravi_data_load_rawframe (frame, current_frameset);
        gravi_data_detector_cleanup (data, parlist);

		/* Option save the preproc file */
		if (gravi_param_get_bool (parlist,"gravity.dfs.bias-subtracted-file")) {
		  
			gravi_data_save_new (data, frameset, NULL, parlist,
								 current_frameset, frame, "gravity_vis",
								 NULL, "BIAS_SUBTRACTED");

			CPLCHECK_CLEAN ("Cannot save the BIAS_SUBTRACTED product");
		}
		
		
		/* Check the shutters */
		if ( !gravi_data_check_shutter_open (data) ) {
		  cpl_msg_warning (cpl_func, "Shutter problem in the OBJECT");
		}

        /* Extract spectrum */
        preproc_data = gravi_extract_spectrum (data, profile_map, dark_map, badpix_map,
                                               sky_maps[isky], parlist);
		CPLCHECK_CLEAN ("Cannot extract spectrum");

		/* Option save the spectrum file */
		if (gravi_param_get_bool (parlist,"gravity.dfs.spectrum-file")) {
			gravi_data_save_new (preproc_data, frameset, NULL, parlist,
								 current_frameset, frame, "gravity_vis",
                                 NULL, GRAVI_SPECTRUM);
			CPLCHECK_CLEAN ("Cannot save the SPECTRUM product");
		}
        
        /* Rescale to common wavelength */
        gravi_align_spectrum (preproc_data, wave_map, p2vm_map);
		CPLCHECK_CLEAN ("Cannot re-interpolate spectrum");

        /* Preproc the Acquisition Camera */
        if (gravi_param_get_bool (parlist,"gravity.test.reduce-acq-cam")) {
            gravi_preproc_acqcam (preproc_data, data, badpix_map);
            CPLCHECK_CLEAN ("Cannot preproc ACQ");
        }

		/* Option save the preproc file */
		if (gravi_param_get_bool (parlist,"gravity.dfs.preproc-file")) {
			gravi_data_save_new (preproc_data, frameset, NULL, parlist,
								 current_frameset, frame, "gravity_vis",
                                 NULL, GRAVI_PREPROC);
			CPLCHECK_CLEAN ("Cannot save the PREPROC product");
		}
        
        /* Move extensions from raw_data and delete it */
        gravi_data_move_ext (preproc_data, data, GRAVI_ARRAY_GEOMETRY_EXT);
        gravi_data_move_ext (preproc_data, data, GRAVI_OPTICAL_TRAIN_EXT);
        gravi_data_move_ext (preproc_data, data, GRAVI_OPDC_EXT);
        gravi_data_move_ext (preproc_data, data, GRAVI_FDDL_EXT);
        gravi_data_move_ext (preproc_data, data, GRAVI_METROLOGY_EXT);
		FREE (gravi_data_delete, data);
		CPLCHECK_CLEAN ("Cannot move ext");

		/* Compute the flux and visibilities for each telescope and
		 * per acquisition with the P2VM applied to preproc_data */
		p2vmred_data = gravi_compute_p2vmred (preproc_data, p2vm_map, mode, parlist);
		CPLCHECK_CLEAN ("Cannot apply p2vm to the preproc data");

        /* Reduce the Acquisition Camera and delete data */
        if (gravi_param_get_bool (parlist,"gravity.test.reduce-acq-cam")) {
            gravi_reduce_acqcam (p2vmred_data, preproc_data);
        }
        
        /* Move extensions and delete preproc */
        gravi_data_move_ext (p2vmred_data, preproc_data, GRAVI_IMAGING_DATA_ACQ_EXT);
        gravi_data_move_ext (p2vmred_data, preproc_data, GRAVI_METROLOGY_EXT);
        gravi_data_move_ext (p2vmred_data, preproc_data, GRAVI_FDDL_EXT);
        gravi_data_move_ext (p2vmred_data, preproc_data, GRAVI_OPDC_EXT);
		FREE (gravi_data_delete, preproc_data);
		CPLCHECK_CLEAN ("Cannot delete preproc");

        /* Reduce the OPDC table */
        gravi_compute_opdc_state (p2vmred_data);
		CPLCHECK_CLEAN ("Cannot reduce OPDC");
        
		/* Reduce the metrology into OI_VIS_MET */
		gravi_metrology_reduce (p2vmred_data, eop_map, diodepos_data, parlist);
		CPLCHECK_CLEAN ("Cannot reduce metrology");

		/* Compute the uv and pointing directions with ERFA */
		gravi_compute_uv (p2vmred_data, eop_map);
		CPLCHECK_CLEAN ("Cannot compute uv");

		gravi_compute_pointing (p2vmred_data, eop_map);
		CPLCHECK_CLEAN ("Cannot compute pointing");

		/* Compute the QC0 about tau0 from piezo signals */
		gravi_compute_tau0 (p2vmred_data);

		/* Compute the SNR_BOOT and GDELAY_BOOT */
		gravi_compute_snr (p2vmred_data, parlist);
		CPLCHECK_MSG ("Cannot compute SNR");

		/* Compute the signals for averaging */
		gravi_compute_signals (p2vmred_data, disp_map, parlist);
		CPLCHECK_MSG ("Cannot compute signals");

		/* Compute rejection flags for averaging */
		gravi_compute_rejection (p2vmred_data, parlist);
		CPLCHECK_MSG ("Cannot compute rejection flags signals");
		
		/* Save the p2vmreduced file */
		if (gravi_param_get_bool (parlist,"gravity.dfs.p2vmred-file")) {
			
			gravi_data_save_new (p2vmred_data, frameset, NULL, parlist,
								 current_frameset, frame, "gravity_vis", NULL, redCatg);

			CPLCHECK_CLEAN ("Cannot save the P2VMREDUCED product");
		}

		/* Visibility and flux are averaged and the followings
		 * are saved in tables VIS, VIS2 and T3 */
		tmpvis_data = gravi_compute_vis (p2vmred_data, parlist);
		CPLCHECK_CLEAN ("Cannot average the P2VMRED frames into VIS");

        /* Merge with already existing */
        if (vis_data == NULL) {
            vis_data = tmpvis_data; tmpvis_data = NULL;
        }
        else {
            cpl_msg_info (cpl_func,"Merge with previous OI_VIS");
            gravi_data_append (vis_data, tmpvis_data, 1);
            FREE (gravi_data_delete, tmpvis_data);
        }

        
		/* Save the astro file, which is a lighter version of the p2vmreduced */
		if (gravi_param_get_bool (parlist,"gravity.dfs.astro-file")) {

		    gravi_data_clean_for_astro (p2vmred_data);
			gravi_data_save_new (p2vmred_data, frameset, NULL, parlist,
								 current_frameset, frame, "gravity_vis",
                                 NULL, GRAVI_ASTROREDUCED);

			CPLCHECK_CLEAN ("Cannot save the ASTROREDUCED product");
		}
		
		cpl_msg_info (cpl_func,"Free the p2vmreduced");
		FREE (cpl_frameset_delete, current_frameset);
		FREE (gravi_data_delete, p2vmred_data);
        
    }
    /* End loop on the input files to reduce */

    /* Compute QC parameters */
    gravi_compute_vis_qc (vis_data);

	/* Compute the QC parameters of the TF 
	 * FIXME: compute QC TF only for CALIB star */
	gravi_compute_tf_qc (vis_data, diamcat_data);

    /* Eventually flatten the OI_FLUX */
    if (gravi_param_get_bool (parlist, "gravity.vis.flat-flux")) {
	  
   	  cpl_msg_info (cpl_func, "Flatten the FLUX with the internal P2VM spectrum");
	  gravi_flat_flux (vis_data, p2vm_map);
	  CPLCHECK_CLEAN ("Cannot flat the OI_FLUX");
	  
	} else {
   	  cpl_msg_info (cpl_func, "Don't flatten the FLUX with the internal P2VM spectrum");
	}

	/* Perform the normalisation of the SC vis2 and visamp
	 * to match those of the FT */
    if (!strcmp (gravi_param_get_string (parlist, "gravity.vis.vis-correction-sc"), "FORCE")) {
	  
	  cpl_msg_info (cpl_func, "Align the SC visibilities on the FT");
	  gravi_normalize_sc_to_ft (vis_data);
	  
	} else {
	  cpl_msg_info (cpl_func, "Don't align the SC visibilities on the FT");
	}

	/* Co-add the observations if requested */
    if (gravi_param_get_bool (parlist, "gravity.postprocess.average-vis")) {
	  gravi_average_vis (vis_data);
	  
	} else {
	  cpl_msg_info (cpl_func, "Don't average the different observation (if any)");
	}

	/* Recompute the TIME column from the MJD column
	 * in all OIFITS tables to follow standard */
	gravi_vis_mjd_to_time (vis_data);
	  
	/* Save the output data file based on the first frame of the frameset */
	cpl_frameset_join (used_frameset, recipe_frameset);
	frame = cpl_frameset_get_position (recipe_frameset, 0);
	
	gravi_data_save_new (vis_data, frameset, NULL, parlist,
			     used_frameset, frame, "gravity_vis", NULL, proCatg);

	CPLCHECK_CLEAN ("Cannot save the VIS product");

	/* Terminate the function */
	goto cleanup;

cleanup:
	/* Deallocation of all variables */
	cpl_msg_info(cpl_func,"Memory cleanup");

	FREELOOP (gravi_data_delete,sky_maps,nb_sky);
	FREE (gravi_data_delete,dark_map);
	FREE (gravi_data_delete,data);
	FREE (gravi_data_delete,preproc_data);
	FREE (gravi_data_delete,profile_map);
	FREE (gravi_data_delete,disp_map);
	FREE (gravi_data_delete,wave_map);
	FREE (gravi_data_delete,badpix_map);
	FREE (gravi_data_delete,p2vm_map);
	FREE (gravi_data_delete,p2vmred_data);
	FREE (gravi_data_delete,vis_data);
	FREE (gravi_data_delete,tmpvis_data);
	FREE (gravi_data_delete,diamcat_data);
	FREE (gravi_data_delete,diodepos_data);
	FREE (gravi_data_delete,eop_map);
	FREE (cpl_frameset_delete,darkcalib_frameset);
	FREE (cpl_frameset_delete,wavecalib_frameset);
	FREE (cpl_frameset_delete,flatcalib_frameset);
	FREE (cpl_frameset_delete,badcalib_frameset);
	FREE (cpl_frameset_delete,p2vmcalib_frameset);
	FREE (cpl_frameset_delete,metpos_frameset);
	FREE (cpl_frameset_delete,dark_frameset);
	FREE (cpl_frameset_delete,diamcat_frameset);
	FREE (cpl_frameset_delete,sky_frameset);
	FREE (cpl_frameset_delete,dispcalib_frameset);
	FREE (cpl_frameset_delete,eop_frameset);
	FREE (cpl_frameset_delete,recipe_frameset);
	FREE (cpl_frameset_delete,current_frameset);
	FREE (cpl_frameset_delete,used_frameset);
    FREE (cpl_free,proCatg);
    FREE (cpl_free,redCatg);
    FREE (cpl_free,skyCatg);
    FREE (cpl_free,mode);
	
	gravi_msg_function_exit(1);
    return (int)cpl_error_get_code();
}
