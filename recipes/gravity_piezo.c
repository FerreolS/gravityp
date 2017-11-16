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
#include "gravi_eop.h"
#include "gravi_metrology.h"

#include "gravi_signal.h"
#include "gravi_vis.h"
#include "gravi_tf.h"

#include "gravi_preproc.h"
// #include "gravi_p2vm.h"

/*-----------------------------------------------------------------------------
                            Private function prototypes
 -----------------------------------------------------------------------------*/

static int gravity_piezo_create(cpl_plugin *);
static int gravity_piezo_exec(cpl_plugin *);
static int gravity_piezo_destroy(cpl_plugin *);
static int gravity_piezo(cpl_frameset *, const cpl_parameterlist *);

/*-----------------------------------------------------------------------------
                            Static variables
 -----------------------------------------------------------------------------*/

static char gravity_piezo_short[] = GRAVI_UNOFFERED"Calibrate the response of the piezo actuators.";
static char gravity_piezo_description[] = GRAVI_UNOFFERED"This recipe compute the response (open loop transfer function) of the piezo actuators used to fringe-track in GRAVITY.\n"
    GRAVI_RECIPE_INPUT"\n"    
    GRAVI_FLAT_MAP"               : flat calibration (PRO.CATG="GRAVI_FLAT_MAP")\n"
    GRAVI_BAD_MAP"                : badpixel calibration (PRO.CATG="GRAVI_BAD_MAP") \n"
    GRAVI_WAVE_MAP"               : wave calibration (PRO.CATG="GRAVI_WAVE_MAP")\n"
    GRAVI_P2VM_MAP"               : p2vm calibration (PRO.CATG="GRAVI_P2VM_MAP")\n"
    GRAVI_DARK_MAP"               : dark calibration  (PRO.CATG="GRAVI_DARK_MAP")\n"
    GRAVI_PIEZOTF_RAW"            : dedicated observations \n"
    GRAVI_RECIPE_OUTPUT"\n"    
    GRAVI_PIEZOTF_MAP"           : Respose of the piezo\n"
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
                    "gravity_piezo",
                    gravity_piezo_short,
                    gravity_piezo_description,
                    "Nabih Azouaoui, Vincent Lapeyrere, JB. Le Bouquin",
                    PACKAGE_BUGREPORT,
                    gravi_get_license(),
                    gravity_piezo_create,
                    gravity_piezo_exec,
                    gravity_piezo_destroy)) {
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
static int gravity_piezo_create(cpl_plugin * plugin)
{
    cpl_recipe    * recipe;                                               
    // cpl_parameter * p;
                                                                       
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
    int isCalib = 1;

    /* Use static names (output_procatg.fits) */
    gravi_parameter_add_static_name (recipe->parameters);

    /* Intermediate files */
    gravi_parameter_add_biassub_file (recipe->parameters);
    gravi_parameter_add_preproc_file (recipe->parameters);

    /* Snr, signal, rejection, vis */
    gravi_parameter_add_compute_snr (recipe->parameters, isCalib);
    gravi_parameter_add_compute_signal (recipe->parameters, isCalib);
    gravi_parameter_add_rejection (recipe->parameters, isCalib);
    gravi_parameter_add_compute_vis (recipe->parameters, isCalib);
	
	return 0;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Execute the plugin instance given by the interface
  @param    plugin  the plugin
  @return   0 if everything is ok
 */
/*----------------------------------------------------------------------------*/
static int gravity_piezo_exec(cpl_plugin * plugin)
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
    recipe_status = gravity_piezo(recipe->frames, recipe->parameters);

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
static int gravity_piezo_destroy(cpl_plugin * plugin)
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
static int gravity_piezo(cpl_frameset * frameset,
					  const cpl_parameterlist * parlist)
{
    cpl_frameset * recipe_frameset=NULL, * wavecalib_frameset=NULL, * dark_frameset=NULL,
	  * darkcalib_frameset=NULL, * flatcalib_frameset=NULL, * p2vmcalib_frameset=NULL,
	  * badcalib_frameset=NULL, *used_frameset=NULL, * current_frameset=NULL;
	
	cpl_frame * frame=NULL;
	
	gravi_data * p2vm_map=NULL, * data=NULL, * wave_map=NULL, * dark_map=NULL,
	  * profile_map=NULL, * badpix_map=NULL, * preproc_data=NULL, * p2vmred_data=NULL, * oi_vis=NULL;
	
	int nb_vis, nb_frame;

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
	
    recipe_frameset = gravi_frameset_extract_piezotf_data (frameset);

	/* To use this recipe the frameset must contain the p2vm, wave and
	 * gain calibration file. */
    if ( cpl_frameset_get_size (p2vmcalib_frameset) !=1 ||
		 cpl_frameset_get_size (wavecalib_frameset) !=1 ||
		 cpl_frameset_get_size (flatcalib_frameset) !=1 ||
		 cpl_frameset_get_size (badcalib_frameset) != 1 ||
		 cpl_frameset_get_size (recipe_frameset) < 1 ||
		 (cpl_frameset_is_empty (dark_frameset) &&
		  cpl_frameset_is_empty (darkcalib_frameset) )) {
	  cpl_error_set_message (cpl_func, CPL_ERROR_ILLEGAL_INPUT,
							 "Illegal number of P2VM, FLAT, WAVE, BAD, DARK, PIEZOTF file on the frameset");
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
							 NULL, frame, "gravity_piezo",
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

	/*
	 * Loop on input RAW frames to be reduced 
	 */
	
	nb_vis = 0;
	oi_vis = gravi_data_new(0);
	
    nb_frame = cpl_frameset_get_size (recipe_frameset);
	
    for (int thread = 0; thread < nb_frame; thread++){
		current_frameset = cpl_frameset_duplicate (used_frameset);

		cpl_msg_info (cpl_func, " ***** OBJECT %d over %d ***** ", thread+1, nb_frame);

		/* 
		 * Reduce the OBJECT
		 */
		
		frame = cpl_frameset_get_position (recipe_frameset, thread);
		data = gravi_data_load_rawframe (frame, current_frameset);
        gravi_data_detector_cleanup (data, parlist);

		/* Option save the preproc file */
		if (gravi_param_get_bool (parlist,"gravity.dfs.bias-subtracted-file")) {
		  
			gravi_data_save_new (data, frameset, NULL, parlist,
								 current_frameset, frame, "gravity_piezo",
								 NULL, "BIAS_SUBTRACTED");

			CPLCHECK_CLEAN ("Cannot save the BIAS_SUBTRACTED product");
		}
		
		
		/* Check the shutters */
		if ( !gravi_data_check_shutter_open (data) ) {
		  cpl_msg_warning (cpl_func, "Shutter problem in the OBJECT");
		}

        /* Extract spectrum */
        preproc_data = gravi_extract_spectrum (data, profile_map, dark_map,
                                               badpix_map, NULL, parlist, 
                                               GRAVI_DET_ALL);
		CPLCHECK_CLEAN ("Cannot extract spectrum");
        
        /* Rescale to common wavelength */
        gravi_align_spectrum (preproc_data, wave_map, p2vm_map, GRAVI_DET_ALL, parlist);
		CPLCHECK_CLEAN ("Cannot re-interpolate spectrum");

        /* Move extensions from raw_data */
        gravi_data_move_ext (preproc_data, data, GRAVI_ARRAY_GEOMETRY_EXT);
        gravi_data_move_ext (preproc_data, data, GRAVI_OPTICAL_TRAIN_EXT);
        gravi_data_move_ext (preproc_data, data, GRAVI_OPDC_EXT);
        gravi_data_move_ext (preproc_data, data, GRAVI_FDDL_EXT);
        gravi_data_move_ext (preproc_data, data, GRAVI_METROLOGY_EXT);
		FREE (gravi_data_delete, data);
		CPLCHECK_CLEAN ("Cannot move ext");

		/* Option save the preproc file */
		if (gravi_param_get_bool (parlist,"gravity.dfs.preproc-file")) {
		  
			gravi_data_save_new (preproc_data, frameset, NULL, parlist,
								 current_frameset, frame, "gravity_piezo",
                                 NULL, GRAVI_PREPROC);

			CPLCHECK_CLEAN ("Cannot save the PREPROC product");
		}

		/* Compute the flux and visibilities for each telescope and
		 * per acquisition with the P2VM applied to preproc_data */
		p2vmred_data = gravi_compute_p2vmred (preproc_data, p2vm_map,
		                                "gravi_single", parlist, GRAVI_DET_ALL);
		CPLCHECK_CLEAN ("Cannot apply p2vm to the preproc data");

        /* Move extensions and delete preproc */
        gravi_data_move_ext (p2vmred_data, preproc_data, GRAVI_METROLOGY_EXT);
        gravi_data_move_ext (p2vmred_data, preproc_data, GRAVI_FDDL_EXT);
        gravi_data_move_ext (p2vmred_data, preproc_data, GRAVI_OPDC_EXT);
        FREE (gravi_data_delete, preproc_data);
		CPLCHECK_CLEAN ("Cannot delete preproc");

        /* Reduce the OPDC table */
        gravi_compute_opdc_state (p2vmred_data);
		CPLCHECK_CLEAN ("Cannot reduce OPDC");
        
		/* Reduce the metrology */
		gravi_metrology_reduce (p2vmred_data, NULL, NULL, parlist);
		CPLCHECK_CLEAN ("Cannot reduce metrology");

		/* Compute the SNR_SMT and GDELAY_SMT columns */
		gravi_compute_snr (p2vmred_data, parlist);
		CPLCHECK_MSG ("Cannot compute SNR");

		/* Compute the additional signals for averaging */
		gravi_compute_signals (p2vmred_data, NULL, parlist);
		CPLCHECK_MSG ("Cannot compute signals");

		/* Compute rejection flags for averaging */
		gravi_compute_rejection (p2vmred_data, parlist);
		CPLCHECK_MSG ("Cannot rejection flags signals");

		/* Add here your computation of QC */
		
		/* Save the PIEZOTF which is in fact a P2VMREDUCED */
		gravi_data_save_new (p2vmred_data, frameset, NULL, parlist,
							 current_frameset, frame, "gravity_piezo", NULL, "PIEZOTF");
		
		CPLCHECK_CLEAN ("Cannot save the PIEZOTF product");

		nb_vis++;
		cpl_msg_info (cpl_func,"Free the p2vmreduced");
		FREE (cpl_frameset_delete, current_frameset);
		FREE (gravi_data_delete, p2vmred_data);
    }
	/* End loop on the input files to reduce */

	/* Terminate the function */
	goto cleanup;

cleanup:
	/* Deallocation of all variables */
	cpl_msg_info(cpl_func,"Memory cleanup");
	
	FREE (gravi_data_delete,dark_map);
	FREE (gravi_data_delete,data);
	FREE (gravi_data_delete,preproc_data);
	FREE (gravi_data_delete,profile_map);
	FREE (gravi_data_delete,wave_map);
	FREE (gravi_data_delete,badpix_map);
	FREE (gravi_data_delete,p2vm_map);
	FREE (gravi_data_delete,p2vmred_data);
	FREE (gravi_data_delete,oi_vis);
	FREE (cpl_frameset_delete,darkcalib_frameset);
	FREE (cpl_frameset_delete,wavecalib_frameset);
	FREE (cpl_frameset_delete,flatcalib_frameset);
	FREE (cpl_frameset_delete,badcalib_frameset);
	FREE (cpl_frameset_delete,p2vmcalib_frameset);
	FREE (cpl_frameset_delete,dark_frameset);
	FREE (cpl_frameset_delete,recipe_frameset);
	FREE (cpl_frameset_delete,current_frameset);
	FREE (cpl_frameset_delete,used_frameset);
	
	gravi_msg_function_exit(1);
    return (int)cpl_error_get_code();
}
