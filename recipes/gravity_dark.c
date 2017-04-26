/* $Id: gravity_dark.c,v 1.29 2011/03/10 09:16:12 nazouaoui Exp $
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
 * $Date: 2011/03/10 09:16:12 $
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

#include "gravi_utils.h"
#include "gravi_pfits.h"
#include "gravi_dfs.h"
#include "gravi_calib.h"

#include "gravi_data.h"

/*-----------------------------------------------------------------------------
                            Private function prototypes
 -----------------------------------------------------------------------------*/

static int gravity_dark_create(cpl_plugin *);
static int gravity_dark_exec(cpl_plugin *);
static int gravity_dark_destroy(cpl_plugin *);
static int gravity_dark(cpl_frameset *, const cpl_parameterlist *);

/*-----------------------------------------------------------------------------
                            Static variables
 -----------------------------------------------------------------------------*/

static char gravity_dark_short[] = "Calibrate the detector noise and background level.";
static char gravity_dark_description[] =
    "This recipe computes the DARK calibration for the SC, the FT and the ACQ detectors. The SC detector is first debiased using the biaspixels, before computing the dark mean and rms. For detectors, the mean dark level of each pixel and the stdev of each pixel are saved in the output product.\n"
    GRAVI_RECIPE_FLOW"\n"
    "* Loop on input dark files and concatenate them\n"
    "* Compute the median and rms of these concatenated files\n"
    "* Save the product (FT, SC, ACQ camera into same product)\n"
    GRAVI_RECIPE_INPUT"\n"    
    GRAVI_DARK_RAW"      : raw dark, all shutters closed (DPR.TYPE=DARK)\n"
    GRAVI_RECIPE_OUTPUT"\n"    
    GRAVI_DARK_MAP"          : dark calibration\n"
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
                    "gravity_dark",
                    gravity_dark_short,
                    gravity_dark_description,
                    "Nabih Azouaoui, Vincent Lapeyrere, JB. Le Bouquin",
                    PACKAGE_BUGREPORT,
                    gravi_get_license(),
                    gravity_dark_create,
                    gravity_dark_exec,
                    gravity_dark_destroy)) {    
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
static int gravity_dark_create(cpl_plugin * plugin)
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

    /* Use static names (output_procatg.fits) */
    gravi_parameter_add_static_name (recipe->parameters);

    /* Bias-method */
    gravi_parameter_add_biasmethod (recipe->parameters);
    
    return 0;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Execute the plugin instance given by the interface
  @param    plugin  the plugin
  @return   0 if everything is ok
 */
/*----------------------------------------------------------------------------*/
static int gravity_dark_exec(cpl_plugin * plugin)
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
    recipe_status = gravity_dark(recipe->frames, recipe->parameters);
                                                                           
    /* Ensure DFS-compliance of the products */                            
    if (cpl_dfs_update_product_header(recipe->frames)) {                   
        if (!recipe_status) recipe_status = (int)cpl_error_get_code();                         
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
static int gravity_dark_destroy(cpl_plugin * plugin)
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
  @brief    Compute the master dark for each dark frames
  @param    frameset   the frames list
  @param    parlist    the parameters list
  @return   0 if everything is ok
 */
/*----------------------------------------------------------------------------*/
static int gravity_dark (cpl_frameset            * frameset,
						   const cpl_parameterlist * parlist)
{
	cpl_frameset * dark_frameset = NULL, * used_frameset = NULL;
	cpl_frame * frame = NULL;
	gravi_data * raw_dark = NULL, * reduced_dark = NULL;
	int comp, nb_frame;

	
	/* Message */
	gravity_print_banner (); 
	cpl_msg_set_time_on();
	cpl_msg_set_component_on();
	gravi_msg_function_start(1);
 
	/* Identify the RAW and CALIB frames in the input frameset */
    cpl_ensure_code (gravi_dfs_set_groups (frameset) == CPL_ERROR_NONE, cpl_error_get_code());

    /* Dispatch the framesets */
    dark_frameset = gravi_frameset_extract_dark_data (frameset);

	/* Check the frameset */
    if (cpl_frameset_is_empty (dark_frameset)) {
	  cpl_error_set_message (cpl_func, CPL_ERROR_ILLEGAL_INPUT,
							 "No DARK_RAW file on the frameset") ;
	  goto cleanup;
    }

    /* Load all DARK in the frameset into a single data */
    nb_frame = cpl_frameset_get_size (dark_frameset);
	used_frameset = cpl_frameset_new ();

    for (comp = 0; comp < nb_frame; comp++){
	    gravi_data * data_tmp;
	    cpl_frame * frame_tmp;

        /* Load this frame */
        frame_tmp = cpl_frameset_get_position (dark_frameset, comp);
        data_tmp  = gravi_data_load_rawframe (frame_tmp, used_frameset);
        gravi_data_detector_cleanup (data_tmp, parlist);

		/* Cleanup unused data */
		gravi_data_erase (data_tmp, GRAVI_METROLOGY_EXT);
		gravi_data_erase (data_tmp, GRAVI_OPDC_EXT);
		gravi_data_erase (data_tmp, GRAVI_FDDL_EXT);

		/* FIXME: shall remove the FT except for the first one */

		if (comp == 0) {
		  /* Use the first frame for merging */
		  frame = frame_tmp;
		  raw_dark = data_tmp; data_tmp = NULL;
		}
		else {
		  /* Merge to first frame */
		  int force = 0;
		  gravi_data_append (raw_dark, data_tmp, force);
		  FREE (gravi_data_delete, data_tmp);
		}
		
		CPLCHECK_CLEAN ("Cannot load all DARK into a single data");
    }

	

	/* Compute the reduced DARK */
	reduced_dark = gravi_compute_dark (raw_dark);
	FREE (gravi_data_delete, raw_dark);
	
	CPLCHECK_CLEAN ("Could not compute the DARK map");
	
	/* Create product frame */
	gravi_data_save_new (reduced_dark, frameset, NULL, parlist,
						 used_frameset, frame, "gravity_dark",
                         NULL, GRAVI_DARK_MAP);
	
	CPLCHECK_CLEAN ("Could not save the DARK map");

	/* Terminate the function */
	goto cleanup;

cleanup:
	/* Deallocation of all variables */
	cpl_msg_info (cpl_func,"Memory cleanup");
	
	FREE (gravi_data_delete, reduced_dark);
	FREE (gravi_data_delete, raw_dark);
    FREE (cpl_frameset_delete, dark_frameset);
	FREE (cpl_frameset_delete, used_frameset);
	
	gravi_msg_function_exit(1);
    return (int)cpl_error_get_code();
}


