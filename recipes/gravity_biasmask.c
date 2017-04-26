/* $Id: gravity_biasmask.c,v 1.29 2009/02/10 09:16:12 llundin Exp $
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


/*-----------------------------------------------------------------------------
                            Private function prototypes
 -----------------------------------------------------------------------------*/

static int gravity_biasmask_create(cpl_plugin *);
static int gravity_biasmask_exec(cpl_plugin *);
static int gravity_biasmask_destroy(cpl_plugin *);
static int gravity_biasmask(cpl_frameset *, const cpl_parameterlist *);

/*-----------------------------------------------------------------------------
                            Static variables
 -----------------------------------------------------------------------------*/

static char gravity_biasmask_short[] = GRAVI_UNOFFERED"Determine which pixels can be used to measure the bias of SC detector.";
static char gravity_biasmask_description[] = GRAVI_UNOFFERED"The recipe creates a binary mask (BIASPIX) indentifying which pixels of the SC detector are not illuminated, and thus could be used as bias-pixels in further processing. The idea would be to input such a mask, as static calibration, in all reductions. However this is not yet implemented, nor demonstrated as necessary.\n"
    GRAVI_RECIPE_FLOW"\n"
    "* Load the input files\n"
    "* Identify the mask\n"
    "* Write product\n"
    GRAVI_RECIPE_INPUT"\n"    
    GRAVI_DARK_RAW"      : raw dark, all shutters closed (DPR.TYPE=DARK)\n"
    GRAVI_FLAT_RAW"  x4  : raw flats, one sutter open (DPR.TYPE=FLAT)\n"
    GRAVI_RECIPE_OUTPUT"\n"
    GRAVI_BIASMASK_MAP"           : biaspixel mask calibration \n"
    "\n";

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
                    "gravity_biasmask",
                    gravity_biasmask_short,
                    gravity_biasmask_description,
                    "Nabih Azouaoui, Vincent Lapeyrere, JB. Le Bouquin",
                    PACKAGE_BUGREPORT,
                    gravi_get_license(),
                    gravity_biasmask_create,
                    gravity_biasmask_exec,
                    gravity_biasmask_destroy)) {
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
static int gravity_biasmask_create(cpl_plugin * plugin)
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

    /* Use static names (output_procatg.fits) */
    gravi_parameter_add_static_name (recipe->parameters);
	
    return 0;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Execute the plugin instance given by the interface
  @param    plugin  the plugin
  @return   0 if everything is ok
 */
/*----------------------------------------------------------------------------*/
static int gravity_biasmask_exec(cpl_plugin * plugin)
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
    recipe_status = gravity_biasmask(recipe->frames, recipe->parameters);
                                                                           
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
static int gravity_biasmask_destroy(cpl_plugin * plugin)
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
  @brief    Compute the DARK, BAD, FLAT, WAVE, P2VM from a list of
            calibration set.
  @param    frameset   the frames list
  @param    parlist    the parameters list
  @return   0 if everything is ok
 */
/*----------------------------------------------------------------------------*/
static int gravity_biasmask(cpl_frameset            * frameset,
		                  const cpl_parameterlist * parlist)
{
	cpl_frameset * dark_frameset=NULL, * flat_frameset=NULL, * used_frameset=NULL;
	
	cpl_frame * frame=NULL;
	
	gravi_data * data = NULL, * dark_map=NULL;
	gravi_data ** raw_data=NULL, * biasmask_map=NULL;

	int nb_frame_gain = 0;
	
	/* Message */
	gravity_print_banner (); 
	cpl_msg_set_time_on();
	cpl_msg_set_component_on();
	gravi_msg_function_start(1);

	/* Get the input frameset */
	cpl_ensure_code(gravi_dfs_set_groups(frameset) == CPL_ERROR_NONE, cpl_error_get_code()) ;

	/* Init the used frameset */
	used_frameset = cpl_frameset_new ();

	/* Extract DARK frameset */
	dark_frameset = gravi_frameset_extract_dark_data (frameset);
	
	/* Extract FLAT frameset */
	flat_frameset = gravi_frameset_extract_flat_data (frameset);

	/* To use this recipe the frameset must contain the p2vm, wave and
	 * gain calibration file. */
    if ( cpl_frameset_is_empty (dark_frameset) ||
		 cpl_frameset_get_size (dark_frameset) != 1 || 
		 cpl_frameset_is_empty (flat_frameset) ||
		 cpl_frameset_get_size (flat_frameset) != 4 ) {
	  cpl_error_set_message (cpl_func, CPL_ERROR_ILLEGAL_INPUT,
							 "Need 1 DARK_RAW and 4 FLAT_RAW");
	  goto cleanup;
    }
	
	
	/*
	 * (1) Identify and extract the dark file
	 */

	cpl_msg_info (cpl_func, " ***** Compute DARK map ***** ");

	/* Load this DARK_RAW */
	frame = cpl_frameset_get_position (dark_frameset, 0);
	data = gravi_data_load_rawframe (frame, used_frameset);
    gravi_data_detector_cleanup (data, parlist);

	/* Compute the dark */
	dark_map = gravi_compute_dark (data);
	FREE (gravi_data_delete, data);

	CPLCHECK_CLEAN ("Cannot compute the DARK map");

	/*
	 * (2) Load the FLAT files
	 */
    
	cpl_msg_info (cpl_func, " ***** Load FLATs ***** ");
	
	/* Identify the flat files */
	nb_frame_gain = cpl_frameset_get_size (flat_frameset);
	raw_data = cpl_calloc (nb_frame_gain, sizeof(gravi_data *));
	
	/* Build the list of FLAT files and output file name */
	for (int i = 0; i < nb_frame_gain; i++) {
		 frame = cpl_frameset_get_position (flat_frameset, i);
		 raw_data[i] = gravi_data_load_rawframe (frame, used_frameset);
         gravi_data_detector_cleanup (raw_data[i], parlist);
	}
    	  
	/*
	 * (2) Compute the BIASMASK from the DARK and FLATs
	 */

	cpl_msg_info (cpl_func, " ***** Compute BIAS_MASK map ***** ");
	biasmask_map = gravi_compute_biasmask (dark_map, raw_data,
                                           nb_frame_gain, parlist);

	CPLCHECK_CLEAN("Cannot compute the BIAS_MASK");
	
	/* Free the list of files */
	FREELOOP (gravi_data_delete, raw_data, nb_frame_gain);
	
	/* Save the BAD */
	frame = cpl_frameset_get_position (dark_frameset, 0);
	gravi_data_save_new (biasmask_map, frameset, NULL, parlist,
						 NULL, frame, "gravity_biasmask",
                         NULL, GRAVI_BIASMASK_MAP);

	CPLCHECK_CLEAN ("Could not save the BAD pixel map");
	
    /* Deallocation of all variables */
 cleanup:
	cpl_msg_info (cpl_func,"Cleanup memory");
	
	FREE (cpl_frameset_delete, dark_frameset);
	FREE (gravi_data_delete, dark_map);
	FREELOOP (gravi_data_delete, raw_data, nb_frame_gain);
	FREE (gravi_data_delete, biasmask_map);
	FREE (cpl_frameset_delete, flat_frameset);
	FREE (cpl_frameset_delete, used_frameset);
	FREE (gravi_data_delete, data);

	/* FIXME: check a *change* of cpl_state instead */
	CPLCHECK_INT ("Could not cleanup memory");

	gravi_msg_function_exit(1);
    return (int)cpl_error_get_code();
}


