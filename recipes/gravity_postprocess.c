/* $Id: gravity_postprocess.c,v 1.29 2011/12/3 09:16:12 nazouaoui Exp $
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

#include "gravi_vis.h"
#include "gravi_tf.h"

#include "gravi_preproc.h"

/*-----------------------------------------------------------------------------
                            Private function prototypes
 -----------------------------------------------------------------------------*/

static int gravity_postprocess_create(cpl_plugin *);
static int gravity_postprocess_exec(cpl_plugin *);
static int gravity_postprocess_destroy(cpl_plugin *);
static int gravity_postprocess(cpl_frameset *, const cpl_parameterlist *);

/*-----------------------------------------------------------------------------
                            Static variables
 -----------------------------------------------------------------------------*/

static char gravity_postprocess_short[] = "Post-process the products, to fine-tune their content.";
static char gravity_postprocess_description[] =
    "This recipe allows to manipulate the product of the GRAVITY pipeline, mostly the VIS. It permits to merge several files together into a single VIS file with all observations; to average the observations of one or several VIS file to increse the SNR; to remove some data (FT, SC); and to resample the SC observation with spectral binning.\n"
    "\n"
    "The list of input files can be P2VMRED, VIS, VIS_CALIBRATED (or even RAW for some parameters). However they should all be compatible in term of setup and observed objets !! Note that the recipe performs only litle checks of the input file content and structure. Thus the user shall ensure the input files are conformable (same polarisation and spectral mode for instante)\n"
    GRAVI_RECIPE_FLOW"\n"
    "* Load the files\n"
    "* Execute request from user\n"
    "* Write product\n"
    GRAVI_RECIPE_INPUT"\n"    
    "Input files    : see above\n"
    GRAVI_RECIPE_OUTPUT"\n"
    "POSTPROCESSED          : Output file\n"
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
                    "gravity_postprocess",
                    gravity_postprocess_short,
                    gravity_postprocess_description,
                    "Nabih Azouaoui, Vincent Lapeyrere, JB. Le Bouquin",
                    PACKAGE_BUGREPORT,
                    gravi_get_license(),
                    gravity_postprocess_create,
                    gravity_postprocess_exec,
                    gravity_postprocess_destroy)) {
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
static int gravity_postprocess_create(cpl_plugin * plugin)
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
    
    /* Averaging */
    gravi_parameter_add_average_vis (recipe->parameters);
    gravi_parameter_add_force_uncertainties (recipe->parameters);

    /* Copy fluxdata */
    gravi_parameter_copy_fluxdata (recipe->parameters);
    
    /* Force */
    p = cpl_parameter_new_value ("gravity.postprocess.force-merge", CPL_TYPE_BOOL,
                                 "Force merging even if inconsistent data",
                                 "gravity.postprocess", FALSE);
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "force-merge");
    cpl_parameter_disable (p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append (recipe->parameters, p);
    
    /* Remove FT */
    p = cpl_parameter_new_value ("gravity.postprocess.remove-ft", CPL_TYPE_BOOL,
                                 "Remove FT extensions",
                                 "gravity.postprocess", FALSE);
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "remove-ft");
    cpl_parameter_disable (p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append (recipe->parameters, p);
    
    /* Remove SC */
    p = cpl_parameter_new_value ("gravity.postprocess.remove-sc", CPL_TYPE_BOOL,
                                 "Remove SC extensions",
                                 "gravity.postprocess", FALSE);
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "remove-sc");
    cpl_parameter_disable (p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append (recipe->parameters, p);
    
    /* Remove OPDC */
    p = cpl_parameter_new_value ("gravity.postprocess.remove-opdc", CPL_TYPE_BOOL,
                                 "Remove OPDC extensions",
                                 "gravity.postprocess", FALSE);
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "remove-opdc");
    cpl_parameter_disable (p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append (recipe->parameters, p);
    
    /* Remove MET */
    p = cpl_parameter_new_value ("gravity.postprocess.remove-met", CPL_TYPE_BOOL,
                                 "Remove METROLOGY related extensions",
                                 "gravity.postprocess", FALSE);
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "remove-met");
    cpl_parameter_disable (p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append (recipe->parameters, p);
    
    /* Resamp SC */
    p = cpl_parameter_new_value ("gravity.postprocess.nbin-lambda-sc", CPL_TYPE_INT,
                                 "Bin SC extensions in spectral dimension",
                                 "gravity.postprocess", 0);
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "nbin-lambda-sc");
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
static int gravity_postprocess_exec(cpl_plugin * plugin)
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
    recipe_status = gravity_postprocess(recipe->frames, recipe->parameters);

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
static int gravity_postprocess_destroy(cpl_plugin * plugin)
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
static int gravity_postprocess(cpl_frameset * frameset,
					  const cpl_parameterlist * parlist)
{
    cpl_frameset *used_frameset=NULL;
	cpl_frame * frame=NULL, *frame_merged=NULL;
	
	gravi_data * data_merged=NULL, * data=NULL;

	/* Message */
	gravity_print_banner (); 
	cpl_msg_set_time_on();
	cpl_msg_set_component_on();
	gravi_msg_function_start(1);

    cpl_ensure_code(gravi_dfs_set_groups(frameset) == CPL_ERROR_NONE,
					cpl_error_get_code());

	cpl_size nframe = cpl_frameset_get_size (frameset);
	int force = gravi_param_get_bool (parlist, "gravity.postprocess.force-merge");

	/* To use this recipe the frameset must not be
	 * empty and have at least two frames */
    if ( nframe<1 ) {
	  cpl_error_set_message (cpl_func, CPL_ERROR_ILLEGAL_INPUT,
							 "Not enough frames in the frameset");
	  goto cleanup;
    }

	/* Check the frame have all the same TAG
	 * (cannot merge different type of frame) */

	/* Insert calibration frame into the used frameset */
	used_frameset = cpl_frameset_new();


	/* Loop on other frames to append them*/
	for (cpl_size f = 0; f < nframe ; f++ ) {

	  /* Load the frame */
	  frame = cpl_frameset_get_position (frameset, f);
	  data  = gravi_data_load_frame (frame, used_frameset);      

	  /* Remove some data */
	  if (gravi_param_get_bool (parlist, "gravity.postprocess.remove-ft")) {
		cpl_msg_info (cpl_func, "Erase FT");
		gravi_data_erase_type (data, "_FT");
		CPLCHECK_CLEAN ("Cannot delete FT");
	  }

	  if (gravi_param_get_bool (parlist, "gravity.postprocess.remove-sc")) {
		cpl_msg_info (cpl_func, "Erase SC");
		gravi_data_erase_type (data, "_SC");
		CPLCHECK_CLEAN ("Cannot delete SC");
	  }

	  if (gravi_param_get_bool (parlist, "gravity.postprocess.remove-opdc")) {
		cpl_msg_info (cpl_func, "Erase OPDC");
		gravi_data_erase_type (data, "OPDC");
		CPLCHECK_CLEAN ("Cannot delete OPDC");
	  }
	  
	  if (gravi_param_get_bool (parlist, "gravity.postprocess.remove-met")) {
		cpl_msg_info (cpl_func, "Erase MET");
		gravi_data_erase_type (data, "METROLOGY");
		gravi_data_erase_type (data, "VIS_MET");
		CPLCHECK_CLEAN ("Cannot delete MET");
	  }

      /* Force uncertainties */
      gravi_force_uncertainties (data, parlist);
      CPLCHECK_CLEAN ("Cannot force uncertainties");
      

	  if (f == 0) {
		/* Use the first frame for merging */
		frame_merged = frame;
		data_merged = data; data = NULL;
	  }
	  else {
		/* Merge for first frame */
		gravi_data_append (data_merged, data, force);
		FREE (gravi_data_delete, data);
	  }
	  
	  CPLCHECK_CLEAN ("Cannot append frames");
	}


	/* Co-add them if required (FIXME: and if VIS) */
    if (gravi_param_get_bool (parlist, "gravity.postprocess.average-vis")) {
	  
      gravi_msg_warning ("FIXME", "Average the different observations = EXPERIMENTAL");
	  gravi_average_vis (data_merged);
	  
	  CPLCHECK_CLEAN ("Cannot average VIS");
	}

	/* Resample them if required (FIXME: and if VIS) */
	cpl_size resamp_sc = gravi_param_get_int (parlist, "gravity.postprocess.nbin-lambda-sc");
    if ( resamp_sc > 1) {
	  
	  gravi_msg_warning ("FIXME", "Resamp the SC data = EXPERIMENTAL");
	  gravi_vis_resamp (data_merged, resamp_sc);
	  
	  CPLCHECK_CLEAN ("Cannot resamp SC");
	}

        /* Add the FLUXDATA column for OIFITS2 standard */
    if (gravi_param_get_bool (parlist, "gravity.postprocess.copy-fluxdata"))
    {
      gravi_vis_copy_fluxdata (data_merged);
    }

	/* Recompute the TIME column from the MJD column
	 * in all OIFITS tables to follow standard */
	gravi_vis_mjd_to_time (data_merged);
	
	/* Save the output data file based on the first frame of the frameset */
	gravi_data_save_new (data_merged, frameset, NULL, NULL, parlist,
						 used_frameset, frame_merged, "gravity_postprocess",
						 NULL, "POSTPROCESSED");

	CPLCHECK_CLEAN ("Cannot save the POSTPROCESSED product");

	/* Terminate the function */
	goto cleanup;

cleanup:
	/* Deallocation of all variables */
	cpl_msg_info(cpl_func,"Memory cleanup");
	
	FREE (gravi_data_delete,data);
	FREE (gravi_data_delete,data_merged);
	FREE (cpl_frameset_delete,used_frameset);
	
	gravi_msg_function_exit(1);
    return (int)cpl_error_get_code();
}
