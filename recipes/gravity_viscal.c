/* $Id: gravity_viscal.c,v 1.29 2011/12/3 09:16:12 nazouaoui Exp $
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
#include <math.h>
#include <time.h>

#include "gravi_utils.h"
#include "gravi_pfits.h"
#include "gravi_data.h"
#include "gravi_dfs.h"
#include "gravi_calib.h"
#include "gravi_p2vmred.h"
#include "gravi_vis.h"
#include "gravi_tf.h"

/*-----------------------------------------------------------------------------
                            Private function prototypes
 -----------------------------------------------------------------------------*/

static int gravity_viscal_create(cpl_plugin *);
static int gravity_viscal_exec(cpl_plugin *);
static int gravity_viscal_destroy(cpl_plugin *);
static int gravity_viscal(cpl_frameset *, const cpl_parameterlist *);

/*-----------------------------------------------------------------------------
                            Static variables
 -----------------------------------------------------------------------------*/

static char gravity_viscal_short[] = "Calibrate visibilities from the transfer function.";
static char gravity_viscal_description[] =
    "This recipe calibrates the visibilities acquired on science target using visibilities acquired on calibrator target. If the DIAMETER_CAT is not provided, the recipe will use the diameter provided in the header to compute the transfer function QC parameters. The corresponding keywords are INS.SOBJ.DIAMETER and FT.ROBJ.DIAMETER. The OI_FLUX data are not yet calibrated."
    "\n"
    "The tag in the DO category can be SINGLE/DUAL and CAL/SCI. They should reflect the mode (SINGLE or DUAL) and the DPR.CATG of the observation (SCIENCE or CALIB). The tag in the PRO.CATG category will be SINGLE/DUAL and CAL/SCI depending on the input tag.\n"
    GRAVI_RECIPE_FLOW"\n"
    "* Loop on all input CALIB files, compute the TF for each of them and write the corresponding product\n"
    "* Loop on all input SCIENCE files, interpolate the TF at that time, calibrate, and write the corresponding product\n"
    GRAVI_RECIPE_INPUT"\n"    
    GRAVI_VIS_SINGLE_SCIENCE" (>=1) : visibilities on sciences\n"
    GRAVI_VIS_SINGLE_CALIB" (>=1) : visibilities on calibrators\n"
    GRAVI_DIAMETER_CAT" (opt)   : catalog of diameter\n"
    GRAVI_RECIPE_OUTPUT"\n"    
    GRAVI_VIS_SINGLE_CALIBRATED" : calibrated science visibilities\n"
    GRAVI_TF_SINGLE_CALIB" :  Transfer Function (TF) estimated on calibrators\n"
    GRAVI_TF_SINGLE_SCIENCE" :  TF interpolated at the time of sciences\n"
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
                    "gravity_viscal",
                    gravity_viscal_short,
                    gravity_viscal_description,
                    "Nabih Azouaoui, Vincent Lapeyrere, JB. Le Bouquin",
                    PACKAGE_BUGREPORT,
                    gravi_get_license(),
                    gravity_viscal_create,
                    gravity_viscal_exec,
                    gravity_viscal_destroy)) {
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
static int gravity_viscal_create(cpl_plugin * plugin)
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
    gravi_parameter_add_static_name (recipe->parameters);
    
    /* delta-time */
    p = cpl_parameter_new_value ("gravity.viscal.delta-time-calib", CPL_TYPE_DOUBLE,
                                 "Delta time to interpolate the TF [s]",
                                 "gravity.viscal", 3600.0);
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "delta-time-calib");
    cpl_parameter_disable (p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append (recipe->parameters, p);
    
    /* force-calib */
    p = cpl_parameter_new_value ("gravity.viscal.force-calib", CPL_TYPE_BOOL,
                                 "Force the calibration, don't check setup",
                                 "gravity.viscal", FALSE);
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "force-calib");
    cpl_parameter_disable (p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append (recipe->parameters, p);

    /* smooth */
    p = cpl_parameter_new_value ("gravity.viscal.nsmooth-tfvis-sc", CPL_TYPE_INT,
                                 "Smooth the TF spectrally by this number of "
                                 "spectral bin, to enhance SNR if needed (only "
                                 "apply to VIS2, VISPHI, VISAMP, T3PHI, T3AMP).",
                                 "gravity.viscal", 0);
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "nsmooth-tfvis-sc");
    cpl_parameter_disable (p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append (recipe->parameters, p);

    p = cpl_parameter_new_value ("gravity.viscal.nsmooth-tfflux-sc", CPL_TYPE_INT,
                                 "Smooth the TF spectrally by this number of "
                                 "spectral bin, to enhance SNR if needed (only "
                                 "apply to FLUX, RVIS, IVIS).",
                                 "gravity.viscal", 0);
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "nsmooth-tfflux-sc");
    cpl_parameter_disable (p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append (recipe->parameters, p);

    p = cpl_parameter_new_value ("gravity.viscal.maxdeg-tfvis-sc", CPL_TYPE_INT,
                                 "Fit the TF spectrally by a polynomial to enhance SNR "
                                 "(only apply to VIS2, VISPHI, VISAMP, T3PHI, T3AMP).",
                                 "gravity.viscal", -1);
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "maxdeg-tfvis-sc");
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
static int gravity_viscal_exec(cpl_plugin * plugin)
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
    recipe_status = gravity_viscal(recipe->frames, recipe->parameters);

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
static int gravity_viscal_destroy(cpl_plugin * plugin)
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
  @brief    Compute the visibilities, and cloture phase and create the io
  	  	  	fits file
  @param    frameset   the frames list
  @param    parlist    the parameters list
  @return   0 if everything is ok
 */
/*----------------------------------------------------------------------------*/
static int gravity_viscal(cpl_frameset            * frameset,
		                  const cpl_parameterlist * parlist)
{
	cpl_frameset * vis_calib_frameset = NULL, * vis_sci_frameset = NULL, *current_frameset = NULL;
	cpl_frameset * tf_calib_frameset = NULL, * used_frameset = NULL, * diamcat_frameset = NULL;
	cpl_frame * frame = NULL;
	
	cpl_propertylist * applist = NULL;

	cpl_errorstate errorstate;
	
	gravi_data ** vis_calibs = NULL, *vis_calib = NULL, * zero_data = NULL, * tf_science = NULL;
	gravi_data * calibrated = NULL, * vis_data = NULL, * diamcat_data = NULL;
	
	int data_mode, nb_frame_tf = 0, nb_frame_calib = 0, nb_frame_sci, i, j, nb_calib = 0;

	/* Message */
	gravity_print_banner (); 
	cpl_msg_set_time_on();
	cpl_msg_set_component_on();
	gravi_msg_function_start(1);

	/* Identify the RAW and CALIB frames in the input frameset */
    cpl_ensure_code(gravi_dfs_set_groups(frameset) == CPL_ERROR_NONE,
					cpl_error_get_code());

    /* Extract a set of vis_calib and vis_sci data frameset */
    vis_calib_frameset = gravi_frameset_extract_vis_calib (frameset);
    vis_sci_frameset = gravi_frameset_extract_vis_science (frameset);
    tf_calib_frameset = gravi_frameset_extract_tf_calib (frameset);
	diamcat_frameset = gravi_frameset_extract_diamcat_map (frameset);


	/* To use this recipe the frameset must contain
	 * at least one VIS_*_CAL frame or TF_*_CAL frame. */
	if ( cpl_frameset_is_empty (vis_calib_frameset) &&
		 cpl_frameset_is_empty (tf_calib_frameset) ) {
	  cpl_error_set_message (cpl_func, CPL_ERROR_ILLEGAL_INPUT,
							"No VIS or TF file on the frameset") ;
	  goto cleanup;
    }
	
    /* Get the number of the frames contained in the frameset */
	nb_frame_tf = cpl_frameset_get_size (tf_calib_frameset);
    nb_frame_calib = cpl_frameset_get_size (vis_calib_frameset);
    nb_frame_sci = cpl_frameset_get_size (vis_sci_frameset);

	/* Init memory */
	vis_calibs = cpl_malloc ( (nb_frame_calib + nb_frame_tf) * sizeof (gravi_data *));
	for (j = 0; j < (nb_frame_calib + nb_frame_tf); j++) vis_calibs[j] = NULL;

	/* 
	 * Load or compute of the transfer function 
	 */

	/* Load the DIAMETER_CAT 
	 * FIXME: how/when install in the used_frameset ? */
	if ( !cpl_frameset_is_empty (diamcat_frameset) ) {
		frame = cpl_frameset_get_position (diamcat_frameset, 0);
		diamcat_data = gravi_data_load_frame (frame, NULL);
	}

	/* Init the frameset for all calibration (TF computed and loaded) */
	used_frameset = cpl_frameset_new();

	/* Loop on the TF to compute */
	for (j = 0; j < nb_frame_calib; j++) {
		errorstate = cpl_errorstate_get();
		
		cpl_msg_info (cpl_func, "*** Compute TF %i over %i ***", j+1, nb_frame_calib);

		/* Load the VIS data and compute TF */
		frame = cpl_frameset_get_position (vis_calib_frameset, j);
		vis_data = gravi_data_load_frame (frame, NULL);
		vis_calib = gravi_compute_tf (vis_data, diamcat_data);

        /* Smooth the TF if required */
        cpl_size smooth_vis_sc = gravi_param_get_int (parlist, "gravity.viscal.nsmooth-tfvis-sc");
        cpl_size smooth_flx_sc = gravi_param_get_int (parlist, "gravity.viscal.nsmooth-tfflux-sc");
        cpl_size maxdeg_sc = gravi_param_get_int (parlist, "gravity.viscal.maxdeg-tfvis-sc");
        gravi_vis_smooth (vis_calib, smooth_vis_sc, smooth_flx_sc, maxdeg_sc);

		/* Check the TF has been computed */
		if (vis_calib == NULL) {
		  cpl_msg_error (cpl_func, "Cannot compute this TF... continue");
		  goto cleanup_rawtf;
		}

		CPLCHECK_GOTO ("Cannot compute the TF", cleanup_rawtf);

		/* Save the TF file */
		data_mode = gravi_data_frame_get_mode (frame);
		
		gravi_data_save_new (vis_calib, frameset, NULL, parlist,
							 NULL, frame, "gravity_vis",
							 NULL, GRAVI_TF_CALIB(data_mode));

		CPLCHECK_GOTO ("Cannot save the TF", cleanup_rawtf);

		/* Store this successfull TF */
		vis_calibs[nb_calib] = vis_calib;
		vis_calib = NULL;
		nb_calib++;

		/* Update the used_frameset -- now used as calibration */
		frame = cpl_frame_duplicate (frame);
		cpl_frame_set_group	(frame, CPL_FRAME_GROUP_CALIB);
		cpl_frameset_insert (used_frameset, frame);
		
		/* Clean memory of the loop */
	cleanup_rawtf:
		FREE (gravi_data_delete, vis_data);
		FREE (gravi_data_delete, vis_calib);
		cpl_errorstate_set (errorstate);
	}
	/* End loop on the TF to compute */

	cpl_msg_info (cpl_func,"*** Load already computed TF ***");
	
	/* Loop on the TF to load */
	for (j = 0; j < nb_frame_tf; j++) {
	    errorstate = cpl_errorstate_get();
	  
		cpl_msg_info (cpl_func," %i over %i", j+1, nb_frame_tf);
		
		frame = cpl_frameset_get_position (tf_calib_frameset, j);
		vis_calib = gravi_data_load_frame (frame, used_frameset);

		CPLCHECK_GOTO("Cannot load the TF", cleanup_caltf);

		/* Store this successfull TF */
		vis_calibs[nb_calib] = vis_calib;
		vis_calib = NULL;
		nb_calib++;

	cleanup_caltf:
		FREE (gravi_data_delete,vis_calib);
		cpl_errorstate_set (errorstate);
	}
	/* End loop on TF to load */
	
	cpl_msg_info (cpl_func,"*** All TF computed or loaded ***");

	cpl_msg_info (cpl_func, "Load or create successfully %i TF over %i input CAL files", nb_calib, nb_frame_calib + nb_frame_tf);
	
	/* 
	 * Compute the zero of the metrology if several TF are availables
	 */

	if ( nb_calib > 1 ) { 
	  errorstate = cpl_errorstate_get();
	  
	  cpl_msg_info (cpl_func, "*** Compute the zero of the metrology -- FIXME: to be done ***");
	  zero_data = gravi_compute_zp (vis_calibs, nb_calib);
	  
	  CPLCHECK_GOTO("Cannot compute ZP", cleanup_zp);
	  
	  gravi_data_save_new (zero_data, frameset, "output.fits", parlist,
						   used_frameset, NULL, "gravity_vis",
                           NULL, GRAVI_ZP_CAL);

	  CPLCHECK_GOTO("Cannot save ZP", cleanup_zp);
	  
	cleanup_zp:
	  FREE (gravi_data_delete,zero_data);
	  cpl_errorstate_set (errorstate);
	}

    	
	/* 
	 * Apply the TF to the SCIENCE files of the frameset
	 */
	
	/* Loop on the SCI files to calibrate */
	for (i = 0; i < nb_frame_sci; i++){
		errorstate = cpl_errorstate_get();
		current_frameset = cpl_frameset_duplicate (used_frameset);
		
		cpl_msg_info (cpl_func, "*** Calibration of file %i over %i ***", i+1, nb_frame_sci);
		
		frame = cpl_frameset_get_position (vis_sci_frameset, i);
		vis_data = gravi_data_load_frame (frame, current_frameset);
		
		tf_science = gravi_data_duplicate (vis_data);
		calibrated = gravi_calibrate_vis (vis_data, vis_calibs, nb_calib, zero_data, tf_science, parlist);
		
		CPLCHECK_GOTO("Cannot calibrate the visibility", cleanup_calib);

		/* Save calibrated visibilities */
		data_mode = gravi_data_frame_get_mode (frame);
		
		gravi_data_save_new (calibrated, frameset, NULL, parlist,
							 current_frameset, frame, "gravity_vis",
							 NULL, GRAVI_VIS_CALIBRATED(data_mode));
		
		CPLCHECK_GOTO("Cannot save the calibrated visibility", cleanup_calib);

		/* Save TF interpolated at the science visibilities */
		data_mode = gravi_data_frame_get_mode (frame);
		
		gravi_data_save_new (tf_science, frameset, NULL, parlist,
							 current_frameset, frame, "gravity_vis", NULL,
							 GRAVI_TF_SCIENCE(data_mode));
		
		CPLCHECK_GOTO("Cannot save the TF interpolated for this visibility", cleanup_calib);

	cleanup_calib:
		FREE (gravi_data_delete,vis_data);
		FREE (gravi_data_delete,tf_science);
		FREE (gravi_data_delete,calibrated);
		FREE (cpl_propertylist_delete,applist);
		FREE (cpl_frameset_delete,current_frameset);
		cpl_errorstate_set (errorstate);
	}
	/* End loop on VIS_*_SCI files to calibrate */
	
	/* Terminate the function */
	goto cleanup;

cleanup:
	/* Deallocation of all variables */
	cpl_msg_info(cpl_func,"Memory cleanup");
	
	FREE (gravi_data_delete,diamcat_data);
	FREE (gravi_data_delete,zero_data);
    FREE (cpl_frameset_delete,tf_calib_frameset);
	FREE (cpl_frameset_delete,vis_calib_frameset);
	FREE (cpl_frameset_delete,vis_sci_frameset);
	FREE (cpl_frameset_delete,diamcat_frameset);
	FREE (cpl_frameset_delete,current_frameset);
	FREE (cpl_frameset_delete,used_frameset);
	FREELOOP (gravi_data_delete,vis_calibs,nb_frame_calib+nb_frame_tf);

	gravi_msg_function_exit(1);
    return (int)cpl_error_get_code();
}


