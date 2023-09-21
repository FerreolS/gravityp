/* $Id: gravity_vis_from_p2vmred.c,v 1.29 2011/12/3 09:16:12 nazouaoui Exp $
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
#include "gravi_p2vmred.h"
#include "gravi_eop.h"
#include "gravi_metrology.h"

#include "gravi_signal.h"
#include "gravi_vis.h"
#include "gravi_tf.h"

#include "gravi_preproc.h"

/*-----------------------------------------------------------------------------
                            Private function prototypes
 -----------------------------------------------------------------------------*/

static int gravity_vis_from_p2vmred_create(cpl_plugin *);
static int gravity_vis_from_p2vmred_exec(cpl_plugin *);
static int gravity_vis_from_p2vmred_destroy(cpl_plugin *);
static int gravity_vis_from_p2vmred(cpl_frameset *, cpl_parameterlist *);

/*-----------------------------------------------------------------------------
                            Static variables
 -----------------------------------------------------------------------------*/
static char gravity_vis_from_p2vmred_short[] = "Compute the visibilities from P2VMRED intermediate product.";
static char gravity_vis_from_p2vmred_description[] =
    "This recipe averages the real-time data of P2VMRED files into a VIS product. It allows to run the reduction with different parameters (for instance for SNR thresholding) without having to re-reduce the files from scratch. Typically the reduction is 4x faster when started from this intermediate product.\n"
    "\n"
    "The tag in the DO category can be SINGLE/DUAL and CAL/SCI. They should reflect the mode (SINGLE or DUAL) and the DPR.CATG of the observation (SCIENCE or CALIB). The tag in the PRO.CATG category will be SINGLE/DUAL and CAL/SCI depending on the input tag.\n"
    GRAVI_RECIPE_FLOW"\n"
    "* Load the input file (loop on input files)\n"
    "* Update the selection flag\n"
    "* Average the real-time visibilities\n"
    "* Write the product\n"
    GRAVI_RECIPE_INPUT"\n"    
    GRAVI_P2VMRED_SINGLE_SCIENCE" : Input intermediate product\n"
    GRAVI_RECIPE_OUTPUT"\n"    
    GRAVI_VIS_SINGLE_SCIENCE"       : OIFITS with uncalibrated visibilities\n"
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
                    "gravity_vis_from_p2vmred",
                    gravity_vis_from_p2vmred_short,
                    gravity_vis_from_p2vmred_description,
                    "JB. Le Bouquin, Vincent Lapeyrere, ",
                    PACKAGE_BUGREPORT,
                    gravi_get_license(),
                    gravity_vis_from_p2vmred_create,
                    gravity_vis_from_p2vmred_exec,
                    gravity_vis_from_p2vmred_destroy)) {
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
static int gravity_vis_from_p2vmred_create(cpl_plugin * plugin)
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

    /* PCA visphi flattening */
    gravi_parameter_add_pca (recipe->parameters);

    /* Averaging */
    gravi_parameter_add_average_vis (recipe->parameters);

    /* Snr */
    gravi_parameter_add_compute_snr (recipe->parameters);
    
    /* Rejection */
    gravi_parameter_add_rejection (recipe->parameters, isCalib);

    /* Parameters for gravi_compute_vis */
    gravi_parameter_add_compute_vis (recipe->parameters, isCalib);

    /* Reduce ACQ_CAM */
    p = cpl_parameter_new_value ("gravity.test.reduce-acq-cam", CPL_TYPE_BOOL,
                                 "If TRUE, reduced ACQ_CAM images",
                                 "gravity.test", FALSE);
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "reduce-acq-cam");
    cpl_parameter_disable (p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append (recipe->parameters, p);

    /* --Use existing selection */
	p = cpl_parameter_new_value ("gravity.signal.use-existing-rejection", CPL_TYPE_BOOL,
                                 "Use existing rejection flags (ignore related options)",
                                 "gravity.signal", FALSE);
	cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "use-existing-rejection");
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
static int gravity_vis_from_p2vmred_exec(cpl_plugin * plugin)
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
    recipe_status = gravity_vis_from_p2vmred(recipe->frames, recipe->parameters);

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
static int gravity_vis_from_p2vmred_destroy(cpl_plugin * plugin)
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
static int gravity_vis_from_p2vmred(cpl_frameset * frameset,
					                cpl_parameterlist * parlist)
{
    cpl_frameset * recipe_frameset=NULL, *pcacalib_frameset=NULL, *used_frameset=NULL;
	
	cpl_frame * frame=NULL;
	
	const char * frame_tag=NULL;
	char * proCatg = NULL, * mode=NULL;
	
	gravi_data * p2vmred_data=NULL, * vis_data=NULL, * tmpvis_data=NULL, * pca_calib_data=NULL;
	
	int nb_frame;

	/* Message */
	gravity_print_banner (); 
	cpl_msg_set_time_on();
	cpl_msg_set_component_on();
	gravi_msg_function_start(1);

    cpl_ensure_code(gravi_dfs_set_groups(frameset) == CPL_ERROR_NONE,
					cpl_error_get_code()) ;

    /* Dispatch the frameset */
    recipe_frameset = gravi_frameset_extract_p2vmred_data (frameset);
    pcacalib_frameset = gravi_frameset_extract_pca_calib (frameset);

	/* To use this recipe the frameset must contain a P2VMREDUCED file. */
    if ( cpl_frameset_get_size (recipe_frameset) < 1 ) {
	  cpl_error_set_message (cpl_func, CPL_ERROR_ILLEGAL_INPUT,
							 "Illegal number of P2VMREDUCED file on the frameset");
	  goto cleanup;
    }

    /* Force some options if phase flattening is to be performed */
    if (gravi_param_get_bool (parlist, "gravity.vis.flatten-visphi")) {
        cpl_parameter *phase_ref = cpl_parameterlist_find (parlist, "gravity.vis.phase-ref-sc");
        cpl_parameter *output_phase = cpl_parameterlist_find (parlist, "gravity.vis.output-phase-sc");

        if (strcmp (cpl_parameter_get_string(phase_ref), "SELF_REF") != 0) {
            cpl_msg_warning (cpl_func, "VISPHI flattening requires phase-ref-sc=SELF_REF, forcing");
            cpl_parameter_set_string (phase_ref, "SELF_REF");
        }

        if (strcmp (cpl_parameter_get_string(output_phase), "SELF_VISPHI") != 0) {
            cpl_msg_warning (cpl_func, "VISPHI flattening requires output-phase-sc=SELF_VISPHI, forcing");
            cpl_parameter_set_string (output_phase, "SELF_VISPHI");
        }
    }

	/* Insert calibration frame into the used frameset */
	used_frameset = cpl_frameset_new();

    if ( !cpl_frameset_is_empty (pcacalib_frameset)) {
        frame = cpl_frameset_get_position (pcacalib_frameset, 0);
        pca_calib_data = gravi_data_load_frame (frame, used_frameset);
    } else
      cpl_msg_info (cpl_func, "There is no PHASE_PCA in the frameset");

	/* 
	 * Select the PRO CATG (based on first frame) 
	 */
   
	frame_tag = cpl_frame_get_tag (cpl_frameset_get_position (recipe_frameset, 0));

	if ((strcmp(frame_tag, GRAVI_P2VMRED_DUAL_CALIB) == 0)) {
	  proCatg = cpl_sprintf (GRAVI_VIS_DUAL_CALIB);
	  mode = cpl_sprintf ("gravi_dual");
	}
	else if ((strcmp(frame_tag, GRAVI_P2VMRED_DUAL_SCIENCE) == 0)) {
	  proCatg = cpl_sprintf (GRAVI_VIS_DUAL_SCIENCE);
	  mode = cpl_sprintf ("gravi_dual"); 
	}
	else if ((strcmp(frame_tag, GRAVI_P2VMRED_SINGLE_CALIB) == 0)) {
	  proCatg = cpl_sprintf (GRAVI_VIS_SINGLE_CALIB);
	  mode = cpl_sprintf ("gravi_single");
	}
	else if ((strcmp(frame_tag, GRAVI_P2VMRED_SINGLE_SCIENCE) == 0)) {
	  proCatg = cpl_sprintf (GRAVI_VIS_SINGLE_SCIENCE);
	  mode = cpl_sprintf ("gravi_single");
	}
	else {
	  proCatg = cpl_sprintf ("UNKNOWN");
	  mode = cpl_sprintf ("gravi_single");
	}

	cpl_msg_info (cpl_func,"Mode of the first frame is: %s (will be used for all frames)", mode);

	/*
	 * Loop on input RAW frames to be reduced 
	 */
	
    nb_frame = cpl_frameset_get_size (recipe_frameset);
    for (int ivis = 0; ivis < nb_frame; ivis++){

		cpl_msg_info (cpl_func, " ***** P2VMREDUCED %d over %d ***** ", ivis+1, nb_frame);

		frame = cpl_frameset_get_position (recipe_frameset, ivis);
		p2vmred_data = gravi_data_load_frame (frame, used_frameset);

		/* Compute rejection flags for averaging */
        if (gravi_param_get_bool (parlist, "gravity.signal.use-existing-rejection")) {
            cpl_msg_info (cpl_func,"Don't recompute SNR and selection, use the existing one");
        } else {
	    /* Find outliers */
            gravi_compute_outliers (p2vmred_data, parlist);
            CPLCHECK_MSG ("Cannot compute outliers");
	    
            /* Compute the SNR/GDELAY */
            gravi_compute_snr (p2vmred_data, parlist);
            CPLCHECK_MSG ("Cannot compute SNR");
            
            /* Compute rejection flags for averaging */
            gravi_compute_rejection (p2vmred_data, parlist);
            CPLCHECK_MSG ("Cannot recompute rejection flags signals");
        }

        /* Loop on the wanted sub-integration */
        cpl_size current_frame = 0;
        while (current_frame >= 0)
        {
            
		/* Visibility and flux are averaged and the followings
		 * are saved in Visibility data in tables VIS, VIS2 and T3 */
		tmpvis_data = gravi_compute_vis (p2vmred_data, parlist, &current_frame);
		CPLCHECK_CLEAN ("Cannot average the P2VMRED frames into VIS");

        /* Set the mean TIME and mean MJD if required */
        if (gravi_param_get_bool (parlist, "gravity.vis.force-same-time") ) {
            cpl_msg_info (cpl_func,"Force same time for all quantities/baselines");
            gravi_vis_force_time (tmpvis_data);
            CPLCHECK_CLEAN ("Cannot average the TIME in OI_VIS");
        }
            
        /* Copy the acquisition camera if requested */
        if (current_frame < 0 && gravi_param_get_bool (parlist, "gravity.test.reduce-acq-cam"))
        {
   	    cpl_msg_info (cpl_func, "Copy ACQ into the VIS file");
            gravi_data_copy_ext_insname (tmpvis_data, p2vmred_data, GRAVI_IMAGING_DATA_ACQ_EXT, INSNAME_ACQ);
        }
        
        /* Merge with already existing */
        if (vis_data == NULL) {
            vis_data = tmpvis_data; tmpvis_data = NULL;
        }
        else {
            cpl_msg_info (cpl_func,"Merge with previous OI_VIS");
            gravi_data_append (vis_data, tmpvis_data, 1);
            FREE (gravi_data_delete, tmpvis_data);
        }
        CPLCHECK_CLEAN ("Cannot merge the visibilities");

        }

		cpl_msg_info (cpl_func,"Free the p2vmreduced");
		FREE (gravi_data_delete, p2vmred_data);
    }
	/* End loop on the input files to reduce */

    /* Use the PCA calibration to flatten the VISPHI */
    if (gravi_param_get_bool (parlist, "gravity.vis.flatten-visphi")) {
        cpl_msg_info (cpl_func, "Flatten VISPHI using PCA");
        gravi_flatten_vis(vis_data, pca_calib_data);
        CPLCHECK_CLEAN ("Cannot apply the VISPHI flattening");
    }

    /* Compute QC parameters */
    gravi_compute_vis_qc (vis_data, 0);
    
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
	  
	  cpl_msg_warning (cpl_func, "Average the different observation (if any) = EXPERIMENTAL");
	  gravi_average_vis (vis_data);
	  
	} else {
	  cpl_msg_info (cpl_func, "Don't average the different observation (if any)");
	}

	/* Recompute the TIME column from the MJD column
	 * in all OIFITS tables to follow standard */
	gravi_vis_mjd_to_time (vis_data);

	/* Save the output data file based on the first frame of the frameset */
	frame = cpl_frameset_get_position (recipe_frameset, 0);
	
	gravi_data_save_new (vis_data, frameset, NULL, NULL, parlist,
                         used_frameset, frame, "gravity_vis_from_p2vmred",
                         NULL, proCatg);

	CPLCHECK_CLEAN ("Cannot save the VIS product");

	/* Terminate the function */
	goto cleanup;

cleanup:
	/* Deallocation of all variables */
	cpl_msg_info(cpl_func,"Memory cleanup");
	
	FREE (gravi_data_delete,p2vmred_data);
    FREE (gravi_data_delete, pca_calib_data);
	FREE (gravi_data_delete,vis_data);
	FREE (gravi_data_delete,tmpvis_data);
	FREE (cpl_frameset_delete,recipe_frameset);
	FREE (cpl_frameset_delete,pcacalib_frameset);
	FREE (cpl_frameset_delete,used_frameset);
    FREE (cpl_free,proCatg);
    FREE (cpl_free,mode);
	
	gravi_msg_function_exit(1);
    return (int)cpl_error_get_code();
}
