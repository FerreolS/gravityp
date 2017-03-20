/* $Id: gravity_disp.c,v 1.29 2015/01/10 09:16:12 nazouaoui Exp $
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
 * $Date: 2015/01/10 09:16:12 $
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
#include "gravi_preproc.h"
#include "gravi_p2vmred.h"
#include "gravi_signal.h"
#include "gravi_vis.h"
#include "gravi_metrology.h"
#include "gravi_disp.h"





/*-----------------------------------------------------------------------------
                            Private function prototypes
 -----------------------------------------------------------------------------*/

static int gravity_disp_create(cpl_plugin *);
static int gravity_disp_exec(cpl_plugin *);
static int gravity_disp_destroy(cpl_plugin *);
static int gravity_disp(cpl_frameset *, const cpl_parameterlist *);

/*-----------------------------------------------------------------------------
                            Static variables
 -----------------------------------------------------------------------------*/

static char gravity_disp_short[] = "Calibrate the linearity and the dispersion of the differential delay lines.";
static char gravity_disp_description[] =
    "This recipe is associated to the template GRAVI_all_disp. It measure the phases obtained on the internal source at the position of the Argon lines and various stretch of the FDDL. It deduces the linearity model and the dispersion model of the differential delay lines. These models are stored as polynomials versus wavelength.\n"
    GRAVI_RECIPE_FLOW"\n"
    "* Reduce all the input DISP files\n"
    "* Compute the dispersion parameters\n"
    "* Write product\n"
    GRAVI_RECIPE_INPUT"\n"    
    GRAVI_FLAT_MAP"               : flat calibration (PRO.CATG="GRAVI_FLAT_MAP")\n"
    GRAVI_BAD_MAP"                : badpixel calibration (PRO.CATG="GRAVI_BAD_MAP") \n"
    GRAVI_WAVE_MAP"               : wave calibration (PRO.CATG="GRAVI_WAVE_MAP")\n"
    GRAVI_P2VM_MAP"               : p2vm calibration (PRO.CATG="GRAVI_P2VM_MAP")\n"
    GRAVI_DARK_MAP"               : dark calibration\n"
    GRAVI_DISP_RAW"  (>50)        : raw dispersion\n"
    GRAVI_RECIPE_OUTPUT"\n"
    GRAVI_DISP_VIS"               : intermediate product\n"
    GRAVI_DISP_MODEL"             : dispersion model of FDDL\n"
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
                    "gravity_disp",
                    gravity_disp_short,
                    gravity_disp_description,
                    "Nabih Azouaoui, Vincent Lapeyrere, JB. Le Bouquin",
                    PACKAGE_BUGREPORT,
                    gravi_get_license(),
                    gravity_disp_create,
                    gravity_disp_exec,
                    gravity_disp_destroy)) {
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
static int gravity_disp_create(cpl_plugin * plugin)
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
    gravi_parameter_add_preproc_file (recipe->parameters);
    gravi_parameter_add_p2vmred_file (recipe->parameters);
    gravi_parameter_add_vis_file (recipe->parameters);

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
static int gravity_disp_exec(cpl_plugin * plugin)
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
    recipe_status = gravity_disp(recipe->frames, recipe->parameters);

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
static int gravity_disp_destroy(cpl_plugin * plugin)
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
  @brief    The perpese of the recipe is to reduce the raw calibration file
  	  	    for dispersion calibration
  @param    frameset   the frames list
  @param    parlist    the parameters list
  @return   0 if everything is ok
 */
/*----------------------------------------------------------------------------*/
static int gravity_disp(cpl_frameset            * frameset,
		                  const cpl_parameterlist * parlist)
{
	cpl_frameset * disp_frameset=NULL, * darkcalib_frameset=NULL,
	  * dark_frameset=NULL, * wavecalib_frameset=NULL,
	  * badcalib_frameset=NULL, * flatcalib_frameset=NULL, * dispvis_frameset=NULL,
	  * p2vmcalib_frameset=NULL, * wavelampcalib_frameset=NULL, *used_frameset=NULL,
	  * current_frameset=NULL;

	cpl_frame * frame=NULL;
	
	gravi_data * disp_map = NULL, * data = NULL, * dark_map = NULL, * wave_map = NULL,
	  * profile_map = NULL, * badpix_map = NULL, * p2vmred_data = NULL, * preproc_data = NULL,
        * p2vm_map = NULL, * vis_data = NULL, * tmpvis_data=NULL, * argon_data=NULL;
	
	int nb_frame;

	/* Message */
	gravity_print_banner (); 
	cpl_msg_set_time_on();
	cpl_msg_set_component_on();
	gravi_msg_function_start(1);

	
	/* Identify the RAW and CALIB frames in the input frameset */
    cpl_ensure_code(gravi_dfs_set_groups(frameset) == CPL_ERROR_NONE,
                        cpl_error_get_code()) ;

	/* Check if a DISP_VIS is already existing */
	dispvis_frameset = gravi_frameset_extract_dispvis_data (frameset);

	/* Insert calibration frame into the used frameset */
	used_frameset = cpl_frameset_new();

	/* No DISP_VIS, reduce all data */
	if (cpl_frameset_is_empty (dispvis_frameset)){

		/* Identify the ARGON, P2VM, DISP, DARK, WAVE, FLAT, BADPIX frames */
		wavelampcalib_frameset = gravi_frameset_extract_wavelamp_map (frameset);
		p2vmcalib_frameset = gravi_frameset_extract_p2vm_map (frameset);
    	wavecalib_frameset = gravi_frameset_extract_wave_map (frameset);
    	flatcalib_frameset = gravi_frameset_extract_flat_map (frameset);
    	badcalib_frameset = gravi_frameset_extract_bad_map (frameset);
		
    	darkcalib_frameset = gravi_frameset_extract_dark_map (frameset);
    	dark_frameset = gravi_frameset_extract_dark_data (frameset);
		
        disp_frameset = gravi_frameset_extract_disp_data (frameset);

		/* Check inputs */
        if ( cpl_frameset_is_empty(p2vmcalib_frameset) ||
			 cpl_frameset_is_empty(wavecalib_frameset) ||
			 cpl_frameset_is_empty(flatcalib_frameset) ||
			 cpl_frameset_is_empty(badcalib_frameset) ||
			 ( cpl_frameset_is_empty(dark_frameset) &&
			   cpl_frameset_is_empty(darkcalib_frameset) ) ||
			 cpl_frameset_is_empty(disp_frameset) ||
			 cpl_frameset_is_empty(wavelampcalib_frameset) ) {

            cpl_error_set_message (cpl_func, CPL_ERROR_ILLEGAL_INPUT,
								   "Missing P2VM, DARK, BAD, WAVE, WAVELAMP, FLAT or DISP frames") ;
            goto cleanup;
        }

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
							   NULL, frame, "gravi_single",
                               NULL, GRAVI_DARK_MAP);
		  
		  CPLCHECK_CLEAN ("Could not save the DARK map");
		}
		else if (!cpl_frameset_is_empty (darkcalib_frameset)) {
		  
		  frame = cpl_frameset_get_position (darkcalib_frameset, 0);
		  dark_map = gravi_data_load_frame (frame, used_frameset);
		  
		  CPLCHECK_CLEAN ("Could not load the DARK map");
		}
		else 
		  cpl_msg_error (cpl_func, "There is no DARK in the frame set");
		
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
		
		CPLCHECK_CLEAN ("Error while loading the calibration map");
		
		
		/* 
		 * Loop on input DISP files 
		 */
		
        nb_frame = cpl_frameset_get_size (disp_frameset);
        for (int ivis = 0; ivis < nb_frame; ivis++) {
			current_frameset = cpl_frameset_duplicate (used_frameset);
		  
		    cpl_msg_info (cpl_func, "*** Processing file %d over %d *** ", ivis+1, nb_frame);
			
        	frame = cpl_frameset_get_position (disp_frameset, ivis);
			data = gravi_data_load_rawframe (frame, current_frameset);
            gravi_data_detector_cleanup (data, parlist);
			
            /* Extract spectrum */
            preproc_data = gravi_extract_spectrum (data, profile_map, dark_map,
                                                   badpix_map, NULL, parlist);
            CPLCHECK_CLEAN ("Cannot extract spectrum");

            /* Rescale to common wavelength */
            gravi_align_spectrum (preproc_data, wave_map, p2vm_map);
            CPLCHECK_CLEAN ("Cannot re-interpolate spectrum");

            /* Move extensions from raw_data */
            gravi_data_move_ext (preproc_data, data, GRAVI_ARRAY_GEOMETRY_EXT);
            gravi_data_move_ext (preproc_data, data, GRAVI_OPTICAL_TRAIN_EXT);
            gravi_data_move_ext (preproc_data, data, GRAVI_OPDC_EXT);
            gravi_data_move_ext (preproc_data, data, GRAVI_FDDL_EXT);
            gravi_data_move_ext (preproc_data, data, GRAVI_METROLOGY_EXT);
            CPLCHECK_CLEAN ("Cannot move ext");
        
            FREE (gravi_data_delete,data);

    		/* Option save the proproc file */
			if (gravi_param_get_bool (parlist,"gravity.dfs.preproc-file")) {
			  
    			gravi_data_save_new (preproc_data, frameset, NULL, parlist,
									 current_frameset, frame, "gravity_disp",
                                     NULL, GRAVI_PREPROC);

    			CPLCHECK_CLEAN ("Could not save the preproc data");
    		}

			/* Compute the flux and visibilities for each telescope and
			 * per acquisition with the P2VM applied to preproc_data */
    		p2vmred_data = gravi_compute_p2vmred (preproc_data, p2vm_map, "gravi_single", parlist);
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
			gravi_metrology_reduce (p2vmred_data);
            CPLCHECK_CLEAN ("Cannot reduce metrology");

			/* Compute the SNR_SMT and GDELAY_SMT columns */
			gravi_compute_snr (p2vmred_data, parlist);
			CPLCHECK_MSG ("Cannot compute SNR");
		
			/* Compute the additional signals for averaging */
			gravi_compute_signals (p2vmred_data, disp_map, parlist);
			CPLCHECK_MSG ("Cannot compute signals");

			/* Compute rejection flags for averaging */
			gravi_compute_rejection (p2vmred_data, parlist);
			CPLCHECK_MSG ("Cannot rejection flags signals");
			
			/* Save the P2VMREDUCED */
			if (gravi_param_get_bool (parlist,"gravity.dfs.p2vmred-file")) {
			  
    			gravi_data_save_new (p2vmred_data, frameset, NULL, parlist,
									 current_frameset, frame, "gravity_disp",
                                     NULL, GRAVI_P2VMRED_SINGLE_CALIB);
				
				CPLCHECK_CLEAN ("Cannot save the P2VMREDUCED product");
			}


			/* Visibility and flux are averaged and the followings
			 * are saved in Visibility data in tables VIS, VIS2 and T3 */
			tmpvis_data = gravi_compute_vis (p2vmred_data, parlist);
			CPLCHECK_CLEAN ("Cannot average the visibilities");

			/* Save the VIS */
			if (gravi_param_get_bool (parlist,"gravity.dfs.vis-file")) {
			  
    			gravi_data_save_new (tmpvis_data, frameset, NULL, parlist,
									 current_frameset, frame, "gravity_disp",
                                     NULL, GRAVI_VIS_SINGLE_CALIB);
				
				CPLCHECK_CLEAN ("Cannot save the VIS product");
			}

            /* Merge with already existing */
            if (ivis == 0) {
                vis_data = tmpvis_data; tmpvis_data = NULL;
            }
            else {
                cpl_msg_info (cpl_func,"Merge with previous OI_VIS");
                gravi_data_append (vis_data, tmpvis_data, 1);
                FREE (gravi_data_delete, tmpvis_data);
            }
			CPLCHECK_CLEAN ("Cannot merge the visibilities");			

			cpl_msg_info (cpl_func,"Free the p2vmreduced");
			FREE (cpl_frameset_delete, current_frameset);
    		FREE (gravi_data_delete, p2vmred_data);
        }
		/* End loop on the input files to reduce */

		/* Recompute the TIME column from the MJD column
		 * in all OIFITS tables to follow standard */
		gravi_vis_mjd_to_time (vis_data);

		/* Identify the WAVELAMP in the input frameset */
    	frame = cpl_frameset_get_position (wavelampcalib_frameset, 0);
		argon_data = gravi_data_load_frame (frame, used_frameset);

		/* Duplicate POS_ARGON into the VIS file */
        gravi_data_copy_ext (vis_data, argon_data, "POS_ARGON");
		
		/* Save the output data file based on the first frame of the frameset */
		cpl_frameset_join (used_frameset, disp_frameset);
		frame = cpl_frameset_get_position (disp_frameset, 0);
		
    	gravi_data_save_new (vis_data, frameset, NULL, parlist,
							 used_frameset, frame, "gravity_disp",
							 NULL, GRAVI_DISP_VIS);
		
		CPLCHECK_CLEAN("Could not save the VIS_SINGLE product");

    	FREE (gravi_data_delete, profile_map);
    	FREE (gravi_data_delete, dark_map);
    	FREE (gravi_data_delete, wave_map);
    	FREE (gravi_data_delete, badpix_map);
    	FREE (gravi_data_delete, p2vm_map);
    	FREE (cpl_frameset_delete, darkcalib_frameset);
    	FREE (cpl_frameset_delete, wavecalib_frameset);
    	FREE (cpl_frameset_delete, dark_frameset);
    	FREE (cpl_frameset_delete, flatcalib_frameset);
    	FREE (cpl_frameset_delete, badcalib_frameset);
        FREE (cpl_frameset_delete, p2vmcalib_frameset);
        FREE (cpl_frameset_delete, wavelampcalib_frameset);
    }
    else {

    	/* Load the DIS_VIS already computed */
		frame = cpl_frameset_get_position (dispvis_frameset, 0);
		vis_data = gravi_data_load_frame (frame, used_frameset);
		
		CPLCHECK_CLEAN ("Cannot load the DISP_VIS file");
		
		disp_frameset = cpl_frameset_duplicate (dispvis_frameset);
    }

	/* 
	 * Compute the dispersion table of a set of disp frames 
	 */

    gravi_disp_cleanup (vis_data);

	disp_map = gravi_compute_disp (vis_data);

	CPLCHECK_CLEAN ("Error while computing the disp_map");

	/* Create product frame */
	frame = cpl_frameset_get_position (disp_frameset, 0);

	/* Save the DISP_MODEL */
	gravi_data_save_new (disp_map, frameset, NULL, parlist,
						 used_frameset, frame, "gravity_disp",
						 NULL, GRAVI_DISP_MODEL);
	
	CPLCHECK_CLEAN("Could not save the DISP_MODEL product");

    /* Deallocation of all variables */
	goto cleanup;

cleanup :
	cpl_msg_info(cpl_func,"Memory cleanup");
	FREE (gravi_data_delete, data);
	FREE (gravi_data_delete, tmpvis_data);
	FREE (gravi_data_delete, vis_data);
	FREE (gravi_data_delete, disp_map);
    FREE (gravi_data_delete, dark_map);
    FREE (gravi_data_delete, wave_map);
    FREE (gravi_data_delete, profile_map);
    FREE (gravi_data_delete, badpix_map);
    FREE (gravi_data_delete, p2vmred_data);
    FREE (gravi_data_delete, p2vm_map);
    FREE (gravi_data_delete, preproc_data);
    FREE (cpl_frameset_delete, disp_frameset);
    FREE (cpl_frameset_delete, dispvis_frameset);
    FREE (cpl_frameset_delete, darkcalib_frameset);
    FREE (cpl_frameset_delete, dark_frameset);
    FREE (cpl_frameset_delete, wavecalib_frameset);
    FREE (cpl_frameset_delete, badcalib_frameset);
    FREE (cpl_frameset_delete, flatcalib_frameset);
    FREE (cpl_frameset_delete, p2vmcalib_frameset);
    FREE (cpl_frameset_delete, wavelampcalib_frameset);
	FREE (cpl_frameset_delete, used_frameset);
	FREE (cpl_frameset_delete, current_frameset);

	gravi_msg_function_exit(1);
    return (int)cpl_error_get_code();
}


