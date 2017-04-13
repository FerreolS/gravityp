/* $Id: gravity_wavelamp.c,v 1.29 2015/01/10 09:16:12 nazouaoui Exp $
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
#include "gravi_wave.h"
#include "gravi_preproc.h"
#include "gravi_disp.h"

/*-----------------------------------------------------------------------------
                            Private function prototypes
 -----------------------------------------------------------------------------*/

static int gravity_wavelamp_create(cpl_plugin *);
static int gravity_wavelamp_exec(cpl_plugin *);
static int gravity_wavelamp_destroy(cpl_plugin *);
static int gravity_wavelamp(cpl_frameset *, const cpl_parameterlist *);

/*-----------------------------------------------------------------------------
                            Static variables
 -----------------------------------------------------------------------------*/

static char gravity_wavelamp_short[] = "Measure the position of the Argon lines in the spectra.";
static char gravity_wavelamp_description[] =
"This recipe is associated to the template gravity_wavelamp.\n"
"It reduces the raw file obtained with the Argon lamp (WAVELAMP) and process it so that it can be used to calibrate the fiber dispersion (recipe gravity_disp).\n"
    GRAVI_RECIPE_FLOW"\n"
    "* Extract the spectrums of the Argon exposure\n"
    "* Interpolate the spectrums into a common wavelength table\n"
    "* Measure the wavelength position of known Argon lines\n"
    "* Write the product\n"
    GRAVI_RECIPE_INPUT"\n"    
    GRAVI_FLAT_MAP"               : flat calibration (PRO.CATG="GRAVI_FLAT_MAP")\n"
    GRAVI_BAD_MAP"                : badpixel calibration (PRO.CATG="GRAVI_BAD_MAP") \n"
    GRAVI_WAVE_MAP"               : wave calibration (PRO.CATG="GRAVI_WAVE_MAP")\n"
    GRAVI_P2VM_MAP"               : p2vm calibration (PRO.CATG="GRAVI_P2VM_MAP")\n"
    GRAVI_WAVELAMP_RAW"           : long exposure of Argon lamp\n"
    GRAVI_DARK_RAW"               : dark of Argon exposure\n"
    GRAVI_RECIPE_OUTPUT"\n"
    GRAVI_WAVELAMP_MAP"           : spectrum of Argon, with position of lines\n"
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
                    "gravity_wavelamp",
                    gravity_wavelamp_short,
                    gravity_wavelamp_description,
                    "Nabih Azouaoui, Vincent Lapeyrere, JB. Le Bouquin",
                    PACKAGE_BUGREPORT,
                    gravi_get_license(),
                    gravity_wavelamp_create,
                    gravity_wavelamp_exec,
                    gravity_wavelamp_destroy)) {
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
static int gravity_wavelamp_create(cpl_plugin * plugin)
{
    cpl_recipe    * recipe;

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

	return 0;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Execute the plugin instance given by the interface
  @param    plugin  the plugin
  @return   0 if everything is ok
 */
/*----------------------------------------------------------------------------*/
static int gravity_wavelamp_exec(cpl_plugin * plugin)
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
    recipe_status = gravity_wavelamp(recipe->frames, recipe->parameters);

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
static int gravity_wavelamp_destroy(cpl_plugin * plugin)
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
static int gravity_wavelamp(cpl_frameset            * frameset,
		                  const cpl_parameterlist * parlist)
{
	cpl_frameset * wavecalib_frameset = NULL,
	  * badcalib_frameset = NULL, * flatcalib_frameset = NULL,
	  * p2vmcalib_frameset = NULL, * dark_frameset = NULL,
	  * wavelamp_frameset = NULL, * used_frameset = NULL,
	  * darkcalib_frameset = NULL;
	
	cpl_frame * frame;
	
	gravi_data * data = NULL, * dark_map = NULL, * wave_map = NULL,
	  * profile_map = NULL, * badpix_map = NULL, * preproc_data = NULL, * p2vm_map = NULL;
	
	/* Message */
	gravity_print_banner (); 
	cpl_msg_set_time_on();
	cpl_msg_set_component_on();
	gravi_msg_function_start(1);
	
    cpl_ensure_code(gravi_dfs_set_groups(frameset) == CPL_ERROR_NONE,
					cpl_error_get_code()) ;

	/* Identify the input framesets */
	wavelamp_frameset = gravi_frameset_extract_wavelamp_data (frameset);
	dark_frameset      = gravi_frameset_extract_dark_data (frameset);
	darkcalib_frameset = gravi_frameset_extract_dark_map (frameset);
	p2vmcalib_frameset = gravi_frameset_extract_p2vm_map (frameset);
	wavecalib_frameset = gravi_frameset_extract_wave_map (frameset);
	flatcalib_frameset = gravi_frameset_extract_flat_map (frameset);
	badcalib_frameset = gravi_frameset_extract_bad_map (frameset);
	
	/* Check input framesets */
	if ( (cpl_frameset_is_empty (dark_frameset) &&
		  cpl_frameset_is_empty (darkcalib_frameset)) ||
		 cpl_frameset_is_empty (wavelamp_frameset) ||
		 cpl_frameset_is_empty (p2vmcalib_frameset) ||
		 cpl_frameset_is_empty (wavecalib_frameset) ||
		 cpl_frameset_is_empty (flatcalib_frameset)  ||
		 cpl_frameset_is_empty (badcalib_frameset) ) {
	  
	  cpl_error_set_message (cpl_func,  CPL_ERROR_ILLEGAL_INPUT,
							 "Mising DARK, WAVELAMP, P2VM, WAVE, FLAT or BAD") ;
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
						   NULL, frame, "gravity_wavelamp",
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

	CPLCHECK_CLEAN ("Error while loading the calibration map");

	
	/* 
	 * Load input WAVELAMP_RAW 
	 */
	frame = cpl_frameset_get_position (wavelamp_frameset, 0);
	data = gravi_data_load_rawframe (frame, used_frameset);
    gravi_data_detector_cleanup (data, parlist);

	/* Create output data (FIXME: probably not necessary) */
	gravi_data * argon_data = gravi_data_new (0);
    cpl_propertylist * argon_header = gravi_data_get_header (argon_data);
	cpl_propertylist_append (argon_header, gravi_data_get_header (data));

    /* Copy the IMAGING_DATA and IMAGING_DETECTOR extensions */
    gravi_data_copy_ext (argon_data, data, GRAVI_IMAGING_DETECTOR_SC_EXT);
    gravi_data_copy_ext (argon_data, data, GRAVI_IMAGING_DATA_SC_EXT);

	FREE (gravi_data_delete, data);

    /* Collapse ARGON */
	cpl_imagelist * imglist = gravi_data_get_cube (argon_data, GRAVI_IMAGING_DATA_SC_EXT);
	cpl_image * img_median = cpl_imagelist_collapse_median_create (imglist);

    /* Replace data in-place */
	cpl_imagelist_empty (imglist);
    cpl_imagelist_set (imglist, img_median, 0);

    /* Extract spectrum */
    preproc_data = gravi_extract_spectrum (argon_data, profile_map, dark_map,
                                           badpix_map, NULL, parlist);
	FREE (gravi_data_delete, argon_data);
    CPLCHECK_CLEAN ("Cannot extract spectrum");
    
	/* Compute the ARGON WAVE */
    //cpl_table * spectrum_table, * argonwave_table;
    //spectrum_table = gravi_data_get_spectrum_data (preproc_data, GRAVI_SC);
	//argonwave_table = gravi_compute_argon_wave (spectrum_table);
    //gravi_data_add_table (preproc_data, NULL, "WAVE_ARGON_RAW", argonwave_table);
    
    /* Rescale to common wavelength */
    gravi_align_spectrum (preproc_data, wave_map, p2vm_map);
    CPLCHECK_CLEAN ("Cannot re-interpolate spectrum");

	/* Compute the ARGON WAVE */
    //spectrum_table = gravi_data_get_spectrum_data (preproc_data, GRAVI_SC);
	//argonwave_table = gravi_compute_argon_wave (spectrum_table);
    //gravi_data_add_table (preproc_data, NULL, "WAVE_ARGON_RESAMP", argonwave_table);
	
	/* Compute position */
	gravi_compute_argon_pos (preproc_data);

	CPLCHECK_CLEAN ("Cannot compute the positions");

	/* Save the output data file */
	gravi_data_save_new (preproc_data, frameset, NULL, parlist,
						 used_frameset, frame, "gravity_wavelamp",
                         NULL, GRAVI_WAVELAMP_MAP);
	
	CPLCHECK_CLEAN("Could not save the WAVELAMP");

    /* Deallocation of all variables */
	goto cleanup;

cleanup :
	cpl_msg_info(cpl_func,"Memory cleanup");
    FREE (cpl_frameset_delete, wavelamp_frameset);
    FREE (cpl_frameset_delete, wavecalib_frameset);
    FREE (cpl_frameset_delete, badcalib_frameset);
    FREE (cpl_frameset_delete, flatcalib_frameset);
    FREE (cpl_frameset_delete, p2vmcalib_frameset);
    FREE (cpl_frameset_delete, used_frameset);
    FREE (cpl_frameset_delete, dark_frameset);
    FREE (cpl_frameset_delete, darkcalib_frameset);
    FREE (gravi_data_delete, dark_map);
    FREE (gravi_data_delete, wave_map);
    FREE (gravi_data_delete, profile_map);
    FREE (gravi_data_delete, badpix_map);
    FREE (gravi_data_delete, p2vm_map);
    FREE (gravi_data_delete, preproc_data);
	
	gravi_msg_function_exit(1);
    return (int)cpl_error_get_code();
}


