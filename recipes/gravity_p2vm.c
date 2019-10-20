/* $Id: gravity_p2vm.c,v 1.29 2009/02/10 09:16:12 llundin Exp $
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
 * $Author: llundin $
 * $Date: 2009/02/10 09:16:12 $
 * $Revision: 1.29 $
 * $Name:  $
 *
 * History :
 *    04/12/2018 use GRAVITY_WAVE.fits calibration file instead of hardcoded values
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
#if defined(__linux__) && defined(__GLIBC__)
#include <malloc.h> //Needed for malloc_trim()
#endif

#include "gravi_data.h"
#include "gravi_pfits.h"
#include "gravi_dfs.h"

#include "gravi_utils.h"

#include "gravi_metrology.h"
#include "gravi_calib.h"
#include "gravi_preproc.h"
#include "gravi_wave.h"
#include "gravi_p2vm.h"

#include "gravi_p2vmred.h"


/*-----------------------------------------------------------------------------
                            Private function prototypes
 -----------------------------------------------------------------------------*/

static int gravity_p2vm_create(cpl_plugin *);
static int gravity_p2vm_exec(cpl_plugin *);
static int gravity_p2vm_destroy(cpl_plugin *);
static int gravity_p2vm(cpl_frameset *, const cpl_parameterlist *);

/*-----------------------------------------------------------------------------
                            Static variables
 -----------------------------------------------------------------------------*/

static char gravity_p2vm_short[] = "Calibrate the instrument bad pixels, wavelength table, interferometric contrast and phase.";
static char gravity_p2vm_description[] =
"This recipe reduces the internal calibrations. As a special sequence of shutter opening is required, it is advised to always build the SOF with a complete sequence of files obtained within a single execution of the p2vm calibration template. However it is still possible to input a SOF with DARK_RAW only, or DARK_RAW and FLAT_RAW only. It is also possible to input a SOF with some already processed calibration (e.g WAVE).\n"
    GRAVI_RECIPE_FLOW"\n"
    "* Compute the dark, write product\n"
    "* Compute the flat, write product\n"
    "* Compute the badpixels, write product\n"
    "* Compute the spectral calibration, write product\n"
    "* Compute the p2vm, write product\n"
    GRAVI_RECIPE_INPUT"\n"    
    GRAVI_DARK_RAW"      : raw dark, all shutters closed (DPR.TYPE=DARK)\n"
    GRAVI_FLAT_RAW"  x4  : raw flats, one shutter open (DPR.TYPE=FLAT)\n"
    GRAVI_P2VM_RAW"  x6  : raw p2vms, two shutters open (DPR.TYPE=P2VM)\n"
    GRAVI_WAVE_RAW"      : raw wavelength calibration for FT (DPR.TYPE=WAVE)\n"
    GRAVI_WAVESC_RAW"    : raw wavelength calibration for SC (DPR.TYPE=WAVE,SC)\n"
    GRAVI_RECIPE_OUTPUT"\n"
    GRAVI_DARK_MAP"          : dark calibration\n"
    GRAVI_FLAT_MAP"          : flat calibration\n"
    GRAVI_BAD_MAP"           : badpixel calibration\n"
    GRAVI_WAVE_MAP"          : wave calibration\n"
    GRAVI_P2VM_MAP"          : p2vm calibration\n"
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
                    "gravity_p2vm",
                    gravity_p2vm_short,
                    gravity_p2vm_description,
                    "Nabih Azouaoui, Vincent Lapeyrere, JB. Le Bouquin",
                    PACKAGE_BUGREPORT,
                    gravi_get_license(),
                    gravity_p2vm_create,
                    gravity_p2vm_exec,
                    gravity_p2vm_destroy)) {
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
static int gravity_p2vm_create(cpl_plugin * plugin)
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

    /* Use static names (output_procatg.fits) */
    gravi_parameter_add_static_name (recipe->parameters);

    /* Debug file */
    gravi_parameter_add_debug_file (recipe->parameters);
    gravi_parameter_add_preproc_file (recipe->parameters);

    /* Bias-method */
    gravi_parameter_add_biasmethod (recipe->parameters);

    /* Badpix and profile */
    gravi_parameter_add_badpix (recipe->parameters);
    gravi_parameter_add_profile (recipe->parameters);
    //gravi_parameter_add_preproc (recipe->parameters);

    /* Wave option */
    gravi_parameter_add_wave (recipe->parameters);

    /* Extraction */
    gravi_parameter_add_extract (recipe->parameters);
    
    /* Phase definition in P2VM */
    p = cpl_parameter_new_enum ("gravity.calib.phase-calibration", CPL_TYPE_STRING,
                                "This option changes the phase reference of the P2VM:\n "
                                "NONE defines phiA(lbd) at zero for all baselines "
                                "(P2VM calibrates only the internal phase-shift of the beam combiner);\n "
                                "CLOSURE defines phiA(lbd) at zero for baselines 01, 02 and 03 "
                                "(P2VM calibrates the phase-shift and the closure-phase of the beam combiner);\n "
                                "DISP defines phiA(lbd) to have zero mean and minimum GD for baselines (01,02,03); "
                                "(P2VM calibrates the phase-shift, the closure-phase and the "
                                "spectral-dispersion of the beam combiner);\n "
                                "FULL defines phiA(lbd) to have zero-GD for baselines (01,02,03)",
                                "(P2VM calibrates the full phase with respect to zero-astrometry);\n "
                                "gravi.p2vm", "FULL", 4, "NONE", "CLOSURE", "DISP", "FULL");
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "phase-calibration");
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
static int gravity_p2vm_exec(cpl_plugin * plugin)
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
    recipe_status = gravity_p2vm(recipe->frames, recipe->parameters);
                                                                           
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
static int gravity_p2vm_destroy(cpl_plugin * plugin)
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
static int gravity_p2vm(cpl_frameset            * frameset,
                        const cpl_parameterlist * parlist)
{
    cpl_frameset * p2vm_frameset=NULL, * wavecalib_frameset=NULL, * darkcalib_frameset=NULL,
            * flatcalib_frameset=NULL, * dark_frameset=NULL, * wave_frameset=NULL, * wavesc_frameset=NULL,
            * badcalib_frameset=NULL, * flat_frameset=NULL, * used_frameset=NULL, * current_frameset=NULL,
			* wave_param_frameset=NULL;

    cpl_frame * frame=NULL, * frame_p2vm=NULL;

    gravi_data * p2vm_map=NULL, * data=NULL, * dark_map=NULL, * wave_map=NULL,
            * profile_map=NULL, * badpix_map=NULL, * wave_param=NULL;
    gravi_data * spectrum_data=NULL;
    gravi_data * preproc_data=NULL;
    gravi_data ** raw_data=NULL;
    gravi_data * p2vmred_data = NULL;

    int nb_frame, nb_frame_gain = 0;

    cpl_propertylist * met_plist=NULL;
    int ** valid_trans = cpl_malloc (2 * sizeof (int*));
    int ** valid_CP = cpl_malloc (2 * sizeof (int*));
    char ext_regexp[500];
    clock_t start;

    for (int i = 0 ; i < 2; i++){
        valid_trans[i] = cpl_calloc (4, sizeof (int));
        valid_CP[i] = cpl_calloc (6, sizeof (int));
    }

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
    darkcalib_frameset = gravi_frameset_extract_dark_map (frameset);
    
    /* Extract BAD frameset */
    badcalib_frameset = gravi_frameset_extract_bad_map (frameset);
    if ( cpl_frameset_is_empty (badcalib_frameset) && cpl_frameset_is_empty (dark_frameset) &&
          cpl_frameset_is_empty (darkcalib_frameset) ) {
        cpl_error_set_message (cpl_func, CPL_ERROR_ILLEGAL_INPUT,"Missing DARK or BAD on the frameset");
        goto cleanup;
    }

    /* Extract FLAT frameset */
    flat_frameset = gravi_frameset_extract_flat_data (frameset);
    flatcalib_frameset = gravi_frameset_extract_flat_map (frameset);
    if ( cpl_frameset_is_empty (flat_frameset) &&
              cpl_frameset_is_empty (flatcalib_frameset) ) {
        cpl_error_set_message (cpl_func, CPL_ERROR_ILLEGAL_INPUT,"Missing FLAT on the frameset");
        goto cleanup;
    }

    /* Extract WAVE frameset */
    wave_frameset = gravi_frameset_extract_wave_data (frameset);
    wavesc_frameset = gravi_frameset_extract_wavesc_data (frameset);
    wavecalib_frameset = gravi_frameset_extract_wave_map (frameset);
    if ( cpl_frameset_is_empty (wave_frameset) &&
              cpl_frameset_is_empty (wavecalib_frameset) ) {
        cpl_error_set_message (cpl_func, CPL_ERROR_ILLEGAL_INPUT,"Missing WAVE on the frameset");
        goto cleanup;
    }

    /* Extract P2VM frameset */
    p2vm_frameset = gravi_frameset_extract_p2vm_data (frameset);
    if ( cpl_frameset_get_size (p2vm_frameset) !=6) {
        cpl_error_set_message (cpl_func, CPL_ERROR_ILLEGAL_INPUT,"Illegal number of P2VM on the frameset");
        goto cleanup;
    }
        
    /* Extract calibration file wave_param frameset */
    wave_param_frameset = gravi_frameset_extract_wave_param (frameset);
        if ( cpl_frameset_is_empty (wave_param_frameset) ) {
        cpl_error_set_message (cpl_func, CPL_ERROR_ILLEGAL_INPUT,"Missing WAVE_PARAM static calibration file");
        goto cleanup;
    }

    /*
     * (1) Identify and extract the dark file
     */

    if (!cpl_frameset_is_empty (dark_frameset)) {
        cpl_msg_info (cpl_func, " ***** Compute DARK map ***** ");

        /* Load this DARK_RAW */
        frame = cpl_frameset_get_position (dark_frameset, 0);
        data = gravi_data_load_rawframe (frame, used_frameset);
        gravi_data_detector_cleanup (data, parlist);

        /* Compute the dark */
        dark_map = gravi_compute_dark (data);
        FREE (gravi_data_delete, data);

        CPLCHECK_CLEAN ("Cannot compute the DARK map");

        /* Save the dark map */
        gravi_data_save_new (dark_map, frameset, NULL, NULL, parlist,
                             NULL, frame, "gravity_p2vm",
                             NULL, GRAVI_DARK_MAP);

        CPLCHECK_CLEAN ("Could not save the DARK map");	  
    }
    else if (!cpl_frameset_is_empty (darkcalib_frameset)) {
        cpl_msg_info (cpl_func, " ***** Get DARK map ***** ");

        /* Load this DARK */
        frame = cpl_frameset_get_position (darkcalib_frameset, 0);
        dark_map = gravi_data_load_frame (frame, used_frameset);

        CPLCHECK_CLEAN ("Could not load the DARK map");

        /* This new DARK may be used latter for BADPIX */
        cpl_frameset_insert (dark_frameset, cpl_frame_duplicate (frame));
    }


    /*
     * (2) Identify and extract the BADPIX file
     */

    if (!cpl_frameset_is_empty (badcalib_frameset)) {
        cpl_msg_info (cpl_func, " ***** Get BAD pixel map ***** ");

        /* The BADPIX is already in the frameset */
        frame = cpl_frameset_get_position (badcalib_frameset, 0);
        badpix_map = gravi_data_load_frame (frame, used_frameset);
    }
    else if (!cpl_frameset_is_empty (dark_frameset) &&
            !cpl_frameset_is_empty (flat_frameset)) {
        cpl_msg_info (cpl_func, " ***** Compute BAD pixel map from DARK and FLAT_RAW ***** ");

        /* Identify the flat files */
        nb_frame_gain = cpl_frameset_get_size (flat_frameset);
        raw_data = cpl_calloc (nb_frame_gain, sizeof(gravi_data *));

        /* Build the list of FLAT files and output file name */
        for (int i = 0; i < nb_frame_gain; i++) {
            frame = cpl_frameset_get_position (flat_frameset, i);
            raw_data[i] = gravi_data_load_rawframe (frame, NULL);
            gravi_data_detector_cleanup (raw_data[i], parlist);
        }

        /* Compute it from the DARK and the FLAT_RAW */
        badpix_map = gravi_compute_badpix (dark_map, raw_data,
                nb_frame_gain, parlist);

        CPLCHECK_CLEAN("Cannot compute the BAD pixel from DARK and FLAT");

        /* Save the BADPIX */
        frame = cpl_frameset_get_position (dark_frameset, 0);
        gravi_data_save_new (badpix_map, frameset, NULL, NULL, parlist,
                             NULL, frame, "gravity_p2vm",
                             NULL, GRAVI_BAD_MAP);

        CPLCHECK_CLEAN ("Could not save the BAD pixel map");

        FREELOOP (gravi_data_delete, raw_data, nb_frame_gain);
    }

    /*
     * (3) Check that the FLAT map in the input frameset
     */

    if (!cpl_frameset_is_empty(flatcalib_frameset)) {
        cpl_msg_info (cpl_func, " ***** Get FLAT map ***** ");

        /* The FLAT is already in the frameset */
        frame = cpl_frameset_get_position (flatcalib_frameset, 0);
        profile_map = gravi_data_load_frame (frame, used_frameset);
    }
    else if (!cpl_frameset_is_empty (flat_frameset)) {
        cpl_msg_info (cpl_func, " ***** Compute FLAT map ***** ");


        /* Check calibs */
        if (!badpix_map || !dark_map) {
            ERROR_CLEAN (CPL_ERROR_ILLEGAL_INPUT, "Missing DARK or BAD in frameset");
        }

        /* Identify the flat files */
        nb_frame_gain = cpl_frameset_get_size (flat_frameset);
        raw_data = cpl_malloc (nb_frame_gain * sizeof(gravi_data *));

        /* Build the list of FLAT files and output file name */
        for (int i = 0; i < nb_frame_gain; i++) {
            frame = cpl_frameset_get_position (flat_frameset, i);
            raw_data[i] = gravi_data_load_rawframe (frame, used_frameset);
            gravi_data_detector_cleanup (raw_data[i], parlist);
        }


        /* Compute the profile from the list of FLAT */
        profile_map = gravi_compute_profile (raw_data, dark_map, badpix_map, nb_frame_gain, parlist);
        CPLCHECK_CLEAN ("Cannot compute the FLAT profile");

        /* Compute the gain QC */
        cpl_propertylist * gain_header;
        gain_header = gravi_compute_gain (raw_data, nb_frame_gain, dark_map);
        CPLCHECK_CLEAN ("Cannot compute the GAIN");

        /* Put the gain QC in profile_map */
        cpl_propertylist * profile_header = gravi_data_get_header (profile_map);
        cpl_propertylist_append (profile_header, gain_header);
        FREE (cpl_propertylist_delete, gain_header);

        /* Save the FLAT map */
        frame = cpl_frameset_get_position (flat_frameset, 0);

        gravi_data_save_new (profile_map, frameset, NULL, NULL, parlist,
                             used_frameset, frame, "gravity_p2vm",
                             NULL, GRAVI_FLAT_MAP);

        CPLCHECK_CLEAN ("Could not save the FLAT profile_map");

        /* Free the list of files */
        FREELOOP (gravi_data_delete, raw_data, nb_frame_gain);
    }

    /*
     * (4) Check and get the WAVE frame
     */

    if (!cpl_frameset_is_empty (wave_param_frameset)) {
        frame = cpl_frameset_get_position (wave_param_frameset, 0);
        wave_param = gravi_data_load_frame (frame, used_frameset);
    }
    else
        cpl_msg_error (cpl_func, "There is no WAVE_PARAM in the frameset (did you forget the static calibration files?)");

    if (!cpl_frameset_is_empty (wavecalib_frameset)) {
        cpl_msg_info (cpl_func, " ***** Get wave map ***** ");

        /* The WAVE in the frame_set is a calibrated wave */
        frame = cpl_frameset_get_position (wavecalib_frameset, 0);
        wave_map = gravi_data_load_frame (frame, used_frameset);
    }
    else if (!cpl_frameset_is_empty (wave_frameset)) {
        cpl_msg_info (cpl_func, " ***** Process the WAVE_RAW ***** ");

        /* Check calibs */
        if ( !badpix_map || !profile_map || !dark_map) {
            ERROR_CLEAN (CPL_ERROR_ILLEGAL_INPUT, "Missing DARK, FLAT, or BAD in frameset");
        }

        /* Create the WAVE product */
        wave_map = gravi_data_new (0);

        /* Get the frame */
        frame = cpl_frameset_get_position (wave_frameset, 0);

        /* Load WAVE_RAW SC */
        snprintf(ext_regexp, 499, "^(%s|%s)$", GRAVI_IMAGING_DATA_SC_EXT, GRAVI_IMAGING_DETECTOR_SC_EXT);
        gravi_data * wave_sc_data = gravi_data_load_rawframe_ext (frame, used_frameset, ext_regexp);
        gravi_data_detector_cleanup (wave_sc_data, parlist);

        /* Reduce WAVE_RAW SC */
        cpl_msg_info (cpl_func, "Extract SC SPECTRUM for WAVE_RAW");
        spectrum_data = gravi_extract_spectrum (wave_sc_data, profile_map, dark_map,
                badpix_map, NULL, parlist, GRAVI_DET_SC);
        FREE (gravi_data_delete, wave_sc_data);

        /* Load WAVE_RAW FT */
        snprintf(ext_regexp, 499, "^(%s|%s)$", GRAVI_IMAGING_DATA_FT_EXT, GRAVI_IMAGING_DETECTOR_FT_EXT);
        gravi_data * wave_ft_data = gravi_data_load_rawframe_ext (frame, NULL, ext_regexp);

        /* Reduce WAVE_RAW FT */
        cpl_msg_info (cpl_func, "Extract FT SPECTRUM for WAVE_RAW");
        gravi_data * ft_spectrum_data = gravi_extract_spectrum (wave_ft_data, profile_map, dark_map,
                badpix_map, NULL, parlist, GRAVI_DET_FT);
        FREE (gravi_data_delete, wave_ft_data);        

        /* Copy FT extensions to spectrum_data */
        gravi_data_move_ext(spectrum_data, ft_spectrum_data, GRAVI_IMAGING_DETECTOR_FT_EXT);
        gravi_data_move_ext(spectrum_data, ft_spectrum_data, GRAVI_SPECTRUM_DATA_FT_EXT);
        cpl_propertylist_copy_property_regexp(gravi_data_get_header(spectrum_data), 
                                              gravi_data_get_header(ft_spectrum_data), 
                                              "^ESO QC ",0);
        FREE (gravi_data_delete, ft_spectrum_data);        

        /* Load WAVE_RAW metrology */
        gravi_data * wave_met_data  = gravi_data_load_rawframe_ext (frame, NULL, GRAVI_METROLOGY_EXT);

        /* Compute the P2VM of the MET from the WAVE_RAW */
        cpl_msg_info (cpl_func, "Compute the P2VM_MET from WAVE_RAW");
        cpl_table * met_table = gravi_data_get_table (wave_met_data, GRAVI_METROLOGY_EXT);
        double lbd_met = gravi_pfits_get_met_wavelength_mean (gravi_data_get_header(wave_met_data), met_table);
        cpl_table * p2vm_met = gravi_metrology_compute_p2vm (met_table, lbd_met);

        /* Set the P2VM_MET in WAVE */
        gravi_data_add_table (wave_map, NULL, GRAVI_P2VM_MET_EXT, p2vm_met);

        CPLCHECK_CLEAN ("Cannot compute P2VM_MET");

        /* Computing OPDs */
        cpl_msg_info (cpl_func, "Compute OPDs for WAVE_RAW");
        gravi_wave_compute_opds (spectrum_data, 
                                 gravi_data_get_table (wave_met_data, GRAVI_METROLOGY_EXT),
                                 GRAVI_DET_ALL);
        FREE (gravi_data_delete, wave_met_data);        

        CPLCHECK_CLEAN ("Cannot compute OPDs");

        /* Save the file with OPD_SC, OPD_FT, OI_VIS_MET */
        if (gravi_param_get_bool (parlist,"gravity.dfs.debug-file")) {
            gravi_data_save_new (spectrum_data, frameset, NULL, NULL, parlist,
                    used_frameset, frame, "gravity_p2vm",
                    NULL, "DEBUG");
        }

        /* Compute wave calibration for FT */
        gravi_compute_wave (wave_map, spectrum_data, GRAVI_FT, parlist, wave_param);

        CPLCHECK_CLEAN ("Cannot compute wave for FT");

        /* Load and process the WAVESC_RAW (if any) */
        if (!cpl_frameset_is_empty (wavesc_frameset)) {
            /* Erase spectrum data */
            cpl_msg_info (cpl_func, "Delete WAVE_RAW data");
            FREE (gravi_data_delete, spectrum_data);

            cpl_msg_info (cpl_func, " ***** Process the WAVESC_RAW ***** ");

            /* Get the frame */
            frame = cpl_frameset_get_position (wavesc_frameset, 0);

            /* Load WAVESC_RAW SC*/
            snprintf(ext_regexp, 499, "^(%s|%s)$", GRAVI_IMAGING_DATA_SC_EXT, GRAVI_IMAGING_DETECTOR_SC_EXT);
            gravi_data * wavesc_sc_data  = gravi_data_load_rawframe_ext (frame, used_frameset, ext_regexp);
            gravi_data_detector_cleanup (wavesc_sc_data, parlist);

            cpl_msg_info (cpl_func, "Extract SC SPECTRUM for WAVESC_RAW");
            spectrum_data = gravi_extract_spectrum (wavesc_sc_data, profile_map, dark_map,
                    badpix_map, NULL, parlist, GRAVI_DET_SC);
            CPLCHECK_CLEAN ("Cannot extract SC spectrum from WAVESC_RAW");
            FREE (gravi_data_delete, wavesc_sc_data);

            /* Load WAVESC_RAW FT*/
            snprintf(ext_regexp, 499, "^(%s|%s)$", GRAVI_IMAGING_DATA_FT_EXT, GRAVI_IMAGING_DETECTOR_FT_EXT);
            gravi_data * wavesc_ft_data  = gravi_data_load_rawframe_ext (frame, used_frameset, ext_regexp);

            cpl_msg_info (cpl_func, "Extract FT SPECTRUM for WAVESC_RAW");
            gravi_data * ft_spectrum_data = gravi_extract_spectrum (wavesc_ft_data, profile_map, dark_map,
                    badpix_map, NULL, parlist, GRAVI_DET_FT);
            CPLCHECK_CLEAN ("Cannot extract FT spectrum from WAVESC_RAW");
            FREE (gravi_data_delete, wavesc_ft_data);

            /* Copy FT extensions to spectrum_data */
            gravi_data_move_ext(spectrum_data, ft_spectrum_data, GRAVI_IMAGING_DETECTOR_FT_EXT);
            gravi_data_move_ext(spectrum_data, ft_spectrum_data, GRAVI_SPECTRUM_DATA_FT_EXT);
            cpl_propertylist_copy_property_regexp(gravi_data_get_header(spectrum_data), 
                                                  gravi_data_get_header(ft_spectrum_data), 
                                                  "^ESO QC ",0);
            FREE (gravi_data_delete, ft_spectrum_data);        

            /* Load WAVESC_RAW metrology */
            gravi_data * wavesc_met_data  = gravi_data_load_rawframe_ext (frame, NULL, GRAVI_METROLOGY_EXT);

            cpl_msg_info (cpl_func, "Compute OPDs for WAVESC_RAW");
            gravi_wave_compute_opds (spectrum_data, 
                                     gravi_data_get_table (wavesc_met_data, GRAVI_METROLOGY_EXT),
                                     GRAVI_DET_SC);
            FREE (gravi_data_delete, wavesc_met_data);

            CPLCHECK_CLEAN ("Cannot process the WAVESC_RAW");

            /* Save the file with OPD_SC, OPD_FT, OI_VIS_MET */
            if (gravi_param_get_bool (parlist,"gravity.dfs.debug-file")) {
                gravi_data_save_new (spectrum_data, frameset, NULL, NULL,
                        parlist, used_frameset, frame, "gravity_p2vm",
                        NULL, "DEBUG");
            }
        }
        else {
            cpl_msg_warning (cpl_func, "No WAVESC_RAW in the SOF:"
                    " SC wavelength will be inaccurate");
        }

        /* Compute wave calibration for SC */
        gravi_compute_wave (wave_map, spectrum_data, GRAVI_SC, parlist, wave_param);

        CPLCHECK_CLEAN ("Cannot compute wave for SC");

        /* Free the  spectrum */
        FREE (gravi_data_delete, spectrum_data);

        /* Compute the QC */
        gravi_wave_qc (wave_map, profile_map);

        /* Save the WAVE map */
        gravi_data_save_new (wave_map, frameset, NULL, NULL, parlist,
                             used_frameset, frame, "gravity_p2vm",
                             NULL, GRAVI_WAVE_MAP);

        CPLCHECK_CLEAN ("Could not save the WAVE map");
    }

    /* 
     * Build the frameset for the P2VM
     */

    /* Check if P2VM are provided */
    if ( cpl_frameset_is_empty (p2vm_frameset) ) {
        cpl_msg_info (cpl_func,"All RAW data reduced... stop recipe");
        goto cleanup;
    }

    /* Check calibs */
    if ( !badpix_map || !profile_map || !wave_map || !dark_map) {
        ERROR_CLEAN (CPL_ERROR_ILLEGAL_INPUT, "Missing DARK, FLAT, BAD, or WAVE in frameset");
    }

    /* Add the FLAT_RAW and WAVE_RAW to the p2vm frameset */
    if ( !cpl_frameset_is_empty (flat_frameset) )
        cpl_frameset_join (p2vm_frameset, flat_frameset);

    if ( !cpl_frameset_is_empty (wave_frameset) )
        cpl_frameset_join (p2vm_frameset, wave_frameset);

    /* Get the number of the p2vm frame contained in the frameset */
    nb_frame = cpl_frameset_get_size (p2vm_frameset);

    /* Check if the 11 files are here */
    if (nb_frame != 11) {
        ERROR_CLEAN (CPL_ERROR_ILLEGAL_INPUT, "Missing P2VM in frameset");
    }

    /*
     * (6) Loop on files of the p2vm frameset
     */

    /* Construction of the p2vm data. */
    /* START EKW 04/12/2018 read wave parameter from calibration file - Load the WAVE_PARAM Parameter */

    cpl_msg_info (cpl_func, " ***** Create the P2VM ***** ");
    p2vm_map = gravi_create_p2vm (wave_map,wave_param);
    CPLCHECK_CLEAN ("Cannot create the P2VM data");

    /* Loop on files */
    int i_wave;
    for (int i = 0; i < nb_frame; i++) {

        cpl_msg_info (cpl_func, " ***** file %d over %d ***** ", i+1, nb_frame );
        current_frameset = cpl_frameset_duplicate (used_frameset);

        /* Load this frame */
        frame = cpl_frameset_get_position (p2vm_frameset, i);
        cpl_frameset_insert (used_frameset, cpl_frame_duplicate (frame));
        snprintf(ext_regexp, 499, "^(%s|%s)$", GRAVI_ARRAY_GEOMETRY_EXT, GRAVI_OPTICAL_TRAIN_EXT);
        gravi_data * hdr_data = gravi_data_load_rawframe_ext (frame, current_frameset, ext_regexp);
        CPLCHECK_CLEAN ("Cannot load data");

        /* Verbose the shutters */
        cpl_msg_info (cpl_func, "Shutters: %d-%d-%d-%d",
                      gravi_data_get_shutter (hdr_data, 0), gravi_data_get_shutter (hdr_data, 1),
                      gravi_data_get_shutter (hdr_data, 2), gravi_data_get_shutter (hdr_data, 3));

        /* Create the product filename as the first P2VM file (1-1-0-0).
         * This is to ensure the run_gravi_reduce.py script check
         * for the correct product */
        if ( frame_p2vm==NULL && gravi_data_check_shutter (hdr_data, 1,1,0,0) )
        {
            cpl_msg_info (cpl_func,"Use this frame for the P2VM product");
            frame_p2vm = frame;
        }

        /*
         * If all shutter open we just continue
         */
        if ( gravi_data_check_shutter (hdr_data, 1,1,1,1) ) {

            /* Here we go to next file */
            cpl_msg_info (cpl_func, "Nothing to be done with the file... yet.");
            i_wave = i;

            FREE (gravi_data_delete, hdr_data);        
            FREE (cpl_frameset_delete, current_frameset);
            continue;
        }
        FREE (gravi_data_delete, hdr_data);        
        /* End if all shutters open */

        /* Load SC */
        snprintf(ext_regexp, 499, "^(%s|%s|%s|%s)$", GRAVI_ARRAY_GEOMETRY_EXT, GRAVI_OPTICAL_TRAIN_EXT, GRAVI_IMAGING_DATA_SC_EXT, GRAVI_IMAGING_DETECTOR_SC_EXT);
        gravi_data * sc_data = gravi_data_load_rawframe_ext (frame, current_frameset, ext_regexp);
        gravi_data_detector_cleanup (sc_data, parlist);

        /* Extract SC spectrum */
        preproc_data = gravi_extract_spectrum (sc_data, profile_map, dark_map,
                badpix_map, NULL, parlist, GRAVI_DET_SC);
        CPLCHECK_CLEAN ("Cannot extract spectrum");
        FREE (gravi_data_delete, sc_data);

        /* Rescale SC to common wavelength */
        gravi_align_spectrum (preproc_data, wave_map, p2vm_map, GRAVI_DET_SC, parlist);
        CPLCHECK_CLEAN ("Cannot re-interpolate spectrum");

        /* Compute the part of the p2vm associated to this file */
        gravi_compute_p2vm (p2vm_map, preproc_data, valid_trans, valid_CP, GRAVI_DET_SC);
        CPLCHECK_CLEAN("Cannot compute the P2VM");

        /* Load FT */
        snprintf(ext_regexp, 499, "^(%s|%s)$", GRAVI_IMAGING_DATA_FT_EXT, GRAVI_IMAGING_DETECTOR_FT_EXT);
        gravi_data * ft_data = gravi_data_load_rawframe_ext (frame, NULL, ext_regexp);

        /* Extract FT spectrum */
        gravi_data * ft_preproc_data = gravi_extract_spectrum (ft_data, profile_map, dark_map,
                badpix_map, NULL, parlist, GRAVI_DET_FT);
        FREE (gravi_data_delete, ft_data);
        CPLCHECK_CLEAN ("Cannot extract spectrum");

        /* Rescale FT to common wavelength */
        gravi_align_spectrum (ft_preproc_data, wave_map, p2vm_map, GRAVI_DET_FT, parlist);
        CPLCHECK_CLEAN ("Cannot re-interpolate spectrum");

        /* Compute the part of the p2vm associated to this file */
        gravi_compute_p2vm (p2vm_map, ft_preproc_data, valid_trans, valid_CP, GRAVI_DET_FT);
        CPLCHECK_CLEAN("Cannot compute the P2VM");

        /* Copy FT extensions to spectrum_data */
        gravi_data_move_ext(preproc_data, ft_preproc_data, GRAVI_IMAGING_DETECTOR_FT_EXT);
        gravi_data_move_ext(preproc_data, ft_preproc_data, GRAVI_SPECTRUM_DATA_FT_EXT);
        FREE (gravi_data_delete, ft_preproc_data);        

        /* Option save the preproc file */
        if (gravi_param_get_bool (parlist,"gravity.dfs.preproc-file")) {

            gravi_data_save_new (preproc_data, frameset, NULL, NULL, parlist,
                    current_frameset, frame, "gravity_p2vm",
                    NULL, GRAVI_PREPROC);
        }

        /* Delete the preproc data */
        start = clock();
        FREE (gravi_data_delete,preproc_data);
        cpl_msg_info(cpl_func, "Execution time to delete preproc_data : %f s",
                     (clock() - start) / (double)CLOCKS_PER_SEC);

        /* End loop on files to build the p2vm */
        FREE (cpl_frameset_delete, current_frameset);
    }


    /* 
     * (7) P2VM normalization
     */

    gravi_p2vm_normalisation (p2vm_map, valid_trans, valid_CP);
    CPLCHECK_CLEAN("Cannot normalise the p2vm_map");


    /* 
     * (8) Analyse the WAVE to get the phase correction 
     *     and the internal spectrum to latter correct.
     */

    cpl_msg_info (cpl_func, " ***** Analyse the WAVE file to calibrate the internal closure and transmission ***** ");

    /* Get the WAVE_SC frame (or WAVE if WAVE_SC does not exist)*/
    if (!cpl_frameset_is_empty (wavesc_frameset)) {
        /* TODO: using wavesc instead of wave, but does not work for now (wrapping issue?) */
        frame = cpl_frameset_get_position (wave_frameset, 0);
    } else {
        frame = cpl_frameset_get_position (wave_frameset, 0);
    }

    /* Load SC WAVE */
    snprintf(ext_regexp, 499, "^(%s|%s|%s|%s)$", GRAVI_ARRAY_GEOMETRY_EXT, GRAVI_OPTICAL_TRAIN_EXT, GRAVI_IMAGING_DATA_SC_EXT, GRAVI_IMAGING_DETECTOR_SC_EXT);
    gravi_data * wave_sc_data = gravi_data_load_rawframe_ext (frame, NULL, ext_regexp);
    gravi_data_detector_cleanup (wave_sc_data, parlist);

    /* Extract SC spectrum */
    preproc_data = gravi_extract_spectrum (wave_sc_data, profile_map, dark_map,
            badpix_map, NULL, parlist, GRAVI_DET_SC);
    CPLCHECK_CLEAN ("Cannot extract SC spectrum");

    /* Move extensions necessary to compute the p2vmred */
    gravi_data_move_ext (preproc_data, wave_sc_data, GRAVI_ARRAY_GEOMETRY_EXT);
    gravi_data_move_ext (preproc_data, wave_sc_data, GRAVI_OPTICAL_TRAIN_EXT);
    FREE (gravi_data_delete, wave_sc_data);
    CPLCHECK_CLEAN ("Cannot move ext");

    /* Rescale to common wavelength */
    gravi_align_spectrum (preproc_data, wave_map, p2vm_map, GRAVI_DET_SC, parlist);
    CPLCHECK_CLEAN ("Cannot re-interpolate SC spectrum");

    /* Compute P2VMRED */
    p2vmred_data = gravi_compute_p2vmred(preproc_data, p2vm_map, "gravi_single", 
                                         parlist, GRAVI_DET_SC);
    FREE (gravi_data_delete, preproc_data);
    CPLCHECK_CLEAN ("Cannot apply p2vm");

    /* Load FT WAVE */
    snprintf(ext_regexp, 499, "^(%s|%s|%s|%s)$", GRAVI_ARRAY_GEOMETRY_EXT, GRAVI_OPTICAL_TRAIN_EXT, GRAVI_IMAGING_DATA_FT_EXT, GRAVI_IMAGING_DETECTOR_FT_EXT);
    gravi_data * wave_ft_data = gravi_data_load_rawframe_ext (frame, NULL, ext_regexp);

    /* Extract FT spectrum */
    gravi_data * ft_preproc_data = gravi_extract_spectrum (wave_ft_data, profile_map, dark_map,
            badpix_map, NULL, parlist, GRAVI_DET_FT);
    CPLCHECK_CLEAN ("Cannot extract FT spectrum");

    /* Move extensions necessary to compute the p2vmred */
    gravi_data_move_ext (ft_preproc_data, wave_ft_data, GRAVI_ARRAY_GEOMETRY_EXT);
    gravi_data_move_ext (ft_preproc_data, wave_ft_data, GRAVI_OPTICAL_TRAIN_EXT);
    FREE (gravi_data_delete, wave_ft_data);
    CPLCHECK_CLEAN ("Cannot move ext");

    /* Rescale to common wavelength */
    gravi_align_spectrum (ft_preproc_data, wave_map, p2vm_map, GRAVI_DET_FT, parlist);
    CPLCHECK_CLEAN ("Cannot re-interpolate FT spectrum");

    /* Compute P2VMRED */
    gravi_data * ft_p2vmred_data = gravi_compute_p2vmred(ft_preproc_data, 
            p2vm_map, "gravi_single", parlist, GRAVI_DET_FT);
    FREE (gravi_data_delete, ft_preproc_data);
    CPLCHECK_CLEAN ("Cannot apply p2vm");

    /* Merge extensions of SC and FT extracted spectra */
    gravi_data_move_ext(p2vmred_data, ft_p2vmred_data, GRAVI_OI_WAVELENGTH_EXT);
    gravi_data_move_ext(p2vmred_data, ft_p2vmred_data, GRAVI_OI_VIS_EXT);
    gravi_data_move_ext(p2vmred_data, ft_p2vmred_data, GRAVI_OI_FLUX_EXT);
    FREE (gravi_data_delete, ft_p2vmred_data);
    CPLCHECK_CLEAN ("Cannot merge SC and FT extensions");

    /* Perform the phase correction */
    if (!strcmp (gravi_param_get_string (parlist, "gravity.calib.phase-calibration"), "CLOSURE")) {
        gravi_p2vm_phase_correction (p2vm_map, p2vmred_data, 0);
    }
    else if (!strcmp (gravi_param_get_string (parlist, "gravity.calib.phase-calibration"), "DISP")) {
        gravi_p2vm_phase_correction (p2vm_map, p2vmred_data, 1);
    }
    else if (!strcmp (gravi_param_get_string (parlist, "gravity.calib.phase-calibration"), "FULL")) {
        gravi_p2vm_phase_correction (p2vm_map, p2vmred_data, 2);
    }
    else {
        cpl_msg_info (cpl_func, "P2VM phases are kept to zero (option phase-calibration=NONE)");
    }
    CPLCHECK_CLEAN ("Cannot recalibrate the P2VM phases");

    /* Add the OI_FLUX to further normalize the flux if needed */
    gravi_p2vm_transmission (p2vm_map, p2vmred_data);
    FREE (gravi_data_delete, p2vmred_data);

    CPLCHECK_CLEAN ("Cannot compute the transmission");

    /* 
     * (9) Create product frame, add DataFlow keywords, save the file, log the
     * saved file in the input frameset. Note that the output filename
     * is already created (from first P2VM file) 
     */

    gravi_data_save_new (p2vm_map, frameset, NULL, NULL, parlist,
                         used_frameset, frame_p2vm, "gravity_p2vm",
                         NULL, GRAVI_P2VM_MAP);

    CPLCHECK_CLEAN("Could not save the P2VM on the output file");

    /* Terminate the function */
    goto cleanup;
    
cleanup:
    /* Deallocation of all variables */
    cpl_msg_info (cpl_func,"Cleanup memory");

    FREE (cpl_frameset_delete, dark_frameset);
    FREE (cpl_frameset_delete, darkcalib_frameset);
    FREE (cpl_frameset_delete, badcalib_frameset);
    FREE (gravi_data_delete, spectrum_data);
    FREE (gravi_data_delete, p2vmred_data);
    FREE (gravi_data_delete, dark_map);
    FREELOOP (gravi_data_delete, raw_data, nb_frame_gain);
    FREE (gravi_data_delete, badpix_map);
    FREE (cpl_propertylist_delete, met_plist);
    FREE (cpl_frameset_delete, wavecalib_frameset);
    FREE (cpl_frameset_delete, flatcalib_frameset);
    FREE (cpl_frameset_delete, flat_frameset);
    FREE (cpl_frameset_delete, wave_frameset);
    FREE (cpl_frameset_delete, wavesc_frameset);
    FREE (cpl_frameset_delete, used_frameset);
    FREE (gravi_data_delete, data);
    FREE (gravi_data_delete, preproc_data);
    FREE (gravi_data_delete, dark_map);
    FREE (gravi_data_delete, wave_map);
    FREE (gravi_data_delete, profile_map);
    FREELOOP (cpl_free, valid_CP, 2);
    FREELOOP (cpl_free, valid_trans, 2);
    FREE (gravi_data_delete, p2vm_map);
    FREE (cpl_frameset_delete, p2vm_frameset);
    FREE (cpl_frameset_delete, current_frameset);
    FREE (cpl_frameset_delete, wave_param_frameset);
    FREE (gravi_data_delete, wave_param);

    //This is a workaround to aliviate PIPE-6316. For some reason the allocation
    //pattern of the recipe causes malloc to keep many pages in the allocation
    //arena which are not returned to the OS (typically calling brk()).
    //This causes issues if there is a fork() call, since then the whole 
    //address space of the parent is duplicated in the child. This is actually
    //the case with the system() call in esorex.
    //malloc_trim() is specific to GLIBC and will ask malloc() to return 
    //as many pages as possible. The code is not portable, hence the guards
#if defined(__linux__) && defined(__GLIBC__)
    malloc_trim(0);
#endif

    gravi_msg_function_exit(1);
    return (int)cpl_error_get_code();
}


