/*
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
#include "gravi_vis.h"
#include "gravi_utils.h"

#include "gravi_astrometry.h"

/*-----------------------------------------------------------------------------
                            Private function prototypes
 -----------------------------------------------------------------------------*/

static int gravity_astrometry_create(cpl_plugin *);
static int gravity_astrometry_exec(cpl_plugin *);
static int gravity_astrometry_destroy(cpl_plugin *);
static int gravity_astrometry(cpl_frameset *, cpl_parameterlist *);

/*-----------------------------------------------------------------------------
                            Static variables
 -----------------------------------------------------------------------------*/
static char gravity_astro_short[] = "Compute astrometric phase reference";
static char gravity_astro_description[] =
    "This recipe computes phase and amplitude referencing for astrometric observations.\n"
    "It supports on- and off-axis observing strategies, as well as the use of swaps.\n"
    GRAVI_RECIPE_FLOW"\n"
    "* If swaps are present: obtain astrometric solution and compute swap phase reference\n"
    "* Compute phase reference for the target.\n"
    "* Write output product with correctly referenced phase.\n"
    GRAVI_RECIPE_INPUT"\n"
    GRAVI_ASTRO_CAL_PHASEREF" : star frames to use for phase referencing\n"
    GRAVI_ASTRO_TARGET":\tplanet frames to be referenced\n"
    GRAVI_ASTRO_SWAP":\talternating star/planet frames for swap observing mode\n"
    GRAVI_RECIPE_OUTPUT"\n"
    GRAVI_ASTRO_PHASE_CALIBRATED" : output astroreduced file with correctly referenced phase\n"
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
                    "gravity_astrometry",
                    gravity_astro_short,
                    gravity_astro_description,
                    "Calvin Sykes, Mathias Nowak, Sebastian Hoenig",
                    PACKAGE_BUGREPORT,
                    gravi_get_license(),
                    gravity_astrometry_create,
                    gravity_astrometry_exec,
                    gravity_astrometry_destroy)) {
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
static int gravity_astrometry_create(cpl_plugin * plugin)
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

    /* Use static names (output_procatg.fits) */
    gravi_parameter_add_static_name (recipe->parameters);

    /* Astrometry parameters */
    gravi_parameter_add_astrometry (recipe->parameters);
    
	return 0;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Execute the plugin instance given by the interface
  @param    plugin  the plugin
  @return   0 if everything is ok
 */
/*----------------------------------------------------------------------------*/
static int gravity_astrometry_exec(cpl_plugin * plugin)
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
    recipe_status = gravity_astrometry(recipe->frames, recipe->parameters);

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
static int gravity_astrometry_destroy(cpl_plugin * plugin)
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

/**
 * @brief Load input astrometric quantities from ASTROREDUCED file(s).
 * 
 * @param frameset Frameset to load from.
 * @param used_frameset Frameset to store used frames in.
 * @param out_data Output array to receive loaded @c gravi_data objects. Must be allocated, or NULL to skip.
 * 
 * @return Array of @c astro_data objects. Will be allocated by this function.
*/
static astro_data** load_data(cpl_frameset *frameset, cpl_frameset *used_frameset, gravi_data **out_data)
{
    cpl_frame *frame = NULL;
    gravi_data *tmp_data = NULL;
    astro_data **data = NULL;

    cpl_size nframes = cpl_frameset_get_size(frameset);

    data = cpl_calloc(nframes, sizeof(astro_data *));
    for (int i = 0; i < nframes; i++) {
        frame = cpl_frameset_get_position(frameset, i);
        tmp_data = gravi_data_load_frame(frame, used_frameset);

        data[i] = gravi_astrometry_load(tmp_data);
        if (out_data)
            out_data[i] = tmp_data;
        else
            FREE(gravi_data_delete, tmp_data);
        CPLCHECK_CLEAN("Could not load data");
    }
    return data;

cleanup:
    FREELOOP(gravi_astrometry_delete, data, nframes);
    return NULL;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    TODO
  @param    frameset   the frames list
  @param    parlist    the parameters list
  @return   0 if everything is ok
 */
/*----------------------------------------------------------------------------*/
static int gravity_astrometry(cpl_frameset * frameset,
                       cpl_parameterlist * parlist)
{    
    cpl_frameset *target_frameset=NULL, *swap_frameset=NULL, *phaseref_frameset=NULL,
                 *used_frameset=NULL;
	
	cpl_frame *frame=NULL;
    cpl_size n_target = 0, n_swap = 0, n_phaseref = 0;
    gravi_data **tgt_data = NULL, **swap_data = NULL;
    astro_data **tgt_astro = NULL, **swap_astro=NULL, **phaseref_astro=NULL;
	
	/* Message */
	gravity_print_banner();
	cpl_msg_set_time_on();
	cpl_msg_set_component_on();
	gravi_msg_function_start(1);

    cpl_boolean debug = cpl_parameter_get_bool(
        cpl_parameterlist_find(parlist, "gravity.astrometry.wait-for-debugger"));

    if (debug) {
        fprintf(stderr, "PID is: %d\n", getpid());
        fprintf(stderr, "Waiting for debugger to attach...\n");
        volatile int attach = 1;
        while (attach) {
            if (attach == 0) break;
        }
    }

    // const char *target = cpl_parameter_get_string(
    //     cpl_parameterlist_find_const(parlist, "gravity.astrometry.target-name")
    // );
    // const char *swap_target = cpl_parameter_get_string(
    //     cpl_parameterlist_find_const(parlist, "gravity.astrometry.swap-target-name")
    // );

    double ft_mean_flux_threshold = cpl_parameter_get_double(
        cpl_parameterlist_find(parlist, "gravity.astrometry.ft-mean-flux"));

    cpl_ensure_code(gravi_dfs_set_groups(frameset) == CPL_ERROR_NONE,
					cpl_error_get_code()) ;

    /* Dispatch the frameset */
    target_frameset = gravi_frameset_extract_astro_target(frameset);
    swap_frameset = gravi_frameset_extract_astro_swap(frameset);
    phaseref_frameset = gravi_frameset_extract_astro_phaseref(frameset);

    // if (cpl_frameset_is_empty(target_frameset)) {
    //     cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
    //                         "No ASTRO_TARGET on the frameset");
    //     goto cleanup;
    // }

    // if ((swap = cpl_frameset_is_empty(swap_frameset))) {
    //     cpl_msg_debug(cpl_func, "No ASTRO_SWAP on the frameset, assuming onaxis");
    // }

    // if (cpl_frameset_is_empty(phaseref_frameset)) {
    //     cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
    //                           "No ASTRO_CAL_PHASEREF on the frameset");
    //     goto cleanup;
    // }

	/* Insert calibration frame into the used frameset */
	used_frameset = cpl_frameset_new();

    /* Target data */
    n_target = cpl_frameset_get_size(target_frameset);
    tgt_data = cpl_malloc(n_target * sizeof(gravi_data*));
    tgt_astro = load_data(target_frameset, used_frameset, tgt_data);
    cpl_msg_debug(cpl_func, "There are %lld ASTRO_TARGET frames", n_target);
    CPLCHECK_CLEAN("Could not load target data");

    /* Swap data */
    n_swap = cpl_frameset_get_size(swap_frameset);
    swap_data = cpl_malloc(n_swap * sizeof(gravi_data*));
    swap_astro = load_data(swap_frameset, used_frameset, swap_data);
    cpl_msg_debug(cpl_func, "There are %lld ASTRO_SWAP frames", n_swap);
    CPLCHECK_CLEAN("Could not load target data");

    /* Data to be used for phase referencing */
    n_phaseref = cpl_frameset_get_size(phaseref_frameset);
    phaseref_astro = load_data(phaseref_frameset, used_frameset, NULL);
    cpl_msg_debug(cpl_func, "There are %lld ASTRO_CAL_PHASEREF frames", n_phaseref);
    CPLCHECK_CLEAN("Could not load phaseref data");

    /* If there are SWAP files, reduce them first */
    /* swap files are modified in-place to store the astrometry */
    if (n_swap > 0)
        gravi_astrometry_reduce_swaps(swap_astro, n_swap, parlist);
    CPLCHECK_CLEAN("Could not reduce swaps");

    /* If there are targets, calculate phase reference */
    /* TODO: if there are no targets, the recipe runs but does nothing */
    /* TODO: it might be desirable to be able to (f.e.) reduce swaps in isolation */
    if (n_target > 0) {
        double ft_mean_flux_tgt = 0.0;
        for (int i = 0; i < n_target; i++)
            ft_mean_flux_tgt += gravi_astrometry_get_mean_ftflux(tgt_astro[i]);
        ft_mean_flux_tgt /= n_target;
        cpl_msg_debug(cpl_func, "ftOnPlanetMeanFlux=%f", ft_mean_flux_tgt);
        
        double ft_mean_flux_ref = 0.0;
        for (int i = 0; i < n_phaseref; i++)
            ft_mean_flux_ref += gravi_astrometry_get_mean_ftflux(phaseref_astro[i]);
        ft_mean_flux_ref /= n_phaseref;
        cpl_msg_debug(cpl_func, "ftOnStarMeanFlux=%f", ft_mean_flux_ref);

        /* Filter based on FT flux threshold, before normalising by coherent FT flux */
        for (int i = 0; i < n_target; i++) {
            gravi_astrometry_filter_ftflux(tgt_astro[i], ft_mean_flux_threshold * ft_mean_flux_tgt);
            gravi_astrometry_normalise_to_ft(tgt_astro[i]);
        }

        for (int i = 0; i < n_phaseref; i++) {
            gravi_astrometry_filter_ftflux(phaseref_astro[i], ft_mean_flux_threshold * ft_mean_flux_ref);
            gravi_astrometry_normalise_to_ft(phaseref_astro[i]);
        }

        cpl_msg_info(cpl_func, "Creating visibility reference from %lld observations", n_phaseref);
        for (int i = 0; i < n_target; i++)
            gravi_astrometry_create_phase_reference(tgt_astro[i], phaseref_astro, n_phaseref, swap_astro, n_swap, parlist);
        CPLCHECK_CLEAN("Could not calculate phase reference");

        for (int i = 0; i < n_target; i++) {
            cpl_table *phaseref_table = gravi_astrometry_get_phase_reference(tgt_astro[i]);

            frame = cpl_frameset_get_position(target_frameset, i);
            gravi_data_add_table(tgt_data[i], NULL, "ASTRO_VISREF", phaseref_table);
            gravi_data_save_new(tgt_data[i], frameset, NULL, NULL, parlist,
                            used_frameset, frame, "gravity_astrometry",
                            NULL, GRAVI_ASTRO_PHASE_CALIBRATED);
            CPLCHECK_CLEAN("Could not save ASTRO_VISREF product");
        }
    }

	/* Terminate the function */
	goto cleanup;

cleanup:
	/* Deallocation of all variables */
	cpl_msg_info(cpl_func, "Memory cleanup");

    FREE(cpl_frameset_delete, target_frameset);
    FREE(cpl_frameset_delete, phaseref_frameset);
    FREE(cpl_frameset_delete, used_frameset);

    FREELOOP(gravi_astrometry_delete, tgt_astro, n_target);
    FREELOOP(gravi_astrometry_delete, swap_astro, n_swap);
    FREELOOP(gravi_astrometry_delete, phaseref_astro, n_phaseref);

	gravi_msg_function_exit(1);
    return (int)cpl_error_get_code();
}
