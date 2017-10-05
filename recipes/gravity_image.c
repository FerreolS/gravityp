/* $Id: gravity_image.c,v 1.29 2009/02/10 09:16:12 llundin Exp $
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
 * $Name: HEAD $
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/*-----------------------------------------------------------------------------
                                Includes
 -----------------------------------------------------------------------------*/

#include <cpl.h>
#include <ctype.h>

#include "gravi_utils.h"
#include "gravi_pfits.h"
#include "gravi_dfs.h"
#include "gravi_image.h"

/*-----------------------------------------------------------------------------
                            Private function prototypes
 -----------------------------------------------------------------------------*/

static int gravity_image_create(cpl_plugin *);
static int gravity_image_exec(cpl_plugin *);
static int gravity_image_destroy(cpl_plugin *);
static int gravity_image(cpl_frameset *, const cpl_parameterlist *);

/*-----------------------------------------------------------------------------
                            Static variables
 -----------------------------------------------------------------------------*/
static char gravity_image_short[] = GRAVI_UNOFFERED"Reconstruct an image from visibilities.";
static char gravity_image_description[] = GRAVI_UNOFFERED
"This recipe launch a Yorick batch executing mira-script.i to process an\n"
"input OIFITS file to produce an image in fits file.\n"
"The input frame is an OIFITS file :\n"
"GRAVI-GRAVI_MIRA-input-file.oifits " GRAVI_MIRA_INPUT_PROCATG "\n"
"and the recipe generates a fits image\n"
"GRAVI-GRAVI_MIRA-image-file.fits " GRAVI_MIRA_OUTPUT_PROCATG "\n";

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
                    "gravity_image",
                    gravity_image_short,
                    gravity_image_description,
                    "Vincent Lapeyrere",
                    PACKAGE_BUGREPORT,
                    gravi_get_license(),
                    gravity_image_create,
                    gravity_image_exec,
                    gravity_image_destroy)) {
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
static int gravity_image_create(cpl_plugin * plugin)
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

    gravi_parameter_add_image (recipe->parameters);

    return 0;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Execute the plugin instance given by the interface
  @param    plugin  the plugin
  @return   0 if everything is ok
 */
/*----------------------------------------------------------------------------*/
static int gravity_image_exec(cpl_plugin * plugin)
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
    recipe_status = gravity_image(recipe->frames, recipe->parameters);
                                                                           
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
static int gravity_image_destroy(cpl_plugin * plugin)
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
  @brief    Interpret the command line options and execute the data processing
  @param    frameset   the frames list
  @param    parlist    the parameters list
  @return   0 if everything is ok
 */
/*----------------------------------------------------------------------------*/
//#define YORICK_BIN 1
static int gravity_image(cpl_frameset            * frameset,
                    const cpl_parameterlist * parlist)
{
#ifdef YORICK_BIN
    const cpl_frame     *   rawframe;
    double                  qc_param = 0.0;
    cpl_propertylist    *   applist;
    cpl_image           *   image;
    cpl_frameset * usedframes;

    /* Use the errorstate to detect an error in a function that does not
       return an error code. */
    cpl_errorstate          prestate = cpl_errorstate_get();

    /* RETRIEVE INPUT PARAMETERS */
    /* No input param to retrieve */
    if (!cpl_errorstate_is_equal(prestate)) {
        return (int)cpl_error_set_message(cpl_func, cpl_error_get_code(),
                                          "Could not retrieve the input "
                                          "parameters");
    }
    
    /* Identify the RAW and CALIB frames in the input frameset */
    cpl_ensure_code(gravi_dfs_set_groups(frameset) == CPL_ERROR_NONE,
                    cpl_error_get_code());

    /* ACCESS INPUT DATA */
    /*  Get the image reconstruction INPUT file  */
    rawframe = cpl_frameset_find_const(frameset, GRAVI_MIRA_INPUT_PROCATG);
    if (rawframe == NULL) {
        /* cpl_frameset_find_const() does not set an error code, when a frame
           is not found, so we will set one here. */
        return (int)cpl_error_set_message(cpl_func, CPL_ERROR_DATA_NOT_FOUND,
                                          "SOF does not have any file tagged "
                                          "with %s", GRAVI_MIRA_INPUT_PROCATG);
    }

    /* Check for a change in the CPL error state */
    /* - if it did change then propagate the error and return */
    cpl_ensure_code(cpl_errorstate_is_equal(prestate), cpl_error_get_code());

    cpl_msg_set_component_on();

    /* get file information */
    const char * filename = cpl_frame_get_filename(rawframe);
    gravi_data * input_data = gravi_data_load_ext(filename, "OI_TARGET");
    cpl_table * oi_target_table = gravi_data_get_table(input_data, "OI_TARGET");
    int n_target = cpl_table_get_nrow(oi_target_table);
    
    /* loop on targets */
    for (int i_target = 0; i_target < n_target; i_target++) {

        /* NOW PERFORMING THE DATA REDUCTION */
        /* Execute the gravi_image function */
        char *target_name = cpl_table_get_string(oi_target_table, "TARGET", i_target);
        //char *target_name = cpl_sprintf("%s", "IRS16C");
        image=gravi_image(rawframe, parlist, target_name);
        if (image == NULL) {
            /* cpl_frameset_find_const() does not set an error code, when a frame
               is not found, so we will set one here. */
            return (int)cpl_error_set_message(cpl_func, CPL_ERROR_NULL_INPUT,
                                              "The gravi_image function return NULL pointer");
        }

        applist = cpl_propertylist_new();

        /* Add the product category  */
        char proCatg[100];
        sprintf(proCatg, "%s", GRAVI_MIRA_OUTPUT_PROCATG);
        cpl_propertylist_append_string(applist, CPL_DFS_PRO_CATG, proCatg);
        /* Select the name extension depending on this catg */
        char catg_ext[100];
        int i, j=0;
        for(i = 0; proCatg[i]; i++) if (proCatg[i]!='_') { catg_ext[j] = tolower(proCatg[i]); j++; }
        catg_ext[j] = '\0';

        /* Add a QC parameter  */
        cpl_propertylist_append_double(applist, "ESO QC QCPARAM", qc_param);

        /* SAVE A DFS-COMPLIANT PRODUCT TO DISK*/
         gravi_dfs_set_groups(frameset);

        /* create product */
        usedframes = cpl_frameset_new();
        cpl_frameset_insert(usedframes, cpl_frame_duplicate(rawframe));

        cpl_msg_info(cpl_func, "Writing image_out.fits");

        /*	if (cpl_dfs_save_propertylist (frameset, NULL, parlist,	usedframes,
                rawframe, "gravity_image", applist, NULL,
                PACKAGE "/" PACKAGE_VERSION , "image_out.fits" ) != CPL_ERROR_NONE){
            cpl_error_set_message(cpl_func, cpl_error_get_code(),
                              "Cannot save the first extension primary header");
            cpl_image_delete(image);
            cpl_propertylist_delete(applist);
            cpl_frameset_delete(usedframes);
            return cpl_error_get_code();
        }
        */
        /* Save the image extensions */
        /*	if (cpl_image_save(image, "image_out.fits", CPL_BPP_IEEE_DOUBLE,
                                applist, CPL_IO_EXTEND) != CPL_ERROR_NONE){
            cpl_error_set_message(cpl_func, cpl_error_get_code(),
                              "Cannot save the image extension");
            cpl_image_delete(image);
            cpl_propertylist_delete(applist);
            cpl_frameset_delete(usedframes);
            return cpl_error_get_code();
        }
        */

        /* Use either the input filename or a name based on recipe */
        char * product_name = NULL;

        if ( cpl_parameterlist_find_const (parlist,"gravity.dfs.static-name") &&
             gravi_param_get_bool (parlist,"gravity.dfs.static-name")) {

          /* Use the recipe name and add the selected catg_ext */
          product_name = cpl_sprintf ("%s_%s_%s.fits", "gravity_image", catg_ext, target_name);
        }
        else {

          /* Remove the extension (last '.') and add the selected catg_ext */
          char * filenoext = cpl_strdup (FILESHORT(filename));
          char * lastdot = strrchr (filenoext, '.');
          if (lastdot != NULL) *lastdot = '\0';
          product_name = cpl_sprintf ("%s_%s_%s.fits", filenoext, catg_ext, target_name);
          FREE (cpl_free, filenoext);
        }

        cpl_dfs_save_image(frameset, NULL, parlist, usedframes,
                rawframe, image, CPL_BPP_IEEE_DOUBLE,
                "gravity_image", applist, NULL,
                PACKAGE "/" PACKAGE_VERSION , product_name);
        cpl_msg_info(cpl_func, "Reconstructed image saved in image_out.fits");

        /* free memory */
        cpl_image_delete(image);
        cpl_propertylist_delete(applist);
        cpl_frameset_delete(usedframes);
        FREE (cpl_free, product_name);
    } // end loop on target
    gravi_data_delete(input_data);


    return (int)cpl_error_get_code();
#else
    return (int)cpl_error_set_message(cpl_func, CPL_ERROR_NULL_INPUT,
                                      "This recipe is only available if "
                                      "the pipeline was compiled with yorick "
                                      "support. Check configure --help");
#endif
}
