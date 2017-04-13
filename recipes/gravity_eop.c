/* 
 * This file is part of the GRAVITY Pipeline
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
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/*-----------------------------------------------------------------------------
                                Includes
 -----------------------------------------------------------------------------*/

#include <cpl.h>
#include <time.h>

#include "gravi_data.h"
#include "gravi_dfs.h"
#include "gravi_pfits.h"

#include "gravi_utils.h"

#include "gravi_eop.h"

/*-----------------------------------------------------------------------------
                            Private function prototypes
 -----------------------------------------------------------------------------*/

static int gravity_eop_create(cpl_plugin *);
static int gravity_eop_exec(cpl_plugin *);
static int gravity_eop_destroy(cpl_plugin *);
static int gravity_eop(cpl_frameset *, const cpl_parameterlist *);
cpl_error_code gravity_eop_compute_qc(cpl_table * eop_table, 
											   cpl_propertylist* header);

/*-----------------------------------------------------------------------------
                            Static variables
 -----------------------------------------------------------------------------*/

static const char gravity_eop_short[] = "Download the last values of the Earth Orientation Parameters and DUT from IERS.";
static const char gravity_eop_description[] =
"This recipe downloads the latest version of the Earth Orientation Parameter \n"
"and DUT from the IERS site. File is created in the current directory. A web connection is required.\n"
        GRAVI_RECIPE_FLOW"\n"
    "* Download the IERS data\n"
    "* Convert into CPL table\n"
    "* Write product\n"
    GRAVI_RECIPE_INPUT"\n"
    "None : No input\n"
    GRAVI_RECIPE_OUTPUT"\n"
    GRAVI_EOP_MAP"           : EOP calibration file (gravity_eop_calib.fits)\n"
    "";

static const char gravity_eop_name[] = "gravity_eop";

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
                    gravity_eop_name,
                    gravity_eop_short,
                    gravity_eop_description,
                    "Cesar Enrique Garcia Dabo",
                    PACKAGE_BUGREPORT,
                    "LL",
                    gravity_eop_create,
                    gravity_eop_exec,
                    gravity_eop_destroy)) {    
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
static int gravity_eop_create(cpl_plugin * plugin)
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

    /* --eop_host */
    p = cpl_parameter_new_value ("gravity.eop.eop_host",
                                CPL_TYPE_STRING, 
                                 "FTP Host to retrieve the EOP from",
                                 "gravity.eop", 
                                 "ftp.iers.org");
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "eop_host");
    cpl_parameter_disable (p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append (recipe->parameters, p);

    /* --eop_urlpath */
    p = cpl_parameter_new_value ("gravity.eop.eop_urlpath",
                                 CPL_TYPE_STRING, 
                                 "FTP URL path of the EOP file to retrieve",
                                 "gravity.eop", 
                                 "/products/eop/rapid/standard/finals2000A.data");
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "eop_urlpath");
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
static int gravity_eop_exec(cpl_plugin * plugin)
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
    recipe_status = gravity_eop(recipe->frames, recipe->parameters);
                                                                           
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
static int gravity_eop_destroy(cpl_plugin * plugin)
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
  @brief    Implement the recipe functionality
  @param    frameset   the frames list
  @param    parlist    the parameters list
  @return   0 if everything is ok
 */
/*----------------------------------------------------------------------------*/
static int gravity_eop(cpl_frameset            * frameset,
                                const cpl_parameterlist * parlist)
{
    const char          *   eop_host;
    const char          *   eop_urlpath;
    cpl_propertylist    *   applist;

	/* Message */
	gravity_print_banner (); 
	cpl_msg_set_time_on();
	cpl_msg_set_component_on();
	gravi_msg_function_start(1);

    /* Use the errorstate to detect an error in a function that does not
       return an error code. */
    cpl_errorstate          prestate = cpl_errorstate_get();

    /* Retrieving eop_host */
	eop_host = gravi_param_get_string (parlist, "gravity.eop.eop_host");

    /* Retrieving eop_urlpath */
	eop_urlpath = gravi_param_get_string (parlist, "gravity.eop.eop_urlpath");

    if (!cpl_errorstate_is_equal(prestate)) {
        return (int)cpl_error_set_message(cpl_func, cpl_error_get_code(),
                                          "Could not retrieve the input "
                                          "parameters");
    }
    

    /* Retrieve EOP file from the site */
    const char * eop_data;
    int          data_length;
    cpl_msg_info (cpl_func, "Retrieving EOP file ");
    eop_data = gravity_eop_download_finals2000A (eop_host, eop_urlpath, &data_length);

    if (eop_data == NULL || !cpl_errorstate_is_equal(prestate)) {
        return (int)cpl_error_set_message(cpl_func, cpl_error_get_code(),
                                        "Could not download data from server");
    }

	/* Convert to a CPL_TABLE */
    cpl_msg_info (cpl_func, "Convert EOP data to cpl_table");
    cpl_table * eop_table = gravity_eop_data_totable (eop_data, data_length);

    /* Check for a change in the CPL error state */
    cpl_ensure_code(cpl_errorstate_is_equal(prestate), cpl_error_get_code());
    
    applist = cpl_propertylist_new();

    /* Add the product category  */
    cpl_propertylist_append_string (applist, CPL_DFS_PRO_CATG, GRAVI_EOP_MAP);
    cpl_propertylist_append_string (applist, "ESO PRO TECH", "CATALOG");
    cpl_propertylist_append_string (applist, "ESO PRO TYPE", "IERS");

    /* Add a QC parameter  */
    gravity_eop_compute_qc (eop_table, applist);
    
    /* Saving the product */
    cpl_table_save (eop_table, applist, NULL, "gravity_eop.fits", CPL_IO_CREATE);

	cpl_msg_info (cpl_func,"Update the frameset");

    /* Updating the frameset */
    cpl_frame * product_frame = cpl_frame_new();
    cpl_frame_set_filename (product_frame, "gravity_eop.fits");
    cpl_frame_set_tag (product_frame, CPL_DFS_PRO_CATG);
    cpl_frame_set_type (product_frame, CPL_FRAME_TYPE_TABLE);
    cpl_frame_set_group (product_frame, CPL_FRAME_GROUP_PRODUCT);
    cpl_frame_set_level (product_frame, CPL_FRAME_LEVEL_FINAL);
    cpl_frameset_insert (frameset, product_frame);
    cpl_ensure_code (cpl_errorstate_is_equal(prestate), cpl_error_get_code());

    /* Cleanup */
    cpl_propertylist_delete (applist);
    cpl_table_delete (eop_table);

	gravi_msg_function_exit(1);
    return (int)cpl_error_get_code();
}

cpl_error_code gravity_eop_compute_qc (cpl_table * eop_table, 
                                       cpl_propertylist* header)
{
    double mjd_start;
    double mjd_lastfinal;
    double mjd_lastprediction;
    int null;

    mjd_start = cpl_table_get_double (eop_table, "MJD", 0, &null);
    for(int i = 0 ; i < cpl_table_get_nrow(eop_table);  i++)
    {
        const char * flag = cpl_table_get_string(eop_table, "FLAG", i);
        if(!strncmp(flag, "I", 1))
            mjd_lastfinal = cpl_table_get_double(eop_table, "MJD", i, &null);
        if(!strncmp(flag, "P", 1))
            mjd_lastprediction = cpl_table_get_double(eop_table, "MJD", i, &null);
    }
	
	cpl_msg_info (cpl_func, "QC EOP MJD START = %.3f", mjd_start);
	cpl_msg_info (cpl_func, "QC EOP MJD LAST FINAL = %.3f", mjd_lastfinal);
	cpl_msg_info (cpl_func, "QC EOP MJD LAST PREDICTION = %.3f", mjd_lastprediction);

    cpl_propertylist_append_double (header, "ESO QC EOP MJD START", mjd_start);
    cpl_propertylist_append_double (header, "ESO QC EOP MJD LAST FINAL", mjd_lastfinal);
    cpl_propertylist_append_double (header, "ESO QC EOP MJD LAST PREDICTION", mjd_lastprediction);
    cpl_propertylist_append_double (header, "MJD-OBS", mjd_lastfinal);
	
	return CPL_ERROR_NONE;
}
