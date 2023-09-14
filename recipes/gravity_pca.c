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
#include "gravi_utils.h"
#include "gravi_calib.h"

/*-----------------------------------------------------------------------------
                            Private function prototypes
 -----------------------------------------------------------------------------*/

static int gravity_pca_create(cpl_plugin *);
static int gravity_pca_exec(cpl_plugin *);
static int gravity_pca_destroy(cpl_plugin *);
static int gravity_pca(cpl_frameset *, const cpl_parameterlist *);

/*-----------------------------------------------------------------------------
                            Static variables
 -----------------------------------------------------------------------------*/
static char gravity_pca_short[] = "Calibrate the phase accuracy with the PCA method.";
static char gravity_pca_description[] = "TODO\n";

/*-----------------------------------------------------------------------------
                                Function code
 -----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
  @brief    Build the list of available plugins, for this module.
  @param    list the plugin list
  @return   0 if everything is ok, 1 otherwise
  @note     Only this function is exported

  Create the recipe instance and make it available to the application using the
  interface.
 */
/*----------------------------------------------------------------------------*/
int cpl_plugin_get_info(cpl_pluginlist *list)
{
    cpl_recipe *recipe = cpl_calloc(1, sizeof *recipe);
    cpl_plugin *plugin = &recipe->interface;

    if (cpl_plugin_init(plugin,
        CPL_PLUGIN_API,
        GRAVI_BINARY_VERSION,
        CPL_PLUGIN_TYPE_RECIPE,
        "gravity_phase_pca",
        gravity_pca_short,
        gravity_pca_description,
        "Calvin Sykes, Shangguan Jinyi, Sebastian Hoenig",
        PACKAGE_BUGREPORT,
        gravi_get_license(),
        gravity_pca_create,
        gravity_pca_exec,
        gravity_pca_destroy
    )) {
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
  @param    plugin the plugin
  @return   0 if everything is ok

  Defining the command-line/configuration parameters for the recipe.
 */
/*----------------------------------------------------------------------------*/
static int gravity_pca_create(cpl_plugin *plugin)
{
    cpl_recipe *recipe;

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
    gravi_parameter_add_static_name(recipe->parameters);

    /* PCA parameters */
    gravi_parameter_add_pca(recipe->parameters);

    return 0;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Execute the plugin instance given by the interface
  @param    plugin the plugin
  @return   0 if everything is ok
 */
/*----------------------------------------------------------------------------*/
static int gravity_pca_exec(cpl_plugin * plugin)
{

    cpl_recipe *recipe;
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
    recipe_status = gravity_pca(recipe->frames, recipe->parameters);

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
static int gravity_pca_destroy(cpl_plugin * plugin)
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
 * @brief Test whether tracking ratio for all baselines exceeds a limit.
 * 
 * @param hdr Header to extract ratios from.
 * @param min_ratio The threshold tracking ratio for acceptance.
 * 
 * @return Nonzero if accepted.
 **/
static int gravi_test_tracking_ratio(const cpl_propertylist *hdr, int min_ratio)
{
    const char *baselines[6] = {"12", "13", "14", "23", "24", "34"};
    const int nbase = 6;
    
    char tr_param_name[32];

    for (int i = 0; i < nbase; i++) {
        sprintf(tr_param_name, "ESO QC TRACKING_RATIO_FT%s", baselines[i]);
        const cpl_property *tr_prop = cpl_propertylist_get_property_const(hdr, tr_param_name);

        double tracking_ratio;
        if (cpl_property_get_type(tr_prop) == CPL_TYPE_DOUBLE)
            tracking_ratio = cpl_property_get_double(tr_prop);
        else if (cpl_property_get_type(tr_prop) == CPL_TYPE_INT)
            tracking_ratio = (double) cpl_property_get_int(tr_prop);
 
        if (tracking_ratio < min_ratio)
            return 0;
    }
    return 1;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    TODO.
  @param    frameset   the frames list
  @param    parlist    the parameters list
  @return   0 if everything is ok
 */
/*----------------------------------------------------------------------------*/
static int gravity_pca(cpl_frameset            * frameset,
		               const cpl_parameterlist * parlist)
{
    cpl_frameset *vis_sci_frameset = NULL;
    cpl_frameset *used_frameset = cpl_frameset_new();
    cpl_frame *frame = NULL;
    
    gravi_data *data_tmp, **data_accepted = NULL;
    cpl_propertylist *hdr = NULL, *hdr0 = NULL, *wave_plist = NULL, *wv_plisti = NULL;
    const char *telescope = NULL, *pola_mode = NULL, *spec_res = NULL;
    int nframes, nwave, nwavei, npol, naccept;
    const int nbase = 6;

    char product_filename[128];

	/* Message */
	gravity_print_banner ();
	cpl_msg_set_time_on();
	cpl_msg_set_component_on();
	gravi_msg_function_start(1);

    int min_tracking_ratio = cpl_parameter_get_int(
        cpl_parameterlist_find_const(parlist, "gravity.calib.pca-tracking-ratio"));

	/* Get the input frameset */
	cpl_ensure_code(gravi_dfs_set_groups(frameset) == CPL_ERROR_NONE,
                    cpl_error_get_code());

    vis_sci_frameset = gravi_frameset_extract_vis_science(frameset);
    if (cpl_frameset_is_empty(vis_sci_frameset)) {
        cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                              "No VIS_SCI on the frameset");
        goto cleanup;
    }

    /* Get header data from first frame */
    frame = cpl_frameset_get_position(vis_sci_frameset, 0);
    data_tmp = gravi_data_load_frame(frame, NULL);
    hdr0 = cpl_propertylist_duplicate(gravi_data_get_header(data_tmp));
    
    telescope = cpl_propertylist_get_string(hdr0, "TELESCOP");
    pola_mode = gravi_pfits_get_pola_mode(hdr0, GRAVI_SC);
    npol = gravi_pfits_get_pola_num(hdr0, GRAVI_SC);
    spec_res = gravi_pfits_get_spec_res(hdr0);
    
    /* Get length of wavelegth axis */
    wave_plist = gravi_data_get_oi_wave_plist(data_tmp, GRAVI_SC, 0, npol);
    nwave = cpl_propertylist_get_int(wave_plist, "NWAVE");

    FREE(gravi_data_delete, data_tmp);
    
    cpl_msg_info(cpl_func, "INS.NAME=%s  POLA.MODE=%s  SPEC.RES=%s  NPOL=%d  NWAVE=%d",
                 telescope, pola_mode, spec_res, npol, nwave);

    /* Check all frames for compatibility and reject if below tracking ratio */
    naccept = 0;
    nframes = cpl_frameset_get_size(vis_sci_frameset);
    data_accepted = cpl_malloc(nframes * sizeof(gravi_data *));
    //memset(data_accepted, 0, nframes * sizeof(gravi_data *));
    for (int n = 0; n < nframes; n++) {
        frame = cpl_frameset_get_position(vis_sci_frameset, n);
        data_tmp = gravi_data_load_frame(frame, NULL);

        CPLCHECK_CLEAN("Could not load input frame");

        hdr = gravi_data_get_header(data_tmp);        
        wv_plisti = gravi_data_get_oi_wave_plist(data_tmp, GRAVI_SC, 0, npol);
        nwavei = cpl_propertylist_get_int(wv_plisti, "NWAVE");
        if (nwave != nwavei)
            cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                  "Input files have inconsistent wavelength axes");
        if (strcmp(pola_mode, gravi_pfits_get_pola_mode(hdr, GRAVI_SC))) {
            cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                  "Input files have multiple INS.POLA.MODE");
            printf("%s\n", pola_mode);
            printf("%s\n", gravi_pfits_get_pola_mode(hdr, GRAVI_SC));
        }
        if (strcmp(spec_res, gravi_pfits_get_spec_res(hdr)))
            cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                  "Input files have multiple INS.SPEC.RES");
        
        /* Check if visibilities are all zeroes (TODO bad data?) */
        cpl_boolean skip = CPL_FALSE;
        cpl_table *vis_tmp = gravi_data_get_oi_vis(data_tmp, GRAVI_SC, n, npol);
        for (int i = 0; i < nbase; i++) {
            const cpl_array *vis_arr = cpl_table_get_array(vis_tmp, "VISPHI", i);
            if ((cpl_array_get_max(vis_arr) < 1e-2) || (cpl_array_get_min(vis_arr) > -1e-2)) {
                cpl_msg_warning(cpl_func, "Input file %s has invalid data and will be skipped [range=%f--%f].",
                    cpl_frame_get_filename(frame), cpl_array_get_max(vis_arr), cpl_array_get_min(vis_arr));
                skip = CPL_TRUE;
                break;
            }
        }

        /* Check if tracking ratio for file exceeds threshold and discard if not */
        if (!skip && gravi_test_tracking_ratio(hdr, min_tracking_ratio)) {
            data_accepted[naccept] = data_tmp;
            cpl_frameset_insert(used_frameset, frame);
            ++naccept;
        } else {
            FREE(gravi_data_delete, data_tmp);
            continue;
        }
    }
    cpl_msg_info(cpl_func, "Accepting %d of %d frames", naccept, nframes);
    data_accepted = cpl_realloc(data_accepted, naccept * sizeof(gravi_data *));

    if (naccept == 0)
        cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
            "None of the input files satisfy the minimum tracking ratio criterion");

    if (cpl_error_get_code())    
        goto cleanup;

    gravi_data *pca_data = gravi_compute_pca(data_accepted, naccept, parlist);
    CPLCHECK_CLEAN("Could not compute PCA decomposition");

    /* Add filenames for accepted files to output */
    cpl_errorstate e_state = cpl_errorstate_get();
    cpl_table *table = gravi_data_get_table(pca_data, GRAVI_PCA_RESID_EXT);
    if (table) {
        cpl_table_new_column(table, "FILES", CPL_TYPE_STRING);
        for (int i = 0; i < naccept; i++) {
            const cpl_frame *frame = cpl_frameset_get_position_const(used_frameset, i);
            cpl_table_set_string(table, "FILES", i, cpl_frame_get_filename(frame));
        }
    } else {
        /* Catch error from accessing null table */
        cpl_errorstate_set(e_state);
    }

    sprintf(product_filename, "GRAVI_PhaseCalib_%s_%s_%s", spec_res, pola_mode, telescope);
    gravi_data_save_new(pca_data, frameset, product_filename, NULL, parlist,
                        used_frameset, NULL, "gravity_phase_pca",
                        NULL, GRAVI_PHASE_PCA);

cleanup:
    FREE(cpl_propertylist_delete, hdr0);
    FREE(cpl_frameset_delete, vis_sci_frameset);
    FREELOOP(gravi_data_delete, data_accepted, naccept);

    return (int)cpl_error_get_code();
}
