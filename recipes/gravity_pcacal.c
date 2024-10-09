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

static int gravity_pcacal_create(cpl_plugin *);
static int gravity_pcacal_exec(cpl_plugin *);
static int gravity_pcacal_destroy(cpl_plugin *);
static int gravity_pcacal(cpl_frameset *, const cpl_parameterlist *);

/*-----------------------------------------------------------------------------
                            Static variables
 -----------------------------------------------------------------------------*/
static char gravity_pcacal_short[] = "Generate static calibration files for flattening phase visibility data using the PCA method.";
static char gravity_pcacal_description[] =
"This recipe produces a PCA calibration file from a set of calibration frames to be used for flattening phase visibility data.\n"
    GRAVI_RECIPE_FLOW"\n"
    "* Select good input frames using tracking ratio criterion.\n"
    "* Compute PCA decomposition for each baseline and polarisation channel\n"
    "* Fit component model and write calibration product\n"
    GRAVI_RECIPE_INPUT"\n"
    GRAVI_VIS_SINGLE_CALIB"  ≥20  : input frames\n"
    GRAVI_RECIPE_OUTPUT"\n"
    GRAVI_PHASE_PCA"      : PCA calibration\n"
    "";

/* Detector "epochs" corresponding to thermal cycles: calibration valid within interval from one date to next. */
#define N_EPOCH 3
static double TIME_MJD_EPOCH_START[N_EPOCH] = {
    57754.000000, // 2017-01-01
    58758.000000, // 2019-10-02
    59178.000000, // 2020-11-25
};

/* Minimum number of valid calibration frames to accept */
static cpl_size MIN_CALIB_FRAMES = 20;

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
        "gravity_pcacal",
        gravity_pcacal_short,
        gravity_pcacal_description,
        "Calvin Sykes, Shangguan Jinyi, Sebastian Hoenig",
        PACKAGE_BUGREPORT,
        gravi_get_license(),
        gravity_pcacal_create,
        gravity_pcacal_exec,
        gravity_pcacal_destroy
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
static int gravity_pcacal_create(cpl_plugin *plugin)
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
    gravi_parameter_add_pcacalib(recipe->parameters);

    return 0;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Execute the plugin instance given by the interface
  @param    plugin the plugin
  @return   0 if everything is ok
 */
/*----------------------------------------------------------------------------*/
static int gravity_pcacal_exec(cpl_plugin * plugin)
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
    recipe_status = gravity_pcacal(recipe->frames, recipe->parameters);

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
static int gravity_pcacal_destroy(cpl_plugin * plugin)
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
 * @brief Test whether tracking ratio for all baselines exceeds a limit.
 * 
 * @param hdr Header to extract ratios from.
 * @param min_ratio The threshold tracking ratio for acceptance.
 * 
 * @return CPL_TRUE if accepted.
 **/
/*----------------------------------------------------------------------------*/
static cpl_boolean gravi_test_tracking_ratio(const cpl_propertylist *hdr, int min_ratio)
{
    const int nbase = 6;    
    char tr_param_name[100];

    for (int i = 0; i < nbase; i++) {
        sprintf(tr_param_name, "ESO QC TRACKING_RATIO_FT%s", GRAVI_BASE_NAME[i]);
        const cpl_property *tr_prop = cpl_propertylist_get_property_const(hdr, tr_param_name);

        int tracking_ratio = 0;
        if (cpl_property_get_type(tr_prop) == CPL_TYPE_INT)
            tracking_ratio = cpl_property_get_int(tr_prop);
        else if (cpl_property_get_type(tr_prop) == CPL_TYPE_DOUBLE)
            tracking_ratio = (int) cpl_property_get_double(tr_prop);
        else
            cpl_msg_error(cpl_func, "Could not get tracking ratio");
 
        if (tracking_ratio < min_ratio)
            return CPL_FALSE;
    }
    return CPL_TRUE;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief    Compute the PCA model from the provided calibration data.
 * @param    frameset   the frames list
 * @param    parlist    the parameters list
 * @return   0 if everything is ok
 */
/*----------------------------------------------------------------------------*/
static int gravity_pcacal(cpl_frameset            * frameset,
		                  const cpl_parameterlist * parlist)
{
    cpl_frameset *vis_cal_frameset = NULL;
    cpl_frameset *used_frameset = cpl_frameset_new();
    cpl_frame *frame = NULL;

    gravi_data *data_tmp = NULL, **data_accepted = NULL, *pca_data = NULL;
    cpl_propertylist *hdr = NULL, *header_first = NULL, *wave_plist = NULL, *wv_plisti = NULL;
    const char *telescope = NULL, *pola_mode = NULL, *spec_res = NULL;
    int nframes, nwave, nwavei, npol, naccept;
    const int nbase = 6;
    char product_filename[100];

	/* Message */
	gravity_print_banner ();
	gravi_msg_function_start(1);

    int min_tracking_ratio = cpl_parameter_get_int(
        cpl_parameterlist_find_const(parlist, "gravity.calib.pca-tracking-ratio"));

	/* Get the input frameset */
	cpl_ensure_code(gravi_dfs_set_groups(frameset) == CPL_ERROR_NONE,
                    cpl_error_get_code());

    vis_cal_frameset = gravi_frameset_extract_vis_calib(frameset);
    if (cpl_frameset_is_empty(vis_cal_frameset)) {
        cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                              "No VIS_CAL on the frameset");
        goto cleanup;
    }

    /* Get header data from first frame */
    frame = cpl_frameset_get_position(vis_cal_frameset, 0);
    data_tmp = gravi_data_load_frame(frame, NULL);
    header_first = cpl_propertylist_duplicate(gravi_data_get_header(data_tmp));
    
    telescope = cpl_propertylist_get_string(header_first, "TELESCOP");
    pola_mode = gravi_pfits_get_pola_mode(header_first, GRAVI_SC);
    npol = gravi_pfits_get_pola_num(header_first, GRAVI_SC);
    spec_res = gravi_pfits_get_spec_res(header_first);

    /* Check on time */
    double time_mjd_obs = cpl_propertylist_get_double(header_first, "MJD-OBS");
    int epoch = -1;
    if (time_mjd_obs < TIME_MJD_EPOCH_START[0]) {
        cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
            "First frame is too old\n");
    }
    for (int i = 0; i < N_EPOCH; i++) {
        /* Select latest epoch date that precedes observation date */
        if (TIME_MJD_EPOCH_START[i] <= time_mjd_obs)
            epoch = i;
    }
    
    /* Get length of wavelength axis */
    wave_plist = gravi_data_get_oi_wave_plist(data_tmp, GRAVI_SC, 0, npol);
    nwave = cpl_propertylist_get_int(wave_plist, "NWAVE");

    FREE(gravi_data_delete, data_tmp);
    
    cpl_msg_info(cpl_func, "INS.NAME=%s  POLA.MODE=%s  SPEC.RES=%s  NPOL=%d  NWAVE=%d",
                 telescope, pola_mode, spec_res, npol, nwave);

    /* Check all frames for compatibility and reject if below tracking ratio */
    naccept = 0;
    nframes = cpl_frameset_get_size(vis_cal_frameset);
    data_accepted = cpl_malloc(nframes * sizeof(gravi_data *));
    for (int n = 0; n < nframes; n++) {
        frame = cpl_frameset_get_position(vis_cal_frameset, n);
        data_tmp = gravi_data_load_frame(frame, NULL);
        CPLCHECK_CLEAN("Could not load input frame");

        /* Check for uniform wavelength axis */
        hdr = gravi_data_get_header(data_tmp);        
        wv_plisti = gravi_data_get_oi_wave_plist(data_tmp, GRAVI_SC, 0, npol);
        nwavei = cpl_propertylist_get_int(wv_plisti, "NWAVE");
        if (nwave != nwavei)
            cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                  "Input files have inconsistent wavelength axes");
        
        /* Check for matching telescope */
        if (strcmp(telescope, cpl_propertylist_get_string(hdr, "TELESCOP"))) {
            cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                  "Input files have multiple TELESCOP");
        }

        /* Check for matching polarisation mode */
        if (strcmp(pola_mode, gravi_pfits_get_pola_mode(hdr, GRAVI_SC))) {
            cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                  "Input files have multiple INS.POLA.MODE");
        }

        /* Check for matching resolution */
        if (strcmp(spec_res, gravi_pfits_get_spec_res(hdr))) {
            cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                  "Input files have multiple INS.SPEC.RES");
        }

        /* Check for matching epoch */
        time_mjd_obs = cpl_propertylist_get_double(hdr, "MJD-OBS");
        if (time_mjd_obs < TIME_MJD_EPOCH_START[epoch] ||
            (epoch < N_EPOCH-1 && time_mjd_obs > TIME_MJD_EPOCH_START[epoch+1])) {
            cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                  "Input files span multiple epochs");
        }
        
        /* Check if visibilities are all zeroes (bad data?) */
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
    
    if (naccept < MIN_CALIB_FRAMES)
        cpl_msg_warning(cpl_func, "Fewer than %lld valid frames provided, calibration may be inaccurate", MIN_CALIB_FRAMES);

    if (cpl_error_get_code())    
        goto cleanup;

    pca_data = gravi_compute_pca(data_accepted, naccept, parlist);
    CPLCHECK_CLEAN("Could not compute PCA decomposition");

    cpl_propertylist *product_header = gravi_data_get_plist(pca_data, GRAVI_PCA_EXT);
    cpl_propertylist_append_double(product_header, "PCA EPOCH BEGIN", TIME_MJD_EPOCH_START[epoch]);
    cpl_propertylist_append_double(product_header, "PCA EPOCH END", epoch < N_EPOCH-1 ? TIME_MJD_EPOCH_START[epoch+1] : 99999);

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

    char* mjd_str = gravi_convert_to_timestamp(TIME_MJD_EPOCH_START[epoch]);
    sprintf(product_filename, "GRAVI.%s.%s_%s_%s.", mjd_str, spec_res, pola_mode, telescope);
    gravi_data_save_new(pca_data, frameset, product_filename, NULL, parlist,
                        used_frameset, NULL, "gravity_phase_pca",
                        NULL, GRAVI_PHASE_PCA);
    // frame = cpl_frameset_get_position(vis_cal_frameset, 0);
    // gravi_data_save_new(pca_data, frameset, NULL, NULL, parlist,
    //                     used_frameset, frame, "gravity_phase_pca",
    //                     product_header, GRAVI_PHASE_PCA);

cleanup:
    FREE(cpl_propertylist_delete, header_first);
    FREE(cpl_frameset_delete, vis_cal_frameset);
    FREELOOP(gravi_data_delete, data_accepted, naccept);
    FREE(gravi_data_delete, pca_data);

    gravi_msg_function_exit(1);
    return (int)cpl_error_get_code();
}
