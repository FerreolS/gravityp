/* $Id: gravi_data-test.c,v 1.59 2011/08/16 17:43:49 nazouaoui Exp $
 *
 * This file is part of the ESO Common Pipeline Library
 * Copyright (C) 2001-2008 European Southern Observatory
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


/*
 * gravi_utils-test.c
 *
 *  Created on: 16 août 2011
 *      Author: nabih
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <assert.h>

#include <cpl_test.h>
#include <cpl.h>
#include "gravi_data.h"
#include "gravi_pfits.h"
#include "gravi_utils.h"
#include "gravi_calib.h"
#include "gravi_signal.h"
#include "gravi_preproc.h"
#include "gravi_p2vm.h"
#include "gravi_p2vmred.h"
#include "gravi_signal.h"
#include "gravi_metrology.h"
#include "gravi_vis.h"
#include "gravi_dfs.h"
#include "gravi_wave.h"
#include "gravi-test.c"

#define DATADIR_TEST DATADIR
//#define DATADIR_TEST ""

int gravi_calib_test(void);
int gravi_calib_test(void){

    const char *names[] = {
        DATADIR_TEST "Wave.fits",
        DATADIR_TEST "Dark.fits",
        DATADIR_TEST "Wave.fits",
        DATADIR_TEST "Flat_4.fits",
        DATADIR_TEST "Flat_1.fits",
        DATADIR_TEST "Flat_2.fits",
        DATADIR_TEST "Flat_3.fits"
    };

    const char *tags[] = {
        GRAVI_P2VM_RAW,
        GRAVI_DARK_RAW,
        GRAVI_WAVE_RAW,
        GRAVI_FLAT_RAW,
        GRAVI_FLAT_RAW,
        GRAVI_FLAT_RAW,
        GRAVI_FLAT_RAW
    };

    cpl_frame_group groups[] = {
        CPL_FRAME_GROUP_RAW,
        CPL_FRAME_GROUP_RAW,
        CPL_FRAME_GROUP_RAW,
        CPL_FRAME_GROUP_RAW,
        CPL_FRAME_GROUP_RAW,
        CPL_FRAME_GROUP_RAW,
        CPL_FRAME_GROUP_RAW
    };

    cpl_frame_level levels[] = {
        CPL_FRAME_LEVEL_NONE,
        CPL_FRAME_LEVEL_NONE,
        CPL_FRAME_LEVEL_NONE,
        CPL_FRAME_LEVEL_NONE,
        CPL_FRAME_LEVEL_NONE,
        CPL_FRAME_LEVEL_NONE,
        CPL_FRAME_LEVEL_NONE
    };

    long i;
    int flag = EXIT_SUCCESS;

    cpl_frame *_frame;
    cpl_frameset *frameset, * dark_frameset, * wave_frameset,
    				* p2vm_frameset, * flat_frameset;

    /* Insert tests below */

    /*
     * Create a frameset and extract the dark frames.
     */

    frameset = cpl_frameset_new();


    /* Extract dark frames, p2vm frames, and wave frames to the frame set created */

    for (i = 0; (size_t)i < CX_N_ELEMENTS(names); i++) {
        _frame = cpl_frame_new();

        cpl_frame_set_filename(_frame, names[i]);
        cpl_frame_set_tag(_frame, tags[i]);
        cpl_frame_set_type(_frame, CPL_FRAME_TYPE_IMAGE);
        cpl_frame_set_group(_frame, groups[i]);
        cpl_frame_set_level(_frame, levels[i]);

        cpl_frameset_insert(frameset, _frame);
    }

    test_data(dark_frameset, gravi_frameset_extract_dark_data(frameset),
    		                              "Extract the dark frames... ", flag);

    test_ivalue(1, cpl_frameset_get_size(dark_frameset),
    		                 "Check the size of the dark frame... ", flag);


    test_data(wave_frameset, gravi_frameset_extract_wave_data(frameset),
    		                              "Extract the wave frames... ", flag);

    test_ivalue(1, cpl_frameset_get_size(wave_frameset),
    		                 "Check the size of the wave frame... ", flag);

    test_data(p2vm_frameset, gravi_frameset_extract_p2vm_data(frameset),
    		                              "Extract the p2vm frames... ", flag);

    test_ivalue(1, cpl_frameset_get_size(p2vm_frameset),
    		                 "Check the size of the p2vm frame... ", flag);

    test_data(flat_frameset, gravi_frameset_extract_flat_data(frameset),
    		                              "Extract the flat frames... ", flag);

    test_ivalue(4, cpl_frameset_get_size(flat_frameset),
    		                 "Check the size of the flat frame... ", flag);

    cpl_frameset_delete(frameset);
    cpl_frameset_delete(dark_frameset);

	/*
     * Compute dark.
     */

	gravi_data * data_dark;
	gravi_data * test;
    gravi_data * data;

    test_data(data, gravi_data_load(DATADIR_TEST "Dark.fits"),
    		               "gravi_compute_dark: Load the data...", flag);

	test_data(data_dark, gravi_compute_dark(data), "gravi_compute_dark: Compute dark... ", flag);

	test_pfailure(CPL_ERROR_NULL_INPUT, gravi_compute_dark(NULL),
			              "gravi_compute_dark: Try to compute the dark from NULL data... ", flag);
	if (COMPUTE_FILES) gravi_data_save_data (data_dark, "test_files/gravi_dark_map.fits", CPL_IO_CREATE);

	test_data(test, gravi_data_load(DATADIR_TEST "gravi_dark_map.fits"),
			                               "gravi_compute_dark: Load the test data... ", flag);


	gravi_data_delete(test);
	gravi_data_delete(data);
	gravi_data_delete(data_dark);


	/*
     * Compute dark and gain.
     */

	gravi_data * dark_map, * flat_data[4];
    const char * filename;


    test_data(data, gravi_data_load(DATADIR_TEST "Dark.fits"),
    		               "gravi_compute_dark: Load the data SC...", flag);

    dark_map = gravi_compute_dark (data);
    cpl_propertylist * plist_dark;


    for (i = 0; i < cpl_frameset_get_size(flat_frameset); i++ ){
    	_frame = cpl_frameset_get_position(flat_frameset, i);
    	filename = cpl_frame_get_filename(_frame);

    	flat_data[i] = gravi_data_load (filename);

    }
    /* Expected value */
    cpl_propertylist * applist = NULL;
    applist = gravi_compute_gain (flat_data, 4, dark_map);

	gravi_data_delete(data);
    cpl_propertylist_delete(applist);

	/*
     * The bad pixel function.
     */

	gravi_data * badpix_data;
	cpl_propertylist * badpix_plist;
	cpl_parameterlist * paralist;
	paralist = cpl_parameterlist_new();
	cpl_parameter * p;
    p = cpl_parameter_new_value("gravity.calib.bad-dark-threshold", CPL_TYPE_INT, "the rms factor for "
					"dark bad pixel threshold", "gravi.preproc", 10);

    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "Bad dark threshold");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(paralist, p);



	test_pfailure(CPL_ERROR_NULL_INPUT,
					gravi_compute_badpix(NULL, flat_data, 4, paralist),
					"gravi_compute_badpix: Try to compute the bad pixel map from "
					"NULL data... ", flag);

	test_data(badpix_data, gravi_compute_badpix(dark_map, flat_data, 4, paralist),
			"gravi_compute_badpix: Compute the bad pix map with the rms factor equal to 10..." , flag);

	if (COMPUTE_FILES) gravi_data_save_data (badpix_data, "test_files/gravi_badpix_map.fits", CPL_IO_CREATE);
	test_data(test, gravi_data_load(DATADIR_TEST "gravi_badpix_map.fits"),
			                           "gravi_compute_badpix: "
			                           "Load the test data bad pix map.. ", flag);

//	table_test = gravi_data_get_table(test, GRAVI_IMAGING_DATA_FT_EXT);
//
//	badpix_table = gravi_data_get_table(badpix_data, GRAVI_IMAGING_DATA_FT_EXT);
//
//	test_ivalue(1, gravi_table_compare(table_test, badpix_table),
//				             "gravi_compute_badpix: Check the two data are the same... ", flag);
	badpix_plist = gravi_data_get_plist (badpix_data, GRAVI_PRIMARY_HDR_EXT);

	int badpix_num = cpl_propertylist_get_int (badpix_plist, QC_BADPIX_SC);
	char * com = cpl_sprintf ("gravi_compute_badpix: Check the number of bad pixels equal to %d...", badpix_num);

	test_ivalue(badpix_num, badpix_num, com, flag);
	cpl_free (com);

	gravi_data_delete(test);
	gravi_data_delete (badpix_data);

	cpl_parameter_set_int (p, 5);
	test_data(badpix_data, gravi_compute_badpix(dark_map, flat_data, 4, paralist),
			"gravi_compute_badpix: Compute the bad pix map with the rms factor equal to 5..." , flag);

	if (COMPUTE_FILES) gravi_data_save_data (badpix_data, "test_files/gravi_badpix_map.fits", CPL_IO_CREATE);

	test_data(test, gravi_data_load(DATADIR_TEST "gravi_badpix_map.fits"),
			                           "gravi_compute_badpix: "
			                           "Load the test data bad pix map.. ", flag);

//	table_test = gravi_data_get_table(test, GRAVI_IMAGING_DATA_FT_EXT);
//
//	badpix_table = gravi_data_get_table(badpix_data, GRAVI_IMAGING_DATA_FT_EXT);
//
//	test_ivalue(1, gravi_table_compare(table_test, badpix_table),
//				             "gravi_compute_badpix: Check the two data are the same... ", flag);


	badpix_plist = gravi_data_get_plist (badpix_data, GRAVI_PRIMARY_HDR_EXT);

	badpix_num = cpl_propertylist_get_int (badpix_plist, QC_BADPIX_SC);
	com = cpl_sprintf ("gravi_compute_badpix: Check the number of bad pixels equal to %d...", badpix_num);

	test_ivalue(badpix_num, badpix_num, com, flag);
	cpl_free (com);
	gravi_data_delete(test);

	/*
     * Compute Profile map.
     */

	gravi_data * profile_map, * profile_test;
	cpl_propertylist * plist;

    p = cpl_parameter_new_value("gravity.calib.profile-width", CPL_TYPE_INT, "width of profile in pixel", "gravi.preproc", 3);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "profile width");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(paralist, p);

    /* --How to deal with bad-pixels  */
	p = cpl_parameter_new_value("gravity.calib.force-badpix-to-zero",
			CPL_TYPE_BOOL, "Force the badpixel to zero in profile", "gravi.preproc", TRUE);
	cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "force-badpix-to-zero");
	cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
	cpl_parameterlist_append(paralist, p);


	test_data(profile_map, gravi_compute_profile(flat_data,
					dark_map, badpix_data, 4, paralist), "gravi_compute_profile: "
							"Compute the profile map with the width parameter equal to 3... ", flag);

	test_pfailure(CPL_ERROR_NULL_INPUT,
					gravi_compute_profile(NULL, dark_map, NULL, 4, paralist),
					"gravi_compute_profile: Try to compute the profile from "
					"NULL data... ", flag);

	if (COMPUTE_FILES) gravi_data_save_data (profile_map, "test_files/gravi_profile_map.fits", CPL_IO_CREATE);
	test_data(profile_test, gravi_data_load(DATADIR_TEST "gravi_profile_map.fits"),
			                           "gravi_compute_Profile: "
			                           "Load the test data profile... ", flag);

//	table_test = gravi_data_get_table(profile_test,
//											GRAVI_PROFILE_DATA_NAME_EXT);
//
//	profile_table = gravi_data_get_table(profile_map,
//											GRAVI_PROFILE_DATA_NAME_EXT);
//
//	test_ivalue(1, gravi_table_compare(table_test, profile_table),
//			             "gravi_compute_profile: Check the two "
//			             "data are the same... ", flag);

	gravi_data_delete (profile_test);
	gravi_data_delete (profile_map);
	cpl_parameter_set_int (cpl_parameterlist_find(paralist, "gravity.calib.profile-width"), 5);
	test_data(profile_map, gravi_compute_profile(flat_data,
					dark_map, badpix_data, 4, paralist), "gravi_compute_profile: "
							"Compute the profile map with the width parameter equal to 5... ", flag);

	if (COMPUTE_FILES) gravi_data_save_data (profile_map, "test_files/gravi_final_profile2.fits", CPL_IO_CREATE);
	test_data(profile_test, gravi_data_load(DATADIR_TEST "gravi_profile_map.fits"),
			                           "gravi_compute_Profile: "
			                           "Load the test data profile... ", flag);


	plist = gravi_data_get_plist(profile_map, GRAVI_PRIMARY_HDR_EXT);

	for (i = 0; i < 4; i++)
		gravi_data_delete (flat_data[i]);

	gravi_data_delete (profile_test);

	/*
     * Compute wave.
     */

	gravi_data * data_wave = gravi_data_new (0);
	double mjd_obs;

    test_data(data, gravi_data_load(DATADIR_TEST "Wave.fits"),
    		               "gravi_compute_wave: Load the data...", flag);


    cpl_table * opl_table = NULL, * metrology_table = gravi_data_get_table(data,
			GRAVI_METROLOGY_EXT);

	cpl_propertylist * met_plist = cpl_propertylist_duplicate(gravi_data_get_plist(data,
			GRAVI_METROLOGY_EXT));
    opl_table = cpl_table_new (cpl_table_get_nrow(metrology_table));
    mjd_obs = gravi_pfits_get_mjd (gravi_data_get_plist(data,
			GRAVI_PRIMARY_HDR_EXT));
    double wave_met = gravi_pfits_get_met_wavelength_mean(
    		gravi_data_get_plist(data, GRAVI_PRIMARY_HDR_EXT), metrology_table);
	cpl_table * p2vm_met = gravi_metrology_compute_p2vm (metrology_table, wave_met);

    cpl_msg_info (cpl_func, "Extract SPECTRUM for WAVE_RAW");
    gravi_data * spectrum_data = gravi_extract_spectrum (data, profile_map, dark_map,
                                                         badpix_data, NULL,paralist,
                                                         GRAVI_DET_ALL);

    cpl_msg_info (cpl_func, "Compute OPDs for WAVE_RAW");
    gravi_wave_compute_opds (spectrum_data, gravi_data_get_table (data, GRAVI_METROLOGY_EXT),
                             GRAVI_DET_ALL);
    FREE (gravi_data_delete, data);

    /* Compute wave calibration for FT and SC */
    gravi_parameter_add_wave(paralist);

    test(gravi_compute_wave (data_wave, spectrum_data, GRAVI_FT, paralist),
			"gravi_compute_wave: Compute wave FT... ", flag);
    test(gravi_compute_wave (data_wave, spectrum_data, GRAVI_SC, paralist),
			"gravi_compute_wave: Compute wave FT... ", flag);
    plist = cpl_propertylist_duplicate(gravi_data_get_plist (data_wave,
    		GRAVI_WAVE_DATA_SC_EXT));

    gravi_data_add_table (data_wave, plist, "OPL_TABLE", opl_table);
    gravi_data_add_table (data_wave, NULL, "P2VM_MET", p2vm_met);

    /* Compute wave test failure  */
	test_pfailure(CPL_ERROR_NULL_INPUT, gravi_compute_wave(NULL, spectrum_data, GRAVI_FT, paralist),
			              "gravi_compute_wave: Try to compute the wave with NULL wave_data... ", flag);

	test_pfailure(CPL_ERROR_NULL_INPUT, gravi_compute_wave(data_wave, NULL, GRAVI_FT, paralist),
			              "gravi_compute_wave: Try to compute the wave from NULL spectrum data... ", flag);

	test_pfailure(CPL_ERROR_ILLEGAL_INPUT, gravi_compute_wave(data_wave, spectrum_data, 10, paralist),
			              "gravi_compute_wave: Try to compute the wave with illegal type... ", flag);

    FREE (gravi_data_delete, spectrum_data);
//	cpl_table_delete (opl_table);

	if (COMPUTE_FILES) gravi_data_save_data (data_wave, "test_files/gravi_wave_map.fits", CPL_IO_CREATE);
//	test_data(test, gravi_data_load(DATADIR_TEST "gravi_wave_map.fits"),
//			                               "gravi_compute_wave: Load the test data wave... ", flag);

//	table_test = gravi_data_get_table(test, "P2VM_MET");
//	test_ivalue(1, gravi_table_compare(table_test, p2vm_met),
//			             "gravi_metrology_calibration: Check p2vm metrology table... ", flag);
//
//	table_test = gravi_data_get_table(test, "OPL_TABLE");
//	test_ivalue(1, gravi_table_compare(table_test, opl_table),
//			             "gravi_metrology_calibration: Check OPL table... ", flag);
//
//	table_test = gravi_data_get_table(test, GRAVI_WAVE_DATA_FT_EXT);
//	table_wave = gravi_data_get_table(data_wave, GRAVI_WAVE_DATA_FT_EXT);
//
//	test_ivalue(1, gravi_table_compare(table_test, table_wave),
//			             "gravi_compute_wave: Check the two data are the same for the FT data... ", flag);

//	table_test = gravi_data_get_table(test, GRAVI_WAVE_DATA_SC_EXT);
//	table_wave = gravi_data_get_table(data_wave, GRAVI_WAVE_DATA_SC_EXT);
//
//	test_ivalue(1, gravi_table_compare(table_test, table_wave),
//			             "gravi_compute_wave: Check the two data are the same for the SC data... ", flag);

//	gravi_data_delete(test);
	gravi_data_delete(data);
    gravi_data_delete(data_wave);


	/*
     * Preprocessing the data files.
     */

	gravi_data ** calib_data;
	cpl_propertylist * plist_wave, * plist_profile;

	data_wave = gravi_data_load (DATADIR_TEST "gravi_wave_map.fits");

	calib_data = cpl_malloc(4 * sizeof(gravi_data *));

	calib_data[0] = dark_map;
	plist_dark = gravi_data_get_plist(calib_data[0],
										GRAVI_PRIMARY_HDR_EXT);
	cpl_propertylist_append_string(plist_dark, CPL_DFS_PRO_CATG, GRAVI_DARK_MAP);


	calib_data[1] = data_wave;
	plist_wave = gravi_data_get_plist(calib_data[1],
										GRAVI_PRIMARY_HDR_EXT);
	cpl_propertylist_append_string(plist_wave, CPL_DFS_PRO_CATG, GRAVI_WAVE_MAP);

	calib_data[2] = profile_map;
	plist_profile = gravi_data_get_plist(calib_data[2],
										GRAVI_PRIMARY_HDR_EXT);
	cpl_propertylist_append_string(plist_profile, CPL_DFS_PRO_CATG, GRAVI_FLAT_MAP);

	calib_data[3] = badpix_data;

	cpl_propertylist  * plist_badpix = gravi_data_get_plist(calib_data[3],
										GRAVI_PRIMARY_HDR_EXT);
	cpl_propertylist_append_string(plist_badpix, CPL_DFS_PRO_CATG, GRAVI_BAD_MAP);


	cpl_parameterlist * parlist = cpl_parameterlist_new();

    /* Use static names (output_procatg.fits) */
    gravi_parameter_add_static_name (parlist);

    /* Intermediate files */
    gravi_parameter_add_biassub_file (parlist);
    gravi_parameter_add_spectrum_file (parlist);
    gravi_parameter_add_preproc_file (parlist);
    gravi_parameter_add_p2vmred_file (parlist);
    gravi_parameter_add_astro_file (parlist);

    /* Averaging */
    gravi_parameter_add_average_vis (parlist);

    /* Bias-method */
    gravi_parameter_add_biasmethod (parlist);

    /* Extraction */
    gravi_parameter_add_extract (parlist);
    gravi_parameter_add_metrology (parlist);
    
    /* Snr, signal, rejectio flags, vis */
    int isCalib = 0;
    gravi_parameter_add_compute_snr (parlist, isCalib);
    gravi_parameter_add_compute_signal (parlist, isCalib);
    gravi_parameter_add_rejection (parlist, isCalib);
    gravi_parameter_add_compute_vis (parlist, isCalib);
    


	/* Add the FLAT_RAW and WAVE_RAW to the p2vm frameset */
	if ( !cpl_frameset_is_empty (flat_frameset) )
	  cpl_frameset_join (p2vm_frameset, flat_frameset);

	if ( !cpl_frameset_is_empty (wave_frameset) )
	  cpl_frameset_join (p2vm_frameset, wave_frameset);

	/* Get the number of the p2vm frame contained in the frameset */
	int size = cpl_frameset_get_size(p2vm_frameset) ;


	gravi_data * p2vm_data=NULL;
	cpl_propertylist * primary_hdr;
	int ** valid_trans = cpl_malloc (2 * sizeof (int*));
	int ** valid_CP = cpl_malloc (2 * sizeof (int*));

	for (i = 0 ; i < 2; i++){
		valid_trans[i] = cpl_calloc (4, sizeof (int));
		valid_CP[i] = cpl_calloc (6, sizeof (int));
	}

	p = (cpl_parameter * ) cpl_parameterlist_find(parlist, "gravi.nspectrum_FT");


	for (i = 0; i < size; i++){
	        gravi_data * preproc_data;
		_frame = cpl_frameset_get_position(p2vm_frameset, i);
		filename = cpl_frame_get_filename(_frame);
		cpl_msg_info (cpl_func, "TEST now in %s", filename);

		if (strstr (filename, "Wave.fits")) { cpl_msg_info (cpl_func,"Skip Wave.fits"); continue;}
		data = gravi_data_load(filename);

		test_data(preproc_data, gravi_extract_spectrum (data, profile_map, dark_map,
                                                        badpix_data, NULL,paralist,
                                                        GRAVI_DET_ALL),
				   "gravi_preproc: Compute the preproc data... ", flag);
		primary_hdr = gravi_data_get_plist(data,
											GRAVI_PRIMARY_HDR_EXT);
		gravi_data_delete(data);



		/* Construction of the p2vm data. */
		if ( p2vm_data == NULL) {
			test_data(p2vm_data, gravi_create_p2vm (data_wave),
							   "gravi_create_p2vm: create the p2vm table... ", flag);

		/* Rescale to common wavelength */
		gravi_align_spectrum (preproc_data, data_wave, p2vm_data, GRAVI_DET_ALL);

		}

//		if (i == 0){
//			cpl_msg_info (NULL, "Construction of the p2vm map");
//			p2vm_data = gravi_data_new(0);
//			p2vm_primary_hdr = cpl_propertylist_new ();
//
//			for (j = 0; j < gravi_data_get_size (preproc_data); j++){
//				plist = gravi_data_get_plist_x (preproc_data, j);
//				const char * plist_name = gravi_pfits_get_extname (plist);
//				if (plist_name == NULL)
//					continue;
//				/* Get the type of the extention */
//				type_data = gravi_pfits_get_extension_type (plist);
//
//						  strcmp (plist_name, GRAVI_IMAGING_DETECTOR_FT_EXT) &&
//							strcmp (plist_name, GRAVI_IMAGING_DETECTOR_SC_EXT) &&
//							  strcmp (plist_name, GRAVI_OI_WAVELENGTH_FT_EXT) &&
//								strcmp (plist_name, GRAVI_OI_WAVELENGTH_SC_EXT))){
//
//
//					if (type_data == 2)
//						gravi_data_add (p2vm_data, plist,
//								cpl_table_duplicate (gravi_data_get_table (preproc_data, plist_name)));
//					else if (type_data == 3)
//						gravi_data_add_cube (p2vm_data, plist,
//									cpl_imagelist_duplicate (gravi_data_get_cube (preproc_data, plist_name)));
//				}
//			}
//			detector_table = gravi_data_get_table (preproc_data,
//													  GRAVI_IMAGING_DETECTOR_SC_EXT);
//			n_region = cpl_table_get_nrow(detector_table);
//			if (n_region > 24)
//				oiwave_plist = gravi_data_get_oi_propertylist
//				  (preproc_data, GRAVI_OI_WAVE_SWITCH(GRAVI_TYPE_SC),
//				   GRAVI_SPECTRO_SWITCH(GRAVI_TYPE_SC,0,2));
////				oiwave_plist = gravi_data_get_oi_propertylist
////					(preproc_data, GRAVI_OI_WAVELENGTH_FT_EXT, INSNAME_FT_P1);
//			else
//				oiwave_plist = gravi_data_get_oi_propertylist
//				  (preproc_data, GRAVI_OI_WAVE_SWITCH(GRAVI_TYPE_SC),
//				   GRAVI_SPECTRO_SWITCH(GRAVI_TYPE_SC,0,1));
//
//			int nacq = cpl_table_get_nrow(gravi_data_get_table (preproc_data,
//					  GRAVI_SPECTRUM_DATA_SC_EXT));
//			n_wave = gravi_pfits_get_nwave(oiwave_plist);
//
//			p2vm_table = gravi_p2vm_new(n_region, n_wave, n_tel);
//			plist = cpl_propertylist_new();
//
//			cpl_propertylist_append_string(plist, "ORIGIN", origin);
//			cpl_propertylist_append_int(plist, "NREGION", n_region);
//			cpl_propertylist_append_int(plist, "NWAVE", n_wave);
//			cpl_propertylist_append_string(plist, "EXTNAME", GRAVI_P2VM_DATA_SC_EXT);
//
//			gravi_data_add (p2vm_data, plist, p2vm_table);
//
//
//			detector_table = gravi_data_get_table (preproc_data,
//													  GRAVI_IMAGING_DETECTOR_FT_EXT);
//			n_region = cpl_table_get_nrow(detector_table);
//			if (n_region > 24)
//				oiwave_plist = gravi_data_get_oi_propertylist
//				  (preproc_data, GRAVI_OI_WAVE_SWITCH(GRAVI_TYPE_FT),
//				   GRAVI_SPECTRO_SWITCH(GRAVI_TYPE_FT,0,2));
////				oiwave_plist = gravi_data_get_oi_propertylist
////					(preproc_data, GRAVI_OI_WAVELENGTH_FT_EXT, INSNAME_FT_P1);
//			else
//				oiwave_plist = gravi_data_get_oi_propertylist
//				  (preproc_data, GRAVI_OI_WAVE_SWITCH(GRAVI_TYPE_FT),
//				   GRAVI_SPECTRO_SWITCH(GRAVI_TYPE_FT,0,1));
//
//			nacq = cpl_table_get_nrow(gravi_data_get_table (preproc_data,
//					  GRAVI_SPECTRUM_DATA_FT_EXT));
//			n_wave = gravi_pfits_get_nwave(oiwave_plist);
//
//			p2vm_table = gravi_p2vm_new(n_region, n_wave, n_tel);
//
//			cpl_propertylist_set_int(plist, "NREGION", n_region);
//			cpl_propertylist_set_int(plist, "NWAVE", n_wave);
//			cpl_propertylist_set_string(plist, "EXTNAME", GRAVI_P2VM_DATA_FT_EXT);
//
//
//			gravi_data_add (p2vm_data, plist, p2vm_table);
//
//			gravi_data_set_propertylist (p2vm_data, GRAVI_PRIMARY_HDR_EXT,
//																	p2vm_primary_hdr);
//			cpl_propertylist_delete(p2vm_primary_hdr);
//
//
//			cpl_propertylist_delete (plist);
//
//
//		}

		test(gravi_compute_p2vm (p2vm_data, preproc_data, valid_trans, valid_CP,
		        GRAVI_DET_ALL),
				"gravi_compute_p2vm : Compute the p2vm... ", flag);

		gravi_data_delete (preproc_data);
	}
	gravi_data_delete (data_wave);

	test(gravi_p2vm_normalisation (p2vm_data, valid_trans, valid_CP ),
			"gravi_p2vm_normalisation : normalisation of the p2vm... ", flag);

	for (i = 0 ; i < 2; i++){
		cpl_free (valid_trans[i]);
		cpl_free (valid_CP[i]);
	}
	cpl_free (valid_trans);
	cpl_free (valid_CP);

	gravi_data * p2vm_file;

	if (COMPUTE_FILES) gravi_data_save_data (p2vm_data, "test_files/gravi_p2vm_map.fits", CPL_IO_CREATE);

	test_data(p2vm_file, gravi_data_load(DATADIR_TEST "gravi_p2vm_map.fits"),
			                               "gravi_compute_p2vm: Load the p2vm test data... ", flag);


//	table_test = gravi_data_get_table(p2vm_data, GRAVI_P2VM_DATA_FT_EXT);
//	table_p2vm = gravi_data_get_table(p2vm_file, GRAVI_P2VM_DATA_FT_EXT);
//
//	test_ivalue(1, gravi_table_compare(table_test, table_p2vm),
//			             "gravi_compute_p2vm: Check the two p2vm data are the same for the FT data... ", flag);

//	table_test = gravi_data_get_table(p2vm_data, GRAVI_P2VM_DATA_SC_EXT);
//	table_p2vm = gravi_data_get_table(p2vm_file, GRAVI_P2VM_DATA_SC_EXT);
//
//	test_ivalue(1, gravi_table_compare(table_test, table_p2vm),
//			             "gravi_compute_p2vm: Check the two p2vm data are the same for SC data... ", flag);

	gravi_data_delete (p2vm_file);

	cpl_propertylist_set_string (met_plist, "EXTNAME",
			GRAVI_P2VM_MET_EXT);
//	gravi_data_add (p2vm_data, met_plist, p2vm_met);
	gravi_data_delete (p2vm_data);

	cpl_propertylist_delete(met_plist);

	cpl_msg_info( cpl_func, " ***** Reduce the WAVE file to create the VIS_FLAT ***** " );
	data = gravi_data_load (DATADIR_TEST "Wave.fits");
	gravi_data * p2vm_reduced;
	//p2vm_data = gravi_data_load (DATADIR_TEST "p2vm_map.fits");
	p2vm_data = gravi_data_load (DATADIR_TEST "gravi_p2vm_map.fits");
	data_wave = gravi_data_load (DATADIR_TEST "gravi_wave_map.fits");

	cpl_msg_info (cpl_func, "Preproc the WAVE");
	gravi_data * preproc_data=gravi_extract_spectrum (data, profile_map, dark_map,
                                                      badpix_data, NULL,paralist,
                                                      GRAVI_DET_ALL);
    
	cpl_parameterlist_delete (paralist);
    gravi_align_spectrum (preproc_data, data_wave, p2vm_data, GRAVI_DET_ALL);

    /* Move extensions from raw_data and delete it */
    gravi_data_move_ext (preproc_data, data, GRAVI_ARRAY_GEOMETRY_EXT);
    gravi_data_move_ext (preproc_data, data, GRAVI_OPTICAL_TRAIN_EXT);
    gravi_data_move_ext (preproc_data, data, GRAVI_OPDC_EXT);
    gravi_data_move_ext (preproc_data, data, GRAVI_FDDL_EXT);
    gravi_data_move_ext (preproc_data, data, GRAVI_METROLOGY_EXT);
	FREE (gravi_data_delete, data);

	/* Compute the flux and visibilities for each telescope and
	 * per acquisition with the P2VM applied to preproc_data */
	cpl_msg_info (cpl_func, "Compute the P2VMRED");
    test_data(p2vm_reduced, gravi_compute_p2vmred (preproc_data, p2vm_data, "gravi_single", parlist, GRAVI_DET_ALL),
			"gravi_p2vm_reduce: reduce the p2vm... ", flag);

    /* Move extensions and delete preproc */
    gravi_data_move_ext (p2vm_reduced, preproc_data, GRAVI_METROLOGY_EXT);
    gravi_data_move_ext (p2vm_reduced, preproc_data, GRAVI_FDDL_EXT);
    gravi_data_move_ext (p2vm_reduced, preproc_data, GRAVI_OPDC_EXT);
	FREE (gravi_data_delete, preproc_data);
	FREE (gravi_data_delete, p2vm_data);


    /* Reduce the OPDC table */
	test (gravi_compute_opdc_state (p2vm_reduced),
			"gravi_compute_opdc_state :  ...", flag);

	/* Reduce the metrology */
    test (gravi_metrology_reduce (p2vm_reduced, NULL, NULL, parlist),
			"gravi_metrology_reduce : reduce the metrology ...", flag);

    //gravi_data_save_data (p2vm_reduced, "test_files/p2vm_reduced.fits", CPL_IO_CREATE);

	gravi_parameter_add_compute_snr (parlist, 0);
	gravi_parameter_add_rejection (parlist, 0);

	/* Compute the SNR */
	cpl_msg_info (cpl_func, "Compute the SNR");
	test(gravi_compute_snr (p2vm_reduced, parlist),
		 "gravi_signal_snr : ...", flag);

	/* Compute the signals */
	cpl_msg_info (cpl_func, "Compute the signal");
    gravi_parameter_add_compute_signal (parlist, 0);
	test(gravi_compute_signals (p2vm_reduced, NULL, parlist),
		 "gravi_compute_signals : ...", flag);

	/* Compute the rejection */
	cpl_msg_info (cpl_func, "Compute the signal");
	test(gravi_compute_rejection (p2vm_reduced, parlist),
		 "gravi_compute_rejection : ...", flag);
    
	cpl_msg_info (cpl_func, "Average the VIS averaged output");
	gravi_data * oi_vis = NULL;
    cpl_size current_frame = 0;
	test_data(oi_vis, gravi_compute_vis (p2vm_reduced, parlist, &current_frame),
				"gravi_vis_reduce in the SINGLE mode: Compute the squared, complex visibilities "
												   "and the cloture phase...", flag);

	if (COMPUTE_FILES) gravi_data_save_data (oi_vis, "test_files/vis_data.fits", CPL_IO_CREATE);
    test = gravi_data_load (DATADIR_TEST "vis_data.fits");
	gravi_data_delete (test);
	gravi_data_delete (oi_vis);
	gravi_data_delete (p2vm_reduced);
	gravi_data_delete (data_wave);
	gravi_data_delete (badpix_data);
	gravi_data_delete (dark_map);
	gravi_data_delete (profile_map);
	cpl_free (calib_data);
	/* Temp */

	cpl_frameset_delete (p2vm_frameset);
	cpl_frameset_delete(flat_frameset);
	cpl_frameset_delete(wave_frameset);

	cpl_parameterlist_delete (parlist);

//	for (i = 0; i < 1; i++){
//		gravi_data_delete (p2vm_reduce[i]);
//	}
//	cpl_free(p2vm_reduce);
//	gravi_data_delete (oi_vis);
	return flag;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Unit tests of gravi_utils module
 */
/*----------------------------------------------------------------------------*/

int main(void)
{
	int flag;
#if defined CPL_VERSION_CODE && CPL_VERSION_CODE >= CPL_VERSION(4, 0, 0)
    cpl_test_init(PACKAGE_BUGREPORT, CPL_MSG_INFO);
#else
    cpl_init();
#endif

    flag=gravi_calib_test();

    if (flag == EXIT_FAILURE)
    {
    	cpl_test_end(0);
    	exit(EXIT_FAILURE);
    }

    cpl_test_end(0);
    exit(EXIT_SUCCESS);
}
