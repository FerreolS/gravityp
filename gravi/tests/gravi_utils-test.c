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
#include "gravi_wave.h"
#include "gravi_cpl.h"
#include "gravi_preproc.h"
#include "gravi_calib.h"
#include "gravi_p2vm.h"
#include "gravi_dfs.h"
#include "gravi-test.h"


/*-----------------------------------------------------------------------------
                              PREPROC Private prototypes
 -----------------------------------------------------------------------------*/
cpl_table * gravi_table_ft_format (cpl_table * table_ft, cpl_table * sky_table_std,
				cpl_table * sky_table_avg, cpl_table * badpix, int n_region, double gain);

cpl_table * gravi_imglist_sc_collapse (cpl_table * profile_table,
                                       cpl_imagelist * raw_imglist,
                                       cpl_imagelist * rawVar_imglist,
                                       cpl_size startx);

cpl_error_code gravi_interpolate_spectrum_table (cpl_table * spectrum_table,
                                                 cpl_table * wave_table,
                                                 cpl_table * oiwave_table,
                                                 cpl_table * specflat_table);


//#include "gravi_data.c"
//#include "gravi_pfits.c"
//#include "gravi_utils.c"
//#include "gravi_wave.c"
//#include "gravi_cpl.c"
//#include "gravi_preproc.c"
//#include "gravi_calib.c"
//#include "gravi_p2vm.c"
//#include "gravi_dfs.c"
//#include "gravi-test.c"


#define DATADIR_TEST DATADIR
//#define DATADIR_TEST ""

int gravi_utils_test(void){

    const char *names[] = {
        "flat1.fits",
        "flat2.fits",
        "flat3.fits",
        "bias1.fits",
        "bias2.fits",
        "bias3.fits",
        "mbias.fits",
        "mflat.fits",
        "science.fits",
        "product.fits"
    };

    const char *tags[] = {
        "FLAT",
        GRAVI_DARK_RAW,
        "FLAT",
        "BIAS",
        GRAVI_DARK_RAW,
        GRAVI_DARK_RAW,
        "MASTER_BIAS",
        GRAVI_DARK_RAW,
        "SCIENCE",
        "SCIENCE_CALIBRATED"
    };

    cpl_frame_group groups[] = {
        CPL_FRAME_GROUP_RAW,
        CPL_FRAME_GROUP_RAW,
        CPL_FRAME_GROUP_RAW,
        CPL_FRAME_GROUP_RAW,
        CPL_FRAME_GROUP_RAW,
        CPL_FRAME_GROUP_RAW,
        CPL_FRAME_GROUP_CALIB,
        CPL_FRAME_GROUP_CALIB,
        CPL_FRAME_GROUP_RAW,
        CPL_FRAME_GROUP_PRODUCT
    };

    cpl_frame_level levels[] = {
        CPL_FRAME_LEVEL_NONE,
        CPL_FRAME_LEVEL_NONE,
        CPL_FRAME_LEVEL_NONE,
        CPL_FRAME_LEVEL_NONE,
        CPL_FRAME_LEVEL_NONE,
        CPL_FRAME_LEVEL_NONE,
        CPL_FRAME_LEVEL_NONE,
        CPL_FRAME_LEVEL_NONE,
        CPL_FRAME_LEVEL_NONE,
        CPL_FRAME_LEVEL_FINAL
    };

    long i;
    int flag=EXIT_SUCCESS;

    cpl_frame *_frame;
    cpl_frameset *frameset, * dark_frameset;

//    cpl_test_init(PACKAGE_BUGREPORT, CPL_MSG_WARNING);

    /* Insert tests below */

    /*
     * Create a frameset and extract the dark frames.
     */

    frameset = cpl_frameset_new();


    /* Add frames to the frame set created */

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

    test_ivalue(4, cpl_frameset_get_size(dark_frameset),
    		                 "Check the size of the dark frame... ", flag);

    cpl_frameset_delete(frameset);
    cpl_frameset_delete(dark_frameset);

    /*
     * Create a an image list from a table of data.
     */

    cpl_imagelist *img_list;
    cpl_table * imaging_data;
    cpl_propertylist * detector_plist;
    gravi_data * data;

    test_data(data, gravi_data_load(DATADIR_TEST "Wave.fits"),
    		               "Load the data...", flag);

    test_data(imaging_data, gravi_data_get_table(data,
    		GRAVI_IMAGING_DATA_FT_EXT), "Get Imaging data table...", flag);

    test_data(detector_plist, gravi_data_get_plist(data,
    	                                     GRAVI_IMAGING_DETECTOR_FT_EXT),
    	                       "Get imaging detector property list... ", flag);

    test_data(img_list, gravi_imagelist_wrap_column(imaging_data,
    		       "PIX"), "Create image list from PIX column... ", flag);
    
    test_pfailure(CPL_ERROR_ILLEGAL_INPUT, gravi_imagelist_wrap_column
    		(imaging_data, "NO_COLUMN"), "Unexisting column", flag);

    test_pfailure(CPL_ERROR_NULL_INPUT, gravi_imagelist_wrap_column(NULL, 0),
    		"Try to create an image list from NULL table... ", flag);



    /*
     * Transform an image to an array.
     */

    cpl_image * img;
    cpl_array * ar;

    test_data(img, cpl_imagelist_get(img_list, 0), "Get an image from "
    		"the image list... ", flag);

	test_data(ar, gravi_array_wrap_image(img), "Get the array from the image... ",
			flag);

	int ar_size = cpl_array_get_size(ar);
	test_ivalue(ar_size, cpl_image_get_size_x(img) * cpl_image_get_size_y(img),
			              "Check the sizes of the two elements... ", flag);

	test_pfailure(CPL_ERROR_NULL_INPUT, gravi_array_wrap_image(NULL),
			                     "Try to create an array from NULL image... ", flag);

	//cpl_imagelist_delete(img_list);
	cpl_array_unwrap(ar);

	/*
	 * Test unit gravi_shutters_check
	 */

	test_pfailure (CPL_ERROR_ILLEGAL_INPUT, gravi_check_shutter (NULL, 1, 1, 1, 1), "gravi_shutters_check : "
			"try the check shutters of a NULL data ...", flag);

	//int shutters;
	test_ivalue (1 , gravi_check_shutter (gravi_data_get_header (data), 1, 1, 1, 1),
			"gravi_shutters_check : check the shutters of an wave data ...", flag);

	//cpl_free (shutters);


	/*
	 * gravi_table_oi_create
	 */
	cpl_table * oi_table;
	test_pfailure (CPL_ERROR_ILLEGAL_INPUT, gravi_table_oi_create(10, 8, "BLABLABLA"), "gravi_table_oi_create : "
			"Try to create and oi table with an illegal name ...", flag);

	test_data (oi_table, gravi_table_oi_create(10, 8, GRAVI_OI_T3_EXT),
			"gravi_table_oi_create : Create an GRAVI_OI_T3_EXT table ...", flag);
	cpl_table_delete (oi_table);

	/*
	 * load test tables
	 */
	gravi_data * dark_map = gravi_data_load (DATADIR_TEST "gravi_dark_map.fits");
	gravi_data * badpix = gravi_data_load (DATADIR_TEST "gravi_badpix_map.fits");
	gravi_data * profile_map = gravi_data_load (DATADIR_TEST "gravi_profile_map.fits");

	/*
	 * gravi_plist_get_qc
	 */

	test_pfailure (CPL_ERROR_NULL_INPUT, gravi_plist_get_qc (NULL), "gravi_plist_get_qc : "
			"try the extract the QC parameters of a NULL data ...", flag);

	cpl_propertylist * qcPlist;

	test_data (qcPlist, gravi_plist_get_qc (gravi_data_get_header (dark_map)), "gravi_plist_get_qc : "
			"Extract the QC parameters of a master dark ...", flag);

	cpl_propertylist_delete (qcPlist);



	/*
	 * gravi_data_remove_badpixel
	 */
    cpl_imagelist * raw_imglist;
    raw_imglist = cpl_imagelist_duplicate (gravi_data_get_cube(data, GRAVI_IMAGING_DATA_SC_EXT));

	test_pfailure (CPL_ERROR_NULL_INPUT, gravi_remove_badpixel_sc(NULL, gravi_data_get_img (badpix, GRAVI_IMAGING_DATA_SC_EXT)),
			"gravi_data_remove_badpixel : "
			"try to remove bad pixels from the input NULL data...", flag);

	test(gravi_remove_badpixel_sc(raw_imglist, gravi_data_get_img (badpix, GRAVI_IMAGING_DATA_SC_EXT)),
			"gravi_data_remove_badpixel : "
			"Remove bad pixels from the input data SC and FT ...", flag);
    cpl_imagelist_delete(raw_imglist);

    /*
	 * gravi_imglist_sc_collapse
	 */
	cpl_propertylist * profile_plist=gravi_data_get_plist(profile_map, GRAVI_PROFILE_DATA_EXT);
	cpl_size startx = 22;//gravi_pfits_get_startx (profile_plist);
	test_pfailure (CPL_ERROR_NULL_INPUT, gravi_imglist_sc_collapse(
			NULL,
			gravi_data_get_cube(data, GRAVI_IMAGING_DATA_SC_EXT),
			gravi_data_get_cube (data, GRAVI_IMAGING_DATA_SC_EXT), startx),
			"gravi_imglist_sc_collapse : "
			"Try to extract the correct format of the images SC using wrong startx ...", flag);

	cpl_table * tableNew;
	test_data (tableNew, gravi_imglist_sc_collapse(
			gravi_data_get_table(profile_map, GRAVI_PROFILE_DATA_EXT),
			gravi_data_get_cube(data, GRAVI_IMAGING_DATA_SC_EXT),
			gravi_data_get_cube (data, GRAVI_IMAGING_DATA_SC_EXT), startx),
			"gravi_imglist_sc_collapse : Compute the correct format of the data SC ...", flag);
	cpl_table_delete (tableNew);


	/*
	 * Test unit of the function that extract spectrum
	 */
	gravi_data * spectrum_data;

	test_pfailure(CPL_ERROR_NULL_INPUT, gravi_extract_spectrum(NULL,
                                                               NULL,
                                                               NULL,  NULL,
                                                               NULL, NULL),
			              "gravi_extract_spectrum: Try to extract the spectrum with a NULL data... ", flag);

    cpl_parameterlist * parlist = cpl_parameterlist_new ();

	test_data (spectrum_data, gravi_extract_spectrum(data,
                                                     profile_map,
                                                     dark_map,  badpix,
                                                     NULL, parlist), "gravi_extract_spectrum: extract the spectrum ...", flag);

    cpl_parameterlist_delete (parlist);




	/*
	 * gravi_create_p2vm_table
	 */

	cpl_table * p2vm;
	cpl_table * detector_table = gravi_data_get_table(profile_map, GRAVI_IMAGING_DETECTOR_SC_EXT);
	test_data (p2vm, gravi_create_p2vm_table (detector_table, 10), "gravi_p2vm_new : "
			"Create a p2vm table with 10 waves canal ...", flag);
	cpl_table_delete (p2vm);

	/*
	 * gravi_table_get_vector
	 */
	char * data_x = "DATA10";
	cpl_vector * vector;
	test_pfailure (CPL_ERROR_NULL_INPUT, gravi_table_get_vector (NULL, 5,
            data_x), "gravi_table_get_vector : Try to get a vector of a column DATA10 from a NULL data ...", flag);

	test_pfailure (CPL_ERROR_ILLEGAL_INPUT, gravi_table_get_vector (gravi_data_get_table (spectrum_data,
			GRAVI_SPECTRUM_DATA_SC_EXT), 5,
            "DATA60"), "gravi_table_get_vector : Try to get a vector of a column who does not exist ...", flag);

	test_data ( vector, gravi_table_get_vector (gravi_data_get_table (spectrum_data,
			GRAVI_SPECTRUM_DATA_SC_EXT), 5,
            data_x), "gravi_table_get_vector : get the vector of the wave cal 5 from the column DATA10 ...", flag);
	cpl_vector_delete (vector);


	/*
	 * gravi_table_ft_format_bis
	 */
	cpl_table * table_ft = gravi_data_get_table (data, GRAVI_IMAGING_DATA_FT_EXT), * tableNew_ft;
	cpl_table * dark_ft = gravi_data_get_table (dark_map, GRAVI_IMAGING_DATA_FT_EXT);
	cpl_table * darkStd_ft = gravi_data_get_table (dark_map, GRAVI_IMAGING_ERR_FT_EXT);
	cpl_table * badpix_ft = gravi_data_get_table (badpix, GRAVI_IMAGING_DATA_FT_EXT);

	test_data (tableNew_ft, gravi_table_ft_format(table_ft, darkStd_ft, dark_ft, badpix_ft, 24, 25),
			"gravi_table_ft_format_bis : "
			"Remove the dark and extract the spectrum FT ...", flag);
	FREE(cpl_table_delete, tableNew_ft);

	FREE(gravi_data_delete, spectrum_data);
	FREE(gravi_data_delete, badpix);
	FREE(gravi_data_delete, dark_map);
	FREE(gravi_data_delete, profile_map);
	FREE(gravi_data_delete, data);

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
    cpl_test_init(PACKAGE_BUGREPORT, CPL_MSG_WARNING);
#else
    cpl_init();
#endif

    flag=gravi_utils_test();

    if (flag == EXIT_FAILURE)
    {
    	cpl_test_end(0);
    	exit(EXIT_FAILURE);
    }
    cpl_test_end(0);
    exit(EXIT_SUCCESS);
}
