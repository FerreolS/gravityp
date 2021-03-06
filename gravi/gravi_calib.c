/* $Id: gravi_calib.c,v 1.10 2012/03/23 15:10:40 nazouaoui Exp $
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

/**
 * @defgroup gravi_calib  Detector calibration
 *
 * This module contains the functions that are used to calibrate the detector.
 * These functions characterize the pixels of the detectors.
 * They are call by the recipe gravity_dark or gravity_p2vm.
 * - function dedicated to dark computation @c gravi_compute_dark(),
 * @c gravi_average_dark()
 * - function dedicated to flat computation : @c gravi_compute_profile(),
 * @c gravi_fit_profile(), @c gravi_create_profile_image()
 * - computation of the gain : @c gravi_compute_gain()
 * - computation of the bad pixels : @c gravi_compute_badpix()
 */
/**@{*/

/*
 * History :
 *               11/01/2019 : Move global parameter from gravi_calib.h
 *                                initialize sigma with 0
 *                                add brackets to 'for' statements
 */
/*----------------------------------------------------------------------------
                                    DEBUG
 -----------------------------------------------------------------------------*/

#define INFO_DEBUG 0
/*-----------------------------------------------------------------------------
                                   Includes
 -----------------------------------------------------------------------------*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cpl.h>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <complex.h>

#include "gravi_data.h"
#include "gravi_dfs.h"
#include "gravi_pfits.h"
#include "gravi_cpl.h"

#include "gravi_utils.h"

#include "gravi_preproc.h"
#include "gravi_calib.h"

/*-----------------------------------------------------------------------------
              ekw 11/01/2019 : Global parameter from gravi_calib.h
 -----------------------------------------------------------------------------*/

static double met_Sep_2016[64] =    {-0.000338233    ,-0.000282107    ,0.000345582    ,0.000381979    ,0.000326764    ,-0.000225664    ,5.35555e-05    ,0.00060938    ,-0.000205619    ,-0.000265577    ,0.00017356    ,0.000265687    ,0.000295987    ,3.05255e-05    ,-0.000498407    ,0.000232075    ,0.00502058    ,0.002045925    ,-0.000110657    ,-0.000403592    ,-0.000065043    ,-0.000433645    ,-0.000626545    ,3.29765e-05    ,-0.001694525    ,-0.0011684    ,6.2795e-06    ,-0.001584115    ,-0.000735553    ,-0.000868538    ,-0.000985087    ,-0.0008204    ,0.000768577    ,0.000848342    ,-0.000134943    ,-0.000385157    ,-8.68822e-05    ,-0.000757366    ,0.000446051    ,-0.000231723    ,0.000790425    ,-0.000638897    ,0.000503496    ,-7.78205e-05    ,0.000287366    ,0.000243789    ,0.000083288    ,-0.000125138    ,-0.000147337    ,4.14224e-05    ,0.000123082    ,-0.00117179    ,4.54785e-05    ,-0.000186707    ,0.000682836    ,0.00090649    ,0.000357256    ,-0.002133845    ,-0.00151895    ,-0.00150048    ,-0.00266423    ,-0.0030716    ,0.000599228    ,-0.001078583};
static double met_Mar_2017[64] =    {0.000481474    ,0.000678173    ,0.00022232    ,9.04735e-05    ,0.001030283    ,0.000400478    ,0.000199579    ,0.000610902    ,0.000870061    ,0.001148845    ,0.000706636    ,0.000491999    ,0.00093158    ,0.001224225    ,0.000652115    ,0.001117025    ,0.00481452    ,0.002359475    ,0.000641491    ,0.000126103    ,0.000158341    ,0.001455785    ,0.000227113    ,0.000366087    ,-0.000372672    ,-0.000814455    ,0.000403834    ,-0.00072791    ,-0.000422227    ,0.000413887    ,-0.000024651    ,0.000106683    ,0.00092596    ,0.000327427    ,0.000775269    ,0.000906505    ,0.000108337    ,-0.000214467    ,0.001249965    ,0.000694693    ,0.000718101    ,0.00083926    ,0.00138818    ,0.00131215    ,0.001113065    ,0.00134113    ,0.000972572    ,0.00073247    ,-3.35943e-05    ,0.000312545    ,0.000365923    ,0.000510906    ,0.000955084    ,0.000808904    ,0.000403238    ,-0.000150186    ,0.000200673    ,-0.00078668    ,-0.00078133    ,-0.00039955    ,-0.00226333    ,-0.00177287    ,0.000435995    ,-0.000397403};
static double met_Jun_2017[64] =    {0.00029923    ,5.42185e-05    ,-0.000023823    ,0.000187444    ,0.000650675    ,0.000135629    ,-0.000047364    ,0.000479489    ,0.000463565    ,0.000536066    ,0.000528847    ,0.000874895    ,0.000702853    ,0.000438988    ,0.000183642    ,0.000536044    ,0.004829985    ,0.00215068    ,5.94134e-05    ,-0.000164468    ,-4.89301e-05    ,0.000400028    ,-0.000277333    ,-0.000169422    ,-0.00136967    ,-0.000998661    ,0.000244959    ,-0.00123726    ,-0.000404182    ,-0.000511517    ,-0.000398515    ,-0.000314167    ,0.000546602    ,0.000547136    ,0.000388259    ,0.000108443    ,-0.000266431    ,-0.000734324    ,0.000946866    ,0.0001005    ,0.000739434    ,0.000530113    ,0.000821013    ,0.000748435    ,0.00094666    ,0.00102861    ,-0.00013052    ,0.000576223    ,-0.000657692    ,-0.000293389    ,0.000146337    ,-0.000338176    ,0.000241386    ,0.00005168    ,0.000495891    ,-1.76273e-05    ,8.84135e-05    ,-0.00153857    ,-0.000850848    ,-0.000986058    ,-0.002700175    ,-0.002620785    ,0.000440797    ,-0.00063843};
static double met_Jul_2017[64] =    {8.19405e-05    ,-0.000061716    ,0.000067827    ,0.000234624    ,0.000521161    ,0.000031451    ,0.000211139    ,0.000578157    ,0.000316435    ,0.000431433    ,0.000431228    ,0.000761259    ,0.000576742    ,0.000299923    ,0.00034315    ,0.000561323    ,0.004834095    ,0.00208259    ,-0.000155192    ,-0.000120373    ,-0.000328756    ,0.000265791    ,-0.000249398    ,-0.000310416    ,-0.001047503    ,-0.000884526    ,0.000431645    ,-0.000890023    ,-0.000330709    ,-0.000540866    ,-0.0002067    ,-0.000334498    ,0.00035945    ,0.000683828    ,0.000464533    ,0.000166613    ,-0.000276387    ,-0.00088093    ,0.001283446    ,0.00016521    ,0.000578248    ,0.000334565    ,0.000762716    ,0.000521964    ,0.000757276    ,0.000959739    ,0.000035001    ,0.000442429    ,-0.0006928    ,-0.000105863    ,-0.00006137    ,-0.000454825    ,-0.000157343    ,-5.68635e-05    ,0.000296553    ,-3.80711e-05    ,0.000384616    ,-0.00149983    ,-0.000650183    ,-0.00076473    ,-0.002494945    ,-0.00246317    ,0.000794264    ,-0.000337915};
static double met_Aug_2017[64] =    {0.000225917    ,0.00033629    ,1.70956e-05    ,0.000178987    ,0.00103744    ,0.000275415    ,0.000107582    ,0.000678724    ,0.000402821    ,0.000713535    ,0.000518528    ,0.000607008    ,0.000506748    ,0.000640514    ,0.000300349    ,0.000522693    ,0.005585285    ,0.00305222    ,0.000771553    ,0.001006485    ,0.000618719    ,0.00130769    ,0.000419028    ,0.000491215    ,-0.000815277    ,-0.000646205    ,0.000640702    ,-0.000881775    ,-0.000200994    ,-0.000258582    ,-0.000465663    ,-0.00004838    ,0.000591627    ,0.000615608    ,0.000197913    ,0.000352664    ,-0.000113565    ,-0.000637229    ,0.001241975    ,-6.659e-06    ,0.000649339    ,0.000641595    ,0.000673658    ,0.000657426    ,0.001263897    ,0.001258565    ,0.000050901    ,0.000318844    ,-9.23866e-05    ,0.00065074    ,0.000793493    ,0.000600095    ,0.001049984    ,0.000858337    ,0.000786737    ,0.000737051    ,0.000576774    ,-0.001384295    ,-0.000642989    ,-0.0008569    ,-0.00208747    ,-0.00247888    ,0.00108834    ,-0.000190064};

/*-----------------------------------------------------------------------------
                                 Private prototypes
 -----------------------------------------------------------------------------*/

cpl_error_code gravi_fit_profile (cpl_vector * values_x0,
                                  cpl_vector * values_y0,
                                  cpl_vector * values_sigma,
                                  cpl_image * mean_img,
                                  int ref_x, int ref_y,
                                  int size_profile,
                                  const char * resolution);

cpl_image * gravi_create_profile_image (cpl_image * mean_img,
                                        cpl_vector * values_x0,
                                        cpl_vector * values_y0,
                                        cpl_vector * values_sigma,
                                        cpl_size iy_min, cpl_size iy_max,
                                        const char * resolution);

/*-----------------------------------------------------------------------------
                             Functions code
 -----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
 * @brief Compute the DARK calibration map
 * 
 * @param raw_data  The input raw dark
 * 
 * @return The output DARK calibration map
 *
 * \exception CPL_ERROR_NULL_INPUT no raw_data as input
 *
 * The dark image of the SC is a saved as full image of the mean dark value
 * and mean dark standard deviation (RON). The dark for the FT is saved
 * into PIX array of the mean and dark standard deviation (RON). And the dark of
 * the metrology is save into the METOLOGY table.
 */
/*----------------------------------------------------------------------------*/

gravi_data * gravi_compute_dark (gravi_data * raw_data)
{
    gravi_msg_function_start(1);
	cpl_ensure (raw_data, CPL_ERROR_NULL_INPUT, NULL);
	
	/* Create the output DARK or SKY map */
	gravi_data * dark_map = gravi_data_new (0);

    /* Dump full header of RAW data */
	cpl_propertylist * dark_header = gravi_data_get_header (dark_map);
	cpl_propertylist * raw_header  = gravi_data_get_header (raw_data);
    cpl_propertylist_append (dark_header, raw_header);

	/* Check if this is a SKY or a DARK */
    const char * dpr_type = gravi_pfits_get_dpr_type (raw_header);
    int isSky = strstr (dpr_type, "SKY")?1:0;

	/* The dark file must contains all the shutter close */
	if ( isSky==0 && !gravi_data_check_shutter_closed (raw_data) ) {
        gravi_pfits_add_check (dark_header, "DARK has some shutter OPEN !!");
	}

	/* The sky file must contains all the shutter open */
	if ( isSky==1 && !gravi_data_check_shutter_open (raw_data) ) {
        gravi_pfits_add_check (dark_header, "SKY has some shutter CLOSED !!");
	}
	
    /*
     * Compute the SC DARK
     */
    if (!gravi_data_has_extension (raw_data, GRAVI_IMAGING_DATA_SC_EXT)) {
        cpl_msg_warning (cpl_func,"The DARK data has no IMAGING_DATA_SC");
    }
    else
    {
	    cpl_msg_info (cpl_func, "Computing the %s of SC",isSky?"SKY":"DARK");
        
        /* Copy IMAGING_DETECTOR into product */
        gravi_data_copy_ext (dark_map, raw_data, GRAVI_IMAGING_DETECTOR_SC_EXT);

		/* Load the IMAGING_DATA table or image list */
		cpl_imagelist * imglist = gravi_data_get_cube (raw_data, GRAVI_IMAGING_DATA_SC_EXT);

		/* Compute the median image of the imagelist */
		cpl_image * median_img = cpl_imagelist_collapse_sigclip_create (imglist, 5, 5, 0.6, CPL_COLLAPSE_MEDIAN_MEAN, NULL);
		CPLCHECK_NUL ("Cannot compute the median dark");

		/* Compute STD. We should see if we use sigma-clipping
		 * or not for the collapse, it may change the bad pixel detection */
		cpl_msg_info (cpl_func,"Compute std with imglist");
		cpl_imagelist * temp_imglist = cpl_imagelist_duplicate (imglist);
		cpl_imagelist_subtract_image (temp_imglist, median_img);
		cpl_imagelist_power (temp_imglist, 2.0);
		// cpl_image * stdev_img = cpl_imagelist_collapse_create (temp_imglist);
		cpl_image * stdev_img = cpl_imagelist_collapse_sigclip_create (temp_imglist, 5, 5, 0.6,
									       CPL_COLLAPSE_MEDIAN_MEAN, NULL);
		FREE (cpl_imagelist_delete, temp_imglist);
		cpl_image_power (stdev_img, 0.5);
		CPLCHECK_NUL ("Cannot compute the STD of the DARK");
		
		/* Compute the QC parameters RMS and MEDIAN */
		cpl_msg_info (cpl_func, "Compute QC parameters");
		double mean_qc = cpl_image_get_median (median_img);
		double darkrms = cpl_image_get_median (stdev_img);
		cpl_propertylist_update_double (dark_header, isSky?QC_MEANSKY_SC:QC_MEANDARK_SC, mean_qc);
		cpl_propertylist_update_double (dark_header, isSky?QC_SKYRMS_SC:QC_DARKRMS_SC, darkrms);

		/* Verbose */
	    cpl_msg_info (cpl_func, "QC_MEDIAN%s_SC = %e",isSky?"SKY":"DARK", mean_qc);
	    cpl_msg_info (cpl_func, "QC_%sRMS_SC = %e",isSky?"SKY":"DARK", darkrms);

		/* Put the data in the output table : dark_map */
		cpl_propertylist * img_plist = gravi_data_get_plist (raw_data, GRAVI_IMAGING_DATA_SC_EXT);
		img_plist = cpl_propertylist_duplicate (img_plist);
		gravi_data_add_img (dark_map, img_plist, GRAVI_IMAGING_DATA_SC_EXT, median_img);
        
		img_plist = cpl_propertylist_duplicate (img_plist);
		gravi_data_add_img (dark_map, img_plist, GRAVI_IMAGING_ERR_SC_EXT, stdev_img);
		CPLCHECK_NUL ("Cannot set the SC data");
        
    } /* End SC case */

    /*
     * Compute the FT DARK
     */
    if (!gravi_data_has_extension (raw_data, GRAVI_IMAGING_DATA_FT_EXT)) {
        cpl_msg_warning (cpl_func,"The DARK data has no IMAGING_DATA_FT");
    }
    else
    {            
	    cpl_msg_info(cpl_func, "Computing the %s of FT",isSky?"SKY":"DARK");

        /* Copy IMAGING_DETECTOR into product */
        gravi_data_copy_ext (dark_map, raw_data, GRAVI_IMAGING_DETECTOR_FT_EXT);
        
		/* Load the IMAGING_DATA table as DOUBLE */
        cpl_msg_info (cpl_func,"Load data");
		cpl_table * table_ft = gravi_data_get_table (raw_data, GRAVI_IMAGING_DATA_FT_EXT);
        cpl_table_cast_column (table_ft, "PIX", "PIX", CPL_TYPE_DOUBLE);
		cpl_imagelist * imglist = gravi_imagelist_wrap_column (table_ft, "PIX");
		CPLCHECK_NUL ("Cannot load the FT data");

		/* Compute the median image of the imagelist */
        cpl_msg_info (cpl_func,"Compute mean and median");
		cpl_image * median_img = cpl_imagelist_collapse_median_create (imglist);
		cpl_image * mean_img   = cpl_imagelist_collapse_create (imglist);
		CPLCHECK_NUL ("Cannot compute the MEAN dark");
        
		/* Compute the std of each pixels */
        cpl_msg_info (cpl_func,"Compute std");
        cpl_imagelist_subtract_image (imglist, mean_img);
        cpl_imagelist_power (imglist, 2.0);
		cpl_image * stdev_img = cpl_imagelist_collapse_create (imglist);
		cpl_image_power (stdev_img, 0.5);

        /* We power back the data, because we are working in-place */
        cpl_imagelist_power (imglist, 0.5);
        cpl_imagelist_add_image (imglist, mean_img);

        /* Delete data */
        cpl_msg_info (cpl_func,"Delete data");
        FREE (gravi_imagelist_unwrap_images, imglist);
		FREE (cpl_image_delete, mean_img);
		CPLCHECK_NUL ("Cannot compute the STD of the DARK");

		/* Compute the QC parameters RMS and MEDIAN */
		cpl_msg_info (cpl_func, "Compute QC parameters");
		double mean_qc = cpl_image_get_mean (median_img);
		double darkrms = cpl_image_get_mean (stdev_img);
		cpl_propertylist_update_double (dark_header, isSky?QC_MEANSKY_FT:QC_MEANDARK_FT, mean_qc);
		cpl_propertylist_update_double (dark_header, isSky?QC_SKYRMS_FT:QC_DARKRMS_FT, darkrms);
		CPLCHECK_NUL ("Cannot compute the QC");

		/* Verbose */
	    cpl_msg_info (cpl_func, "QC_MEDIAN%s_FT = %e",isSky?"SKY":"DARK", mean_qc);
	    cpl_msg_info (cpl_func, "QC_%sRMS_FT = %e",isSky?"SKY":"DARK", darkrms);

        /* Create the output DARK table, with a single row */
		cpl_table * median_table = cpl_table_extract (table_ft, 0, 1);
        cpl_array * median_array = gravi_array_wrap_image (median_img);
        cpl_table_set_array (median_table, "PIX", 0, median_array);
        FREE (cpl_array_unwrap, median_array);
        FREE (cpl_image_delete, median_img);
		CPLCHECK_NUL("Cannot set median in table");

		/* Put median dark in the output gravi_data */
		gravi_data_add_table (dark_map, NULL, GRAVI_IMAGING_DATA_FT_EXT, median_table);
		CPLCHECK_NUL("Cannot save median in gravi_data");
        
        /* Create the output DARK RMS table, with a single row */
		cpl_table * stdev_table = cpl_table_extract (table_ft, 0, 1);
        cpl_array * stdev_array = gravi_array_wrap_image (stdev_img);
        cpl_table_set_array (stdev_table, "PIX", 0, stdev_array);
        FREE (cpl_array_unwrap, stdev_array);
        FREE (cpl_image_delete, stdev_img);
		CPLCHECK_NUL("Cannot set rms in table");

		/* Put median dark in the output gravi_data */
		gravi_data_add_table (dark_map, NULL, GRAVI_IMAGING_ERR_FT_EXT, stdev_table);
		CPLCHECK_NUL("Cannot save median in gravi_data");
        
    } /* End case FT */
    
    /*
     * Compute the METROLOGY DARK
     */
    if (( isSky==1 )||(!gravi_data_has_extension (raw_data, GRAVI_METROLOGY_EXT))) {
        
        if ( isSky==0 )
            cpl_msg_warning (cpl_func,"The DARK data has no METROLOGY");

    }
    else
    {
        
        cpl_msg_info(cpl_func, "Computing the %s of METROLOGY",isSky?"SKY":"DARK");
        
        /* Load the IMAGING_DATA table as DOUBLE */
        cpl_msg_info (cpl_func,"Load data");
        cpl_table * table_met = gravi_data_get_table (raw_data, GRAVI_METROLOGY_EXT);
        cpl_table_cast_column (table_met, "VOLT", "VOLT", CPL_TYPE_DOUBLE);
        cpl_imagelist * imglist = gravi_imagelist_wrap_column (table_met, "VOLT");
        CPLCHECK_NUL ("Cannot load the VOLT data");
        
        /* Compute the median image of the imagelist */
        cpl_msg_info (cpl_func,"Compute mean and median");
        cpl_image * median_img = cpl_imagelist_collapse_median_create (imglist);
        cpl_image * mean_img   = cpl_imagelist_collapse_create (imglist);
        CPLCHECK_NUL ("Cannot compute the MEAN dark");
        
        /* Compute the std of each pixels */
        cpl_msg_info (cpl_func,"Compute std");
        cpl_imagelist_subtract_image (imglist, mean_img);
        cpl_imagelist_power (imglist, 2.0);
        cpl_image * stdev_img = cpl_imagelist_collapse_create (imglist);
        cpl_image_power (stdev_img, 0.5);
        
        /* We power back the data, because we are working in-place */
        cpl_imagelist_power (imglist, 0.5);
        cpl_imagelist_add_image (imglist, mean_img);
        
        /* Delete data */
        cpl_msg_info (cpl_func,"Delete data");
        FREE (gravi_imagelist_unwrap_images, imglist);
        FREE (cpl_image_delete, mean_img);
        CPLCHECK_NUL ("Cannot compute the STD of the DARK");
        
        /* Compute the QC parameters RMS and MEDIAN */
        cpl_msg_info (cpl_func, "Compute QC parameters");
        double mean_qc = cpl_image_get_mean_window (median_img,1,1,64,1);
        double darkrms = cpl_image_get_mean_window (stdev_img,1,1,64,1);
        double darkmin = cpl_image_get_min_window (median_img,1,1,64,1);
        double darkmax = cpl_image_get_max_window (median_img,1,1,64,1);
        if (darkmin>0) darkmin=0;
        if (darkmax<0) darkmax=0;
        cpl_propertylist_update_double (dark_header, QC_MEANDARK_MET, mean_qc);
        cpl_propertylist_update_double (dark_header, QC_DARKRMS_MET, darkrms);
        cpl_propertylist_update_double (dark_header, QC_DARKRANGE_MET, darkmax-darkmin);
        CPLCHECK_NUL ("Cannot compute the QC");
        
        /* Verbose */
        cpl_msg_info (cpl_func, "QC_MEDIAN%s_MET = %e","DARK", mean_qc);
        cpl_msg_info (cpl_func, "QC_%sRMS_MET = %e","DARK", darkrms);
        cpl_msg_info (cpl_func, "QC_%sRANGE_MET = %e","DARK", darkmax-darkmin);
        
        /* Create the output DARK table, with a single row */
        cpl_table * median_table = cpl_table_extract (table_met, 0, 1);
        cpl_array * median_array = gravi_array_wrap_image (median_img);
        
        /* Put to zero the dark of the fiber coupler unit */
        for (cpl_size diode = 64; diode < 80; diode++){
            cpl_array_set(median_array, diode, 0);
        }
        
        /* Check on time */
        double time_mjd_obs= cpl_propertylist_get_double (raw_header, "MJD-OBS");
        
        /* If the data is too old, we are using an old dataset */
        for (cpl_size diode = 0; diode < 64; diode++)
        {
            if (time_mjd_obs < 57747) cpl_array_set(median_array, diode, met_Sep_2016[diode]); /* 25 decembre 2016 */
            else if (time_mjd_obs < 57851) cpl_array_set(median_array, diode, met_Mar_2017[diode]); /* 8 April */
            else if (time_mjd_obs < 57924) cpl_array_set(median_array, diode, met_Jun_2017[diode]); /* 20 Juin */
            else if (time_mjd_obs < 57990) cpl_array_set(median_array, diode, met_Jul_2017[diode]); /* 25 Juillet */
            else if (time_mjd_obs < 58208.01) cpl_array_set(median_array, diode, met_Aug_2017[diode]);
        }
        
        cpl_table_set_array (median_table, "VOLT", 0, median_array);
        FREE (cpl_array_unwrap, median_array);
        FREE (cpl_image_delete, median_img);
        CPLCHECK_NUL("Cannot set median in table");
        
        /* Put median dark in the output gravi_data */
        gravi_data_add_table (dark_map, NULL, GRAVI_METROLOGY_EXT, median_table);
        CPLCHECK_NUL("Cannot save median in gravi_data");
        
        /* Create the output DARK RMS table, with a single row */
        cpl_table * stdev_table = cpl_table_extract (table_met, 0, 1);
        cpl_array * stdev_array = gravi_array_wrap_image (stdev_img);
        cpl_table_set_array (stdev_table, "VOLT", 0, stdev_array);
        FREE (cpl_array_unwrap, stdev_array);
        FREE (cpl_image_delete, stdev_img);
        CPLCHECK_NUL("Cannot set rms in table");
        
        /* Put median dark in the output gravi_data */
        gravi_data_add_table (dark_map, NULL, GRAVI_METROLOGY_ERR_EXT, stdev_table);
        CPLCHECK_NUL("Cannot save median in gravi_data");
        
    }  /* End case METROLOGY */

    /*
     * Compute the ACQ DARK
     */
    if (!gravi_data_has_extension (raw_data, GRAVI_IMAGING_DATA_ACQ_EXT)) {
        cpl_msg_warning (cpl_func,"The DARK data has no IMAGING_DATA_ACQ");
    }
    else
    {
	    cpl_msg_info (cpl_func, "Computing the %s of ACQ",isSky?"SKY":"DARK");
        
		/* Load the IMAGING_DATA table or image list */
		cpl_imagelist * imglist = gravi_data_get_cube (raw_data, GRAVI_IMAGING_DATA_ACQ_EXT);

		/* Compute the median image of the imagelist */
		cpl_image * median_img = cpl_imagelist_collapse_median_create (imglist);
		CPLCHECK_NUL ("Cannot compute the median dark");

		/* Compute the QC parameters RMS and MEDIAN */
		cpl_msg_info (cpl_func, "Compute QC parameters");
		double mean_qc = cpl_image_get_median (median_img);
        
		/* Verbose */
	    cpl_msg_info (cpl_func, "QC_MEDIAN%s_ACQ = %e",isSky?"SKY":"DARK", mean_qc);
		cpl_propertylist_update_double (dark_header, isSky?"ESO QC MEDIANSKY ACQ":"ESO QC MEDIANDARK ACQ",
                                        mean_qc);

		/* Put the data in the output table : dark_map */
		cpl_propertylist * img_plist = gravi_data_get_plist (raw_data, GRAVI_IMAGING_DATA_ACQ_EXT);
		img_plist = cpl_propertylist_duplicate (img_plist);
		gravi_data_add_img (dark_map, img_plist, GRAVI_IMAGING_DATA_ACQ_EXT, median_img);
		CPLCHECK_NUL("Cannot save median in gravi_data");
    } 

    gravi_msg_function_exit(1);
    return dark_map;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Average several DARK calibration map
 *
 * @param data     List of input DARK map
 * @param ndata    Number of DARK in the list
 *
 * @return A newly allocated DARK calibration map
 *
 * \exception CPL_ERROR_NULL_INPUT no data as input
 *
 */
/*---------------------------------------------------------------------------*/

gravi_data * gravi_average_dark (gravi_data ** data, cpl_size ndata)
{
  gravi_msg_function_start(0);
  cpl_ensure (data, CPL_ERROR_NULL_INPUT, NULL);
  
  gravi_msg_warning ("FIXME","Averaging DARK or SKY is experimental");

  gravi_data * output_data = gravi_data_duplicate (data[0]);

  /* Average the IMAGING_DATA of SC */
  cpl_msg_info (cpl_func,"Average IMAGING_DATA of SC");
  
  cpl_image * darksc_image = gravi_data_get_img (output_data, GRAVI_IMAGING_DATA_SC_EXT);
  for (int file = 1; file < ndata; file++) {
      cpl_image_add (darksc_image, gravi_data_get_img (data[file],  GRAVI_IMAGING_DATA_SC_EXT));
  }
  cpl_image_divide_scalar (darksc_image, ndata);
  CPLCHECK_NUL ("Cannot average DARK of SC");

  /* Average the IMAGING_ERR of SC */
  cpl_msg_info (cpl_func,"Average IMAGING_ERR of SC");
  
  cpl_image * stdevsc_image = gravi_data_get_img (output_data, GRAVI_IMAGING_ERR_SC_EXT);
  for (int file = 1; file < ndata; file++) {
      cpl_image_add (stdevsc_image, gravi_data_get_img (data[file],  GRAVI_IMAGING_ERR_SC_EXT));
  }
  cpl_image_divide_scalar (stdevsc_image, ndata);
  CPLCHECK_NUL ("Cannot average DARKERR of SC");

  /* Average the IMAGING_DATA of FT */
  cpl_msg_info (cpl_func,"Average IMAGING_DATA of FT");
  
  cpl_array * darkft_array;
  darkft_array = cpl_table_get_data_array (gravi_data_get_table (output_data, GRAVI_IMAGING_DATA_FT_EXT), "PIX")[0];
  for (int file = 1; file < ndata; file++) {
      cpl_array_add (darkft_array, cpl_table_get_data_array (gravi_data_get_table (data[file], GRAVI_IMAGING_DATA_FT_EXT), "PIX")[0]);
  }
  cpl_array_divide_scalar (darkft_array, ndata);
  CPLCHECK_NUL ("Cannot average DARK of FT");
  
  /* Average the IMAGING_ERR of FT */
  cpl_msg_info (cpl_func,"Average IMAGING_ERR of FT");
  
  cpl_array * errft_array;
  errft_array = cpl_table_get_data_array (gravi_data_get_table (output_data, GRAVI_IMAGING_ERR_FT_EXT), "PIX")[0];
  for (int file = 1; file < ndata; file++) {
      cpl_array_add (errft_array, cpl_table_get_data_array (gravi_data_get_table (data[file], GRAVI_IMAGING_ERR_FT_EXT), "PIX")[0]);
  }
  cpl_array_divide_scalar (errft_array, ndata);
  CPLCHECK_NUL ("Cannot average DARK of FT");
  
  gravi_msg_function_exit(0);
  return output_data;
}


/*----------------------------------------------------------------------------*/
/**
 * @brief Fit the profile parameters in an image
 * 
 * @param values_x0     Output vector of spectral channels (0..nx-1)
 * @param values_y0     Output vector of profile center in spatial direction
 * @param values_sigma  Output vector of profile width in spatial direction
 * @param mean_img      Input image (nx,ny)
 * @param ref_x, ref_y  Coordinate of a known pixel inside the spectra
 * @param size_profile  Spatial extend of the data to extract and fit
 * @param resolution    Spectral resolution (LOW, MED, HIGH)
 *
 * \exception CPL_ERROR_NULL_INPUT input data is missing
 * \exception CPL_ERROR_ILLEGAL_INPUT ref_x or ref_y are outside boundaries
 *
 * values_x0, values_y0, values_sigma are of size nx and unit [pixel] in
 * the FITS convention (1..nx). The function starts to fit the profile
 * toward the right from ref_x, then toward the left from ref_x.
 * In resolution HIGH, the ref_y is updated at each spectral channel to
 * follow the curvature when extracting the data.
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_fit_profile (cpl_vector * values_x0,
                                  cpl_vector * values_y0,
                                  cpl_vector * values_sigma,
                                  cpl_image * mean_img,
                                  int ref_x, int ref_y,
                                  int size_profile,
                                  const char * resolution)
{
    gravi_msg_function_start(0);
    cpl_ensure_code (values_x0,    CPL_ERROR_NULL_INPUT);
    cpl_ensure_code (values_y0,    CPL_ERROR_NULL_INPUT);
    cpl_ensure_code (values_sigma, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code (mean_img,     CPL_ERROR_NULL_INPUT);
    cpl_ensure_code (resolution,   CPL_ERROR_NULL_INPUT);

    /* Compute middle of size_profile */
    int middle = floor (size_profile / 2);

    /* Get the image dimensions */
    int nx = cpl_image_get_size_x (mean_img);
    int ny = cpl_image_get_size_y (mean_img);

    cpl_ensure_code (ref_x > 0 && ref_x < nx, CPL_ERROR_ILLEGAL_INPUT);
    cpl_ensure_code (ref_y > 0 && ref_y < ny, CPL_ERROR_ILLEGAL_INPUT);

    /* y0: best fit profile position (y0_ is of previous channel)
     * sigma: best fit profile width (sigma_ is of previous channel) */
    double y0, y0_ = ref_y;
    double sigma, sigma_ = 0;
    
    /* Loop on spectral channels. Here fit all the
     * points under the reference point */
    int coord_x = ref_x;
    while (coord_x <= nx) {
        
        /* Extract the corresponding column in the 2D image */
        cpl_vector * line_data;
        line_data = cpl_vector_new_from_image_column (mean_img,
                                                      coord_x);
        
        /* Construction of the vector of references points y_data and
         * its values z_data for the gaussian fit */
        cpl_vector * y_data, * z_data;
        y_data = cpl_vector_new (size_profile);
        
        if (! (strcmp(resolution, "LOW") && strcmp(resolution, "MED")) )	{
            /* Case LOW and MED */
            for (cpl_size i = 0; i < size_profile; i++){
                cpl_vector_set (y_data, i, ref_y - size_profile/2 + i);
            }
            z_data = cpl_vector_extract (line_data,
                                         ref_y - size_profile/2,
                                         ref_y - size_profile/2+size_profile-1, 1);
        }
        else {
            /* Case HIGH */
            if ((ref_y - middle) <= 0 ) {
                for (cpl_size i = 0; i < size_profile; i++){
                    cpl_vector_set (y_data, i, i);
                }
                z_data = cpl_vector_extract (line_data, 0,
                                             size_profile - 1, 1);
            }
            else if ((ref_y - middle +  size_profile) >= ny ){
                for (cpl_size i = 0; i < size_profile; i++){
                    cpl_vector_set (y_data, i, ny - size_profile + i);
                }
                z_data = cpl_vector_extract (line_data,
                                             ny - size_profile,
                                             ny - 1, 1);
            }
            else {
                for (cpl_size i = 0; i < size_profile; i++){
                    cpl_vector_set(y_data, i, ref_y + (i - middle));
                }
                z_data = cpl_vector_extract (line_data, ref_y - middle,
                                             ref_y - middle + size_profile - 1, 1);
            }
        }
        
        /* Fit a Gaussian to the extracted data */
        cpl_errorstate prestate = cpl_errorstate_get();
        double area, mse, offset = 0;
        cpl_vector_fit_gaussian (y_data, NULL, z_data, NULL,
                                 CPL_FIT_ALL, &y0, &sigma, &area,
                                 &offset, &mse, NULL, NULL);
        
        /* If the fit fail, we keep the value of the
         * previous spectral channel */
        if (cpl_error_get_code() == CPL_ERROR_CONTINUE){
            cpl_errorstate_set (prestate);
            if (coord_x == ref_x) {
                y0 =  ref_y;
                sigma = 1; 
            }else{
                y0 =  y0_; 
                sigma = sigma_;
            }
            cpl_msg_warning (cpl_func, "Cannot fit profile of channel %d, new y0=%e", coord_x, y0);
        }
        
        CPLCHECK_MSG ("Error during the gaussian fit");
        
        /* If more than 1 pixel shift, we keep the value of the
         * previous spectral channel */
        if ((fabs(y0 - y0_) >= 1) && (coord_x != ref_x)) {
            cpl_msg_warning (cpl_func, "Too much difference of channel %d with previous (%e pixels)", coord_x, y0 - y0_);
            y0 = y0_;
            sigma = sigma_;
        }
        
        /* In HIGH Resolution, we update ref_y in order to
         * follow the curvature when extracting the data */
        if (! (strcmp(resolution, "HIGH"))) {
            ref_y = floor(y0);
        }
        
        /* Record x0 and sigma to compare to next channel */
        y0_ = y0;
        sigma_ = sigma;
        
        /* Save best fit y0 and sigma into output vectors */
        cpl_vector_set (values_x0, (coord_x - 1), coord_x - 1);
        cpl_vector_set (values_y0, (coord_x - 1), y0);
        cpl_vector_set (values_sigma, (coord_x - 1), sigma);
        
        /* Increment spectral channel number */
        coord_x ++;
        
        FREE (cpl_vector_delete, line_data);
        FREE (cpl_vector_delete, z_data);
        FREE (cpl_vector_delete, y_data);
        
    } /* End fit the column ref_x -> nx */
    
    
    /* Reset coord_x to the middle of the spectra */
    coord_x = ref_x - 1;
    
    /* In HIGH, we use the best-fit position for this channel */
    if (! (strcmp(resolution, "HIGH"))) {
        ref_y = cpl_vector_get (values_y0, (ref_x - 1));
    }
    
    
    /* Loop on spectral channels. Here fit all the
     * points under the reference point */
    while (coord_x > 0) {
        
        /* Extract the corresponding column in the 2D image */
        cpl_vector * line_data;
        line_data = cpl_vector_new_from_image_column (mean_img,
                                                      coord_x);
        
        /* Construction of the vector of references points y_data and
         * its values z_data for the gaussian fit */
        cpl_vector * y_data, * z_data;
        y_data = cpl_vector_new(size_profile);
        
        if (! (strcmp(resolution, "LOW") && strcmp(resolution, "MED")) )	{
            /* Case LOW and MED */
            for (cpl_size i = 0; i < size_profile; i++){
                cpl_vector_set(y_data, i, ref_y - size_profile/2 + i);
            }
            z_data = cpl_vector_extract (line_data,
                                         ref_y - size_profile/2,
                                         ref_y - size_profile/2 +size_profile-1, 1);
        }
        else {
            /* Case HIGH */
            if ((ref_y - middle) <= 0 ){
                for (cpl_size i = 0; i < size_profile; i++){
                    cpl_vector_set(y_data, i, i);
                }
                z_data = cpl_vector_extract (line_data, 0,
                                             size_profile - 1, 1);
            }
            else if ((ref_y - middle + size_profile) >= ny ){
                for (cpl_size i = 0; i < size_profile; i++){
                    cpl_vector_set(y_data, i, ny - size_profile + i);
                }
                z_data = cpl_vector_extract (line_data,
                                             ny - size_profile,
                                             ny - 1, 1);
            }
            else{
                for (cpl_size i = 0; i < size_profile; i++){
                    cpl_vector_set(y_data, i, ref_y + (i - middle));
                }
                z_data = cpl_vector_extract (line_data, ref_y - middle,
                                             ref_y - middle + size_profile - 1, 1);
            }
        }
        
        /* Fit a Gaussian to the extracted data */
        cpl_errorstate prestate = cpl_errorstate_get();
        double area, mse, offset = 0;
        cpl_vector_fit_gaussian (y_data, NULL, z_data, NULL,
                                 CPL_FIT_ALL, &y0, &sigma, &area,
                                 &offset, &mse, NULL, NULL);
        
        /* If the fit fail, we keep the value of the
         * previous spectral channel */
        if (cpl_error_get_code() == CPL_ERROR_CONTINUE){
            cpl_errorstate_set (prestate);
            if (coord_x == ref_x-1){
                y0 =  ref_y;
                sigma = 1;
            } else {
                y0 = y0_; 
                sigma = sigma_;
            }
            cpl_msg_warning (cpl_func, "Cannot fit profile of channel %d, new y0=%e", coord_x, y0);
        }
        
        CPLCHECK_MSG ("Error during the gaussian fit");
        
        /* If more than 1 pixel shift, we keep the value of the
         * previous spectral channel */
        if ((fabs(y0 - y0_) >= 1) && (coord_x != ref_x-1)) {
            cpl_msg_warning (cpl_func, "Too much difference of channel %d with previous (%e pixels)", coord_x, y0 - y0_);
            y0 = y0_;
            sigma = sigma_;
        }
        
        /* In HIGH Resolution, we update ref_y in order to
         * follow the curvature when extracting the data */
        if (! (strcmp(resolution, "HIGH"))) {
            ref_y = floor(y0);
        }
        
        /* Record x0 and sigma to compare to next channel */
        y0_ = y0;
        sigma_ = sigma;
        
        /* Save best fit y0 and sigma into output vectors */
        cpl_vector_set (values_x0, (coord_x - 1), coord_x);
        cpl_vector_set (values_y0, (coord_x - 1), y0 );
        cpl_vector_set (values_sigma, (coord_x - 1), sigma);
        
        /* Increment spectral channel number */
        coord_x --;
        
        FREE (cpl_vector_delete, line_data);
        FREE (cpl_vector_delete, z_data);
        FREE (cpl_vector_delete, y_data);
        
    } /* End fit the column ref_x -> 0 */
    

    gravi_msg_function_exit(0);
    return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Create the profile image from image and/or fitted params
 * 
 * @param mean_img      Input image (nx,ny)
 * @param values_x0     Input vector of spectral channels (0..nx-1)
 * @param values_y0     Input vector of profile center in spatial direction
 * @param values_sigma  Input vector of profile width in spatial direction
 * @param iy_min        Min coordinate where to fill the profile (0..ny-1)
 * @param iy_max        Max coordinate where to fill the profile (0..ny-1)
 * @param mode          "BOX", "PROFILE", "GAUSS"
 * 
 * @return The image of the profile
 *
 * \exception CPL_ERROR_NULL_INPUT input data is missing
 * \exception CPL_ERROR_ILLEGAL_INPUT iy_min or iy_max are outside boundaries
 *
 * The vectors values_x0, values_y0, values_sigma are first smoothed
 * (duplicated) with a polynomial interpolation before being used.
 * The image of the profile is built either from mean_img (LOW, MED) or
 * is a boxcard of 5 pixels (HIGH).
 */
/*----------------------------------------------------------------------------*/

cpl_image * gravi_create_profile_image (cpl_image * mean_img,
                                        cpl_vector * values_x0,
                                        cpl_vector * values_y0,
                                        cpl_vector * values_sigma,
                                        cpl_size iy_min,
                                        cpl_size iy_max,
                                        const char * mode)
{
    int nv = 0;
	gravi_msg_function_start(0);
    cpl_ensure (values_x0,    CPL_ERROR_NULL_INPUT, NULL);
    cpl_ensure (values_y0,    CPL_ERROR_NULL_INPUT, NULL);
    cpl_ensure (values_sigma, CPL_ERROR_NULL_INPUT, NULL);
    cpl_ensure (mean_img,     CPL_ERROR_NULL_INPUT, NULL);
    cpl_ensure (mode,         CPL_ERROR_NULL_INPUT, NULL);

    /* Get the image dimensions */
    cpl_size nx = cpl_image_get_size_x (mean_img);
    cpl_size ny = cpl_image_get_size_y (mean_img);

    cpl_msg_info (cpl_func, "iy_min = %lld, iy_max = %lld, ny = %lld",
                  iy_min, iy_max, ny);
    
    cpl_ensure (iy_min >= 0 && iy_min < ny, CPL_ERROR_ILLEGAL_INPUT, NULL);
    cpl_ensure (iy_max >= 0 && iy_max < ny, CPL_ERROR_ILLEGAL_INPUT, NULL);

    /* Filter the profile params with a polynomial fit in spectral direction
     * FIXME: could be a running median instead of fit, surely more stable.
     * However the values contains some 'accidents' or modulations that
     * may not be real and are efficiently removed by polynomial fit. */
    cpl_matrix * matrix = cpl_matrix_wrap (1, nx, cpl_vector_get_data(values_x0));
    cpl_size power;
    
    /* Fit the y0    --  apply a median filter first */
    cpl_polynomial * y0_poly = cpl_polynomial_new(1);
    power = 4;
    cpl_vector * temp_median = cpl_vector_filter_median_create(values_y0,3);
    cpl_polynomial_fit (y0_poly, matrix, NULL, temp_median, NULL,
                        CPL_FALSE, NULL, &power);
    cpl_vector_delete (temp_median);
    CPLCHECK_NUL ("Cannot fit the y0");
    
    /* Fit the sigma */
    cpl_polynomial * sigma_poly = cpl_polynomial_new(1);
    power = 5;
    cpl_polynomial_fit (sigma_poly, matrix, NULL, values_sigma,
                        NULL, CPL_FALSE, NULL, &power);
    CPLCHECK_NUL ("Cannot fit the sigma");
    
    /* Compute the new y0 and sigma from these fits */
    cpl_vector * valuesfit_y0 = cpl_vector_new(nx);
    cpl_vector * valuesfit_sig = cpl_vector_new(nx);
    
    for (cpl_size ix = 0; ix < nx; ix++){
        double result;
        result = cpl_polynomial_eval_1d (y0_poly, cpl_vector_get(values_x0, ix), NULL);
        cpl_vector_set(valuesfit_y0, ix, result);
        result = cpl_polynomial_eval_1d (sigma_poly, cpl_vector_get(values_x0, ix), NULL);
        cpl_vector_set(valuesfit_sig, ix, result);
    }
    FREE (cpl_polynomial_delete, y0_poly);
    FREE (cpl_polynomial_delete, sigma_poly);
    cpl_matrix_unwrap (matrix);
    
    /* 
     * Allocate image profile
     */
    
    cpl_image * region_img;
    region_img = cpl_image_new (nx, ny, CPL_TYPE_DOUBLE);
    cpl_image_fill_window (region_img, 1,1,nx,ny, 0.0);
	
    /* Loop on spectral direction */
    for (cpl_size ix = 0; ix < nx; ix++){
        
        double sum_flux = 0, sum_flux2 = 0;
        
        /* Loop on spatial direction */
        for (cpl_size iy = iy_min; iy <= iy_max; iy++ ){
            
            double result;
            
            /* We use the measured profile */
            if (!strcmp (mode, "PROFILE")) {
                result = cpl_image_get (mean_img, ix+1, iy+1, &nv);
            }
            /* We use a fixed 5 pixel boxcar profile */
            else if (!strcmp (mode, "BOX")) {
                result = ( (fabs(iy - cpl_vector_get (valuesfit_y0, ix)) < 3 ) ? 1.0 : 0.0);
            }
            /* We use a Gaussian profile */
            else if (!strcmp (mode, "GAUSS")) {
                result = (iy - cpl_vector_get (valuesfit_y0, ix)) /
                          cpl_vector_get (valuesfit_sig, ix);
                result = exp( - pow (result, 2) / 2);
            } else {
                cpl_msg_error (cpl_func, "BUG, report to DRS team");
                return NULL;
            }
            
            /* Fill the profile image of this region */
            cpl_image_set (region_img, ix + 1, iy + 1, result);
            
            /* Compute normalization coefficients */
            sum_flux  += result;
            sum_flux2 += result * result;
        } /* End loop on spatial direction */
        
        /* Keep only effective part of the profile 
         * Force pixel <1e-7 to zero */
        for (cpl_size iy = iy_min; iy <= iy_max; iy++ ) {
            double img_j = cpl_image_get (region_img, ix+1, iy+1, &nv);
            if (img_j / sum_flux < 1e-7) 
                cpl_image_set (region_img, ix+1, iy+1, 0.0);
        }
        
        /* When we use a complex profile, we have to normalize
         * it to make it flux conservative at extraction */
        if (!strcmp (mode, "PROFILE") || !strcmp (mode, "GAUSS") ) {
            
            double sum_flux = 0.0;
            double sum_flux2 = 0.0;
            for (cpl_size iy = iy_min; iy <= iy_max; iy++ ) {
                double img_j = cpl_image_get (region_img, ix+1, iy+1, &nv);
                sum_flux  += img_j;
                sum_flux2 += img_j * img_j;
            }
            
            for (cpl_size iy = iy_min; iy <= iy_max; iy++ ) {
                double img_j = cpl_image_get (region_img, ix+1, iy+1, &nv);
                cpl_image_set (region_img, ix+1, iy+1, img_j * sum_flux / sum_flux2);
            }
        }
		
    } /* End loop on columns = spectral direction */

    /* Deallocated of the fitted variables */
    FREE (cpl_vector_delete, valuesfit_y0);
    FREE (cpl_vector_delete, valuesfit_sig);

	gravi_msg_function_exit(0);
    return region_img;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Computes the spatial profile of each spectrum for
 * 	  	  optimal extraction purpose.
 * 
 * @param flats_data	The FLAT_RAW datas. Each of the flats has
 *                      the shutter of only one telescope open
 * @param dark_map	 	The DARK calibration map
 * @param bad_map	    BAD calibration map
 * @param nflat  	    The number of FLAT file inputs
 * @param params	  	Input parameter list with :
 *                      - profile-width : Width of the detector window extracted
 *                       around the default position of each spectrum, and on
 *                       which the profile will be applied to perform the extraction.
 *                      - force-badpix-to-zero : Force the badpixel to zero in profile.
 *                      - profile-mode : Method to compute the extraction profile.
 *                       PROFILE corresponds to the pixel intensities measured in the
 *                       FLAT files (Gaussian like with FWHM of approx 1.5 pixel).
 *                       This is the AUTO option for the Low and Med spectral
 *                       resolution. GAUSS corresponds to a Gaussian fit of the
 *                       (non-zero) pixel intensities measured in the FLAT files.
 *                        BOX corresponds to a box-card of 6 pixels centered on the
 *                        spectra measured in the FLAT files. This is the AUTO option
 *                        for High spectral resolution.
 * @return The gravi_data with FLATs and PROFILE maps
 * 
 * \exception CPL_ERROR_NULL_INPUT input data is missing
 * \exception CPL_ERROR_ILLEGAL_INPUT not all shutter are opened or profile_width
 * option not > 0 or missing table in the input data
 *
 * For each region defined in the IMAGING_DETECTOR_SC table this function
 * retrieves the profile of the spectrum passing by the reference point
 * defined in the table IMAGING_DETECTOR_SC. The FT data have no profile.
 */
/*----------------------------------------------------------------------------*/

gravi_data * gravi_compute_profile(gravi_data ** flats_data,
								   gravi_data * dark_map, gravi_data * bad_map,
								   int nflat, const cpl_parameterlist * params)
{
	/* Verbose */
	gravi_msg_function_start(1);
	cpl_ensure (flats_data, CPL_ERROR_NULL_INPUT, NULL);
	cpl_ensure (dark_map,   CPL_ERROR_NULL_INPUT, NULL);

    /* Ensure all beam open once */
    cpl_ensure (gravi_data_check_shutter_beam (flats_data, nflat),
                CPL_ERROR_ILLEGAL_INPUT, NULL);

	/* Construct the product */
	gravi_data * out_data = gravi_data_new (0);
    cpl_propertylist * out_header = gravi_data_get_header (out_data);
    
    /* Dump full header of first RAW data */
	gravi_data * raw0 = flats_data[0];
    cpl_propertylist_append (out_header, gravi_data_get_header (flats_data[0]));

	/* Copy IMAGING_DETECTOR into product */
	gravi_data_copy_ext (out_data, raw0, GRAVI_IMAGING_DETECTOR_FT_EXT);
	gravi_data_copy_ext (out_data, raw0, GRAVI_IMAGING_DETECTOR_SC_EXT);

	/* Create the FLUX array*/
	cpl_array * array_flux_FT = cpl_array_new(nflat, CPL_TYPE_DOUBLE);
	cpl_array * array_flux_SC = cpl_array_new(nflat, CPL_TYPE_DOUBLE);

    /*
     * (1) Compute the FLAT of the FT
     */
    if (!gravi_data_has_extension (flats_data[0], GRAVI_IMAGING_DATA_FT_EXT)) {
        cpl_msg_warning (cpl_func,"The FLAT data has no IMAGING_DATA_FT");
    }
    else
    {
	    cpl_msg_info (cpl_func, "Computing the FLAT of FT");

        /* Get the dark of FT */
        cpl_table * darkft_table;
        darkft_table = gravi_data_get_table (dark_map,  GRAVI_IMAGING_DATA_FT_EXT);
        CPLCHECK_NUL ("Cannot get DARK data of FT");

        /* Get the DARK as an image */
        cpl_image * darkft_img;
        darkft_img = gravi_image_from_column (darkft_table, "PIX", 0);

        /* Collapse each FLAT file */
        cpl_imagelist *imglist_ft = cpl_imagelist_new ();
        for (int file = 0; file < nflat; file++) {

            cpl_table * dataft_table;
            dataft_table = gravi_data_get_table (flats_data[file],  GRAVI_IMAGING_DATA_FT_EXT);
            CPLCHECK_NUL ("Cannot get data");

            cpl_imagelist * imglistft_tmp;
            imglistft_tmp = gravi_imagelist_wrap_column (dataft_table, "PIX");

            /* Remove DARK */
            cpl_imagelist_subtract_image (imglistft_tmp, darkft_img);

            /* Collapse the DITs of this FLAT */
            cpl_imagelist_set (imglist_ft, cpl_imagelist_collapse_create (imglistft_tmp), file);
            gravi_imagelist_unwrap_images (imglistft_tmp);

            /* compute the flux of this flat*/
            cpl_array_set(array_flux_FT, file, cpl_image_get_flux(cpl_imagelist_get(imglist_ft, file)));
        }

        /* Collapse the FLAT files for all telescopes together */
        cpl_image * flatft_img = cpl_imagelist_collapse_create (imglist_ft);
        cpl_imagelist_delete (imglist_ft);

        /* Create the flat_table */
        cpl_table * flatft_table = cpl_table_extract (gravi_data_get_table (flats_data[nflat-1],
                                   GRAVI_IMAGING_DATA_FT_EXT), 0, 1);
        
        cpl_array * flatft_array = gravi_array_wrap_image (flatft_img);
        cpl_table_set_array (flatft_table, "PIX", 0, flatft_array);

        /* Remove median image and array */
        FREE (cpl_array_unwrap, flatft_array);
        FREE (cpl_image_delete, flatft_img);
        FREE (cpl_image_delete, darkft_img);

        /* Set the FLAT as IMAGING_DATA into product */
        gravi_data_add_table (out_data, NULL, GRAVI_IMAGING_DATA_FT_EXT, flatft_table);
        
    } /* End FLAT of FT */

    /* 
     * (2) General data used by FLAT and PROFILE for SC
     */
    
	/* Get that the detector_table */
    cpl_table * detector_table;
	detector_table = gravi_data_get_table (flats_data[0], GRAVI_IMAGING_DETECTOR_SC_EXT);
	cpl_ensure (detector_table, CPL_ERROR_ILLEGAL_INPUT, NULL);

	/* Get that the header of first file and of dark */
    cpl_propertylist * dark_header = gravi_data_get_header (dark_map);
	cpl_propertylist * flat0_header = gravi_data_get_header (flats_data[0]);

    /* Get necessary information */
    int det_startx = gravi_pfits_get_window_start (flat0_header);
	const char * resolution = gravi_pfits_get_resolution (flat0_header);
	const char * pola_mode = gravi_pfits_get_pola_mode(flat0_header, GRAVI_SC);
	int nb_region = cpl_table_get_nrow (detector_table);
    
	/* Get the DARK and BAD data */
    cpl_image * dark_img, * bad_img;
	dark_img   = gravi_data_get_img (dark_map, GRAVI_IMAGING_DATA_SC_EXT);
	bad_img = gravi_data_get_img (bad_map, GRAVI_IMAGING_DATA_SC_EXT);
    cpl_ensure (dark_img && bad_img, CPL_ERROR_ILLEGAL_INPUT, NULL);
    
    /* 
     * (3) Compute the FLAT for the SC 
     */
    cpl_msg_info (cpl_func, "Computing the FLAT of SC");
    
    /* Imagelist to store the collapsed FLATs with all telescopes */
	cpl_imagelist * temp_imglist = cpl_imagelist_new ();

	for (int file = 0; file < nflat; file++){
        
		/* Extract necessary parameters and construct the output table */
        cpl_imagelist * data_imglist;
		data_imglist = gravi_data_get_cube (flats_data[file], GRAVI_IMAGING_DATA_SC_EXT);

		/* Extract data with DARK, BADPIX  in  [ADU] */
        data_imglist = cpl_imagelist_duplicate (data_imglist);
        cpl_imagelist_subtract_image (data_imglist, dark_img);
        gravi_remove_badpixel_sc (data_imglist, bad_img);

        /* Collapse the DITs of this FLAT */
        cpl_image * collapsed_img = cpl_imagelist_collapse_sigclip_create (data_imglist, 5, 5, 0.6, CPL_COLLAPSE_MEDIAN_MEAN, NULL);
		FREE (cpl_imagelist_delete, data_imglist);

        /* Save this FLAT in the imagelist to collapse them */
        cpl_imagelist_set (temp_imglist, collapsed_img,
                           cpl_imagelist_get_size (temp_imglist));

        /* compute the flux of this flat*/
        cpl_array_set(array_flux_SC, file, cpl_image_get_flux(collapsed_img));

		CPLCHECK_NUL ("Error");
	} /* End loop on FLATs*/

    /* Collapse the FLATs for all telescopes together */
    cpl_image * allflat_img;
	allflat_img = cpl_imagelist_collapse_create (temp_imglist);
	FREE (cpl_imagelist_delete, temp_imglist);

	CPLCHECK_NUL ("Cannot collapse FLATs");

    /* Get the extension (illumated part) of the combined FLAT */
    int * ext_dim = gravi_image_extract_dimension (allflat_img);
    int fullstartx = ext_dim[0] + det_startx - 2;

    /* Crop the FLAT into these new dimensions 
     * FIXME: would be better to no crop FLAT */
    cpl_image * flatsc_img;
    flatsc_img = cpl_image_extract (allflat_img, ext_dim[0], 1,
                                    ext_dim[0] + ext_dim[1] - 1,
                                    cpl_image_get_size_y (allflat_img));
	FREE (cpl_image_delete, allflat_img);

    /* Set the cropped FLAT of SC in output data 
     * STARTX is in FITS convention (start 1) while
     * FULLSTARTX seems in C convention (start 0) */
	cpl_propertylist * flat_plist = cpl_propertylist_new ();
	cpl_propertylist_update_int (flat_plist, PROFILE_FULLSTARTX, fullstartx);
	cpl_propertylist_update_int (flat_plist, PROFILE_STARTX, ext_dim[0]);
	cpl_propertylist_update_int (flat_plist, PROFILE_NX, ext_dim[1]);
    
	gravi_data_add_img (out_data, flat_plist, GRAVI_IMAGING_DATA_SC_EXT, flatsc_img);
	CPLCHECK_NUL ("Cannot set FLAT");


    /* 
     * (4) Compute the PROFILE for the SC 
     */
    cpl_msg_info (cpl_func, "Computing the PROFILE of SC");


    /* List of image to store each FLAT */
	cpl_image ** median_img = cpl_calloc (4, sizeof(cpl_image *));

	/* kernel for median filter. Use more 
     * pixels (~15) in HIGH, to skip the cluster */
    cpl_size size_kernel = !strcmp (resolution, "HIGH") ? 15 : 5;
    cpl_msg_info (cpl_func, "Median filtering over %lld spectral pixels", size_kernel);
    cpl_mask * kernel = cpl_mask_new (size_kernel, 1);
    cpl_mask_not (kernel);

	for (int file = 0; file < nflat; file++){
		/* Get the primary header of each file */
        cpl_propertylist * flat_header;
		flat_header = gravi_data_get_header (flats_data[file]);

		/* Extract necessary parameters and construct the output table */
		cpl_imagelist * data_imglist;
        data_imglist = gravi_data_get_cube (flats_data[file], GRAVI_IMAGING_DATA_SC_EXT);
		cpl_ensure (flat_header && data_imglist, CPL_ERROR_NULL_INPUT, NULL);

		/* Extract data with DARK, BADPIX  in  [ADU] */
        data_imglist = cpl_imagelist_duplicate (data_imglist);
        cpl_imagelist_subtract_image (data_imglist, dark_img);
        gravi_remove_badpixel_sc (data_imglist, bad_img);
        
        /* Collapse the DITs of this FLAT */
        cpl_image * collapsed_img = cpl_imagelist_collapse_sigclip_create (data_imglist, 5, 5, 0.6, CPL_COLLAPSE_MEDIAN_MEAN, NULL);
        
		FREE (cpl_imagelist_delete, data_imglist);

        /* Create a filtered version of this FLAT */
        cpl_image * filtered_img = cpl_image_duplicate (collapsed_img);
        cpl_image_filter_mask (filtered_img, collapsed_img, kernel,
                               CPL_FILTER_MEDIAN, CPL_BORDER_FILTER);
        FREE (cpl_image_delete, collapsed_img);

        /* Crop it */
        cpl_image * crop_img;
        crop_img = cpl_image_extract (filtered_img, ext_dim[0], 1,
                                      ext_dim[0] + ext_dim[1] - 1,
                                      cpl_image_get_size_y (filtered_img));
        FREE (cpl_image_delete, filtered_img);
        
		/* Save this filtered version in the median_img[] 
         * According to the beam open */
        int id = gravi_get_shutter_id (flat_header);
        median_img[id] = crop_img;

		CPLCHECK_NUL ("Error");
	} /* End loop on FLATs*/
    
	FREE (cpl_mask_delete, kernel);
    
	/* Get the image dimensions after crop */
	int nx = cpl_image_get_size_x (median_img[0]);
	int ny = cpl_image_get_size_y (median_img[0]);

	/* Get the profile width parameter */
	int profile_width = gravi_param_get_int (params, "gravity.calib.profile-width");
    cpl_ensure (profile_width > 0, CPL_ERROR_ILLEGAL_INPUT, NULL);
    
	/* Compute the size_profile. For LOW and MEDIUM it is automatic.
     * For HIGH it comes from the option */
    int n_darkline;
    cpl_propertylist * flat_header;
    int window_mode;
	flat_header = gravi_data_get_header (flats_data[0]);
    int size_profile;
    /* if the windowing mode if after the 05/12/2016 */
    if ( cpl_propertylist_get_float(flat_header, "MJD-OBS") > 57728 ){
        n_darkline= 0;
        window_mode = 1;
        cpl_msg_info(cpl_func, "Windowing after MJD = 57728");
		if (! (strcmp(resolution, "LOW") && strcmp(resolution, "MED")) &&  !strcmp(pola_mode, "COMB"))	{
			size_profile = ny/nb_region;
			cpl_msg_info(cpl_func, "Use a computed size_profile of %d", size_profile);
		}
		else{
			size_profile = (fmod(profile_width, 2) == 0) ? profile_width + 1 : profile_width;
			cpl_msg_info(cpl_func, "Use a given size profile of %d", size_profile);
		}
    }
    /* else keep old windowing for backward compatibility */
    else {
        n_darkline= 1;
        window_mode = 0;
		if (! (strcmp(resolution, "LOW") && strcmp(resolution, "MED")))	{
			size_profile = (ny-(nb_region+1)*n_darkline)/nb_region;
			cpl_msg_info(cpl_func, "Use a computed size_profile of %d", size_profile);
		}
		else{
			size_profile = (fmod(profile_width, 2) == 0) ? profile_width + 1 : profile_width;
			cpl_msg_info(cpl_func, "Use a given size profile of %d", size_profile);
		}
    }

	/* Construction of the PROFILE_DATA table */
	cpl_table * profile_table = cpl_table_new (1);

	/* Create the DATA# columns with their dimension */
	cpl_array * dimension = cpl_array_new (2, CPL_TYPE_INT);
	cpl_array_set (dimension, 0, nx);
	cpl_array_set (dimension, 1, ny);
	for (int region = 0; region < nb_region; region ++){
        const char * data = GRAVI_DATA[region];
		cpl_table_new_column_array (profile_table, data,
									CPL_TYPE_DOUBLE, nx * ny);
		cpl_table_set_column_dimensions (profile_table, data, dimension);
	}
	FREE (cpl_array_delete, dimension);

	/* Construction of the PROFILE_PARAMS table */
	cpl_table * params_table = cpl_table_new (nb_region);
	cpl_table_new_column_array (params_table, "CENTERY", CPL_TYPE_DOUBLE, nx);
	cpl_table_new_column_array (params_table, "WIDTH", CPL_TYPE_DOUBLE, nx);

    
    /* Construction of a map used to zero the badpixels
     * Mask is:     good pixels = 1  ;  bad pixels = 0    */
    cpl_image * mask_img = NULL;

    if ( !gravi_param_get_bool (params,
          "gravity.calib.force-badpix-to-zero") )
    {
        cpl_msg_info (cpl_func, "Bad pixels are *not*"
                      "forced to zero in profiles");
        cpl_image_fill_window (mask_img, 1, 1, nx, ny, 1.0);
    }
    else {
        cpl_msg_info (cpl_func, "Bad pixels are "
                      "forced to zero in profiles");
        mask_img = cpl_image_extract (bad_img, ext_dim[0], 1,
                                      ext_dim[0] + ext_dim[1] - 1, ny);
        cpl_image_threshold (mask_img, 0.5, 0.5, 1.0, 0.0);
    }

    /* Update mask from FLAT itself in LOW,
     * FIXME: to decide what we zero exactly. */
	if ( !strcmp (resolution, "LOW") )
    {
        cpl_msg_info (cpl_func, "Pixels with low FLAT values"
                      " are forced to zero in profiles.");
        
        double threshold = cpl_propertylist_get_double (dark_header, QC_DARKRMS_SC);
        cpl_image * mask2_img = cpl_image_duplicate (flatsc_img);
        cpl_image_threshold (mask2_img, threshold, threshold, 0.0, 1.0);
        cpl_image_multiply (mask_img, mask2_img);
        FREE (cpl_image_delete, mask2_img);
    }

    /* Define the mode */
    const char * mode = gravi_param_get_string_default (params,
                        "gravity.calib.profile-mode", "AUTO");
    
    if (!strcmp(mode, "AUTO")) {
        if (!strcmp(resolution, "LOW"))  mode = "PROFILE";
        if (!strcmp(resolution, "MED"))  mode = "PROFILE";
        if (!strcmp(resolution, "HIGH")) mode = "BOX";
    }

    cpl_msg_info (cpl_func, "Profile computed with mode: %s  (%s)", mode, resolution);
    
    /* Loop on regions */
    for (int region = 0; region < nb_region ; region++) {

        int tel_1 = gravi_region_get_tel (detector_table, region, 0);
        int tel_2 = gravi_region_get_tel (detector_table, region, 1);
        CPLCHECK_NUL ("Cannot get the telescope from region");
        
        /* Compute the mean image between the first
         * and second telescope */
        cpl_image * mean_img;
        mean_img = cpl_image_add_create (median_img[tel_1],
                                         median_img[tel_2]);
        cpl_image_divide_scalar (mean_img, 2);

        /* Get the reference coordinates of this region. That is
         * a point where the spectra is supposed to be for sure 
         * In unit of the full detector window */
        int ref0_x = gravi_table_get_value (detector_table, "CENTER", region, 0);
        int ref0_y = gravi_table_get_value (detector_table, "CENTER", region, 1);

        /* Convert the reference coordinates in unit of the cropped data 
         * Actually we force ref_x and ref_y in LOW and MED. */
        int ref_x, ref_y;

        /* if the windowing mode if after the 05/12/2016 */
        if ( window_mode == 1 ){
            if (! (strcmp(resolution, "LOW") && strcmp(resolution, "MED")) &&  !strcmp(pola_mode, "COMB") )	{
                ref_x = nx / 2;
                ref_y = (region+1)*(size_profile)-size_profile/2;
            } else  {
                //ref_x = ref0_x - det_startx - (1 + ext_dim[0]);
                ref_x = ref0_x - (1 + ext_dim[0]);
                ref_y = ref0_y;
            }
        }
        else {
            if (! (strcmp(resolution, "LOW") && strcmp(resolution, "MED")) )	{
                ref_x = nx / 2;
                ref_y = (region+1)*(size_profile+n_darkline)-size_profile/2;
            } else  {
                //ref_x = ref0_x - det_startx - (1 + ext_dim[0]);
                ref_x = ref0_x - (1 + ext_dim[0]);
                ref_y = ref0_y;
            }
        }

        /*
         * Fit Gaussian to profile of each spectral channel
         */
        cpl_vector * values_x0 = cpl_vector_new (nx);
        cpl_vector * values_y0 = cpl_vector_new (nx);
        cpl_vector * values_sigma = cpl_vector_new (nx);
        
        gravi_fit_profile (values_x0, values_y0, values_sigma,
                           mean_img, ref_x, ref_y, size_profile,
                           resolution);
        CPLCHECK_NUL ("Cannot fit data into profile params");
        
        /* Fill the PROFILE_PARAM table for this region */
        cpl_array * values_arr;
        values_arr = cpl_array_wrap_double (cpl_vector_get_data (values_y0), nx);
        cpl_table_set_array (params_table, "CENTERY", region, values_arr);
        cpl_array_unwrap (values_arr);
        
        values_arr = cpl_array_wrap_double (cpl_vector_get_data(values_sigma), nx);
        cpl_table_set_array (params_table, "WIDTH", region, values_arr);
        cpl_array_unwrap (values_arr);

        
        /*
         * Create and compute the profile image of this region
         */
        
        /* Define the spatial pixels over which we fill the profile */
        cpl_size iy_min = 0, iy_max = ny-1;
        if (! (strcmp(resolution, "LOW") && strcmp(resolution, "MED")) ){
            iy_min = ref_y - size_profile/2;
            iy_max = ref_y + size_profile/2-1;
        }

        /* Fill the profile image */
        cpl_image * region_img;
        region_img = gravi_create_profile_image (mean_img, values_x0,
                                                 values_y0, values_sigma,
                                                 iy_min, iy_max, mode);
        FREE (cpl_vector_delete, values_x0);
        FREE (cpl_vector_delete, values_y0);
        FREE (cpl_vector_delete, values_sigma);
        CPLCHECK_NUL ("Cannot build profile image");

        /* Force some pixel to zero with computed mask */
        cpl_image_multiply (region_img, mask_img);
        
        /* Fill the profile_table for this region */
        cpl_array * array = gravi_array_wrap_image (region_img);
        cpl_table_set_array (profile_table, GRAVI_DATA[region], 0, array);
        FREE (cpl_array_unwrap, array);
        FREE (cpl_image_delete, region_img);

        /* Free the mean image of the 2 flats of this region */
        FREE (cpl_image_delete, mean_img);
        
	} /* End loop on regions */

	/* Save the PROFILE_DATA table */
	cpl_propertylist * profile_plist = cpl_propertylist_new ();
    cpl_propertylist_copy_property (profile_plist, flat_plist, PROFILE_FULLSTARTX);
    cpl_propertylist_copy_property (profile_plist, flat_plist, PROFILE_STARTX);
    cpl_propertylist_copy_property (profile_plist, flat_plist, PROFILE_NX);
	gravi_data_add_table (out_data, profile_plist,
                          GRAVI_PROFILE_DATA_EXT, profile_table);
	
	/* Save the PROFILE_PARAMS table */
	cpl_propertylist * params_plist = cpl_propertylist_duplicate (profile_plist);
	gravi_data_add_table (out_data, params_plist,
                          GRAVI_PROFILE_PARAMS_EXT, params_table);
    

	/* 
     * (5) Add the QC parameter of the lateral positioning of the first region 
     */
    
	cpl_propertylist * main_header = gravi_data_get_header (out_data);
	double qc_value;
	char qc_name[100];
	cpl_size idx;

	for (int file = 0; file < nflat; file++){
		  sprintf (qc_name, "ESO QC FLATFLUX FT%i", file+1);
		  qc_value = cpl_array_get(array_flux_FT, file, NULL);
		  cpl_propertylist_append_double (main_header, qc_name, qc_value);
		  cpl_propertylist_set_comment (main_header, qc_name, "[ADU] Total flux per DIT");
		  cpl_msg_info (cpl_func,"%s = %f [ADU]", qc_name, qc_value);
	}
	for (int file = 0; file < nflat; file++){
		  sprintf (qc_name, "ESO QC FLATFLUX SC%i", file+1);
		  qc_value = cpl_array_get(array_flux_SC, file, NULL);
		  cpl_propertylist_append_double (main_header, qc_name, qc_value);
		  cpl_propertylist_set_comment (main_header, qc_name, "[ADU] Total flux per DIT");
		  cpl_msg_info (cpl_func,"%s = %f [ADU]", qc_name, qc_value);
	}

	for (int reg = 0; reg < nb_region; reg += 12) {

	  /* Median lateral position */
	  qc_value = cpl_array_get_median (cpl_table_get_array (params_table, "CENTERY", reg));
	  qc_value = floor (qc_value * 1e6) * 1e-6;
	  sprintf (qc_name, "ESO QC PROFILE_CENTER SC%i MED", reg+1);
	  cpl_propertylist_append_double (main_header, qc_name, qc_value);
	  cpl_propertylist_set_comment (main_header, qc_name, "[pixel] position of region");
	  cpl_msg_info (cpl_func,"%s = %f [pixel]", qc_name, qc_value);

	  /* Median width */
	  qc_value = cpl_array_get_median (cpl_table_get_array (params_table, "WIDTH", reg));
	  qc_value = floor (qc_value * 1e6) * 1e-6;
	  sprintf (qc_name, "ESO QC PROFILE_WIDTH SC%i MED", reg+1);
	  cpl_propertylist_append_double (main_header, qc_name, qc_value);
	  cpl_propertylist_set_comment (main_header, qc_name, "[pixel] width of region");
	  cpl_msg_info (cpl_func,"%s = %f [pixel]", qc_name, qc_value);

	  /* Lateral position in the left */
	  idx = 1 * ext_dim[1] / 6;
	  qc_value = cpl_array_get (cpl_table_get_array (params_table, "CENTERY", reg), idx, NULL);
	  qc_value = floor (qc_value * 1e6) * 1e-6;
	  sprintf (qc_name, "ESO QC PROFILE_CENTER SC%i LEFT", reg+1);
	  cpl_propertylist_append_double (main_header, qc_name, qc_value);
	  cpl_propertylist_set_comment (main_header, qc_name, "[pixel] at STARTX+NX/6");
	  cpl_msg_info (cpl_func,"%s = %f [pixel] for x=%lld", qc_name, qc_value, idx);

	  /* Lateral position in the right */
	  idx = 5 * ext_dim[1] / 6;
	  qc_value = cpl_array_get (cpl_table_get_array (params_table, "CENTERY", reg), idx, NULL);
	  qc_value = floor (qc_value * 1e6) * 1e-6;
	  sprintf (qc_name, "ESO QC PROFILE_CENTER SC%i RIGHT", reg+1);
	  cpl_propertylist_append_double (main_header, qc_name, qc_value);
	  cpl_propertylist_set_comment (main_header, qc_name, "[pixel] at STARTX+5*NX/6");
	  cpl_msg_info (cpl_func,"%s = %f [pixel] for x=%lld", qc_name, qc_value, idx);

	  CPLCHECK_NUL ("Cannot compute QC parameters");
	}

	/* Deallocation of the variables */
    FREE (cpl_image_delete, mask_img);    
    FREELOOP (cpl_image_delete, median_img, 4);
    FREE (cpl_free, ext_dim);
    cpl_array_delete(array_flux_FT);
    cpl_array_delete(array_flux_SC);

	/* Verbose */
	gravi_msg_function_exit(1);
	return out_data;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Compute mean detector gain
 * 
 * @param flats_data    The input raw FLAT data, one per beam
 * @param nrawgain      4
 * @param dark_map      The input DARK calibration map
 *
 * \exception CPL_ERROR_NULL_INPUT input data is missing
 * \exception CPL_ERROR_ILLEGAL_INPUT not all shutter are opened
 *
 * It returs a propertylist with the QC value of the mean SC and mean FT
 * gain, in [ADU/e].
 */
/*----------------------------------------------------------------------------*/

cpl_propertylist * gravi_compute_gain (gravi_data ** flats_data,
                                       int nrawgain,
                                       gravi_data * dark_map)
{
	int nv;
	const cpl_size maxdeg = 1, mindeg = 0;
	const cpl_size slope_deg = 1;

	/* Verbose */
	gravi_msg_function_start(1);
	cpl_ensure (flats_data, CPL_ERROR_NULL_INPUT, NULL);
	cpl_ensure (dark_map,   CPL_ERROR_NULL_INPUT, NULL);
    
    /* Ensure all beam open once */
    cpl_ensure (gravi_data_check_shutter_beam (flats_data, nrawgain),
                CPL_ERROR_ILLEGAL_INPUT, NULL);
	
	/* Get a sequence of 4 single-shutter open files 
     * FIXME: remove this part isn't it ?? */
	cpl_msg_info (cpl_func, "Search for 4-shutter sequence");
	
	gravi_data ** gain_file = cpl_calloc(4,sizeof(gravi_data*));
    
	for (int i = 0; i < nrawgain; i++){
		cpl_propertylist * flat_header = gravi_data_get_header (flats_data[i]);

		if (gravi_check_shutter (flat_header, 1,0,0,0)) {
		  gain_file[0] = (flats_data[i]);
		}
		if (gravi_check_shutter (flat_header, 0,1,0,0)) {
		  gain_file[1] = (flats_data[i]);
		}
		if (gravi_check_shutter (flat_header, 0,0,1,0)) {
		  gain_file[2] = (flats_data[i]);
		}
		if (gravi_check_shutter (flat_header, 0,0,0,1)) {
		  gain_file[3] = (flats_data[i]);
		}
	}

    /* 
     * Create the output. It is a plist only since we
     * don't create a gravi_data
     */
    cpl_propertylist * output_plist = cpl_propertylist_new();

    
    /*
     *
     * Compute the gain for SC
     *
     */
    if (!gravi_data_has_extension (gain_file[0], GRAVI_IMAGING_DATA_SC_EXT)) {
        cpl_msg_warning (cpl_func,"The FLAT data has no IMAGING_DATA_SC");
    }
    else
    {
	    cpl_msg_info (cpl_func, "Computing the gain of SC");

        /* Get the size of the image */
		cpl_image * image = gravi_data_get_img (gain_file[0], GRAVI_IMAGING_DATA_SC_EXT);
		cpl_size size = cpl_image_get_size_x (image) * cpl_image_get_size_y (image);

        /* Build matrix and vector to fit var = f(mean) */
		cpl_vector * vector_var = cpl_vector_new (4*size);
		cpl_matrix * matrix_mean = cpl_matrix_new (1, 4*size);

		/* Get a pointer to the dark image */
		cpl_msg_info (cpl_func,"DARK and BADPIX of SC are not used");

		/* Loop on files to fill the variance and mean vectors */
		for (int file = 0; file < 4; file++) {
		  
			/* Get the imagelist and property lists */
            cpl_imagelist * data_imglist;
			data_imglist = gravi_data_get_cube (gain_file[file], GRAVI_IMAGING_DATA_SC_EXT);

            /* Compute <image> and <image^2> */
			cpl_imagelist * imglist1 = cpl_imagelist_new ();
			cpl_imagelist * imglist2 = cpl_imagelist_new ();
			for (cpl_size i = 0; i < cpl_imagelist_get_size (data_imglist); i ++){
				image = cpl_image_cast (cpl_imagelist_get (data_imglist, i), CPL_TYPE_DOUBLE);
				cpl_imagelist_set (imglist1, image, i);
				cpl_imagelist_set (imglist2, cpl_image_power_create (image, 2), i);
			}
	        cpl_image * image1_mean = cpl_imagelist_collapse_create (imglist1);
	        cpl_image * image2_mean = cpl_imagelist_collapse_create (imglist2);
			FREE (cpl_imagelist_delete, imglist1);
            FREE (cpl_imagelist_delete, imglist2);

			/* Loop on pixels to fill the vector and matrix structures
			 * used latter in the PTC fit */
            cpl_size nx = cpl_image_get_size_x (image1_mean);
            cpl_size ny = cpl_image_get_size_y (image1_mean);
			for (cpl_size i = 0; i <  nx; i ++) {
				for (cpl_size j = 0; j <  ny; j ++) {
					cpl_vector_set (vector_var, file * size + (j + ny*i),
                                    cpl_image_get (image2_mean, i+1, j+1, &nv) -
                                    pow (cpl_image_get (image1_mean, i+1, j+1, &nv), 2));
					cpl_matrix_set (matrix_mean, 0, file * size + (j + ny*i),
                                    cpl_image_get (image1_mean, i+1, j+1, &nv));
				}
			}
			/* End loop on pixels */

			/* Deallocation of variables */
			FREE (cpl_image_delete, image1_mean);
			FREE (cpl_image_delete, image2_mean);
		}
		/* End loop on files */

        /* Fit the curve variance versus the median,
           thus the output gains are in ADU/e */
        cpl_polynomial * fit_slope = cpl_polynomial_new (1);
        cpl_polynomial_fit (fit_slope, matrix_mean, NULL, vector_var, NULL,
                           CPL_FALSE, &mindeg, &maxdeg);

        /* Get the slope */
        double slope = cpl_polynomial_get_coeff (fit_slope, &slope_deg);
        cpl_msg_info (cpl_func, "mean gain SC = %.4f [adu/e-] Mean gain of detector", slope);
        cpl_propertylist_append_double (output_plist, QC_MEANGAIN_SC, slope);
        cpl_propertylist_set_comment (output_plist, QC_MEANGAIN_SC, "[adu/e-] Mean gain of SC detector" );
        
        /* Delete */
        cpl_vector_delete (vector_var);
        cpl_matrix_delete (matrix_mean);
        cpl_polynomial_delete (fit_slope);
	}
	/* End SC case */

    /*
     *
     * Compute the gain for FT
     *
     */
    if (!gravi_data_has_extension (gain_file[0], GRAVI_IMAGING_DATA_SC_EXT)) {
        cpl_msg_warning (cpl_func,"The FLAT data has no IMAGING_DATA_FT");
    }
    else
    {
	    cpl_msg_info (cpl_func, "Computing the gain of FT");

        /* Get DARK image */
        cpl_table * dark_table;
		dark_table = gravi_data_get_table (dark_map, GRAVI_IMAGING_DATA_FT_EXT);
        cpl_image * dark_img = gravi_image_from_column (dark_table, "PIX", 0);

        /* Get the size of the image */
		cpl_size size = cpl_table_get_column_depth (dark_table, "PIX");

        /* Build matrix and vector to fit var = f(mean) */
		cpl_vector * vector_var = cpl_vector_new (4 * size);
		cpl_matrix * matrix_mean = cpl_matrix_new (1, 4 * size);

		/* Compute the gain for each file */
		for (int file = 0; file < 4; file++) {
            
			/* Get the tables */
            cpl_table * data_table;
            data_table = gravi_data_get_table (gain_file[file], GRAVI_IMAGING_DATA_FT_EXT);
            cpl_table_cast_column (data_table, "PIX", "PIX", CPL_TYPE_DOUBLE);
            cpl_imagelist * data_imglist = gravi_imagelist_from_column (data_table,"PIX");
            CPLCHECK_NUL ("Cannot get data");

            /* Remove DARK */
	        cpl_imagelist_subtract_image (data_imglist, dark_img);

            /* Compute MEAN */
	        cpl_image * mean_img = cpl_imagelist_collapse_create (data_imglist);

            /* Compute VARIANCE */
	        cpl_imagelist_subtract_image (data_imglist, mean_img);
            cpl_imagelist_power (data_imglist, 2.0);
            cpl_image * var_img = cpl_imagelist_collapse_create (data_imglist);
            
            FREE (cpl_imagelist_delete, data_imglist);

			/* Loop on pixels to fill the vector and matrix structures
			 * used latter in the PTC fit */
            cpl_size nx = cpl_image_get_size_x (mean_img);
            cpl_size ny = cpl_image_get_size_y (mean_img);
			for (cpl_size i = 0; i <  nx; i ++) {
				for (cpl_size j = 0; j <  ny; j ++) {
					cpl_vector_set (vector_var, file * size + (j + ny*i),
                                    cpl_image_get (var_img, i+1, j+1, &nv));
					cpl_matrix_set (matrix_mean, 0, file * size + (j + ny*i),
                                    cpl_image_get (mean_img, i+1, j+1, &nv));
				}
			}
			/* End loop on pixels */
            
			/* Deallocation of variables */
			FREE (cpl_image_delete, var_img);
			FREE (cpl_image_delete, mean_img);
		} /* End loop on files */

        /* Fit the curve variance versus the median,
           thus the output gains are in ADU/e */
        cpl_polynomial * fit_slope = cpl_polynomial_new (1);
        cpl_polynomial_fit (fit_slope, matrix_mean, NULL, vector_var, NULL,
                            CPL_FALSE, &mindeg, &maxdeg);

        /* Get the slope */
        double slope = cpl_polynomial_get_coeff (fit_slope, &slope_deg);
        cpl_msg_info (cpl_func, "mean gain FT = %.4f [adu/e-] Mean gain of detector", slope);
        cpl_propertylist_append_double (output_plist, QC_MEANGAIN_FT, slope);
        cpl_propertylist_set_comment (output_plist, QC_MEANGAIN_FT, "[adu/e-] Mean gain of FT detector");
        
        /* Delete */
        FREE (cpl_vector_delete, vector_var);
        FREE (cpl_matrix_delete, matrix_mean);
        FREE (cpl_polynomial_delete, fit_slope);
        FREE (cpl_image_delete, dark_img);
        
	} /* End FT case*/

    cpl_free (gain_file);

	/* Verbose */
	gravi_msg_function_exit(1);
	return output_plist;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Identify the bad pixels in the DARK map and create the BAD map.
 * 
 * @param dark_map 	   The input dark map calibration
 * @param flats_data   The input raw flats (optional)
 * @param nflat        The number of input flats (shall be 4)
 * @param params       Input parameter list with :
 *                     - bad-dark-threshold : the rms factor for dark bad
 *                     pixel threshold.
 * 
 * @return The BAD map with the detected bad pixels.
 * 
 * \exception CPL_ERROR_NULL_INPUT input data is missing
 * \exception CPL_ERROR_ILLEGAL_INPUT The number of input flats is not 4
 *
 * Pixel with dark value or read-out-noise out of a specified
 * range are declared as bad pixels.
 *
 * If the flats_data are provided, they are collapsed together.
 * Pixel with mean flat value lower than a horizontal median filtering
 * are declared as bad pixels. Only pixel with a mean flat value higher
 * than 100 adu are inspected. The mask of inspected pixels
 * is saved in IMAGING_MASK_SC. This can be applied with normal
 * *or* defocused FLATs. Still to be validated in LOW.
 * 
 * FIXME: in term of implementation, this function is weird as
 * some of the computation are somehow done in-place after duplication
 * of the data...
 */
/*----------------------------------------------------------------------------*/

gravi_data * gravi_compute_badpix (gravi_data * dark_map,
                                   gravi_data ** flats_data,
                                   int nflat,
                                   const cpl_parameterlist * params)
{
	/* Verbose */
	gravi_msg_function_start(1);
	cpl_ensure (dark_map,   CPL_ERROR_NULL_INPUT, NULL);
	cpl_ensure (params,     CPL_ERROR_NULL_INPUT, NULL);
    cpl_ensure (nflat==4 || nflat==0, CPL_ERROR_ILLEGAL_INPUT, NULL);

	/* Construction of the bad pixels map */
    gravi_data * bad_map = gravi_data_new(0);
	cpl_propertylist * bad_header = gravi_data_get_header (bad_map);

    /* Dump full header into product */
	cpl_propertylist * dark_header = gravi_data_get_header (dark_map);
    cpl_propertylist_append (bad_header, dark_header);

    /* 
     * Compute bad pixels of FT 
     */

    if (!gravi_data_has_extension (dark_map, GRAVI_IMAGING_DATA_FT_EXT)) {
        cpl_msg_warning (cpl_func,"The DARK map has no IMAGING_DATA_FT");
    }
    else
    {
	    cpl_msg_info (cpl_func,"Compute BADPIXEL of FT");

        /* Copy necessary tables */
        gravi_data_copy_ext (bad_map, dark_map, GRAVI_IMAGING_DETECTOR_FT_EXT);
	  
	    /* This is the FT */
		cpl_table * dark_table = gravi_data_get_table (dark_map, GRAVI_IMAGING_DATA_FT_EXT);
        cpl_table * std_table  = gravi_data_get_table (dark_map, GRAVI_IMAGING_ERR_FT_EXT);

        /* Verify it has one single row */
        cpl_size n_row = cpl_table_get_nrow (dark_table);
        cpl_ensure (n_row == 1, CPL_ERROR_ILLEGAL_INPUT, NULL);
        cpl_ensure (n_row == cpl_table_get_nrow (std_table), CPL_ERROR_ILLEGAL_INPUT, NULL);
        
		/* Get the rms factor for dark bad pixel threshold */
		int bad_dark_factor = gravi_param_get_int (params, "gravity.calib.bad-dark-threshold");

		CPLCHECK_NUL ("Cannot get data");

		/* Get the dark rms and the mean dark */
		double dark_rms = cpl_propertylist_get_double (dark_header, QC_DARKRMS_FT);
		double dark_mean = cpl_propertylist_get_double (dark_header, QC_MEANDARK_FT);
        CPLCHECK_NUL ("Cannot get QC");
        

        /* Compute the specified range are declared as bad pixels */
        double range_max = dark_mean + bad_dark_factor * dark_rms;
        cpl_msg_info (cpl_func,"FT threshold on mean value = %f",range_max);
        
        /* Get the min std to declare pixel as bad */
        double std_min = 1e-3;
        
        /* Get the max std to declare pixel as bad */
        cpl_size npix = cpl_table_get_column_depth (dark_table, "PIX");
        cpl_array * std_array  = cpl_table_get_data_array (std_table, "PIX")[0];
        cpl_vector * std_vector = cpl_vector_new (npix);
        for (cpl_size pix = 0; pix < npix; pix++) {
            cpl_vector_set (std_vector, pix, cpl_array_get (std_array, pix, NULL));
        }
        cpl_vector_sort (std_vector, CPL_SORT_ASCENDING);
        /* use 21/24 percentil to include the rms of the metrology*/
        cpl_size percentil = (int) npix*0.875;
        double std_max = cpl_vector_get (std_vector, percentil)*(1+bad_dark_factor/8.75);
        cpl_msg_info (cpl_func,"FT threshold on std value = %f",std_max);

		/* Compute the number of bad pixel */
		int count_bp_dark = 0;

        /* Create output table */
		cpl_table * bad_table = cpl_table_extract (dark_table, 0, 1);

        /* Loop on pixels */
        cpl_array * dark_array = cpl_table_get_data_array (dark_table, "PIX")[0];
        cpl_array * bad_array  = cpl_table_get_data_array (bad_table, "PIX")[0];

        for (cpl_size pix = 0; pix < npix; pix++) {
            cpl_array_set (bad_array, pix, 0);
            /* flag on the mean value of the dark */
            if (cpl_array_get (dark_array, pix, NULL) > range_max) {
                cpl_array_set (bad_array, pix, BADPIX_DARK);
                count_bp_dark ++;
                cpl_msg_info(cpl_func,"Detected a bad FT pixel at position %lli based on its mean flux: %f ADU", pix, cpl_array_get (dark_array, pix, NULL));
            } else if (cpl_array_get (std_array, pix, NULL) < std_min) {
                cpl_array_set (bad_array, pix, BADPIX_DARK);
                count_bp_dark ++;
            } else if (cpl_array_get (std_array, pix, NULL) > std_max) {
                cpl_array_set (bad_array, pix, BADPIX_DARK);
                count_bp_dark ++;
                cpl_msg_warning(cpl_func,"Detected a bad FT pixel at position %lli based on its rms: %f ADU", pix, cpl_array_get (std_array, pix, NULL));
            }
                
        }

        /* Set QC parameter */
		cpl_propertylist_append_int (bad_header, QC_BADPIX_FT, count_bp_dark);
		cpl_msg_info (cpl_func, "QC_BADPIX_FT = %d", count_bp_dark);

        /* Set the badpixel map of FT in output data */
		gravi_data_add_table (bad_map, NULL, GRAVI_IMAGING_DATA_FT_EXT, bad_table);
    
      FREE (cpl_vector_delete, std_vector);
	} /* End FT case */
    
    /* 
     * Compute bad pixels of SC
     */
    
    if (!gravi_data_has_extension (dark_map, GRAVI_IMAGING_DATA_SC_EXT)) {
        cpl_msg_warning (cpl_func,"The DARK map has no IMAGING_DATA_SC");
    }
    else
    {
	    cpl_msg_info (cpl_func,"Compute BADPIXEL of SC");
        
        /* Copy necessary tables */
        gravi_data_copy_ext (bad_map, dark_map, GRAVI_IMAGING_DETECTOR_SC_EXT);
	  
	    /* This is the SC */
		cpl_image * dark_img = gravi_data_get_img (dark_map, GRAVI_IMAGING_DATA_SC_EXT);
        cpl_image * std_img  = gravi_data_get_img (dark_map, GRAVI_IMAGING_ERR_SC_EXT);
		cpl_size nx = cpl_image_get_size_x (dark_img);
		cpl_size ny = cpl_image_get_size_y (dark_img);
        CPLCHECK_NUL ("Cannot get the SC data");        
		
		/* Compute an image with only the high-frequency of dark */
		cpl_image * darkhf_img = cpl_image_cast (dark_img, CPL_TYPE_DOUBLE);
		cpl_mask * kernel = cpl_mask_new (9, 9);
		cpl_mask_not (kernel);
		cpl_image_filter_mask (darkhf_img, dark_img, kernel, CPL_FILTER_MEDIAN, CPL_BORDER_FILTER);
        FREE (cpl_mask_delete, kernel);
        
		cpl_image_subtract (darkhf_img, dark_img);
        cpl_image_multiply_scalar (darkhf_img, -1.0);
        CPLCHECK_NUL ("Cannot create darkhf");

		/* Get the rms factor for dark bad pixel threshold */
		int bad_dark_factor = gravi_param_get_int (params, "gravity.calib.bad-dark-threshold");

		/* Accepted range for DARK mean */
		double dark_rms  = cpl_propertylist_get_double (dark_header, QC_DARKRMS_SC);
		double dark_mean = cpl_image_get_mean (darkhf_img);
		double dark_max = dark_mean + 2*bad_dark_factor * dark_rms;
		double dark_min = dark_mean - 2*bad_dark_factor * dark_rms;

        /* Accepted range for dark STD */
		double dark_rms_std = cpl_image_get_stdev (std_img);
		double std_max = bad_dark_factor * dark_rms_std;
		double std_min = 0.05 * dark_rms_std;

		/* Cannot detect bad pix on LOW mode */
		if (nx < 60) {
            cpl_msg_warning (cpl_func, "Don't detect SC bapixels in mode LOW");
			dark_max = 50000; //7; can't find a good param for slpit and combined
			dark_min = -5000;
            std_max *= 5;
            std_min *= 5; // FIXME: weird ??
		}

        /* If we provide FLATs, we also use them. flat_img 
         * will be filled with a filtered version of the FLAT */
        cpl_image * flat_img = NULL, * flatmed_img = NULL, * flatmask_img = NULL;
        
        if (!nflat) {
            cpl_msg_info (cpl_func,"No FLATs provided, detect only DARK");
        }
        else {
            cpl_ensure (flats_data, CPL_ERROR_NULL_INPUT, NULL);
            cpl_msg_info (cpl_func,"FLATs used to detect badpix");

            /* Init co-add FLAT image */
            flat_img = cpl_image_new (nx,ny,CPL_TYPE_DOUBLE);
            cpl_image_fill_window (flat_img, 1, 1, nx, ny, 0.0);
    
            /* Co-add FLATs, remove a self-estimate bias value
             * (the normal bias-pixels maybe wrong because of defocus), 
             * assuming only half of the detector is illuminated */
            for (int f = 0; f<nflat; f++) {
                cpl_image * img = cpl_imagelist_collapse_median_create (
                                   gravi_data_get_cube (flats_data[f],
                                   GRAVI_IMAGING_DATA_SC_EXT));
                cpl_image_subtract (img, dark_img);
        
                double bias = gravi_image_get_quantile (img, 0.25);
                cpl_image_subtract_scalar (img, bias);
                cpl_image_add (flat_img, img);
                cpl_image_delete (img);
                CPLCHECK_NUL ("Cannot add flats");
            } /* End co-add flats */
    
            /* Create kernel for horizontal median filtering 
             * Note that the large cluster is ~11 pixel wide */
            cpl_size kernel_x = 5;
            if (nx > 100) kernel_x = 11;
            if (nx > 1000) kernel_x = 31;
            cpl_msg_info (cpl_func,"Kernel of (%lld,%i) pixels for median filtering", kernel_x, 1);
            cpl_mask * kernel = cpl_mask_new (kernel_x, 1);
            cpl_mask_not (kernel);
    
            /* Run the median filter */
            flatmed_img = cpl_image_duplicate (flat_img);
            cpl_image_filter_mask (flatmed_img, flat_img, kernel, CPL_FILTER_MEDIAN, CPL_BORDER_FILTER);
            FREE (cpl_mask_delete, kernel);

            /* Compute a map of illuminated pixels at 100 [adu] */
            double threshold = 100.0;
            flatmask_img = cpl_image_duplicate (flatmed_img);
            cpl_image_threshold (flatmask_img, threshold, threshold, 0, 1.0);
            
            double fraction = cpl_image_get_flux (flatmask_img) / (nx*ny) * 100.0;
            cpl_msg_info (cpl_func, "Fraction of detector illuminated = %.2f%%  (>%.1f adu)",
                          fraction, threshold);

            /* Add this image of illuminated part in BAD map */
            gravi_data_add_img (bad_map, NULL, GRAVI_IMAGING_MASK_SC_EXT,
                                flatmask_img);

        } /* End case FLAT provided */

        
		/* Compute the number of bad pixel */
		int count_bp_dark = 0;
		int count_bp_rms = 0;
		int count_bp_flat = 0;
		int count_bp = 0;

        /* Create the badpixel image */
		cpl_image * bad_img = cpl_image_new (nx, ny, CPL_TYPE_INT);

		/* Loop on pixels */
		for (cpl_size i = 0; i < nx; i++){
			for (cpl_size j = 0; j < ny; j++){
              int nv = 0, is_bad = 0, flag = 0;

			  /* Tag the bad pixel on its mean value ... */
              double dij = cpl_image_get (darkhf_img, i + 1, j + 1, &nv);
			  if ( (dij > dark_max)||
                   (dij < dark_min)) {
				  flag += BADPIX_DARK;
				  count_bp_dark ++;
				  is_bad = 1;
			  }
			  /* ... and its variance */
			  double stdij = cpl_image_get (std_img, i + 1, j + 1, &nv);
			  if ( (stdij > std_max) ||
                   (stdij < std_min) ) {
				  flag += BADPIX_RMS;
				  count_bp_rms ++;
				  is_bad = 1;
			  }
              /* ... and the flat */
              if (flat_img) {
              double flat = cpl_image_get (flat_img, i + 1, j + 1, &nv);
              double med  = cpl_image_get (flatmed_img, i + 1, j + 1, &nv);
              double mask = cpl_image_get (flatmask_img, i + 1, j + 1, &nv);
              if ( flat < 0.5 * med && mask && i>4 && i<nx-4 ) {
                  flag += BADPIX_FLAT;
                  count_bp_flat ++;
                  is_bad = 1;
              } }

              /* Set this tag */
			  count_bp += is_bad;
			  cpl_image_set (bad_img, i + 1, j + 1, flag);

			  CPLCHECK_NUL ("Cannot compute bad pixel");
			}
		} /* End loop on pixels */

        /* Set the badpixel map of SC in IMAGING_DATA_SC */
		gravi_data_add_img (bad_map, NULL, GRAVI_IMAGING_DATA_SC_EXT, bad_img);

        /* Update QC parameter of main header */
		cpl_propertylist_append_int (bad_header, QC_BADPIX_SC, count_bp);
		cpl_msg_info (cpl_func, "QC_BADPIX_SC (total) = %d (%.2f%%)",
			      count_bp, (100.0 * count_bp) / (nx*ny));
		cpl_propertylist_append_int (bad_header, QC_BADPIX_DARK_SC, count_bp_dark);
		cpl_msg_info (cpl_func, "QC_BADPIX_DARK_SC = %d", count_bp_dark);
		cpl_propertylist_append_int (bad_header, QC_BADPIX_RMS_SC, count_bp_rms);
		cpl_msg_info (cpl_func, "QC_BADPIX_RMS_SC = %d", count_bp_rms);
        cpl_propertylist_append_int (bad_header, QC_BADPIX_FLAT_SC, count_bp_flat);
        cpl_msg_info (cpl_func, "QC_BADPIX_FLAT_SC = %d", count_bp_flat);
        
		FREE (cpl_image_delete, darkhf_img);
		FREE (cpl_image_delete, flat_img);
        FREE (cpl_image_delete, flatmed_img);
	} /* End SC case*/

    /* 
     * Compute bad pixels of ACQ
     */
    
    if (!gravi_data_has_extension (dark_map, GRAVI_IMAGING_DATA_ACQ_EXT)) {
        cpl_msg_warning (cpl_func,"The DARK map has no IMAGING_DATA_ACQ");
    }
    else
    {
	    cpl_msg_info (cpl_func,"Compute BADPIXEL of ACQ");
                
	    /* This is the SC */
		cpl_image * dark_img = gravi_data_get_img (dark_map, GRAVI_IMAGING_DATA_ACQ_EXT);
        CPLCHECK_NUL ("Cannot get dark");

		/* Compute an image with only the high-frequency of dark */
		cpl_image * darkhf_img = cpl_image_cast (dark_img, CPL_TYPE_DOUBLE);
		cpl_mask * kernel = cpl_mask_new (21, 21);
		cpl_mask_not (kernel);
		cpl_image_filter_mask (darkhf_img, dark_img, kernel, CPL_FILTER_MEDIAN, CPL_BORDER_FILTER);
        FREE (cpl_mask_delete, kernel);
        
		cpl_image_subtract (darkhf_img, dark_img);
        cpl_image_multiply_scalar (darkhf_img, -1.0);
        cpl_image_abs (darkhf_img);
        CPLCHECK_NUL ("Cannot create darkhf");

        /* Get the STD of the dark_hf */
        double dark_std = cpl_image_get_median (darkhf_img);
        cpl_msg_info (cpl_func, "DARK_ACQ_STD = %e", dark_std);

        /* Detect bad pixel as >5 STD */
        double threshold = 20 * dark_std + 1e-10;
        cpl_image_threshold (darkhf_img, threshold, threshold, 0, 1);
        cpl_image * bad_img = cpl_image_cast (darkhf_img, CPL_TYPE_INT);

        /* Remove pixel in the middle of a square pattern */
        for (cpl_size x=2; x<cpl_image_get_size_x (bad_img);x++) {
            for (cpl_size y=2; y<cpl_image_get_size_y (bad_img);y++) {
                int nv = 0;
                if (cpl_image_get (bad_img,x,y-1,&nv) &&
                    cpl_image_get (bad_img,x,y+1,&nv) &&
                    cpl_image_get (bad_img,x-1,y,&nv) &&
                    cpl_image_get (bad_img,x+1,y,&nv)) {
                    cpl_image_set (bad_img,x,y, 1);
                }
            }
        }

        /* Set QC parameter */
        cpl_size count_bp = 0;
        int nv = 0;

        /* EKW 10/01/2019 Attention, I added brackets to avoid . [-Wmisleading-indentat */
        for (cpl_size x=0; x<cpl_image_get_size_x (bad_img);x++) {
            for (cpl_size y=0; y<cpl_image_get_size_y (bad_img);y++) {
                if (cpl_image_get (bad_img,x+1,y+1,&nv)) count_bp++;
            }
        }
        
		cpl_propertylist_append_int (bad_header, "ESO QC BADPIX ACQ", count_bp);
        
        /* Set the badpixel map of ACQ in IMAGING_DATA_ACQ */
		gravi_data_add_img (bad_map, NULL, GRAVI_IMAGING_DATA_ACQ_EXT, bad_img);

		FREE (cpl_image_delete, darkhf_img);
    }
    
	/* Verbose */
 	gravi_msg_function_exit(1);
	return bad_map;
}


/*----------------------------------------------------------------------------*/
/**
 * @brief Create BIASMASK for SC from raw FLATs and raw DARK
 * 
 * @param dark_map 	   The input dark map calibration
 * @param flats_data   The input raw flats (optional)
 * @param nflat        The number of input flats (shall be 4)
 * @param params       The parameter list (no parameter used)
 * 
 * @return The BIASMASK map with the un-illuminated pixels.
 * 
 * \exception CPL_ERROR_NULL_INPUT input data is missing
 * \exception CPL_ERROR_ILLEGAL_INPUT The number of input flats is not 4
 *
 * Pixel with values lower than 100 adu in the collapsed FLAT
 * are flag with 1 in the BIASMASK, while illuminated pixels
 * (>100adu) are flag with 0 in BIASMASK.
 *
 * FIXME: improve the way we define the threshold. Use an 'ouverture'
 * filtering maybe.
 */
/*----------------------------------------------------------------------------*/

gravi_data * gravi_compute_biasmask (gravi_data * dark_map,
                                     gravi_data ** flats_data,
                                     int nflat,
                                     const cpl_parameterlist * params)
{
	/* Verbose */
	gravi_msg_function_start(1);
	cpl_ensure (dark_map,   CPL_ERROR_NULL_INPUT, NULL);
	cpl_ensure (params,     CPL_ERROR_NULL_INPUT, NULL);
    cpl_ensure (nflat==4 || nflat==0, CPL_ERROR_ILLEGAL_INPUT, NULL);

	/* Construction of the bad pixels map */
    gravi_data * biasmask_map = gravi_data_new(0);
	cpl_propertylist * biasmask_header = gravi_data_get_header (biasmask_map);

    /* Dump full header into product */
	cpl_propertylist * dark_header = gravi_data_get_header (dark_map);
    cpl_propertylist_append (biasmask_header, dark_header);
    
    /* Copy necessary tables */
    gravi_data_copy_ext (biasmask_map, dark_map, GRAVI_IMAGING_DETECTOR_SC_EXT);
	
	/* This is the SC */
	cpl_image * dark_img = gravi_data_get_img (dark_map, GRAVI_IMAGING_DATA_SC_EXT);
	cpl_size nx = cpl_image_get_size_x (dark_img);
	cpl_size ny = cpl_image_get_size_y (dark_img);
    CPLCHECK_NUL ("Cannot get the SC data");
    
    /* Init co-add FLAT image */
    cpl_image * flat_img = cpl_image_new (nx,ny,CPL_TYPE_DOUBLE);
    cpl_image_fill_window (flat_img, 1, 1, nx, ny, 0.0);
    
    /* Co-add FLATs, remove a self-estimate bias value
     * (the normal bias-pixels maybe wrong because of defocus), 
     * assuming only half of the detector is illuminated */
    for (int f = 0; f<nflat; f++) {
        cpl_image * img = cpl_imagelist_collapse_median_create (
                           gravi_data_get_cube (flats_data[f],
                           GRAVI_IMAGING_DATA_SC_EXT));
        cpl_image_subtract (img, dark_img);
    
        double bias = gravi_image_get_quantile (img, 0.25);
        cpl_image_subtract_scalar (img, bias);
        cpl_image_add (flat_img, img);
        cpl_image_delete (img);
        CPLCHECK_NUL ("Cannot add flats");
    } /* End co-add flats */
    
    /* Create kernel for horizontal median filtering 
     * Note that the large cluster is ~11 pixel wide */
    cpl_size kernel_x = 5;
    if (nx > 100) kernel_x = 11;
    if (nx > 1000) kernel_x = 31;
    cpl_msg_info (cpl_func,"Kernel of (%lld,%i) pixels for median filtering", kernel_x, 1);
    cpl_mask * kernel = cpl_mask_new (kernel_x, 1);
    cpl_mask_not (kernel);
    
    /* Run the median filter */
    cpl_image * flatmed_img = cpl_image_duplicate (flat_img);
    cpl_image_filter_mask (flatmed_img, flat_img, kernel, CPL_FILTER_MEDIAN, CPL_BORDER_FILTER);
    FREE (cpl_mask_delete, kernel);

    /* Compute a map of illuminated pixels at 100 [adu] */
    double threshold = 100.0;
    cpl_image * biasmask_img = cpl_image_duplicate (flatmed_img);
    cpl_image_threshold (biasmask_img, threshold, threshold, 1.0, 0.0);
    
    double fraction = cpl_image_get_flux (biasmask_img) / (nx*ny) * 100.0;
    cpl_msg_info (cpl_func, "Fraction of detector in bias mask = %.2f%%  (<%.1f adu)",
                  fraction, threshold);

    /* Add this image of illuminated part in BAD map */
    gravi_data_add_img (biasmask_map, NULL, GRAVI_BIAS_MASK_SC_EXT,
                        biasmask_img);
    
	FREE (cpl_image_delete, flat_img);
    FREE (cpl_image_delete, flatmed_img);
    
	/* Verbose */
 	gravi_msg_function_exit(1);
	return biasmask_map;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Create piezo transfer function for Kalman Calibration & monitoring
 *
 * @param data         The input raw data
 * @param params       The parameter list (no parameter used)
 *
 * @return The table of the OPDC data with the correct QC insrted
 *
 * \exception CPL_ERROR_NULL_INPUT input data is missing
 *
 * 1- Read the piezo commands (in volts) from the OPDC table
 * 2- Read the OPD measured by the fringe tracker (in radians)
 * 3- Using an SVD inversion, get the 20 parameters of the transfer function:
 *      OPD(n)=a_i*piezo(n-1)+b_i*piezo(n-2)+c_i*piezo(n-3)+d_i*piezo(n-4)+e_i*piezo(n-5)
 * where i is the piezo number (between 1 and 4), and a,b,c,d,e the 5 values of the
 * autoregressive function of degree 5 (AR5).
 * 4- Compute and store the QC parameters: residual errors, delay, gain, etc...
 *
 */
/*----------------------------------------------------------------------------*/

gravi_data * gravi_compute_piezotf (gravi_data * data,
                                    const cpl_parameterlist * params)
{
    /* Verbose */
    gravi_msg_function_start(1);
    cpl_ensure (data,   CPL_ERROR_NULL_INPUT, NULL);
    cpl_ensure (params,     CPL_ERROR_NULL_INPUT, NULL);

    /* Construction of the piezo tf product */
    gravi_data * piezo_tf = gravi_data_new(0);
    cpl_propertylist * piezotf_header = gravi_data_get_header (piezo_tf);

    /* Dump full header into product */
    cpl_propertylist * data_header = gravi_data_get_header (data);
    cpl_propertylist_append (piezotf_header, data_header);

    /* Copy necessary tables */
    gravi_data_copy_ext (piezo_tf, data, GRAVI_OPDC_EXT);
    
    /* Get the OPDC table */
    cpl_table * opdc = gravi_data_get_table (piezo_tf, GRAVI_OPDC_EXT);

    /* Get the necessary variables */
    int nbase = 6;
    int ntel = 4;
    int nresp = 5;
    int nsmooth = 201;
    char qc_name[100];
    char name[100];
    double phase;
    int base1,base2,base3,sign1,sign2,sign3;
    int ndit_small, dit_matrix;
    double twoPi = 2.0 * CPL_MATH_PI;
    
    
    /* DO THE COMPUTATION*/
    
    //  Determine the actuator to use, and associated gain, based on ESO.DPR.TYPE
    
    char * actuator;
    const char * dpr_type = cpl_propertylist_get_string(data_header, "ESO DPR TYPE");
    if (strcmp(dpr_type, "PIEZOTF") == 0)
        actuator = cpl_strdup("PIEZO_DL_OFFSET");
    else if (strcmp(dpr_type, "VLTITF") == 0)
        actuator = cpl_strdup("VLTI_DL_OFFSET");
    else
    {
        cpl_msg_error(cpl_func, "ESO DPR TYPE = %s is not supported!!!", dpr_type);
        return NULL;
    }
    cpl_msg_info(cpl_func, "Using %s actuator (DPR.TYPE=%s)", actuator, dpr_type);
    
    //  Read data array from OPDC table
    
    cpl_size ndit   = cpl_table_get_nrow (opdc);
    ndit_small=ndit-nsmooth-(nresp+1);
    cpl_msg_info (cpl_func, "Preparing matrix inversion with NDIT = %lld",ndit);
    cpl_array** opd   = cpl_table_get_data_array (opdc, "OPD");
    cpl_array** piezo = cpl_table_get_data_array (opdc, actuator);
    CPLCHECK_NUL ("Cannot read the OPDC data");
    
    // create temporary arrays
    
    cpl_array * phase_array = cpl_array_new(nbase,  CPL_TYPE_DOUBLE);
    cpl_array * piezo_array = cpl_array_new(ntel,  CPL_TYPE_DOUBLE);
    cpl_matrix * phase_matrix = cpl_matrix_new (ndit_small*nbase,1);
    cpl_matrix * piezo_matrix = cpl_matrix_new (ndit_small*nbase, nresp*ntel);
    cpl_matrix * piezo_header_resp = cpl_matrix_new (nresp*ntel,1);
    
    // Unwrap phase of OPD
    for (cpl_size dit = 0 ; dit < ndit-1 ; dit ++)
        for (cpl_size base = 0 ; base < nbase; base ++)
        {
            phase=cpl_array_get(opd[dit+1],base, NULL)-cpl_array_get(opd[dit],base, NULL);
            phase-= twoPi * floor( phase / twoPi );
            if (phase > CPL_MATH_PI) phase -= twoPi;
            cpl_array_set(opd[dit+1],base, phase+cpl_array_get(opd[dit],base, NULL));
        }
    
    // filter OPDs, Commands and create Matrixes
    
    for (cpl_size dit = 0 ; dit < ndit-nsmooth; dit ++)
    {
        // Filter OPDs
        cpl_array_fill_window (phase_array, 0, nbase, 0.0);
        for (cpl_size smooth = 0 ; smooth < nsmooth; smooth ++)
            cpl_array_add( phase_array, opd[smooth+dit] );
        cpl_array_multiply_scalar (phase_array, -1.0 / nsmooth);
        cpl_array_add( phase_array, opd[ (int) (nsmooth/2+dit) ] );
        
        // Filter PIEZO commands
        cpl_array_fill_window_double (piezo_array, 0, ntel, 0.0);
        for (cpl_size smooth = 0 ; smooth < nsmooth; smooth ++)
            cpl_array_add( piezo_array, piezo[smooth+dit] );
        cpl_array_multiply_scalar (piezo_array, -1.0 / nsmooth);
        cpl_array_add( piezo_array, piezo[ (int) (nsmooth/2+dit) ] );
        
        // Store values into phase matrix for SVD inversion
        if (( dit-(nresp+1) >= 0 )&( dit-(nresp+1) < ndit_small))
            for (cpl_size base = 0 ; base < nbase; base ++)
            {
                int dit_matrix =dit-(nresp+1)+base*ndit_small;
                cpl_matrix_set(phase_matrix,dit_matrix,0,cpl_array_get(phase_array, base, NULL));
            }
        
        // Store values into piezo matrix for SVD inversion
        for (cpl_size tel = 0 ; tel < ntel; tel ++)
        {
            switch (tel)
            {
                default:
                case 0:
                    base1=0;
                    sign1=-1;
                    base2=1;
                    sign2=-1;
                    base3=2;
                    sign3=-1;
                    break;
                case 1:
                    base1=0;
                    sign1=1;
                    base2=3;
                    sign2=-1;
                    base3=4;
                    sign3=-1;
                    break;
                case 2:
                    base1=1;
                    sign1=1;
                    base2=3;
                    sign2=1;
                    base3=5;
                    sign3=-1;
                    break;
                case 3:
                    base1=2;
                    sign1=1;
                    base2=4;
                    sign2=1;
                    base3=5;
                    sign3=1;
                    break;
            }
            
            for (cpl_size resp = 0 ; resp < nresp; resp ++)
                if (( dit-(nresp+1)+(1+resp) >= 0 )&( dit-(nresp+1)+(1+resp) < ndit_small))
                {
                    dit_matrix =dit-(nresp+1)+(1+resp)+base1*ndit_small;
                    cpl_matrix_set(piezo_matrix,dit_matrix,resp*ntel+tel,sign1*cpl_array_get(piezo_array, tel, NULL));
                    dit_matrix =dit-(nresp+1)+(1+resp)+base2*ndit_small;
                    cpl_matrix_set(piezo_matrix,dit_matrix,resp*ntel+tel,sign2*cpl_array_get(piezo_array, tel, NULL));
                    dit_matrix =dit-(nresp+1)+(1+resp)+base3*ndit_small;
                    cpl_matrix_set(piezo_matrix,dit_matrix,resp*ntel+tel,sign3*cpl_array_get(piezo_array, tel, NULL));
                }
        }
    }
    CPLCHECK_NUL ("Cannot create matrix for SVD inversion");
    
    // resolve SVD
    cpl_msg_info (cpl_func, "Doing SVD inversion" );
    cpl_matrix * piezo_resp = cpl_matrix_solve_normal(piezo_matrix,phase_matrix); // coef_vis is 20x1
    cpl_matrix * residuals_fit = cpl_matrix_product_create(piezo_matrix,piezo_resp);
    cpl_matrix_subtract (residuals_fit,phase_matrix); //residuals_fit is (ndit-5)*nbase
    CPLCHECK_NUL ("Failed to do SVD inversion");
    
    // Scale DL response to arbitrary PIEZO response of 17.4 rad/Volt
    if (strcmp(dpr_type, "VLTITF") == 0)
    {
        for (cpl_size tel = 0 ; tel < ntel; tel ++)
        {
            double gain = 0.0;
            for (cpl_size resp = 0 ; resp < nresp; resp ++)
                gain += cpl_matrix_get (piezo_resp, resp*ntel+tel, 0);
            gain /= 17.4;
            for (cpl_size resp = 0 ; resp < nresp; resp ++)
                cpl_matrix_set (piezo_resp, resp*ntel+tel, 0, cpl_matrix_get (piezo_resp, resp*ntel+tel, 0) / gain );
        }
    }
    
    // Get piezo response from header
    
    for (cpl_size tel = 0 ; tel < ntel; tel ++)
        for (cpl_size resp = 0 ; resp < nresp; resp ++)
        {
            sprintf (name, "ESO FT KAL P%lld_RESP%lld", tel+1, resp+1);
            if (cpl_propertylist_has (piezotf_header, name))
                cpl_matrix_set( piezo_header_resp, resp*ntel+ tel, 0, cpl_propertylist_get_double (piezotf_header, name));
            else
                cpl_matrix_set( piezo_header_resp, resp*ntel+ tel, 0, 0.0);
        }
    
    // output QC parameters
    
    sprintf (qc_name, "ESO QC FT KAL P_FIT");
    cpl_propertylist_update_double (piezotf_header, qc_name, cpl_matrix_get_stdev( residuals_fit ) );
    cpl_propertylist_set_comment (piezotf_header, qc_name, "Fitting standard deviation [rad]");
    cpl_msg_info (cpl_func, "Fit standard deviation = %e [rad]", cpl_matrix_get_stdev( residuals_fit ) );
    
    
    sprintf (name, "ESO FT RATE");
    double sampling = cpl_propertylist_get_double (piezotf_header, name);
    
    for (cpl_size tel = 0 ; tel < ntel; tel ++)
    {
        for (cpl_size resp = 0 ; resp < nresp; resp ++)
            {
            sprintf (qc_name, "ESO QC FT KAL P%lld_RESP%lld", tel+1, resp+1);
            cpl_propertylist_update_double (piezotf_header, qc_name, cpl_matrix_get( piezo_resp, resp*ntel+ tel,0 ) );
            cpl_propertylist_set_comment (piezotf_header, qc_name, "Kalman piezo response");
            
            cpl_msg_info (cpl_func, "QC FT KAL P%lld_RESP%lld = %5.5g [rad/Volts]", tel+1, resp+1, cpl_matrix_get( piezo_resp, resp*ntel+ tel,0 ));
            
            }
        
        // Get gain in rad/Volts
        
        double QC_gain=0.0;
        for (cpl_size resp = 0 ; resp < nresp; resp ++)
            QC_gain+=cpl_matrix_get( piezo_resp, resp*ntel + tel,0 );
        
        sprintf (qc_name, "ESO QC FT KAL P%lld_GAIN", tel+1);
        cpl_propertylist_update_double (piezotf_header, qc_name, QC_gain );
        cpl_propertylist_set_comment (piezotf_header, qc_name, "Open loop gain [rad/Volts]");
        
        // Get pur delay in milliseconds

        double QC_delay=0.0;
        for (cpl_size resp = 0 ; resp < nresp; resp ++)
            QC_delay+=(resp+1)*cpl_matrix_get( piezo_resp, resp*ntel + tel,0 )*sampling;
        QC_delay/=QC_gain;
        
        sprintf (qc_name, "ESO QC FT KAL P%lld_DELAY", tel+1);
        cpl_propertylist_update_double (piezotf_header, qc_name, QC_delay );
        cpl_propertylist_set_comment (piezotf_header, qc_name, "Open loop latency [ms]");
        
        // Get standard deviation
        double QC_std=0.0;
        for (cpl_size resp = 0 ; resp < nresp; resp ++)
            QC_std+=(cpl_matrix_get( piezo_resp, resp*ntel + tel,0 ) - cpl_matrix_get( piezo_header_resp, resp*ntel+ tel, 0))
            * (cpl_matrix_get( piezo_resp, resp*ntel + tel,0 ) - cpl_matrix_get( piezo_header_resp, resp*ntel+ tel, 0));
        
        sprintf (qc_name, "ESO QC FT KAL P%lld_STDEV", tel+1);
        cpl_propertylist_update_double (piezotf_header, qc_name, sqrt(QC_std)/sqrt(nresp) );
        cpl_propertylist_set_comment (piezotf_header, qc_name, "Stdev of RTC [radians]");
        
    }
    CPLCHECK_NUL ("Failed to generate and store QC parameters");

    // delete disposable arrays
    cpl_array_delete(phase_array);
    cpl_array_delete(piezo_array);
    cpl_matrix_delete(phase_matrix);
    cpl_matrix_delete(piezo_matrix);
    cpl_matrix_delete(piezo_header_resp);
    cpl_matrix_delete(piezo_resp);
    cpl_matrix_delete(residuals_fit);
    cpl_free(actuator);

    /* Verbose */
    gravi_msg_function_exit(1);
    return piezo_tf;
}

/*----------------------------------------------------------------------------*/


/**@}*/
