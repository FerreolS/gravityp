/* $Id: gravi_vis.c,v 1.10 2014/11/12 15:10:40 nazouaoui Exp $
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
 * @defgroup gravity_acqcam  TBD
 */
/**@{*/

/*-----------------------------------------------------------------------------
                                   Includes
 -----------------------------------------------------------------------------*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#define _XOPEN_SOURCE
#include <cpl.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "gravi_dfs.h"
#include "gravi_data.h"
#include "gravi_pfits.h"

#include "gravi_cpl.h"
#include "gravi_utils.h"

#include "gravi_acqcam.h"

/*-----------------------------------------------------------------------------
                               Private prototypes
 -----------------------------------------------------------------------------*/

cpl_table * gravi_acqcam_table_create (cpl_imagelist * acqcam_imglist,
                                       cpl_propertylist * header);

cpl_error_code gravi_acqcam_set_calibration (cpl_frameset * frameset);
cpl_error_code gravi_acqcam_set_parameters (cpl_propertylist * header);
cpl_imagelist * gravi_acqcam_convert (cpl_imagelist * input_imglist, cpl_propertylist * header);
int gravi_acqcam_isblink (cpl_imagelist * imglist, cpl_size pos);


/*-----------------------------------------------------------------------------
                             Functions code
 -----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
 * @brief Reduce the ACQ camera images
 *  
 * @param data        The gravi_data input/output
 * 
 * The routine read the images from the ACQ camera, reduce these images
 * with the XXXXX algorithm. The resulting tables is append into
 * data.
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_reduce_acqcam (gravi_data * data, cpl_frameset * frameset)
{
    gravi_msg_function_start(1);
    cpl_ensure_code (data,     CPL_ERROR_NULL_INPUT);
    cpl_ensure_code (frameset, CPL_ERROR_NULL_INPUT);

    /* Check if extension exist */
    if (!gravi_data_has_extension (data, GRAVI_IMAGING_DATA_ACQ_EXT)) {
        cpl_msg_warning (cpl_func, "Cannot reduce the ACQCAM, not data");
        return CPL_ERROR_NONE;
    }


    cpl_propertylist * data_header;
    data_header = gravi_data_get_header (data);
    CPLCHECK_MSG ("Cannot get data or header");

    /* Load calibrations from input frameset */
    gravi_acqcam_set_calibration (frameset);

    /* Load calibrations from input frameset */
    gravi_acqcam_set_parameters (data_header);

    /* Get the data and header. (descramble the 16 sub-windows)  */
    cpl_imagelist * acqcam_imglist;
    acqcam_imglist = gravi_data_get_cube (data, GRAVI_IMAGING_DATA_ACQ_EXT);
    acqcam_imglist = gravi_acqcam_convert (acqcam_imglist, data_header);

    /* Create the output table */
    cpl_table * acqcam_table;
    acqcam_table = gravi_acqcam_table_create (acqcam_imglist, data_header);
	CPLCHECK_MSG ("Cannot compute acqcam_table");

    /* Delete data */
    FREE (cpl_imagelist_delete, acqcam_imglist);

    /* Add this output table in the gravi_data */
	gravi_data_add_table (data, NULL, GRAVI_OI_VIS_ACQ_EXT, acqcam_table);
	CPLCHECK_MSG ("Cannot add acqcam_table in data");
    
    gravi_msg_function_exit(1);
    return CPL_ERROR_NONE;   
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Reduce the ACQ camera images
 *  
 * @param acqcam_imglist     The imagelist with the ACQ images 
 * @param header             The header corresponding to these images
 * 
 * The routine creates a cpl_table with NDIT x 4 rows, ordered
 * T1,T2,T3,T4,T1,T2... and with the following columns:
 * 
 * TIME [us]
 * FIELD_X, FIELD_Y  [???]
 * PUPIL_X, PUPIL_Y, PUPIL_Z  [???]
 * ABERATION, array of 27 zernick coefficients [???]
 */
/*----------------------------------------------------------------------------*/

cpl_table * gravi_acqcam_table_create (cpl_imagelist * acqcam_imglist,
                                       cpl_propertylist * header)
{
    gravi_msg_function_start(1);
    
    cpl_ensure (acqcam_imglist, CPL_ERROR_NULL_INPUT, NULL);
    cpl_ensure (header,         CPL_ERROR_NULL_INPUT, NULL);

    cpl_size ntel = 4, nzernick = 27;
    cpl_size nrow = cpl_imagelist_get_size (acqcam_imglist);

    /* Here we shall build the table */
    cpl_table * output_table;
    output_table = cpl_table_new (nrow * ntel);

    /* Time column shall contain the time from PCR.ACQ.START in [us] */
    cpl_table_new_column (output_table, "TIME", CPL_TYPE_INT);
    cpl_table_set_column_unit (output_table, "TIME", "us");

    /* Field positions (or we use array of 2) */
    cpl_table_new_column (output_table, "FIELD_SC_X", CPL_TYPE_DOUBLE);
    cpl_table_new_column (output_table, "FIELD_SC_Y", CPL_TYPE_DOUBLE);
    cpl_table_new_column (output_table, "FIELD_FT_X", CPL_TYPE_DOUBLE);
    cpl_table_new_column (output_table, "FIELD_FT_Y", CPL_TYPE_DOUBLE);

    /* Pupil positions (or we use array of 3)  */
    cpl_table_new_column (output_table, "PUPIL_X", CPL_TYPE_DOUBLE);
    cpl_table_new_column (output_table, "PUPIL_Y", CPL_TYPE_DOUBLE);
    cpl_table_new_column (output_table, "PUPIL_Z", CPL_TYPE_DOUBLE);

    /* Aberation as arrays of 27 coefficients */
    cpl_table_new_column_array (output_table, "ABERATION", CPL_TYPE_DOUBLE, nzernick);

    /* Get if the first is blinked */
    int blinkOFF = gravi_acqcam_isblink (acqcam_imglist, 0);

    /* Loop on DIT in cube */
    for (cpl_size row = 0; row < nrow; row++) {

        /* Fill the TIME column (same value for all beams) */
        double time = gravi_pfits_get_time_acqcam (header, row);
        for (int tel = 0; tel < ntel; tel ++)
            cpl_table_set (output_table, "TIME", row*ntel+tel, time);

        /* Compute the FIELD, PUPIL and ABERR for this image */
        cpl_image * img = cpl_imagelist_get (acqcam_imglist,row);
        // gvoacqProcessImageAC (img, blinkOFF, ACMODE_FIELD, &allDataVars);
        // gvoacqProcessImageAC (img, blinkOFF, ACMODE_PUPIL, &allDataVars);
        // gvoacqProcessImageAC (img, blinkOFF, ACMODE_ABERR, &allDataVars);

        /* Use the blinking to update the pupil sky */
        // storeSky (img, STORSKY_PUPIL);
  
        /* Fill the columns */
        // int tel=0;
        // cpl_table_set (output_table, "FIELD_SC_X", row*ntel+tel, (allDataVars.acFiPos)->value[0].arm1);
        // cpl_table_set (output_table, "FIELD_SC_Y", row*ntel+tel, (allDataVars.acFiPos)->value[1].arm1);
        // cpl_table_set (output_table, "FIELD_FT_X", row*ntel+tel, (allDataVars.acFiPos2)->value[0].arm1);
        // cpl_table_set (output_table, "FIELD_FT_Y", row*ntel+tel, (allDataVars.acFiPos2)->value[1].arm1);
        // cpl_table_set (output_table, "PUPIL_X", row*ntel+tel, (allDataVars.acPtPos)->value[0].arm1);
        // cpl_table_set (output_table, "PUPIL_Y", row*ntel+tel, (allDataVars.acPtPos)->value[1].arm1);
        // cpl_table_set (output_table, "PUPIL_Z", row*ntel+tel, (allDataVars.acPtPos)->value[2].arm1);
        // 
        // tel=1;
        // cpl_table_set (output_table, "FIELD_SC_X", row*ntel+tel, (allDataVars.acFiPos)->value[0].arm2);
        // cpl_table_set (output_table, "FIELD_SC_Y", row*ntel+tel, (allDataVars.acFiPos)->value[1].arm2);
        // cpl_table_set (output_table, "FIELD_FT_X", row*ntel+tel, (allDataVars.acFiPos2)->value[0].arm2);
        // cpl_table_set (output_table, "FIELD_FT_Y", row*ntel+tel, (allDataVars.acFiPos2)->value[1].arm2);
        // cpl_table_set (output_table, "PUPIL_X", row*ntel+tel, (allDataVars.acPtPos)->value[0].arm2);
        // cpl_table_set (output_table, "PUPIL_Y", row*ntel+tel, (allDataVars.acPtPos)->value[1].arm2);
        // cpl_table_set (output_table, "PUPIL_Z", row*ntel+tel, (allDataVars.acPtPos)->value[2].arm2);
        // 
        // tel=2;
        // cpl_table_set (output_table, "FIELD_SC_X", row*ntel+tel, (allDataVars.acFiPos)->value[0].arm3);
        // cpl_table_set (output_table, "FIELD_SC_Y", row*ntel+tel, (allDataVars.acFiPos)->value[1].arm3);
        // cpl_table_set (output_table, "FIELD_FT_X", row*ntel+tel, (allDataVars.acFiPos2)->value[0].arm3);
        // cpl_table_set (output_table, "FIELD_FT_Y", row*ntel+tel, (allDataVars.acFiPos2)->value[1].arm3);
        // cpl_table_set (output_table, "PUPIL_X", row*ntel+tel, (allDataVars.acPtPos)->value[0].arm3);
        // cpl_table_set (output_table, "PUPIL_Y", row*ntel+tel, (allDataVars.acPtPos)->value[1].arm3);
        // cpl_table_set (output_table, "PUPIL_Z", row*ntel+tel, (allDataVars.acPtPos)->value[2].arm3);
        // 
        // tel=3;
        // cpl_table_set (output_table, "FIELD_SC_X", row*ntel+tel, (allDataVars.acFiPos)->value[0].arm4);
        // cpl_table_set (output_table, "FIELD_SC_Y", row*ntel+tel, (allDataVars.acFiPos)->value[1].arm4);
        // cpl_table_set (output_table, "FIELD_FT_X", row*ntel+tel, (allDataVars.acFiPos2)->value[0].arm4);
        // cpl_table_set (output_table, "FIELD_FT_Y", row*ntel+tel, (allDataVars.acFiPos2)->value[1].arm4);
        // cpl_table_set (output_table, "PUPIL_X", row*ntel+tel, (allDataVars.acPtPos)->value[0].arm4);
        // cpl_table_set (output_table, "PUPIL_Y", row*ntel+tel, (allDataVars.acPtPos)->value[1].arm4);
        // cpl_table_set (output_table, "PUPIL_Z", row*ntel+tel, (allDataVars.acPtPos)->value[2].arm4);

        /* Blink */
        if (blinkOFF) blinkOFF = 0;
        else blinkOFF = 1;
        
    } /* End loop on DIT in cube */

    gravi_msg_function_exit(1);
    return output_table;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Load the calibration for the ACQ_CAM camera
 * 
 * @param frameset : shall include FLAT_AQC and BAD_ACQ
 * 
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_acqcam_set_calibration (cpl_frameset * frameset)
{
    gravi_msg_function_start(1);
    cpl_ensure_code (frameset, CPL_ERROR_NULL_INPUT);
    
    cpl_frame * frame = NULL;

    /* FLAT_ACQ */
    frame = cpl_frameset_find (frameset, GRAVI_FLAT_ACQ_MAP);
    cpl_ensure_code (frame, CPL_ERROR_ILLEGAL_INPUT);

	cpl_msg_info (cpl_func, "FLAT_ACQ file as  %s\n",
                  cpl_frame_get_filename (frame));

    /* Set in global variable */
    // MasterFlat = cpl_image_load (cpl_frame_get_filename (frame), CPL_TYPE_DOUBLE, 0, 0);
    // cpl_ensure_code (MasterFlat, CPL_ERROR_ILLEGAL_INPUT);
    
    
    /* BAD_ACQ */
    frame = cpl_frameset_find (frameset, GRAVI_BAD_ACQ_MAP);
    cpl_ensure_code (frame, CPL_ERROR_ILLEGAL_INPUT);

	cpl_msg_info (cpl_func, "BAD_ACQ file as  %s\n",
                  cpl_frame_get_filename (frame));
    
    /* Set in global variable */
    cpl_image * bad_img = cpl_image_load (cpl_frame_get_filename (frame), CPL_TYPE_DOUBLE, 0, 0);
    // DeadPixelMap = cpl_mask_threshold_image_create (bad_img, 0, 2);
    // cpl_mask_not (DeadPixelMap);
    // cpl_ensure_code (DeadPixelMap, CPL_ERROR_ILLEGAL_INPUT);
    FREE (cpl_image_delete, bad_img);

    // ...

    gravi_msg_function_exit(1);
    return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Load the parameters for the ACQ_CAM reduction
 * 
 * @param header:
 * 
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_acqcam_set_parameters (cpl_propertylist * header)
{
    gravi_msg_function_start(1);
    
    double FIfieldWindowSizeAndSigma[2] = { 110, 3.0 };
    // setFIfitWindowSizeAndSigma(FIfieldWindowSizeAndSigma);

    int coordinates[8] = { 287, 803, 1257, 1699, 228, 222, 204, 214 };
    // setCoordinates(coordinates);

    double FIminFWHM     = 1.5;
    // setFIminFWHM(FIminFWHM);

    int    FIwindowGauss = 16;
    // setFIwindowGauss(FIwindowGauss);

    int    ABSLensletSizeIN       = 10;
    // setABSLensletSize(ABSLensletSizeIN);

    int    RefWindowPosPT[8]      = { 292, 759, 1230, 1698, 1376, 1371, 1375, 1383 };
    // setRefWindowPosPT(RefWindowPosPT);

    int    PTImageWindowSize      = 220; /* pupil tracker image window */
    // setPTImageWindowSize(PTImageWindowSize);
    

    // ...

    gravi_msg_function_exit(1);
    return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Descramble, if necessary, the 16 sub-windows into a larger image
 *        Data are duplicated.
 * @param input_img
 * 
 */
/*----------------------------------------------------------------------------*/
cpl_imagelist * gravi_acqcam_convert (cpl_imagelist * input_imglist, cpl_propertylist * header)
{
    gravi_msg_function_start(1);
    cpl_ensure (input_imglist, CPL_ERROR_NULL_INPUT, NULL);
    int nv = 0;

    cpl_imagelist * output_imglist = NULL;
    
    cpl_image * output_img = NULL;
    cpl_image * input_img = cpl_imagelist_get (input_imglist, 0);
    
    cpl_size nx = cpl_image_get_size_x (input_img);
    cpl_size ny = cpl_image_get_size_y (input_img);

    cpl_size Nx = 2048;
    cpl_size Ny = 1536;


    if (nx == Nx && ny == Ny) {
        /* Create Image */
        output_imglist = cpl_imagelist_duplicate (input_imglist);
        
    } else {
        
        /* Get size of each sub-window */
        cpl_size sizex = cpl_propertylist_get_int (header, "ESO DET1 FRAMES NX");
        cpl_size sizey = cpl_propertylist_get_int (header, "ESO DET1 FRAMES NY");
        cpl_ensure (sizex > 1, CPL_ERROR_ILLEGAL_INPUT, NULL);
        cpl_ensure (sizey > 1, CPL_ERROR_ILLEGAL_INPUT, NULL);
        
        /* Loop on images */
        cpl_size nimg = cpl_imagelist_get_size (input_imglist);
        for (cpl_size img = 0; img < nimg; img++) {

            /* Create Image */
            input_img  = cpl_imagelist_get (input_imglist, img);
            output_img = cpl_image_new (Nx,Ny,cpl_image_get_type (input_img));
            
            /* The image is cut into 16 sub-frames */
            for (int frame = 0; frame < 16; frame++) {
                cpl_size framex = frame-(frame/4)*4;
                cpl_size framey = frame/4;
                
                /* Read the start of this frame */
                char name[90];
                sprintf (name,"ESO DET1 FRAM%d STRX",frame+1);
                cpl_size strx = cpl_propertylist_get_int (header, name);
                cpl_ensure (strx > 1, CPL_ERROR_ILLEGAL_INPUT, NULL);
                sprintf (name,"ESO DET1 FRAM%d STRY",frame+1);
                cpl_size stry = cpl_propertylist_get_int (header, name);
                cpl_ensure (stry > 1, CPL_ERROR_ILLEGAL_INPUT, NULL);
           
                for (cpl_size i=1; i<=sizex; i++)
                    for (cpl_size j=1; j<=sizey; j++) {
                        double content = cpl_image_get (input_img, i+framex*sizex, j+framey*sizey, &nv);
                        cpl_image_set (output_img, i+strx, j+stry, content);
                    }
           
            } /* End loop on sub-windows */

            cpl_imagelist_set (output_imglist, output_img, img);
        } /* End loop on images */
        
    } /* End case image is split in 16 */

    gravi_msg_function_exit(1);
    return output_imglist;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Return 0 if flux(pos) > flux(pos+1)
 * 
 */
/*----------------------------------------------------------------------------*/
int gravi_acqcam_isblink (cpl_imagelist * imglist, cpl_size pos)
{
    gravi_msg_function_start(1);
    cpl_ensure (imglist, CPL_ERROR_NULL_INPUT, -1);
    
    int blinkOFF = 1;

    /* Ensure two images */
    cpl_size nrow = cpl_imagelist_get_size (imglist);
    if (nrow == 1) return blinkOFF;

    /* Part of the image to analyse */
    cpl_size llx =  0; //RefWindowPosPT[0] - PTImageWindowSize / 2;
    cpl_size urx =  0; //RefWindowPosPT[3] + PTImageWindowSize / 2;
    cpl_size lly =  0; //RefWindowPosPT[4] - PTImageWindowSize / 2;
    cpl_size ury =  0; //RefWindowPosPT[7] + PTImageWindowSize / 2;

    /* tell if the first image is blink on or off */
    double flux0 = cpl_image_get_flux_window (cpl_imagelist_get (imglist,pos), llx, lly, urx, ury);
    double flux1 = cpl_image_get_flux_window (cpl_imagelist_get (imglist,pos+1), llx, lly, urx, ury);
    if (flux0 > flux1) blinkOFF = 0;

    gravi_msg_function_exit(1);
    return blinkOFF;
}



/**@}*/
