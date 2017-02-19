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

cpl_imagelist * gravi_acqcam_convert (cpl_imagelist * input_imglist, cpl_propertylist * header);

int gravi_acqcam_isblink (cpl_imagelist * imglist, cpl_size pos);


/*-----------------------------------------------------------------------------
                             Functions code
 -----------------------------------------------------------------------------*/

cpl_error_code gravi_preproc_acqcam (gravi_data *output_data,
                                     gravi_data *input_data,
                                     gravi_data *bad_map)
{
    gravi_msg_function_start(1);
    cpl_ensure_code (output_data, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code (input_data,  CPL_ERROR_NULL_INPUT);

    /* Check if extension exist */
    if (!gravi_data_has_extension (input_data, GRAVI_IMAGING_DATA_ACQ_EXT)) {
        cpl_msg_warning (cpl_func, "Cannot preproc the ACQCAM, not data");
        return CPL_ERROR_NONE;
    }

    cpl_propertylist * data_header;
    data_header = gravi_data_get_header (input_data);
    CPLCHECK_MSG ("Cannot get data or header");

    /* Construct a mask of badpixels */
    cpl_image * badpix_img = gravi_data_get_img (bad_map, GRAVI_IMAGING_DATA_ACQ_EXT);
    cpl_mask * badpix_mask = cpl_mask_threshold_image_create (badpix_img, 0.5, 10000);
    CPLCHECK_MSG ("Cannot get BAD map for ACQ");

    /* Get the imagelist */
    cpl_imagelist * imglist;
    imglist = gravi_data_get_cube (input_data, GRAVI_IMAGING_DATA_ACQ_EXT);
    CPLCHECK_MSG ("Cannot get image for ACQ");
    
    /* 
     * Loop on images to cleanup-badpixels 
     */
    cpl_size nrow = cpl_imagelist_get_size (imglist);
    for (cpl_size row = 0; row < nrow; row++) {

        /* Get image */
        cpl_image * img = cpl_imagelist_get (imglist, row);

        /* Cleanup-badpixel */
        cpl_image_reject_from_mask (img, badpix_mask);
        cpl_detector_interpolate_rejected (img);
        CPLCHECK_MSG ("Cannot clean-up badpixel");
    }
    
    FREE (cpl_mask_delete, badpix_mask);

    /* 
     * Convert to old format
     */

    /* Allocate new memory */
    imglist = gravi_acqcam_convert (imglist, data_header);
    CPLCHECK_MSG ("Cannot convert ACQ");
        
    /* Get the size */
    cpl_image * img = cpl_imagelist_get (imglist, 0);
    cpl_size nx = cpl_image_get_size_x (img);
    cpl_size ny = cpl_image_get_size_y (img);

    /* 
     * Remove the pupil background by the mean of blinking
     */

    if (nrow == 1) {
        gravi_msg_warning ("FIXME","Cannot remove blinked pupil (no blink)");
        
    } else {
    
        int blink = gravi_acqcam_isblink (imglist, 0) == 1 ? 0 : 1;
        cpl_size nblink = 0;
        
        /* Coadd the blinked image */
        cpl_image * blink_img = cpl_image_new (nx,ny, CPL_TYPE_DOUBLE);
        for (cpl_size row = blink; row < nrow; row +=2) {
            cpl_image_add (blink_img, cpl_imagelist_get (imglist, row));
            nblink ++;
        }
        cpl_image_divide_scalar (blink_img, nblink);
        cpl_image_fill_window (blink_img, 1, 1, nx, 1200, 0.0);
        
        gravi_msg_warning ("FIXME","Coadd the pupil blink so far.");
        gravi_msg_warning ("FIXME","Can use a vertical median to remove background");
        
        /* Remove the blink only for pupil */
        for (cpl_size row = 0; row < nrow; row ++) {
            cpl_image_subtract (cpl_imagelist_get (imglist, row), blink_img);
            CPLCHECK_MSG ("Cannot remove blinked pupil");
        }
    
        FREE (cpl_image_delete, blink_img);
    }

    /* Set in output */
    gravi_data_add_cube (output_data, NULL, GRAVI_IMAGING_DATA_ACQ_EXT, imglist);
    
    gravi_msg_function_exit(1);
    return CPL_ERROR_NONE;   
}

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

cpl_error_code gravi_reduce_acqcam (gravi_data * output_data,
                                    gravi_data * input_data)
{
    gravi_msg_function_start(1);
    cpl_ensure_code (output_data, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code (input_data,  CPL_ERROR_NULL_INPUT);

    /* Check if extension exist */
    if (!gravi_data_has_extension (input_data, GRAVI_IMAGING_DATA_ACQ_EXT)) {
        cpl_msg_warning (cpl_func, "Cannot reduce the ACQCAM, not data");
        return CPL_ERROR_NONE;
    }

    /* Get the header */
    cpl_propertylist * data_header;
    data_header = gravi_data_get_header (input_data);
    CPLCHECK_MSG ("Cannot get data or header");

    /* Get the data and header */
    cpl_imagelist * acqcam_imglist;
    acqcam_imglist = gravi_data_get_cube (input_data, GRAVI_IMAGING_DATA_ACQ_EXT);

    /* Create the output table */
    cpl_table * acqcam_table;
    acqcam_table = gravi_acqcam_table_create (acqcam_imglist, data_header);
	CPLCHECK_MSG ("Cannot compute acqcam_table");

    /* Add this output table in the gravi_data */
	gravi_data_add_table (output_data, NULL, GRAVI_OI_VIS_ACQ_EXT, acqcam_table);
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

    /* 
     * Compute TIME column
     */

    /* Time column shall contain the time from PCR.ACQ.START in [us] */
    cpl_table_new_column (output_table, "TIME", CPL_TYPE_INT);
    cpl_table_set_column_unit (output_table, "TIME", "us");

    /* Loop on DIT in cube */
    for (cpl_size row = 0; row < nrow; row++) {
        
        /* Fill the TIME column (same value for all beams) */
        double time = gravi_pfits_get_time_acqcam (header, row);
        for (int tel = 0; tel < ntel; tel ++)
            cpl_table_set (output_table, "TIME", row*ntel+tel, time);
    }
    

//    /* 
//     * Compute PUPIL column position
//     */
//
//    /* Read the references for pupil sensor */
//    double Spot5RefData[32];
//    for (int i=0; i< 4; i++ ) {
//        for (int j=0; j<4 ; j++ ) {
//            Spot5RefData[j+4*i] = gravi_pfits_get_ptfc_acqcam (header, i+4*j + 1);
//            Spot5RefData[j+4*i+16] = gravi_pfits_get_ptfc_acqcam (header, i+4*j+16 + 1);
//        }
//    }
//    
//    /* Define parameters for pupil sensor */
//    double FIfwhm_minthreshold = 1.99; /* pixels */
//    double PTFitStateErrorThreshold = 10;
//
//    /* Allocate memory for pupil sensor outputs */
//    double PupilPosition[12];      // 4 tel x XYZ
//    double PupilPositionError[12]; // 4 tel x XYZ
//    double barycentre[32];         // 4 tel x 4 sub x XY
//    double barycentreError[32];    // 4 tel x 4 sub x XY
//    double PTRefPosition[128];     // 4 tel x 16 spots x XY
//    double PTRefPositionError[128];// 4 tel x 16 spots x XY
//    double SpotsFlux[64];          // 4 tel x 16 spots
//    int    PTErrorStatus[4];       // 4 tel
//
//    /* Pupil positions (or we use array of 3)  */
//    cpl_table_new_column (output_table, "PUPIL_X", CPL_TYPE_DOUBLE);
//    cpl_table_new_column (output_table, "PUPIL_Y", CPL_TYPE_DOUBLE);
//    cpl_table_new_column (output_table, "PUPIL_Z", CPL_TYPE_DOUBLE);
//    
//    /* Loop on DIT in cube */
//    for (cpl_size row = 0; row < nrow; row++) {
//
//        /* Call the PUPIL algorithm */
//        cpl_image * img = cpl_imagelist_get (acqcam_imglist,row);
//        gvacqPupilTracker (img, Spot5RefData,
//                           FIfwhm_minthreshold, PTFitStateErrorThreshold,
//                           PupilPosition, PupilPositionError,
//                           barycentre, barycentreError,
//                           PTRefPosition, PTRefPositionError,
//                           SpotsFlux, PTErrorStatus);
//
//        /* Fill the columns */
//        for (int tel=0; tel<4; tel++) {
//            cpl_table_set (output_table, "PUPIL_X", row*ntel+tel, PupilPosition[tel*3+0]);
//            cpl_table_set (output_table, "PUPIL_Y", row*ntel+tel, PupilPosition[tel*3+1]);
//            cpl_table_set (output_table, "PUPIL_Z", row*ntel+tel, PupilPosition[tel*3+2]);
//        }
//
//    } /* End loop on DIT in cube */
//
    
    /* 
     * Other column
     */
    
    /* Field positions (or we use array of 2) */
    // cpl_table_new_column (output_table, "FIELD_SC_X", CPL_TYPE_DOUBLE);
    // cpl_table_new_column (output_table, "FIELD_SC_Y", CPL_TYPE_DOUBLE);
    // cpl_table_new_column (output_table, "FIELD_FT_X", CPL_TYPE_DOUBLE);
    // cpl_table_new_column (output_table, "FIELD_FT_Y", CPL_TYPE_DOUBLE);

    /* Aberation as arrays of 27 coefficients */
    // cpl_table_new_column_array (output_table, "ABERATION", CPL_TYPE_DOUBLE, nzernick);
    

    gravi_msg_function_exit(1);
    return output_table;
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
        output_imglist = cpl_imagelist_new ();
        
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
            CPLCHECK_NUL ("Error when convert ACQ");
        } /* End loop on images */
        
    } /* End case image is split in 16 */

    gravi_msg_function_exit(1);
    return output_imglist;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Return 1 if flux(pos) < flux(pos+1)
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
    cpl_size llx =  1;
    cpl_size urx =  2048;
    cpl_size lly =  1200;
    cpl_size ury =  1536;

    /* tell if the first image is blink on or off */
    double flux0 = cpl_image_get_flux_window (cpl_imagelist_get (imglist,pos), llx, lly, urx, ury);
    double flux1 = cpl_image_get_flux_window (cpl_imagelist_get (imglist,pos+1), llx, lly, urx, ury);
    if (flux0 > flux1) blinkOFF = 0;

    gravi_msg_function_exit(1);
    return blinkOFF;
}



/**@}*/

