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

cpl_error_code gravi_reduce_acqcam (gravi_data * data)
{
    gravi_msg_function_start(1);
    cpl_ensure_code (data, CPL_ERROR_NULL_INPUT);

    /* Check if extension exist */
    if (!gravi_data_has_extension (data, GRAVI_IMAGING_DATA_ACQ_EXT)) {
        cpl_msg_warning (cpl_func, "Cannot reduce the ACQCAM, not data");
        return CPL_ERROR_NONE;
    }

    /* Get the data and header */
    cpl_imagelist * acqcam_imglist;
    acqcam_imglist = gravi_data_get_cube (data, GRAVI_IMAGING_DATA_ACQ_EXT);

    cpl_propertylist * data_header;
    data_header = gravi_data_get_header (data);
    CPLCHECK_MSG ("Cannot get data or header");

    /* Create the output table */
    cpl_table * acqcam_table;
    acqcam_table = gravi_acqcam_table_create (acqcam_imglist, data_header);
	CPLCHECK_MSG ("Cannot compute acqcam_table");

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
    cpl_table_new_column (output_table, "FIELD_X", CPL_TYPE_DOUBLE);
    cpl_table_new_column (output_table, "FIELD_Y", CPL_TYPE_DOUBLE);

    /* Pupil positions (or we use array of 3)  */
    cpl_table_new_column (output_table, "PUPIL_X", CPL_TYPE_DOUBLE);
    cpl_table_new_column (output_table, "PUPIL_Y", CPL_TYPE_DOUBLE);
    cpl_table_new_column (output_table, "PUPIL_Z", CPL_TYPE_DOUBLE);

    /* Aberation as arrays of 27 coefficients */
    cpl_table_new_column_array (output_table, "ABERATION", CPL_TYPE_DOUBLE, nzernick);


    /* Loop on DIT in cube */
    for (cpl_size row = 0; row < nrow; row++) {

        /* Fill the TIME column (same value for all beams) */
        double time = gravi_pfits_get_time_acqcam (header, row);
        for (int tel = 0; tel < ntel; tel ++)
            cpl_table_set (output_table, "TIME", row*ntel+tel, time);
        
        /* 
         * FIXME: fill the table with measurments done in imagelist !!
         */
        if (row == 0)
            gravi_msg_warning ("FIXME", "ACQ_CAM images not reduced yet!!");
        

    } /* End loop on DIT in cube */

    gravi_msg_function_exit(1);
    return output_table;
}


/**@}*/
