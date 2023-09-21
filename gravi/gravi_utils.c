/* $Id: gravi_utils.c,v 1.10 2011/05/31 06:10:40 nazouaoui Exp $
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
 * @defgroup gravi_utils     Miscellaneous Utilities
 *
 * This module contains some utilities to access to specific quantities.
 *
 */
/**@{*/

/*-----------------------------------------------------------------------------
                                   Includes
 -----------------------------------------------------------------------------*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cpl.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "gravi_data.h"

#include "gravi_dfs.h"
#include "gravi_pfits.h"
#include "gravi_cpl.h"

#include "gravi_utils.h"

/*-----------------------------------------------------------------------------
                                   Gloval variable
 -----------------------------------------------------------------------------*/

/* Baseline: [6 bases][tel1, tel2] */
int GRAVI_BASE_TEL[6][2] = { {0,1}, {0,2}, {0,3}, {1,2}, {1,3}, {2,3} };
char GRAVI_BASE_NAME[6][3] = {"12","13","14","23","24","34"};

/* Closing triangle: [6 bases][2 closing triangle][base1, base2] */
int GRAVI_TRI_BASE[6][2][2] = { {{1,3},{2,4}}, {{0,3},{2,5}}, {{0,4},{1,5}}, {{0,1},{4,5}}, {{0,2},{3,5}}, {{1,2},{3,4}}};
int GRAVI_TRI_SIGN[6][2][2] = { {{1,-1},{1,-1}},  {{1,1},{1,-1}}, {{1,1},{1,1}}, {{-1,1},{1,-1}}, {{-1,1},{1,1}}, {{-1,1},{-1,1}} };

/* Closing telescope: [3 closures][tel1, tel2, tel3] */
int GRAVI_CLO_TEL[4][3] = {{0,1,2}, {0,1,3}, {0,2,3}, {1,2,3}};
char GRAVI_CLO_NAME[4][4] = {"123", "124", "134", "234"};

/* Baseline in closure: [3 closures][base1, base2, -base3] */
int GRAVI_CLO_BASE[4][3] = {{0,3,1}, {0,4,2}, {1,5,2}, {3,5,4}};

/* DATA and DATAERR region names */
char GRAVI_DATA[50][7] =
  {"DATA1", "DATA2", "DATA3", "DATA4", "DATA5", "DATA6", "DATA7", "DATA8", "DATA9", "DATA10",
   "DATA11","DATA12","DATA13","DATA14","DATA15","DATA16","DATA17","DATA18","DATA19","DATA20",
   "DATA21","DATA22","DATA23","DATA24","DATA25","DATA26","DATA27","DATA28","DATA29","DATA30",
   "DATA31","DATA32","DATA33","DATA34","DATA35","DATA36","DATA37","DATA38","DATA39","DATA40",
   "DATA41","DATA42","DATA43","DATA44","DATA45","DATA46","DATA47","DATA48","DATA49","DATA50"};

/* DATA and DATAERR region names */
char GRAVI_DATAERR[50][10] =
  {"DATAERR1", "DATAERR2", "DATAERR3", "DATAERR4", "DATAERR5", "DATAERR6", "DATAERR7", "DATAERR8", "DATAERR9", "DATAERR10",
   "DATAERR11","DATAERR12","DATAERR13","DATAERR14","DATAERR15","DATAERR16","DATAERR17","DATAERR18","DATAERR19","DATAERR20",
   "DATAERR21","DATAERR22","DATAERR23","DATAERR24","DATAERR25","DATAERR26","DATAERR27","DATAERR28","DATAERR29","DATAERR30",
   "DATAERR31","DATAERR32","DATAERR33","DATAERR34","DATAERR35","DATAERR36","DATAERR37","DATAERR38","DATAERR39","DATAERR40",
   "DATAERR41","DATAERR42","DATAERR43","DATAERR44","DATAERR45","DATAERR46","DATAERR47","DATAERR48","DATAERR49","DATAERR50"};


/* IP of GRAVITY beams */
int GRAVI_LABINPUT[4] = {7,5,3,1};

/*-----------------------------------------------------------------------------
                                   Function code
 -----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
 * @brief    Get the pipeline copyright and license
 * @return   The copyright and license string
 * 
 * The function returns a pointer to the statically allocated license string.
 * This string should not be modified using the returned pointer.
 */
/*----------------------------------------------------------------------------*/

const char * gravi_get_license(void)
{
    const char  *   gravi_license = 
        "This file is part of the GRAVI Instrument Pipeline\n"
        "Copyright (C) 2002,2003 European Southern Observatory\n"
        "\n"
        "This program is free software; you can redistribute it and/or modify\n"
        "it under the terms of the GNU General Public License as published by\n"
        "the Free Software Foundation; either version 2 of the License, or\n"
        "(at your option) any later version.\n"
        "\n"
        "This program is distributed in the hope that it will be useful,\n"
        "but WITHOUT ANY WARRANTY; without even the implied warranty of\n"
        "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n"
        "GNU General Public License for more details.\n"
        "\n"
        "You should have received a copy of the GNU General Public License\n"
        "along with this program; if not, write to the Free Software\n"
        "Foundation, Inc., 59 Temple Place, Suite 330, Boston, \n"
        "MA  02111-1307  USA" ;
    return gravi_license ;
}

cpl_error_code gravi_msg_warning (const char * component, const char * msg)
{
    cpl_msg_warning (component,"***********************************************");
    cpl_msg_warning (component,"                                               ");
    cpl_msg_warning (component," %s ", msg);
    cpl_msg_warning (component,"                                               ");
    cpl_msg_warning (component,"***********************************************");
    return CPL_ERROR_NONE;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief    Create the oi table (oi_vis, oi_vis2, oi_t3)
 * 
 * @param  nwave 		 	The size of the wavelength
 * @param  nrow			    The number of acquisitions or of the file using
 * 	  	  	  	  	  	    to compute the visibilities
 * @param  oi_name			The name of the table to create
 * 
 * @return   The oi table needed with the number of acquisitions and size of the
 * wavelength of each spectrum
 * 
 * The function returns a cpl_table
 */
/*---------------------------------------------------------------------------*/

cpl_table * gravi_table_oi_create (int nwave, int nrow, const char * oi_name)
{
    gravi_msg_function_start(0);
	cpl_ensure (nwave>0, CPL_ERROR_ILLEGAL_INPUT, NULL);
	cpl_ensure (nrow>0,  CPL_ERROR_ILLEGAL_INPUT, NULL);
	cpl_ensure (oi_name, CPL_ERROR_NULL_INPUT, NULL);

    /* We init most of the column */
    int init = 1;

	cpl_table * oi_table = NULL;

	/* Creating the columns */

	if (!strcmp(oi_name, GRAVI_OI_VIS_EXT)){
		oi_table = cpl_table_new(nrow * 6);

		cpl_table_new_column (oi_table, "TARGET_ID", CPL_TYPE_INT);
		cpl_table_set_column_savetype(oi_table, "TARGET_ID", CPL_TYPE_SHORT);
		if (init) cpl_table_fill_column_window_int (oi_table, "TARGET_ID", 0, nrow * 6, 0);

		cpl_table_new_column (oi_table, "TIME", CPL_TYPE_DOUBLE);
		if (init) cpl_table_fill_column_window_double (oi_table, "TIME", 0, nrow * 6, 0);
		cpl_table_set_column_unit (oi_table, "TIME", "s");

		cpl_table_new_column (oi_table, "MJD", CPL_TYPE_DOUBLE);
		if (init) cpl_table_fill_column_window_double (oi_table, "MJD", 0, nrow * 6, 0.0);
		cpl_table_set_column_unit (oi_table, "MJD", "d");
		
		cpl_table_new_column (oi_table, "INT_TIME", CPL_TYPE_DOUBLE);
		if(init) cpl_table_fill_column_window_double (oi_table, "INT_TIME", 0, nrow * 6, 0.0);
		cpl_table_set_column_unit (oi_table, "INT_TIME", "s");

		cpl_table_new_column_array (oi_table, "VISDATA", CPL_TYPE_DOUBLE_COMPLEX, nwave);
		cpl_table_set_column_unit (oi_table, "VISDATA", "e");
		cpl_table_new_column_array (oi_table, "VISERR", CPL_TYPE_DOUBLE_COMPLEX, nwave);
		cpl_table_set_column_unit (oi_table, "VISERR", "e");
		
		cpl_table_new_column_array (oi_table, "VISAMP", CPL_TYPE_DOUBLE, nwave);
		cpl_table_new_column_array (oi_table, "VISAMPERR", CPL_TYPE_DOUBLE, nwave);
		
		cpl_table_new_column_array (oi_table, "VISPHI", CPL_TYPE_DOUBLE, nwave);
		cpl_table_set_column_unit (oi_table, "VISPHI", "deg");
		
		cpl_table_new_column_array (oi_table, "VISPHIERR", CPL_TYPE_DOUBLE, nwave);
		cpl_table_set_column_unit (oi_table, "VISPHIERR", "deg");

		cpl_table_new_column_array (oi_table, "RVIS", CPL_TYPE_DOUBLE, nwave);
		cpl_table_set_column_unit (oi_table, "RVIS", "e");
        
		cpl_table_new_column_array (oi_table, "RVISERR", CPL_TYPE_DOUBLE, nwave);
		cpl_table_set_column_unit (oi_table, "RVISERR", "e");

		cpl_table_new_column_array (oi_table, "IVIS", CPL_TYPE_DOUBLE, nwave);
		cpl_table_set_column_unit (oi_table, "IVIS", "e");
        
		cpl_table_new_column_array (oi_table, "IVISERR", CPL_TYPE_DOUBLE, nwave);
		cpl_table_set_column_unit (oi_table, "IVISERR", "e");
		
		cpl_table_new_column (oi_table, "UCOORD", CPL_TYPE_DOUBLE);
		if(init) cpl_table_fill_column_window_double (oi_table, "UCOORD", 0, nrow * 6, 0.0);
		cpl_table_set_column_unit (oi_table, "UCOORD", "m");
		
		cpl_table_new_column (oi_table, "VCOORD", CPL_TYPE_DOUBLE);
		if(init) cpl_table_fill_column_window_double (oi_table, "VCOORD", 0, nrow * 6, 0.0);
		cpl_table_set_column_unit (oi_table, "VCOORD", "m");
		
		cpl_table_new_column_array (oi_table, "STA_INDEX", CPL_TYPE_INT, 2);
		cpl_table_set_column_savetype(oi_table, "STA_INDEX", CPL_TYPE_SHORT);
	}
	else if (!strcmp(oi_name, GRAVI_OI_VIS2_EXT)){
		oi_table = cpl_table_new(nrow * 6);
		cpl_table_new_column (oi_table, "TARGET_ID", CPL_TYPE_INT);
		cpl_table_set_column_savetype(oi_table, "TARGET_ID", CPL_TYPE_SHORT);
		if(init) cpl_table_fill_column_window_int (oi_table, "TARGET_ID", 0, nrow * 6, 0);
		cpl_table_new_column (oi_table, "TIME", CPL_TYPE_DOUBLE);
		if(init) cpl_table_fill_column_window_double (oi_table, "TIME", 0, nrow * 6, 0);
		cpl_table_set_column_unit (oi_table, "TIME", "s");
		cpl_table_new_column (oi_table, "MJD", CPL_TYPE_DOUBLE);
		if(init) cpl_table_fill_column_window_double (oi_table, "MJD", 0, nrow * 6, 0.0);
		cpl_table_set_column_unit (oi_table, "MJD", "d");
		cpl_table_new_column (oi_table, "INT_TIME", CPL_TYPE_DOUBLE);
		if(init) cpl_table_fill_column_window_double (oi_table, "INT_TIME",
															  0, nrow * 6, 0.0);
		cpl_table_set_column_unit (oi_table, "INT_TIME", "s");
		cpl_table_new_column_array (oi_table, "VIS2DATA",
														CPL_TYPE_DOUBLE, nwave);

		cpl_table_new_column_array (oi_table, "VIS2ERR",
														CPL_TYPE_DOUBLE, nwave);
		cpl_table_new_column (oi_table, "UCOORD", CPL_TYPE_DOUBLE);
		if(init) cpl_table_fill_column_window_double (oi_table, "UCOORD",
															  0, nrow * 6, 0.0);
		cpl_table_set_column_unit (oi_table, "UCOORD", "m");
		cpl_table_new_column (oi_table, "VCOORD", CPL_TYPE_DOUBLE);
		if(init) cpl_table_fill_column_window_double (oi_table, "VCOORD",
															  0, nrow * 6, 0.0);
		cpl_table_set_column_unit (oi_table, "VCOORD", "m");

		cpl_table_new_column_array (oi_table, "STA_INDEX", CPL_TYPE_INT, 2);
		cpl_table_set_column_savetype(oi_table, "STA_INDEX", CPL_TYPE_SHORT);
	}
	else if (!strcmp(oi_name, GRAVI_OI_T3_EXT)){
		oi_table = cpl_table_new(nrow * 4);
		cpl_table_new_column (oi_table, "TARGET_ID", CPL_TYPE_INT);
		cpl_table_set_column_savetype(oi_table, "TARGET_ID", CPL_TYPE_SHORT);
		if(init) cpl_table_fill_column_window_int (oi_table, "TARGET_ID",
															  0, nrow * 4, 0);
		cpl_table_new_column (oi_table, "TIME", CPL_TYPE_DOUBLE);
		if(init) cpl_table_fill_column_window_double (oi_table, "TIME", 0, nrow * 4, 0);
		cpl_table_set_column_unit (oi_table, "TIME", "s");
		cpl_table_new_column (oi_table, "MJD", CPL_TYPE_DOUBLE);
		if(init) cpl_table_fill_column_window_double (oi_table, "MJD", 0, nrow * 4, 0.0);
		cpl_table_set_column_unit (oi_table, "MJD", "d");
		cpl_table_new_column (oi_table, "INT_TIME", CPL_TYPE_DOUBLE);
		cpl_table_set_column_unit (oi_table, "INT_TIME", "s");
		if(init) cpl_table_fill_column_window_double (oi_table, "INT_TIME",
															0, nrow * 4, 0.0);
		cpl_table_new_column_array (oi_table, "T3AMP", CPL_TYPE_DOUBLE, nwave);
		cpl_table_new_column_array (oi_table, "T3AMPERR", CPL_TYPE_DOUBLE, nwave);
		cpl_table_new_column_array (oi_table, "T3PHI", CPL_TYPE_DOUBLE, nwave);
		cpl_table_set_column_unit (oi_table, "T3PHI", "deg");
		cpl_table_new_column_array (oi_table, "T3PHIERR",
													  CPL_TYPE_DOUBLE, nwave);
		cpl_table_set_column_unit (oi_table, "T3PHIERR", "deg");
		cpl_table_new_column (oi_table, "U1COORD", CPL_TYPE_DOUBLE);
		if(init) cpl_table_fill_column_window_double (oi_table, "U1COORD",
															0, nrow * 4, 0.0);
		cpl_table_set_column_unit (oi_table, "U1COORD", "m");
		cpl_table_new_column (oi_table, "V1COORD", CPL_TYPE_DOUBLE);
		if(init) cpl_table_fill_column_window_double (oi_table, "V1COORD",
															0, nrow * 4, 0.0);
		cpl_table_set_column_unit (oi_table, "V1COORD", "m");
		cpl_table_new_column (oi_table, "U2COORD", CPL_TYPE_DOUBLE);
		if(init) cpl_table_fill_column_window_double (oi_table, "U2COORD",
															0, nrow * 4, 0.0);
		cpl_table_set_column_unit (oi_table, "U2COORD", "m");
		cpl_table_new_column (oi_table, "V2COORD", CPL_TYPE_DOUBLE);
		if(init) cpl_table_fill_column_window_double (oi_table, "V2COORD",
															0, nrow * 4, 0.0);
		cpl_table_set_column_unit (oi_table, "V2COORD", "m");
		cpl_table_new_column_array (oi_table, "STA_INDEX", CPL_TYPE_INT, 3);
		cpl_table_set_column_savetype(oi_table, "STA_INDEX", CPL_TYPE_SHORT);
	}
	else if (!strcmp(oi_name, GRAVI_OI_FLUX_EXT)) {
		oi_table = cpl_table_new(nrow * 4);
		cpl_table_new_column (oi_table, "TARGET_ID", CPL_TYPE_INT);
		cpl_table_set_column_savetype(oi_table, "TARGET_ID", CPL_TYPE_SHORT);
		if(init) cpl_table_fill_column_window_int (oi_table, "TARGET_ID",
															0, nrow * 4, 0);
		cpl_table_new_column (oi_table, "TIME", CPL_TYPE_DOUBLE);
		cpl_table_set_column_unit (oi_table, "TIME", "s");
		if(init) cpl_table_fill_column_window_double (oi_table, "TIME", 0, nrow * 4, 0);
		cpl_table_new_column (oi_table, "MJD", CPL_TYPE_DOUBLE);
		cpl_table_set_column_unit (oi_table, "MJD", "d");
		if(init) cpl_table_fill_column_window_double (oi_table, "MJD", 0, nrow * 4, 0.0);
		cpl_table_new_column (oi_table, "INT_TIME", CPL_TYPE_DOUBLE);
		cpl_table_set_column_unit (oi_table, "INT_TIME", "s");
		if(init) cpl_table_fill_column_window_double (oi_table, "INT_TIME",
														  0, nrow * 4, 0.0);
		cpl_table_new_column_array (oi_table, "FLUX",CPL_TYPE_DOUBLE, nwave);
		cpl_table_set_column_unit (oi_table, "FLUX", "e");
		cpl_table_new_column_array (oi_table, "FLUXERR",CPL_TYPE_DOUBLE, nwave);
		cpl_table_set_column_unit (oi_table, "FLUXERR", "e");
		cpl_table_new_column (oi_table, "STA_INDEX", CPL_TYPE_INT);
		cpl_table_set_column_savetype(oi_table, "STA_INDEX", CPL_TYPE_SHORT);
		if(init) cpl_table_fill_column_window_int (oi_table, "STA_INDEX",
												0, nrow * 4, 0);
	}
	else {
		cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
								"The name of the table is incorrect");
		return NULL;
	}

	/* Set the flags */
	cpl_table_new_column_array (oi_table, "FLAG", CPL_TYPE_INT, nwave);
	cpl_table_set_column_savetype(oi_table, "FLAG", CPL_TYPE_BOOL);
	cpl_array * bool_array = cpl_array_new(nwave, CPL_TYPE_INT);
	cpl_array_fill_window_int(bool_array, 0, nwave, CPL_FALSE);
	cpl_table_fill_column_window_array (oi_table, "FLAG",
					0, cpl_table_get_nrow(oi_table), bool_array);
	cpl_array_delete (bool_array);

	CPLCHECK_NUL ("Cannot create oi_table");
	
    gravi_msg_function_exit(0);
	return oi_table;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Return the base of a region
 * 
 * @param  imaging_detector    The IMAGING_DETECTOR table
 * @param  region 		 	   The region (0..47)
 * 
 * @return The base corresponding to this region
 *
 * The function returns a integer from 0 to 5 (or -1 on failure)
 * If REGNAME is 12-A-C, the function returns 0 
 * If REGNAME is 21-A-C, the function returns 0 
 * If REGNAME is 23-A-C, the function returns 3
 * If REGNAME is 34-A-C, the function returns 5
 * ...
 */
/*---------------------------------------------------------------------------*/

int gravi_region_get_base (cpl_table *imaging_detector, int region)
{
	cpl_ensure (imaging_detector, CPL_ERROR_NULL_INPUT, -1);
	cpl_ensure (region>=0 && region<48, CPL_ERROR_ILLEGAL_INPUT, -1);

	const char * regname_i;
	int base = 0;
	const char * bases[6]={"12","13","14","23","24","34"};
	const char * bases_inv[6]={"21","31","41","32","42","43"};


	/* get the regname of region */
	regname_i = cpl_table_get_string(imaging_detector, "REGNAME", region);

	for (base=0; base<6; base++){
		if (!(strncmp(bases[base], regname_i, 2)) || !(strncmp(bases_inv[base], regname_i, 2))){
			return base;
		}
	}
	cpl_msg_error (cpl_func, "Cannot find the baseline of region %s", regname_i);
	return -1;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Return the sign of a base by looking at the PORT order
 * 
 * @param  imaging_detector    The IMAGING_DETECTOR table
 * @param  base  		 	   The base (0..5)
 * 
 * @return The sign of this base (-1/+1), or 0 or error
 *
 * The sign is defined by the order of the beam in the PORT column.
 * If the first PORT is the largest, the sign is -1. If the first
 * port is the smallest, the sign is +1.
 */
/*---------------------------------------------------------------------------*/

int gravi_region_get_base_sign (cpl_table *imaging_detector, int base)
{
	cpl_ensure (imaging_detector,   CPL_ERROR_NULL_INPUT, 0);
	cpl_ensure (base>=0 && base<=5, CPL_ERROR_ILLEGAL_INPUT, 0);

	cpl_size n_region = cpl_table_get_nrow (imaging_detector);
    const char * basename = GRAVI_BASE_NAME[base];

    /* Loop on region */
    for (cpl_size reg = 0; reg < n_region; reg ++) {
        const char * regname = cpl_table_get_string (imaging_detector, "REGNAME", reg);

        /* Check if this region match the requested base */
        if ( (regname[0] == basename[0] && regname[1] == basename[1]) ||
             (regname[0] == basename[1] && regname[1] == basename[0]) ) {

            /* Get the port array */
            const cpl_array * port = cpl_table_get_array (imaging_detector, "PORTS", reg);
            cpl_ensure (port, CPL_ERROR_ILLEGAL_INPUT, 0);

            /* Check the port order */
            if (cpl_array_get_int (port, 0, NULL) < cpl_array_get_int (port, 1, NULL)) return  1;
            if (cpl_array_get_int (port, 0, NULL) > cpl_array_get_int (port, 1, NULL)) return -1;
        }
    }

    return 0;
}


/*---------------------------------------------------------------------------*/
/**
 * @brief Return the polarisation id of a region
 * 
 * @param  imaging_detector    The IMAGING_DETECTOR table
 * @param  region 		 	   The region (0..47)
 * 
 * @return The polarisation id corresponding to this region
 *
 * The function returns 0 or 1 (or -1 on failure).
 * If REGNAME is 12-A-C, the function returns 0 
 * If REGNAME is 12-A-S, the function returns 0 
 * If REGNAME is 12-A-P, the function returns 1
 */
/*---------------------------------------------------------------------------*/

int gravi_region_get_pol (cpl_table *imaging_detector, int region)
{
	cpl_ensure (imaging_detector, CPL_ERROR_NULL_INPUT, -1);
	cpl_ensure (region>=0 && region<48, CPL_ERROR_ILLEGAL_INPUT, -1);

	/* get the regname of region */
	const char * regname_i;
	regname_i = cpl_table_get_string (imaging_detector, "REGNAME", region);

	/* check its name */
	if (regname_i[5] == 'S' || regname_i[5] == 'C') return 0;
	if (regname_i[5] == 'P') return 1;
	
	cpl_msg_error (cpl_func, "Cannot find the polarisation of region %s", regname_i);
	return -1;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Return the phase id of a region
 * 
 * @param  imaging_detector    The IMAGING_DETECTOR table
 * @param  region 		 	   The region (0..47)
 * 
 * @return The phase id corresponding to this region
 *
 * The function returns 0 to 3 (or -1 on failure).
 * If REGNAME is 12-A-C, the function returns 0 
 * If REGNAME is 12-B-C, the function returns 1 
 * If REGNAME is 12-C-C, the function returns 2
 * If REGNAME is 12-D-C, the function returns 3
 */
/*---------------------------------------------------------------------------*/

int gravi_region_get_phaseid (cpl_table *imaging_detector, int region)
{
	cpl_ensure (imaging_detector, CPL_ERROR_NULL_INPUT, -1);
	cpl_ensure (region>=0 && region<48, CPL_ERROR_ILLEGAL_INPUT, -1);

	/* get the regname of region */
	const char * regname_i;
	regname_i = cpl_table_get_string (imaging_detector, "REGNAME", region);

	/* check its name */
	if (regname_i[3] == 'A') return 0;
	if (regname_i[3] == 'B') return 1;
	if (regname_i[3] == 'C') return 2;
	if (regname_i[3] == 'D') return 3;
	
	cpl_msg_error (cpl_func, "Cannot find the phaseid of region %s", regname_i);
	return -1;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Return the telescope id (0,1,2,3) in a beam of a region
 *
 * @param  imaging_detector    The table IMAGING_DETECTOR
 * @param  region 		 	   The region to query
 * @param  beam 		 	   The beam to query (0 or 1)
 *
 * If REGNAME is 34-A-C, the function return 2 for beam 0 and 3 for beam 1 
 * If REGNAME is 21-B-S, the function return 1 for beam 0 and 0 for beam 1 
 * ...
 */
/*---------------------------------------------------------------------------*/

int gravi_region_get_tel (cpl_table *imaging_detector, int region, int beam)
{
	cpl_ensure (imaging_detector,       CPL_ERROR_NULL_INPUT, -1);
	cpl_ensure (region>=0 && region<48, CPL_ERROR_ILLEGAL_INPUT, -1);
	cpl_ensure (beam==0 || beam==1,     CPL_ERROR_ILLEGAL_INPUT, -1);

	/* get the regname of region */
	const char * regname_i;
	regname_i = cpl_table_get_string (imaging_detector, "REGNAME", region);
    cpl_ensure (regname_i, CPL_ERROR_ILLEGAL_INPUT, -1);
    
    /* Memory to read */
    char *ptr, str[2] = "0";
    str[0] = regname_i[beam];

    /* Read the telescope */
    long tel = strtol (str, &ptr, 10) - 1;

	cpl_ensure (tel>=0 && tel<4, CPL_ERROR_ILLEGAL_INPUT, -1);
	return (int)tel;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Return the phase character of a region
 * 
 * @param  imaging_detector    The IMAGING_DETECTOR table
 * @param  region 		 	   The region (0..47)
 * 
 * @return The phase char corresponding to this region
 *
 * The function returns 'A', 'B', 'C' or 'D'
 * If REGNAME is 12-A-C, the function returns 'A' 
 * If REGNAME is 12-B-C, the function returns 'B' 
 * If REGNAME is 12-C-C, the function returns 'C'
 * If REGNAME is 12-D-C, the function returns 'D'
 */
/*---------------------------------------------------------------------------*/

char gravi_region_get_phase (cpl_table *imaging_detector, int region)
{
	cpl_ensure (imaging_detector, CPL_ERROR_NULL_INPUT, -1);
	cpl_ensure (region>=0 && region<48, CPL_ERROR_ILLEGAL_INPUT, -1);

	/* get the regname of region */
	const char * regname_i;
	regname_i = cpl_table_get_string (imaging_detector, "REGNAME", region);

	/* check its name */
	// return regname_i[3];
	if (regname_i[3] == 'A') return 'A';
	if (regname_i[3] == 'B') return 'B';
	if (regname_i[3] == 'C') return 'C';
	if (regname_i[3] == 'D') return 'D';
	
	cpl_msg_error (cpl_func, "Cannot find the phase of region %s", regname_i);
	return -1;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Find the region matching base, phase and pol
 * 
 * @param  imaging_detector    The IMAGING_DETECTOR table
 * @param  base 		 	   The requested base (0..5)
 * @param  phase 		 	   The requested phase ('A','B','C' or 'D')
 * @param  pol   		 	   The requested pol (0 or 1)
 * 
 * @return The region id (or -1 on failure)
 */
/*---------------------------------------------------------------------------*/

int gravi_get_region (cpl_table *img_det, int base, char phase, int pol)
{
  cpl_ensure (img_det, CPL_ERROR_NULL_INPUT, -1);

  for (cpl_size reg = 0 ; reg < cpl_table_get_nrow (img_det); reg++) {
	if ( gravi_region_get_pol (img_det, reg) == pol &&
		 gravi_region_get_phase (img_det, reg) == phase &&
		 gravi_region_get_base (img_det, reg) == base ) return reg;
  }
  
  return -1;
}


/*---------------------------------------------------------------------------*/
/**
 * @brief    Get the number of spectral element between lambdamin et lambdamax
 * @param  wave_table The table loaded from the IMAGING_DETECTOR field
 * @param  lambda_min 		 	   The number of the first telescope
 * @param  lambda_max 		 	   The number of the second telescope
 * 
 * @return   the number of spectral element
 * 
 * The function returns a integer
 */
/*---------------------------------------------------------------------------*/

int gravi_wave_get_nlambda(cpl_table *wave_data, double lambda_min, double lambda_max)
{
	cpl_ensure (wave_data, CPL_ERROR_NULL_INPUT, -1);
	cpl_ensure (lambda_max>lambda_min, CPL_ERROR_ILLEGAL_INPUT, -1);
	cpl_ensure (lambda_min>1e-6 && lambda_min<3e-6, CPL_ERROR_ILLEGAL_INPUT, -1);
	cpl_ensure (lambda_max>1e-6 && lambda_max<3e-6, CPL_ERROR_ILLEGAL_INPUT, -1);
	
	int n_wave = cpl_table_get_column_depth (wave_data, "DATA1");
	int n_region = cpl_table_get_ncol (wave_data);
    int n_element;
	double wave, max_wave=0., min_wave = 3*pow(10, -6);
    double sc_medium_res, sc_high_res, sc_low_res;
	cpl_size ind_min=0, ind_max=n_wave-1;
	char* regname;

	/* find first pixel */
	while ((lambda_min-max_wave >= 0.0001*pow(10, -6)) && (ind_min < n_wave)) {
		for (int region=0; region < n_region; region++){
			regname = GRAVI_DATA[region];
			wave=cpl_array_get(cpl_table_get_array(wave_data, regname, 0), ind_min, NULL);
			if (wave >=max_wave) max_wave=wave;
		}
		ind_min++;
	}
	ind_min--;

	CPLCHECK_INT ("Cannot get first pixel");

	/* find last pixel */
	while ((lambda_max-min_wave <= -0.0001*pow(10, -6)) && (ind_max >= 0)) {
		for (int region=0; region < n_region; region++){
			regname = GRAVI_DATA[region];
			wave=cpl_array_get(cpl_table_get_array(wave_data, regname, 0), ind_max, NULL);
			if (wave <=min_wave) min_wave=wave;
		}
		ind_max--;
	}
	ind_max++;
	
	CPLCHECK_INT ("Cannot get last pixel");

	/* Compute the mean spectral resolution */
	cpl_array * all_steps = cpl_array_new (n_region * (n_wave-1), CPL_TYPE_DOUBLE);
	for (int region=0; region < n_region; region++) {
	  const cpl_array * data = cpl_table_get_array (wave_data, GRAVI_DATA[region], 0);
	  for (wave=0; wave < n_wave-1; wave ++)
		cpl_array_set (all_steps, region*(n_wave-1)+wave, fabs (cpl_array_get (data, wave+1, NULL) - cpl_array_get (data, wave, NULL) ));
	}
    
	double res = cpl_array_get_median (all_steps);
	FREE (cpl_array_delete, all_steps);
    
    /* n_element valid for the FT case or default if SC spectral channel change too much from nominal */
    n_element = ind_max-ind_min+1;

    /* Constant wavelength channel width for FT - pkervell 16/03/2016 */
    /* Nominal spectral channel width for FT in meters */
    /* ft_res = 0.082308 * pow(10, -6);*/
    /* Allow for a +/- 5% variation margin around nominal value
     if (( res >= 0.95*ft_res) && (res <= 1.05*ft_res)) {
     res = ft_res;
     int n_element = ind_max-ind_min+1;
     cpl_msg_info (cpl_func, "FT spectral resolution res = %e m, n_element = %i", res, n_element);
     }*/
    
    /* Constant wavelength channel width for SC - pkervell 16/03/2016 */
    /* Nominal spectral channel width for SC in meters */
    sc_low_res = 0.0418181818 * pow(10, -6);
    sc_medium_res = 0.0022009569 * pow(10, -6);
    sc_high_res =   0.0002642248 * pow(10, -6);

    /* SC low resolution case */
    /* Not implemented in low resolution due to tighter wavelength limits than 1.99-2.45 microns */
    if (( res >= 0.95*sc_low_res) && (res <= 1.05*sc_low_res)) {
        n_element = round(((lambda_max-lambda_min)/sc_low_res) + 2); // plus two, because it is always better to have a little margin at low resolution, especially with the 3 pixels interpolation
        cpl_msg_info (cpl_func, "Chosen SC LOW spectral resolution element = %e m, n_element = %i", sc_low_res, n_element);
    }
    
    /* SC medium resolution case */
    /* Allow for a +/- 5% variation margin around nominal value */
    if ((res >= 0.95*sc_medium_res ) && (res <= 1.05*sc_medium_res)) {
        n_element = round(((lambda_max-lambda_min)/sc_medium_res) + 1);
        cpl_msg_info (cpl_func, "Chosen SC MEDIUM fixed spectral element = %e m, n_element = %i", sc_medium_res, n_element);
    }

    
    /* SC high resolution case */
    /* Allow for a +/- 10% variation margin around nominal value */
    if ((res >= 0.9*sc_high_res ) && (res <= 1.1*sc_high_res)) {
        n_element = round(((lambda_max-lambda_min)/sc_high_res) + 1);
        cpl_msg_info (cpl_func, "Chosen SC HIGH fixed spectral element = %e m, n_element = %i", sc_high_res, n_element);
    }
	
	CPLCHECK_INT ("Cannot compute nlambda");
    
    cpl_msg_info (cpl_func, "min=%e, max=%e, n=%f, res=%e, n=%i", lambda_min, lambda_max, (lambda_max-lambda_min)/res, res, n_element);

    // add security here, or the interpolation will crash sometimes
    if (n_element > n_wave)
        n_element = n_wave;
    
    return n_element;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Compute startx and nx of the illuminated part of the image
 * @param img_profile     The input image to look at
 * @return A pointer to int[2] with startx and nx
 *
 * The function collapse the vertical (y) dimension and run a median
 * filter over x. Then the limits are defined by the region when the 
 * the median filter is larger than 0.1 x MAX
 */
/*---------------------------------------------------------------------------*/

int * gravi_image_extract_dimension (cpl_image * img_profile)
{
	cpl_ensure (img_profile, CPL_ERROR_NULL_INPUT, NULL);

//    /* Collapse */
//    cpl_image * collapse_img = cpl_image_collapse_create (img_profile, 0);
//    cpl_size nx = cpl_image_get_size_x (collapse_img);
//
//    /* Median filter */
//    cpl_size size = nx < 60 ? 2 : 5;
//    cpl_mask * kernel = cpl_mask_new (size, 1); 
//    cpl_mask_not (kernel);
//
//    cpl_image * filt_img = cpl_image_duplicate (collapse_img);
//    cpl_image_filter_mask (filt_img, collapse_img, kernel,
//                           CPL_FILTER_MEDIAN, CPL_BORDER_FILTER);
//
//    /* Search for limits */
//    double max = cpl_image_get_max (filt_img) * 0.1;
//    double imin = nx, imax = 0;
//    for (cpl_size i = 1; i <= nx; nx++) {
//        if (cpl_image_get (filt_img, i, 1) > max) {
//            if (i < imin) imin = i;
//            if (i > imax) imax = i;
//        }
//    }
//    
//    /* Return */
//    int dim = cpl_malloc (2 * sizeof (int));
//    dim[0] = imin;
//    dim[1] = imax - imin + 1;
//        
	double sig;
	cpl_vector * vector, * vect, * vect_mean;

	cpl_size nx = cpl_image_get_size_x (img_profile);
	cpl_size ny = cpl_image_get_size_y (img_profile);
    
	vect_mean = cpl_vector_new_from_image_row (img_profile, 1);
	for (cpl_size i = 2; i <= ny; i++) {
		vect = cpl_vector_new_from_image_row (img_profile, i);
		cpl_vector_add (vect_mean, vect);
		cpl_vector_delete (vect);
	}

	cpl_vector_divide_scalar (vect_mean, ny);
	if (nx < 60) { // case for LOW res
		vector = cpl_vector_filter_median_create (vect_mean, 2);
	}
	else{ // case for MED and HIGH
		vector = cpl_vector_filter_median_create (vect_mean, 5);
	}

	double max = cpl_vector_get_max  (vector);
	sig = max * 0.10; // cut the edge to 10 % of the flux
    
	double sum = 0;
	int i_2 = nx, i_1 = 0;
	for (cpl_size i = 1; i < nx - 1; i++){
		if (cpl_vector_get (vector, i) > sig){
			sum ++;
			if ((cpl_vector_get (vector, i - 1) > sig) && (cpl_vector_get (vector, i + 1) < sig))
				i_2 = i;
		}
	}

	for (cpl_size i = 1; i < nx - 1; i++){
		if (cpl_vector_get (vector, nx - i) > sig){
			sum ++;
			if ((cpl_vector_get (vector, nx - i - 1) < sig) && (cpl_vector_get (vector, nx - i + 1) > sig))
				i_1 = nx - i;
		}
	}

	/* increase by 1 pixels */
	i_1 -= 1;
	i_2 += 2;
	if (i_1 < 0) i_1 = 0;
	if (i_2 >= nx) i_2 = nx-1;

    /* Fill output */
	int * dim = cpl_malloc (2 * sizeof (int));
	dim [0] = i_1 + 1;
	dim [1] = i_2 - i_1;

	cpl_vector_delete (vect_mean);
	cpl_vector_delete (vector);

	return dim;
}


/*----------------------------------------------------------------------------*/
/**
 * @brief Retrieve STA_INDEX corresponding to a given input
 * @param gravi_input		GRAVITY input index (1..4)
 * 
 * @return   STA_INDEX
 */
/*----------------------------------------------------------------------------*/

short gravi_sta_index(int gravi_input, cpl_table * optical_train_table, cpl_table * array_geometry_table)
{
  gravi_msg_function_start(0);
  cpl_ensure (optical_train_table,  CPL_ERROR_NULL_INPUT, 0);
  cpl_ensure (array_geometry_table, CPL_ERROR_NULL_INPUT, 0);
  
  short lab_input=9;
  switch (gravi_input) {
  case 1:
    lab_input = GRAVI_LABINPUT_1;
    break;
  case 2:
    lab_input = GRAVI_LABINPUT_2;
    break;
  case 3:
    lab_input = GRAVI_LABINPUT_3;
    break;
  case 4:
    lab_input = GRAVI_LABINPUT_4;
    break;
  default:
    cpl_error_set_message (cpl_func, CPL_ERROR_ILLEGAL_INPUT,
						   "impossible GRAVITY input");
    return 0;
  }

  cpl_msg_debug(cpl_func, "Start function");
  
  /* Get tel_name corresponding to lab_input in optical_train_table */
  int row, null;
  for (row=0; row < cpl_table_get_nrow(optical_train_table); ++row) {
    if (cpl_table_get_int(optical_train_table, "VALUE2", row, &null) ==
	lab_input) {
      break;
    }
  }
  if (row >= cpl_table_get_nrow(optical_train_table)) {
    cpl_error_set_message(cpl_func, CPL_ERROR_DATA_NOT_FOUND,
			  "lab input not found in optical train table");
    return 0;
  }

  char * tel_name;
  tel_name = cpl_sprintf ("%s",cpl_table_get_string(optical_train_table, "TEL_NAME", row));
  
  /* trim tel_name */
  {
	int cur = strlen (tel_name);
	while (cur) {
	  if (tel_name[cur-1] == ' ') {
		--cur;
		tel_name[cur]=0;
	  } else break;
	}
  }
  
  /* Get STA_INDEX corresponding to tel_name in array_geometry_table */
  char * cur_name = NULL;
  int found=0;
  
  for (row = 0; (!found) && (row < cpl_table_get_nrow(optical_train_table)); ++row) {
	
      cur_name = cpl_sprintf ("%s",cpl_table_get_string(array_geometry_table, "TEL_NAME", row));
	
    /* trim cur_name */
    int cur = strlen (cur_name);
    while (cur) {
      if (cur_name[cur-1] == ' ') {
	--cur;
	cur_name[cur]=0;
      } else break;
    }
    if (!strcmp (tel_name, cur_name)) found=1;
    FREE (cpl_free, cur_name);
    if (found) break;
  }

  FREE (cpl_free, tel_name);
 
  if (!found) {
    cpl_error_set_message(cpl_func, CPL_ERROR_DATA_NOT_FOUND,
			  "telescope not found in array geometry");
    return 0;
  }


  short res = cpl_table_get_int(array_geometry_table, "STA_INDEX", row, &null);

  cpl_msg_debug(cpl_func, "gravi_input=%i; lab_input=%i; sta_index=%i", gravi_input, lab_input, res);
  cpl_msg_debug(cpl_func, "End function");
  
  gravi_msg_function_exit(0);
  return res;
}

int gravi_get_shutter (cpl_propertylist * header, int tel)
{
  cpl_ensure (header, CPL_ERROR_ILLEGAL_INPUT, -1);
  
  char key[90];
  sprintf (key, SHUTTER_KEY"%d ST", tel+1);

  /* Try the old and the new shutter name */
  if ( !cpl_propertylist_has (header, key) ) {
	sprintf (key, GRAVI_SHUTTER_KEY"1%d ST", tel+1);
	
	if ( !cpl_propertylist_has (header, key) ) {
	  cpl_error_set_message (cpl_func, CPL_ERROR_ILLEGAL_INPUT, "Missing shutter");
	  return -1;
	}
  }

  return cpl_propertylist_get_bool (header, key);
}

int gravi_check_shutter (cpl_propertylist * header, int t0, int t1, int t2, int t3)
{
  cpl_ensure (header, CPL_ERROR_ILLEGAL_INPUT, -1);

  if ( gravi_get_shutter (header, 0) == t0 &&
	   gravi_get_shutter (header, 1) == t1 &&
	   gravi_get_shutter (header, 2) == t2 &&
	   gravi_get_shutter (header, 3) == t3 )
	return 1;

  return 0;
}

int gravi_get_shutter_id (cpl_propertylist * header)
{
    cpl_ensure (header, CPL_ERROR_ILLEGAL_INPUT, -1);
    
    if (gravi_check_shutter (header,1,0,0,0)) return 0;
    if (gravi_check_shutter (header,0,1,0,0)) return 1;
    if (gravi_check_shutter (header,0,0,1,0)) return 2;
    if (gravi_check_shutter (header,0,0,0,1)) return 3;

    cpl_error_set_message (cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                           "More than one shutter open, or none");
    return -1;
}

int gravi_get_shutter_baseid (cpl_propertylist * header)
{
    cpl_ensure (header, CPL_ERROR_ILLEGAL_INPUT, -1);
    
    if (gravi_check_shutter (header,1,1,0,0)) return 0;
    if (gravi_check_shutter (header,1,0,1,0)) return 1;
    if (gravi_check_shutter (header,1,0,0,1)) return 2;
    if (gravi_check_shutter (header,0,1,1,0)) return 3;
    if (gravi_check_shutter (header,0,1,0,1)) return 4;
    if (gravi_check_shutter (header,0,0,1,1)) return 5;

    cpl_error_set_message (cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                           "More than two shutter open, or less");
    return -1;
}

int gravi_data_check_shutter_beam (gravi_data ** datas, int nb_datas)
{
    cpl_ensure (datas, CPL_ERROR_ILLEGAL_INPUT, -1);

    /* Wrong number of flats */
    if ( nb_datas != 4 ) return 0;

    /* Count the single-shutter open for each tel */
    int is_open[] = {0,0,0,0};
    for (cpl_size f = 0 ; f < 4 ; f++) {
        cpl_propertylist * header = gravi_data_get_header (datas[f]);
        is_open[0] += gravi_check_shutter (header,1,0,0,0);
        is_open[1] += gravi_check_shutter (header,0,1,0,0);
        is_open[2] += gravi_check_shutter (header,0,0,1,0);
        is_open[3] += gravi_check_shutter (header,0,0,0,1);
    }

    cpl_msg_info (cpl_func,"is_open = %i %i %i %i", is_open[0], is_open[1], is_open[2], is_open[3]);

    /* Check all open, and once only */
    int check = 1;
    for (cpl_size f = 0 ; f < 4 ; f++) if (is_open[f] != 1) check *= 0;

    return check;
}


cpl_size gravi_spectrum_get_nwave (const cpl_table * table)
{
    return cpl_table_get_column_depth (table,"DATA1");
}

cpl_size gravi_spectrum_get_nregion (const cpl_table * table)
{
    cpl_ensure (table, CPL_ERROR_NULL_INPUT, -1);

    cpl_size nregion = 0;
    for (cpl_size col = 0; col < 48 ; col++) {
        nregion += cpl_table_has_column (table, GRAVI_DATA[col]);
    }

    return nregion;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Return the total flux in DATA# regions
 *
 * @param table     The input table
 * @return The flux
 *
 * Compute the total flux on DATA# regions by looping on
 * regions and rows.
 */
/*----------------------------------------------------------------------------*/

double gravi_spectrum_get_flux (const cpl_table * table)
{
    cpl_ensure (table, CPL_ERROR_NULL_INPUT, 0);
    
    cpl_size n_reg = gravi_spectrum_get_nregion (table);
    cpl_size n_row = cpl_table_get_nrow (table);
    double flux = 0.0;
    
    for (cpl_size reg = 0; reg < n_reg ; reg++) {
        
        cpl_size size = cpl_table_get_column_depth (table, GRAVI_DATA[reg]);
        cpl_array ** arrays = cpl_table_get_data_array ((cpl_table *)table,GRAVI_DATA[reg]);
        
        for (cpl_size row=0; row < n_row ; row++)
            flux += cpl_array_get_mean (arrays[row]) * size;
    }
    
    return flux;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Compute the envelope value
 * 
 * @param opd        Input vector with GD in [m]
 * @param wave       Input wave channel
 * @param nwave     Number of wave channel (define spectral resolution)
 *
 * @return A newly allocated vector with the computed envelope
 *
 * The vector is filled with the envelope amplitude at each OPD_GD position
 * for the spectral channel wave, assuming the K-band is split into n_wave
 * channel. Thus n_wave defines the spectral resolution (coherence length).
 * The enveloppe is modeled with a Gaussian function.
 *
 * Assume wave are ordered increasingly.
 */
/*---------------------------------------------------------------------------*/

cpl_vector * gravi_compute_envelope (const cpl_vector *opd, int wave, int nwave)
{
    cpl_ensure (opd, CPL_ERROR_NULL_INPUT, NULL);
    cpl_ensure (wave>=0 && wave<nwave, CPL_ERROR_ILLEGAL_INPUT, NULL);
    
	/* Compute delta_lambda and lambda from experience */
	double delta_lambda = (nwave > GRAVI_LBD_FTSC) ? 0.45 / nwave * 3 : 0.13;
    double lambda = 2.0 + 0.45 / nwave * wave;
    
	 /* Compute coherent length */
	double coh_len= (lambda*lambda) / delta_lambda * 1.e-6;

	long nrow = cpl_vector_get_size (opd);
	cpl_vector * envelope = cpl_vector_new (nrow);

	/* Gaussian enveloppe */
	for (long row = 0; row < nrow; row++){
		double value = cpl_vector_get (opd, row);
        cpl_vector_set (envelope, row, exp(-1*(value*value)/(coh_len*coh_len/2.)));
        CPLCHECK_NUL ("Cannot compute envelope");
	}

	return envelope;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Return the total flux in imagelist
 *
 * @param  imglist     The input imagelist
 * @return The flux
 *
 * Compute the total flux by looping on image.
 */
/*----------------------------------------------------------------------------*/

double gravi_imagelist_get_flux (const cpl_imagelist * imglist)
{
    cpl_ensure (imglist, CPL_ERROR_NULL_INPUT, 0);

    cpl_size size = cpl_imagelist_get_size (imglist);
    double flux = 0.0;
    
    for (cpl_size row = 0; row < size; row ++)
        flux += cpl_image_get_flux (cpl_imagelist_get_const (imglist, row));

    return flux;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Return the longuest sequence with constant LKDT
 *
 * @param  oi_table   The table to explore
 * @param  ntel       The number of sample per obs (4 or 6)
 * @param  first      Return the first obs of the sequence 
 * @param  nobs       Return the nb of obs in the sequence
 *
 * Search for the longuest sequence for which the LKDT_MET_FC
 * are the same for the ntel beams. If several are of the same
 * length, the first is returned.
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_lkdt_get_sequence (cpl_table * oi_table,
                                        cpl_size ntel,
                                        cpl_size *first, cpl_size *nobs)
{
    cpl_ensure_code (oi_table, CPL_ERROR_NULL_INPUT);

    cpl_size cfirst = 0, cnobs = 0;
    cpl_size nrow = cpl_table_get_nrow (oi_table) / ntel;
    double * lkdt = cpl_table_get_data_double (oi_table, "LKDT_MET_FC");
    cpl_ensure_code (lkdt, CPL_ERROR_ILLEGAL_INPUT); 

    for (cpl_size row = 1; row < nrow; row++) {
        
        /* If change of lockdate, we start counting again */
        for (cpl_size tel = 0; tel < ntel; tel++)
            if (lkdt[row*ntel+tel] != lkdt[(row-1)*ntel+tel]) cfirst = row;

        /* If more than before, we save this sequence */
        cnobs = row - cfirst + 1;
        if (cnobs > *nobs) {*nobs = cnobs; *first = cfirst;}
    }
    
    return CPL_ERROR_NONE;
}

int gravi_conf_get_iss (int gravi_beam, cpl_propertylist * header)
{
    cpl_ensure (gravi_beam >=0 && gravi_beam < 4, CPL_ERROR_ILLEGAL_INPUT, 0);
    cpl_ensure (header, CPL_ERROR_NULL_INPUT, 0);

    int iss = 0;
    char name[100];

    for (iss = 1; iss<=4; iss++) {
        sprintf (name, "ESO ISS CONF INPUT%i", iss);
        if ( cpl_propertylist_get_int (header, name) ==
             GRAVI_LABINPUT[gravi_beam]) return iss;
    }
    
    cpl_error_set_message (cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                          "Could not find CONF INPUT for beam %i",
                          gravi_beam);

    return 0;
}

const char * gravi_conf_get_telname (int gravi_beam, cpl_propertylist * header)
{
    cpl_ensure (gravi_beam >=0 && gravi_beam < 4, CPL_ERROR_ILLEGAL_INPUT, NULL);
    cpl_ensure (header, CPL_ERROR_NULL_INPUT, NULL);
    char name[100];

    for (int iss = 1; iss<=4; iss++) {
        
        sprintf (name, "ESO ISS CONF INPUT%i", iss);
        if ( !cpl_propertylist_has (header, name)) continue;
        
        if ( cpl_propertylist_get_int (header, name) ==
             GRAVI_LABINPUT[gravi_beam]) {

            sprintf (name, "ESO ISS CONF T%iNAME", iss);
            if ( !cpl_propertylist_has (header, name)) continue;
            
            return cpl_propertylist_get_string (header, name);
        }
    }
    
    cpl_msg_debug (cpl_func,"Telscope name cannot be find in header");
    return NULL;
}

                                                                                               
cpl_error_code gravi_dump_the_boss (double ra, double dec)
{
    /* GC coordinates in [rad] */
    double c_ra  = gravi_ra_to_rad ("17:45:40.03599");
    double c_dec = gravi_dec_to_rad ("-29:00:28.1699");

    /* Distance in [rad] */
    double dist = acos ( sin (c_dec) * sin (dec) + cos (c_dec) * cos (dec) * cos (ra - c_ra) );

    /* Distance in [as] */
    dist = dist / CPL_MATH_PI * 180 * 3600;
    cpl_msg_info (cpl_func, "Distance from GC: %g [as]", dist);
    
    if (dist > 20.0) return CPL_ERROR_NONE;
    
    cpl_msg_info (cpl_func, "*****************************************************************************************");
    cpl_msg_info (cpl_func, "                                                                                         ");
    cpl_msg_info (cpl_func, "                You are reducing data on the Galactic Center !!!!                        ");
    cpl_msg_info (cpl_func, "                                                                                         ");
    cpl_msg_info (cpl_func, "                                      MAY THE FRINGES BE WITH YOU !!!!                   ");
    cpl_msg_info (cpl_func, "                                                                                         ");
    cpl_msg_info (cpl_func, "                                                                                         ");
    cpl_msg_info (cpl_func, "                                                                                         ");
    cpl_msg_info (cpl_func, "                                `...`  ```  ` ```` ````                                  ");
    cpl_msg_info (cpl_func, "                          `....:/:::--://+/:--:/::--.-``````                             ");
    cpl_msg_info (cpl_func, "                      `.-://::.-----::/+oohyso/++:/:::.-.-://-.`                         ");
    cpl_msg_info (cpl_func, "                 ```.:/:+o+/::/+ooossoo+oosso++::+/:://:+:osso//-``                      ");
    cpl_msg_info (cpl_func, "               ```.:+ooos+:/+shddhhhysosooso+o++ooo/+sso+////:///:-`                     ");
    cpl_msg_info (cpl_func, "              ``.-:+yso/::+shdhhddhysssyhyhyyyysoyso+shys+o////::++:.`                   ");
    cpl_msg_info (cpl_func, "            `.-/+++syo//+syhdddddhhyyhddhhyyyhysyhyysyhdss+o+/o++:/+/:-.`                ");
    cpl_msg_info (cpl_func, "           `:/ss++oo/osyhddddmmmdddhddmdhhyhhhhhhhhyssyyoo/s+:+hoo/:://-.`               ");
    cpl_msg_info (cpl_func, "          `-+ss+oo++oyddmmmmmNNmddddddddhhhdhhhhhhyyyyy+o++so/ssyys+/+o+-.`              ");
    cpl_msg_info (cpl_func, "         `.+oysssosshdmmmmmmmdmddddhhddddddysssssyysyso+///+osysyoyyo+os+:.              ");
    cpl_msg_info (cpl_func, "         `-/syosyyydmmdddhhhyhhhhhhyhhddhso////+ssoo+/:----/:+oshs+ss++sso/`             ");
    cpl_msg_info (cpl_func, "         `:oyyohdhdmmdyyyysssssooo+oo+++/:----:/+::::-...---../ohho/so/osoo:`            ");
    cpl_msg_info (cpl_func, "        `./yyyhmdNNmhyssssoo+++//:::::-------------.......--..-/yyy/+ssoooo+-`           ");
    cpl_msg_info (cpl_func, "        `.ohshmmmNNdyoooooo++///::-:------------..........-----:shss+oooo+so/.           ");
    cpl_msg_info (cpl_func, "        `-syshmmNNNdsooooo+++//:::--------------...........----:+hhsss/++++o/.           ");
    cpl_msg_info (cpl_func, "         `+syhdmNNNhsoooo++++/::::-----.-------..............---/dhyso++++oo:.           ");
    cpl_msg_info (cpl_func, "          -yhhmmNNNhyooo+++++//::::--....-------...-........---:/ddhs/+/o+o+:``          ");
    cpl_msg_info (cpl_func, "          `shdmNMNmdysoo+++++/////:::------------------......---:smh+/+/++o/.`           ");
    cpl_msg_info (cpl_func, "          `oydmNMMmhyo+++++ooo+++++///::---:::::://///////:...--:sNd++o:/o+:.            ");
    cpl_msg_info (cpl_func, "          `+ydmmNMNhy++ooyhhhhdddddhyo++/:-://+shdddddhhyys+-..-:sNdoss//o/-`            ");
    cpl_msg_info (cpl_func, "           osmdmNMmhsoyyddmmmdddmmmmdhyy+:::/oshdmmmddhhhyso+/--:sNmhhy++o:``            ");
    cpl_msg_info (cpl_func, "           +hmmNNNmyohhhhhhhyyhhhdhhdhdh+:--:oyhddhdddhhyysso/:::smmdhhoo/-`             ");
    cpl_msg_info (cpl_func, "          `yhmNNmmhsyhyhhdmmmNmdyhdddddhs:--+ydyhhdhdmmmddhso+/y/+yhhhs+//+`             ");
    cpl_msg_info (cpl_func, "          `yhhdmNNdoshshdddhhddhssyhydmd+-.-/ydysyysyhhyosyhy+:s-/dmh+///ss.             ");
    cpl_msg_info (cpl_func, "          .syyhdNNh+oysssyyyyyssssyyyhhy/-.-:ossoosyysso++//////-:sds/:--/o.             ");
    cpl_msg_info (cpl_func, "          -sssyhddo+oooosssssssssssssyys/-..-:///++osoo++///::----/s/+/---+`             ");
    cpl_msg_info (cpl_func, "          `oosshyy++ooo+++ooooo++++ossso+--.-::::-::/+///::-------:+/so+/--`             ");
    cpl_msg_info (cpl_func, "           :+ossshsosooo++////:://+osso+/:-..-:::-----------------//-:+/:-.              ");
    cpl_msg_info (cpl_func, "           -+ososhhsyssoo+///::::/+sso++:--...-::---------:::--::/o:--://-`              ");
    cpl_msg_info (cpl_func, "           -++osshmhyyysoo+//////+osso++/---...-::---------::::::sy/::::--`              ");
    cpl_msg_info (cpl_func, "           `o++sshdhyyyyysooo++/+osyssso///:--:/:+/---:---::::::/ys/:/:.--`              ");
    cpl_msg_info (cpl_func, "            :/osyhhyyyyyyssooo++oyyyhdmhysssoohy++so/:::::/:::::::+o+::...               ");
    cpl_msg_info (cpl_func, "            `+oossyysyyhyyssooosyssshdmmddddyso/::/oyo/////////::://::---`               ");
    cpl_msg_info (cpl_func, "             .+sssoosyyyyyysoosyssssyhhhhhhhs++/////oys++++++//::/+///:.`                ");
    cpl_msg_info (cpl_func, "               `..`.syhyyyysssyyssssssssooss+//::////ohs++++++/::/.```                   ");
    cpl_msg_info (cpl_func, "                    +yyyhyysssssyyssssyssoooooo++///+osyo++++//:-:`                      ");
    cpl_msg_info (cpl_func, "                    -syyhyysossyhdmmdhhyyyyysyshddhyssss+++++/::::                       ");
    cpl_msg_info (cpl_func, "                    `syyyyysosssyhddyso+//+++o+oso+/++o+//+++/:::.                       ");
    cpl_msg_info (cpl_func, "                     :yyyyysossssssssso++++++/++//////++/+++//::.                        ");
    cpl_msg_info (cpl_func, "                     `+yyyyysssssssyssssssosso+///////++++++///-`                        ");
    cpl_msg_info (cpl_func, "                      `syyyhyysssssyyyysssssso++///////++/++//:.                         ");
    cpl_msg_info (cpl_func, "                       .syhhhyysyyssssooooo+++///////////+++//:.                         ");
    cpl_msg_info (cpl_func, "                        .shhyyyssssooo++///:::::/::/://++o+//:-.                         ");
    cpl_msg_info (cpl_func, "                         `+yyyyssooooo+/////:::::::://++oo+/::-.                         ");
    cpl_msg_info (cpl_func, "                           -yhyyssssso++/+++++//////++osso//::--`                        ");
    cpl_msg_info (cpl_func, "                            -hhhhyyyyysooosssooo++oosyyso+/:::--.                        ");
    cpl_msg_info (cpl_func, "                             /hddddddhhhyyyyyyyyyyyhhyys+//:::---.                       ");
    cpl_msg_info (cpl_func, "                             `ohdmmmmmmmmmmmmdddddhhyyso+//:::----.                      ");
    cpl_msg_info (cpl_func, "                              .yhddmmmmmmmmdddhyyyyssso+///:::----.`                     ");
    cpl_msg_info (cpl_func, "                               :hhdddddddhhhyyyssoooooo++/:::::--....`                   ");
    cpl_msg_info (cpl_func, "                               `ohdhhhhhhyyyysssooooooo++//::::---.--.`                  ");
    cpl_msg_info (cpl_func, "                                .yhhhyyyyyyyyyyssssssoo++//:::::--......`                ");
    cpl_msg_info (cpl_func, "                             `.-:hhdhhyyyyhhyyyyssssoo+++//:::::---.......-+ssoo-        ");
    cpl_msg_info (cpl_func, "                      ``.-/oyhmsmmhdddhyhyyyyysssssooo+++///:::::---......:+dNNNd.       ");
    cpl_msg_info (cpl_func, "              ```.-/+syhmNNNNNydMdhhhdhyyssssssoooooooo+++///::::------:+ymNNNmdy/````   ");
    cpl_msg_info (cpl_func, "     ``.--:/+oyhdmNNNNNNNNNNNNoMMhhhhhhyssssssooosooo+++++////:::---:+ymNNNNmy/-:oyhho/- ");
    cpl_msg_info (cpl_func, "``-/sydmmmmNNNNNNNNNNNNNNNNNNdyMNhyhhhhyssooooooooooo++++++/////::+ymNNNNNh+-:odmmmmmmmm ");
    cpl_msg_info (cpl_func, "sdmmmmNNNNNNNNNNNNNNNNNNNNNNNsdMNyyhhhhhysssooooooso+++++++++//+ydNMMMNms:-odNNNNNNmmmmm ");
    cpl_msg_info (cpl_func, "mNNNNNNNNNNNNNNNNNNNNNNNNNNMNomMNyyhhhhhhyyyssooooooo++++++++ydNMMMMNh+-+hNNNNNNNNNNNmmm ");
    cpl_msg_info (cpl_func, "NNNNNNNNNNNNNNNNNNNNNNNNNNMMM+mMMmhhhhhhhyyysssooooooo+++oydNMMMMNms:/ymNNNNNNNNNNNNNNNm ");
    cpl_msg_info (cpl_func, "NNNNNNNNNNMMMNNNNNNNNNNNNNMMMshMMMmdhhhhhyyyyssoooooooshmNMMMMMNh+:smNNNNNNNNNNNNNNNNNNN ");
    cpl_msg_info (cpl_func, "MMNNNNNNMMMMMMNNNNNNNNNNNMMMMd+MMMMmdhhhhyyyysoooooshmMMMMMMNh+:odNNNNNNNNNNNNNNNNNNNNNN ");
    cpl_msg_info (cpl_func, "MMMMMMMMMMMMMMMMMMMNNNNMMMMMMN:mMMMNmhhdhhhyysooydNMMMMMMNy+:odNMMNNNNNNNNNNNNNNNNNNNNNN ");
    cpl_msg_info (cpl_func, "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMy/MMMMMdhhhyyyyhdNMMMMMNms//sdNMMMNNMNNNNNNNNNNNNNNNNNNNNN ");
    cpl_msg_info (cpl_func, "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMN/oMMMMMmdhhmNMMMMMMNho/+yNMMMMMMMMMNNNNNNNNNNNNNNNNNNNNNN ");
    cpl_msg_info (cpl_func, "                                                  coutesis: The GRAVITY pipeline team    ");
    cpl_msg_info (cpl_func, "*****************************************************************************************");
    return CPL_ERROR_NONE;
}

/**@}*/
