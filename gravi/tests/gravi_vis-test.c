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
 *
 *  History 
 *  ekw  14/11/2018 correct unused variable / forward declaration of gravi_average_self_visphi / cmin-cmax pointer
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
#include <cpl_table.h>
#include <cpl_array.h>
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

#define STR(x) x

int gravi_vis_test(void);

cpl_error_code gravi_average_self_visphi(cpl_table * oi_vis_avg, cpl_table * oi_vis,
               cpl_array * wavenumber, const char * phase_ref, int* cmin, int* cmax, int nrange);

int gravi_vis_test(void){

	int flag = EXIT_SUCCESS;
	int nrange=0;
	int cmin, cmax;
	char filename[128];

	sprintf(filename, "%s%s", STR(DATADIR), "/oi_vis_SC.fits");
	cpl_table* oi_vis_SC = cpl_table_load 	(filename,1,0 );

	sprintf(filename, "%s%s", STR(DATADIR), "/vis_SC.fits");
	cpl_table* vis_SC = cpl_table_load 	(filename,1,0);

	sprintf(filename, "%s%s", STR(DATADIR), "/wavenumber_sc_table.fits");
	cpl_table* wavenumber_sc_tab = cpl_table_load (filename,1,0);

	double * waven_double = cpl_table_get_data_double(wavenumber_sc_tab, "TEST");
	cpl_array *wavenumber_sc = cpl_array_wrap_double(waven_double,14);

	flag = gravi_average_self_visphi(oi_vis_SC, vis_SC, wavenumber_sc, "SELF_REF", &cmin, &cmax, nrange);
	
        FREE(cpl_table_delete,  oi_vis_SC);
	
	FREE(cpl_table_delete,  vis_SC);
	
	cpl_array_unwrap(wavenumber_sc);
	FREE(cpl_table_delete,  wavenumber_sc_tab );

	
 return flag;
}
/*----------------------------------------------------------------------------*/
/**
  @brief    Unit tests of gravi_vis module
 */
/*----------------------------------------------------------------------------*/

int main(void)
{
    int flag;

    cpl_test_init(PACKAGE_BUGREPORT, CPL_MSG_INFO);

    flag=gravi_vis_test();

    if (flag == EXIT_FAILURE)
    {
    	cpl_test_end(0);
    	exit(EXIT_FAILURE);
    }

    cpl_test_end(0);
    exit(EXIT_SUCCESS);
}
