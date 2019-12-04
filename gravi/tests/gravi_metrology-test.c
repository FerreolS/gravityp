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
 *  History :
 *  14/11/2018  correct unused variable warnings
 *  12/11/2018  add gravi_data *static_param_data
 *  15/01/2019  unused parameter telname
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

/* declaration of functions to test */
double gravi_metrology_get_fc_focus (cpl_propertylist * header, int gv, gravi_data* gravi_data_focus);
double gravi_metrology_get_fc_shift (cpl_propertylist * header, int gv, gravi_data* gravi_data_focus);
cpl_error_code gravi_metrology_get_astig (cpl_propertylist * header, int gv, double * amplitude, double * angle, double * radius);

int gravi_metrology_test(void);

#define STR(x) x

int gravi_metrology_test(void){
	cpl_msg_info (cpl_func,"EKKI ENTER gravity_metrology-test\n");
    int flag = EXIT_SUCCESS;

    /* Prepare gravi_data *static_param_data with the default focus data */
    char filename[128];
    sprintf(filename, "%s%s", STR(DATADIR), "/GRAVI_STATIC_CALIB.fits");
    //sprintf(filename, "%s","/home/grav/pipeline_ekki/execute/GRAVI_STATIC_CALIB.fits");
    gravi_data * gravi_data_focus = gravi_data_load_ext(filename,"FOCUSPAR");

    /* Treat header keywords */
    cpl_propertylist *propertylist = cpl_propertylist_new();
    cpl_propertylist_prepend_int (propertylist,"ESO ISS CONF INPUT1", 1);
    cpl_propertylist_prepend_int (propertylist,"ESO ISS CONF INPUT2", 3);
    cpl_propertylist_prepend_int (propertylist,"ESO ISS CONF INPUT3", 5);
    cpl_propertylist_prepend_int (propertylist,"ESO ISS CONF INPUT4", 7);

    cpl_propertylist_prepend_string (propertylist,"ESO ISS CONF T1NAME", "AT1");
    cpl_propertylist_prepend_string (propertylist,"ESO ISS CONF T2NAME", "AT2");
    cpl_propertylist_prepend_string (propertylist,"ESO ISS CONF T3NAME", "AT3");
    cpl_propertylist_prepend_string (propertylist,"ESO ISS CONF T4NAME", "AT4");
    /* -------------------------------------------------------------------------------------------------------- */

    cpl_msg_info (cpl_func,"gravi_metrology_get_fc_focus \n---------------------------------------------------\n");
    double defocus[4];

    defocus[0] = gravi_metrology_get_fc_focus (propertylist, 0, gravi_data_focus);
    defocus[1] = gravi_metrology_get_fc_focus (propertylist, 1,  gravi_data_focus);
    defocus[2] = gravi_metrology_get_fc_focus (propertylist, 2,  gravi_data_focus);
    defocus[3] = gravi_metrology_get_fc_focus (propertylist, 3,  gravi_data_focus);
    cpl_msg_info (cpl_func,"Defocus : %e %e %e %e \n", defocus[0], defocus[1], defocus[2], defocus[3]);

    cpl_propertylist_prepend_float (propertylist,"ESO MET GV1 AT FC FOCUS", -70.0);
    cpl_propertylist_prepend_float (propertylist,"ESO MET GV2 AT FC FOCUS", -90.0);
    cpl_propertylist_prepend_float (propertylist,"ESO MET GV3 AT FC FOCUS",  30.0);
    cpl_propertylist_prepend_float (propertylist,"ESO MET GV4 AT FC FOCUS", -70.0);

    /* const char * telname = gravi_conf_get_telname(1, propertylist ); */

    defocus[0] = gravi_metrology_get_fc_focus (propertylist, 0,  gravi_data_focus);
    defocus[1] = gravi_metrology_get_fc_focus (propertylist, 1,  gravi_data_focus);
    defocus[2] = gravi_metrology_get_fc_focus (propertylist, 2,  gravi_data_focus);
    defocus[3] = gravi_metrology_get_fc_focus (propertylist, 3,  gravi_data_focus);

    cpl_msg_info (cpl_func,"Defocus with header : %e %e %e %e\n", defocus[0], defocus[1], defocus[2], defocus[3]);

    /* -------------------------------------------------------------------------------------------------------- */
    cpl_msg_info (cpl_func,"Test : gravi_metrology_get_fc_shift \n-------------------------------------------------------\n");

    /* test without header keywords for default values */
    double shift[4];

    shift[0] = gravi_metrology_get_fc_shift (propertylist, 0,  gravi_data_focus);
    shift[1] = gravi_metrology_get_fc_shift (propertylist, 1,  gravi_data_focus);
    shift[2] = gravi_metrology_get_fc_shift (propertylist, 2,  gravi_data_focus);
    shift[3] = gravi_metrology_get_fc_shift (propertylist, 3,  gravi_data_focus);

    cpl_msg_info (cpl_func,"Shift  : %e  %e %e %e \n", shift[0], shift[1], shift[2], shift[3]);

    /* Now add header keywords with different values */
    cpl_propertylist_prepend_float (propertylist,"ESO MET GV1 UT FC SHIFT", -450.0);
    cpl_propertylist_prepend_float (propertylist,"ESO MET GV2 UT FC SHIFT", -350.0);
    cpl_propertylist_prepend_float (propertylist,"ESO MET GV3 UT FC SHIFT", -50.0);
    cpl_propertylist_prepend_float (propertylist,"ESO MET GV4 UT FC SHIFT", -525.0);

    cpl_propertylist_prepend_float (propertylist,"ESO MET GV1 AT FC SHIFT",  -25.0 *1.8/8.0);
    cpl_propertylist_prepend_float (propertylist,"ESO MET GV2 AT FC SHIFT",  -30.0 *1.8/8.0);
    cpl_propertylist_prepend_float (propertylist,"ESO MET GV3 AT FC SHIFT",   45.0 *1.8/8.0);
    cpl_propertylist_prepend_float (propertylist,"ESO MET GV4 AT FC SHIFT",  -50.0 *1.8/8.0);

    shift[0] = gravi_metrology_get_fc_shift (propertylist, 0,  gravi_data_focus);
    shift[1] = gravi_metrology_get_fc_shift (propertylist, 1,  gravi_data_focus);
    shift[2] = gravi_metrology_get_fc_shift (propertylist, 2,  gravi_data_focus);
    shift[3] = gravi_metrology_get_fc_shift (propertylist, 3,  gravi_data_focus);

    cpl_msg_info (cpl_func,"Shift with header %e %e %e %e \n", shift[0],shift[1], shift[2],shift[3]);

    /* -------------------------------------------------------------------------------------------------------- */
    cpl_msg_info (cpl_func,"Test : gravi_metrology_get_astig \n-------------------------------------------------------\n");

    /* test without header keywords for default values */
    double AstigmAmplitude; // [nm]
    double AstigmTheta; // [rad]
    double AstigmRadius; // [m]
    cpl_error_code err;

    err = gravi_metrology_get_astig (propertylist, 0, &AstigmAmplitude, &AstigmTheta, &AstigmRadius);
    cpl_msg_info (cpl_func,"0: AstigmAmplitude %e  - AstigmTheta %e  - AstigmRadius %e \n",AstigmAmplitude,AstigmTheta, AstigmRadius);
    err = gravi_metrology_get_astig (propertylist, 1, &AstigmAmplitude, &AstigmTheta, &AstigmRadius);
    cpl_msg_info (cpl_func,"1: AstigmAmplitude %e  - AstigmTheta %e  - AstigmRadius %e \n",AstigmAmplitude,AstigmTheta, AstigmRadius);
    err = gravi_metrology_get_astig (propertylist, 2, &AstigmAmplitude, &AstigmTheta, &AstigmRadius);
    cpl_msg_info (cpl_func,"2: AstigmAmplitude %e  - AstigmTheta %e  - AstigmRadius %e \n",AstigmAmplitude,AstigmTheta, AstigmRadius);
    err = gravi_metrology_get_astig (propertylist, 3, &AstigmAmplitude, &AstigmTheta, &AstigmRadius);
    cpl_msg_info (cpl_func,"3: AstigmAmplitude %e  - AstigmTheta %e  - AstigmRadius %e \n",AstigmAmplitude,AstigmTheta, AstigmRadius);

    cpl_propertylist_prepend_float (propertylist,"ESO MET GV1 UT ASTIG AMP", 182.0);
    cpl_propertylist_prepend_float (propertylist,"ESO MET GV1 UT ASTIG ANG",  -2.0);
    cpl_propertylist_prepend_float (propertylist,"ESO MET GV2 UT ASTIG AMP", 185.0);
    cpl_propertylist_prepend_float (propertylist,"ESO MET GV2 UT ASTIG ANG",  18.0);
    cpl_propertylist_prepend_float (propertylist,"ESO MET GV3 UT ASTIG AMP", 113.0);
    cpl_propertylist_prepend_float (propertylist,"ESO MET GV3 UT ASTIG ANG",  20.0);
    cpl_propertylist_prepend_float (propertylist,"ESO MET GV4 UT ASTIG AMP", 242.0);
    cpl_propertylist_prepend_float (propertylist,"ESO MET GV4 UT ASTIG ANG",  19.0);

    cpl_propertylist_prepend_float (propertylist,"ESO MET GV1 AT ASTIG AMP", 164.0);
    cpl_propertylist_prepend_float (propertylist,"ESO MET GV1 AT ASTIG ANG",   1.0);
    cpl_propertylist_prepend_float (propertylist,"ESO MET GV2 AT ASTIG AMP", 166.0);
    cpl_propertylist_prepend_float (propertylist,"ESO MET GV2 AT ASTIG ANG",  28.0);
    cpl_propertylist_prepend_float (propertylist,"ESO MET GV3 AT ASTIG AMP",  99.0);
    cpl_propertylist_prepend_float (propertylist,"ESO MET GV3 AT ASTIG ANG",   0.4);
    cpl_propertylist_prepend_float (propertylist,"ESO MET GV4 AT ASTIG AMP", 266.0);
    cpl_propertylist_prepend_float (propertylist,"ESO MET GV4 AT ASTIG ANG",  25.0);

    gravi_metrology_get_astig (propertylist, 0, &AstigmAmplitude, &AstigmTheta, &AstigmRadius);
    printf ("0 with header AstigmAmplitude %e  - AstigmTheta %e  - AstigmRadius %e \n",AstigmAmplitude,AstigmTheta, AstigmRadius);
    gravi_metrology_get_astig (propertylist, 1, &AstigmAmplitude, &AstigmTheta, &AstigmRadius);
    printf ("1  with header AstigmAmplitude %e  - AstigmTheta %e  - AstigmRadius %e \n",AstigmAmplitude,AstigmTheta, AstigmRadius);
    gravi_metrology_get_astig (propertylist, 2, &AstigmAmplitude, &AstigmTheta, &AstigmRadius);
    printf ("2  with header AstigmAmplitude %e  - AstigmTheta %e  - AstigmRadius %e \n",AstigmAmplitude,AstigmTheta, AstigmRadius);
    gravi_metrology_get_astig (propertylist, 3, &AstigmAmplitude, &AstigmTheta, &AstigmRadius);
    printf ("3  with header AstigmAmplitude %e  - AstigmTheta %e  - AstigmRadius %e \n",AstigmAmplitude,AstigmTheta, AstigmRadius);


    /* gravi_pfits_get_fangle_acqcam */

    FREE(cpl_propertylist_delete, propertylist  );
    gravi_data_delete(gravi_data_focus);

 return flag;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Unit tests of gravi_metrology module
 */
/*----------------------------------------------------------------------------*/

int main(void)
{
	int flag;

    cpl_test_init(PACKAGE_BUGREPORT, CPL_MSG_INFO);

    flag=gravi_metrology_test();

    if (flag == EXIT_FAILURE)
    {
    	cpl_test_end(0);
    	exit(EXIT_FAILURE);
    }

    cpl_test_end(0);
    exit(EXIT_SUCCESS);
}
