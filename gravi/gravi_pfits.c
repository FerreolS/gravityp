/* $Id: gravi_pfits.c,v 1.12 2011/04/31 06:10:40 nazouaoui Exp $
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
 * @defgroup gravi_pfits  Gravity property list
 *
 * This module implement the function to access the a given property of the header.
 * Some of the are specific for one keywords like @c gravi_pfits_det_starx()
 * (not all of them are documented), and some are more generic like
 * @c gravi_pfits_get_double()
 *
 */
/**@{*/

/*
 * History :
 * 10/01/2019  EKW fix Warning : unused parameter : power
 * 06/11/2019  EkW commit for FE : make use of dedicated keyword introduced Oct 2019
 * 08/11/2019  EkW correct uninitilized parameter (PIPE-8072)
 */
/*-----------------------------------------------------------------------------
                                   Includes
 -----------------------------------------------------------------------------*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cpl.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <erfa.h>

#include "gravi_pfits.h"
#include "gravi_utils.h"

/*-----------------------------------------------------------------------------
                              Function codes
 -----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
 * @brief    find out the name of the propertylist
 * @param    plist       property list to read from
 * @return   The requested value or the pointer to the string
 */
/*----------------------------------------------------------------------------*/

int gravi_pfits_get_startx (const cpl_propertylist * plist)
{
    cpl_errorstate prestate = cpl_errorstate_get();
    int value = cpl_propertylist_get_int(plist, PROFILE_STARTX);
    cpl_ensure (cpl_errorstate_is_equal(prestate), cpl_error_get_code(), 0);
    return value;
}

int gravi_pfits_get_fullstartx (const cpl_propertylist * plist)
{
    cpl_errorstate prestate = cpl_errorstate_get();
    int value = cpl_propertylist_get_int(plist, PROFILE_FULLSTARTX);
    cpl_ensure (cpl_errorstate_is_equal(prestate), cpl_error_get_code(), 0);
    return value;
}

int gravi_pfits_get_nx (const cpl_propertylist * plist)
{
    cpl_errorstate prestate = cpl_errorstate_get();
    int value = cpl_propertylist_get_int(plist, PROFILE_NX);
    cpl_ensure (cpl_errorstate_is_equal(prestate), cpl_error_get_code(), 0);
    return value;
}

int gravi_pfits_get_window_start (const cpl_propertylist * plist)
{
    cpl_errorstate prestate = cpl_errorstate_get();
    int value = cpl_propertylist_get_int(plist, "ESO DET2 FRAM STRX");
    cpl_ensure (cpl_errorstate_is_equal(prestate), cpl_error_get_code(), 0.0);
    return value;
}

double gravi_pfits_get_diameter (const cpl_propertylist * plist, int type_data)
{
    cpl_errorstate prestate = cpl_errorstate_get();
    double value = 0.0;

    const char * keyword = ( type_data == GRAVI_FT ? "ESO FT ROBJ "
            "DIAMETER" : "ESO INS SOBJ DIAMETER" );
    if (cpl_propertylist_has (plist, keyword)) {
        value =  cpl_propertylist_get_double( plist,  keyword );
    }
    else if ( !gravi_pfits_is_calib (plist)) {
        cpl_msg_warning (cpl_func, "The keyword %s does not "
                "exist in the propertylist", keyword);
    }

    cpl_ensure(cpl_errorstate_is_equal(prestate), cpl_error_get_code(), 0.0);
    return value;
}

double gravi_pfits_get_pmra (const cpl_propertylist * plist)
{
    const char * keyword = "ESO FT ROBJ PMA";
    double pma = gravi_pfits_get_double_silentdefault (plist, keyword, 0.);
    return pma; // [as/year]
}

double gravi_pfits_get_pmdec (const cpl_propertylist * plist)
{
    const char * keyword = "ESO FT ROBJ PMD";
    double pmd = gravi_pfits_get_double_silentdefault (plist, keyword, 0.);
    return pmd; // [as/year]
}

double gravi_pfits_get_plx (const cpl_propertylist * plist)
{
    const char * keyword = "ESO FT ROBJ PARALLAX";
    double plx = gravi_pfits_get_double_silentdefault (plist, keyword, 0.);
    return plx; // [as]
}

const char * gravi_pfits_get_extname (const cpl_propertylist * plist)
{
    cpl_errorstate prestate = cpl_errorstate_get();
    const char * value = cpl_propertylist_get_string(plist, "EXTNAME");
    cpl_ensure (cpl_errorstate_is_equal(prestate), cpl_error_get_code(), "");
    return value;
}

const char * gravi_pfits_get_dpr_type (const cpl_propertylist * plist)
{
    const char * value = cpl_propertylist_get_string(plist, "ESO DPR TYPE");
    cpl_ensure (value != NULL, cpl_error_get_code(), NULL);
    return value;
}

const char * gravi_pfits_get_resolution (const cpl_propertylist * plist)
{
    const char * value = cpl_propertylist_get_string(plist, "ESO INS FILT2 NAME");
    cpl_ensure(value != NULL, cpl_error_get_code(), NULL);
    return value;
}

const char * gravi_pfits_get_spec_res (const cpl_propertylist * plist)
{
    const char * value = cpl_propertylist_get_string(plist, "ESO INS SPEC RES");
    cpl_ensure(value != NULL, cpl_error_get_code(), NULL);
    return value;
}

const char * gravi_pfits_get_pola_mode (const cpl_propertylist * plist, int type_data )
{
    const char * keyword = ( type_data == GRAVI_FT ? "ESO FT POLA MODE" : "ESO INS POLA MODE");
    const char * value = cpl_propertylist_get_string(plist, keyword);
    cpl_ensure (value != NULL, cpl_error_get_code(), NULL);
    return value;
}

int gravi_data_frame_get_mode (const cpl_frame * frame)
{
    cpl_ensure (frame, CPL_ERROR_NULL_INPUT, 0);
    const char * value = cpl_frame_get_tag (frame);

    cpl_ensure (value != NULL, cpl_error_get_code(), 0);
    if (strstr(value,"SINGLE")) {
        return MODE_SINGLE;
    }
    if (strstr(value, "DUAL")) {
        return MODE_DUAL;
    }

    return 0;
}

int gravi_pfits_get_mode (const cpl_propertylist * plist)
{
    const char * type;

    if (cpl_propertylist_has (plist, "ESO PRO CATG"))
        type = cpl_propertylist_get_string (plist,"ESO PRO CATG");
    else if (cpl_propertylist_has (plist, "ESO DPR TYPE"))
        type = cpl_propertylist_get_string (plist,"ESO DPR TYPE");
    else return -1;

    if (strstr (type,"SINGLE")) {
        return MODE_SINGLE;
    }
    if (strstr (type, "DUAL")) {
        return MODE_DUAL;
    }

    return 0;
}

int gravi_pfits_get_axis (const cpl_propertylist * plist)
{
    const char * type=NULL;

    /* FE 2019-10-31 now making use of dedicated keyword introduced Oct 2019
       if not in header, use historic method (which is wrong, because the named
       position of onaxis and offaxis are identical, and therefore ICS always
       returns offaxis */

    if (cpl_propertylist_has (plist, "ESO INS ROOF POS")) {
        type = cpl_propertylist_get_string (plist,"ESO INS ROOF POS");
    } else {
      if (cpl_propertylist_has (plist, "ESO INS OPTI11 ID")) {
	    type = cpl_propertylist_get_string (plist,"ESO INS OPTI11 ID");
	    cpl_msg_warning (cpl_func,"Using ESO INS OPTI11 ID (making use of historic method)");
      } else return -1;
    }

    if (type == NULL) return -1;
    if (strstr (type,"OFFAXIS")) {
       return MODE_OFFAXIS;
    }
    if (strstr (type, "ONAXIS")) {
       return MODE_ONAXIS;
    }

    return 0;
}

const char * gravi_pfits_get_mode_name (const cpl_propertylist * plist)
{
    cpl_ensure (plist, CPL_ERROR_NULL_INPUT, NULL);
    const char * type;

    if (cpl_propertylist_has (plist, "ESO PRO CATG"))
        type = cpl_propertylist_get_string (plist,"ESO PRO CATG");
    else if (cpl_propertylist_has (plist, "ESO DPR TYPE"))
        type = cpl_propertylist_get_string (plist,"ESO DPR TYPE");
    else return "UNKNOWN";

    if (strstr (type,"SINGLE")) {
        return "SINGLE";
    }
    if (strstr (type, "DUAL")) {
        return "DUAL";
    }

    return "UNKNOWN";
}

int gravi_pfits_get_pola_num (const cpl_propertylist * plist, int type_data )
{
    const char * keyword = ( type_data == GRAVI_FT ? "ESO FT POLA MODE" : "ESO INS POLA MODE");
    const char * value = cpl_propertylist_get_string(plist, keyword);
    cpl_ensure (value != NULL, cpl_error_get_code(), 0);
    return ( !(strcmp(value,"SPLIT")) ? 2 : 1 );
}

int gravi_pfits_get_extension_type (const cpl_propertylist * plist)
{
    cpl_errorstate prestate = cpl_errorstate_get();

    const char * value = cpl_propertylist_get_string (plist, "XTENSION");
    cpl_ensure (value, CPL_ERROR_ILLEGAL_INPUT, 0);

    int ext_type;

    if (! strcmp(value, "IMAGE"))
        ext_type = 3;
    else if (! strcmp(value, "BINTABLE"))
        ext_type = 2;
    else if (! strcmp(value, "TABLE"))
        ext_type = 2;
    else {
        cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                "Type extension is incorrect");
        ext_type = 0;
    }
    cpl_ensure(cpl_errorstate_is_equal(prestate), cpl_error_get_code(), 0.0);
    return ext_type;
}

double gravi_pfits_get_metfc_lockmjd (const cpl_propertylist * plist, int tel)
{
    cpl_errorstate prestate = cpl_errorstate_get();

    char *lock = cpl_sprintf ("ESO OCS MET LKDT_FC%i", tel+1);
    double value = gravi_convert_to_mjd (cpl_propertylist_get_string (plist, lock));

    if (!cpl_errorstate_is_equal(prestate)) {
        cpl_errorstate_set (prestate);
        cpl_msg_warning (cpl_func, "Cannot read %s", lock);
        value = 0.0;
    }

    cpl_free (lock);
    return value;
}

double gravi_pfits_get_met_wavelength_mean (const cpl_propertylist * plist, cpl_table * met_table)
{
    double met_wavelength_ref = gravi_pfits_get_met_wavelength (plist)*1e-9;
    double lambda_met_mean = met_wavelength_ref +
                cpl_table_get_column_mean(met_table, "LAMBDA_LASER")*1e-9;

    return (lambda_met_mean);
}

double gravi_pfits_get_met_wavelength (const cpl_propertylist * plist)
{
    cpl_errorstate prestate = cpl_errorstate_get();

    /* double power; */
    double wavelength = 0.0;

    /*
     * If cannot read MLC laser power
     */
    if (cpl_propertylist_has (plist, "ESO INS MLC POWER") == FALSE) {
        wavelength = gravi_pfits_get_double(plist, "ESO INS MLC WAVELENG");
        cpl_msg_warning (cpl_func, "Cannot read the laser POWER : use MLC wavelength : %f ", wavelength);
    }
    /*
     * If MLC laser ON
     */
    else if (gravi_pfits_get_double(plist, "ESO INS MLC POWER") != 0) {
        wavelength = gravi_pfits_get_double(plist, "ESO INS MLC WAVELENG");
        cpl_msg_info(cpl_func, "Using laser MLC wavelength : %f ", wavelength);
    }
    /*
     * If MLAS laser ON
     */
    else if (gravi_pfits_get_double(plist, "ESO INS MLAS LPOW") != 0) {
        /*  power=gravi_pfits_get_double(plist, "ESO INS MLAS LPOW"); */
        wavelength = cpl_propertylist_get_double(plist, "ESO INS MLAS LWAV");
        cpl_msg_info(cpl_func, "Using laser MLAS wavelength : %f ", wavelength);
    }

    if (!cpl_errorstate_is_equal(prestate)) {
        cpl_msg_warning (cpl_func, "Cannot read the laser wavelength in the header : %s", cpl_error_get_message());
        cpl_errorstate_set (prestate);
    }

    return wavelength;
}

int gravi_pfits_get_met_mode (const cpl_propertylist * plist)
{
    const char * value = gravi_pfits_get_string_default (plist, "ESO INS MET MODE","ON");

    if ( !strcmp(value, "FAINT") ) {
        return MET_FAINT_HEADER;
    }
    else{
        return MET_NORMAL_HEADER;
    }

}

/*-----------------------------------------------------------------------------*/
/**
 * @brief Extract metrology faint settings
 * @param plist    The input header
 * @return A new cpl_vector with the parameters
 */
/*-----------------------------------------------------------------------------*/

cpl_vector * gravi_pfits_get_met_faint_params (const cpl_propertylist * plist)
{
    const char * faint_keys[]={"ESO INS ANLO3 RATE1",\
                               "ESO INS ANLO3 REPEAT1",\
                               "ESO INS ANLO3 TIMER1",\
                               "ESO INS ANLO3 RATE2",\
                               "ESO INS ANLO3 REPEAT2",\
                               "ESO INS ANLO3 TIMER2"};

    cpl_vector * faint_params = cpl_vector_new(sizeof(faint_keys)/sizeof(faint_keys[0]));

    cpl_vector_fill(faint_params, 0.0);

    for(size_t i=0;i<sizeof(faint_keys)/sizeof(faint_keys[0]);i++)
    {
        cpl_errorstate prestate = cpl_errorstate_get();
        cpl_type par_type = cpl_propertylist_get_type(plist, faint_keys[i]);
        if(par_type==CPL_TYPE_INT)
        {
            cpl_vector_set(faint_params,i,(double) cpl_propertylist_get_int(plist, faint_keys[i]));
        }
        else
        {
            cpl_vector_set(faint_params,i,cpl_propertylist_get_double(plist, faint_keys[i]));
        }
        cpl_ensure(cpl_errorstate_is_equal(prestate), cpl_error_get_code(), faint_params);
    }

    return faint_params;
}

double gravi_pfits_get_geolat (const cpl_propertylist * plist)
{
    cpl_errorstate prestate = cpl_errorstate_get();
    double value = cpl_propertylist_get_double(plist, "ESO ISS GEOLAT");
    cpl_ensure(cpl_errorstate_is_equal(prestate), cpl_error_get_code(), 0.0);
    return value;
}

double gravi_pfits_get_geolon (const cpl_propertylist * plist)
{
    cpl_errorstate prestate = cpl_errorstate_get();
    double value = cpl_propertylist_get_double(plist, "ESO ISS GEOLON");
    cpl_ensure(cpl_errorstate_is_equal(prestate), cpl_error_get_code(), 0.0);
    return value;
}

double gravi_pfits_get_geoelev (const cpl_propertylist * plist)
{
    cpl_errorstate prestate = cpl_errorstate_get();
    double value = cpl_propertylist_get_double (plist, "ESO ISS GEOELEV");
    cpl_ensure(cpl_errorstate_is_equal(prestate), cpl_error_get_code(), 0.0);
    return value;
}

const char * gravi_pfits_get_feed (const cpl_propertylist * plist)
{
    const char * value = gravi_pfits_get_string_default (plist, "ESO ISS FEED MODE", "SINGLE_STS");
    return value;
}

const char * gravi_pfits_get_sobj (const cpl_propertylist * plist)
{
    const char * value = gravi_pfits_get_string_default (plist, "ESO INS SOBJ NAME", "UNKNOWN_SC");
    return value;
}

double gravi_pfits_get_sobj_diam (const cpl_propertylist * plist)
{
    cpl_errorstate prestate = cpl_errorstate_get();
    double value = cpl_propertylist_get_double (plist, "ESO INS SOBJ DIAMETER");
    cpl_ensure(cpl_errorstate_is_equal(prestate), cpl_error_get_code(), 0.0);
    return value;
}

double gravi_pfits_get_sobj_x (const cpl_propertylist * plist)
{
    double value = gravi_pfits_get_double_default (plist, "ESO INS SOBJ X", 0.0);
    return value;
}

double gravi_pfits_get_sobj_y (const cpl_propertylist * plist)
{
    double value = gravi_pfits_get_double_default (plist, "ESO INS SOBJ Y", 0.0);
    return value;
}

double gravi_pfits_get_sobj_offx (const cpl_propertylist * plist)
{
    double value = gravi_pfits_get_double_default (plist, "ESO INS SOBJ OFFX", 0.0);
    return value;
}

double gravi_pfits_get_sobj_offy (const cpl_propertylist * plist)
{
    double value = gravi_pfits_get_double_default (plist, "ESO INS SOBJ OFFY", 0.0);
    return value;
}

double gravi_pfits_get_gvctu_x (const cpl_propertylist * plist)
{
    double value_default = gravi_pfits_get_sobj_x(plist);
    double value1 = gravi_pfits_get_double_default (plist, "ESO FT KAL SEPCTUX1", value_default);
    double value2 = gravi_pfits_get_double_default (plist, "ESO FT KAL SEPCTUX2", value_default);
    double value3 = gravi_pfits_get_double_default (plist, "ESO FT KAL SEPCTUX3", value_default);
    double value4 = gravi_pfits_get_double_default (plist, "ESO FT KAL SEPCTUX4", value_default);
    return (value1+value2+value3+value4)*0.25;
}

double gravi_pfits_get_gvctu_y (const cpl_propertylist * plist)
{
    double value_default = gravi_pfits_get_sobj_y(plist);
    double value1 = gravi_pfits_get_double_default (plist, "ESO FT KAL SEPCTUY1", value_default);
    double value2 = gravi_pfits_get_double_default (plist, "ESO FT KAL SEPCTUY2", value_default);
    double value3 = gravi_pfits_get_double_default (plist, "ESO FT KAL SEPCTUY3", value_default);
    double value4 = gravi_pfits_get_double_default (plist, "ESO FT KAL SEPCTUY4", value_default);
    return (value1+value2+value3+value4)*0.25;
}

const char * gravi_pfits_get_robj (const cpl_propertylist * plist)
{
    const char * value = gravi_pfits_get_string_default (plist, "ESO FT ROBJ NAME", "UNKNOWN_FT");
    return value;
}

double gravi_pfits_get_robj_diam (const cpl_propertylist * plist)
{
    cpl_errorstate prestate = cpl_errorstate_get();
    double value = cpl_propertylist_get_double(plist, "ESO FT ROBJ DIAMETER");
    cpl_ensure(cpl_errorstate_is_equal(prestate), cpl_error_get_code(), 0.0);
    return value;
}


/*-----------------------------------------------------------------------------
                                 TIMING
 -----------------------------------------------------------------------------*/

double gravi_pfits_get_fddlwindow (const cpl_propertylist * plist)
{
    cpl_errorstate prestate = cpl_errorstate_get();
    double value = cpl_propertylist_get_double(plist, "ESO INS FDDL WINDOW");
    cpl_ensure(cpl_errorstate_is_equal(prestate), cpl_error_get_code(), 0.0);
    return value;
}

double gravi_pfits_get_mjd (const cpl_propertylist * plist)
{
    cpl_errorstate prestate = cpl_errorstate_get();
    double value = cpl_propertylist_get_double(plist, "MJD-OBS");
    cpl_ensure (cpl_errorstate_is_equal(prestate), cpl_error_get_code(), 0.0);
    return value;
}

const char * gravi_pfits_get_met_ph (const cpl_propertylist * plist)
{
    const char * tim1_start = cpl_propertylist_get_string(plist, "ESO OCS MET PH_DATE");
    cpl_ensure (tim1_start != NULL, cpl_error_get_code(), NULL);
    return tim1_start;
}

int gravi_pfits_has_gdzero (const cpl_propertylist * plist,
                                       int tel)
{
    cpl_ensure (plist, CPL_ERROR_NULL_INPUT, 0);

    int return_val = 0;
    char name[100];

    sprintf (name, "ESO QC MET GD_ZERO_FC%i", tel);
    if (cpl_propertylist_has(plist, name))
    {
        return_val = 1;
    } else {
        sprintf (name, "ESO PRO MET GD_ZERO_FC%i", tel);
        if (cpl_propertylist_has(plist, name))
        return_val = 1;
    }

    return return_val;
}

double gravi_pfits_get_gdzero (const cpl_propertylist * plist,
                                       int tel)
{
    cpl_ensure (plist, CPL_ERROR_NULL_INPUT, 0);

    double output;
    char name[100];

    sprintf (name, "ESO QC MET GD_ZERO_FC%i", tel);
    if (cpl_propertylist_has(plist, name)) {
        /* Try to read this keyword as a double */
        output = gravi_pfits_get_double(plist, name);
    } else {
        sprintf (name, "ESO PRO MET GD_ZERO_FC%i", tel);
        if (cpl_propertylist_has(plist, name)) {
            /* Try to read this keyword as a double */
            output = gravi_pfits_get_double(plist, name);
        } else
            output = 0;
    }

    return output;
}

int gravi_pfits_has_oplzero (const cpl_propertylist * plist,
                                       int tel)
{
    cpl_ensure (plist, CPL_ERROR_NULL_INPUT, 0);

    int return_val = 0;
    char name[100];

    sprintf (name, "ESO QC MET OPL_ZERO_FC%i", tel);
    if (cpl_propertylist_has(plist, name)) {
        return_val = 1;
    } else {
        sprintf (name,"ESO PRO MET OPL_ZERO_FC%i", tel);
        if (cpl_propertylist_has(plist, name))
        return_val = 1;
    }

    return return_val;
}

double gravi_pfits_get_oplzero (const cpl_propertylist * plist,
                                       int tel)
{
    cpl_ensure (plist, CPL_ERROR_NULL_INPUT, 0);

    double output;
    char name[100];

    sprintf (name, "ESO QC MET OPL_ZERO_FC%i", tel);
    if (cpl_propertylist_has(plist, name)) {
        /* Try to read this keyword as a double */
        output = gravi_pfits_get_double(plist, name);
    } else {
        sprintf (name, "ESO PRO MET OPL_ZERO_FC%i", tel);
        if (cpl_propertylist_has(plist, name)) {
            /* Try to read this keyword as a double */
            output = gravi_pfits_get_double(plist, name);
        } else
            output = 0;
    }

    return output;
}


const char * gravi_pfits_get_start_sc (const cpl_propertylist * plist)
{
    /* SC = TIM1 = DET2 */
    const char * tim1_start = cpl_propertylist_get_string(plist, "ESO INS TIM1 START");
    cpl_ensure (tim1_start != NULL, cpl_error_get_code(), NULL);
    return tim1_start;
}

const char * gravi_pfits_get_start_acqcam (const cpl_propertylist * plist)
{
    /* ACQCAM = TIM2 = DET1 */
    const char * tim1_start = cpl_propertylist_get_string(plist, "ESO INS TIM2 START");
    cpl_ensure (tim1_start != NULL, cpl_error_get_code(), NULL);
    return tim1_start;
}

const char * gravi_pfits_get_start_prcacq (const cpl_propertylist * plist)
{
    /* This is the recording start of RMN tables */
    const char * acq_start = cpl_propertylist_get_string(plist, "ESO PCR ACQ START");
    cpl_ensure (acq_start != NULL, cpl_error_get_code(), NULL);
    return acq_start;
}

double gravi_pfits_get_dit_acqcam (const cpl_propertylist * plist)
{
    /* ACQCAM = TIM2 = DET1 */
    cpl_errorstate prestate = cpl_errorstate_get();
    double value = cpl_propertylist_get_double(plist, "ESO DET1 SEQ1 DIT");
    cpl_ensure(cpl_errorstate_is_equal(prestate), cpl_error_get_code(), 0.0);
    return value;
}

double gravi_pfits_get_dit_sc (const cpl_propertylist * plist)
{
    /* SC = TIM1 = DET2 */
    cpl_errorstate prestate = cpl_errorstate_get();
    double value = cpl_propertylist_get_double(plist, "ESO DET2 SEQ1 DIT");
    cpl_ensure(cpl_errorstate_is_equal(prestate), cpl_error_get_code(), 0.0);
    return value;
}

double gravi_pfits_get_dit_ft (const cpl_propertylist * plist)
{
    /* FT = PCRACQ = DET3 */
    cpl_errorstate prestate = cpl_errorstate_get();
    double value = cpl_propertylist_get_double(plist, "ESO DET3 SEQ1 DIT");
    cpl_ensure(cpl_errorstate_is_equal(prestate), cpl_error_get_code(), 0.0);
    return value;
}

double gravi_pfits_get_period_sc (const cpl_propertylist * plist)
{
    /* SC = TIM1 = DET2 */
    cpl_errorstate prestate = cpl_errorstate_get();
    double value = cpl_propertylist_get_double(plist, "ESO INS TIM1 PERIOD");
    cpl_ensure(cpl_errorstate_is_equal(prestate), cpl_error_get_code(), 0.0);
    return value;
}

double gravi_pfits_get_period_acqcam (const cpl_propertylist * plist)
{
    /* ACQCAM = TIM2 = DET1 */
    cpl_errorstate prestate = cpl_errorstate_get();
    double value = cpl_propertylist_get_double(plist, "ESO INS TIM2 PERIOD");
    cpl_ensure(cpl_errorstate_is_equal(prestate), cpl_error_get_code(), 0.0);
    return value;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Time of the middle of the SC exposure row in [us],
 *        counted from PRC.ACQ.START
 * @param header   The main header
 * @param row      The frame number (0..NDIT-1)
 * @return The time value
 */
/*---------------------------------------------------------------------------*/

double gravi_pfits_get_time_sc (const cpl_propertylist * header, cpl_size row)
{
    gravi_msg_function_start(0);
    cpl_ensure (header, CPL_ERROR_NULL_INPUT, 0.0);
    cpl_errorstate prestate = cpl_errorstate_get();

    /* Time of the middle of the exposure in
     * [us] with respect to PRC.ACQ.START.  */
    double window_time = gravi_pfits_get_fddlwindow (header);
    double period    = gravi_pfits_get_period_sc (header);

    double extra_time_faint = 0.0; 
    if (gravi_pfits_get_met_mode(header)==MET_FAINT_HEADER) {
        extra_time_faint = 1.0;
    }

    double time = 86400 * 1e6 *
        (gravi_convert_to_mjd (gravi_pfits_get_start_sc (header)) -
         gravi_convert_to_mjd (gravi_pfits_get_start_prcacq (header)) ) +
        (period-window_time-extra_time_faint)/2. * 1e6 +
        row * period * 1e6;

    /* Return 0 if error */
    cpl_ensure (cpl_errorstate_is_equal(prestate), cpl_error_get_code(), 0.0);

    gravi_msg_function_exit(0);
    return time;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Time of the middle of the ACQCAM exposure row in [us],
 *        counted from PRC.ACQ.START
 * @param header   The main header
 * @param row      The frame number (0..NDIT-1)
 * @return The time value
 */
/*---------------------------------------------------------------------------*/

double gravi_pfits_get_time_acqcam (const cpl_propertylist * header, cpl_size row)
{
    gravi_msg_function_start(0);
    cpl_ensure (header, CPL_ERROR_NULL_INPUT, 0.0);
    cpl_ensure (row>=0, CPL_ERROR_ILLEGAL_INPUT, 0.0);
    cpl_errorstate prestate = cpl_errorstate_get();

    /* Time of the middle of the exposure in
     * [us] with respect to PRC.ACQ.START.  */
    double period    = gravi_pfits_get_period_acqcam (header);

    double time = 86400 * 1e6 *
        (gravi_convert_to_mjd (gravi_pfits_get_start_acqcam (header)) -
         gravi_convert_to_mjd (gravi_pfits_get_start_prcacq (header)) ) +
        period/2. * 1e6 +
        row * period * 1e6;

    /* Return 0 if error */
    cpl_ensure (cpl_errorstate_is_equal(prestate), cpl_error_get_code(), 0.0);

    gravi_msg_function_exit(0);
    return time;
}

double gravi_pfits_get_ptfc_acqcam (const cpl_propertylist * plist, int spot)
{
    cpl_errorstate prestate = cpl_errorstate_get();
    char name[90];
    sprintf (name, "ESO ACQ PTFC REFPOS%i", spot);
    double value = cpl_propertylist_get_double (plist, name);
    cpl_ensure (cpl_errorstate_is_equal(prestate), cpl_error_get_code(), 0.0);
    return value;
}

double gravi_pfits_get_drotoff (const cpl_propertylist * plist, int tel)
{
    char name[90];

    /* Angle of binary on ACQ image */
    sprintf (name, "ESO INS DROTOFF%i", tel+1);
    double drottoff = cpl_propertylist_get_double (plist, name);

    return drottoff;
}

double gravi_pfits_get_northangle_acqcam (const cpl_propertylist * plist, int tel)
{
    double fangle = 0.0;

    char name[90];
    sprintf (name, "ESO INS DROTOFF%i", tel+1);

    if (cpl_propertylist_has(plist, "ESO ACQ NORTHANG"))
    {
        fangle = cpl_propertylist_get_double (plist, "ESO ACQ NORTHANG");
        cpl_msg_warning( cpl_func, "Using ACQ Camera North angle from header" );
        cpl_msg_warning (cpl_func, "North angle = %.2f [deg] / NorthACQ in Y to X", fangle);
    }
    else
    {
    if (cpl_propertylist_has (plist, "ESO INS SOBJ X") &&
        cpl_propertylist_has (plist, "ESO INS SOBJ Y") &&
        cpl_propertylist_has (plist, name)) {

        /* Angle of binary on ACQ image */
        double drottoff = cpl_propertylist_get_double (plist, name);

        /* Position angle of binary */
        double dx = cpl_propertylist_get_double (plist, "ESO INS SOBJ X");
        double dy = cpl_propertylist_get_double (plist, "ESO INS SOBJ Y");

        /* If in a mapping template, remove dithering offset */
        double posangle = 0.0;
        if (cpl_propertylist_has(plist, "ESO INS SOBJ OFFX") &&
            cpl_propertylist_has(plist, "ESO INS SOBJ OFFY")) {

            dx -= gravi_pfits_get_double_default (plist, "ESO INS SOBJ OFFX", 0.0);
            dy -= gravi_pfits_get_double_default (plist, "ESO INS SOBJ OFFY", 0.0);
        }
        if ((fabs(dx)>0.0) || (fabs(dy)>0.0))
            posangle = atan2 (dx, dy) * CPL_MATH_DEG_RAD;

        /* Angle of North in ACQ image, from vertical (y+) to right (x+) */
        fangle = - posangle - drottoff + 270;
        if (fangle >= 180) fangle -= 360.0;
        if (fangle < -180) fangle += 360.0;
        if (fangle >= 180) fangle -= 360.0;
        if (fangle < -180) fangle += 360.0;

        cpl_msg_info (cpl_func, "Acquisition camera North angle (tel=%i) = %.2f [deg] / NorthACQ in Y to X", (tel+1), fangle);
    }
    else
    {
        cpl_msg_warning (cpl_func, "Cannot compute North angle: fangle = 0.0");
    }
    }

    return fangle;
}


double gravi_pfits_get_zenithangle_beamb_acqcam (const cpl_propertylist * plist, int tel, int n, int nrow)
{
    double zangle = 0.0;

    char pupilrot_start_name[90];
    char pupilrot_end_name[90];
    char drot_enc_start_name[90];
    char drot_enc_end_name[90];
    char kmiroff_name[90];
    char rotoff_name[90];

    sprintf (pupilrot_start_name, "ESO ISS IP%i PUPILROT START", 7-tel*2);
    sprintf (pupilrot_end_name, "ESO ISS IP%i PUPILROT END", 7-tel*2);
    sprintf (drot_enc_start_name, "ESO INS DROT%i ENC START", tel+1);
    sprintf (drot_enc_end_name, "ESO INS DROT%i ENC END", tel+1);
    sprintf (kmiroff_name,  "ESO INS KMIRROROFFSET%i", tel+1);
    sprintf (rotoff_name, "ESO INS ROTATIONOFFSET%i", tel+1);


    if (cpl_propertylist_has (plist, pupilrot_start_name) &&
        cpl_propertylist_has (plist, pupilrot_end_name) &&
        cpl_propertylist_has (plist, drot_enc_start_name) &&
        cpl_propertylist_has (plist, drot_enc_end_name) &&
        cpl_propertylist_has (plist, kmiroff_name) &&
        cpl_propertylist_has (plist, rotoff_name)) {

        /* VLTI Lab image rotation */
        double pupilrot1 = cpl_propertylist_get_double(plist, pupilrot_start_name);
        double pupilrot2 = cpl_propertylist_get_double(plist, pupilrot_end_name);

        double prdif = pupilrot2 - pupilrot1;
        if (prdif > 180)
            prdif -= 360;
        if (prdif < -180)
            prdif += 360;
        double pupilrot = pupilrot1 + prdif*n/(nrow-1);

        /* Derotator encoder values */
        int drotenc1 = cpl_propertylist_get_int(plist, drot_enc_start_name);
        int drotenc2 = cpl_propertylist_get_int(plist, drot_enc_end_name);

        double drotencdif = (double)drotenc2 - (double)drotenc1;
        double posenc = (double)drotenc1 + drotencdif*n/(nrow-1);

        /* K-mirror offset angle */
        double kmirroroffset = cpl_propertylist_get_double(plist, kmiroff_name);

        /* ISS offset angle */
        double rotoff = cpl_propertylist_get_double(plist, rotoff_name);

        zangle = pupilrot + 90. -  posenc / 200. - kmirroroffset - rotoff + 180.;
        if (zangle >= 180) zangle -= 360.0;
        if (zangle < -180) zangle += 360.0;

        cpl_msg_info (cpl_func, "Acquisition camera Zenith angle Beam B (tel=%i) = %.2f [deg] / ZenithACQ in Y to X", (tel+1), zangle);
    }
    else
    {
        cpl_msg_warning (cpl_func, "Cannot compute Zenith angle for Beam B: zangle = 0.0");
    }

    return zangle;
}


const char * gravi_pfits_get_insname (const cpl_propertylist * plist)
{
    const char * value = cpl_propertylist_get_string(plist, "INSNAME");
    cpl_ensure (value != NULL, cpl_error_get_code(), NULL);
    return value;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief FT system gain in [ADU/e]
 * @param plist    The main header
 * @return The gain value
 */
/*---------------------------------------------------------------------------*/

double gravi_pfits_get_ft_gain (const cpl_propertylist * plist)
{
    double gain = 0.0;

    if ( cpl_propertylist_has (plist, "MJD-OBS") &&
         cpl_propertylist_has (plist, "ESO INS DET3 GAIN") &&
         cpl_propertylist_get_double (plist, "MJD-OBS") >= 57584.83 &&
         cpl_propertylist_get_double (plist, "MJD-OBS") <  57617.0) {

        /* Between 57584.83 and 57617.0, the keyword is valid but written
         * in [adu/e] and the value 3.0 shall be replaced by 1.7 [adu/e] */
        gain = cpl_propertylist_get_double (plist, "ESO INS DET3 GAIN");
        if (gain == 3.0) gain = 1.7;
        cpl_msg_warning (cpl_func,"Use FT gain of %.3f [adu/e] (wrong units and value in header)", gain);
    } else if ( cpl_propertylist_has (plist, "MJD-OBS") &&
                cpl_propertylist_has (plist, "ESO INS DET3 GAIN") &&
                cpl_propertylist_get_double (plist, "MJD-OBS") >= 57617.0) {

        /* After 57617.0, the keyword is valid and written in [e/adu] */
        gain = 1. / cpl_propertylist_get_double (plist, "ESO INS DET3 GAIN");
        cpl_msg_info (cpl_func,"Use FT gain of %.3f [adu/e] from header", gain);
    } else {
        /* The keyword is not valid at all */
        gain = 25.0;
        cpl_msg_warning (cpl_func,"Force FT gain to %.3f [adu/e] (wrong value or no value in header)", gain);
    }

    /* gain in ADU/e */
    return gain;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief SC system gain in [ADU/e]
 * @param plist    The main header
 * @return The gain value
 */
/*---------------------------------------------------------------------------*/

double gravi_pfits_get_sc_gain (const cpl_propertylist * plist)
{
    double gain = 0.0;

    if ( cpl_propertylist_has (plist, "MJD-OBS") &&
         cpl_propertylist_has (plist, "ESO INS DET2 GAIN") &&
         cpl_propertylist_get_double (plist, "MJD-OBS") > 57617.0) {

        /* After 57617.0, the keyword is valid and written in [e/adu] */
        gain = 1. / cpl_propertylist_get_double (plist, "ESO INS DET2 GAIN");
        cpl_msg_info (cpl_func,"Use SC gain of %.3f [adu/e] from header", gain);
    } else {
        /* The keyword is not valid at all */
        gain = 0.5;
        cpl_msg_warning (cpl_func,"Force SC gain to %.3f [adu/e]", gain);
    }

    /* gain in ADU/e */
    return gain;
}



/*-----------------------------------------------------------------------------
                    Extract sub-set of propertylist
 -----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------*/
/**
 * @brief Create OIFITS keywords to satisfy standar
 * @param header    The input header
 * @return A new cpl_parameterlist with the keywords
 */
/*-----------------------------------------------------------------------------*/

cpl_propertylist * gravi_plist_get_oifits_keywords (cpl_propertylist * header)
{
    gravi_msg_function_start(1);

    cpl_propertylist * o_plist = cpl_propertylist_new ();

    /* CONTENT */
    cpl_propertylist_update_string (o_plist, "CONTENT", "OIFITS2");
    cpl_propertylist_update_string (o_plist, "REFERENC", "2001PASP..112.1133P");

    /* OBSERVER */
    const char * observer = gravi_pfits_get_string_default (header, "ESO ISS OPER", "Unknown");
    cpl_propertylist_update_string (o_plist, "OBSERVER", observer);

    /* PROG_ID */
    const char * prog_id = gravi_pfits_get_string_default (header, "ESO OBS PROG ID", "Unknown");
    cpl_propertylist_update_string (o_plist, "PROG_ID", prog_id);

    /* PROCSOFT */
    char * procsoft = cpl_sprintf("GRAVITY pipeline %s", PACKAGE_VERSION);
    cpl_propertylist_update_string (o_plist, "PROCSOFT", procsoft);
    FREE (cpl_free, procsoft);

    /* OBJECT -- may be overwritten by ESO esorex mechanism */
    const char * sc_object = gravi_pfits_get_sobj (header);
    const char * ft_object = gravi_pfits_get_robj (header);
    char * object = cpl_sprintf ("%s,%s", sc_object, ft_object);
    cpl_propertylist_update_string (o_plist, "OBJECT", object);
    FREE (cpl_free, object);

    /* INSMODE */
    char * mode;
    mode = cpl_sprintf ("%s,%s,%s,%s",
                        gravi_pfits_get_mode_name (header),
                        gravi_pfits_get_spec_res (header),
                        gravi_pfits_get_pola_mode (header, GRAVI_SC),
                        gravi_pfits_get_pola_mode (header, GRAVI_FT));
    cpl_propertylist_update_string (o_plist, "INSMODE", mode);
    FREE (cpl_free, mode);

    CPLCHECK_NUL ("Cannot fill the OIFITS2 specific keywords");

    gravi_msg_function_exit(1);
    return o_plist;
}

/*-----------------------------------------------------------------------------*/
/**
 * @brief Extract parameters from a product header
 * @param header    The input header
 * @param param     The parameter name
 * @return A new cpl_parameter with the parameter
 */
/*-----------------------------------------------------------------------------*/
cpl_parameter * gravi_pfits_get_extrapixel_param(const cpl_propertylist * header)
{
    /* Check inputs */
    cpl_ensure (header, CPL_ERROR_NULL_INPUT, NULL);

    int key_index = 1;
    const char*  param_value=NULL;
    const char * param_name = "extra-pixel-ft";
    char key[100];

    sprintf (key, "ESO PRO REC1 PARAM%d NAME",key_index);
    while (cpl_propertylist_has(header, key)){
        if (!strcmp(cpl_propertylist_get_string(header, key), param_name))
        {
	    sprintf (key, "ESO PRO REC1 PARAM%d VALUE",key_index);
            param_value = cpl_propertylist_get_string(header, key);
            cpl_msg_info (cpl_func, "Retrieve the extra-pixel-ft option from the p2vm: [%s] ", param_value);
        }
        sprintf (key, "ESO PRO REC1 PARAM%d NAME", ++key_index);
    }

    cpl_parameter * p = cpl_parameter_new_value ("gravity.preproc.extra-pixel-ft", CPL_TYPE_BOOL,
                                 "Include the 6th pixels ot the FT",
                                 "gravity.preproc", TRUE);
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "extra-pixel-ft");

    if (param_value != NULL && !strcmp(param_value, "true"))
        cpl_parameter_set_bool(p, 1);
    else
        cpl_parameter_set_bool(p, 0);

    cpl_msg_info(cpl_func, "Set extra-pixel-ft option : [%d] ", cpl_parameter_get_bool(p));
    return p;
}

/*-----------------------------------------------------------------------------*/
/**
 * @brief Extract QC parameters
 * @param header    The input header
 * @return A new cpl_parameterlist with the QC keywords
 */
/*-----------------------------------------------------------------------------*/

cpl_propertylist *  gravi_plist_get_qc (cpl_propertylist * header)
{
    const char * qc = " QC ";

    /* Check inputs */
    cpl_ensure (header, CPL_ERROR_NULL_INPUT, NULL);

    /* Research all the qc parameter inside the primary header of the data */
    cpl_propertylist * applist = cpl_propertylist_new();
    cpl_size size = cpl_propertylist_get_size (header);
    for (cpl_size i = 0; i < size; i++) {
        cpl_property * p = cpl_propertylist_get (header, i);
        const char * p_name = cpl_property_get_name (p);
        if ( (strstr(p_name, qc) != NULL) ||
             (strstr(p_name, GRAVI_NIGHT_OBS) != NULL) ) {
            cpl_type type_qc = cpl_property_get_type (p);
            switch (type_qc) {
                case CPL_TYPE_FLOAT :
                    if (cpl_property_get_float (p) != NAN)
                        cpl_propertylist_append_property(applist, p);
                    else
                        cpl_msg_warning(cpl_func, "The parameter %s is not correct", cpl_property_get_name (p));
                    break;
                case CPL_TYPE_DOUBLE :
                    if (cpl_property_get_double (p) != NAN)
                        cpl_propertylist_append_property(applist, p);
                    else
                        cpl_msg_warning(cpl_func, "The parameter %s is not correct", cpl_property_get_name (p));
                    break;
                case CPL_TYPE_INT :
                    if (cpl_property_get_int (p) != NAN)
                        cpl_propertylist_append_property(applist, p);
                    else
                        cpl_msg_warning(cpl_func, "The parameter %s is not correct", cpl_property_get_name (p));
                    break;
                case CPL_TYPE_STRING :
                    if (cpl_property_get_string (p))
                        cpl_propertylist_append_property(applist, p);
                    else
                        cpl_msg_warning(cpl_func, "The parameter %s is not correct", cpl_property_get_name (p));
                    break;
                case CPL_TYPE_CHAR :
                    if (cpl_property_get_char (p))
                        cpl_propertylist_append_property(applist, p);
                    else
                        cpl_msg_warning(cpl_func, "The parameter %s is not correct", cpl_property_get_name (p));
                    break;
                case CPL_TYPE_BOOL :
                    cpl_propertylist_append_property (applist, p);
                    break;
                default :
                    cpl_msg_error (cpl_func,"'%s' is an invalid type of property",p_name);
                    cpl_error_set_message(cpl_func, CPL_ERROR_INVALID_TYPE,
                               "invalid type of property");
                    return NULL;
            }
        }
    }
    return applist;
}

/*-----------------------------------------------------------------------------
                    Function to convert string
 -----------------------------------------------------------------------------*/

double gravi_convert_to_mjd (const char * start)
{
    cpl_ensure (start, CPL_ERROR_NULL_INPUT, 0.0);

    /* Cut the string: 2015-04-01T20:36:59.380228
    * and convert each in int or double */
    char * str = cpl_strdup (start);

    str[4] = '\n';
    str[7] = '\n';
    str[10] = '\n';
    str[13] = '\n';
    str[16] = '\n';

    int iy = atoi (str+0); // YYYY
    int im = atoi (str+5); // MM
    int id = atoi (str+8); // DD
    double dhr  = atof (str+11); // hh
    double dmin = atof (str+14); // mm
    double dsec = atof (str+17); // ss.ss
    cpl_free (str);

    /* Get the MJD at 00:00 with ERFA [d] */
    double dmjd0, dmjd;
    eraCal2jd (iy, im, id, &dmjd0, &dmjd);

    /* Return the full MJD [d] */
    return dmjd + (dhr + (dmin + dsec/60.0)/60.0)/24.0;
}

char * gravi_convert_to_timestamp (double mjd)
{
    int year, month, day, hour, minute;
    double fraction, second;

    eraJd2cal (2400000.5, mjd, &year, &month, &day, &fraction);

    hour = (int)(fraction/3600.0);
    fraction -= hour*3600.0;
    minute = (int)(fraction/60.0);
    fraction -= minute*60.0;
    second = fraction;

    return cpl_sprintf ("%04i-%02i-%02iT%02i:%02i:%06.3f", year, month, day, hour, minute, second);
}

double gravi_pfits_get_decep (const cpl_propertylist * plist, double coef)
{
    cpl_ensure (plist, CPL_ERROR_NULL_INPUT, 0.0);
    cpl_ensure (coef >= 0 && coef <= 1, CPL_ERROR_ILLEGAL_INPUT, 0.0);

    /* Get RA and DEC in [rad] */
    // double raz  = gravi_pfits_get_robj_raep (plist);
    double decz = gravi_pfits_get_robj_decep (plist);

    /* Get and convert the offsets from [mas] to [rad] */
    double xi  = gravi_pfits_get_sobj_x (plist) / 3600000. * CPL_MATH_RAD_DEG * coef;
    double eta = gravi_pfits_get_sobj_y (plist) / 3600000. * CPL_MATH_RAD_DEG * coef;

    double sdecz = sin(decz);
    double cdecz = cos(decz);
    double denom = cdecz - eta * sdecz;

    /* Compute new coordinates in [rad] */
    double dec = atan2 (sdecz+eta*cdecz, sqrt(xi*xi + denom*denom));

    /* Return in [rad] */
    return dec;
}

double gravi_ra_to_rad (const char *stri)
{
    cpl_ensure (stri, CPL_ERROR_NULL_INPUT, 0.0);

    /* Assume the format is HH MM SS.SSSS */
    char * str = cpl_strdup (stri);
    str[2] = '\0';
    str[5] = '\0';

    /* Ra in [hours] */
    double out = 0.0;
    out += atof(str+0);
    out += atof(str+3) / 60;
    out += atof(str+6) / 3600;

    cpl_free (str);
    return out / 12 * CPL_MATH_PI;
}

double gravi_dec_to_rad (const char *stri)
{
    cpl_ensure (stri, CPL_ERROR_NULL_INPUT, 0.0);

    /* Assume the format is +DD MM SS.SSSS */
    char * str = cpl_strdup (stri);
    str[3] = '\0';
    str[6] = '\0';

    /* Ra in [hours] */
    double out = 0.0;
    out += atof(str+1);
    out += atof(str+4) / 60;
    out += atof(str+7) / 3600;

    out *= (str[0]=='-' ? -1.0 : 1.0);

    cpl_free (str);
    return out / 180 * CPL_MATH_PI;
}


/*-----------------------------------------------------------------------------
                    Function to read the target coordinates
 -----------------------------------------------------------------------------*/

double gravi_pfits_get_raep(const cpl_propertylist * plist, double coef)
{
    cpl_ensure (plist, CPL_ERROR_NULL_INPUT, 0.0);
    cpl_ensure (coef >= 0 && coef <= 1, CPL_ERROR_ILLEGAL_INPUT, 0.0);

    /* Get RA and DEC in [rad] */
    double raz  = gravi_pfits_get_robj_raep (plist);
    double decz = gravi_pfits_get_robj_decep (plist);

    /* Get and convert the offsets from [mas] to [rad] */
    double xi  = gravi_pfits_get_sobj_x (plist) / 3600000. * CPL_MATH_RAD_DEG * coef;
    double eta = gravi_pfits_get_sobj_y (plist) / 3600000. * CPL_MATH_RAD_DEG * coef;

    double sdecz = sin(decz);
    double cdecz = cos(decz);
    double denom = cdecz - eta * sdecz;

    /* Compute new coordinates in [rad] */
    double ra  = atan2 (xi,denom) + raz;

    /* Make ra within 0-2pi */
    ra = fmod (ra, CPL_MATH_2PI);
    if ( ra < 0.0 ) ra += CPL_MATH_2PI;

    /* Return in [rad] */
    return ra;
}

double gravi_pfits_get_robj_raep(const cpl_propertylist * plist)
{
    cpl_ensure (plist, CPL_ERROR_NULL_INPUT, 99.);
    cpl_errorstate prestate = cpl_errorstate_get();
    char tmp[15];

    /* Read the parameter as string */
    const char * str_value;
    str_value = gravi_pfits_get_string_default (plist, "ESO FT ROBJ ALPHA", "000000.00");

    double value = 99.;

    /* Read of the form '043555.24' */
    cpl_msg_debug( cpl_func, "Found '%s'", str_value );

    /* Read 2 first chars */
    strncpy(tmp, str_value, 2);
    tmp[2] = '\0';
    cpl_msg_debug( cpl_func, "Found tmp '%s' -> %f", tmp, atof(tmp) );
    value = atof(tmp);

    /* Read 2 more chars */
    strncpy(tmp, str_value+2, 4);
    tmp[2] = '\0';
    cpl_msg_debug( cpl_func, "Found tmp '%s' -> %f", tmp, atof(tmp) );
    value += atof(tmp) / 60.0;

    /* Read remaining data */
    strcpy(tmp, str_value+4);
    cpl_msg_debug( cpl_func, "Found tmp '%s' -> %f", tmp, atof(tmp) );
    value += atof(tmp) / 3600.0;

    /* Convert hours to deg */
    value *= 360. / 24;

    /* Verbose */
    cpl_msg_debug( cpl_func, "Convert RA='%s' into RA=%fdeg", str_value, value);

    if (cpl_error_get_code() == CPL_ERROR_DATA_NOT_FOUND) {
        cpl_errorstate_set (prestate);
        cpl_msg_warning(cpl_func, "rarp doesn't exist in this file.");
        value=0;
    }

    /* Check for a change in the CPL error state */
    /* - if it did change then propagate the error and return */
    cpl_ensure(cpl_errorstate_is_equal(prestate), cpl_error_get_code(), 0.0);

    /* Convert in [rad] */
    return value * CPL_MATH_RAD_DEG;
}

double gravi_pfits_get_robj_decep (const cpl_propertylist * plist)
{
    cpl_ensure (plist, CPL_ERROR_NULL_INPUT, 99.);

    cpl_errorstate prestate = cpl_errorstate_get();
    double value = 99.;
    char tmp[15];
    double sign;

    const char * str_value;
    str_value = gravi_pfits_get_string_default (plist, "ESO FT ROBJ DELTA", "+000000.00");

    /* Read of the form '163033.49' -- issue with +/- */
    cpl_msg_debug( cpl_func, "Found '%s'", str_value );

    /* Check the first digit is a sign + or -
       and discard it */
    if ( str_value[0] == '-' ) {
        sign = -1.0;
        str_value = str_value + 1;
    } else if ( str_value[0] == '+' ) {
        sign = +1.0;
        str_value = str_value + 1;
    } else {
        sign = +1.0;
    }

    /* Read 2 first chars */
    strncpy(tmp, str_value, 2);
    tmp[2] = '\0';
    cpl_msg_debug( cpl_func, "Found tmp '%s' -> %f", tmp, atof(tmp) );
    value = atof(tmp);

    /* Read 2 more chars */
    strncpy(tmp, str_value+2, 4);
    tmp[2] = '\0';
    cpl_msg_debug( cpl_func, "Found tmp '%s' -> %f", tmp, atof(tmp) );
    value += atof(tmp) / 60.0;

    /* Read remaining data */
    strcpy(tmp, str_value+4);
    cpl_msg_debug( cpl_func, "Found tmp '%s' -> %f", tmp, atof(tmp) );
    value += atof(tmp) / 3600.0;

    /* Apply sign */
    value *= sign;

    /* Verbose */
    cpl_msg_debug( cpl_func, "Convert DEC='%s' into DEC=%fdeg", str_value, value);

    if (cpl_error_get_code() == CPL_ERROR_DATA_NOT_FOUND) {
        cpl_errorstate_set (prestate);
        cpl_msg_warning(cpl_func, "decep doesn't exist in this file.");

        value=0;
    }

    /* Check for a change in the CPL error state */
    /* - if it did change then propagate the error and return */
    cpl_ensure(cpl_errorstate_is_equal(prestate), cpl_error_get_code(), 0.0);

    /* Convert in [rad] */
    return value * CPL_MATH_RAD_DEG;
}

/*-----------------------------------------------------------------------------
                              Test function codes
 -----------------------------------------------------------------------------*/

int gravi_pfits_is_calib (const cpl_propertylist * plist)
{
    cpl_ensure (plist, CPL_ERROR_NULL_INPUT, 0);
    const char * opt1 = "ESO INS OPTI1 NAME";

    if ( !cpl_propertylist_has (plist, opt1) ) return 0;

    const char * value = cpl_propertylist_get_string (plist, opt1);
    if ( !strcmp (value,"CALIB") ) return 1;

    return 0;
}


/*----------------------------------------------------------------------------
                Generic functions to read and set keyword
  ----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/
/**
 * @brief    Get the double value of the given property list entry.
 * @param    self   A property list.
 * @param    name   The property name to look up.
 * @return   The double value stored in the list entry. The function returns 0 if
 *           an error occurs and an appropriate error code is set.
 *
 * Unlike cpl_propertylist_get_double(), this variant will happily convert from
 * integer to double.
 */
/*----------------------------------------------------------------------------*/

double gravi_pfits_get_double (const cpl_propertylist * self, const char * name)
{
    cpl_ensure (self, CPL_ERROR_NULL_INPUT, 0);
    cpl_ensure (name, CPL_ERROR_NULL_INPUT, 0);
    cpl_ensure (cpl_propertylist_has (self, name), CPL_ERROR_DATA_NOT_FOUND, 0);

    cpl_type type = cpl_propertylist_get_type (self, name);
    switch (type) {
        case CPL_TYPE_CHAR:
        case CPL_TYPE_UCHAR:
        case CPL_TYPE_BOOL:
        case CPL_TYPE_SHORT:
        case CPL_TYPE_USHORT:
        case CPL_TYPE_INT:
        case CPL_TYPE_UINT:
        case CPL_TYPE_LONG:
        case CPL_TYPE_LONG_LONG:
        case CPL_TYPE_ULONG:
        case CPL_TYPE_SIZE:
            return cpl_propertylist_get_long(self, name);
        case CPL_TYPE_FLOAT:
        case CPL_TYPE_DOUBLE:
            return cpl_propertylist_get_double(self, name);
        case CPL_TYPE_STRING:
            cpl_msg_debug (cpl_func,"FITS card %s is string '%s' to double %f",
                           name,cpl_propertylist_get_string(self, name),atof(cpl_propertylist_get_string(self, name)));
            return atof(cpl_propertylist_get_string(self, name));
        default:
            cpl_error_set(cpl_func, CPL_ERROR_TYPE_MISMATCH);
            return 0.;
    }
    return 0.;
}

cpl_error_code gravi_pfits_ensure_double (cpl_propertylist * self, const char * name)
{
    cpl_ensure_code (self, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code (name, CPL_ERROR_NULL_INPUT);

    /* Get value and comment */
    double value = gravi_pfits_get_double (self, name);
    char * comment = cpl_sprintf ("%s", cpl_propertylist_get_comment (self, name));

    /* Delete card and set back to double */
    cpl_propertylist_erase (self, name);
    cpl_propertylist_append_double (self, name, value );
    cpl_propertylist_set_comment (self, name, comment);
    cpl_free (comment);

    return cpl_error_get_code();
}

/*
 * Set the header keyword NAME+EXT to value
 * protected from nan
 */
cpl_error_code gravi_pfits_update_double (cpl_propertylist * plist,
                                          const char * full_name, double value)
{
    cpl_error_code code;
    cpl_ensure_code (plist,     CPL_ERROR_NULL_INPUT);
    cpl_ensure_code (full_name, CPL_ERROR_NULL_INPUT);

    /* Set value with check for nan */
    if ( isnan(value) )
        cpl_propertylist_update_double (plist, full_name, GRAVI_NAN_DOUBLE);
    else
        cpl_propertylist_update_double (plist, full_name, value);

    /* Check */
    if ( (code=cpl_error_get_code()) ) {
        cpl_msg_warning (cpl_func, "Cannot set keyword: %s", full_name);
        return cpl_error_set_message(cpl_func, code, "Cannot set keyword: %s", full_name);
    }

    return CPL_ERROR_NONE;
}

cpl_error_code gravi_pfits_update_int (cpl_propertylist * plist,
                                       const char * full_name, int value)
{
    cpl_error_code code;
    cpl_ensure_code (plist,     CPL_ERROR_NULL_INPUT);
    cpl_ensure_code (full_name, CPL_ERROR_NULL_INPUT);

    /* Set value with check for nan */
    if ( isnan((double)value) )
        cpl_propertylist_update_int (plist, full_name, GRAVI_NAN_INT);
    else
        cpl_propertylist_update_int (plist, full_name, value);

    /* Check */
    if ( (code=cpl_error_get_code()) ) {
        cpl_msg_warning (cpl_func, "Cannot set keyword: %s", full_name);
        return cpl_error_set_message(cpl_func, code, "Cannot set keyword: %s", full_name);
    }

    return CPL_ERROR_NONE;
}

const char * gravi_pfits_get_string_default (const cpl_propertylist * plist,
                                             const char *name,
                                             const char *def)
{
    cpl_ensure (plist, CPL_ERROR_NULL_INPUT, def);
    cpl_ensure (name,  CPL_ERROR_NULL_INPUT, def);

    const char *output;

    if (cpl_propertylist_has(plist, name)) {
        /* Try to read this keyword */
        output = cpl_propertylist_get_string(plist, name);
    } else {
        /* Get the default and add warning only if not CALIB */
        output = def;

        if (!gravi_pfits_is_calib (plist))
            cpl_msg_warning (cpl_func, "Can't find keyword %s (use '%s')", name, def);
        else
            cpl_msg_info (cpl_func, "Can't find keyword %s (use '%s')", name, def);
    }

    return output;
}

double gravi_pfits_get_double_default (const cpl_propertylist * plist,
                                       const char *name,
                                       double def)
{
    cpl_ensure (plist, CPL_ERROR_NULL_INPUT, def);
    cpl_ensure (name,  CPL_ERROR_NULL_INPUT, def);

    double output;

    if (cpl_propertylist_has(plist, name)) {
        /* Try to read this keyword as a double */
        output = gravi_pfits_get_double(plist, name);
    } else {
        /* Get the default and add warning only if not CALIB */
        output = def;

        if (!gravi_pfits_is_calib (plist))
            cpl_msg_warning (cpl_func, "Can't find keyword %s (use '%f')", name, def);
        else
            cpl_msg_info (cpl_func, "Can't find keyword %s (use '%f')", name, def);
    }

    return output;
}

double gravi_pfits_get_double_silentdefault (const cpl_propertylist * plist,
                                             const char *name,
                                             double def)
{
    cpl_ensure (plist, CPL_ERROR_NULL_INPUT, def);
    cpl_ensure (name,  CPL_ERROR_NULL_INPUT, def);
    double output;

    if (cpl_propertylist_has(plist, name))
        output = gravi_pfits_get_double (plist, name);
    else
        output = def;

    return output;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Add a QC.CHECK keyword to the header
 *
 * @param header: the propertylist to update (in-place)
 * @param msg: the string of the message
 *
 * Add a QC.CHECK keyword to the header, and dump a warning
 * message. The function increments the QC.CHECK.FLAGS
 * counter and create a new string keyword QC.CHECK.MSGi
 * This function allows to propagate in HEADER the most critical
 * warning generated by the pipeline.
 */
/*---------------------------------------------------------------------------*/

cpl_error_code gravi_pfits_add_check (cpl_propertylist * header, const char *msg)
{
    gravi_msg_function_start(0);
    cpl_ensure_code (header, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code (msg,    CPL_ERROR_NULL_INPUT);

    /* Increment counter */
    int i = 1;
    if (cpl_propertylist_has (header, "ESO QC CHECK FLAGS"))
        i = cpl_propertylist_get_int (header, "ESO QC CHECK FLAGS") + 1;

    cpl_propertylist_update_int (header, "ESO QC CHECK FLAGS", i);
    /* CPLCHECK_MSG ("Cannot add check flags2..."); */

    /* Set message */
    char qc_name[80];
    sprintf (qc_name, "ESO QC CHECK MSG%i", i);

    cpl_msg_warning (cpl_func, "%s = '%s'", qc_name, msg);
    cpl_propertylist_append_string (header, qc_name, msg);

    CPLCHECK_MSG ("Cannot add check msg...");

    gravi_msg_function_exit(0);
    return CPL_ERROR_NONE;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Add the ESO PRO REC# PIPE LAST_BUILD in header
 *
 * @param header: the propertylist to update (in-place)
 *
 * The header is updated with a string keyword 'ESO PRO REC# PIPE LAST_BUILD'
 * where # is incremented to avoid overwriting the same keyword. The
 * value is set to __DATE__ __TIME__, which are compiler-macro with the last
 * time of full rebuilt.
 */
/*---------------------------------------------------------------------------*/

cpl_error_code gravi_pfits_add_pipe_build (cpl_propertylist * header)
{
    gravi_msg_function_start(0);
    cpl_ensure_code (header, CPL_ERROR_NULL_INPUT);

    char name[100];
    char value[100];

    /* Increment counter until empty slot */
    int i = 0;
    do {
        i++;
        sprintf (name, "ESO PRO REC%i PIPE LAST_BUILD", i);
    } while (cpl_propertylist_has (header, name));

    /* Define the last build string */
    sprintf (value, "%s %s", __DATE__,__TIME__);
    cpl_msg_info (cpl_func, "%s = '%s'", name, value);

    /* Write into header */
    cpl_propertylist_update_string (header, name, value);
    cpl_propertylist_set_comment (header, name, "Last 'make clean all install'");
    CPLCHECK_MSG ("Cannot add PIPE LAST_BUILD...");

    gravi_msg_function_exit(0);
    return CPL_ERROR_NONE;
}


/**@}*/
