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

#ifndef GRAVI_UTILS_H
#define GRAVI_UTILS_H

/*-----------------------------------------------------------------------------
                                   Includes
 -----------------------------------------------------------------------------*/

#include <cpl.h>
#include"gravi_data.h"

/*-----------------------------------------------------------------------------
                                Macros
 -----------------------------------------------------------------------------*/

/* Manipulate bits */
#define gravi_bit_set(number, pos) do{(number) |= (1 << ((int)pos));}while(0)
#define gravi_bit_clear(number, pos)  do{(number) &= (~(1 << ((int)pos)));}while(0)
#define gravi_bit_get(number, pos) (1 & ((number) >> ((int)pos)))


/* Check the CPL errorstate and rise a new error */
#define CPLCHECK(msg) do{int code; if( (code=cpl_error_get_code ()) ) { cpl_msg_error(cpl_func,msg); cpl_error_set_message(cpl_func, code, msg); }} while (0)

/* Check the CPL errorstate and goto exit the cpl error code  */
#define CPLCHECK_MSG(msg) do{int code; if( (code=cpl_error_get_code ()) ) { cpl_msg_error(cpl_func,msg); return cpl_error_set_message(cpl_func, code, msg); }} while (0)

/* Check the CPL errorstate and goto exit as NULL */
#define CPLCHECK_NUL(msg) do{int code; if( (code=cpl_error_get_code ()) ) { cpl_msg_error(cpl_func,msg); cpl_error_set_message(cpl_func, code, msg);  return NULL; }} while (0)

/* Check the CPL errorstate and goto exit the cpl error code as int */
#define CPLCHECK_INT(msg) do{int code; if( (code=cpl_error_get_code ()) ) { cpl_msg_error(cpl_func,msg); cpl_error_set_message(cpl_func, code, msg); return code; }} while (0)

/* Check the CPL errorstate and goto 'clean' */
#define CPLCHECK_CLEAN(msg) do{int code; if( (code=cpl_error_get_code ()) ) { cpl_msg_error(cpl_func,msg);  cpl_error_set_message(cpl_func, code, msg); goto cleanup; }} while (0)

/* Check the CPL errorstate and goto 'tag' */
#define CPLCHECK_GOTO(msg,tag) do{int code; if( (code=cpl_error_get_code ()) ) { cpl_msg_error(cpl_func,msg);  cpl_error_set_message(cpl_func, code, msg); goto tag; }} while (0)

/* Check the CPL errorstate and goto 'clean' */
#define ERROR_CLEAN(code,msg) do{cpl_msg_error(cpl_func,msg);  cpl_error_set_message(cpl_func, code, msg); goto cleanup;} while (0)

/* Check the CPL errorstate and goto exit as NULL */
#define CHECK_NUL(flag,msg) do{ if( flag ) { cpl_msg_error(cpl_func,msg); cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT, msg);  return NULL; }} while (0)

/* Check the CPL errorstate and goto exit the cpl error code  */
#define CHECK_MSG(flag,msg) do{ if( flag ) { cpl_msg_error(cpl_func,msg); return cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT, msg); }} while (0)

/* Protected memory free */
#define FREE(function,variable) do{ if (variable) { function(variable); variable=NULL;}  } while (0)

/* Protected memory free of an array of variable */
#define FREELOOP(function,variable,n) do{ if (variable) { for (int iloopfree = 0 ; iloopfree < n; iloopfree++) { if (variable[iloopfree]) { function( variable[iloopfree] ); variable[iloopfree] = NULL; } } if (variable) { cpl_free (variable); variable = NULL; } } } while (0)

/* Get a the short filename */
#define FILESHORT(file) (strrchr(file, '/') ? strrchr(file, '/') + 1 : file)

#define TEST_MESSAGE  cpl_msg_info (cpl_func,">>>>>>> TEST <<<<<<<")


/* Much faster than pow (data, 2) on my compilator */
#define gravi_pow2(data) data*data

/* Timer function to enter and exit function */
# define gravi_msg_function_start(flag) clock_t timer_function = clock();  do{ if(flag) cpl_msg_info(cpl_func,"Start function %s", cpl_func); }while(0)
# define gravi_msg_function_exit(flag) do{ if(flag) cpl_msg_info(cpl_func,"Exit function %s (%.6f s)",cpl_func,(double)(clock()-timer_function)/(double)CLOCKS_PER_SEC); }while(0)

/* Test shutters all open or all closed */
# define gravi_data_get_shutter(data,tel) gravi_get_shutter (gravi_data_get_header(data),tel)
# define gravi_data_check_shutter(data,t0,t1,t2,t3) gravi_check_shutter (gravi_data_get_header(data), t0,t1,t2,t3)
# define gravi_data_check_shutter_closed(data) gravi_check_shutter (gravi_data_get_header(data), 0,0,0,0)
# define gravi_data_check_shutter_open(data) gravi_check_shutter (gravi_data_get_header(data), 1,1,1,1)
# define gravi_spectrum_get_npol(table) (gravi_spectrum_get_nregion (table) > 24 ? 2 : 1)


/*-----------------------------------------------------------------------------
                        Defines and Static variables
 -----------------------------------------------------------------------------*/

/* 
 * Baseline and telescope order.
 * Also define closing triangles of each baseline
 * and closure phases 
 */

extern int GRAVI_BASE_TEL[6][2];
extern char GRAVI_BASE_NAME[6][3];
extern int GRAVI_TRI_BASE[6][2][2];
extern int GRAVI_TRI_SIGN[6][2][2];
extern int GRAVI_CLO_BASE[4][3];
extern int GRAVI_CLO_TEL[4][3];
extern char GRAVI_CLO_NAME[4][4];
extern char GRAVI_DATA[50][7];
extern char GRAVI_DATAERR[50][10];

/*
 * Lab input connected to each GRAVITY input window
 */
#define GRAVI_LABINPUT_1 7
#define GRAVI_LABINPUT_2 5
#define GRAVI_LABINPUT_3 3
#define GRAVI_LABINPUT_4 1

#define SHUTTER_KEY "IPAG INS SHUT" //"GRAVITY SHUT"
#define SHUTTER_KEY1 "IPAG INS SHUT1 ST" //"GRAVITY SHUT 1"
#define SHUTTER_KEY2 "IPAG INS SHUT2 ST" //"GRAVITY SHUT 2"
#define SHUTTER_KEY3 "IPAG INS SHUT3 ST" //"GRAVITY SHUT 3"
#define SHUTTER_KEY4 "IPAG INS SHUT4 ST" //"GRAVITY SHUT 4"
#define GRAVI_SHUTTER_KEY "ESO INS SHUT"
#define GRAVI_SHUTTER_KEY1 "ESO INS SHUT 1"
#define GRAVI_SHUTTER_KEY2 "ESO INS SHUT 2"
#define GRAVI_SHUTTER_KEY3 "ESO INS SHUT 3"
#define GRAVI_SHUTTER_KEY4 "ESO INS SHUT 4"
#define POLAR_1 "S"
#define POLAR_2 "P"
#define POLAR_3 "C"
#define GRAVI_POLAR(pol,npol) (npol==1 ? "C" : (pol==0 ? "S" : "P") )
#define PHASE_1 "A"
#define PHASE_2 "B"
#define PHASE_3 "C"
#define PHASE_4 "D"
#define COHERENCE "COHERENCE"
#define TRANSMISSION "TRANSMISSION"
#define PHASE "PHASE"

#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define IMIN(a,b) (((a) < (b)) ? (a) : (b))

/*-----------------------------------------------------------------------------
                              Public prototypes
 -----------------------------------------------------------------------------*/

const char * gravi_get_license (void) ;
cpl_error_code gravi_msg_fixme (const char * msg);
cpl_error_code gravi_dump_the_boss (double ra, double dec);

int gravi_get_shutter (cpl_propertylist *, int);
int gravi_check_shutter (cpl_propertylist *, int, int, int, int);
int gravi_data_check_shutter_beam (gravi_data ** datas, int nb_datas);
int gravi_get_shutter_id (cpl_propertylist * header);
int gravi_get_shutter_baseid (cpl_propertylist * header);

int gravi_region_get_base(cpl_table *imaging_detector, int region);
int gravi_region_get_pol (cpl_table *imaging_detector, int region);
int gravi_get_region (cpl_table *img_det, int base, char phase, int pol);
int gravi_region_get_tel (cpl_table *imaging_detector, int region, int beam);
int gravi_region_get_phaseid (cpl_table *imaging_detector, int region);
char gravi_region_get_phase (cpl_table *imaging_detector, int region);
int gravi_region_get_base_sign (cpl_table *imaging_detector, int base);

int gravi_wave_get_nlambda(cpl_table *wave_data, double lambda_min, double lambda_max);

cpl_table * gravi_table_oi_create (int , int , const char * );

int * gravi_image_extract_dimension (cpl_image * );

short gravi_sta_index (int gravi_input, cpl_table * optical_train_table, cpl_table * array_geometry_table);

cpl_size gravi_spectrum_get_nregion (const cpl_table * table);
cpl_size gravi_spectrum_get_nwave (const cpl_table * table);

double gravi_spectrum_get_flux (const cpl_table * table);
double gravi_imagelist_get_flux (const cpl_imagelist * imglist);

cpl_error_code gravi_lkdt_get_sequence (cpl_table * oi_table,
                                        cpl_size ntel,
                                        cpl_size *first,
                                        cpl_size *nobs);

cpl_vector * gravi_compute_envelope (const cpl_vector * opd, int wave, int n_wave);


#endif
