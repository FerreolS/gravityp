/* $Id: gravi_dfs.h,v 1.9 2011/04/31 06:10:40 nazouaoui Exp $
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

/*
 * History
 *     12/11/2018 add cpl_frameset * gravi_frameset_extract_static_param
 *                   add STATIC_PARAM
 *     04/12/2018 add cpl_frameset * gravi_frameset_extract_wave_param (cpl_frameset * frameset);
 *                   add WAVE_PARAM
 */
#ifndef GRAVI_DFS_H
#define GRAVI_DFS_H

/*-----------------------------------------------------------------------------
                                   Define
 -----------------------------------------------------------------------------*/

/* Flag for recipe helps */
#define GRAVI_UNOFFERED      "*UNOFFERED* "
#define GRAVI_RECIPE_FLOW    "\nReduction steps:\n---------------------"
#define GRAVI_RECIPE_INPUT   "\nDO.CATG in input SoF:\n---------------------"
#define GRAVI_RECIPE_OUTPUT  "\nPRO.CATG of products:\n---------------------"


/* DO.CATG for RAW data */
#define GRAVI_PIEZOTF_RAW                   "PIEZOTF_RAW"
#define GRAVI_P2VM_RAW                      "P2VM_RAW"
#define GRAVI_DISP_RAW                      "DISP_RAW"
#define GRAVI_DARK_RAW                      "DARK_RAW"
#define GRAVI_WAVE_RAW                      "WAVE_RAW"
#define GRAVI_WAVESC_RAW                    "WAVESC_RAW"
#define GRAVI_WAVELAMP_RAW                  "WAVELAMP_RAW"
#define GRAVI_FLAT_RAW                      "FLAT_RAW"
#define GRAVI_SINGLE_CALIB_RAW              "SINGLE_CAL_RAW"
#define GRAVI_SINGLE_SCIENCE_RAW            "SINGLE_SCI_RAW"
#define GRAVI_DUAL_CALIB_RAW                "DUAL_CAL_RAW"
#define GRAVI_DUAL_SCIENCE_RAW              "DUAL_SCI_RAW"
#define GRAVI_DUAL_SKY_RAW                  "DUAL_SKY_RAW"
#define GRAVI_SINGLE_SKY_RAW                "SINGLE_SKY_RAW"

/* PRO.CATG / DO.CATG for intermediate product */
#define GRAVI_DISP_VIS                      "DISP_VIS"
#define GRAVI_PREPROC                       "PREPROC"
#define GRAVI_SPECTRUM                      "SPECTRUM"
#define GRAVI_P2VMRED_SINGLE_CALIB          "SINGLE_CAL_P2VMRED"
#define GRAVI_P2VMRED_SINGLE_SCIENCE        "SINGLE_SCI_P2VMRED"
#define GRAVI_P2VMRED_DUAL_CALIB            "DUAL_CAL_P2VMRED"
#define GRAVI_P2VMRED_DUAL_SCIENCE          "DUAL_SCI_P2VMRED"
#define GRAVI_ASTROREDUCED                  "ASTROREDUCED"
#define GRAVI_KEY_PATCH                     "KEY_PATCH"
#define GRAVI_STATIC_PARAM                  "STATIC_PARAM"
#define GRAVI_WAVE_PARAM                    "WAVE_PARAM"

/* PRO.CATG / DO.CATG for calibration product */
#define GRAVI_BAD_MAP                       "BAD"
#define GRAVI_FLAT_MAP                      "FLAT"
#define GRAVI_WAVE_MAP                      "WAVE"
#define GRAVI_P2VM_MAP                      "P2VM"
#define GRAVI_DARK_MAP                      "DARK"
#define GRAVI_EOP_MAP                       "EOP_PARAM"
#define GRAVI_DISP_MODEL                    "DISP_MODEL"
#define GRAVI_DIAMETER_CAT                  "DIAMETER_CAT"
#define GRAVI_DIODE_POSITION		        "DIODE_POSITION"
#define GRAVI_NAB_CAL                       "NAB_CAL"
#define GRAVI_DEBUG_MAP                     "DEBUG"
#define GRAVI_BIASMASK_MAP                  "BIASMASK"
#define GRAVI_WAVELAMP_MAP                  "WAVELAMP"
#define GRAVI_PIEZOTF_MAP                   "PIEZOTF"

#define GRAVI_SINGLE_SKY_MAP                "SINGLE_SKY"
#define GRAVI_DUAL_SKY_MAP                  "DUAL_SKY"

#define GRAVI_FLAT_ACQ_MAP                  "FLAT_ACQ"
#define GRAVI_BAD_ACQ_MAP                   "BAD_ACQ"

/* PRO.CATG / DO.CATG for Visibility product */
#define GRAVI_VIS_SINGLE_SCIENCE            "SINGLE_SCI_VIS"
#define GRAVI_VIS_SINGLE_CALIB              "SINGLE_CAL_VIS"
#define GRAVI_VIS_DUAL_SCIENCE              "DUAL_SCI_VIS"
#define GRAVI_VIS_DUAL_CALIB                "DUAL_CAL_VIS"

#define GRAVI_VIS_SINGLE_CALIBRATED         "SINGLE_SCI_VIS_CALIBRATED"
#define GRAVI_VIS_DUAL_CALIBRATED           "DUAL_SCI_VIS_CALIBRATED"
#define GRAVI_VIS_CALIBRATED(data_mode) (data_mode==MODE_DUAL?GRAVI_VIS_DUAL_CALIBRATED:GRAVI_VIS_SINGLE_CALIBRATED)

/* PRO.CATG / DO.CATG for Transfer Function product */
#define GRAVI_TF_SINGLE_CALIB               "SINGLE_CAL_TF"
#define GRAVI_TF_SINGLE_SCIENCE             "SINGLE_SCI_TF"
#define GRAVI_TF_SCIENCE(data_mode) (data_mode==MODE_DUAL?GRAVI_TF_DUAL_SCIENCE:GRAVI_TF_SINGLE_SCIENCE)

#define GRAVI_TF_DUAL_SCIENCE               "DUAL_SCI_TF"
#define GRAVI_TF_DUAL_CALIB                 "DUAL_CAL_TF"
#define GRAVI_TF_CALIB(data_mode) (data_mode==MODE_DUAL?GRAVI_TF_DUAL_CALIB:GRAVI_TF_SINGLE_CALIB)

/* Still unsuported data */
#define GRAVI_ZP_CAL                        "ZP_CAL"
#define GRAVI_MIRA_INPUT_PROCATG            "VIS_CALIBRATED"
#define GRAVI_MIRA_OUTPUT_PROCATG           "IMAGE"

/*-----------------------------------------------------------------------------
                                Public prototypes
 -----------------------------------------------------------------------------*/

void gravity_print_banner (void);

cpl_error_code gravi_dfs_set_groups(cpl_frameset *);

cpl_error_code gravi_parameter_disable (cpl_parameter * p);

cpl_parameter * gravi_parameter_add_badpix (cpl_parameterlist *self);
cpl_parameter * gravi_parameter_add_profile (cpl_parameterlist *self);
cpl_parameter * gravi_parameter_add_preproc (cpl_parameterlist *self);
cpl_parameter * gravi_parameter_add_wave (cpl_parameterlist *self);
cpl_parameter * gravi_parameter_add_metrology (cpl_parameterlist *self);

cpl_parameter * gravi_parameter_add_static_name (cpl_parameterlist *self);
cpl_parameter * gravi_parameter_add_debug_file (cpl_parameterlist *self);
cpl_parameter * gravi_parameter_add_biassub_file (cpl_parameterlist *self);
cpl_parameter * gravi_parameter_add_spectrum_file (cpl_parameterlist *self);
cpl_parameter * gravi_parameter_add_preproc_file (cpl_parameterlist *self);
cpl_parameter * gravi_parameter_add_p2vmred_file (cpl_parameterlist *self);
cpl_parameter * gravi_parameter_add_vis_file (cpl_parameterlist *self);
cpl_parameter * gravi_parameter_add_astro_file (cpl_parameterlist *self);

cpl_parameter * gravi_parameter_add_biasmethod (cpl_parameterlist *self);
cpl_parameter * gravi_parameter_add_extract (cpl_parameterlist *self);

cpl_parameter * gravi_parameter_add_average_vis (cpl_parameterlist *self);
cpl_parameter * gravi_parameter_copy_fluxdata (cpl_parameterlist *self);
cpl_parameter * gravi_parameter_add_force_uncertainties (cpl_parameterlist *self);

cpl_error_code gravi_parameter_add_compute_snr (cpl_parameterlist *self, int isCalib);
cpl_error_code gravi_parameter_add_compute_signal (cpl_parameterlist *self, int isCalib);
cpl_error_code gravi_parameter_add_rejection (cpl_parameterlist *self, int iscalib);
cpl_error_code gravi_parameter_add_compute_vis (cpl_parameterlist *self, int iscalib);
cpl_error_code gravi_parameter_add_image (cpl_parameterlist *self);

cpl_frameset * gravi_frameset_extract (cpl_frameset * frameset, const char ** frame_tags, int nb_tabs);
cpl_frameset * gravi_frameset_extract_wave_map(cpl_frameset * );
cpl_frameset * gravi_frameset_extract_bad_map(cpl_frameset * );
cpl_frameset * gravi_frameset_extract_biasmask_map(cpl_frameset * );
cpl_frameset * gravi_frameset_extract_dark_data (cpl_frameset * );
cpl_frameset * gravi_frameset_extract_dark_map(cpl_frameset * );
cpl_frameset * gravi_frameset_extract_p2vm_data(cpl_frameset * );
cpl_frameset * gravi_frameset_extract_disp_data(cpl_frameset * );
cpl_frameset * gravi_frameset_extract_wave_data(cpl_frameset * );
cpl_frameset * gravi_frameset_extract_wavesc_data(cpl_frameset * );
cpl_frameset * gravi_frameset_extract_flat_data(cpl_frameset * );
cpl_frameset * gravi_frameset_extract_flat_map(cpl_frameset * );
cpl_frameset * gravi_frameset_extract_p2vm_map(cpl_frameset * );
cpl_frameset * gravi_frameset_extract_tf_calib (cpl_frameset * );
cpl_frameset * gravi_frameset_extract_vis_calib (cpl_frameset * );
cpl_frameset * gravi_frameset_extract_vis_science (cpl_frameset * );
cpl_frameset * gravi_frameset_extract_sky_data(cpl_frameset * );
cpl_frameset * gravi_frameset_extract_wavelamp_data(cpl_frameset * );
cpl_frameset * gravi_frameset_extract_dispvis_data(cpl_frameset * );
cpl_frameset * gravi_frameset_extract_wavelamp_map(cpl_frameset * frameset);
cpl_frameset * gravi_frameset_extract_disp_map (cpl_frameset * frameset);
cpl_frameset * gravi_frameset_extract_met_pos (cpl_frameset * frameset);
cpl_frameset * gravi_frameset_extract_fringe_data (cpl_frameset * frameset);
cpl_frameset * gravi_frameset_extract_p2vmred_data (cpl_frameset * frameset);
cpl_frameset * gravi_frameset_extract_piezotf_data (cpl_frameset * frameset);
cpl_frameset * gravi_frameset_extract_diamcat_map (cpl_frameset * frameset);
cpl_frameset * gravi_frameset_extract_eop_map (cpl_frameset * frameset);
cpl_frameset * gravi_frameset_extract_patch (cpl_frameset * frameset);
cpl_frameset * gravi_frameset_extract_static_param (cpl_frameset * frameset);
cpl_frameset * gravi_frameset_extract_wave_param (cpl_frameset * frameset);

const char * gravi_param_get_string (const cpl_parameterlist * parlist, const char * name);
double gravi_param_get_double (const cpl_parameterlist *, const char *);
int gravi_param_get_bool (const cpl_parameterlist *, const char *);
int gravi_param_get_int (const cpl_parameterlist *, const char *);

double gravi_param_get_double_default (const cpl_parameterlist *, const char *, double);
int gravi_param_get_bool_default (const cpl_parameterlist *, const char *, int);
int gravi_param_get_int_default (const cpl_parameterlist *, const char *, int);
const char * gravi_param_get_string_default (const cpl_parameterlist * parlist, const char * name, const char * def);

cpl_error_code gravi_check_frameset (cpl_frameset *frameset, const char * tag, int min, int max);

#endif
