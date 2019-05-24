/* $Id: gravi_pfits.h,v 1.8 2007/07/31 06:10:40 llundin Exp $
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

#ifndef GRAVI_PFITS_H
#define GRAVI_PFITS_H

/*-----------------------------------------------------------------------------
                                   Includes
 -----------------------------------------------------------------------------*/

#include <cpl.h>

/*-----------------------------------------------------------------------------
                                   Defines
 -----------------------------------------------------------------------------*/

#define DET_DIT "IPAG DET DIT"
#define GRAVI_DET_DIT "ESO DET2 SEQ1 DIT"
#define GRAVI_NIGHT_OBS "ESO PRO NIGHT OBS"
#define DPR_TYPE	"ESO DPR TYPE"

#define GRAVI_PRIMARY_HDR_EXT "PRIMARY_HDR"

#define GRAVI_IMAGING_DATA_ACQ_EXT "IMAGING_DATA_ACQ"
#define GRAVI_IMAGING_DATA_FT_EXT "IMAGING_DATA_FT"
#define GRAVI_IMAGING_DATA_SC_EXT "IMAGING_DATA_SC"
#define GRAVI_IMAGING_ERR_SC_EXT "IMAGING_ERR_SC"
#define GRAVI_IMAGING_ERR_FT_EXT "IMAGING_ERR_FT"
#define GRAVI_IMAGING_MASK_SC_EXT "IMAGING_MASK_SC"

#define GRAVI_IMAGING_DETECTOR_SC_EXT "IMAGING_DETECTOR_SC"
#define GRAVI_IMAGING_DETECTOR_FT_EXT "IMAGING_DETECTOR_FT"
#define GRAVI_IMAGING_DETECTOR_EXT(type) (type==GRAVI_SC?GRAVI_IMAGING_DETECTOR_SC_EXT:GRAVI_IMAGING_DETECTOR_FT_EXT)

#define GRAVI_SPECTRUM_DATA_SC_EXT "SPECTRUM_DATA_SC"
#define GRAVI_SPECTRUM_DATA_FT_EXT "SPECTRUM_DATA_FT"
#define GRAVI_SPECTRUM_DATA_EXT(type) (type==GRAVI_SC?GRAVI_SPECTRUM_DATA_SC_EXT:GRAVI_SPECTRUM_DATA_FT_EXT)
#define GRAVI_SPECTRUMFLAT_DATA_SC_EXT "SPECTRUMFLAT_DATA_SC"

#define GRAVI_BIAS_MASK_SC_EXT "BIAS_MASK_SC"

#define GRAVI_METROLOGY_EXT "METROLOGY"
#define GRAVI_METROLOGY_ERR_EXT "METROLOGY_ERR"
#define GRAVI_OPDC_EXT "OPDC"
#define GRAVI_WAVE_ARGON_EXT "WAVE_ARGON"
#define GRAVI_WAVE_DATA_FT_EXT "WAVE_DATA_FT"
#define GRAVI_WAVE_DATA_SC_EXT "WAVE_DATA_SC"
#define GRAVI_WAVE_DATA_FT_EXT "WAVE_DATA_FT"
#define GRAVI_WAVE_DATA_EXT(type) (type==GRAVI_SC?GRAVI_WAVE_DATA_SC_EXT:GRAVI_WAVE_DATA_FT_EXT)
#define GRAVI_WAVE_FIBRE_FT_EXT "WAVE_FIBRE_FT"
#define GRAVI_WAVE_FIBRE_SC_EXT "WAVE_FIBRE_SC"
#define GRAVI_WAVE_FIBRE_EXT(type) (type==GRAVI_SC?GRAVI_WAVE_FIBRE_SC_EXT:GRAVI_WAVE_FIBRE_FT_EXT)
#define GRAVI_P2VM_MET_EXT "P2VM_MET"
#define GRAVI_P2VM_DATA_SC_EXT "P2VM_SC"
#define GRAVI_P2VM_DATA_FT_EXT "P2VM_FT"
#define GRAVI_P2VM_DATA_EXT(type) (type==GRAVI_SC?GRAVI_P2VM_DATA_SC_EXT:GRAVI_P2VM_DATA_FT_EXT)
#define GRAVI_FDDL_EXT "FDDL"

#define GRAVI_PROFILE_DATA_EXT "PROFILE_DATA"
#define GRAVI_PROFILE_PARAMS_EXT "PROFILE_PARAMS"
#define GRAVI_IMAGING_DETECTOR_SC_EXT "IMAGING_DETECTOR_SC"
#define GRAVI_IMAGING_DETECTOR_FT_EXT "IMAGING_DETECTOR_FT"
#define GRAVI_OI_ARRAY_EXT "OI_ARRAY"
#define GRAVI_ARRAY_GEOMETRY_EXT "ARRAY_GEOMETRY"
#define GRAVI_OPTICAL_TRAIN_EXT "OPTICAL_TRAIN"
#define GRAVI_OI_TARGET_EXT "OI_TARGET"

#define GRAVI_OI_VIS_MET_EXT "OI_VIS_MET"
#define GRAVI_OI_VIS_ACQ_EXT "OI_VIS_ACQ"

#define GRAVI_OI_WAVELENGTH_EXT "OI_WAVELENGTH"
#define GRAVI_OI_VIS_EXT "OI_VIS"
#define GRAVI_OI_FLUX_EXT "OI_FLUX"
#define GRAVI_OI_T3_EXT "OI_T3"
#define GRAVI_OI_VIS2_EXT "OI_VIS2"

#define GRAVI_NAN_DOUBLE -99.
#define GRAVI_NAN_FLOAT -99.
#define GRAVI_NAN_INT -99

//#define LAMBDA_MET 0.000001908287 /* updated with actual valeur - pkervell - 19aug15 */
#define LAMBDA_MET 0.000001908254

/* QC WAVE*/
#define QC_PHASECHI2 "ESO QC PHASE_CALIBRATION_CHI2"
#define QC_CHI2WAVE(type) (type==GRAVI_SC?"ESO QC METFITRMS WAVESC":"ESO QC METFITRMS WAVEFT")
#define QC_MINWAVE(type) (type==GRAVI_SC?"ESO QC MINWAVE SC":"ESO QC MINWAVE FT")
#define QC_MAXWAVE(type) (type==GRAVI_SC?"ESO QC MAXWAVE SC":"ESO QC MAXWAVE FT")
#define QC_RMS_RESIDUALS(type) (type==GRAVI_SC?"ESO QC RMSWAVE SC":"ESO QC RMSWAVE FT")
#define OPD_COEFF_SIGN(type) (type==GRAVI_SC?"ESO QC OPD_COEFF_SIGN SC":"ESO QC OPD_COEFF_SIGN FT")

/* QC DARK */
#define QC_MEANDARK_SC "ESO QC MEDIANDARK SC"
#define QC_DARKRMS_SC "ESO QC DARKRMS SC"
#define QC_MEANSKY_SC "ESO QC MEDIANSKY SC"
#define QC_SKYRMS_SC "ESO QC SKYRMS SC"

#define QC_MEANDARK_FT "ESO QC MEANDARK FT"
#define QC_DARKRMS_FT "ESO QC DARKRMS FT"
#define QC_MEANSKY_FT "ESO QC MEANSKY FT"
#define QC_SKYRMS_FT "ESO QC SKYRMS FT"

#define QC_MEANDARK_MET "ESO QC MEANDARK MET"
#define QC_DARKRMS_MET "ESO QC DARKRMS MET"

#define QC_MEANDARK "ESO QC MEANDARK"
#define QC_DARKRMS "ESO QC DARKRMS"

/* QC FLAT */
#define QC_MEANGAIN_SC "ESO QC MEANGAIN SC"
#define QC_BADPIX_SC "ESO QC BADPIX SC"
#define QC_BADPIX_DARK_SC "ESO QC BADPIX_DARK SC"
#define QC_BADPIX_RMS_SC "ESO QC BADPIX_RMS SC"
#define QC_BADPIX_FLAT_SC "ESO QC BADPIX_FLAT SC"
#define QC_MEANGAIN_FT "ESO QC MEANGAIN FT"
#define QC_BADPIX_FT "ESO QC BADPIX FT"
#define PROFILE_STARTX "ESO PRO PROFILE STARTX"
#define PROFILE_FULLSTARTX "ESO PRO PROFILE FULLSTARTX"
#define PROFILE_NX "ESO PRO PROFILE NX"
    

/* QC P2VM */
#define QC_MEANCOH_SC 		"ESO QC P2VM_COHERENCE_AVG_SC"
#define QC_RMSCOH_SC 		"ESO QC P2VM_COHERENCERMS_AVG_SC"
#define QC_RMSPHASE_SC 		"ESO QC P2VM_PHASERMS_AVG_SC"
#define QC_MEANCOH_FT 		"ESO QC P2VM_COHERENCE_AVG_FT"
#define QC_RMSCOH_FT		"ESO QC P2VM_COHERENCERMS_AVG_FT"
#define QC_RMSPHASE_FT		"ESO QC P2VM_PHASERMS_AVG_FT"


/* Type of data */
#define GRAVI_SC 0 
#define GRAVI_FT 1
#define GRAVI_TYPE(type) (type==GRAVI_SC?"SC":"FT")

enum gravi_detector_type
{
    GRAVI_DET_FT,
    GRAVI_DET_SC,
    GRAVI_DET_ALL
};

/* Modes */
#define MODE_SINGLE 1
#define MODE_DUAL   2
#define MODE_ONAXIS  1
#define MODE_OFFAXIS 2

/* INSNAME_SC or INSNAME_FT. Deal with polar */
#define INSNAME_FT_P1 "GRAVITY_FT_P1"
#define INSNAME_FT_P2 "GRAVITY_FT_P2"
#define INSNAME_SC_P1 "GRAVITY_SC_P1"
#define INSNAME_SC_P2 "GRAVITY_SC_P2"
#define INSNAME_FT "GRAVITY_FT"
#define INSNAME_SC "GRAVITY_SC"
#define GRAVI_INSNAME(type,pol,npol) (type==GRAVI_SC ? (npol==1?INSNAME_SC:(pol==0?INSNAME_SC_P1:INSNAME_SC_P2) ) : (npol==1?INSNAME_FT:(pol==0?INSNAME_FT_P1:INSNAME_FT_P2) ) )

/* EXTVER = {10,11,12} for SC
   EXTVER = {20,21,22} for FT
   0,1,2 beeing the combined mode, the first polar, and the second polar */
#define GRAVI_EXTVER(type,pol,npol) ( (type+1)*10 + (npol==1?0:(pol+1)) )

#define gravi_pfits_get_sobj_decep(plist) gravi_pfits_get_decep(plist, 1.0)
#define gravi_pfits_get_sobj_raep(plist) gravi_pfits_get_raep(plist, 1.0)
#define gravi_pfits_get_mid_decep(plist) gravi_pfits_get_decep(plist, 0.5)
#define gravi_pfits_get_mid_raep(plist) gravi_pfits_get_raep(plist, 0.5)
#define gravi_pfits_get_type_decep(plist,type) (type==GRAVI_SC?gravi_pfits_get_sobj_decep(plist):gravi_pfits_get_robj_decep(plist))
#define gravi_pfits_get_type_raep(plist,type) (type==GRAVI_SC?gravi_pfits_get_sobj_raep(plist):gravi_pfits_get_robj_raep(plist))

#define gravi_pfits_get_dit(plist, type) (type==GRAVI_SC ? gravi_pfits_get_dit_sc(plist) : gravi_pfits_get_dit_ft(plist))

/*-----------------------------------------------------------------------------
                              Private prototypes
 ----------------------------------------------------------------------------- */


const char * gravi_pfits_get_met_ph(const cpl_propertylist * );
int gravi_pfits_has_gdzero (const cpl_propertylist * plist, int tel);
double gravi_pfits_get_gdzero (const cpl_propertylist * plist, int tel);
int gravi_pfits_has_oplzero (const cpl_propertylist * plist, int tel);
double gravi_pfits_get_oplzero (const cpl_propertylist * plist, int tel);
double gravi_pfits_get_metfc_lockmjd (const cpl_propertylist * plist, int tel);
double gravi_pfits_get_met_wavelength (const cpl_propertylist * plist);
double gravi_pfits_get_met_wavelength_mean (const cpl_propertylist * plist, cpl_table * met_table);

const char * gravi_pfits_get_start_sc (const cpl_propertylist * plist);
const char * gravi_pfits_get_start_acqcam (const cpl_propertylist * plist);
const char * gravi_pfits_get_start_prcacq (const cpl_propertylist * plist);

double gravi_pfits_get_fddlwindow (const cpl_propertylist * plist);

double gravi_pfits_get_period_sc (const cpl_propertylist * plist);
double gravi_pfits_get_period_acqcam (const cpl_propertylist * plist);
double gravi_pfits_get_dit_ft (const cpl_propertylist * plist);
double gravi_pfits_get_dit_sc (const cpl_propertylist * plist);
double gravi_pfits_get_dit_acqcam (const cpl_propertylist * plist);

double gravi_pfits_get_time_sc (const cpl_propertylist * header, cpl_size row);
double gravi_pfits_get_time_acqcam (const cpl_propertylist * header, cpl_size row);

double gravi_pfits_get_mjd (const cpl_propertylist * plist);

double gravi_pfits_get_ft_gain (const cpl_propertylist * plist);
double gravi_pfits_get_sc_gain (const cpl_propertylist * plist);

const char * gravi_pfits_get_extname(const cpl_propertylist * );
int gravi_pfits_get_extension_type(const cpl_propertylist * plist);

double gravi_pfits_get_sobj_diam (const cpl_propertylist * plist);
double gravi_pfits_get_robj_diam (const cpl_propertylist * plist);
double gravi_pfits_get_diameter(const cpl_propertylist * plist, int type_data);

double gravi_pfits_get_ptfc_acqcam (const cpl_propertylist * plist, int spot);
double gravi_pfits_get_drotoff (const cpl_propertylist * plist, int tel);
double gravi_pfits_get_fangle_acqcam (const cpl_propertylist * plist, int tel);

const char * gravi_pfits_get_resolution(const cpl_propertylist * plist);
const char * gravi_pfits_get_dpr_type(const cpl_propertylist * plist);
const char * gravi_pfits_get_insname(const cpl_propertylist * plist);
const char * gravi_pfits_get_spec_res(const cpl_propertylist * plist);
const char * gravi_pfits_get_pola_mode(const cpl_propertylist * plist, int type_data);

int gravi_pfits_is_calib (const cpl_propertylist * plist);
int gravi_pfits_get_pola_num(const cpl_propertylist * plist, int type_data );
int gravi_pfits_get_mode (const cpl_propertylist * plist);
int gravi_pfits_get_axis (const cpl_propertylist * plist);
const char * gravi_pfits_get_mode_name (const cpl_propertylist * plist);
int gravi_data_frame_get_mode(const cpl_frame * frame);

int gravi_pfits_get_window_start (const cpl_propertylist * plist);
int gravi_pfits_get_startx (const cpl_propertylist * plist);
int gravi_pfits_get_fullstartx (const cpl_propertylist * plist);
int gravi_pfits_get_nx (const cpl_propertylist * plist);

const char * gravi_pfits_get_robj(const cpl_propertylist * plist);
const char * gravi_pfits_get_sobj(const cpl_propertylist * plist);

double gravi_pfits_get_decep(const cpl_propertylist * plist, double coef);
double gravi_pfits_get_raep(const cpl_propertylist * plist, double coef);
double gravi_pfits_get_robj_decep(const cpl_propertylist * plist);
double gravi_pfits_get_robj_raep(const cpl_propertylist * plist);
double gravi_pfits_get_sobj_x(const cpl_propertylist * plist);
double gravi_pfits_get_sobj_y(const cpl_propertylist * plist);

double gravi_pfits_get_plx(const cpl_propertylist * plist);
double gravi_pfits_get_pmra(const cpl_propertylist * plist);
double gravi_pfits_get_pmdec(const cpl_propertylist * plist);

double gravi_pfits_get_geoelev(const cpl_propertylist * plist);
double gravi_pfits_get_geolat(const cpl_propertylist * plist);
double gravi_pfits_get_geolon(const cpl_propertylist * plist);

double gravi_ra_to_rad (const char *stri);
double gravi_dec_to_rad (const char *stri);
double gravi_convert_to_mjd (const char * );
char * gravi_convert_to_timestamp (double mjd);
cpl_error_code gravi_pfits_ensure_double(cpl_propertylist * self, const char * name);

cpl_propertylist *  gravi_plist_get_qc (cpl_propertylist * );
cpl_propertylist *  gravi_plist_get_oifits_keywords (cpl_propertylist * header);

double gravi_pfits_get_double(const cpl_propertylist * self, const char * name);
double gravi_pfits_get_double_default(const cpl_propertylist * plist, const char *name,double def);
double gravi_pfits_get_double_silentdefault(const cpl_propertylist * plist, const char *name,double def);
const char * gravi_pfits_get_string_default (const cpl_propertylist * plist, const char *name, const char *def);

cpl_error_code gravi_pfits_add_check (cpl_propertylist * header, char *msg);
cpl_error_code gravi_pfits_add_pipe_build (cpl_propertylist * header);

cpl_error_code gravi_pfits_update_double (cpl_propertylist * plist, const char * name, double value);
cpl_error_code gravi_pfits_update_int (cpl_propertylist * plist, const char * name, int value);


#endif
