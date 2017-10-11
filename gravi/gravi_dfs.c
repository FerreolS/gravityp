/* $Id: gravi_dfs.c,v 1.6 2011/04/31 06:10:40 llundin Exp $
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
 * @defgroup gravi_dfs  DFS related functions
 */
/**@{*/

/*-----------------------------------------------------------------------------
                                   Includes
 -----------------------------------------------------------------------------*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "time.h"
#include <string.h>
#include <math.h>
#include <cpl.h>

#include "gravi_dfs.h"
#include "gravi_utils.h"

/*-----------------------------------------------------------------------------
                             Functions code
 -----------------------------------------------------------------------------*/

void gravity_print_banner (void) 
{ 
    cpl_msg_info(__func__, "**********************************************"); 
    cpl_msg_info(__func__, "  Welcome to GRAVITY Pipeline release %s", 
                 PACKAGE_VERSION);
    cpl_msg_info(__func__, "  Last rebuilt at %s %s",__DATE__,__TIME__);
    cpl_msg_info(__func__, "**********************************************"); 
} 

/*----------------------------------------------------------------------------*/
/**
 * @brief    Set the group as RAW or CALIB in a frameset
 * @param    set     the input frameset
 * @return   CPL_ERROR_NONE iff OK
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_dfs_set_groups(cpl_frameset * set)
{
    cpl_ensure_code (set, CPL_ERROR_NULL_INPUT);
  
    cpl_errorstate prestate = cpl_errorstate_get();
    cpl_frame * frame = NULL;
    int         i, nb_frame;
    nb_frame = cpl_frameset_get_size(set);


    /* Loop on frames */
    for (i = 0; i < nb_frame; i++) {

    	frame = cpl_frameset_get_position(set, i);
        const char * tag = cpl_frame_get_tag(frame);

        if (tag == NULL) {
            cpl_msg_warning(cpl_func, "Frame %d has no tag", i);
        } else if ((!strcmp(tag, GRAVI_DARK_RAW)) ||
				   (!strcmp(tag, GRAVI_FLAT_RAW)) ||
				   (!strcmp(tag, GRAVI_WAVE_RAW)) ||
				   (!strcmp(tag, GRAVI_WAVELAMP_RAW)) ||
				   (!strcmp(tag, GRAVI_P2VM_RAW)) ||
				   (!strcmp(tag, GRAVI_SINGLE_SCIENCE_RAW)) ||
				   (!strcmp(tag, GRAVI_SINGLE_CALIB_RAW)) ||
				   (!strcmp(tag, GRAVI_DUAL_SCIENCE_RAW)) ||
				   (!strcmp(tag, GRAVI_DUAL_CALIB_RAW))||
				   (!strcmp(tag, GRAVI_DUAL_SKY_RAW)) ||
				   (!strcmp(tag, GRAVI_SINGLE_SKY_RAW))||
                   (!strcmp(tag, GRAVI_VIS_SINGLE_CALIB)) ||
				   (!strcmp(tag, GRAVI_VIS_SINGLE_SCIENCE)) ||
                   (!strcmp(tag, GRAVI_VIS_DUAL_CALIB)) ||
                   (!strcmp(tag, GRAVI_VIS_DUAL_SCIENCE)) ||
				   (!strcmp(tag, GRAVI_P2VMRED_SINGLE_CALIB)) ||
				   (!strcmp(tag, GRAVI_P2VMRED_SINGLE_SCIENCE)) ||
				   (!strcmp(tag, GRAVI_P2VMRED_DUAL_CALIB)) ||
				   (!strcmp(tag, GRAVI_P2VMRED_DUAL_SCIENCE)) ||
				   (!strcmp(tag, GRAVI_MIRA_INPUT_PROCATG))||
				   (!strcmp(tag, GRAVI_VIS_SINGLE_CALIBRATED)) ||
				   (!strcmp(tag, GRAVI_VIS_DUAL_CALIBRATED)) ||
				   (!strcmp(tag, GRAVI_DISP_RAW))  ){
		  /* RAW frames */
		  cpl_frame_set_group(frame, CPL_FRAME_GROUP_RAW);
        }else if ((!strcmp(tag, GRAVI_DARK_MAP)) ||
				  (!strcmp(tag, GRAVI_FLAT_MAP)) ||
				  (!strcmp(tag, GRAVI_WAVE_MAP)) ||
				  (!strcmp(tag, GRAVI_P2VM_MAP)) ||
				  (!strcmp(tag, GRAVI_BAD_MAP)) ||
				  (!strcmp(tag, GRAVI_BIASMASK_MAP)) ||
				  (!strcmp(tag, GRAVI_PREPROC)) ||
				  (!strcmp(tag, GRAVI_TF_SINGLE_SCIENCE)) ||
				  (!strcmp(tag, GRAVI_TF_SINGLE_CALIB)) ||
				  (!strcmp(tag, GRAVI_WAVELAMP_MAP)) ||
				  (!strcmp(tag, GRAVI_TF_DUAL_SCIENCE))  ||
				  (!strcmp(tag, GRAVI_TF_DUAL_CALIB)) || 
				  (!strcmp(tag, GRAVI_ZP_CAL)) ||
				  (!strcmp(tag, GRAVI_DISP_VIS)) ||
				  (!strcmp(tag, GRAVI_DIAMETER_CAT)) ||
				  (!strcmp(tag, GRAVI_DISP_MODEL)) ||
				  (!strcmp(tag, GRAVI_DIODE_POSITION))||
	              (!strcmp(tag, GRAVI_KEY_PATCH))){
        	/* CALIB frames */
        	cpl_frame_set_group(frame, CPL_FRAME_GROUP_CALIB);
        }else if (
        		(!strcmp(tag, GRAVI_MIRA_OUTPUT_PROCATG)) ||
                (!strcmp(tag, GRAVI_NAB_CAL)) ){
        	/* PRODUCT frames */
        	cpl_frame_set_group(frame, CPL_FRAME_GROUP_PRODUCT);
        }

    }

    if (!cpl_errorstate_is_equal(prestate)) {
        return cpl_error_set_message(cpl_func, cpl_error_get_code(),
                                     "Could not identify RAW and CALIB "
                                     "frames");
    }

    return CPL_ERROR_NONE;
}

/*---------------------------------------------------------------------------*/

cpl_error_code gravi_parameter_disable (cpl_parameter * p)
{
    cpl_ensure_code (p, CPL_ERROR_NULL_INPUT);
    cpl_parameter_disable (p, CPL_PARAMETER_MODE_CLI);
    cpl_parameter_disable (p, CPL_PARAMETER_MODE_ENV);
    cpl_parameter_disable (p, CPL_PARAMETER_MODE_CFG);
    return CPL_ERROR_NONE;
}

/*---------------------------------------------------------------------------*/

cpl_parameter * gravi_parameter_add_badpix (cpl_parameterlist *self)
{
    cpl_ensure (self, CPL_ERROR_NULL_INPUT, NULL);
    
    cpl_parameter *p;
    p = cpl_parameter_new_value ("gravity.calib.bad-dark-threshold", CPL_TYPE_INT,
                                 "the rms factor for "
                                 "dark bad pixel threshold",
                                 "gravity.calib", 10);
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "bad-dark-threshold");
    cpl_parameter_disable (p, CPL_PARAMETER_MODE_ENV);
	cpl_parameterlist_append (self, p);

    return p;
}

cpl_parameter * gravi_parameter_add_profile (cpl_parameterlist *self)
{
    cpl_ensure (self, CPL_ERROR_NULL_INPUT, NULL);
    
    cpl_parameter *p;
    
    /* Method for profile */
    p = cpl_parameter_new_enum ("gravity.calib.profile-mode", CPL_TYPE_STRING,
                                "Method to compute the extraction profile. "
                                "PROFILE corresponds to the pixel intensities measured in the "
                                "FLAT files (Gaussian like with FWHM of approx 1.5 pixel). "
                                "This is the AUTO option for the Low and Med spectral resolution. "
                                "GAUSS corresponds to a Gaussian fit of the (non-zero) pixel intensities measured "
                                "in the FLAT files. BOX corresponds to a box-card of 6 pixels centered "
                                "on the spectra measured in the FLAT files. This is the AUTO option for High "
                                "spectral resolution",
                                "gravity.calib", "AUTO",
                                4, "AUTO", "PROFILE", "GAUSS", "BOX");
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "profile-mode");
    cpl_parameter_disable (p, CPL_PARAMETER_MODE_ENV);
	cpl_parameterlist_append (self, p);

    /* How to deal with bad-pixels  */
	p = cpl_parameter_new_value ("gravity.calib.force-badpix-to-zero", CPL_TYPE_BOOL,
                                 "Force the badpixel to zero in profile",
                                 "gravity.calib", TRUE);
	cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "force-badpix-to-zero");
	cpl_parameter_disable (p, CPL_PARAMETER_MODE_ENV);
	cpl_parameterlist_append (self, p);

	/* The width of the profile element */
    p = cpl_parameter_new_value ("gravity.calib.profile-width", CPL_TYPE_INT,
                                 "Width of the detector window extracted around the default "
                                 "position of each spectrum, and on which the profile "
                                 "will be applied to perform the extraction.",
                                 "gravity.calib", 6);
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "profile-width");
    cpl_parameter_disable (p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append (self, p);
    

    return p;
}

cpl_parameter * gravi_parameter_add_wave (cpl_parameterlist *self)
{
    cpl_ensure (self, CPL_ERROR_NULL_INPUT, NULL);

    cpl_parameter *p;

    /* Method for profile */
    p = cpl_parameter_new_value ("gravity.calib.force-wave-ft-equal", CPL_TYPE_BOOL,
                                "Force the spatial order of the wavelength 2D fit for FT to "
                                 "zero (so all region share the same calibration). "
                                 "This is used to build the P2VM calibration of the TAC "
                                 "real-time code running on the instrument ifself.",
                                 "gravity.calib", FALSE);
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "force-wave-ft-equal");
    cpl_parameter_disable (p, CPL_PARAMETER_MODE_ENV);
	cpl_parameterlist_append (self, p);

    return p;
}

cpl_parameter * gravi_parameter_add_static_name (cpl_parameterlist *self)
{
    cpl_ensure (self, CPL_ERROR_NULL_INPUT, NULL);
    
    cpl_parameter *p;
	p = cpl_parameter_new_value ("gravity.dfs.static-name", CPL_TYPE_BOOL,
                                 "Use static names for the products (for ESO)",
                                 "gravity.dfs", FALSE);
	cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "static-name");
	cpl_parameter_disable (p, CPL_PARAMETER_MODE_ENV);
	cpl_parameterlist_append (self, p);

    return p;
}

cpl_parameter * gravi_parameter_add_debug_file (cpl_parameterlist *self)
{
    cpl_ensure (self, CPL_ERROR_NULL_INPUT, NULL);
    
    cpl_parameter *p;
	p = cpl_parameter_new_value ("gravity.dfs.debug-file", CPL_TYPE_BOOL,
                                 "Save additional debug file(s)",
                                 "gravity.dfs", FALSE);
	cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "debug-file");
	cpl_parameter_disable (p, CPL_PARAMETER_MODE_ENV);
	cpl_parameterlist_append (self, p);

    return p;
}

cpl_parameter * gravi_parameter_add_biassub_file (cpl_parameterlist *self)
{
    cpl_ensure (self, CPL_ERROR_NULL_INPUT, NULL);
    
    cpl_parameter *p;
	p = cpl_parameter_new_value ("gravity.dfs.bias-subtracted-file", CPL_TYPE_BOOL,
                                 "Save the BIAS_SUBTRACTED intermediate product",
                                 "gravity.dfs", FALSE);
	cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "bias-subtracted-file");
	cpl_parameter_disable (p, CPL_PARAMETER_MODE_ENV);
	cpl_parameterlist_append (self, p);
    
    return p;
}

cpl_parameter * gravi_parameter_add_spectrum_file (cpl_parameterlist *self)
{
    cpl_ensure (self, CPL_ERROR_NULL_INPUT, NULL);
    
    cpl_parameter *p;
	p = cpl_parameter_new_value ("gravity.dfs.spectrum-file", CPL_TYPE_BOOL,
                                 "Save the SPECTRUM intermediate product",
                                 "gravity.dfs", FALSE);
	cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "spectrum-file");
	cpl_parameter_disable (p, CPL_PARAMETER_MODE_ENV);
	cpl_parameterlist_append (self, p);
    
    return p;
}

cpl_parameter * gravi_parameter_add_preproc_file (cpl_parameterlist *self)
{
    cpl_ensure (self, CPL_ERROR_NULL_INPUT, NULL);
    
    cpl_parameter *p;
	p = cpl_parameter_new_value ("gravity.dfs.preproc-file", CPL_TYPE_BOOL,
                                 "Save the PREPROC intermediate product",
                                 "gravity.dfs", FALSE);
	cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "preproc-file");
	cpl_parameter_disable (p, CPL_PARAMETER_MODE_ENV);
	cpl_parameterlist_append (self, p);
    
    return p;
}

cpl_parameter * gravi_parameter_add_p2vmred_file (cpl_parameterlist *self)
{
    cpl_ensure (self, CPL_ERROR_NULL_INPUT, NULL);
    
    cpl_parameter *p;
	p = cpl_parameter_new_value ("gravity.dfs.p2vmred-file", CPL_TYPE_BOOL,
                                 "Save the P2VMRED intermediate product",
                                 "gravity.dfs", FALSE);
	cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "p2vmreduced-file");
	cpl_parameter_disable (p, CPL_PARAMETER_MODE_ENV);
	cpl_parameterlist_append (self, p);
    
    return p;
}

cpl_parameter * gravi_parameter_add_vis_file (cpl_parameterlist *self)
{
    cpl_ensure (self, CPL_ERROR_NULL_INPUT, NULL);
    
    cpl_parameter *p;
	p = cpl_parameter_new_value ("gravity.dfs.vis-file", CPL_TYPE_BOOL,
                                 "Save the VIS intermediate product",
                                 "gravity.dfs", FALSE);
	cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "vis-file");
	cpl_parameter_disable (p, CPL_PARAMETER_MODE_ENV);
	cpl_parameterlist_append (self, p);
    
    return p;
}

cpl_parameter * gravi_parameter_add_astro_file (cpl_parameterlist *self)
{
    cpl_ensure (self, CPL_ERROR_NULL_INPUT, NULL);
    
    cpl_parameter *p;
	p = cpl_parameter_new_value ("gravity.dfs.astro-file", CPL_TYPE_BOOL,
                                 "Save the ASTROREDUCED intermediate product",
                                 "gravity.dfs", FALSE);
	cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "astro-file");
	cpl_parameter_disable (p, CPL_PARAMETER_MODE_ENV);
	cpl_parameterlist_append (self, p);
    
    return p;
}

cpl_parameter * gravi_parameter_add_biasmethod (cpl_parameterlist *self)
{
    cpl_ensure (self, CPL_ERROR_NULL_INPUT, NULL);
    
    cpl_parameter *p;
    p = cpl_parameter_new_enum ("gravity.preproc.bias-method", CPL_TYPE_STRING,
                                "Method to average the biaspixels when cleaning-up\n "
                                "the SC detector (only applied to MED and LOW). Ideally\n "
                                "the same value shall be used when reducing the DARK\n "
                                "with gravity_dark and the OBJECT with gravity_vis.",
                                "gravity.preproc", "MEDIAN",
                                2, "MEDIAN", "MEDIAN_PER_COLUMN");
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "bias-method");
    cpl_parameter_disable (p, CPL_PARAMETER_MODE_ENV);
	cpl_parameterlist_append (self, p);

    return p;
}

cpl_parameter * gravi_parameter_add_metrology (cpl_parameterlist *self)
{
    cpl_ensure (self, CPL_ERROR_NULL_INPUT, NULL);
    
    cpl_parameter *p;
	p = cpl_parameter_new_value ("gravity.metrology.acq-correction-delay",
                                 CPL_TYPE_DOUBLE,
                                 "Delay between the end of ACQ frame and correction\n "
                                 "offset seen by the metrology diodes, in seconds.",
                                 "gravity.metrology", 0.1);
	cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "acq-correction-delay");
	cpl_parameter_disable (p, CPL_PARAMETER_MODE_ENV);
	cpl_parameterlist_append (self, p);
    
    return p;
}
    
cpl_parameter * gravi_parameter_add_extract (cpl_parameterlist *self)
{
    cpl_ensure (self, CPL_ERROR_NULL_INPUT, NULL);
    
    cpl_parameter *p;
    p = cpl_parameter_new_value ("gravity.preproc.ditshift-sc", CPL_TYPE_INT,
                                 "Shift the time of SC DITs by an integer value to\n "
                                 "account for lost frames in exposure (issue on the\n "
                                 "instrument side, report to instrument team). The\n "
                                 "time of all DITs in exposure are increased by\n "
                                 "ditshift x PERIOD. ditshift can be 0,\n " 
                                 "positive (system has lost one SC DIT), or negative "
                                 "(SC desynchronized).",
                                 "gravity.preproc",0);
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "ditshift-sc");
    cpl_parameter_disable (p, CPL_PARAMETER_MODE_ENV);
	cpl_parameterlist_append (self, p);

    return p;
}

cpl_parameter * gravi_parameter_add_average_vis (cpl_parameterlist *self)
{
    cpl_ensure (self, CPL_ERROR_NULL_INPUT, NULL);
    cpl_parameter *p;
	p = cpl_parameter_new_value ("gravity.postprocess.average-vis", CPL_TYPE_BOOL,
                                 "Average the results from the different input files (if any)\n "
                                 "in the output product, instead of simply appending them.",
                                 "gravity.postprocess", FALSE);
	cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "average-vis");
	cpl_parameter_disable (p, CPL_PARAMETER_MODE_ENV);
	cpl_parameterlist_append (self, p);
    
    return p;
}

cpl_parameter * gravi_parameter_add_force_uncertainties (cpl_parameterlist *self)
{
    cpl_ensure (self, CPL_ERROR_NULL_INPUT, NULL);
    
    cpl_parameter *p;
	p = cpl_parameter_new_value ("gravity.postprocess.fluxerr-sc", CPL_TYPE_DOUBLE,
                                 "Force the uncertainty in FLUX of SC",
                                 "gravity.postprocess", 0.0);
	cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "fluxerr-sc");
	cpl_parameter_disable (p, CPL_PARAMETER_MODE_ENV);
	cpl_parameterlist_append (self, p);

	p = cpl_parameter_new_value ("gravity.postprocess.visamperr-sc", CPL_TYPE_DOUBLE,
                                 "Force the uncertainty in VISAMP of SC",
                                 "gravity.postprocess", 0.0);
	cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "visamperr-sc");
	cpl_parameter_disable (p, CPL_PARAMETER_MODE_ENV);
	cpl_parameterlist_append (self, p);

	p = cpl_parameter_new_value ("gravity.postprocess.visphierr-sc", CPL_TYPE_DOUBLE,
                                 "Force the uncertainty in VISPHI of SC",
                                 "gravity.postprocess", 0.0);
	cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "visphierr-sc");
	cpl_parameter_disable (p, CPL_PARAMETER_MODE_ENV);
	cpl_parameterlist_append (self, p);
    
	p = cpl_parameter_new_value ("gravity.postprocess.vis2err-sc", CPL_TYPE_DOUBLE,
                                 "Force the uncertainty in VIS2 of SC",
                                 "gravity.postprocess", 0.0);
	cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "vis2err-sc");
	cpl_parameter_disable (p, CPL_PARAMETER_MODE_ENV);
	cpl_parameterlist_append (self, p);
    
    return p;
}

cpl_error_code gravi_parameter_add_compute_snr (cpl_parameterlist *self, int isCalib)
{
    cpl_ensure_code (self, CPL_ERROR_NULL_INPUT);
    cpl_parameter *p;
    
    /* Number of FT samples to average to compute SNR and GDELAY */
	p = cpl_parameter_new_value ("gravity.signal.nsmooth-snr-ft", CPL_TYPE_INT,
                                 "Number of samples to average coherently when computing\n "
                                 "the real-time SNR and GDELAY of the FT (shall correspond\n "
                                 "to the atmospheric coherence time). The integration\n "
                                 "window runs from -nsmooth -> +nsmooth.",
                                 "gravity.signal", 5);
	cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "nsmooth-snr-ft");
	cpl_parameter_disable (p, CPL_PARAMETER_MODE_ENV);
	cpl_parameterlist_append (self, p);

    return CPL_ERROR_NONE;
}

cpl_error_code gravi_parameter_add_compute_signal (cpl_parameterlist *self, int isCalib)
{
    cpl_ensure_code (self, CPL_ERROR_NULL_INPUT);
    return CPL_ERROR_NONE;
}

cpl_error_code gravi_parameter_add_rejection (cpl_parameterlist *self, int isCalib)
{
    cpl_ensure_code (self, CPL_ERROR_NULL_INPUT);
    
    cpl_parameter *p;
	p = cpl_parameter_new_value ("gravity.signal.snr-min-ft", CPL_TYPE_DOUBLE,
                                 "SNR threshold to accept FT frames (>0). It raises the first bit (<<0)\n "
                                 "of column REJECTION_FLAG of FT.",
                                 "gravity.signal", isCalib ? 30.0 : 3.0);
	cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "snr-min-ft");
	cpl_parameter_disable (p, CPL_PARAMETER_MODE_ENV);
	cpl_parameterlist_append (self, p);

	/* OPDC_STATE threshold for fringe DET in FT */
	p = cpl_parameter_new_value ("gravity.signal.global-state-min-ft", CPL_TYPE_DOUBLE,
                                 "Minimum OPDC state to accept FT frames (>=0) It raises the second bit\n "
                                 "(<<1) of column REJECTION_FLAG of FT.",
                                 "gravity.signal", 2.0);
	cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "global-state-min-ft");
	cpl_parameter_disable (p, CPL_PARAMETER_MODE_ENV);
	cpl_parameterlist_append (self, p);
    
	p = cpl_parameter_new_value ("gravity.signal.global-state-max-ft", CPL_TYPE_DOUBLE,
                                 "Maximum OPDC state to accept FT frames (>=0) It raises the second bit\n "
                                 "(<<1) of column REJECTION_FLAG of FT.",
                                 "gravity.signal", 4.0);
	cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "global-state-max-ft");
	cpl_parameter_disable (p, CPL_PARAMETER_MODE_ENV);
	cpl_parameterlist_append (self, p);

	/* STATE threshold for fringe DET in FT */
	p = cpl_parameter_new_value ("gravity.signal.state-min-ft", CPL_TYPE_DOUBLE,
                                 "Minimum OPDC state per baseline to accept FT frames (>=0) It raises\n "
                                 "the second bit (<<1) of column REJECTION_FLAG of FT.",
                                 "gravity.signal", 1.0);
	cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "state-min-ft");
	cpl_parameter_disable (p, CPL_PARAMETER_MODE_ENV);
	cpl_parameterlist_append (self, p);
	
	/* Minimum detection ratio to accept SC frame */
	p = cpl_parameter_new_value ("gravity.signal.tracking-min-sc", CPL_TYPE_DOUBLE,
                                 "Minimum ratio of accepted FT frames in order to accept a SC frames (0..1),\n "
                                 "that is, for each SC DIT, the fraction of the time the\n "
                                 "REJECTION_FLAG of the FT is not 0.\n "
                                 "It raises the first bit (<<0) of column REJECTION_FLAG of SC",
                                 "gravity.signal", 0.8);
	cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "tracking-min-sc");
	cpl_parameter_disable (p, CPL_PARAMETER_MODE_ENV);
	cpl_parameterlist_append (self, p);

	/* vFactor threshold to accept SC frame */
	p = cpl_parameter_new_value ("gravity.signal.vfactor-min-sc", CPL_TYPE_DOUBLE,
                                 "vFactor threshold to accept SC frame (0..1).\n ",
                                 "It raises the second bit (<<1) of column REJECTION_FLAG of SC",
                                 "gravity.signal", isCalib ? 0.8 : 0.1);
	cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "vfactor-min-sc");
	cpl_parameter_disable (p, CPL_PARAMETER_MODE_ENV);
	cpl_parameterlist_append (self, p);
    
    /* Flag to activate metrology zero calculation */
    p = cpl_parameter_new_value ("gravity.signal.use-met-zero", CPL_TYPE_BOOL,
                                 "Flag to activate metrology zero calculation in OPD_DISP.\n ",
                                 "If disabled all metrology zeroes are set to zero. ",
                                 "gravity.signal", TRUE);
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "use-met-zero");
    cpl_parameter_disable (p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append (self, p);
    
    return CPL_ERROR_NONE;
}

cpl_error_code gravi_parameter_add_compute_vis (cpl_parameterlist *self, int isCalib)
{
    cpl_ensure_code (self, CPL_ERROR_NULL_INPUT);
    cpl_parameter *p;
    
    /* Max-frame */
	p = cpl_parameter_new_value ("gravity.vis.max-frame", CPL_TYPE_INT,
                                 "Maximum number of frames to integrate \n "
                                 "coherently into an OIFITS entry",
                                 "gravity.vis", 10000);
	cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "max-frame");
	cpl_parameter_disable (p, CPL_PARAMETER_MODE_ENV);
	cpl_parameterlist_append (self, p);
    
    /* Force same time for all baselines */
	p = cpl_parameter_new_value ("gravity.vis.force-same-time", CPL_TYPE_BOOL,
                                 "Force all baseline/quantities to have\n "
                                 "strickly the same TIME and MJD columns",
                                 "gravity.vis", FALSE);
	cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "force-same-time");
	cpl_parameter_disable (p, CPL_PARAMETER_MODE_ENV);
	cpl_parameterlist_append (self, p);
    
    /* Debias SC VIS2 */
	p = cpl_parameter_new_value ("gravity.vis.debias-sc", CPL_TYPE_BOOL,
                                 "Subtract the V2 bias from SC",
                                 "gravity.vis", TRUE);
	cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "debias-sc");
	cpl_parameter_disable (p, CPL_PARAMETER_MODE_ENV);
	cpl_parameterlist_append (self, p);

    /* Debias FT VIS2 */
	p = cpl_parameter_new_value ("gravity.vis.debias-ft", CPL_TYPE_BOOL,
                                 "Subtract the V2 bias from FT",
                                 "gravity.vis", TRUE);
	cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "debias-ft");
	cpl_parameter_disable (p, CPL_PARAMETER_MODE_ENV);
	cpl_parameterlist_append (self, p);

	/* Number of bootstrap */
	p = cpl_parameter_new_value ("gravity.vis.nboot", CPL_TYPE_INT,
                                 "Number of bootstraps to compute error (1..100)",
                                 "gravity.vis", isCalib ? 1 : 20);
	cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "nboot");
	cpl_parameter_disable (p, CPL_PARAMETER_MODE_ENV);
	cpl_parameterlist_append (self, p);

    /* Visibility correction */
	p = cpl_parameter_new_enum ("gravity.vis.vis-correction-sc", CPL_TYPE_STRING,
                                "Correction of SC visibility from losses due to long integration,\n "
                                "using the measured visibility losses with the FT (VFACTOR\n "
                                "and/or PFACTOR) or by forcing\n "
                                "the SC visibilities to match those of the FT (FORCE). Possible\n "
                                "choices are:",
                                "gravity.vis",
								isCalib ? "NONE" : "VFACTOR", 5, "VFACTOR", "PFACTOR",
                                "VFACTOR_PFACTOR","FORCE", "NONE");
	cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "vis-correction-sc");
	cpl_parameter_disable (p, CPL_PARAMETER_MODE_ENV);
	cpl_parameterlist_append (self, p);

    /* Phase referencing */
    p = cpl_parameter_new_enum ("gravity.vis.phase-ref-sc", CPL_TYPE_STRING,
                                "Reference phase used to integrate the SC frames.\n "
                                "Use a self-estimate of the phase, fitted by poly. (SELF_REF)\n "
                                "Use the FT phase only, interpolated in lbd (PHASE_REF)\n "
                                "Use the FT+MET-SEP.UV phase (IMAGING_REF).",
                                "gravity.vis", "AUTO", 5,
                                "SELF_REF","PHASE_REF","IMAGING_REF","AUTO","NONE");
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "phase-ref-sc");
    cpl_parameter_disable (p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append (self, p);
  
    /* Phase cleaning in the final OIFITS */
    p = cpl_parameter_new_enum ("gravity.vis.output-phase-sc", CPL_TYPE_STRING,
                                "With DIFFERENTIAL, the mean group-delay and mean\n "
                                "phases are removed from the output VISPHI in the\n "
                                "final OIFITS file. With ABSOLUTE, the VISPHI is\n "
                                "kept unmodified.",
                                "gravity.vis", "AUTO", 3,
                                "DIFFERENTIAL","ABSOLUTE","AUTO");
    cpl_parameter_set_alias (p, CPL_PARAMETER_MODE_CLI, "output-phase-sc");
    cpl_parameter_disable (p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append (self, p);
    
    return CPL_ERROR_NONE;
}

cpl_error_code gravi_parameter_add_image (cpl_parameterlist *self)
{
    cpl_ensure_code (self, CPL_ERROR_NULL_INPUT);
    cpl_parameter *p;

    /* Fill the parameters list */
    /* --isotropic */
/*    p = cpl_parameter_new_value("gravi.gravity_image.isotropic_option",
            CPL_TYPE_BOOL, "a flag", "gravi.gravity_image", FALSE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "isotropic");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);
*/

    /* --pixelsize */
    p = cpl_parameter_new_value("gravi.gravity_image.pixelsize",
            CPL_TYPE_DOUBLE, "size of the pixel (milliarcseconds)",
            "gravi.gravity_image", 0.2);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "pixelsize");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(self, p);

    /* --dim */
    p = cpl_parameter_new_value("gravi.gravity_image.dim",
            CPL_TYPE_INT, "number of pixels per side of the image",
            "gravi.gravity_image", 100);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "dim");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(self, p);

    /* --regul */
    p = cpl_parameter_new_value("gravi.gravity_image.regul",
            CPL_TYPE_STRING, "name of regularization method",
            "gravi.gravity_image", "totvar");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "regul");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(self, p);

    /* --regul_mu */
    p = cpl_parameter_new_value("gravi.gravity_image.regul_mu",
            CPL_TYPE_DOUBLE, "global regularization weight",
            "gravi.gravity_image", 1E4);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "regul_mu");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(self, p);

    /* --maxeval */
    p = cpl_parameter_new_value("gravi.gravity_image.maxeval",
            CPL_TYPE_INT, "maximum number of evaluations of the objective function",
            "gravi.gravity_image", 2000);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "maxeval");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(self, p);

    /* --timeout */
    p = cpl_parameter_new_value("gravi.gravity_image.timeout",
            CPL_TYPE_DOUBLE, "Maximum execution time of Mira process (s)",
            "gravi.gravity_image", 60.);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "timeout");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(self, p);


    return CPL_ERROR_NONE;
}

/*---------------------------------------------------------------------------*/

cpl_frameset * gravi_frameset_extract (cpl_frameset * frameset,
									   const char ** frame_tags,
									   int nb_tags)
{
    cpl_ensure (frameset,   CPL_ERROR_NULL_INPUT, NULL);
    cpl_ensure (frame_tags, CPL_ERROR_NULL_INPUT, NULL);
	cpl_ensure (nb_tags>0,  CPL_ERROR_ILLEGAL_INPUT, NULL);

	int nb_frame = cpl_frameset_get_size (frameset);
	cpl_frameset * output_frameset = cpl_frameset_new();

	/* Loop on frames in the frameset */
	for (int i = 0; i < nb_frame; i++){
	  
		cpl_frame * frame = cpl_frameset_get_position (frameset, i);
		const char * frame_tag = cpl_frame_get_tag (frame) ;

		/* Loop on requested tags */
		for (int j = 0; j < nb_tags; j++) {
		  if (strcmp(frame_tag, frame_tags[j]) == 0) {
			cpl_frameset_insert (output_frameset, cpl_frame_duplicate(frame));
			break;
		  }
		}
		
	} /* End loop on frames*/

	return output_frameset;
}

cpl_frameset * gravi_frameset_extract_p2vm_data (cpl_frameset * frameset) {
  const char *tags[] = {GRAVI_P2VM_RAW};
  return gravi_frameset_extract (frameset, tags, 1);
}
cpl_frameset * gravi_frameset_extract_disp_data (cpl_frameset * frameset) {
  const char *tags[] = {GRAVI_DISP_RAW};
  return gravi_frameset_extract (frameset, tags, 1);
}
cpl_frameset * gravi_frameset_extract_dark_data (cpl_frameset * frameset) {
  const char *tags[] = {GRAVI_DARK_RAW};
  return gravi_frameset_extract (frameset, tags, 1);
}
cpl_frameset * gravi_frameset_extract_flat_data (cpl_frameset * frameset) {
  const char *tags[] = {GRAVI_FLAT_RAW};
  return gravi_frameset_extract (frameset, tags, 1);
}
cpl_frameset * gravi_frameset_extract_diamcat_map (cpl_frameset * frameset) {
  const char *tags[] = {GRAVI_DIAMETER_CAT};
  return gravi_frameset_extract (frameset, tags, 1);
}
cpl_frameset * gravi_frameset_extract_fringe_data (cpl_frameset * frameset) {
  const char *tags[] = {GRAVI_SINGLE_CALIB_RAW, GRAVI_SINGLE_SCIENCE_RAW, GRAVI_DUAL_CALIB_RAW, GRAVI_DUAL_SCIENCE_RAW};
  return gravi_frameset_extract (frameset, tags, 4);
}
cpl_frameset * gravi_frameset_extract_p2vmred_data (cpl_frameset * frameset) {
  const char *tags[] = {GRAVI_P2VMRED_SINGLE_CALIB, GRAVI_P2VMRED_SINGLE_SCIENCE,
						GRAVI_P2VMRED_DUAL_CALIB, GRAVI_P2VMRED_DUAL_SCIENCE};
  return gravi_frameset_extract (frameset, tags, 4);
}
cpl_frameset * gravi_frameset_extract_piezotf_data (cpl_frameset * frameset) {
  const char *tags[] = {GRAVI_PIEZOTF_RAW};
  return gravi_frameset_extract (frameset, tags, 1);
}
cpl_frameset * gravi_frameset_extract_sky_data (cpl_frameset * frameset) {
  const char *tags[] = {GRAVI_DUAL_SKY_RAW, GRAVI_SINGLE_SKY_RAW};
  return gravi_frameset_extract (frameset, tags, 2);
}
cpl_frameset * gravi_frameset_extract_wave_data (cpl_frameset * frameset) {
  const char *tags[] = {GRAVI_WAVE_RAW};
  return gravi_frameset_extract (frameset, tags, 1);
}
cpl_frameset * gravi_frameset_extract_wavesc_data (cpl_frameset * frameset) {
  const char *tags[] = {GRAVI_WAVESC_RAW};
  return gravi_frameset_extract (frameset, tags, 1);
}
cpl_frameset * gravi_frameset_extract_dispvis_data (cpl_frameset * frameset) {
  const char *tags[] = {GRAVI_DISP_VIS};
  return gravi_frameset_extract (frameset, tags, 1);
}
cpl_frameset * gravi_frameset_extract_disp_map (cpl_frameset * frameset) {
  const char *tags[] = {GRAVI_DISP_MODEL};
  return gravi_frameset_extract (frameset, tags, 1);
}
cpl_frameset * gravi_frameset_extract_met_pos (cpl_frameset * frameset) {
  const char *tags[] = {GRAVI_DIODE_POSITION};
  return gravi_frameset_extract (frameset, tags, 1);
}
cpl_frameset * gravi_frameset_extract_wavelamp_map (cpl_frameset * frameset) {
  const char *tags[] = {GRAVI_WAVELAMP_MAP};
  return gravi_frameset_extract (frameset, tags, 1);
}
cpl_frameset * gravi_frameset_extract_wavelamp_data (cpl_frameset * frameset) {
  const char *tags[] = {GRAVI_WAVELAMP_RAW};
  return gravi_frameset_extract (frameset, tags, 1);
}
cpl_frameset * gravi_frameset_extract_tf_calib (cpl_frameset * frameset) {
  const char *tags[] = {GRAVI_TF_SINGLE_CALIB, GRAVI_TF_DUAL_CALIB, GRAVI_ZP_CAL};
  return gravi_frameset_extract (frameset, tags, 3);
}
cpl_frameset * gravi_frameset_extract_vis_calib (cpl_frameset * frameset) {
  const char *tags[] = {GRAVI_VIS_SINGLE_CALIB, GRAVI_VIS_DUAL_CALIB};
  return gravi_frameset_extract (frameset, tags, 2);
}
cpl_frameset * gravi_frameset_extract_vis_science (cpl_frameset * frameset) {
  const char *tags[] = {GRAVI_VIS_SINGLE_SCIENCE, GRAVI_VIS_DUAL_SCIENCE};
  return gravi_frameset_extract (frameset, tags, 2);
}
cpl_frameset * gravi_frameset_extract_p2vm_map (cpl_frameset * frameset) {
  const char *tags[] = {GRAVI_P2VM_MAP};
  return gravi_frameset_extract (frameset, tags, 1);
}
cpl_frameset * gravi_frameset_extract_flat_map (cpl_frameset * frameset){
  const char *tags[] = {GRAVI_FLAT_MAP};
  return gravi_frameset_extract (frameset, tags, 1);
}
cpl_frameset * gravi_frameset_extract_dark_map (cpl_frameset * frameset) {
  const char *tags[] = {GRAVI_DARK_MAP};
  return gravi_frameset_extract (frameset, tags, 1);
}
cpl_frameset * gravi_frameset_extract_wave_map (cpl_frameset * frameset) {
  const char *tags[] = {GRAVI_WAVE_MAP};
  return gravi_frameset_extract (frameset, tags, 1);
}
cpl_frameset * gravi_frameset_extract_bad_map (cpl_frameset * frameset) {
  const char *tags[] = {GRAVI_BAD_MAP};
  return gravi_frameset_extract (frameset, tags, 1);
}
cpl_frameset * gravi_frameset_extract_biasmask_map (cpl_frameset * frameset) {
  const char *tags[] = {GRAVI_BIASMASK_MAP};
  return gravi_frameset_extract (frameset, tags, 1);
}
cpl_frameset * gravi_frameset_extract_eop_map (cpl_frameset * frameset) {
  const char *tags[] = {GRAVI_EOP_MAP};
  return gravi_frameset_extract (frameset, tags, 1);
}
cpl_frameset * gravi_frameset_extract_patch (cpl_frameset * frameset) {
  const char *tags[] = {GRAVI_KEY_PATCH};
  return gravi_frameset_extract (frameset, tags, 1);
}

/*---------------------------------------------------------------------------*/
/* 
 * Get the parameter from the list. Provide a default in case the parameter
 * is NOT in the list.
 */
/*---------------------------------------------------------------------------*/

double gravi_param_get_double_default (const cpl_parameterlist * parlist, const char * name, double def)
{
  const cpl_parameter * tmp = cpl_parameterlist_find_const(parlist, name);
  
  if (tmp) {
    return cpl_parameter_get_double (tmp);
  } else {
    cpl_msg_info (cpl_func, "Could not find the parameter '%s':, use %f", name, def);
    return def;
  }
}

int gravi_param_get_int_default (const cpl_parameterlist * parlist, const char * name, int def)
{
  const cpl_parameter * tmp = cpl_parameterlist_find_const(parlist, name);
  
  if (tmp) {
    return cpl_parameter_get_int (tmp);
  } else {
    cpl_msg_info (cpl_func, "Could not find the parameter '%s': use %i", name, def);
    return def;
  }
}

int gravi_param_get_bool_default (const cpl_parameterlist * parlist, const char * name, int def)
{
  const cpl_parameter * tmp = cpl_parameterlist_find_const(parlist, name);
  
  if (tmp) {
    return cpl_parameter_get_bool (tmp);
  } else {
    cpl_msg_info (cpl_func, "Could not find the boolean parameter '%s': use %s", name, (def==0?"FALSE":"TRUE"));
    return def;
  }
}

const char * gravi_param_get_string_default (const cpl_parameterlist * parlist, const char * name, const char * def)
{
  const cpl_parameter * tmp = cpl_parameterlist_find_const(parlist, name);
  
  if (tmp) {
    return cpl_parameter_get_string (tmp);
  } else {
    cpl_msg_info (cpl_func, "Could not find the string parameter '%s': use %s", name, def);
    return def;
  }
}

double gravi_param_get_double (const cpl_parameterlist * parlist, const char * name)
{
  const cpl_parameter * tmp = cpl_parameterlist_find_const(parlist, name);
  double def = 0.0;
  
  if (tmp) {
    return cpl_parameter_get_double (tmp);
  } else {
    cpl_error_set_message (cpl_func, CPL_ERROR_ILLEGAL_INPUT, "Could not find the parameter '%s':, use %f", name, def);
    return def;
  }
}

int gravi_param_get_int (const cpl_parameterlist * parlist, const char * name)
{
  const cpl_parameter * tmp = cpl_parameterlist_find_const(parlist, name);
  int def = 0;

  if (tmp) {
    return cpl_parameter_get_int (tmp);
  } else {
    cpl_error_set_message (cpl_func, CPL_ERROR_ILLEGAL_INPUT, "Could not find the parameter '%s': use %i", name, def);
    return def;
  }
}

int gravi_param_get_bool (const cpl_parameterlist * parlist, const char * name)
{
  const cpl_parameter * tmp = cpl_parameterlist_find_const(parlist, name);
  int def = 0;
  
  if (tmp) {
    return cpl_parameter_get_bool (tmp);
  } else {
    cpl_error_set_message (cpl_func, CPL_ERROR_ILLEGAL_INPUT, "Could not find the boolean parameter '%s': use %s", name, (def==0?"FALSE":"TRUE"));
    return def;
  }
}

const char * gravi_param_get_string (const cpl_parameterlist * parlist, const char * name)
{
  const cpl_parameter * tmp = cpl_parameterlist_find_const(parlist, name);
  
  if (tmp) {
    return cpl_parameter_get_string (tmp);
  } else {
    cpl_error_set_message (cpl_func, CPL_ERROR_ILLEGAL_INPUT, "Could not find the string parameter '%s': use UNKNOWN", name);
    return "UNKNOWN";
  }
}

cpl_error_code gravi_check_frameset (cpl_frameset *frameset, const char * tag, int min, int max)
{
  int flag = 0;
  int nf = cpl_frameset_count_tags (frameset, tag);
  char * msg = cpl_sprintf ("Need %i<#<%i '%s' in frameset (%i provided)", min, max, tag, nf);
  
  if (nf < min || nf > max) {
	cpl_msg_error (cpl_func, "%s",msg);
    cpl_error_set_message (cpl_func, CPL_ERROR_ILLEGAL_INPUT, "%s", msg);
	flag ++;
  } else {
	cpl_msg_info (cpl_func, "%s", msg);
  }

  cpl_free (msg);
  return (flag) ? CPL_ERROR_ILLEGAL_INPUT : CPL_ERROR_NONE;
}

/**@}*/
