/*
 * This file is part of the GRAVI Pipeline
 * Copyright (C) 2022 European Southern Observatory
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

/*-----------------------------------------------------------------------------
                                   Includes
 -----------------------------------------------------------------------------*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include <string.h>
#include "gravi_utils.h"
#include "gravi_idp.h"
#include "gravi_pfits.h"
#include "gravi_dfs.h"
#include "gravi_cpl.h"
#include "gravi_vis.h"

/*-----------------------------------------------------------------------------*/
/**
 * @brief Create IDP keywords to satisfy standard
 * @param oi_vis2_SC_allpol Merged (all pol) squared science visibilities
 * @param oi_T3_SC_allpol Merged (all pol) science T3
 * @param oi_wave_SC_allpol Merged (all pol) science wavelengths
 * @param header    The input header
 * @return A new cpl_parameterlist with the keywords
 */
/*-----------------------------------------------------------------------------*/
cpl_propertylist * gravi_idp_compute (gravi_data * vis_data,
                                      cpl_propertylist * header,
                                      cpl_frameset * frameset)
{
    cpl_propertylist * idp_plist = cpl_propertylist_new ();
    int nbase = 6;

    int npol_sc = gravi_pfits_get_pola_num (header, GRAVI_SC);
    char qc_name[100];

    /* Scratch tables to merge polarization tables */
    cpl_table * oi_vis2_SC_allpol = NULL;
    cpl_table * oi_T3_SC_allpol = NULL;
    cpl_table * oi_wave_SC_allpol = NULL;
 
    /* There are products like the p2vmred that do not have all columns */
    if (gravi_data_has_extension(vis_data, GRAVI_OI_VIS2_EXT))
    {
        double vis2err = 0;
        double visphierr = 0;
        double t3phierr = 0;
        oi_vis2_SC_allpol = cpl_table_duplicate(gravi_data_get_oi_vis2 (vis_data, GRAVI_SC, 0, npol_sc));
        oi_T3_SC_allpol = cpl_table_duplicate(gravi_data_get_oi_t3 (vis_data, GRAVI_SC, 0, npol_sc));
        oi_wave_SC_allpol = cpl_table_duplicate(gravi_data_get_oi_wave (vis_data, GRAVI_SC, 0, npol_sc));
        for (int pol = 1; pol < npol_sc; pol++) {
            cpl_table_insert(oi_vis2_SC_allpol, gravi_data_get_oi_vis2 (vis_data, GRAVI_SC, pol, npol_sc), cpl_table_get_nrow(oi_vis2_SC_allpol));
            cpl_table_insert(oi_T3_SC_allpol, gravi_data_get_oi_t3 (vis_data, GRAVI_SC, pol, npol_sc), cpl_table_get_nrow(oi_T3_SC_allpol));
            cpl_table_insert(oi_wave_SC_allpol, gravi_data_get_oi_wave (vis_data, GRAVI_SC, pol, npol_sc), cpl_table_get_nrow(oi_wave_SC_allpol));
        }
        
        for (int pol = 0; pol < npol_sc; pol++) {
            double this_pol_vis2err = gravi_table_get_column_flagged_mean(gravi_data_get_oi_vis2 (vis_data, GRAVI_SC, pol, npol_sc), "VIS2ERR");
            vis2err += this_pol_vis2err * this_pol_vis2err;
            double this_pol_visphierr = gravi_table_get_column_flagged_mean(gravi_data_get_oi_vis (vis_data, GRAVI_SC, pol, npol_sc), "VISPHIERR");
            visphierr += this_pol_visphierr * this_pol_visphierr;
            double this_pol_t3phierr = gravi_table_get_column_flagged_mean(gravi_data_get_oi_t3 (vis_data, GRAVI_SC, pol, npol_sc), "T3PHIERR");
            t3phierr += this_pol_t3phierr * this_pol_t3phierr;
        }

        /* Compute QC values that aggregate on polarizations */
        sprintf (qc_name, "VIS2ERR");
        cpl_propertylist_update_double (idp_plist, qc_name, sqrt(vis2err));
        cpl_propertylist_set_comment (idp_plist, qc_name, "Representative squared visibility error [%]");

        sprintf (qc_name, "VISPHERR");
        cpl_propertylist_update_double (idp_plist, qc_name, sqrt(visphierr));
        cpl_propertylist_set_comment (idp_plist, qc_name, "Representative visibility phase error [%]");

        sprintf (qc_name, "T3PHIERR");
        cpl_propertylist_update_double (idp_plist, qc_name, sqrt(t3phierr));
        cpl_propertylist_set_comment (idp_plist, qc_name, "Representative closure phase error [deg]");

        double min_uvcoord, max_uvcoord;
        gravi_data_get_minmax_uvcoord(oi_vis2_SC_allpol, &min_uvcoord, &max_uvcoord);

        double min_eff_wave = cpl_table_get_column_min(oi_wave_SC_allpol, "EFF_WAVE");
        double max_eff_wave = cpl_table_get_column_max(oi_wave_SC_allpol, "EFF_WAVE");

        double base_max = max_uvcoord / min_eff_wave;
        double base_min = min_uvcoord / max_eff_wave;
        if(isnan(base_max) || isinf(base_max))
            base_max = 0;
        if(isnan(base_min) || isinf(base_min))
            base_min = 0;

        sprintf (qc_name, "BASE_MAX");
        cpl_propertylist_update_double (idp_plist, qc_name, base_max);
        cpl_propertylist_set_comment (idp_plist, qc_name, "Maximum baseline / Minimum effective wavelenth");

        sprintf (qc_name, "BASE_MIN");
        cpl_propertylist_update_double (idp_plist, qc_name, base_min);
        cpl_propertylist_set_comment (idp_plist, qc_name, "Minimum baseline / Maximum effective wavelenth");

        /* The rows in oi_wave_SC_allpol contain each wavelenght npol_sc times,
           since it has been aggregated. Therefore dividing by npol_sc times */
        sprintf (qc_name, "NUM_CHAN");
        cpl_propertylist_update_int (idp_plist, qc_name, cpl_table_get_nrow(oi_wave_SC_allpol) / npol_sc );
        cpl_propertylist_set_comment (idp_plist, qc_name, "Number of wavelength channels");

        sprintf (qc_name, "WAVELMAX");
        double avg_band = cpl_table_get_column_mean(oi_wave_SC_allpol, "EFF_BAND");
        cpl_propertylist_update_double (idp_plist, qc_name, (max_eff_wave - avg_band / 2.) * 1e9);
        cpl_propertylist_set_comment (idp_plist, qc_name, "[nm] Maximum wavelength");

        sprintf (qc_name, "WAVELMIN");
        cpl_propertylist_update_double (idp_plist, qc_name, (min_eff_wave + avg_band / 2.) * 1e9);
        cpl_propertylist_set_comment (idp_plist, qc_name, "[nm] Minimum wavelength");

        cpl_table_duplicate_column(oi_wave_SC_allpol, "SPEC_RES", oi_wave_SC_allpol, "EFF_WAVE");
        cpl_table_divide_columns(oi_wave_SC_allpol,"SPEC_RES", "EFF_BAND");

        sprintf (qc_name, "SPEC_RES");
        cpl_propertylist_update_double (idp_plist, qc_name, cpl_table_get_column_mean(oi_wave_SC_allpol, "SPEC_RES"));
        cpl_propertylist_set_comment (idp_plist, qc_name, "Spectral resolution");

        /* This is the mean INT_TIME, which includes duplicated entries due to the several polarizations.
           Since it is a mean, the final value should be the same */
        sprintf (qc_name, "EXPTIME");
        double mean_int_time = gravi_table_get_column_flagged_mean(oi_vis2_SC_allpol, "INT_TIME");
        cpl_propertylist_update_double (idp_plist, qc_name, mean_int_time);
        cpl_propertylist_set_comment (idp_plist, qc_name, "Exposure time");
         /* This is the sum of all INT_TIME divided by the number of baselines. The mean is multiplied
           by the number of rows and then divided by the number of polarizations since they have been
           aggregated. Finally divided by the number of baselines as specified in PIPE-9900 */
        if(!cpl_propertylist_has(header, "TEXPTIME"))
        {
            cpl_size nsets = cpl_table_get_nrow(oi_vis2_SC_allpol) / GRAVI_NBASE ;
            double texptime = 0;
            for (cpl_size  set = 0 ; set < nsets ; set++)
            {
                cpl_table * int_time_this_set = cpl_table_extract(oi_vis2_SC_allpol, GRAVI_NBASE * set, GRAVI_NBASE);
                double max_int_time = gravi_table_get_column_flagged_max(int_time_this_set, "INT_TIME");
                texptime += max_int_time;
                cpl_table_delete(int_time_this_set);
            }
            
            sprintf (qc_name, "TEXPTIME");
            cpl_propertylist_update_double (idp_plist, qc_name, texptime);
            cpl_propertylist_set_comment (idp_plist, qc_name, "Total exposure time");
        }
        else
        {
            cpl_propertylist_update_double (idp_plist, "TEXPTIME",
                    cpl_propertylist_get_double(header, "TEXPTIME") );
            cpl_propertylist_set_comment (idp_plist, "TEXPTIME", "Total exposure time");
        }
    }

    /* PRODCATG */
    cpl_propertylist_update_string (idp_plist, "PRODCATG", "SCIENCE.VISIBILITY.UNCALIBRATED");
    cpl_propertylist_set_comment (idp_plist, "PRODCATG", "Data product category");

    /* MJD-OBS */
    double mjd_obs_first = DBL_MAX;
    if(frameset != NULL)
    {
        const cpl_frame *frame;
        cpl_frameset * science_frames = gravi_frameset_extract_fringe_data(frameset);
        cpl_frameset_iterator *it = cpl_frameset_iterator_new(science_frames);
        while ((frame = cpl_frameset_iterator_get(it)) != NULL) {
            cpl_propertylist * this_frame_header = cpl_propertylist_load(cpl_frame_get_filename(frame), 0);
            double mjd_obs = gravi_pfits_get_mjd(this_frame_header);
            if (mjd_obs < mjd_obs_first)
                mjd_obs_first = mjd_obs;
            cpl_frameset_iterator_advance(it, 1);
            cpl_propertylist_delete(this_frame_header);
        }
        cpl_frameset_delete(science_frames);
        cpl_frameset_iterator_delete(it);
    }
    if (mjd_obs_first == DBL_MAX)
        mjd_obs_first = gravi_pfits_get_mjd(header);
    cpl_propertylist_update_double (idp_plist, "MJD-OBS", mjd_obs_first);
    cpl_propertylist_set_comment (idp_plist, "MJD-OBS", "Start of observation");
    /* Delete MJD-OBS from the main header, since it will take precedence over what is computed here*/
    cpl_propertylist_erase (header, "MJD-OBS");

    /* MJD-END */
    double mjd_obs_last = 0;
    double exptime_last = 0;
    if(frameset != NULL)
    {
        const cpl_frame *frame;
        cpl_frameset * science_frames = gravi_frameset_extract_fringe_data(frameset);
        cpl_frameset_iterator *it = cpl_frameset_iterator_new(science_frames);
        while ((frame = cpl_frameset_iterator_get(it)) != NULL) {
            cpl_propertylist * this_frame_header = cpl_propertylist_load(cpl_frame_get_filename(frame), 0);
            double mjd_obs = gravi_pfits_get_mjd(this_frame_header);
            if (mjd_obs > mjd_obs_last)
            {
                mjd_obs_last = mjd_obs;
                exptime_last =  cpl_propertylist_get_double(this_frame_header, "EXPTIME");
            }   
            cpl_frameset_iterator_advance(it, 1);
            cpl_propertylist_delete(this_frame_header);
        }
        cpl_frameset_delete(science_frames);
        cpl_frameset_iterator_delete(it);
    }
    if (mjd_obs_last == 0)
        mjd_obs_last = gravi_pfits_get_mjd(header);
    if (exptime_last == 0)
    {
        if ( cpl_propertylist_has(header, "EXPTIME") )
            exptime_last = cpl_propertylist_get_double(header, "EXPTIME");
        else
            exptime_last = cpl_propertylist_get_double(idp_plist, "EXPTIME");
    }

    cpl_propertylist_update_double (idp_plist, "MJD-END",
        mjd_obs_last + exptime_last / 86400.);
    cpl_propertylist_set_comment (idp_plist, "MJD-END", "End of observation");

    /* OBID */
    cpl_propertylist_update_int (idp_plist, "OBID1",
            cpl_propertylist_get_int(header, "ESO OBS ID"));
    cpl_propertylist_set_comment (idp_plist, "OBID1", "Obseration Block ID");

    /* NCOMBINE */
    if(frameset != NULL)
    {
        cpl_frameset * science_frames = gravi_frameset_extract_fringe_data(frameset);
        cpl_size nscience = cpl_frameset_get_size(science_frames);
        cpl_frameset_delete(science_frames);
        if(!cpl_propertylist_has(header, "NCOMBINE"))
        {
            if (nscience != 0)
            {
                cpl_propertylist_update_int (idp_plist, "NCOMBINE", nscience);
                cpl_propertylist_set_comment (idp_plist, "NCOMBINE", "Number of raw science combined");
            }
        }
        else
        {
            cpl_propertylist_update_int (idp_plist, "NCOMBINE",
                    cpl_propertylist_get_int(header, "NCOMBINE") );
            cpl_propertylist_set_comment (idp_plist, "NCOMBINE", "Number of raw science combined");
        }
    }
    /* OBSTECH */
    // Only create OBSTECH if it does not exist yet.
    // This is needed for For gravity_viscal which starts from
    // products and does not have a ESO DPR TECH anymore
    if(cpl_propertylist_has(header, "ESO DPR TECH"))
    {
        cpl_propertylist_update_string (idp_plist, "OBSTECH",
                cpl_propertylist_get_string(header, "ESO DPR TECH") );
        cpl_propertylist_set_comment (idp_plist, "OBSTECH", "Observation technique");
    }
    else if(cpl_propertylist_has(header, "OBSTECH"))
    {
        cpl_propertylist_update_string (idp_plist, "OBSTECH",
                cpl_propertylist_get_string(header, "OBSTECH") );
        cpl_propertylist_set_comment (idp_plist, "OBSTECH", "Observation technique");
    }

    /* SPECSYS */
    cpl_propertylist_update_string (idp_plist, "SPECSYS", "TOPOCENT");
    cpl_propertylist_set_comment (idp_plist, "SPECSYS", "Frame of reference for spectral coordinates");

    /* TIMESYS */
    cpl_propertylist_update_string (idp_plist, "TIMESYS", "UTC");
    cpl_propertylist_set_comment (idp_plist, "TIMESYS", "Time system");

    /* SPEC_ERR */
    /* According to https://jira.eso.org/browse/PIPE-9900 this
       is hard-coded depending on resolution */
    double spec_err = 0;
    const char * resolution = gravi_pfits_get_resolution (header);
    if ( !strcmp (resolution, "HIGH") )
        spec_err = 0.28;
    if ( !strcmp (resolution, "MED") )
        spec_err = 2.2;
    if ( !strcmp (resolution, "LOW") )
        spec_err = 50;
    cpl_propertylist_update_double (idp_plist, "SPEC_ERR", spec_err);
    cpl_propertylist_set_comment (idp_plist, "SPEC_ERR", "Statistical error in spectral coordinate");

    /* SPEC_SYE */
    /* Hard-coded to the values in Sanchez-Bermudez et al. 2017 */
    cpl_propertylist_update_double (idp_plist, "SPEC_SYE", 0.1);
    cpl_propertylist_set_comment (idp_plist, "SPEC_SYE", "Systematic error in spectral coordinate");

    /* PROV keywords */
    if(frameset != NULL)
    {
        const cpl_frame *frame;
        size_t i_prov = 1;
        char prov_keyword[8];
        cpl_frameset_iterator *it = cpl_frameset_iterator_new(frameset);
        while ((frame = cpl_frameset_iterator_get(it)) != NULL) {
            if (strcmp(cpl_frame_get_tag(frame), GRAVI_SINGLE_SCIENCE_RAW) == 0 || 
                strcmp(cpl_frame_get_tag(frame), GRAVI_DUAL_SCIENCE_RAW) == 0)
            {
                snprintf(prov_keyword, 7, "PROV%zu",i_prov);
                const char * filename = cpl_frame_get_filename(frame);
                const char * filename_no_path = strrchr(filename, '/');
                if (filename_no_path == NULL)
                    filename_no_path = filename;
                else
                    filename_no_path += 1;
                cpl_propertylist_update_string(idp_plist, prov_keyword, filename_no_path);
                i_prov++;
            }
            cpl_frameset_iterator_advance(it, 1);
        }
        cpl_frameset_iterator_delete(it);
    }
        
    if (gravi_data_has_extension(vis_data, GRAVI_OI_ARRAY_EXT))
    {
        cpl_table * oi_array = gravi_data_get_table (vis_data, GRAVI_OI_ARRAY_EXT);
        cpl_table_new_column(oi_array, "FOV", CPL_TYPE_DOUBLE);
        cpl_table_new_column(oi_array, "FOVTYPE", CPL_TYPE_STRING);
        double fov;
        const char * telname = gravi_conf_get_telname (0, header);
        if (telname == NULL) {
            cpl_msg_warning(cpl_func, "Cannot get TELNAME, FOV is not determined");
            fov = 0.0;
        } else {
            if (telname[0] == 'U')
                fov = 0.03;  // Hard-coded UT FOV
            else
                fov = 0.14;  // Hard-coded AT FOV
        }

        cpl_table_fill_column_window_double(oi_array, "FOV",  0, cpl_table_get_nrow(oi_array), fov);
        cpl_table_fill_column_window_string(oi_array, "FOVTYPE",  0, cpl_table_get_nrow(oi_array),"RADIUS");
    }

    /* Delete scratch tables */
    cpl_table_delete(oi_vis2_SC_allpol);
    cpl_table_delete(oi_T3_SC_allpol);
    cpl_table_delete(oi_wave_SC_allpol);

    return idp_plist;
}

