/*
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
 * @defgroup gravi_eop  Earth Orientation Parameters
 *
 * This module implements the function link the the Earth Orientation Parameters.
 * It contains the function called by the recipe @c gravity_eop to generate the static
 * calibration file : @c gravity_eop_download_finals2000A() and
 * @c gravity_eop_data_totable()
 *
 * It also implements the computation the UV coordinates making use of this EOP
 * calibration file : @c gravi_compute_pointing_uv()
 *
 */
/**@{*/

/*-----------------------------------------------------------------------------
                                   Includes
 -----------------------------------------------------------------------------*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <unistd.h>

#include <erfa.h>

#include "gravi_eop.h"
#include "gravi_cpl.h"
#include "gravi_pfits.h"

#include "gravi_utils.h"

/*-----------------------------------------------------------------------------
                             Private prototypes and define
 -----------------------------------------------------------------------------*/


cpl_error_code gravi_eop_interpolate (cpl_size n, double *mjd,
                                      double *pmx, double *pmy,
                                      double *dut, cpl_table * eop_table,
                                      cpl_propertylist * header);

void eraAtboq (double rc, double dc, eraASTROM *astrom, double enuob[3]);
void eraAtcoq (double rc, double dc, double pmr, double pmd, double px,
               double rv, eraASTROM *astrom, double enuob[3]);
void dtp2s (double xi, double eta, double raz,
            double decz, double *ra, double *dec);
void rotate_vector (double in[3], double angle, double axis[3], double out[3]);
void difference (double x[3], double y[3], double z[3]);
void multiply (double xyz[3], double factor);
void normalize (double xyz[3]);
void cross (double x[3], double y[3], double z[3]);

/*-----------------------------------------------------------------------------
                             Functions code
 -----------------------------------------------------------------------------*/

/* 
 * Private function to interpolate EOP parameter from EOP_PARAM
 * and manipulate coordinates and 3D vectors
 */

cpl_error_code gravi_eop_interpolate (cpl_size n, double *mjd,
                                      double *pmx, double *pmy,
                                      double *dut,
                                      cpl_table * eop_table,
                                      cpl_propertylist * header)
{  
  gravi_msg_function_start(1);
  cpl_ensure_code (eop_table, CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (header, CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (n>0, CPL_ERROR_ILLEGAL_INPUT);
  cpl_ensure_code (mjd, CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (pmx, CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (pmy, CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (dut, CPL_ERROR_NULL_INPUT);

  /* Wrap inputs in vectors */
  cpl_vector *vmjd = cpl_vector_wrap (n, mjd);
  cpl_vector *vpmx = cpl_vector_wrap (n, pmx);
  cpl_vector *vpmy = cpl_vector_wrap (n, pmy);
  cpl_vector *vdut = cpl_vector_wrap (n, dut);

  /* Check validity of input mjd with table */
  if (cpl_vector_get_min (vmjd) < cpl_table_get_column_min (eop_table, "MJD") ||
      cpl_vector_get_max (vmjd) > cpl_table_get_column_max (eop_table, "MJD"))
  {
    cpl_msg_warning (cpl_func, "Some MJD are outside the EOP_PARAM range (MJD=%.2f, %.2f..%.2f). Use EOP and DUT=0.0s",
                     cpl_vector_get_mean (vmjd),
                     cpl_table_get_column_min (eop_table, "MJD"),
                     cpl_table_get_column_max (eop_table, "MJD"));
    cpl_vector_unwrap (vmjd);
    cpl_vector_unwrap (vpmx);
    cpl_vector_unwrap (vpmy);
    cpl_vector_unwrap (vdut);
    gravi_msg_function_exit(1);
    return CPL_ERROR_NONE;
  }

  /* Check if within the prediction range */
  if(cpl_vector_get_max (vmjd) > cpl_propertylist_get_double (header, "ESO QC EOP MJD LAST FINAL"))
  {
    cpl_msg_warning (cpl_func, "Some MJD are inside the EOP_PARAM prediction range...");
  }

  /* Wrap output vectors */
  cpl_size nrow = cpl_table_get_nrow (eop_table);
  cpl_vector * eop_mjd = cpl_vector_wrap (nrow, cpl_table_get_data_double (eop_table, "MJD"));
  cpl_vector * eop_pmx = cpl_vector_wrap (nrow, cpl_table_get_data_double (eop_table, "PMX"));
  cpl_vector * eop_pmy = cpl_vector_wrap (nrow, cpl_table_get_data_double (eop_table, "PMY"));
  cpl_vector * eop_dut = cpl_vector_wrap (nrow, cpl_table_get_data_double (eop_table, "DUT"));
  
  /* Perform interpolations */
  cpl_bivector *mjd_pmx = cpl_bivector_wrap_vectors (vmjd, vpmx);
  cpl_bivector *eop_mjd_pmx = cpl_bivector_wrap_vectors (eop_mjd, eop_pmx);
  cpl_bivector_interpolate_linear (mjd_pmx, eop_mjd_pmx);
  cpl_bivector_unwrap_vectors (mjd_pmx);
  cpl_bivector_unwrap_vectors (eop_mjd_pmx);

  cpl_bivector *mjd_pmy = cpl_bivector_wrap_vectors (vmjd, vpmy);
  cpl_bivector *eop_mjd_pmy = cpl_bivector_wrap_vectors (eop_mjd, eop_pmy);
  cpl_bivector_interpolate_linear (mjd_pmy, eop_mjd_pmy);
  cpl_bivector_unwrap_vectors (mjd_pmy);
  cpl_bivector_unwrap_vectors (eop_mjd_pmy);

  cpl_bivector *mjd_dut = cpl_bivector_wrap_vectors (vmjd, vdut);
  cpl_bivector *eop_mjd_dut = cpl_bivector_wrap_vectors (eop_mjd, eop_dut);
  cpl_bivector_interpolate_linear (mjd_dut, eop_mjd_dut);
  cpl_bivector_unwrap_vectors (mjd_dut);
  cpl_bivector_unwrap_vectors (eop_mjd_dut);
  
  /* Print diagnostics */
  cpl_msg_info(cpl_func, "EOP averages: MJD=%.2f DUT1=%.3f PMX=%.3farcsec PMY=%.3farcsec",
               cpl_vector_get_mean (vmjd), cpl_vector_get_mean (vdut),
               cpl_vector_get_mean (vpmx), cpl_vector_get_mean (vpmy));

  /* Unwrap vectors */
  cpl_vector_unwrap (vmjd);
  cpl_vector_unwrap (vpmx);
  cpl_vector_unwrap (vpmy);
  cpl_vector_unwrap (vdut);
  cpl_vector_unwrap (eop_mjd);
  cpl_vector_unwrap (eop_pmx);
  cpl_vector_unwrap (eop_pmy);
  cpl_vector_unwrap (eop_dut);
  
  CPLCHECK_MSG ("Cannot interpolate EOP");

  gravi_msg_function_exit(1);
  return CPL_ERROR_NONE;
}

/*
 * Quick transform from BCRS to Observed
 *     rc,dc  double   ICRS right ascension at J2000.0 (radians)
 */
void eraAtboq (double rc, double dc, eraASTROM *astrom, double enuob[3])
{
    double ri, di;
    double aob, zob, hob, dob, rob;

    /* Transform from BCRS to Intermediate */
    eraAtciq (rc, dc, 0.0, 0.0, 0.0, 0.0, astrom, &ri, &di);

    /* Transform from Intermediate to Observed. */
    eraAtioq (ri, di, astrom, &aob, &zob, &hob, &dob, &rob);

    /* Convert from equatorial to cartesian */
    eraS2c(CPL_MATH_PI/2.0-aob, CPL_MATH_PI/2.0-zob, enuob);
}

void rotate_vector (double in[3], double angle, double axis[3], double out[3])
{
    double rv[3];
    double rm[3][3];
    eraSxp(angle, axis, rv);
    eraRv2m(rv, rm);
    eraRxp(rm, in, out);
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Compute the pointing directions and projected baselines
 * 
 * @param input_table    input/output data
 * @param header         input header
 * @param eop_table:     table containing the EOP pameters
 * @param eop_header:    header of EOP table
 * @param save_pointing: save (E_U,E_V,E_V,E_AZ,E_ZD) is specified
 * @param array_table:   OI_ARRAY table (optional)
 * 
 * For each DIT of the input table, compute [E_U,E_V,E_W,E_AZ,E_ZD]_Obs as the
 * transformation into Observed reference frame of the orthonormal
 * [E_U,E_V,E_W,E_AZ,E_ZD]_ICRS defined in the ICRS.
 * This way, the real-time projected baseline can be recomputed easily off-line.
 * [E_U,E_V,E_W,E_AZ,E_ZD]_Obs does not form an orthonormal basis,
 * due to the effects of precession, nutation, aberration...
 * The quantities [E_U,E_V,E_W,E_AZ,E_ZD]_Obs are stored as new columns
 * in the input table, is save_pointing is specified.
 * If the array_table is specified, the projected baseline [UCOORD,VCOORD]
 * is calculated.
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_eop_pointing_uv (cpl_table * input_table,
                                      cpl_propertylist * header,
                                      cpl_table * eop_table,
                                      cpl_propertylist * eop_header,
                                      int save_pointing,
                                      cpl_table * array_table)
{
    gravi_msg_function_start(1);
    cpl_ensure_code (input_table, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code (header,      CPL_ERROR_NULL_INPUT);
    
    /* Check and announce the optional parameters */
    if (save_pointing > 0) {
        cpl_msg_info (cpl_func, "Will save [E_U,E_V,E_W,E_AZ,E_ZD]");
    }
    int compute_uv;
    if (array_table != NULL) {
        cpl_msg_info (cpl_func, "Will compute [UCOORD,VCOORD]");
        compute_uv = 1;
    } else {
        compute_uv = 0;
    }
    
    double t_skip = 1./24/3600, mjd0 = -1.0, mjd1 = -1.0;
    cpl_msg_info (cpl_func, "Compute pointing with full ERFA every %.2f s",t_skip*24*3600.0);
    
    /* If requested, loop on bases to compute the physical 3D baseline as T2-T1
    * Note: The telescope locations STAXYZ are given in the [West,South,Up]
    *   frame (ESO convention), whereas the UVW, calculated later on, are in
    *   the [East,North,Up] frame. Hence the sign changes below on X and Y. */
    double baseline[GRAVI_NBASE][3];
    double * uCoord;
    double * vCoord;
    if (compute_uv) {
        for ( int base = 0; base < GRAVI_NBASE; base ++) {
           int tel1=0; while ( cpl_table_get (array_table, "STA_INDEX", tel1, NULL) != gravi_table_get_value (input_table, "STA_INDEX", base, 0) ) tel1++;
           int tel2=0; while ( cpl_table_get (array_table, "STA_INDEX", tel2, NULL) != gravi_table_get_value (input_table, "STA_INDEX", base, 1) ) tel2++;
           baseline[base][0] = -(gravi_table_get_value (array_table, "STAXYZ", tel2, 0) - gravi_table_get_value (array_table, "STAXYZ", tel1, 0));
           baseline[base][1] = -(gravi_table_get_value (array_table, "STAXYZ", tel2, 1) - gravi_table_get_value (array_table, "STAXYZ", tel1, 1));
           baseline[base][2] = +(gravi_table_get_value (array_table, "STAXYZ", tel2, 2) - gravi_table_get_value (array_table, "STAXYZ", tel1, 2));
        }
        uCoord = cpl_table_get_data_double (input_table, "UCOORD");
        vCoord = cpl_table_get_data_double (input_table, "VCOORD");
    }
    
    /* Read UTC as MJD */
    double *mjd = cpl_table_get_data_double (input_table, "MJD");
    double mean_mjd = cpl_table_get_column_mean (input_table, "MJD");
    
    /* If requested, create [E_U,E_V,E_V,E_AZ,E_ZD] columns */
    cpl_array ** p_u;
    cpl_array ** p_v;
    cpl_array ** p_w;
    cpl_array ** p_az;
    cpl_array ** p_zd;
    if (save_pointing) {
        cpl_msg_info (cpl_func, "Saving E_U, E_V, E_W, E_AZ, E_ZD");
        cpl_table_new_column_array (input_table, "E_U", CPL_TYPE_DOUBLE, 3);
        cpl_table_new_column_array (input_table, "E_V", CPL_TYPE_DOUBLE, 3);
        cpl_table_new_column_array (input_table, "E_W", CPL_TYPE_DOUBLE, 3);
        cpl_table_new_column_array (input_table, "E_AZ", CPL_TYPE_DOUBLE, 3);
        cpl_table_new_column_array (input_table, "E_ZD", CPL_TYPE_DOUBLE, 3);
        p_u = cpl_table_get_data_array (input_table,"E_U");
        p_v = cpl_table_get_data_array (input_table,"E_V");
        p_w = cpl_table_get_data_array (input_table,"E_W");
        p_az = cpl_table_get_data_array (input_table,"E_AZ");
        p_zd = cpl_table_get_data_array (input_table,"E_ZD");
        CPLCHECK_MSG ("Cannot create [E_U, E_V, E_W, E_AZ, E_ZD]");
    }
    
    /* Get the location of the observer from the header */
    double elev = gravi_pfits_get_geoelev (header); // Height in [m]
    double lon = gravi_pfits_get_geolon (header) * CPL_MATH_RAD_DEG; // Lat in [rad]
    double lat = gravi_pfits_get_geolat (header) * CPL_MATH_RAD_DEG; // Lon in [rad], East positive
    CPLCHECK_MSG ("Cannot get the observer location");
    
    /* If EOP are supplied, interpolate to the mean MJD */
    double dut1 = 0, pmx = 0, pmy = 0;
    if (eop_table != NULL) {
        gravi_eop_interpolate (1, &mean_mjd, &pmx, &pmy, &dut1,
                            eop_table, eop_header);
                            CPLCHECK_MSG ("Cannot interpolate");
    } else {
        cpl_msg_warning (cpl_func, "No EOP_PARAM. Use EOP and DUT=0.0s");
    }
    
    /* We use the mid-point between FT and SC to compute the UV coordinates */
    double rc  = gravi_pfits_get_mid_raep  (header); // [rad]
    double dc  = gravi_pfits_get_mid_decep (header); // [rad]
    double pmr = gravi_pfits_get_pmra  (header) * CPL_MATH_RAD_DEG / 3600.0 / cos(dc); // dRA/dt, not cos(Dec)xdRA/dt [rad/year]
    double pmd = gravi_pfits_get_pmdec (header) * CPL_MATH_RAD_DEG / 3600.0; // dDec/dt [rad/year]
    double parallax = gravi_pfits_get_plx (header); // [as]
    double sysvel = 0.0;
    CPLCHECK_MSG ("Cannot get the header data");
    
    /* Prepare the following loop computations */
    eraASTROM astrom;
    double eo;
    double eps = 10.0 / 3600.0 * CPL_MATH_RAD_DEG; // 10 arcsec has been chosen for optimal accuracy
    double rb, db;
    double eUb[3], eVb[3], eWb[3], eZb[3];
    double eWo_up[3], eWo_um[3], eWo_vp[3], eWo_vm[3];
    double eWb_up[3], eWb_um[3], eWb_vp[3], eWb_vm[3];
    double rb_up, db_up, rb_um, db_um, rb_vp, db_vp, rb_vm, db_vm;
    double eUo[3], eVo[3], eWo[3], eAZo[3], eZDo[3];
    double ez[3] = {0.0, 0.0, 1.0}; // Zenith direction in ENU frame
    double norm;
    double pressure = 0.0; // Pressure at zero to disable atmospheric refraction
    double temperature = 0.0;
    double humidity = 0.0;
    double wavelength = 0.0;
    
    /* Loop on rows. The baselines are not assumed to share the
     * same MJD since this function is called after the averaging */
    cpl_size n_row = cpl_table_get_nrow (input_table);
    for (cpl_size row = 0 ; row < n_row ; row++) {
        
        /* If the time is strickly the same as previous row, we skip this computation */
        if (mjd[row] != mjd1 ) {
            
            /* Full transformation from ICRS to Observed
             * update every t_skip, otherwise rotate earth only */
            if ( fabs (mjd[row]-mjd0) > t_skip ) {
                eraApco13 (2400000.5, mjd[row], dut1, lon, lat, elev,
                            pmx/3600.0*CPL_MATH_RAD_DEG, pmy/3600.0*CPL_MATH_RAD_DEG,
                            pressure, temperature, humidity, wavelength,
                            &astrom, &eo);
                mjd0 = mjd[row];
            }
            else
                eraAper13 (2400000.5, mjd[row] + dut1/(24.0*3600.0), &astrom);
            
            /* Apply proper motion (ICRS->BCRS) */
            eraPmpx(rc, dc, pmr, pmd, parallax, sysvel, astrom.pmt, astrom.eb, eWb);
            eraC2s(eWb, &rb, &db);
            
            /* Create (eU,eV,eW) in the BCRS */
            eraS2c(0.0, CPL_MATH_PI/2.0, eZb);  // eZc is the unit vector to the ICRS pole
            eraPxp(eZb, eWb, eUb);  // eUc is the cross product of eZc and eWc
            eraPn(eUb, &norm, eUb);  // eUc is normalized to a unit vector
            eraPxp(eWb, eUb, eVb);  // eWc is the cross product of eWc and eUc
            
            /* Create a 10 arcsec cardinal asterism around eWb in the BCRS */
            rotate_vector(eWb, -eps, eVb, eWb_up);
            rotate_vector(eWb, +eps, eVb, eWb_um);
            rotate_vector(eWb, +eps, eUb, eWb_vp);
            rotate_vector(eWb, -eps, eUb, eWb_vm);
            eraC2s(eWb_up, &rb_up, &db_up);
            eraC2s(eWb_um, &rb_um, &db_um);
            eraC2s(eWb_vp, &rb_vp, &db_vp);
            eraC2s(eWb_vm, &rb_vm, &db_vm);
            
            /* Transform the pointing direction and the cardinal asterism from BCRS to observed */
            eraAtboq (rb   , db   , &astrom, eWo   );
            eraAtboq (rb_up, db_up, &astrom, eWo_up);
            eraAtboq (rb_um, db_um, &astrom, eWo_um);
            eraAtboq (rb_vp, db_vp, &astrom, eWo_vp);
            eraAtboq (rb_vm, db_vm, &astrom, eWo_vm);
            
            /* Compute the observed (eUo,eVo,eWo) reference frame */
            eraPxp(eWo_up, eWo_um, eUo);
            eraPxp(eWo, eUo, eUo);
            eraSxp(1./(2.0*eps), eUo, eUo);
            eraPxp(eWo_vp, eWo_vm, eVo);
            eraPxp(eWo, eVo, eVo);
            eraSxp(1./(2.0*eps), eVo, eVo);
            
            /* Using eWo and zenith directions, compute eAz */
            eraPxp(eWo, ez, eAZo);
            eraPn(eAZo, &norm, eAZo);
            
            /* Using eWo and azimuth directions, compute eZd */
            eraPxp(eWo, eAZo, eZDo);
            eraPn(eZDo, &norm, eZDo);
        }
    
        /* If requested, store [E_U,E_V,E_V,E_AZ,E_ZD] columns */
        if (save_pointing) {
            if (mjd[row] != mjd1 ) {
                double * eU = cpl_malloc (sizeof(double) * 3);
                double * eV = cpl_malloc (sizeof(double) * 3);
                double * eW = cpl_malloc (sizeof(double) * 3);
                double * eAZ = cpl_malloc (sizeof(double) * 3);
                double * eZD = cpl_malloc (sizeof(double) * 3);
                for ( cpl_size c = 0; c < 3; c++) {
                    eU[c] = eUo[c];
                    eV[c] = eVo[c];
                    eW[c] = eWo[c];
                    eAZ[c] = eAZo[c];
                    eZD[c] = eZDo[c];
                }
                /* Wrap into vectors. It makes the data valid, and is the fastest
                * This and the duplication take most of the time of this function */
                p_u[row]  = cpl_array_wrap_double (eU,  3);
                p_v[row]  = cpl_array_wrap_double (eV,  3);
                p_w[row]  = cpl_array_wrap_double (eW,  3);
                p_az[row] = cpl_array_wrap_double (eAZ, 3);
                p_zd[row] = cpl_array_wrap_double (eZD, 3);
            } else {
                p_u[row]  = cpl_array_duplicate (p_u [row-1]);
                p_v[row]  = cpl_array_duplicate (p_v [row-1]);
                p_w[row]  = cpl_array_duplicate (p_w [row-1]);
                p_az[row] = cpl_array_duplicate (p_az[row-1]);
                p_zd[row] = cpl_array_duplicate (p_zd[row-1]);
            }
        }
        
        /* If requested, compute the projected baseline [UCOORD,VCOORD]
        * Note: The baseline length is corrected by the air refractive index
        *       at the HeNe wavelength and Paranal pressure. We want to use
        *       the vacuum baseline. */
        double n_air = 1.0002028;
        if (compute_uv) {
            int base = row % GRAVI_NBASE;
            uCoord[row] = (eUo[0] * baseline[base][0] + eUo[1] * baseline[base][1] + eUo[2] * baseline[base][2])/n_air;
            vCoord[row] = (eVo[0] * baseline[base][0] + eVo[1] * baseline[base][1] + eVo[2] * baseline[base][2])/n_air;
        }
        
        mjd1 = mjd[row];
        CPLCHECK_MSG ("Cannot run the ERFA transform");
    } /* End loop on rows */
    
    gravi_msg_function_exit(1);
    return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Compute the pointing directions and projected baselines in OI_VIS
 * 
 * @param p2vmred_data:   input/output data
 * @param eop_data:       data containing the EOP pameters
 * 
 * Compute the projected baselines [UCOORD,VCOORD] for FT and SC
 * Save the pointing directions [E_U,E_V,E_V,E_AZ,E_ZD] for the SC only
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_compute_pointing_uv (gravi_data * p2vmred_data, gravi_data * eop_data)
{
  gravi_msg_function_start(1);
  cpl_ensure_code (p2vmred_data, CPL_ERROR_NULL_INPUT);

  /* Get the header */
  cpl_propertylist * hdr_data = gravi_data_get_header (p2vmred_data);

  /* Get the OI_ARRAY table */
  cpl_table * oi_array = gravi_data_get_table (p2vmred_data, GRAVI_OI_ARRAY_EXT);

  /* For each type of data SC / FT */
  int type_data, ntype_data = 2;
  for (type_data = 0; type_data < ntype_data ; type_data ++) {

	int npol = gravi_pfits_get_pola_num (hdr_data, type_data);
	cpl_table * oi_vis =  gravi_data_get_oi_vis (p2vmred_data, type_data, 0, npol);

	/* Verbose */
	cpl_msg_info(cpl_func, "Compute pointing for %s ",type_data==GRAVI_FT?"FT":"SC");

    /* Compute for this polarisation */
    gravi_eop_pointing_uv (oi_vis, hdr_data,
                          (eop_data ? gravi_data_get_table_x (eop_data, 0) : NULL),
                          (eop_data ? gravi_data_get_header (eop_data) : NULL),
					      type_data==GRAVI_FT?0:1,
					      oi_array);
	CPLCHECK_MSG ("Cannot compute pointing");
    
	/* If second polarisation, for SC only, duplicate [E_U,E_V,E_V,E_AZ,E_ZD] */
	if ((npol > 1) && (type_data==GRAVI_SC)) {
	  cpl_msg_debug (cpl_func,"Duplicate in the 2nd polarisation");

	  cpl_table * oi_vis_1 =  gravi_data_get_oi_vis (p2vmred_data, type_data, 1, npol);
	  cpl_table_duplicate_column (oi_vis_1, "E_U", oi_vis, "E_U");
	  cpl_table_duplicate_column (oi_vis_1, "E_V", oi_vis, "E_V");
	  cpl_table_duplicate_column (oi_vis_1, "E_W", oi_vis, "E_W");
	  cpl_table_duplicate_column (oi_vis_1, "E_AZ", oi_vis, "E_AZ");
	  cpl_table_duplicate_column (oi_vis_1, "E_ZD", oi_vis, "E_ZD");

	  CPLCHECK_MSG ("Cannot duplicate");
	}

	/* If second polarisation, duplicate [UCOORD,VCOORD] */
	if (npol > 1) {
	  cpl_msg_info (cpl_func,"Duplicate in the 2nd polarisation");

	  cpl_table * oi_vis_1 =  gravi_data_get_oi_vis (p2vmred_data, type_data, 1, npol);
	  cpl_table_copy_data_double (oi_vis_1, "UCOORD", cpl_table_get_data_double (oi_vis, "UCOORD"));
	  cpl_table_copy_data_double (oi_vis_1, "VCOORD", cpl_table_get_data_double (oi_vis, "VCOORD"));
	}
	CPLCHECK_MSG ("Cannot duplicate");

  } /* End loop on FT/SC */

  gravi_msg_function_exit(1);
  return CPL_ERROR_NONE;
}

/**@}*/
