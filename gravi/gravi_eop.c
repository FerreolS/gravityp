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

#define _POSIX_C_SOURCE 200809L

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#include <unistd.h>
#include <fcntl.h>

#include <erfa.h>

#include "gravi_eop.h"
#include "gravi_cpl.h"
#include "gravi_pfits.h"

#include "gravi_utils.h"

/*-----------------------------------------------------------------------------
                             Private prototypes and define
 -----------------------------------------------------------------------------*/

#define LINE_SIZE 188
#define BUFFER_LENGTH 256

char * gravity_eop_get_ftp_file (int socketfd, int * data_length);
int gravity_get_socket_connection (const char * host, const char * port);
int gravity_eop_send_ftpcmd (int sockfd, const char *cmd);
int gravity_eop_send_pasv (int sockfd, const char *cmd);
int gravity_eop_ftp_reply (int sockfd, char ** message);
int gravity_eop_verify_ftp_code (char * msg, int length);

cpl_error_code gravi_eop_interpolate (cpl_size n, double *mjd,
									  double *pmx, double *pmy,
									  double *dut,
                                      cpl_table * eop_table,
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
    if (array_table > 0) {
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
    double baseline[6][3];
    double * uCoord;
    double * vCoord;
    if (compute_uv) {
        for ( int base = 0; base < 6; base ++) {
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
            int base = row % 6;
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

/*----------------------------------------------------------------------------*/
/**
 * @brief    Retrieve the Earth Orientation Parameters computed by IERS
 * 
 * @param    eop_host     The FTP host to retrieve the data from
 * @param    eop_urlpath  The full path to the data file
 * @param    data_length  The total size of the data buffer returned (returned)
 * @return   A string buffer with the full contents of the data
 * 
 * This function will connect to a given FTP host specified in eop_host
 * and the given eop_urlpath and retrieve the ascii file with the EOP data.
 * 
 * Possible #_cpl_error_code_ set in this function:
 * - CPL_ERROR_NULL_INPUT if eop_host, data_length or eop_urlpath are NULL
 * - CPL_ERROR_DATA_NOT_FOUND if the connection to the host cannot be 
 *   successfully established.
 * - CPL_ERROR_DATA_NOT_FOUND if the FTP transaction cannot be fullfilled.
 */
/*----------------------------------------------------------------------------*/

char * gravity_eop_download_finals2000A (const char * eop_host,
                                         const char * eop_urlpath,
                                         int * data_length)
{
    gravi_msg_function_start(1);
	
    const char ftp_port[] = "21";
    int cmd_socket, data_socket;

    /* Check and dump the input */
    cpl_ensure (eop_host,    CPL_ERROR_NULL_INPUT, NULL);
    cpl_ensure (eop_urlpath, CPL_ERROR_NULL_INPUT, NULL);
    cpl_ensure (data_length, CPL_ERROR_NULL_INPUT, NULL);
    cpl_msg_debug (cpl_func, "Using URL ftp://%s%s", eop_host, eop_urlpath);

    /* Getting the communication socket. use non-blocking mode for it */
    cmd_socket = gravity_get_socket_connection(eop_host, ftp_port);
    if (cmd_socket == 0) 
    {
        cpl_error_set_message(cpl_func, CPL_ERROR_DATA_NOT_FOUND,
            "Couldn't connect to the host");
        return NULL;
    }
    int flags = fcntl(cmd_socket, F_GETFL, 0);
    fcntl(cmd_socket, F_SETFL, flags | O_NONBLOCK);

    if(!gravity_eop_ftp_reply(cmd_socket, NULL))
    {
        cpl_error_set_message(cpl_func, CPL_ERROR_DATA_NOT_FOUND,
            "FTP server didn't reply");
        close(cmd_socket);
        return NULL;
    }
    cpl_msg_debug(cpl_func, "SEND");
    if(!gravity_eop_send_ftpcmd(cmd_socket, "USER anonymous\n"))
    {
        close(cmd_socket);
        return NULL;
    }
    if(!gravity_eop_send_ftpcmd(cmd_socket, "PASS ftp@eso.org\n"))
    {
        close(cmd_socket);
        return NULL;
    }
    int data_port = gravity_eop_send_pasv(cmd_socket, "PASV\n");
    if(!data_port)
    {
        close(cmd_socket);
        return NULL;
    }

    /* Getting the data socket in passive mode */
    char data_port_s[256];
    snprintf(data_port_s, 255, "%d", data_port);
    data_socket = gravity_get_socket_connection(eop_host, data_port_s);
    
    /* Retrieving the file */
    if(!gravity_eop_send_ftpcmd(cmd_socket, "TYPE I\n"))
        return NULL;
    char * retr_command = cpl_malloc(strlen(eop_urlpath) + 7);
    snprintf(retr_command, strlen(eop_urlpath) + 7, "RETR %s\n", eop_urlpath);
    if(!gravity_eop_send_ftpcmd(cmd_socket, retr_command))
    {
        close(cmd_socket);
        cpl_free(retr_command);
        return NULL;
    }
    cpl_free(retr_command);

    char * data_eop = gravity_eop_get_ftp_file(data_socket, data_length);

    /* Close connection and free resources */
    close(cmd_socket);
    close(data_socket);

    gravi_msg_function_exit(1);
    return data_eop;
}

int gravity_get_socket_connection (const char * host, const char * port)
{
    int sockfd;
    gravi_msg_function_start(0);

    /* IP Name resolution */

    /* Set the hints. First we initialize to 0 the structure.
       Only retrieve IPv4 or IPv6 if configured in the system */
    struct addrinfo hints;
    memset(&hints, 0, sizeof(hints));
    hints.ai_flags = AI_ADDRCONFIG;
    hints.ai_socktype = SOCK_STREAM;
    /* Getting the list of IP addresses */
    cpl_msg_debug(cpl_func, "Getting IP");
    struct addrinfo * addr_list ;
    if (getaddrinfo(host, port, &hints, &addr_list) != 0) {
        cpl_error_set_message(cpl_func, CPL_ERROR_DATA_NOT_FOUND,
            "Couldn't get address for host");
        return 0;
    }

    /* Connecting to the server for the FTP commands. 
       The first address to which we can connect will be used.
       addr_list is a linked list */
    cpl_msg_debug(cpl_func, "Connecting to server");
    struct addrinfo *this_addr;
    for(this_addr = addr_list; this_addr != NULL; this_addr = this_addr->ai_next)
    {
        /* Opening the socket */
        if ((sockfd = socket(this_addr->ai_family, this_addr->ai_socktype,
                this_addr->ai_protocol)) == -1) {
            continue;
        }

        if (connect(sockfd, this_addr->ai_addr, this_addr->ai_addrlen) == -1) {
            close(sockfd);
            continue;
        }
        cpl_msg_debug(cpl_func, "Connection established");
        break;
    }

    if (this_addr == NULL)
    {
        cpl_error_set_message(cpl_func, CPL_ERROR_DATA_NOT_FOUND,
            "Couldn't connect to the host");
        return 0;
    }

    freeaddrinfo(addr_list);

    gravi_msg_function_exit(0);
    return sockfd;
}

int gravity_eop_send_ftpcmd (int sockfd, const char *cmd)
{
    gravi_msg_function_start(0);
    cpl_msg_debug(cpl_func, "Sending FTP command <<%s>>", cmd);

    if(write(sockfd, cmd, strlen(cmd))==0)
    {
        cpl_error_set_message(cpl_func, CPL_ERROR_DATA_NOT_FOUND,
            "Problem during FTP transaction");
        return 0;
    }

    char * msg;
    if(!gravity_eop_ftp_reply(sockfd, &msg))
    {
        cpl_error_set_message(cpl_func, CPL_ERROR_DATA_NOT_FOUND,
            "Problem during FTP transaction");
        return 0;
    }
    if (msg != NULL)
        free(msg);
	
    gravi_msg_function_exit(0);
    return 1;
}


int gravity_eop_send_pasv (int sockfd, const char *cmd)
{
    gravi_msg_function_start(0);

    if(write(sockfd, cmd, strlen(cmd))==0)
    {
        cpl_error_set_message(cpl_func, CPL_ERROR_DATA_NOT_FOUND,
            "Problem during FTP transaction");
        return 0;
    }

    char * msg;
    char * new_con;
    unsigned int v[6];

    if(!gravity_eop_ftp_reply(sockfd, &msg))
    {
        cpl_error_set_message(cpl_func, CPL_ERROR_DATA_NOT_FOUND,
            "Problem during FTP transaction");
        return 0;
    }
    new_con = strchr(msg,'(');
    if (new_con == NULL)
        return 0;
    sscanf(new_con+1,"%u,%u,%u,%u,%u,%u",&v[2],&v[3],&v[4],&v[5],&v[0],&v[1]);

    //Get the new port connection of passive mode
    //This is coded in the reply of PASV
    //See http://www.freefire.org/articles/ftpexample.php
    int data_port = v[0] * 256 + v[1]; 

    if (msg != NULL)
        free(msg);
	
    gravi_msg_function_exit(0);
    return data_port;
}

int gravity_eop_ftp_reply (int sockfd, char ** message)
{
    gravi_msg_function_start(0);
    int  length, n;
    char buffer[BUFFER_LENGTH];
    char * msg = NULL;
   
    length = 0;
    if(message != NULL)
        *message = NULL;

    /* Wait for the server to fill the reply, since we are using non-blocking recv */
    cpl_msg_debug (cpl_func, "sleep for 5s...");
    sleep(5);
    cpl_msg_debug (cpl_func, "... done.");
	
    while( ( n = recv(sockfd, buffer, BUFFER_LENGTH - 1, 0) ) > 0)
    {
        if(msg == NULL)
            msg = strndup(buffer, n);
        else 
        {
            msg = realloc(msg, length + n + 1);
            strncpy(msg + length, buffer, n);
        }
        length += n;
    }
    if (errno == EAGAIN || errno == EWOULDBLOCK) // No messages were available
        errno= 0;
    if(length == 0)
    {
        free(msg);
        return 0;
    }
    msg[length] = '\0';

    cpl_msg_debug(cpl_func,"FTP reply: <<%s>>", msg);

    /* verify */
    int verify = gravity_eop_verify_ftp_code(msg, length + 1);

    if(message != NULL && verify)
        *message = msg;
    else
        free(msg);
	
    gravi_msg_function_exit(0);
    return verify;
}

char * gravity_eop_get_ftp_file(int sockfd, int * data_length)
{
    gravi_msg_function_start(1);
    int  length, n;
    char buffer[BUFFER_LENGTH];
    char * msg = NULL;

    length = 0;

    /* Sleep */
    cpl_msg_debug (cpl_func, "sleep 10s...");
    sleep(10);
    cpl_msg_debug (cpl_func, "... done.");
	
    /* Get the data */
    cpl_msg_info (cpl_func, "Get the data");
    while( ( n = recv(sockfd, buffer, BUFFER_LENGTH - 1, 0) ) > 0)
    {
        if(msg == NULL)
            msg = strndup(buffer, n);
        else 
        {
            msg = realloc(msg, length + n + 1);
            strncpy(msg + length, buffer, n);
        }
        length += n;
    }
    if(length == 0)
        return 0; 
    msg[length] = '\0';
    *data_length = length+1;

    gravi_msg_function_exit(1);
    return msg;
}

int gravity_eop_verify_ftp_code (char * msg, int length)
{
    gravi_msg_function_start(0);
    char * line = msg;

    //FTP protocol specifies that lines starting with 2xx codes are ok
    //Starting with 3xx are ok but the server expects some extra input
    //Starting with 4xx, 5xx or 6xx it denotes an error.
    while(line[0] == '1' || line[0] == '2' || line[0] == '3')
    {
        line = strchr(line, '\n');
        if(line == NULL || line - msg + 2 == length)
            return 1;
        line = line + 1;
    }
	
    gravi_msg_function_exit(0);
    return 0;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief    Export a raw string buffer containing EOP data to a CPL table
 * 
 * @param    eop_data     The string buffer with the data
 * @param    data_length  The total size of the data buffer 
 * @return   A string buffer with the full contents of the data
 * 
 * This function convert the ascii file retrieve from FTP
 * and convert it to a CPL table.
 * 
 * Possible #_cpl_error_code_ set in this function:
 * - CPL_ERROR_NULL_INPUT if eop_data is NULL
 * - CPL_ERROR_NULL_INPUT if data_length doesn't correspond with the expected 
 *   EOP records length
 */
/*----------------------------------------------------------------------------*/

cpl_table * gravity_eop_data_totable (const char * eop_data, int data_length)
{
    gravi_msg_function_start(1);
    cpl_ensure (eop_data, CPL_ERROR_NULL_INPUT, NULL);

    if(!((data_length - 1) % LINE_SIZE == 0))
    {
        cpl_error_set_message(cpl_func, CPL_ERROR_NULL_INPUT,
            "Raw data doesn't have a fixed record width");
        return 0;
    }

    /* Create tables */
    cpl_size n_entries = (data_length - 1 ) / LINE_SIZE;
    cpl_table * eop_table = cpl_table_new (n_entries);
    cpl_msg_info (cpl_func, " EOP data has a total of %"CPL_SIZE_FORMAT" entries", n_entries);

    /* Create columns */
    cpl_table_new_column (eop_table, "MJD",  CPL_TYPE_DOUBLE);
    cpl_table_new_column (eop_table, "PMX",  CPL_TYPE_DOUBLE);
    cpl_table_new_column (eop_table, "PMY",  CPL_TYPE_DOUBLE);
    cpl_table_new_column (eop_table, "DUT",  CPL_TYPE_DOUBLE);
    cpl_table_new_column (eop_table, "FLAG", CPL_TYPE_STRING);

    /* Set units */
    cpl_table_set_column_unit (eop_table, "MJD", "d");
    cpl_table_set_column_unit (eop_table, "PMX", "arcsec");
    cpl_table_set_column_unit (eop_table, "PMY", "arcsec");
    cpl_table_set_column_unit (eop_table, "DUT", "s");

    /* Fill the columns from the string buffer */
    for(cpl_size i=0; i<n_entries; i++)
    {
        char flag[3];
        strncpy(flag, eop_data+i*LINE_SIZE+16, 1);
        flag[2] = '\0';
        cpl_table_set_string(eop_table, "FLAG", i, flag);

        cpl_table_set_double(eop_table, "MJD", i, atof(eop_data+i*LINE_SIZE+7));
        if(!strncmp(flag, "I", 1) || !strncmp(flag, "P", 1))
        {
            cpl_table_set_double(eop_table, "PMX", i, atof(eop_data+i*LINE_SIZE+18));
            cpl_table_set_double(eop_table, "PMY", i, atof(eop_data+i*LINE_SIZE+37));
            cpl_table_set_double(eop_table, "DUT", i, atof(eop_data+i*LINE_SIZE+58));
        }
    }

    /* Remove the NULL columns */
    cpl_table_unselect_all (eop_table);
    cpl_table_or_selected_invalid (eop_table, "PMX");
    cpl_table_or_selected_invalid (eop_table, "PMY");
    cpl_table_or_selected_invalid (eop_table, "DUT");
    cpl_msg_info (cpl_func,"Found %lld invalid", cpl_table_count_selected (eop_table));
    cpl_table_erase_selected (eop_table);

    gravi_msg_function_exit(1);
    return eop_table;
}

/**@}*/
