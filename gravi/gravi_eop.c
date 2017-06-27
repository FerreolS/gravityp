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
									  double *dut, gravi_data * eop_data);

void eraAtcoq (double rc, double dc, double pmr, double pmd, double px,
			   double rv, eraASTROM *astrom, double enuob[3]);
void dtp2s (double xi, double eta, double raz,
			double decz, double *ra, double *dec);

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

cpl_error_code gravi_eop_interpolate (cpl_size n, double *mjd, double *pmx, double *pmy, double *dut, gravi_data * eop_data)
{  
  gravi_msg_function_start(1);
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
  cpl_vector_fill (vpmx, 0.0);
  cpl_vector_fill (vpmy, 0.0);
  cpl_vector_fill (vdut, 0.0);

  /* Check if data available. If not, we return zero */
  if ( eop_data == NULL )
  {
	cpl_msg_warning (cpl_func, "No EOP_PARAM. Use EOP and DUT=0.0s");
    cpl_vector_unwrap (vmjd);
    cpl_vector_unwrap (vpmx);
    cpl_vector_unwrap (vpmy);
    cpl_vector_unwrap (vdut);
    return CPL_ERROR_NONE;
  }
  
  cpl_table * eop_table = gravi_data_get_table_x (eop_data, 0);
  cpl_propertylist * header = gravi_data_get_header (eop_data);

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
 * Quick transform from Celestial to Observed
 *     rc,dc  double   ICRS right ascension at J2000.0 (radians)
 *     pr     double   RA proper motion dRA/dt, not cos(Dec) x dRA/dt (radians/year)
 *     pd     double   Dec proper motion dDec/dt (radians/year)
 *     px     double   parallax (arcsec)
 *     rv     double   radial velocity (km/s, +ve if receding)
 */
void eraAtcoq(double rc, double dc, double pmr, double pmd, double px, double rv, eraASTROM *astrom, double enuob[3])
{
    double ri, di;
    double aob, zob, hob, dob, rob;

    /* Transform from Celestial to Intermediate */
    eraAtciq (rc, dc, pmr, pmd, px, rv, astrom, &ri, &di);

    /* Transform Intermediate to Observed. */
    eraAtioq (ri, di, astrom, &aob, &zob, &hob, &dob, &rob);

    /* Convert from equatorial to cartesian */
    enuob[0] = sin(zob) * sin(aob);
    enuob[1] = sin(zob) * cos(aob);
    enuob[2] = cos(zob);
}

void dtp2s(double xi, double eta, double raz, double decz, double *ra, double *dec)
{
    double sdecz = sin(decz);
    double cdecz = cos(decz);
    double denom = cdecz - eta * sdecz;
    double d = atan2(xi, denom) + raz;
    *ra = fmod(d, 2.0 * CPL_MATH_PI);
    *dec = atan2(sdecz + eta * cdecz, sqrt(xi * xi + denom * denom));
}

void difference(double x[3], double y[3], double z[3])
{
    for (int i=0; i<3; i++) {
	  z[i] = x[i] - y[i];
    }
}

void multiply(double xyz[3], double factor)
{
    xyz[0] *= factor;
    xyz[1] *= factor;
    xyz[2] *= factor;
}

void normalize(double xyz[3])
{
    double norm = sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1] + xyz[2]*xyz[2]);
    multiply (xyz, 1.0/norm);
}

void cross(double x[3], double y[3], double z[3])
{
    z[0] = x[1] * y[2] - x[2] * y[1];
    z[1] = x[2] * y[0] - x[0] * y[2];
    z[2] = x[0] * y[1] - x[1] * y[0];
}



cpl_error_code gravi_eop_pointing (cpl_table * input_table,
                                   cpl_propertylist * header,
                                   gravi_data * eop_data)
{
    gravi_msg_function_start(1);
    cpl_ensure_code (input_table, CPL_ERROR_NULL_INPUT);
    
    double t_skip = 1./24/3600, mjd0 = -1.0, mjd1 = -1.0;
    cpl_msg_info (cpl_func, "Compute pointing with full ERFA every %.2f s",t_skip*24*3600.0);
    
    /* Read UTC as MJD */
    double *mjd = cpl_table_get_data_double (input_table, "MJD");
    double mean_mjd = cpl_table_get_column_mean (input_table, "MJD");
    
	/* Create new columns */
	cpl_table_new_column_array (input_table, "E_U", CPL_TYPE_DOUBLE, 3);
	cpl_table_new_column_array (input_table, "E_V", CPL_TYPE_DOUBLE, 3);
	cpl_table_new_column_array (input_table, "E_W", CPL_TYPE_DOUBLE, 3);
	cpl_table_new_column_array (input_table, "E_AZ", CPL_TYPE_DOUBLE, 3);
	cpl_table_new_column_array (input_table, "E_ZD", CPL_TYPE_DOUBLE, 3);
	cpl_array ** p_u = cpl_table_get_data_array (input_table,"E_U");
	cpl_array ** p_v = cpl_table_get_data_array (input_table,"E_V");
	cpl_array ** p_w = cpl_table_get_data_array (input_table,"E_W");
	cpl_array ** p_az = cpl_table_get_data_array (input_table,"E_AZ");
	cpl_array ** p_zd = cpl_table_get_data_array (input_table,"E_ZD");
    CPLCHECK_MSG ("Cannot create data");

	/* Checked by J. Woillez */
	double elev = gravi_pfits_get_geoelev (header); // Height in [m]
	double lon = gravi_pfits_get_geolon (header) * CPL_MATH_RAD_DEG; // Lat in [rad]
	double lat = gravi_pfits_get_geolat (header) * CPL_MATH_RAD_DEG; // Lon in [rad], East positive
	CPLCHECK_MSG ("Cannot get the data");

	/* Compute Earth Orientation Paramters at given MJDs */
	double dut1, pmx, pmy;
	gravi_eop_interpolate (1, &mean_mjd, &pmx, &pmy, &dut1, eop_data);
	CPLCHECK_MSG ("Cannot interpolate");

	/* We use the FT coordinate for all uv coordinates */
	double raep = gravi_pfits_get_mid_raep (header); // [rad]
	double decp = gravi_pfits_get_mid_decep (header); // [rad]
	double pmra   = gravi_pfits_get_pmra (header) * CPL_MATH_RAD_DEG / 3600.0 / cos(decp); // dRA/dt, not cos(Dec)xdRA/dt [rad/year]
	double pmdec  = gravi_pfits_get_pmdec (header) * CPL_MATH_RAD_DEG / 3600.0; // dDec/dt [rad/year]
	double parallax = gravi_pfits_get_plx (header); // [as]
	double sysvel = 0.0;
	CPLCHECK_MSG ("Cannot get the header data");

	/* Allocate memory for the tmp computations */
	eraASTROM astrom;
	double eo, rcUp, dcUp, rcUm, dcUm, rcVp, dcVp, rcVm, dcVm;
	double enuobUp[3], enuobUm[3], enuobVp[3], enuobVm[3];
	double ez[3] = {0.0, 0.0, 1.0}; // Zenith direction in ENU frame

	/* Step for finite difference, 10 arcsec chosen for optimal accuracy */
	double eps = 10.0 / 3600.0 * CPL_MATH_RAD_DEG; // [rad]

	/* Prepare centered finite differences
	 * eU corresponds to +RA
	 * eV corresponds to +DEC */
	dtp2s (+eps, 0.0, raep, decp, &rcUp, &dcUp);
	dtp2s (-eps, 0.0, raep, decp, &rcUm, &dcUm);
	dtp2s (0.0, +eps, raep, decp, &rcVp, &dcVp);
	dtp2s (0.0, -eps, raep, decp, &rcVm, &dcVm);

	/* Loop on rows. the baselines are not assumed to share the
	 * same MJD since this function is called after the averaging */
	cpl_size n_row = cpl_table_get_nrow (input_table);
	for (cpl_size row = 0 ; row < n_row ; row++) {

	  /* Transformation from ICRS to Observed (neglect refraction, all at zero)
	   * Update precession/nudation every t_skip, otherwise rotate earth only
	   * If the time is strickly the same as previous row, we skip this computation */

	  if (mjd[row] != mjd1 ) {

		/* Allocate memory */
		double * eu = cpl_malloc (sizeof(double) * 3);
		double * ev = cpl_malloc (sizeof(double) * 3);
		double * ew = cpl_malloc (sizeof(double) * 3);
		double * eaz = cpl_malloc (sizeof(double) * 3);
		double * ezd = cpl_malloc (sizeof(double) * 3);

		if ( fabs (mjd[row]-mjd0) > t_skip ) {
		  eraApco13 (2400000.5, mjd[row], dut1, lon, lat, elev, pmx/3600.0*CPL_MATH_RAD_DEG, pmy/3600.0*CPL_MATH_RAD_DEG, 0.0, 0.0, 0.0, 0.0, &astrom, &eo);
		  mjd0 = mjd[row];
		}
		else
		  eraAper13 (2400000.5, mjd[row] + dut1/(24.0*3600.0), &astrom);

		/* Transform from celestial to intermediate, compute eU */
		eraAtcoq (rcUp, dcUp, pmra, pmdec, parallax, sysvel, &astrom, enuobUp);
		eraAtcoq (rcUm, dcUm, pmra, pmdec, parallax, sysvel, &astrom, enuobUm);
		difference (enuobUp, enuobUm, eu);
		normalize (eu);

		/* Transform from celestial to intermediate, compute eV */
		eraAtcoq (rcVp, dcVp, pmra, pmdec, parallax, sysvel, &astrom, enuobVp);
		eraAtcoq (rcVm, dcVm, pmra, pmdec, parallax, sysvel, &astrom, enuobVm);
		difference (enuobVp, enuobVm, ev);
		normalize (ev);

		/* Transform from celestial to intermediate, compute eW */
		eraAtcoq (raep, decp, pmra, pmdec, parallax, sysvel, &astrom, ew);
		multiply (ew, -1.0);

		/* Using pointing and zenith directions, compute eAz */
		cross (ew, ez, eaz);
		normalize (eaz);

		/* Using pointing and azimuth directions, compute eZd */
		cross (ew, eaz, ezd);
		normalize (ezd);

		/* Wrap into vectors. It makes the data valid, and is the fastest
		 * This and the duplication take most of the time of this function */
		p_u[row] = cpl_array_wrap_double (eu, 3);
		p_v[row] = cpl_array_wrap_double (ev, 3);
		p_w[row] = cpl_array_wrap_double (ew, 3);
		p_az[row] = cpl_array_wrap_double (eaz, 3);
		p_zd[row] = cpl_array_wrap_double (ezd, 3);

	  } else {
		p_u[row] = cpl_array_duplicate (p_u[row-1]);
		p_v[row] = cpl_array_duplicate (p_v[row-1]);
		p_w[row] = cpl_array_duplicate (p_w[row-1]);
		p_az[row] = cpl_array_duplicate (p_az[row-1]);
		p_zd[row] = cpl_array_duplicate (p_zd[row-1]);
	  }

	  mjd1 = mjd[row];
	  CPLCHECK_MSG ("Cannot run the ERFA transform");
	} /* End loop on rows */
  
  gravi_msg_function_exit(1);
  return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Compute the pointing direction
 * 
 * @param p2vmred_data:   input/output data
 * @param eop_data:       data containing the EOP pameters
 * 
 * Compute the pointing direction in Observed coordinate for every frames of
 * the SC, so that the real-time projected baseline can be recomputed easily
 * off-line. [e_u, e_v, e_w] are the unitary vectors of the uv-plane
 * on the local coordinates. These quantities are stored as new column in the
 * OI_VIS tables of the SC.
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_compute_pointing (gravi_data * p2vmred_data, gravi_data * eop_data)
{
  gravi_msg_function_start(1);
  cpl_ensure_code (p2vmred_data, CPL_ERROR_NULL_INPUT);

  /* Get the header */
  cpl_propertylist * hdr_data = gravi_data_get_header (p2vmred_data);

  /* For each type of data SC / FT */
  int type_data, ntype_data = 2;
  for (type_data = 0; type_data < ntype_data ; type_data ++) {

	int pol = 0, npol = gravi_pfits_get_pola_num (hdr_data, type_data);
	cpl_table * oi_vis =  gravi_data_get_oi_vis (p2vmred_data, type_data, pol, npol);
	double t_skip = 1./24/3600, mjd0 = -1.0, mjd1 = -1.0;

	/* Verbose */
	if (type_data == GRAVI_FT) {
	  cpl_msg_debug (cpl_func, "Don't compute ERFA pointing for FT"); continue;
	}
	cpl_msg_info(cpl_func, "Compute pointing for %s (full ERFA every %.2f s)",type_data==GRAVI_FT?"FT":"SC",t_skip*24*3600.0);

	/* Read UTC as MJD */
	double *mjd = cpl_table_get_data_double (oi_vis, "MJD");
	double mean_mjd = cpl_table_get_column_mean (oi_vis, "MJD");

	/* Create new columns */
	cpl_table_new_column_array (oi_vis, "E_U", CPL_TYPE_DOUBLE, 3);
	cpl_table_new_column_array (oi_vis, "E_V", CPL_TYPE_DOUBLE, 3);
	cpl_table_new_column_array (oi_vis, "E_W", CPL_TYPE_DOUBLE, 3);
	cpl_table_new_column_array (oi_vis, "E_AZ", CPL_TYPE_DOUBLE, 3);
	cpl_table_new_column_array (oi_vis, "E_ZD", CPL_TYPE_DOUBLE, 3);
	cpl_array ** p_u = cpl_table_get_data_array (oi_vis,"E_U");
	cpl_array ** p_v = cpl_table_get_data_array (oi_vis,"E_V");
	cpl_array ** p_w = cpl_table_get_data_array (oi_vis,"E_W");
	cpl_array ** p_az = cpl_table_get_data_array (oi_vis,"E_AZ");
	cpl_array ** p_zd = cpl_table_get_data_array (oi_vis,"E_ZD");

	/* Checked by J. Woillez */
	double elev = gravi_pfits_get_geoelev (hdr_data); // Height in [m]
	double lon = gravi_pfits_get_geolon (hdr_data) * CPL_MATH_RAD_DEG; // Lat in [rad]
	double lat = gravi_pfits_get_geolat (hdr_data) * CPL_MATH_RAD_DEG; // Lon in [rad], East positive
	CPLCHECK_MSG ("Cannot get the data");

	/* Compute Earth Orientation Paramters at given MJDs */
	double dut1, pmx, pmy;
	gravi_eop_interpolate (1, &mean_mjd, &pmx, &pmy, &dut1, eop_data);
	CPLCHECK_MSG ("Cannot interpolate");

	/* We use the FT coordinate for all uv coordinates */
	double raep = gravi_pfits_get_mid_raep (hdr_data); // [rad]
	double decp = gravi_pfits_get_mid_decep (hdr_data); // [rad]
	double pmra   = gravi_pfits_get_pmra (hdr_data) * CPL_MATH_RAD_DEG / 3600.0 / cos(decp); // dRA/dt, not cos(Dec)xdRA/dt [rad/year]
	double pmdec  = gravi_pfits_get_pmdec (hdr_data) * CPL_MATH_RAD_DEG / 3600.0; // dDec/dt [rad/year]
	double parallax = gravi_pfits_get_plx (hdr_data); // [as]
	double sysvel = 0.0;
	CPLCHECK_MSG ("Cannot get the header data");

	/* Allocate memory for the tmp computations */
	eraASTROM astrom;
	double eo, rcUp, dcUp, rcUm, dcUm, rcVp, dcVp, rcVm, dcVm;
	double enuobUp[3], enuobUm[3], enuobVp[3], enuobVm[3];
	double ez[3] = {0.0, 0.0, 1.0}; // Zenith direction in ENU frame

	/* Step for finite difference, 10 arcsec chosen for optimal accuracy */
	double eps = 10.0 / 3600.0 * CPL_MATH_RAD_DEG; // [rad]

	/* Prepare centered finite differences
	 * eU corresponds to +RA
	 * eV corresponds to +DEC */
	dtp2s (+eps, 0.0, raep, decp, &rcUp, &dcUp);
	dtp2s (-eps, 0.0, raep, decp, &rcUm, &dcUm);
	dtp2s (0.0, +eps, raep, decp, &rcVp, &dcVp);
	dtp2s (0.0, -eps, raep, decp, &rcVm, &dcVm);

	/* Loop on rows. the baselines are not assumed to share the
	 * same MJD since this function is called after the averaging */
	cpl_size n_row = cpl_table_get_nrow (oi_vis);
	for (cpl_size row = 0 ; row < n_row ; row++) {

	  /* Transformation from ICRS to Observed (neglect refraction, all at zero)
	   * Update precession/nudation every t_skip, otherwise rotate earth only
	   * If the time is strickly the same as previous row, we skip this computation */

	  if (mjd[row] != mjd1 ) {

		/* Allocate memory */
		double * eu = cpl_malloc (sizeof(double) * 3);
		double * ev = cpl_malloc (sizeof(double) * 3);
		double * ew = cpl_malloc (sizeof(double) * 3);
		double * eaz = cpl_malloc (sizeof(double) * 3);
		double * ezd = cpl_malloc (sizeof(double) * 3);

		if ( fabs (mjd[row]-mjd0) > t_skip ) {
		  eraApco13 (2400000.5, mjd[row], dut1, lon, lat, elev, pmx/3600.0*CPL_MATH_RAD_DEG, pmy/3600.0*CPL_MATH_RAD_DEG, 0.0, 0.0, 0.0, 0.0, &astrom, &eo);
		  mjd0 = mjd[row];
		}
		else
		  eraAper13 (2400000.5, mjd[row] + dut1/(24.0*3600.0), &astrom);

		/* Transform from celestial to intermediate, compute eU */
		eraAtcoq (rcUp, dcUp, pmra, pmdec, parallax, sysvel, &astrom, enuobUp);
		eraAtcoq (rcUm, dcUm, pmra, pmdec, parallax, sysvel, &astrom, enuobUm);
		difference (enuobUp, enuobUm, eu);
		normalize (eu);

		/* Transform from celestial to intermediate, compute eV */
		eraAtcoq (rcVp, dcVp, pmra, pmdec, parallax, sysvel, &astrom, enuobVp);
		eraAtcoq (rcVm, dcVm, pmra, pmdec, parallax, sysvel, &astrom, enuobVm);
		difference (enuobVp, enuobVm, ev);
		normalize (ev);

		/* Transform from celestial to intermediate, compute eW */
		eraAtcoq (raep, decp, pmra, pmdec, parallax, sysvel, &astrom, ew);
		multiply (ew, -1.0);

		/* Using pointing and zenith directions, compute eAz */
		cross (ew, ez, eaz);
		normalize (eaz);

		/* Using pointing and azimuth directions, compute eZd */
		cross (ew, eaz, ezd);
		normalize (ezd);

		/* Wrap into vectors. It makes the data valid, and is the fastest
		 * This and the duplication take most of the time of this function */
		p_u[row] = cpl_array_wrap_double (eu, 3);
		p_v[row] = cpl_array_wrap_double (ev, 3);
		p_w[row] = cpl_array_wrap_double (ew, 3);
		p_az[row] = cpl_array_wrap_double (eaz, 3);
		p_zd[row] = cpl_array_wrap_double (ezd, 3);

	  } else {
		p_u[row] = cpl_array_duplicate (p_u[row-1]);
		p_v[row] = cpl_array_duplicate (p_v[row-1]);
		p_w[row] = cpl_array_duplicate (p_w[row-1]);
		p_az[row] = cpl_array_duplicate (p_az[row-1]);
		p_zd[row] = cpl_array_duplicate (p_zd[row-1]);
	  }

	  mjd1 = mjd[row];
	  CPLCHECK_MSG ("Cannot run the ERFA transform");
	} /* End loop on rows */

	/* Fill second polarisation */
	if ( npol > 1) {
	  cpl_msg_debug (cpl_func,"Duplicate in the 2sd polarisation");

	  cpl_table * oi_vis_1 =  gravi_data_get_oi_vis (p2vmred_data, type_data, 1, npol);
	  cpl_table_duplicate_column (oi_vis_1, "E_U", oi_vis, "E_U");
	  cpl_table_duplicate_column (oi_vis_1, "E_V", oi_vis, "E_V");
	  cpl_table_duplicate_column (oi_vis_1, "E_W", oi_vis, "E_W");
	  cpl_table_duplicate_column (oi_vis_1, "E_AZ", oi_vis, "E_AZ");
	  cpl_table_duplicate_column (oi_vis_1, "E_ZD", oi_vis, "E_ZD");

	  CPLCHECK_MSG ("Cannot duplicate");
	}

  } /* End loop on FT/SC */

  gravi_msg_function_exit(1);
  return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Compute the UCOORD and VCOORD (uv)
 * 
 * @param p2vmred_data:   input/output data
 * @param eop_data:       data containing the EOP pameters
 * 
 * This function re-computes the UCOORD and VCOORD of the OI_VIS
 * tables using the HEADER coordinates, the OI_ARRAY baseline and
 * the MJD time of the OI_VIS table. It is based on ERFA library.
 * The SC and FT OI_VIS tables are updated.
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_compute_uv (gravi_data * p2vmred_data, gravi_data * eop_data)
{
  gravi_msg_function_start(1);
  cpl_ensure_code (p2vmred_data, CPL_ERROR_NULL_INPUT);

  /* Get the header */
  cpl_propertylist * hdr_data = gravi_data_get_header (p2vmred_data);
  cpl_table * oi_array = gravi_data_get_table (p2vmred_data,GRAVI_OI_ARRAY_EXT);

  /* For each type of data SC / FT */
  int type_data, ntype_data = 2;
  for (type_data = 0; type_data < ntype_data ; type_data ++) {

	int pol = 0, npol = gravi_pfits_get_pola_num (hdr_data, type_data);
	cpl_table * oi_vis =  gravi_data_get_oi_vis (p2vmred_data, type_data, pol, npol);
	double t_skip = 2./24/3600, mjd0 = -1.0, mjd1 = -1.0;

	cpl_msg_info(cpl_func, "Compute uv for %s (full ERFA every %.2f s)",type_data==GRAVI_FT?"FT":"SC",t_skip*24*3600.0);

	/* Read UTC as MJD */
	double *mjd = cpl_table_get_data_double (oi_vis, "MJD");
	double mean_mjd = cpl_table_get_column_mean (oi_vis, "MJD");

	/* Checked by J. Woillez */
	double elev = gravi_pfits_get_geoelev (hdr_data); // Height in [m]
	double lon = gravi_pfits_get_geolon (hdr_data) * CPL_MATH_RAD_DEG; // Lat in [rad]
	double lat = gravi_pfits_get_geolat (hdr_data) * CPL_MATH_RAD_DEG; // Lon in [rad], East positive

	/* Compute Earth Orientation Paramters at mean MJD */
	double dut1, pmx, pmy;
	gravi_eop_interpolate (1, &mean_mjd, &pmx, &pmy, &dut1, eop_data);
    
	/* We use the FT coordinate for all uv coordinates */
	double raep = gravi_pfits_get_mid_raep (hdr_data); // [rad]
	double decp = gravi_pfits_get_mid_decep (hdr_data); // [rad]
	double pmra   = gravi_pfits_get_pmra (hdr_data) * CPL_MATH_RAD_DEG / 3600.0 / cos(decp); // dRA/dt, not cos(Dec)xdRA/dt [rad/year]
	double pmdec  = gravi_pfits_get_pmdec (hdr_data) * CPL_MATH_RAD_DEG / 3600.0; // dDec/dt [rad/year]
	double parallax = gravi_pfits_get_plx (hdr_data); // [as]
	double sysvel = 0.0;

	/* Allocate memory for the tmp computations */
	eraASTROM astrom;
	double eo, rcUp, dcUp, rcUm, dcUm, rcVp, dcVp, rcVm, dcVm;
	double enuobUp[3], enuobUm[3], enuobVp[3], enuobVm[3], eu[3], ev[3];

	CPLCHECK_MSG ("Cannot get the data");

	/* Step for finite difference, 10 arcsec chosen for optimal accuracy */
	double eps = 10.0 / 3600.0 * CPL_MATH_RAD_DEG; // [rad]

	/* Prepare centered finite differences
	 * eU corresponds to +RA
	 * eV corresponds to +DEC */
	dtp2s (+eps, 0.0, raep, decp, &rcUp, &dcUp);
	dtp2s (-eps, 0.0, raep, decp, &rcUm, &dcUm);
	dtp2s (0.0, +eps, raep, decp, &rcVp, &dcVp);
	dtp2s (0.0, -eps, raep, decp, &rcVm, &dcVm);

	/* Loop on bases to compute the physical 3D baseline as T2-T1 [m]
	 * Note: The telescope locations STAXYZ are given in the [West,South,Up]
	 * frame (ESO convention), whereas the UVW, calculated later on, are in
	 * the [East,North,Up] frame. Hence the sign changes below on X and Y. */
	double baseline[6][3];
	for ( int base = 0; base < 6; base ++) {
	  int tel1=0; while ( cpl_table_get (oi_array, "STA_INDEX", tel1, NULL) != gravi_table_get_value (oi_vis, "STA_INDEX", base, 0) ) tel1++;
	  int tel2=0; while ( cpl_table_get (oi_array, "STA_INDEX", tel2, NULL) != gravi_table_get_value (oi_vis, "STA_INDEX", base, 1) ) tel2++;
	  baseline[base][0] = -(gravi_table_get_value (oi_array, "STAXYZ", tel2, 0) - gravi_table_get_value (oi_array, "STAXYZ", tel1, 0));
	  baseline[base][1] = -(gravi_table_get_value (oi_array, "STAXYZ", tel2, 1) - gravi_table_get_value (oi_array, "STAXYZ", tel1, 1));
	  baseline[base][2] = +(gravi_table_get_value (oi_array, "STAXYZ", tel2, 2) - gravi_table_get_value (oi_array, "STAXYZ", tel1, 2));
	}

	double * uCoord = cpl_table_get_data_double (oi_vis, "UCOORD");
	double * vCoord = cpl_table_get_data_double (oi_vis, "VCOORD");

	/* Loop on rows. */
	cpl_size n_row = cpl_table_get_nrow (oi_vis);
	for (cpl_size row = 0 ; row < n_row ; row++) {

	  /* Transformation from ICRS to Observed (neglect refraction, all at zero)
	   * Update precession/nudation every t_skip, otherwise rotate earth only
	   * If the time is strickly the same as previous row, we skip this computation */

	  if (mjd[row] != mjd1 ) {

		if ( fabs (mjd[row]-mjd0) > t_skip ) {
		  eraApco13 (2400000.5, mjd[row], dut1, lon, lat, elev, pmx/3600.0*CPL_MATH_RAD_DEG, pmy/3600.0*CPL_MATH_RAD_DEG, 0.0, 0.0, 0.0, 0.0, &astrom, &eo);
		  mjd0 = mjd[row]; }
		else
		  eraAper13 (2400000.5, mjd[row] + dut1/(24.0*3600.0), &astrom);

		/* Transform from celestial to intermediate, compute eU */
		eraAtcoq (rcUp, dcUp, pmra, pmdec, parallax, sysvel, &astrom, enuobUp);
		eraAtcoq (rcUm, dcUm, pmra, pmdec, parallax, sysvel, &astrom, enuobUm);
		difference (enuobUp, enuobUm, eu);
		normalize (eu);

		/* Transform from celestial to intermediate, compute eV */
		eraAtcoq (rcVp, dcVp, pmra, pmdec, parallax, sysvel, &astrom, enuobVp);
		eraAtcoq (rcVm, dcVm, pmra, pmdec, parallax, sysvel, &astrom, enuobVm);
		difference (enuobVp, enuobVm, ev);
		normalize (ev);
	  }

 	  /* Project physical baseline into u,v */
	  int base = row % 6;
	  uCoord[row] = eu[0] * baseline[base][0] + eu[1] * baseline[base][1] + eu[2] * baseline[base][2];
	  vCoord[row] = ev[0] * baseline[base][0] + ev[1] * baseline[base][1] + ev[2] * baseline[base][2];

	  mjd1 = mjd[row];
	  CPLCHECK_MSG ("Cannot compute the uv");
	} /* End loop on rows */

	/* Copy second polarisation. Assume they have same uv-plane */
	if (npol > 1) {
	  cpl_msg_info (cpl_func,"Duplicate in the 2sd polarisation");

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
