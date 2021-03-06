/* $Id: gravi_vis.c,v 1.10 2014/11/12 15:10:40 nazouaoui Exp $
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
 * @defgroup gravi_vis  Averaging the individual DITs into a final OIFITS
 *
 * This module processes the p2vmreduced file to produce the final product of the
 * @c gravity_vis recipe. The functions @c gravi_compute_vis()  computes the averaged
 * quantities, main steps are :
 * - @c gravi_vis_average_bootstrap() to average the visibilities (see
 * Algorithms/Average complex visibilities, Algorithms/Average squared visibilities)
 * - @c gravi_flux_average_bootstrap() to average the fluxes (see Algorithms/AverageFlux)
 * - @c gravi_t3_average_bootstrap() (see Algorithms/Average closure-phase)
 *
 */
/**@{*/

/*
 * History
 *    21.11.2018   memory leak in gravi_average_self_visphi
 *                      changes marked as 'EKW'
 *    11/01/2019   Fix Warning parameter 'ret'
 */
/*-----------------------------------------------------------------------------
                                   Includes
 -----------------------------------------------------------------------------*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cpl.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <complex.h>
#include <string.h>

#include "gravi_data.h"
#include "gravi_dfs.h"
#include "gravi_pfits.h"
#include "gravi_cpl.h"

#include "gravi_utils.h"

#include "gravi_vis.h"
#include "gravi_eop.h"
#include "gravi_tf.h"

/*-----------------------------------------------------------------------------
                              Private prototypes
 -----------------------------------------------------------------------------*/

double gravi_randn (void);

cpl_error_code gravi_array_online_variance(cpl_array * data, cpl_array * mean, cpl_array * variance, int n);
cpl_error_code gravi_array_online_variance_res(cpl_array ** data,
											   int n, int rephase);
cpl_error_code gravi_flux_average_bootstrap(cpl_table * oi_flux_avg,
											cpl_table * oi_flux,
											int nboot,
					    double outlier_threshold);
cpl_error_code gravi_t3_average_bootstrap(cpl_table * oi_t3_avg,
										  cpl_table * oi_vis,
										  cpl_table * oi_flux,
										  int nboot,
										  int use_vFactor,
                                          int use_pFactor,
					  double outlier_threshold);
cpl_error_code gravi_vis_average_bootstrap (cpl_table * oi_vis_avg,
											cpl_table * oi_vis2_avg,
											cpl_table * oi_vis,
											int nboot, 
											const char * phase_ref,
											int use_vFactor,
                                            int use_pFactor,
					    int use_debiasing,
					    double outlier_threshold);

cpl_error_code gravi_vis_flag_nan (cpl_table * oi_table);

cpl_error_code gravi_vis_average_amp (cpl_table *oi_table, const char *name,  const char *err, int nbase);
cpl_error_code gravi_vis_average_phi (cpl_table *oi_table, const char *name,  const char *err, int nbase);
cpl_error_code gravi_vis_average_value (cpl_table *oi_table, const char *name,  const char *err, int nbase);
cpl_error_code gravi_vis_resamp_amp (cpl_table * oi_table, const char * name, const char * err,
									 cpl_size nsamp, cpl_size nwave_new);
cpl_error_code gravi_vis_resamp_phi (cpl_table * oi_table, const char * name, const char * err,
									 cpl_size nsamp, cpl_size nwave_new);
cpl_error_code gravi_vis_smooth_amp (cpl_table * oi_table, const char * name, const char * err,
									 cpl_size nsamp);
cpl_error_code gravi_vis_smooth_phi (cpl_table * oi_table, const char * name, const char * err,
									 cpl_size nsamp);

cpl_error_code gravi_vis_fit_amp (cpl_table * oi_table, const char * name,
                                  const char * err, cpl_size maxdeg);

cpl_error_code gravi_vis_compute_column_mean (cpl_table * out_table,
                                              cpl_table * in_table,
                                              const char * name, int ntel);

cpl_error_code gravi_vis_flag_median (cpl_table * oi_table, const char * data, const char *flag, double value);

cpl_error_code gravi_average_self_visphi(cpl_table * oi_vis_avg,
                                         cpl_table * oi_vis,
                                         cpl_array * wavenumber,
                                         const char * phase_ref, int* cmin, int* cmax, int nrange);

double gdAbacusErrPhi(double x);

/*-----------------------------------------------------------------------------
                                  Function code
 -----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
 * @brief Normal distribution pseudo-random generator
 */
/*----------------------------------------------------------------------------*/

double gravi_randn (void)
{
  int nsamp = 50;
  double samp[] = {0.97446388,  0.78404357,  2.24226141,  1.85507201,  1.10792943,
				   1.34028771, -0.15399594,  0.07316682,  1.60898976,  0.33111245,
				   0.76767625, -2.1252529 ,  0.3898138 ,  2.1218198 ,  0.51703696,
				   0.38451722,  1.07581416, -0.61435275,  1.91926679,  1.10020069,
				   1.82407999,  1.07367663,  0.46105875,  0.45497282,  1.65549611,
				   1.21647974, -0.32725523, -0.36477508,  0.43947414,  1.0242778 ,
				   2.05617949,  1.06163165,  1.24564147,  2.36249995,  0.20676319,
				   1.30886256,  0.7122533 ,  2.28503709,  0.7134141 , -0.19104819,
				   2.9925884 ,  0.95761567,  2.11770457,  0.34763896,  0.30040327,
				   2.3535165 ,  1.65839907,  1.89819461,  1.67480833,  1.11174145};

  /* FIXME: build a better normal random generator !! */
  return samp[rand()%nsamp];
}


/*-----------------------------------------------------------------------------*/

cpl_error_code gravi_array_online_variance(cpl_array * data, cpl_array * mean, cpl_array * variance, int n)
{
  cpl_ensure_code (data,     CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (mean,     CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (variance, CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (n>=0, CPL_ERROR_ILLEGAL_INPUT);
  
  double delta  = 0.0;
  double rdata  = 0.0;
  double rmean  = 0.0;
  double n1 = 1. / ( (double)n + 1.0 );
  double n2 = (double)n / ( (double)n + 1.0 );
  
  cpl_size size = cpl_array_get_size (data);

  /* delta = (x - mean)/(n+1)
	 mean = mean + delta
	 M2 = M2*n/(n+1) + delta*(x - mean) */

  cpl_size w;
  int nv = 0.0;
  for (w = 0; w < size; w ++) {
	/* delta = (x - mean)/(n+1) */
	rdata = cpl_array_get (data, w, &nv);
	rmean = cpl_array_get (mean, w, &nv);
	delta = ( rdata - rmean ) * n1;
	/* mean = mean + delta */
	rmean = rmean + delta;
	cpl_array_set (mean, w, rmean);
	/* M2 = M2*n/(n+1) + delta*(x - mean) */
	cpl_array_set (variance, w,
				   cpl_array_get (variance, w, &nv) * n2 + delta * ( rdata - rmean ));
  }

  int code; 
  if ( (code=cpl_error_get_code()) ) {
	return cpl_error_set_message(cpl_func, code, "Cannot do online variance");
  }

  return CPL_ERROR_NONE;
}

/*-----------------------------------------------------------------------------*/
/**
 * @brief On-line variance of arrays
 * 
 * @param data     list of 4 arrays, all initialized:
 *                 0: current_sample,
 *                 1: running_mean,
 *                 2: first_boot (to keep track if necessary),
 *                 3: running_variance
 * @param n        the sample number
 * @param rephase  if 1, the value are consired being phases
 * 
 * Compute the online-variance of a sample of 'data' with increasing number
 * of sample (n), in-place. Mean and variance should be init to zero at n=0 
 * See https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance 
 */
/*-----------------------------------------------------------------------------*/

cpl_error_code gravi_array_online_variance_res(cpl_array ** data,
											   int n, int rephase)
{
  gravi_msg_function_start(0);
  cpl_ensure_code (data, CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (n>=0, CPL_ERROR_ILLEGAL_INPUT);

  cpl_size size = cpl_array_get_size (data[0]);
  cpl_msg_debug(cpl_func,"Start function");

  /* If first boot, we init the runnning mean and running variance
   * and we store in the first array */
  if ( n == 0 ) {
  	cpl_array_add (data[2], data[0]);
  }

  /* Recenter phase around the mean phase before computing its VARIANCE */
  if (rephase) {
	for (int w=0; w<size; w++ ) {
	  cpl_array_set (data[0], w,
					 carg( cexp (1*I* (cpl_array_get(data[0], w, NULL) -
									   cpl_array_get(data[2], w, NULL))) ) );
	}
  }

  /* Run the gravi_online */
  gravi_array_online_variance(data[0], data[1], data[3], n);

  /* Free the current boot to prepare for next integration */
  cpl_array_fill_window (data[0], 0, size, 0.0);

  CPLCHECK_MSG ("Error in online variance");
  
  gravi_msg_function_exit(0);
  return CPL_ERROR_NONE;
}

/*-----------------------------------------------------------------------------*/
/**
 * @brief Average the flux of all DITs into a final, averaged value
 * 
 * @param oi_flux_avg:  already allocated OI_FLUX table to be filled
 * 
 * @param oi_flux:      input OI_FLUX with all individual DITs
 * @param nseg:         number of segment for the boostrap
 * @param nboot:        number of boostrap
 * @param tel:          beam to consider (0..3)
 *
 * Average the flux and compute the final flux.
 * The errors are computed with the boostraping method.
 */
/*-----------------------------------------------------------------------------*/

cpl_error_code gravi_flux_average_bootstrap(cpl_table * oi_flux_avg,
											cpl_table * oi_flux,
											int nboot,
					    double outlier_threshold)
{
   gravi_msg_function_start(0);
   cpl_ensure_code (oi_flux_avg, CPL_ERROR_ILLEGAL_OUTPUT);
   cpl_ensure_code (oi_flux,     CPL_ERROR_NULL_INPUT);
   cpl_ensure_code (nboot>0, CPL_ERROR_ILLEGAL_INPUT);

   /* parameters */
  int nv = 0, ntel = 4;
  cpl_size nrow = cpl_table_get_nrow (oi_flux) / ntel;
  cpl_size nwave = cpl_table_get_column_depth (oi_flux, "FLUX");

  /* Pointer to columns, to speed-up */
  cpl_array ** pFLUX     = cpl_table_get_data_array (oi_flux, "FLUX");
  cpl_array ** pFLUXERR  = cpl_table_get_data_array (oi_flux, "FLUXERR");
  cpl_array ** pFLAG     = cpl_table_get_data_array (oi_flux, "FLAG");
  double * pINTTIME = cpl_table_get_data_double (oi_flux, "INT_TIME");
  double * pMJD     = cpl_table_get_data_double (oi_flux, "MJD");
  CPLCHECK_MSG ("Cannot get the data");

  /* Loop on tel */
  for (cpl_size tel = 0; tel < ntel; tel++) {
      
  /* Tel for base and base for closure */
  cpl_size nvalid = 0;

  /* Arrays to store the final, integrated quantities 
   * 0: current boot, 1: running_mean, 2: first boot, 3: variance */
  cpl_array **flux_res = gravi_array_new_list (4, CPL_TYPE_DOUBLE, nwave);
  double total_exptime = 0.0, mjd_avg = 0.0;

  /* 
   * (0) Optimize the number of segment
   */
  
  /* Get the number of non-rejected frames */
  int * flag = cpl_malloc( sizeof(int) * nrow );
  for ( int row=0 ; row<nrow; row++ ) {
	flag[row] = 0; nvalid++;
  }
  
  /* Build an optimal number of segment and nrow_per_segment */
  cpl_size nrow_per_seg = CPL_MAX(nvalid / CPL_MIN (nrow, 100), 1);
  cpl_size nseg = nvalid / nrow_per_seg;

  /* Ensure there are at least 5 samples to bootstrap on,
   * if no add montecarlo samples */
  cpl_size nsamp = 5, nmontecarlo = CPL_MAX (nsamp - nseg, 0);

  cpl_msg_info ("Stat", "%6lld valid frames over %6lld (%5.1f%%), make %4lld seg. of %5lld (miss %lld), add %lld MonteCarlo",
				nvalid, nrow, (double)nvalid/(double)nrow*100.0,
				nseg, nrow_per_seg, nvalid - nseg*nrow_per_seg, nmontecarlo);
  
  /* Case we have at least one valid frame */
  if ( nvalid > 0 ) {

    /*
     * (1) Pre-integration over segment, to bootstrap on less statistic 
     */
	
	cpl_array **flux = gravi_array_new_list (nseg + nmontecarlo, CPL_TYPE_DOUBLE, nwave);
  
	/* Loop on segment */
	cpl_size row = -1;
	for ( int seg = 0 ; seg < nseg + nmontecarlo; seg ++ ) {
	  cpl_msg_debug(cpl_func,"pre integration of seg %d start with row %lld", seg, row);

	  /* Find nrow_per_seg valid frame to integrate in this segment */
	  cpl_size r = 0;
	  while ( r < nrow_per_seg ) {
		row = (row + 1) % nrow;
		if ( flag[row] ) {continue;} else {r++;}

        /* Get indices */
        cpl_size rtel = row * ntel + tel;
        
		/* Compute the total integration time.
		 * Do not integrate for the MonteCarlo samples */
		if (seg < nseg) {
		  total_exptime += pINTTIME[rtel];
		  mjd_avg += pMJD[rtel] * pINTTIME[rtel];
		}
		
        /* fast-no-CPL integration: get pointers on data */
		double * tflux = cpl_array_get_data_double (flux[seg]);

		/* Loop on wave */
		for ( int w=0; w<nwave; w++ ) {
          double FLUX    = cpl_array_get (pFLUX[rtel], w, NULL);
          double FLUXERR = cpl_array_get (pFLUXERR[rtel], w, NULL);
	  int outlier_flag = cpl_array_get (pFLAG[rtel], w, NULL);
            CPLCHECK_MSG ("Cannot get data");

	    /* Reject outlier */
	    if (outlier_flag) {
	      FLUX = 0.0;
	    }
		  /* Add noise if this is a Monte Carlo sample. */
		  if ( seg > nseg-1 ) {
			FLUX += 2 * FLUXERR * gravi_randn();
		  }
          
		  /* flux = < FLUX > over the frames in [e] */
		  tflux[w] += FLUX;

		} /* End loop on wave */
	  } /* End loop on rows in this segment */
    }/* End loop on segments */

    /*
     * (2) Compute the variance by bootstraping on the segments
     */
	
    /* Loop on bootstramp to compute the avg and the 
     * variance by the bootstraping methode */
    srand(1);
    
    for ( int boot = 0 ; boot < nboot ; boot ++ ) {
    	cpl_msg_debug(cpl_func,"Bootstrap %d over %d", boot+1, nboot);
			  
		/* Integrate nseg segment randomly selected.
		 * This loop is vectorialized in spectral direction */
		for (int rowb = 0; rowb < nseg; rowb ++){
				  
    	  /* For the first bootstrap, we use all observed samples
		   * For the others, we also includes the possible montecarlo
    	   * FIXME: verify the uniformity of rand for small nrows */
		  int rows;
		  if (boot == 0 )  rows = rowb;
		  else             rows = rand()%(nseg+nmontecarlo);
				  
		  /* Integrate the flux of selected segments */
		  cpl_array_add (flux_res[0], flux[rows]);
		}
		/* End loop on selected segments */

		/* Compute the VARIANCE over the bootstraped samples with the 'online_variance' algorithm.
		 * See https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance */
		gravi_array_online_variance_res (flux_res, boot, 0);
		
    	CPLCHECK_MSG("while computing the variances over the bootstrap");
	}
	/* End Loop on bootstrap */
			  
    /* Free list of segments and running mean */
	FREELOOP (cpl_array_delete, flux, nseg + nmontecarlo);

	/* Convert variance from bootstrap to RMS */
	cpl_msg_debug(cpl_func,"Put the RMS over bootstrap");	
	cpl_array_power (flux_res[3], 0.5);
	CPLCHECK_MSG("while converting variance -> rms");

	/* Normalise integration */
	mjd_avg /= total_exptime;
	
  }

  /* 
   * Step (1,2) in case there are no valid frames at all
   */
  if (nvalid == 0) {
	cpl_msg_debug (cpl_func,"Not valid frames, force zero and infinit RMS");
	cpl_array_fill_window (flux_res[3], 0, nwave, 1e10);
	mjd_avg = cpl_table_get_column_mean (oi_flux, "MJD");
  }

  /* Renormalise by the number of outliers.
     Also flag channels with too much outliers. */
  cpl_array ** pflag  = cpl_table_get_data_array (oi_flux_avg,  "FLAG");
  cpl_array * array = gravi_table_get_column_sum_array (oi_flux, "FLAG", tel, ntel);
  CPLCHECK_MSG("cannot get data");
  
  for (int wave = 0; wave < nwave; wave++) {
    double value = cpl_array_get (array, wave, NULL) / nrow;
    if (value < 1.0) {
      cpl_array_set (flux_res[2], wave, cpl_array_get (flux_res[2], wave, NULL) / (1.0 - value));
    }
    if (value > outlier_threshold) {
      cpl_array_set (pflag[tel], wave, 1);
    }
  }
  FREE (cpl_array_delete, array);

  /* 
   * (3) Save the results on the oi_flux_avg tables 
   */

  cpl_table_set_array (oi_flux_avg, "FLUX", tel, flux_res[2]);
  cpl_table_set_array (oi_flux_avg, "FLUXERR", tel, flux_res[3]);
  CPLCHECK_MSG("filling FLUX and FLUXERR");

  /* Flag the data with >100% error */
  gravi_vis_flag_relative_threshold (oi_flux_avg, "FLUXERR", "FLUX", "FLAG", 1.0);
  gravi_vis_flag_lower (oi_flux_avg, "FLUXERR", "FLAG", 0.0);
  CPLCHECK_MSG("cannot flag baddata data");

  /* Compute the total integration time */
  cpl_msg_debug(cpl_func,"Total integration time = %.3f s", total_exptime);
  cpl_table_set_double (oi_flux_avg, "INT_TIME", tel, total_exptime);
  cpl_table_set_double (oi_flux_avg, "MJD", tel, mjd_avg);
  cpl_table_set (oi_flux_avg, "NVALID", tel, nvalid);
  cpl_table_set (oi_flux_avg, "NDIT", tel, nrow);

  /* Set the TARGET_ID and STA_INDEX */
  cpl_table_set_int (oi_flux_avg, "TARGET_ID", tel, cpl_table_get_int (oi_flux, "TARGET_ID", tel, &nv));
  cpl_table_set_int (oi_flux_avg, "STA_INDEX", tel, cpl_table_get_int (oi_flux, "STA_INDEX", tel, &nv));

  FREELOOP (cpl_array_delete, flux_res, 4);
  cpl_free(flag);

  } /* End loop on tel */
  
  gravi_msg_function_exit(0);
  return CPL_ERROR_NONE;
}

/*-----------------------------------------------------------------------------*/
/**
 * @brief Average the closure-phase of all DITs into a final, averaged value
 * 
 * @param oi_t3_avg:  already allocated OI_T3 table to be filled
 * 
 * @param oi_vis:     input OI_VIS with all individual DITs
 * @param oi_flux:    input OI_FLUX with all individual DITs
 * @param nseg:       number of segment for the boostrap
 * @param nboot:      number of boostrap
 * @param closure:    triplet to consider (0..3)
 *
 * Average the bispectrum and compute the final closure-phase, considering the 
 * rejection flags. The errors are computed with the boostraping method.
 */
/*-----------------------------------------------------------------------------*/

cpl_error_code gravi_t3_average_bootstrap(cpl_table * oi_t3_avg,
										  cpl_table * oi_vis,
										  cpl_table * oi_flux,
										  int nboot,
										  int use_vFactor,
                                          int use_pFactor,
					  double outlier_threshold)
{
  gravi_msg_function_start(0);
  cpl_ensure_code (oi_t3_avg, CPL_ERROR_ILLEGAL_OUTPUT);
  cpl_ensure_code (oi_vis,    CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (oi_flux,   CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (nboot>0,    CPL_ERROR_ILLEGAL_INPUT);

  /* Tel for base and base for closure 
   * Assume all the observations are made with the same array geometry.
   * sta_index_t3[0] <-> tel1[clo[j][0]] = tel1[clo[j][2]]
   * sta_index_t3[1] <-> tel1[clo[j][1]] = tel2[clo[j][0]]
   * sta_index_t3[2] <-> tel2[clo[j][2]] = tel2[clo[j][1]]
   */
  int nv = 0, nbase = 6, ntel = 4, nclo = 4;

  cpl_size nrow = cpl_table_get_nrow (oi_vis) / nbase;
  cpl_size nwave = cpl_table_get_column_depth (oi_vis, "VISDATA");

  /* Pointer to column, to speed-up */
  cpl_array ** pVISDATA = cpl_table_get_data_array (oi_vis, "VISDATA");
  cpl_array ** pVISERR  = cpl_table_get_data_array (oi_vis, "VISERR");
  cpl_array ** pFLUX    = cpl_table_get_data_array (oi_flux, "FLUX");
  cpl_array ** pFLAG    = cpl_table_get_data_array (oi_vis, "FLAG");
  double * pINTTIME = cpl_table_get_data_double (oi_vis, "INT_TIME");
  double * pMJD     = cpl_table_get_data_double (oi_vis, "MJD");
  double * pUCOORD  = cpl_table_get_data_double (oi_vis, "UCOORD");
  double * pVCOORD  = cpl_table_get_data_double (oi_vis, "VCOORD");
  cpl_array ** pVFACTOR = use_vFactor?cpl_table_get_data_array (oi_vis, "V_FACTOR"):NULL;
  double * pPFACTOR = use_pFactor?cpl_table_get_data_double (oi_vis, "P_FACTOR"):NULL;
  CPLCHECK_MSG ("Cannot get the data");

  /* Loop on closure */
  for (cpl_size closure = 0; closure < nclo; closure++) {

  cpl_size nvalid = 0;
  int base0 = GRAVI_CLO_BASE[closure][0];
  int base1 = GRAVI_CLO_BASE[closure][1];
  int base2 = GRAVI_CLO_BASE[closure][2];
  int ctel0 = GRAVI_CLO_TEL[closure][0];
  int ctel1 = GRAVI_CLO_TEL[closure][1];
  int ctel2 = GRAVI_CLO_TEL[closure][2];
  
  /* Arrays to store the final, integrated quantities 
   * 0: current boot, 1: running_mean, 2: first boot, 3: variance */
  cpl_array **t3Amp_res = gravi_array_new_list (4, CPL_TYPE_DOUBLE, nwave);
  cpl_array **t3Phi_res = gravi_array_new_list (4, CPL_TYPE_DOUBLE, nwave);
  double total_exptime = 0.0, mjd_avg = 0.0;
  double u1Coord = 0.0, v1Coord = 0.0, u2Coord = 0.0, v2Coord = 0.0;
  
  /* 
   * (0) Optimize the number of segment
   */
  
  /* Get the number of non-rejected frames */
  int * flag = cpl_table_get_data_int (oi_vis, "REJECTION_FLAG");
  cpl_ensure_code (flag, CPL_ERROR_ILLEGAL_INPUT);
  int * flagclo = cpl_malloc( sizeof(int) * nrow );
  for ( int row=0 ; row<nrow; row++ ) {
	flagclo[row] = flag[row * nbase + base0] + flag[row * nbase + base1] + flag[row * nbase + base2];
	if ( flagclo[row] == 0 ) nvalid++;
  }
  
  /* Build an optimal number of segment and nrow_per_segment */
  cpl_size nrow_per_seg = CPL_MAX (nvalid / CPL_MIN (nrow, 100), 1);
  cpl_size nseg = nvalid / nrow_per_seg;

  /* Ensure there are at least 5 samples to bootstrap on,
   * if no add montecarlo samples */
  cpl_size nsamp = 5, nmontecarlo = CPL_MAX (nsamp - nseg, 0);

  cpl_msg_info ("Stat", "%6lld valid frames over %6lld (%5.1f%%), make %4lld seg. of %5lld (miss %lld), add %lld MonteCarlo",
				nvalid, nrow, (double)nvalid/(double)nrow*100.0,
				nseg, nrow_per_seg, nvalid - nseg*nrow_per_seg, nmontecarlo);
  
  /* Case we have at least one valid frame */
  if ( nvalid > 0 ) {

    /*
     * (1) Pre-integration over segment, to bootstrap on less statistic 
     */
	
	cpl_array **bisp = gravi_array_new_list (nseg + nmontecarlo, CPL_TYPE_DOUBLE_COMPLEX, nwave);
    cpl_array **F012 = gravi_array_new_list (nseg + nmontecarlo, CPL_TYPE_DOUBLE, nwave);

	/* Loop on segment */
	cpl_size row = -1;
	for ( int seg = 0 ; seg < nseg + nmontecarlo ; seg ++ ) {
	  cpl_msg_debug(cpl_func,"pre integration of seg %d start with row %lld", seg, row);

	  /* Find nrow_per_seg valid frame to integrate in this segment */
	  cpl_size r = 0;
	  while ( r < nrow_per_seg ) {
		row = (row + 1) % nrow;
		if ( flagclo[row] ) {continue;} else {r++;}

        /* Get indices */
        cpl_size rbase0 = row * nbase + base0;
        cpl_size rbase1 = row * nbase + base1;
        cpl_size rbase2 = row * nbase + base2;
    
		/* Compute the total integration time.
		 * Do not integrate for the MonteCarlo samples */
		if (seg < nseg) {
		  total_exptime += pINTTIME[rbase0];
		  mjd_avg += pMJD[rbase0] * pINTTIME[rbase0];
		  u1Coord += pUCOORD[rbase0] * pINTTIME[rbase0];
		  v1Coord += pVCOORD[rbase0] * pINTTIME[rbase0];
		  u2Coord += pUCOORD[rbase1] * pINTTIME[rbase0];
		  v2Coord += pVCOORD[rbase1] * pINTTIME[rbase0];
		}
		
        /* fast-no-CPL integration: get pointers on data */
		double complex * tbisp = cpl_array_get_data_double_complex (bisp[seg]);
		double *tF012 = cpl_array_get_data_double (F012[seg]);
		CPLCHECK_MSG ("Cannot get data");
        
        double PFACTOR0 = (use_pFactor?pPFACTOR[rbase0]:1.0);
        double PFACTOR1 = (use_pFactor?pPFACTOR[rbase1]:1.0);
        double PFACTOR2 = (use_pFactor?pPFACTOR[rbase2]:1.0);
		CPLCHECK_MSG ("Cannot get FACTOR data");

		/* Loop on wave */
		for ( int w=0; w<nwave; w++ ) {
          double complex Vis0 = cpl_array_get_complex (pVISDATA[rbase0], w, NULL);
          double complex Vis1 = cpl_array_get_complex (pVISDATA[rbase1], w, NULL);
          double complex Vis2 = cpl_array_get_complex (pVISDATA[rbase2], w, NULL);
          double complex VisErr0 = cpl_array_get_complex (pVISERR[rbase0], w, NULL);
          double complex VisErr1 = cpl_array_get_complex (pVISERR[rbase1], w, NULL);
          double complex VisErr2 = cpl_array_get_complex (pVISERR[rbase2], w, NULL);
          double F0 = cpl_array_get (pFLUX[row * ntel + ctel0], w, NULL);
          double F1 = cpl_array_get (pFLUX[row * ntel + ctel1], w, NULL);
          double F2 = cpl_array_get (pFLUX[row * ntel + ctel2], w, NULL);
          double VFACTOR0 = (use_vFactor?cpl_array_get (pVFACTOR[rbase0], w, NULL):1.0);
          double VFACTOR1 = (use_vFactor?cpl_array_get (pVFACTOR[rbase1], w, NULL):1.0);
          double VFACTOR2 = (use_vFactor?cpl_array_get (pVFACTOR[rbase2], w, NULL):1.0);
	  int outlier_flag0 = cpl_array_get (pFLAG[rbase0], w, NULL);
	  int outlier_flag1 = cpl_array_get (pFLAG[rbase1], w, NULL);
	  int outlier_flag2 = cpl_array_get (pFLAG[rbase2], w, NULL);

	  /* Reject outlier */
	  if (outlier_flag0 || outlier_flag1 || outlier_flag2) {
	    Vis0 = 0.0+0.0*I; 
	    Vis1 = 0.0+0.0*I;
	    Vis2 = 0.0+0.0*I;
	    VisErr0 = 0.0+0.0*I;
	    VisErr1 = 0.0+0.0*I;
	    VisErr2 = 0.0+0.0*I;
	    F0 = 0.0;
	    F1 = 0.0;
	    F2 = 0.0;
	  }
	  
		  /* Add noise if this is a Monte Carlo sample.
		   * APPROX: Noise is only added to the coherent fluxes */
		  if ( seg > nseg-1 ) {
			Vis0 += 1.*I * cimag(VisErr0) * gravi_randn();
			Vis0 +=   1. * creal(VisErr0) * gravi_randn();
			Vis1 += 1.*I * cimag(VisErr1) * gravi_randn();
			Vis1 +=   1. * creal(VisErr1) * gravi_randn();
			Vis2 += 1.*I * cimag(VisErr2) * gravi_randn();
			Vis2 +=   1. * creal(VisErr2) * gravi_randn();
		  }
		  
		  /* bisp = < v1*v2*conj(v3) > over the frames in [e^3] */
		  tbisp[w] += Vis0 * Vis1 * conj (Vis2);

		  /* F012 = < F0*F1*F2 > over the frames in [e^3] 
           * corrected from expected visibility losses */
          tF012[w] += F0 * F1 * F2 *
                      sqrt (CPL_MAX (VFACTOR0 * VFACTOR1 * VFACTOR2 *
                                     PFACTOR0 * PFACTOR1 * PFACTOR2, 0.0));
		  
		} /* End loop on wave */
	  } /* End loop on rows in this segment */
    }/* End loop on segments */

    /*
     * (2) Compute the variance by bootstraping on the segments
     */
	
    /* Loop on bootstramp to compute the avg and the 
     * variance by the bootstraping methode */
    srand(1);
    
    for ( int boot = 0 ; boot < nboot ; boot ++ ) {
    	cpl_msg_debug(cpl_func,"Bootstrap %d over %d", boot+1, nboot);
			  
		/* Init the integration of nseg segment randomly selected */
		cpl_array * bisp_boot = gravi_array_init_double_complex (nwave, 0.0 + I*0.0);
		cpl_array * f012_boot = gravi_array_init_double (nwave, 0.0);
				
		/* Integrate nseg segment randomly selected.
		 * This loop is vectorialized in spectral direction */
		for (int rowb = 0; rowb < nseg; rowb ++){
		  
    	  /* For the first bootstrap, we use all observed samples
		   * For the others, we also includes the possible montecarlo
    	   * FIXME: verify the uniformity of rand for small nrows */
		  int rows;
		  if (boot == 0 )  rows = rowb;
		  else             rows = rand()%(nseg+nmontecarlo);
				  
		  /* Integrate the bispectre and flux of selected segments */
		  cpl_array_add (bisp_boot, bisp[rows]);
		  cpl_array_add (f012_boot, F012[rows]);
		}
		/* End loop on selected segments */

        /* Make sure the geometric flux is not null */
        gravi_array_threshold_min (f012_boot, 1e-15);
        
		/* Compute the argument and the module of the bispectrum */
		FREE (cpl_array_delete, t3Amp_res[0]);
		t3Amp_res[0] = cpl_array_duplicate(bisp_boot);
		cpl_array_abs (t3Amp_res[0]);
		cpl_array_divide (t3Amp_res[0], f012_boot);
		
		FREE (cpl_array_delete, t3Phi_res[0]);
		t3Phi_res[0] = cpl_array_duplicate(bisp_boot);
		cpl_array_arg (t3Phi_res[0]);

		/* Compute the VARIANCE over the bootstraped samples with the 'online_variance' algorithm.
		 * See https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance */
		gravi_array_online_variance_res (t3Amp_res, boot, 0);
		gravi_array_online_variance_res (t3Phi_res, boot, 1);
		
		FREE (cpl_array_delete, bisp_boot);
		FREE (cpl_array_delete, f012_boot);
    	CPLCHECK_MSG("while computing the variances over the bootstrap");
	}
	/* End Loop on bootstrap */
			  
    /* Free list of segments and running mean */
	FREELOOP (cpl_array_delete, bisp, nseg + nmontecarlo);
	FREELOOP (cpl_array_delete, F012, nseg + nmontecarlo);

	/* Convert variance from bootstrap to RMS */
	cpl_msg_debug(cpl_func,"Put the RMS over bootstrap");
	cpl_array_power (t3Phi_res[3], 0.5);
	cpl_array_power (t3Amp_res[3], 0.5);
	CPLCHECK_MSG("while converting variance -> rms");

	/* Normalizase integration */
	mjd_avg /= total_exptime;
	u1Coord /= total_exptime;
	v1Coord /= total_exptime;
	u2Coord /= total_exptime;
	v2Coord /= total_exptime;
	
  } 

  /* 
   * Step (1,2) in case there are no valid frames at all
   */
  if (nvalid == 0) {
	cpl_msg_debug (cpl_func,"Not valid frames, force zero and infinit RMS");
	cpl_array_fill_window (t3Amp_res[3], 0, nwave, 1e10);
	cpl_array_fill_window (t3Phi_res[3], 0, nwave, 1e10);
	mjd_avg = cpl_table_get_column_mean (oi_vis, "MJD");
  }

  /* Flag channels with too much outliers. */
  cpl_array ** pflag  = cpl_table_get_data_array (oi_t3_avg,  "FLAG");
  CPLCHECK_MSG("cannot get data");

  cpl_array * array = cpl_array_new (nwave, CPL_TYPE_INT);
  cpl_array_fill_window (array, 0, nwave, 0.0);
  int * sum_flag = cpl_array_get_data_int (array);
  
  for (int row=0; row < nrow; row++) {
      cpl_array * arr0 = pFLAG[row * nbase + base0];
      cpl_array * arr1 = pFLAG[row * nbase + base1];
      cpl_array * arr2 = pFLAG[row * nbase + base2];
      for (int wave = 0; wave < nwave; wave++) {
	if (cpl_array_get (arr0, wave, NULL) ||
	    cpl_array_get (arr1, wave, NULL) ||
	    cpl_array_get (arr2, wave, NULL) )
	  sum_flag[wave] ++;
      }
  }
  
  for (int wave = 0; wave < nwave; wave++) {
    double value = cpl_array_get (array, wave, NULL) / nrow;
    if (value > outlier_threshold) {
      cpl_array_set (pflag[closure],  wave, 1);
    }
  }
  FREE (cpl_array_delete, array);

  /* 
   * (3) Save the results on the oi_t3_avg tables 
   */
	
  /* Save the cloture amplitude on the oi_T3 tables  */
  cpl_table_set_array (oi_t3_avg, "T3AMP", closure, t3Amp_res[2]);
  cpl_table_set_array (oi_t3_avg, "T3AMPERR", closure, t3Amp_res[3]);
  CPLCHECK_MSG("filling T3AMP");
			  
  /* Save the cloture phase on the oi_T3 tables  */
  gravi_table_set_array_phase (oi_t3_avg, "T3PHI", closure, t3Phi_res[2]);
  gravi_table_set_array_phase (oi_t3_avg, "T3PHIERR", closure, t3Phi_res[3]);
  CPLCHECK_MSG("filling T3PHI");

  /* Flag the data with >100% error or >60deg error or negative errors */
  gravi_vis_flag_threshold (oi_t3_avg, "T3PHIERR", "FLAG", 60.0);
  gravi_vis_flag_threshold (oi_t3_avg, "T3AMPERR", "FLAG", 1.0);
  gravi_vis_flag_lower (oi_t3_avg, "T3AMPERR", "FLAG", 0.0);
  gravi_vis_flag_median (oi_t3_avg, "T3PHIERR", "FLAG", 5.0);
  CPLCHECK_MSG("cannot flag baddata data");
  
  /* Compute the total integration time and MJD */
  cpl_msg_debug(cpl_func,"Total integration time = %.3f s", total_exptime);
  cpl_table_set_double (oi_t3_avg, "INT_TIME", closure, total_exptime);
  cpl_table_set_double (oi_t3_avg, "MJD", closure, mjd_avg);
  cpl_table_set_double (oi_t3_avg, "U1COORD", closure, u1Coord);
  cpl_table_set_double (oi_t3_avg, "V1COORD", closure, v1Coord);
  cpl_table_set_double (oi_t3_avg, "U2COORD", closure, u2Coord);
  cpl_table_set_double (oi_t3_avg, "V2COORD", closure, v2Coord);
  cpl_table_set (oi_t3_avg, "NVALID", closure, nvalid);
  cpl_table_set (oi_t3_avg, "NDIT", closure, nrow);
  
  /* Set the TARGET_ID */
  cpl_table_set_int (oi_t3_avg, "TARGET_ID", closure, cpl_table_get_int (oi_vis, "TARGET_ID", base0, &nv));

  /* Set STA_INDEX */
  cpl_array * sta_index = cpl_array_new (3, CPL_TYPE_INT);
  cpl_array_set_int (sta_index, 0, cpl_table_get_int (oi_flux,"STA_INDEX", ctel0, &nv));
  cpl_array_set_int (sta_index, 1, cpl_table_get_int (oi_flux,"STA_INDEX", ctel1, &nv));
  cpl_array_set_int (sta_index, 2, cpl_table_get_int (oi_flux,"STA_INDEX", ctel2, &nv));
  cpl_table_set_array (oi_t3_avg, "STA_INDEX", closure, sta_index);
  FREE (cpl_array_delete, sta_index);
  
  /* Free the aggregate flags for closures */
  FREELOOP (cpl_array_delete, t3Phi_res, 4);
  FREELOOP (cpl_array_delete, t3Amp_res, 4);
  FREE (cpl_free, flagclo);
  
  } /* End loop on closure */
  
  gravi_msg_function_exit(0);
  return CPL_ERROR_NONE;
}

/*-----------------------------------------------------------------------------*/
/**
 * @brief Average the visibility of all DITs into a final, averaged value
 * 
 * @param oi_vis_avg:     already allocated OI_VIS table to be filled
 * @param oi_vis2_avg:    already allocated OI_VIS2 table to be filled
 * 
 * @param oi_vis:         input OI_VIS with all individual DITs
 * @param oi_flux:        input OI_FLUX with all individual DITs
 * @param nseg:           number of segment for the boostrap
 * @param nboot:          number of boostrap
 * @param base:           base to consider (0..5)
 * @param phase_ref:      phase used to rephase the DITs
 * @param use_vFactor:    use the VFACTOR to compensate for vis. losses
 *
 * Average the coherent flux and photometric fluex, then normalize and
 * compute the final visibilities, considering the rejection flags.
 * The errors are computed with the boostraping method.
 */
/*-----------------------------------------------------------------------------*/

cpl_error_code gravi_vis_average_bootstrap (cpl_table * oi_vis_avg,
											cpl_table * oi_vis2_avg,
											cpl_table * oi_vis,
											int nboot,
											const char * phase_ref,
											int use_vFactor,
                                            int use_pFactor,
					    int use_debiasing,
					    double outlier_threshold)
{
  gravi_msg_function_start(0);
  cpl_ensure_code (oi_vis_avg,  CPL_ERROR_ILLEGAL_OUTPUT);
  cpl_ensure_code (oi_vis2_avg, CPL_ERROR_ILLEGAL_OUTPUT);
  cpl_ensure_code (oi_vis,      CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (nboot>0, CPL_ERROR_ILLEGAL_INPUT);

  /* Parameters */
  int nv = 0, nbase = 6;
  cpl_size nrow = cpl_table_get_nrow (oi_vis) / nbase;
  cpl_size nwave = cpl_table_get_column_depth (oi_vis, "VISDATA");

  /* Pointer to columns, to speed-up */
  cpl_array ** pVISDATA = cpl_table_get_data_array (oi_vis, "VISDATA");
  cpl_array ** pVISERR  = cpl_table_get_data_array (oi_vis, "VISERR");
  cpl_array ** pFNORM   = cpl_table_get_data_array (oi_vis, "F1F2");
  cpl_array ** pFLAG    = cpl_table_get_data_array (oi_vis, "FLAG");
  double * pINTTIME = cpl_table_get_data_double (oi_vis, "INT_TIME");
  double * pMJD     = cpl_table_get_data_double (oi_vis, "MJD");
  double * pUCOORD  = cpl_table_get_data_double (oi_vis, "UCOORD");
  double * pVCOORD  = cpl_table_get_data_double (oi_vis, "VCOORD");
  cpl_array ** pVFACTOR = use_vFactor?cpl_table_get_data_array (oi_vis, "V_FACTOR"):NULL;
  double * pPFACTOR = use_pFactor?cpl_table_get_data_double (oi_vis, "P_FACTOR"):NULL;
  CPLCHECK_MSG ("Cannot get the data");

  /* Get the reference phase */
  cpl_array ** pPHASEREF = NULL;
  if (phase_ref && strcmp (phase_ref, "NONE"))
      pPHASEREF = cpl_table_get_data_array (oi_vis, phase_ref);
  CPLCHECK_MSG ("Cannot get the reference phase data (column missing?)");
  
  /* Loop on base */
  for (cpl_size base = 0; base < nbase; base ++) {

  /* Tel for base */
  cpl_size nvalid = 0;

  /* Arrays to store the final, integrated quantities 
   * 0: current boot, 1: running_mean, 2: first boot, 3: variance */
  cpl_array **visR_res = gravi_array_new_list (4, CPL_TYPE_DOUBLE, nwave);
  cpl_array **visI_res = gravi_array_new_list (4, CPL_TYPE_DOUBLE, nwave);
  cpl_array **vis2_res = gravi_array_new_list (4, CPL_TYPE_DOUBLE, nwave);
  cpl_array **visAmp_res = gravi_array_new_list (4, CPL_TYPE_DOUBLE, nwave);
  cpl_array **visPhi_res = gravi_array_new_list (4, CPL_TYPE_DOUBLE, nwave);
  double total_exptime = 0.0, mjd_avg = 0.0;
  double uCoord = 0.0, vCoord = 0.0; 

  /* 
   * (0) Optimize the number of segment
   */
  
  /* Get the number of non-rejected frames */
  int * flag = cpl_table_get_data_int (oi_vis, "REJECTION_FLAG");
  cpl_ensure_code (flag, CPL_ERROR_ILLEGAL_INPUT);
  for ( int row=0 ; row<nrow; row++ ) if ( flag[row * nbase + base] == 0 ) nvalid++;

  /* Build an optimal number of segment and nrow_per_segment */
  cpl_size nrow_per_seg = CPL_MAX (nvalid / CPL_MIN (nrow, 100), 1);
  cpl_size nseg = nvalid / nrow_per_seg;

  /* Ensure there are at least 5 samples to bootstrap on,
   * if no add montecarlo samples */
  cpl_size nsamp = 5, nmontecarlo = CPL_MAX (nsamp - nseg, 0);

  cpl_msg_info ("Stat", "%6lld valid frames over %6lld (%5.1f%%), make %4lld seg. of %5lld (miss %lld), add %lld MonteCarlo",
				nvalid, nrow, (double)nvalid/(double)nrow*100.0,
				nseg, nrow_per_seg, nvalid - nseg*nrow_per_seg, nmontecarlo);

  /* Case we have at least one valid frame */
  if ( nvalid > 0 ) {

    /*
     * (1) Pre-integration over segment, to bootstrap on less statistic 
     */
    cpl_array **visR  = gravi_array_new_list (nseg + nmontecarlo, CPL_TYPE_DOUBLE, nwave);
    cpl_array **visI  = gravi_array_new_list (nseg + nmontecarlo, CPL_TYPE_DOUBLE, nwave);			
    cpl_array **POWER = gravi_array_new_list (nseg + nmontecarlo, CPL_TYPE_DOUBLE, nwave);			
    cpl_array **F12   = gravi_array_new_list (nseg + nmontecarlo, CPL_TYPE_DOUBLE, nwave);			
    cpl_array **F1F2  = gravi_array_new_list (nseg + nmontecarlo, CPL_TYPE_DOUBLE, nwave);

    /* Loop on segment */
    cpl_size row = -1;
    for ( int seg = 0 ; seg < nseg + nmontecarlo ; seg ++ ) {
	    cpl_msg_debug(cpl_func,"pre integration of seg %d start with row %lld", seg, row);
	  
    	/* Find nrow_per_seg valid frame to integrate in this segment */
    	cpl_size r = 0;
    	while ( r < nrow_per_seg ) {
    	  row = (row + 1) % nrow;
    	  if ( flag[row * nbase + base] ) {continue;} else {r++;}

          /* Get indices */
          cpl_size rbase = row * nbase + base;
    
    	  /* Compute the total integration time.
		   * Do not integrate for the MonteCarlo samples */
		  if (seg < nseg) {
			total_exptime += pINTTIME[rbase];
			mjd_avg += pMJD[rbase] * pINTTIME[rbase];
			uCoord += pUCOORD[rbase] * pINTTIME[rbase];
			vCoord += pVCOORD[rbase] * pINTTIME[rbase];
		  }
		  
    	  /* fast-no-CPL integration: get pointers on data */
    	  double *tR  = cpl_array_get_data_double (visR[seg]);
    	  double *tI  = cpl_array_get_data_double (visI[seg]);
    	  double *tP  = cpl_array_get_data_double (POWER[seg]);
    	  double *tF1F2 = cpl_array_get_data_double (F1F2[seg]);
    	  double *tF12  = cpl_array_get_data_double (F12[seg]);
		  CPLCHECK_MSG ("Cannot get data");

          double  PFACTOR = (use_pFactor?pPFACTOR[rbase]:1.0);
          CPLCHECK_MSG ("Cannot get FACTOR data");
          
    	  /* Loop on wave */
    	  for (int w = 0; w < nwave; w++) {
            double VFACTOR = use_vFactor?cpl_array_get (pVFACTOR[rbase], w, NULL):1.0;
            double PHASEREF = pPHASEREF?cpl_array_get (pPHASEREF[rbase], w, NULL):0.0;
            double FNORM = cpl_array_get (pFNORM[rbase], w, NULL);
            double complex Vis  = cpl_array_get_complex (pVISDATA[rbase], w, NULL);
            double complex VErr = cpl_array_get_complex (pVISERR[rbase], w, NULL);
	    double mR  = creal (Vis);
	    double mI  = cimag (Vis);
	    double eR = creal (VErr);
	    double eI = cimag (VErr);
            int outlier_flag = cpl_array_get (pFLAG[rbase], w, NULL);
            CPLCHECK_MSG ("Cannot get data");

	    /* Reject outlier */
	    if (outlier_flag) {
	      mR = 0.0; mI = 0.0;
	      eR = 0.0; eI = 0.0;
	      FNORM = 0.0;
	    }

			/* Add noise if this is a Monte Carlo sample.
			 * APPROX: Noise is only added to the coherent flux */
			if ( seg > nseg-1 ) {
			  mR += 2 * eR * gravi_randn();
			  mI += 2 * eI * gravi_randn();
			}
            if (PHASEREF) {
                /* Integrate <R> and <I> rephased */
                tR[w] += cos(PHASEREF) * mR - sin(PHASEREF) * mI;
                tI[w] += cos(PHASEREF) * mI + sin(PHASEREF) * mR;
            } else {
                /* Integrate directly without rephasing */
                tR[w] += mR;
                tI[w] += mI;
            }
    
    		/* Compute the flux <F1F2> and <sqrt(|F1F2|)> x sign(F1F2)
    		 * corrected by the vFactor if needed */
            tF1F2[w] += FNORM * VFACTOR * PFACTOR;
            tF12[w]  += sqrt( CPL_MAX (FNORM * VFACTOR * PFACTOR, 0.0) );
			
    		/* Integrate < R2 + I2 - sR2 - sI2 > */
    		if (use_debiasing) {
    		  tP[w] += mR*mR + mI*mI - eR*eR - eI*eI;
    		} else {
    		  tP[w] += mR*mR + mI*mI;
    		}
    	  } /* End loop on wave */
    	} /* End loop on rows in this segment */
    }/* End loop on segments */

    
    /*
     * (2) Compute the variance by bootstraping on the segments
     */
    
    /* Loop on bootstramp to compute the avg and the 
     * variance by the bootstraping methode */
    srand(1);
    
    for ( int boot = 0 ; boot < nboot ; boot ++ ) {
    	cpl_msg_debug (cpl_func,"Bootstrap %d over %d", boot+1, nboot);
    
    	/* Init the itegration of nseg segments */
    	cpl_array * F12_boot   = gravi_array_init_double (nwave, 0.0);
    	cpl_array * F1F2_boot  = gravi_array_init_double (nwave, 0.0);
    			  
    	/* Integrate nseg segments randomly selected.
    	 * This loop is vectorialized in spectral direction */
    	for (int rowb = 0; rowb < nseg; rowb ++) {
    			  	
    	  /* For the first bootstrap, we use all observed samples
		   * For the others, we also includes the possible montecarlo
    	   * FIXME: verify the uniformity of rand for small nrows */
		  int rows;
    	  if (boot == 0 )  rows = rowb;
    	  else             rows = rand()%(nseg+nmontecarlo);
    
    	  /* Integrate the selected segments */
    	  cpl_array_add (visR_res[0], visR[rows]);
    	  cpl_array_add (visI_res[0], visI[rows]);
    	  cpl_array_add (vis2_res[0], POWER[rows]);
    	  cpl_array_add (F12_boot, F12[rows]);
    	  cpl_array_add (F1F2_boot, F1F2[rows]);
    	}
    	/* End loop to integrate nseg segments randomly selected 
    	 * We now have a random realisation of a high SNR dataset 
    	 * to which we can apply the estimators */
  
        /* Make sure the geometric flux is not null */
        gravi_array_threshold_min (F1F2_boot, 1e-15);
    			  
    	/* Energy estimator
    	   vis2 = POWER /  F1F2 */
    	cpl_array_divide (vis2_res[0], F1F2_boot);
    	CPLCHECK_MSG("while computing the energie integration");
    			  
    	/* Norm of coherent integration
    	   visAmp = sqrt(visR^2 + visR^2) / F12 */
	FREE (cpl_array_delete, visAmp_res[0]);
    	visAmp_res[0] = gravi_array_compute_norm2 (visR_res[0], visI_res[0]);
    	cpl_array_power (visAmp_res[0], 0.5);
    	cpl_array_divide (visAmp_res[0], F12_boot);
    	CPLCHECK_MSG("while computing the norm of the coherent integration");
    	
    	/* Phase of coherent integration 
    	   visPhi = arctan( visI, visR ) */
		FREE (cpl_array_delete, visPhi_res[0]);
    	visPhi_res[0] = gravi_array_wrap_complex (visR_res[0], visI_res[0]);
		cpl_array_arg (visPhi_res[0]);
    			
    	/* Compute the VARIANCE over the bootstraped samples with the 'online_variance' algorithm.
    	 * See https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance */
     	gravi_array_online_variance_res (vis2_res, boot, 0);
    	gravi_array_online_variance_res (visAmp_res, boot, 0);
    	gravi_array_online_variance_res (visR_res, boot, 0);
    	gravi_array_online_variance_res (visI_res, boot, 0);
    	gravi_array_online_variance_res (visPhi_res, boot, 1);
    	CPLCHECK_MSG("while computing the variances over the bootstrap");

    	/* Delete these temporary integrations */
    	FREE (cpl_array_delete, F12_boot);
    	FREE (cpl_array_delete, F1F2_boot);
    	CPLCHECK_MSG("while freeing during bootstrap");
    }
    /* End bootstrap */
    
    /* Free list of segments */
    FREELOOP (cpl_array_delete, visR, nseg + nmontecarlo);
    FREELOOP (cpl_array_delete, visI, nseg + nmontecarlo);
    FREELOOP (cpl_array_delete, POWER, nseg + nmontecarlo);
    FREELOOP (cpl_array_delete, F12, nseg + nmontecarlo);
    FREELOOP (cpl_array_delete, F1F2, nseg + nmontecarlo);

	/* Convert variance from bootstrap to RMS */
	cpl_msg_debug (cpl_func,"Put the RMS over bootstrap");	
	cpl_array_power (vis2_res[3], 0.5);
	cpl_array_power (visAmp_res[3], 0.5);
	cpl_array_power (visPhi_res[3], 0.5);
	cpl_array_power (visR_res[3], 0.5);
	cpl_array_power (visI_res[3], 0.5);
	CPLCHECK_MSG ("while converting variance -> rms");

	/* Normalize integration */
	mjd_avg /= total_exptime;
	uCoord /= total_exptime;
	vCoord /= total_exptime;
	
  } 

  /* 
   * Step (1,2) in case there are no valid frames at all
   */
  if (nvalid == 0) {
	cpl_msg_debug (cpl_func,"Not valid frames, force zero and infinit RMS");
	cpl_array_fill_window (vis2_res[3], 0, nwave, 1e10);
	cpl_array_fill_window (visR_res[3], 0, nwave, 1e10);
	cpl_array_fill_window (visI_res[3], 0, nwave, 1e10);
	cpl_array_fill_window (visAmp_res[3], 0, nwave, 1e10);
	cpl_array_fill_window (visPhi_res[3], 0, nwave, 1e10);
	mjd_avg = cpl_table_get_column_mean (oi_vis, "MJD");
  }

  /* Renormalise by the number of outliers.
     Also flag channels with too much outliers. */
  cpl_array ** pflag_vis  = cpl_table_get_data_array (oi_vis_avg,  "FLAG");
  cpl_array ** pflag_vis2 = cpl_table_get_data_array (oi_vis2_avg, "FLAG");
  cpl_array * array = gravi_table_get_column_sum_array (oi_vis, "FLAG", base, nbase);
  CPLCHECK_MSG("cannot get data");
  
  for (int wave = 0; wave < nwave; wave++) {
    double value = cpl_array_get (array, wave, NULL) / nrow;
    if (value < 1.0) {
      cpl_array_set (visR_res[2], wave, cpl_array_get (visR_res[2], wave, NULL) / (1.0 - value));
      cpl_array_set (visI_res[2], wave, cpl_array_get (visI_res[2], wave, NULL) / (1.0 - value));
    }
    if (value > outlier_threshold) {
      cpl_array_set (pflag_vis[base],  wave, 1);
      cpl_array_set (pflag_vis2[base], wave, 1);
    }
  }
  FREE (cpl_array_delete, array);
    
  /* 
   * (3) Save the results on the oi_vis_avg tables 
   */
  
  cpl_table_set_array (oi_vis2_avg, "VIS2DATA", base, vis2_res[2]);
  cpl_table_set_array (oi_vis2_avg, "VIS2ERR", base, vis2_res[3]);
  CPLCHECK_MSG("filling VIS2");
  
  gravi_table_set_array_double_complex (oi_vis_avg, "VISDATA", base, visR_res[2], visI_res[2]);
  gravi_table_set_array_double_complex (oi_vis_avg, "VISERR", base, visR_res[3], visI_res[3]);
  CPLCHECK_MSG("filling VISDATA");

  cpl_table_set_array (oi_vis_avg, "RVIS", base, visR_res[2]);
  cpl_table_set_array (oi_vis_avg, "RVISERR", base, visR_res[3]);
  CPLCHECK_MSG("filling RVIS");

  cpl_table_set_array (oi_vis_avg, "IVIS", base, visI_res[2]);
  cpl_table_set_array (oi_vis_avg, "IVISERR", base, visI_res[3]);
  CPLCHECK_MSG("filling IVIS");
  
  cpl_table_set_array (oi_vis_avg, "VISAMP", base, visAmp_res[2]);
  cpl_table_set_array (oi_vis_avg, "VISAMPERR", base, visAmp_res[3]);
  CPLCHECK_MSG("filling VISAMP");
  
  gravi_table_set_array_phase (oi_vis_avg, "VISPHI", base, visPhi_res[2]);
  gravi_table_set_array_phase (oi_vis_avg, "VISPHIERR", base, visPhi_res[3]);
  CPLCHECK_MSG("filling VISPHI");

  /* Flag the data with >100% error or >60deg error */
  gravi_vis_flag_threshold (oi_vis2_avg, "VIS2ERR", "FLAG", 1.0);
  gravi_vis_flag_lower (oi_vis2_avg, "VIS2ERR", "FLAG", 0.0);
  gravi_vis_flag_median (oi_vis2_avg, "VIS2ERR", "FLAG", 5.0);
  gravi_vis_flag_threshold (oi_vis_avg,  "VISAMPERR", "FLAG", 1.0);
  gravi_vis_flag_lower (oi_vis_avg,  "VISAMPERR", "FLAG", 0.0);
  gravi_vis_flag_median (oi_vis_avg, "VISPHIERR", "FLAG", 5.0);
  CPLCHECK_MSG("cannot flag baddata data");
  
  /* Compute the total integration time */
  cpl_table_set (oi_vis_avg,  "INT_TIME", base, total_exptime);
  cpl_table_set (oi_vis2_avg, "INT_TIME", base, total_exptime);
  cpl_table_set (oi_vis_avg,  "MJD", base, mjd_avg);
  cpl_table_set (oi_vis2_avg, "MJD", base, mjd_avg);
  cpl_table_set (oi_vis_avg,  "UCOORD", base, uCoord);
  cpl_table_set (oi_vis_avg,  "VCOORD", base, vCoord);
  cpl_table_set (oi_vis2_avg, "UCOORD", base, uCoord);
  cpl_table_set (oi_vis2_avg, "VCOORD", base, vCoord);
  CPLCHECK_MSG("cannot fill time");

  /* Set some statistics */
  cpl_table_set (oi_vis2_avg, "NVALID", base, nvalid);
  cpl_table_set (oi_vis2_avg, "NDIT", base, nrow);
  cpl_table_set (oi_vis_avg, "NVALID", base, nvalid);
  cpl_table_set (oi_vis_avg, "NDIT", base, nrow);
  CPLCHECK_MSG("cannot fill nvalid");
  
  /* Set the TARGET_ID and STA_INDEX */
  cpl_table_set_int (oi_vis_avg,  "TARGET_ID", base, cpl_table_get_int (oi_vis, "TARGET_ID", base, &nv));
  cpl_table_set_int (oi_vis2_avg, "TARGET_ID", base, cpl_table_get_int (oi_vis, "TARGET_ID", base, &nv));
  cpl_table_set_array (oi_vis_avg,  "STA_INDEX", base, cpl_table_get_array (oi_vis, "STA_INDEX", base));
  cpl_table_set_array (oi_vis2_avg, "STA_INDEX", base, cpl_table_get_array (oi_vis, "STA_INDEX", base));

  /* Free variance */
  FREELOOP (cpl_array_delete, vis2_res, 4);
  FREELOOP (cpl_array_delete, visAmp_res, 4);
  FREELOOP (cpl_array_delete, visPhi_res, 4);
  FREELOOP (cpl_array_delete, visR_res, 4);
  FREELOOP (cpl_array_delete, visI_res, 4);

  }
  /* End loop on bases */
  
  gravi_msg_function_exit(0);
  return CPL_ERROR_NONE;
}
 
double gdAbacusErrPhi(double x)
/**
 * Estimate true phase rms from the cross-spectrum variance.
 * see Petrov, Roddier and Aime, JOSAA vol 3, N.5, may 1986 p 634.
 * and Petrov's Thesis, p. 50 ff.
 * I replace the piecewise interpolation usually used by a polynomial
 * approximation of the function:
 * if z=log10(y), 
 * then z=(C[1]*X^7+C[2]*X^6+C[3]*X^5+C[4]*X^4+C[5]*X^3+C[6]*X^2+C[7]*X+C[8])
 * and y=10^z.
 * where
 * C[01]= 2.71918080109099
 * C[02]=-17.1901043936273
 * C[03]= 45.0654103760899
 * C[04]=-63.4441678243197
 * C[05]= 52.3098941426378
 * C[06]=-25.8090699917488
 * C[07]= 7.84352873962491
 * C[08]=-1.57308595820081
 * This interpolation is valid in the range x=[0.1,1.74413]. 
 * Error is 1% everywhere except above x=1.73 where it is 10%
 * Below x=0.1: y=x
 * Above x=M_PI/sqrt(3.0), y = blanking (impossible value)
 * Above x=1.74413, we take: y=0.691/(pi/sqrt(3.0)-x)
 */
{
  const double Asymptot = CPL_MATH_PI / sqrt(3.0);
  double c[8] = {2.7191808010909,
    -17.1901043936273,
    45.0654103760899,
    -63.4441678243197,
    52.3098941426378,
    -25.8090699917488,
    7.84352873962491,
    -1.57308595820081};

  double x2, x3, x4, x5, x6, x7, z;
  if (x > Asymptot) {
    return (1E10);
  }
  if (x > 1.74413) {
    return (0.691 / (Asymptot - x));
  }
  if (x < 0.1) {
    return (x);
  }
  x2 = x*x;
  x3 = x2*x;
  x4 = x2*x2;
  x5 = x3*x2;
  x6 = x3*x3;
  x7 = x6*x;
  z = c[0] * x7 + c[1] * x6 + c[2] * x5 + c[3] * x4 + c[4] * x3 + c[5] * x2 + c[6] * x + c[7];
  return pow(10, z);
}

/*-----------------------------------------------------------------------------*/
/**
 * @brief Compute Averaged VISPHI in the manner described, e.g., in F. Millour's thesis.
 * 
 * @param oi_vis_avg:     already allocated OI_VIS table to be filled
 * 
 * @param oi_vis:         input OI_VIS with all individual DITs
 * @param base:           base to consider (0..5)
 * @param phase_ref:      phase used to rephase the DITs
 *
 */
/*-----------------------------------------------------------------------------*/

cpl_error_code gravi_average_self_visphi(cpl_table * oi_vis_avg,
    cpl_table * oi_vis,
    cpl_array * wavenumber,
    const char * phase_ref, int* cmin, int* cmax, int nrange)
{
  gravi_msg_function_start(0);
  cpl_ensure_code(oi_vis_avg, CPL_ERROR_ILLEGAL_OUTPUT);
  cpl_ensure_code(oi_vis, CPL_ERROR_NULL_INPUT);

  /* Parameters */
  int nbase = 6;
  cpl_size nrow = cpl_table_get_nrow(oi_vis) / nbase;
  cpl_size nwave = cpl_table_get_column_depth(oi_vis, "VISDATA");

  int use_crange = (nrange > 0);
  if (use_crange) {
    for (int i = 0; i < nrange; ++i) {
      if (cmax[i] - cmin[i] < 1) {
        use_crange = 0;
        break;
      }
      if (cmax[i] > nwave - 1) {
        use_crange = 0;
        break;
      }
      if (cmin[i] < 0) {
        use_crange = 0;
        break;
      }
    }
    if (use_crange) {
      for (int i = 0; i < nrange; ++i) cpl_msg_info("Reference Channel", "Part %02d [%3d:%3d]", i + 1, cmin[i], cmax[i]);
    } else {
      cpl_msg_info("Warning (SELF_VISPHI)", "Invalid Ranges Found, continuing with default Method.");
    }
  }

  /* Pointer to columns, to speed-up */
  cpl_array ** pVISDATA = cpl_table_get_data_array(oi_vis, "VISDATA");
  cpl_array ** pVISERR = cpl_table_get_data_array(oi_vis, "VISERR");
  CPLCHECK_MSG("Cannot get the data");

  /* Get the reference phase unless it has already been removed*/
  cpl_array ** pPHASEREF = NULL;
  if (phase_ref && strcmp(phase_ref, "NONE"))
    pPHASEREF = cpl_table_get_data_array(oi_vis, phase_ref);
  CPLCHECK_MSG("Cannot get the reference phase data (column missing?)");

  /* Loop on base */
  for (cpl_size base = 0; base < nbase; base++) {

    /* number of valid rows for base */
    cpl_size nvalid = 0;

    /* Arrays to store the final, integrated quantities 
     * 0: running_mean, 2: variance */
    cpl_array **visPhi_res = gravi_array_new_list(2, CPL_TYPE_DOUBLE, nwave);

    /* Get the number of non-rejected frames */
    int * flag = cpl_table_get_data_int(oi_vis, "REJECTION_FLAG");
    cpl_ensure_code(flag, CPL_ERROR_ILLEGAL_INPUT);
    for (int row = 0; row < nrow; row++) if (flag[row * nbase + base] == 0) nvalid++;

    cpl_msg_info("Stat (SELF_VISPHI)", "%6lld valid frames over %6lld (%5.1f%%)",
        nvalid, nrow, (double) nvalid / (double) nrow * 100.0);

    /* 
     * return in case there are no valid frames at all
     */
    if (nvalid == 0) {
      cpl_msg_debug(cpl_func, "No valid frames, force zero and infinite RMS");
      cpl_array_fill_window(visPhi_res[1], 0, nwave, 1e10);
      gravi_table_set_array_phase(oi_vis_avg, "VISPHI", base, visPhi_res[0]);
      gravi_table_set_array_phase(oi_vis_avg, "VISPHIERR", base, visPhi_res[1]);
      CPLCHECK_MSG("filling VISPHI");
      /* Free variance */
      FREELOOP(cpl_array_delete, visPhi_res, 2);
      break;
    }

    cpl_array **Vis = gravi_array_new_list(nvalid, CPL_TYPE_DOUBLE_COMPLEX, nwave);
    cpl_array **EVis = gravi_array_new_list(nvalid, CPL_TYPE_DOUBLE_COMPLEX, nwave);

    cpl_array **W1 = gravi_array_new_list(nvalid, CPL_TYPE_DOUBLE_COMPLEX, nwave);
    cpl_array **EW1 = gravi_array_new_list(nvalid, CPL_TYPE_DOUBLE_COMPLEX, nwave);
    /* Loop on row */

    for (cpl_size currentRow = 0, validRowIndex = -1; currentRow < nrow; currentRow++) {
      if (flag[currentRow * nbase + base]) {
        continue;
      } else {
        /* Get indices */
        cpl_size rbase = currentRow * nbase + base;
        validRowIndex++;

        /* fast-no-CPL integration: get pointers on data */
        double complex *ptrC = cpl_array_get_data_double_complex(Vis[validRowIndex]);
        double complex *ptrEC = cpl_array_get_data_double_complex(EVis[validRowIndex]);
        CPLCHECK_MSG("Cannot get data");

        /* Loop on wave */
        for (int w = 0; w < nwave; w++) {
          double PHASEREF = pPHASEREF ? cpl_array_get(pPHASEREF[rbase], w, NULL) : 0.0;
          double complex vis = cpl_array_get_double_complex(pVISDATA[rbase], w, NULL);
          double complex viserr = cpl_array_get_double_complex(pVISERR[rbase], w, NULL);
          CPLCHECK_MSG("Cannot get data");
          if (PHASEREF) {
            /* rephase <R> and <I> */
            ptrC[w] = (cos(PHASEREF) * creal(vis) - sin(PHASEREF) * cimag(vis)) +
                I * (cos(PHASEREF) * cimag(vis) + sin(PHASEREF) * creal(vis));
          } else {
            /* no rephasing */
            ptrC[w] = vis;
          }
          ptrEC[w] = viserr;
        } /* End loop on wave */
      } /* End frame not flagged */
    } /* End loop on rows */

    for (int irow = 0; irow < nvalid; irow++) {

      /* Normalize the phasor to avoid bad pixels */
      /*
            cpl_array * norm = cpl_array_duplicate(Vis[irow]);
            cpl_array_abs(norm);
            cpl_array_divide(Vis[irow], norm);
       */

      /* Compute and remove the mean group delay in [m] */
      double mean_delay = 0.0;
      gravi_array_get_group_delay_loop(&Vis[irow], NULL, wavenumber, &mean_delay, 1, 2e-3, CPL_FALSE);
      gravi_array_multiply_phasor(Vis[irow], -2 * I * CPL_MATH_PI * mean_delay, wavenumber);

      /* Compute and remove the mean phase in [rad] */
      double mean_phase = carg(cpl_array_get_mean_complex(Vis[irow]));
      cpl_array_multiply_scalar_complex(Vis[irow], cexp(-I * mean_phase));
    }

    cpl_array *CRef = cpl_array_new(nwave, CPL_TYPE_DOUBLE_COMPLEX);
    cpl_array *ECRef = cpl_array_new(nwave, CPL_TYPE_DOUBLE_COMPLEX);
    double complex *pCRef = cpl_array_get_data_double_complex(CRef);
    double complex *pECRef = cpl_array_get_data_double_complex(ECRef);

    /* Production of the reference channel: mean of all R and I parts
     * of all channels except the one considered ( eq 2.3) */
    for (int irow = 0; irow < nvalid; irow++) {
      /*reset sum to 0 */
      double complex totalVis = 0.0 + I * 0.0;
      double complex totalEVis = 0.0 + I * 0.0;
      /* fast-no-CPL integration: get pointers on data */
      double complex *pVis = cpl_array_get_data_double_complex(Vis[irow]);
      double complex *pEVis = cpl_array_get_data_double_complex(EVis[irow]);

      double complex *pW1 = cpl_array_get_data_double_complex(W1[irow]);
      double complex *pEW1 = cpl_array_get_data_double_complex(EW1[irow]);

      CPLCHECK_MSG("Cannot get data");
      if (use_crange) {
        /* sum all Vis for this row */
        int nchans = 0;
        for (int i = 0; i < nrange; ++i) {
          for (int w = cmin[i]; w < cmax[i]; ++w) {
            totalVis += pVis[w];
            totalEVis += pEVis[w];
            nchans++;
          }
        }
        for (int w = 0; w < nwave; w++) {
          cpl_array_set_double_complex(CRef, w, totalVis / nchans);
          cpl_array_set_double_complex(ECRef, w, totalEVis / nchans);
        }
      } else {
        /* sum all Vis for this row */
        for (int w = 0; w < nwave; w++) {
          totalVis += pVis[w];
          totalEVis += pEVis[w];
        }
        /* then construct Cref by substracting current R and I 
         * at that Wlen and make the arithmetic mean. The code permits
         * to avoid not only the channel itself removed but the 2*radius
         * around (not activated). */
        int iw = 0;
        int radius = 0;
        int divider = nwave - (2 * radius) - 1;
        for (; iw < radius; iw++) {
          cpl_array_set_double_complex(CRef, iw, totalVis / nwave);
          cpl_array_set_double_complex(ECRef, iw, totalEVis / nwave);
        }
        for (; iw < nwave - radius; iw++) {
          double complex tmp = 0.0 + I * 0.0;
          double complex Etmp = 0.0 + I * 0.0;
          for (int j = iw; j < iw + 2 * radius + 1; ++j) tmp += pVis[iw];
          cpl_array_set_double_complex(CRef, iw, (totalVis - tmp) / divider);
          for (int j = iw; j < iw + 2 * radius + 1; ++j) Etmp += pEVis[iw];
          cpl_array_set_double_complex(ECRef, iw, (totalEVis - Etmp) / divider);
        }
        for (; iw < radius; iw++) {
          cpl_array_set_double_complex(CRef, iw, totalVis / nwave);
          cpl_array_set_double_complex(ECRef, iw, totalEVis / nwave);
        }
      }
      /* Now the interspectrum is C*~C_Ref. Store in w1. */
      for (int w = 0; w < nwave; w++) {
        pW1[w] = pVis[w] * conj(pCRef[w]);

        /* Please have a look to the F. Millour thesis
           (http://tel.archives-ouvertes.fr/tel-00134268),
           pp.91-92 (eq. 4.55 to 4.58) */
        pEW1[w] = (creal(pEVis[w]) * pow(creal(pCRef[w]), 2) +
            creal(pECRef[w]) * pow(creal(pVis[w]), 2) +
            cimag(pVis[w]) * pow(cimag(pCRef[w]), 2) +
            cimag(pECRef[w]) + pow(cimag(pVis[w]), 2)) +

            I * (
            cimag(pVis[w]) * pow(creal(pCRef[w]), 2) +
            cimag(pECRef[w]) * pow(creal(pVis[w]), 2) +
            creal(pVis[w]) * pow(cimag(pCRef[w]), 2) +
            creal(pECRef[w]) + pow(cimag(pVis[w]), 2)
            );
      }
    }

    FREE(cpl_array_delete, CRef);
    FREE(cpl_array_delete, ECRef);
    FREELOOP(cpl_array_delete, Vis, nvalid);
    FREELOOP(cpl_array_delete, EVis, nvalid);

    double *pPhi = cpl_array_get_data_double(visPhi_res[0]);
    double *pPhiErr = cpl_array_get_data_double(visPhi_res[1]);

    /* Compute mean VisPhi as average of selected frames. */
    for (int w = 0; w < nwave; w++) {
      cpl_array *cpxVisVect = cpl_array_new(nvalid, CPL_TYPE_DOUBLE_COMPLEX);

      /* The W1 vector */
      for (int irow = 0; irow < nvalid; irow++) {
        const double complex *pW1 = cpl_array_get_data_double_complex_const(W1[irow]);
        cpl_array_set_double_complex(cpxVisVect, irow, pW1[w]);
      }
      /* The Phase Herself */
      /* average re and im */
      double complex w1Avg;
      w1Avg = cpl_array_get_mean_complex(cpxVisVect);

      /* store */
      pPhi[w] = atan2(cimag(w1Avg), creal(w1Avg));
      /* WE USE THE STATISTICAL ERROR FOR BINNING */
      w1Avg = conj(w1Avg);
      cpl_array *Vect = cpl_array_new(nvalid, CPL_TYPE_DOUBLE);
      for (int irow = 0; irow < nvalid; irow++) {
        const double complex *tW1 = cpl_array_get_data_double_complex_const(W1[irow]);
        /* add w1*conj(w1Avg) to vector*/
        cpl_array_set_double(Vect, irow, atan2(cimag(tW1[w] * w1Avg), creal(tW1[w] * w1Avg)));
      }
      double x = cpl_array_get_stdev(Vect);
      /* Err on Phi must be corrected with an abacus*/
      pPhiErr[w] = gdAbacusErrPhi(x / sqrt(nvalid));
      /* START EKW 21/11/2018 */
      FREE(cpl_array_delete,  cpxVisVect);
      FREE(cpl_array_delete,  Vect);
      /* END EKW 21/11/2018 */
    }

    gravi_table_set_array_phase(oi_vis_avg, "VISPHI", base, visPhi_res[0]);
    gravi_table_set_array_phase(oi_vis_avg, "VISPHIERR", base, visPhi_res[1]);
    CPLCHECK_MSG("filling VISPHI");
    /* Free variance */
    FREELOOP(cpl_array_delete, visPhi_res, 2);
    /* START EKW 21/11/2018 */
    FREELOOP(cpl_array_delete, EW1 , nvalid);
    FREELOOP(cpl_array_delete, W1  , nvalid);
    /* END EKW 21/11/2018 */

  }

  /* End loop on bases */

  gravi_msg_function_exit(0);
  return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief The function average the individual frames of a P2VMREDUCED file
 *        into a final, single observation per base and per tel.
 * 
 * @param p2vmred_data  P2VMREDUCED data containing OI_FLUX and OI_VIS tables
 * 	  	                coming from the @c gravi_compute_p2vmred
 * @param parlist       Parameters of the recipe
 */
/*----------------------------------------------------------------------------*/

gravi_data * gravi_compute_vis (gravi_data * p2vmred_data,
                                const cpl_parameterlist * parlist,
                                cpl_size * current_frame)
{
	gravi_msg_function_start(1);
	cpl_ensure (p2vmred_data, CPL_ERROR_NULL_INPUT, NULL);
	cpl_ensure (parlist,      CPL_ERROR_NULL_INPUT, NULL);
	
	int nv, nbase = 6, ntel = 4;

    /* 
     * Compute the limit of integration
     */

    /* Get the current position and the maximum nb to integrate */
    cpl_size max_frame = gravi_param_get_int (parlist, "gravity.vis.max-frame");
    cpl_msg_info (cpl_func,"Average %lli frames starting from %lli",
                  max_frame, *current_frame);
    
    

    /* Get the SC DIT */
	cpl_propertylist * p2vmred_header = gravi_data_get_header (p2vmred_data);
    double dit_sc = gravi_pfits_get_dit_sc (p2vmred_header) * 1e6;
    
    /* Get the vis_SC table for first polarisation */
    int npol_sc = gravi_pfits_get_pola_num (p2vmred_header, GRAVI_SC);
    cpl_table * vis_SC = gravi_data_get_oi_vis (p2vmred_data, GRAVI_SC, 0, npol_sc);
    cpl_size nrow = cpl_table_get_nrow (vis_SC) / nbase;
    
    /* Get first and last frame for this integration */
    cpl_size sframe = *current_frame;
    cpl_size eframe = CPL_MIN (*current_frame + max_frame - 1, nrow-1);
    
    /* Check if we reached the end of the file, or increment */
    if (eframe >= nrow-1)  *current_frame = -1;
    else                   *current_frame += max_frame;
    
    /* Compute start and end-time */
    double start_time, end_time;
    start_time  = cpl_table_get (vis_SC, "TIME", sframe*nbase, &nv) - dit_sc/2;
    end_time    = cpl_table_get (vis_SC, "TIME", eframe*nbase, &nv) + dit_sc/2;
    
    /* Compute verbose */
    cpl_msg_info (cpl_func,"Integrate frames: first = %lli  last = %lli", sframe, eframe);
    cpl_msg_info (cpl_func,"start = %f  end = %f [s]", start_time*1e-6, end_time*1e-6);
    

	/* 
	 * Prepare the output 
	 */
	
    cpl_msg_info(cpl_func, "Construction of the averaged output data");

    gravi_data * vis_data = gravi_data_new (0);
    cpl_propertylist * vis_header = gravi_data_get_header (vis_data);
    cpl_propertylist_append (vis_header, p2vmred_header);
		
    /* Copy the oi tables needed in output data.
     * This will duplicate all OI_WAVELENGTH tables */
    gravi_data_copy_ext (vis_data, p2vmred_data, GRAVI_OI_ARRAY_EXT);
    gravi_data_copy_ext (vis_data, p2vmred_data, GRAVI_OI_TARGET_EXT);
    gravi_data_copy_ext (vis_data, p2vmred_data, GRAVI_OI_WAVELENGTH_EXT);
    
    CPLCHECK_NUL ("Cannot get tables for output data");
    
    /* 
     * Start with FT 
     */
    
    if (gravi_data_has_type (p2vmred_data, "_FT") <= 0 ) {
		cpl_msg_info (cpl_func, "P2VMRED data has no FT extensions");
    }
    else {
		/* Reduction parameters */
		int v_factor_flag_ft = 0;
		int p_factor_flag_ft = 0;
		int debiasing_flag_ft = gravi_param_get_bool (parlist, "gravity.vis.debias-ft");
        int nboot_ft = gravi_param_get_int (parlist, "gravity.vis.nboot");
        const char * phase_ref_ft = "SELF_REF";
	double outlier_threshold_ft = 0.0;
		CPLCHECK_NUL("Cannot get parameters");

        cpl_msg_info (cpl_func, "Bias subtraction of V2 for FT is %s",debiasing_flag_ft?"ENABLE":"DISABLE");
		cpl_msg_info (cpl_func, "Reference phase for FT is %s",phase_ref_ft);


        /*
         * Loop on polarisations 
         */
        int npol_ft = gravi_pfits_get_pola_num (p2vmred_header, GRAVI_FT);
        for (int pol = 0; pol < npol_ft; pol++) {
            cpl_msg_info (cpl_func, "Start FT polarisation %d over %d",pol+1, npol_ft);
            
            /* Get the input table of FT */
            cpl_table * vis_FT = gravi_data_get_oi_vis (p2vmred_data, GRAVI_FT, pol, npol_ft);
            cpl_table * flux_FT = gravi_data_get_oi_flux (p2vmred_data, GRAVI_FT, pol, npol_ft);
            int nwave_ft = cpl_table_get_column_depth (vis_FT, "VISDATA");
            CPLCHECK_NUL ("Cannot get data");

            /* Create averated product */
            cpl_table * oi_vis2_FT = gravi_table_oi_create (nwave_ft, 1, GRAVI_OI_VIS2_EXT);
            gravi_table_new_column (oi_vis2_FT, "NDIT", NULL, CPL_TYPE_INT);
            gravi_table_new_column (oi_vis2_FT, "NVALID", NULL, CPL_TYPE_INT);
            
            cpl_table * oi_vis_FT = gravi_table_oi_create (nwave_ft, 1, GRAVI_OI_VIS_EXT);
            gravi_table_new_column (oi_vis_FT, "NDIT", NULL, CPL_TYPE_INT);
            gravi_table_new_column (oi_vis_FT, "NVALID", NULL, CPL_TYPE_INT);

            cpl_table * oi_T3_FT = gravi_table_oi_create (nwave_ft, 1, GRAVI_OI_T3_EXT);
            gravi_table_new_column (oi_T3_FT, "NDIT", NULL, CPL_TYPE_INT);
            gravi_table_new_column (oi_T3_FT, "NVALID", NULL, CPL_TYPE_INT);

            cpl_table * oi_flux_FT = gravi_table_oi_create (nwave_ft, 1, GRAVI_OI_FLUX_EXT);
            gravi_table_new_column (oi_flux_FT, "NDIT", NULL, CPL_TYPE_INT);
            gravi_table_new_column (oi_flux_FT, "NVALID", NULL, CPL_TYPE_INT);
            CPLCHECK_NUL ("Cannot create product");

            /* Keep only selected rows */
            vis_FT  = gravi_table_extract_time_interval (vis_FT, start_time, end_time);
            flux_FT = gravi_table_extract_time_interval (flux_FT, start_time, end_time);

            /* 
             * Compute OIVIS2 and OIVIS for FT
             */
            cpl_msg_info (cpl_func, "Compute OIVIS2 and OIVIS for FT");
            
            gravi_vis_average_bootstrap (oi_vis_FT, oi_vis2_FT, vis_FT,
                                         nboot_ft,
                                         phase_ref_ft,
                                         v_factor_flag_ft,
                                         p_factor_flag_ft,
                                         debiasing_flag_ft,
					 outlier_threshold_ft);
            CPLCHECK_NUL("Cannot average the FT frames");

            /* 
             * Compute OIT3 for FT
             */
            cpl_msg_info (cpl_func, "Compute OIT3 for FT");

            gravi_t3_average_bootstrap (oi_T3_FT, vis_FT, flux_FT,
                                        nboot_ft,
                                        v_factor_flag_ft,
                                        p_factor_flag_ft,
					outlier_threshold_ft);
            CPLCHECK_NUL("Cannot average t3 of FT");
            
            /* 
             * Compute OI_FLUX for FT
             */
            cpl_msg_info (cpl_func, "Compute OI_FLUX for FT");
		
            gravi_flux_average_bootstrap (oi_flux_FT, flux_FT,
                                          nboot_ft,
					  outlier_threshold_ft);
            CPLCHECK_NUL("Cannot average flux of FT");

            /* 
             * Add tables in the vis_data
             */
            cpl_propertylist * vis_plist   = gravi_data_get_oi_vis_plist (p2vmred_data, GRAVI_FT, pol, npol_ft);
            cpl_propertylist * oivis_plist = cpl_propertylist_new();
            cpl_propertylist_copy_property (oivis_plist, vis_plist, "DATE-OBS");
            cpl_propertylist_copy_property (oivis_plist, vis_plist, "EXTVER");
            cpl_propertylist_copy_property (oivis_plist, vis_plist, "ARRNAME");
            cpl_propertylist_copy_property (oivis_plist, vis_plist, "INSNAME");
            cpl_propertylist_update_string (oivis_plist, "AMPTYP","absolute");
            cpl_propertylist_update_string (oivis_plist, "PHITYP","differential");
            cpl_propertylist_update_int (oivis_plist, "PHIORDER",1);
            gravi_data_add_table (vis_data, oivis_plist, GRAVI_OI_VIS_EXT, oi_vis_FT);
            
            cpl_propertylist * oivis2_plist = cpl_propertylist_new();
            cpl_propertylist_copy_property (oivis2_plist, vis_plist, "DATE-OBS");
            cpl_propertylist_copy_property (oivis2_plist, vis_plist, "EXTVER");
            cpl_propertylist_copy_property (oivis2_plist, vis_plist, "ARRNAME");
            cpl_propertylist_copy_property (oivis2_plist, vis_plist, "INSNAME");
            gravi_data_add_table (vis_data, oivis2_plist, GRAVI_OI_VIS2_EXT, oi_vis2_FT);
            
            cpl_propertylist * oit3_plist = cpl_propertylist_new();
            cpl_propertylist_copy_property (oit3_plist, vis_plist, "DATE-OBS");
            cpl_propertylist_copy_property (oit3_plist, vis_plist, "EXTVER");
            cpl_propertylist_copy_property (oit3_plist, vis_plist, "ARRNAME");
            cpl_propertylist_copy_property (oit3_plist, vis_plist, "INSNAME");
            gravi_data_add_table (vis_data, oit3_plist, GRAVI_OI_T3_EXT, oi_T3_FT);
            
            cpl_propertylist * oiflux_plist = cpl_propertylist_new();
            cpl_propertylist_copy_property (oiflux_plist, vis_plist, "DATE-OBS");
            cpl_propertylist_copy_property (oiflux_plist, vis_plist, "EXTVER");
            cpl_propertylist_copy_property (oiflux_plist, vis_plist, "ARRNAME");
            cpl_propertylist_copy_property (oiflux_plist, vis_plist, "INSNAME");
            cpl_propertylist_update_string (oiflux_plist, "CALSTAT", "U");
            cpl_propertylist_set_comment (oiflux_plist, "CALSTAT", "Uncalibrated flux per telescope");
            gravi_data_add_table (vis_data, oiflux_plist, GRAVI_OI_FLUX_EXT, oi_flux_FT);
            
            CPLCHECK_NUL ("Cannot add tables");

            FREE (cpl_table_delete, vis_FT);
            FREE (cpl_table_delete, flux_FT);
            
        } /* end loop on pol */
    } /* End FT */

    /* 
     * Then with SC 
     */
    
    if (gravi_data_has_type (p2vmred_data, "_SC") <= 0 ) {
		cpl_msg_info (cpl_func, "P2VMRED data has no SC extensions");
    }
    else {

        /* Read the object separation. Actually should
         * be the mode maybe more than the separation */
        double SOBJ_X = gravi_pfits_get_sobj_x (p2vmred_header);
        double SOBJ_Y = gravi_pfits_get_sobj_y (p2vmred_header);
        double SOBJ_R = sqrt(SOBJ_X*SOBJ_X + SOBJ_Y*SOBJ_Y);

        /* Reference phase requested by user, deal with AUTO */
        const char * phase_ref_sc = gravi_param_get_string (parlist, "gravity.vis.phase-ref-sc");

        if ( !strcmp (phase_ref_sc,"AUTO") && SOBJ_R > 0.001)
        {
            phase_ref_sc = "IMAGING_REF";
            if (cpl_table_has_column (vis_SC, phase_ref_sc) == 0) {
                phase_ref_sc = "PHASE_REF";
                cpl_msg_warning (cpl_func, "No table 'IMAGING_REF', changing mode to phase_ref_sc='PHASE_REF'");
                CPLCHECK_NUL("tototo....");
            }
        }
        else if ( !strcmp (phase_ref_sc,"AUTO") && SOBJ_R < 0.001) 
            phase_ref_sc = "PHASE_REF";
        
        cpl_msg_info (cpl_func, "Reference phase for SC is %s",phase_ref_sc);

        /* Output phase requested by user, deal with AUTO */
        const char * output_phase_sc = gravi_param_get_string (parlist, "gravity.vis.output-phase-sc");
        
        if ( !strcmp (output_phase_sc,"AUTO") && SOBJ_R > 0.001)
            output_phase_sc = "ABSOLUTE";
        else if ( !strcmp (output_phase_sc,"AUTO") && SOBJ_R < 0.001) 
            output_phase_sc = "DIFFERENTIAL";
        
        /*
          force SELF_REF if SELF_VISPHI unless NONE has been selected (for what purpose?)
        */
        if (!strcmp(output_phase_sc, "SELF_VISPHI")) {
          if (strcmp(phase_ref_sc, "NONE")) {
            phase_ref_sc = "SELF_REF"; /*we are robust to phase_ref_sc=NONE but SELF_REF is better */
            cpl_msg_info(cpl_func, "Reference phase for SC forced to %s due to option SELF_VISPHI", phase_ref_sc);
           } else cpl_msg_info(cpl_func, "Reference phase for SC is %s", phase_ref_sc);
        } else cpl_msg_info(cpl_func, "Reference phase for SC is %s", phase_ref_sc);
        
        cpl_msg_info (cpl_func, "Output phase for SC is %s",output_phase_sc);
        
		/* Other reduction parameters */
		int v_factor_flag_sc = strstr (gravi_param_get_string (parlist, "gravity.vis.vis-correction-sc"),"VFACTOR") ? 1 : 0;
		int p_factor_flag_sc = strstr (gravi_param_get_string (parlist, "gravity.vis.vis-correction-sc"),"PFACTOR") ? 1 : 0;
		int debiasing_flag_sc = gravi_param_get_bool (parlist, "gravity.vis.debias-sc");
		int nboot_sc = gravi_param_get_int (parlist, "gravity.vis.nboot");
        const char* rangeString = gravi_param_get_string (parlist, "gravity.vis.output-phase-channels");
	double outlier_threshold_sc = gravi_param_get_double (parlist, "gravity.vis.outlier-fraction-threshold");
	cpl_msg_info (cpl_func, "Bias subtraction of V2 for SC is %s",debiasing_flag_sc?"ENABLE":"DISABLE");
	cpl_msg_info (cpl_func, "Threshold for fraction of outlier is %.3f",outlier_threshold_sc);
        /*
          force VFACTOR and PFACTOR to 0 for SELF_VISPHI by precaution.
        */
        if (!strcmp(output_phase_sc, "SELF_VISPHI")) {
          v_factor_flag_sc= 0;
          p_factor_flag_sc= 0;
          cpl_msg_info(cpl_func, "vFactor correction for SC is %s due to option SELF_VISPHI", v_factor_flag_sc ? "ENABLE" : "DISABLE");
          cpl_msg_info(cpl_func, "pFactor correction for SC is %s due to option SELF_VISPHI", p_factor_flag_sc ? "ENABLE" : "DISABLE");
        } else {
		cpl_msg_info (cpl_func, "vFactor correction for SC is %s",v_factor_flag_sc?"ENABLE":"DISABLE");
		cpl_msg_info (cpl_func, "pFactor correction for SC is %s",p_factor_flag_sc?"ENABLE":"DISABLE");
        }

		CPLCHECK_NUL("Cannot get parameters");

        /* 
         * Loop on polarisations
         */
        for (int pol = 0; pol < npol_sc; pol++) {
            cpl_msg_info (cpl_func, "Start SC polarisation %d over %d",pol+1, npol_sc);
            
            /* Get the input table of SC */
            cpl_table * vis_SC = gravi_data_get_oi_vis (p2vmred_data, GRAVI_SC, pol, npol_sc);
            cpl_table * flux_SC = gravi_data_get_oi_flux (p2vmred_data, GRAVI_SC, pol, npol_sc);
            cpl_table * oi_wavelengthsc = gravi_data_get_oi_wave (p2vmred_data, GRAVI_SC, pol, npol_sc);
            CPLCHECK_NUL ("Cannot get data");

            /* Compute the wavenumber for SC in [m^-1] */
            cpl_array * wavenumber_sc;
            int nwave_sc = cpl_table_get_column_depth (vis_SC, "VISDATA");
            wavenumber_sc = cpl_array_new (nwave_sc, CPL_TYPE_DOUBLE);
            for (cpl_size wave = 0; wave < nwave_sc; wave ++){
                cpl_array_set (wavenumber_sc, wave, 1./cpl_table_get (oi_wavelengthsc, "EFF_WAVE", wave, &nv));
            }

            CPLCHECK_NUL ("Cannot build the wave and wavenumber");

            /* Create averaged tables */
            cpl_table * oi_vis2_SC = gravi_table_oi_create (nwave_sc, 1, GRAVI_OI_VIS2_EXT);
            gravi_table_new_column (oi_vis2_SC, "NDIT", NULL, CPL_TYPE_INT);
            gravi_table_new_column (oi_vis2_SC, "NVALID", NULL, CPL_TYPE_INT);
            
            cpl_table * oi_vis_SC = gravi_table_oi_create (nwave_sc, 1, GRAVI_OI_VIS_EXT);
            gravi_table_new_column (oi_vis_SC, "NDIT", NULL, CPL_TYPE_INT);
            gravi_table_new_column (oi_vis_SC, "NVALID", NULL, CPL_TYPE_INT);
            gravi_table_new_column (oi_vis_SC, "GDELAY", "m", CPL_TYPE_DOUBLE);
            gravi_table_new_column (oi_vis_SC, "PHASE", "rad", CPL_TYPE_DOUBLE);
            
            cpl_table * oi_T3_SC = gravi_table_oi_create (nwave_sc, 1, GRAVI_OI_T3_EXT);
            gravi_table_new_column (oi_T3_SC, "NDIT", NULL, CPL_TYPE_INT);
            gravi_table_new_column (oi_T3_SC, "NVALID", NULL, CPL_TYPE_INT);
            
            cpl_table * oi_flux_SC = gravi_table_oi_create (nwave_sc, 1, GRAVI_OI_FLUX_EXT);
            gravi_table_new_column (oi_flux_SC, "NDIT", NULL, CPL_TYPE_INT);
            gravi_table_new_column (oi_flux_SC, "NVALID", NULL, CPL_TYPE_INT);
            gravi_table_new_column (oi_flux_SC, "LKDT_MET_FC", "mjd", CPL_TYPE_DOUBLE);

            CPLCHECK_NUL("Cannot create columns in averaged OIFITS...");

            /* Keep only selected rows */
            vis_SC  = gravi_table_extract_time_interval (vis_SC, start_time, end_time);
            flux_SC = gravi_table_extract_time_interval (flux_SC, start_time, end_time);
                        
            /* 
             * Compute OIVIS2 and OIVIS for SC
             */
            cpl_msg_info (cpl_func, "Compute OIVIS2 and OIVIS for SC");
            
            gravi_vis_average_bootstrap (oi_vis_SC, oi_vis2_SC, vis_SC,
                                         nboot_sc,
                                         phase_ref_sc,
                                         v_factor_flag_sc,
                                         p_factor_flag_sc,
                                         debiasing_flag_sc,
					 outlier_threshold_sc);            
            CPLCHECK_NUL("Cannot average the SC frames");

            /* Compute other columns */
            gravi_vis_compute_column_mean (oi_vis_SC, vis_SC, "OPD_MET_FC", 6);
            gravi_vis_compute_column_mean (oi_vis_SC, vis_SC, "PHASE_REF_COEFF", 6);
            gravi_vis_compute_column_mean (oi_vis_SC, vis_SC, "E_U", 6);
            gravi_vis_compute_column_mean (oi_vis_SC, vis_SC, "E_V", 6);
            gravi_vis_compute_column_mean (oi_vis_SC, vis_SC, "E_W", 6);
            gravi_vis_compute_column_mean (oi_vis_SC, vis_SC, "E_AZ", 6);
            gravi_vis_compute_column_mean (oi_vis_SC, vis_SC, "E_ZD", 6);
            CPLCHECK_NUL("Cannot compute means.");
                
            /* Re compute the astrometric phase from VISDATA to deal with absolute phase
             * VISDATA as well as (R,I) remains unchanged. The goal is to split
             * the differential phase, calibratable so far, from the absolute phase */
            if ( !strcmp (output_phase_sc,"DIFFERENTIAL") ) {
            for (int base = 0; base < nbase; base++) {
                
                /* We duplicate VISDATA, to keep the value untouched in the table */
                cpl_array * visData_sc, * visErr_sc;
                visData_sc = cpl_array_cast (cpl_table_get_array (oi_vis_SC, "VISDATA", base),
                                             CPL_TYPE_DOUBLE_COMPLEX);
                
                /* Normalize the phasor by error, to discard crazy points */
                visErr_sc = cpl_array_duplicate (visData_sc);
                cpl_array_abs (visErr_sc);
                cpl_array_divide (visData_sc, visErr_sc);

		/* Get the flag */
		cpl_array * flag_sc;
                flag_sc = (cpl_array *)cpl_table_get_array (oi_vis_SC, "FLAG", base);
                
                /* Compute and remove the mean group delay in [m] */
                double mean_delay = 0.0;
                gravi_array_get_group_delay_loop (&visData_sc, &flag_sc, wavenumber_sc, &mean_delay, 1, 2e-3, CPL_FALSE);
                gravi_array_multiply_phasor (visData_sc, - 2*I*CPL_MATH_PI * mean_delay, wavenumber_sc);
                
                /* Save this delay [m] */
                cpl_table_set (oi_vis_SC, "GDELAY", base, mean_delay);
                
                /* Compute and remove the mean phase in [rad] */
                double mean_phase = carg (cpl_array_get_mean_complex (visData_sc));
                cpl_array_multiply_scalar_complex (visData_sc, cexp(- I * mean_phase));
                
                /* Save this phase [rad] */
                cpl_table_set (oi_vis_SC, "PHASE", base, mean_phase);
                
                /* Set back the phase in [deg] */
                cpl_array_arg (visData_sc);
                gravi_table_set_array_phase (oi_vis_SC, "VISPHI", base, visData_sc);
                cpl_array_delete (visData_sc);
                cpl_array_delete (visErr_sc);
                
                CPLCHECK_NUL("when computing the astrometric phase");

            } /* End loop on base */
            }
            if (!strcmp(output_phase_sc, "SELF_VISPHI")) {
              int* cmin=NULL;
              int* cmax=NULL;
              int nrange=0;
              /* rather complex C way to analyse strings defining clusters of wavelengths like " [ 12345 : 12346 , 678 : 680 , 822:864]" */
              if (strcmp(rangeString, "UNKNOWN")) {
              /*find number of ranges (they are separated by ',')*/
                char *str, *str1 ;
                int l=strlen(rangeString)+1;
                int j=0;
                int i=0;
                str=(char*)malloc(l);
                strncpy(str,rangeString,sizeof(str)-1); /*for future use*/
                for (i = 0; i<l; ++i) if (str[i]!=' ' && str[i]!='\t') str[j++]=str[i]; //remove ALL blanks.
                str[j]='\0';
                l=strlen(str)+1;
                str1=(char*)malloc(l);
                strncpy(str1,str,sizeof(str1)-1); /*for future use, as strtok destroys its arguments*/
            
                char *token;
                token = strtok(str, "[,]");
                while (token) {
                  nrange++;
                  token = strtok(NULL, "[,]");
                }
                if (nrange > 1) {
                  cmin=(int*) calloc(nrange,sizeof(int));
                  cmax=(int*) calloc(nrange,sizeof(int));

                  char *str2, *subtoken;
                  char *saveptr1, *saveptr2;
                  for (j = 0; ; j++, str1 = NULL) {
                      token = strtok_r(str1, "[,]" , &saveptr1);
                      if (token == NULL)
                         break;
                      for (str2 = token, i=0; i<2 ; str2 = NULL, ++i) {
                         subtoken = strtok_r(str2, ":", &saveptr2);
                         if (subtoken == NULL)
                             break;
                         /* int ret=sscanf(subtoken,"%d", (i==0)?&(cmin[j]):&(cmax[j])); */
                      }
                  }
                }
               }

               gravi_average_self_visphi(oi_vis_SC, vis_SC, wavenumber_sc, phase_ref_sc, cmin, cmax, nrange);
             }
            /* 
             * Compute OIT3 for SC
             */
            cpl_msg_info (cpl_func, "Compute OIT3 for SC");

            gravi_t3_average_bootstrap (oi_T3_SC, vis_SC, flux_SC,
                                        nboot_sc,
                                        v_factor_flag_sc,
                                        p_factor_flag_sc,
					outlier_threshold_sc);
            CPLCHECK_NUL("Cannot average t3 of SC");
            
            /* 
             * Compute OI_FLUX for SC
             */
            cpl_msg_info (cpl_func, "Compute OI_FLUX for SC");

            gravi_flux_average_bootstrap (oi_flux_SC, flux_SC,
                                          nboot_sc, outlier_threshold_sc);
            CPLCHECK_NUL("Cannot average flux of SC");
            
            /* Compute other columns */
            gravi_vis_compute_column_mean (oi_flux_SC, flux_SC, "OPD_MET_FC", 4);
            gravi_vis_compute_column_mean (oi_flux_SC, flux_SC, "FT_POS", 4);
            gravi_vis_compute_column_mean (oi_flux_SC, flux_SC, "SC_POS", 4);
            gravi_vis_compute_column_mean (oi_flux_SC, flux_SC, "OPL_AIR", 4);
            CPLCHECK_NUL ("Cannot compute mean columns");
            
            /* Save the FC metrology lock date */
            for (int tel = 0; tel < ntel; tel++){
                double lockdate = gravi_pfits_get_metfc_lockmjd (p2vmred_header, tel);
                cpl_table_set (oi_flux_SC, "LKDT_MET_FC", tel, lockdate);
            }

            /* 
             * Add tables in the vis_data
             */
            cpl_propertylist * vis_plist   = gravi_data_get_oi_vis_plist (p2vmred_data, GRAVI_SC, pol, npol_sc);
            cpl_propertylist * oivis_plist = cpl_propertylist_new();
            cpl_propertylist_copy_property (oivis_plist, vis_plist, "DATE-OBS");
            cpl_propertylist_copy_property (oivis_plist, vis_plist, "EXTVER");
            cpl_propertylist_copy_property (oivis_plist, vis_plist, "ARRNAME");
            cpl_propertylist_copy_property (oivis_plist, vis_plist, "INSNAME");
            cpl_propertylist_update_string (oivis_plist, "AMPTYP","absolute");
            if ( !strcmp (output_phase_sc,"DIFFERENTIAL") )
                cpl_propertylist_update_string (oivis_plist, "PHITYP","differential");
            if ( !strcmp (output_phase_sc,"ABSOLUTE") )
                cpl_propertylist_update_string (oivis_plist, "PHITYP","absolute");
            cpl_propertylist_update_int (oivis_plist, "PHIORDER",1);
            gravi_data_add_table (vis_data, oivis_plist, GRAVI_OI_VIS_EXT, oi_vis_SC);
            
            cpl_propertylist * oivis2_plist = cpl_propertylist_new();
            cpl_propertylist_copy_property (oivis2_plist, vis_plist, "DATE-OBS");
            cpl_propertylist_copy_property (oivis2_plist, vis_plist, "EXTVER");
            cpl_propertylist_copy_property (oivis2_plist, vis_plist, "ARRNAME");
            cpl_propertylist_copy_property (oivis2_plist, vis_plist, "INSNAME");
            gravi_data_add_table (vis_data, oivis2_plist, GRAVI_OI_VIS2_EXT, oi_vis2_SC);
            
            cpl_propertylist * oit3_plist = cpl_propertylist_new();
            cpl_propertylist_copy_property (oit3_plist, vis_plist, "DATE-OBS");
            cpl_propertylist_copy_property (oit3_plist, vis_plist, "EXTVER");
            cpl_propertylist_copy_property (oit3_plist, vis_plist, "ARRNAME");
            cpl_propertylist_copy_property (oit3_plist, vis_plist, "INSNAME");
            gravi_data_add_table (vis_data, oit3_plist, GRAVI_OI_T3_EXT, oi_T3_SC);
            
            cpl_propertylist * oiflux_plist = cpl_propertylist_new();
            cpl_propertylist_copy_property (oiflux_plist, vis_plist, "DATE-OBS");
            cpl_propertylist_copy_property (oiflux_plist, vis_plist, "EXTVER");
            cpl_propertylist_copy_property (oiflux_plist, vis_plist, "ARRNAME");
            cpl_propertylist_copy_property (oiflux_plist, vis_plist, "INSNAME");
            cpl_propertylist_update_string (oiflux_plist, "CALSTAT", "U");
            cpl_propertylist_set_comment (oiflux_plist, "CALSTAT", "Uncalibrated flux per telescope");
            gravi_data_add_table (vis_data, oiflux_plist, GRAVI_OI_FLUX_EXT, oi_flux_SC);
            
            CPLCHECK_NUL ("Cannot add tables");
            
            /* Delete waves */
            FREE (cpl_array_delete, wavenumber_sc);
            CPLCHECK_NUL ("Cannot delete wavenumber");

            FREE (cpl_table_delete, vis_SC);
            FREE (cpl_table_delete, flux_SC);
            
        } /* end loop on pol */
    } /* End SC */
    
        
    /* Add the night-obs in the main HEADER. This is the MJD of the begining
       of the night, similar for all files from noon to noon (~local time) */
    double mjd_obs = cpl_table_get_column_mean (gravi_data_get_table (vis_data, GRAVI_OI_VIS_EXT), "MJD");
    cpl_propertylist_update_int (vis_header, GRAVI_NIGHT_OBS, (int)floor(mjd_obs-0.625));

	gravi_msg_function_exit(1);
	return vis_data;
}


/*----------------------------------------------------------------------------*/
/**
 * @brief The function compute the QC parameters for a VIS (averaged) data
 * 
 * @param vis_data  VIS data containing OI_FLUX and OI_VIS tables
 * 	  	            comming from the @c gravi_compute_vis
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_compute_vis_qc (gravi_data * vis_data)
{
	gravi_msg_function_start(1);
	cpl_ensure_code (vis_data, CPL_ERROR_NULL_INPUT);
	
	int nv, nbase = 6, ntel=4, nclo=4;
	char qc_name[100];

	/* 
	 * Prepare the output 
	 */
	
    cpl_propertylist * vis_header = gravi_data_get_header (vis_data);
    cpl_propertylist * plist = gravi_data_get_header (vis_data);

    
    /* 
     * Start with FT 
     */
    if (gravi_data_has_type (vis_data, "_FT") <= 0 ) {
		cpl_msg_info (cpl_func, "VIS data has no FT extensions");
    }
    else {

        /* Loop on polarisations */
        int npol_ft = gravi_pfits_get_pola_num (vis_header, GRAVI_FT);
        for (int pol = 0; pol < npol_ft; pol++) {
            cpl_msg_info (cpl_func, "Start FT polarisation %d over %d",pol+1, npol_ft);
            
            /* Loop on bases to compute OIVIS2 and OIVIS for FT
             */
            cpl_msg_info (cpl_func, "Compute QC OIVIS2 and OIVIS for FT");
            
            cpl_table * oi_vis2_FT = gravi_data_get_oi_vis2 (vis_data, GRAVI_FT, pol, npol_ft);
            cpl_table * oi_vis_FT = gravi_data_get_oi_vis (vis_data, GRAVI_FT, pol, npol_ft);

            for (int base = 0; base < nbase; base++) {
                
                /* Add the QC parameters for FT */
                
                sprintf (qc_name, "ESO QC VISPHIERR_FT%s_P%d AVG", GRAVI_BASE_NAME[base], pol+1);
                cpl_propertylist_update_double (plist, qc_name, gravi_table_get_column_mean (oi_vis_FT, "VISPHIERR", base, nbase));
                cpl_propertylist_set_comment (plist, qc_name, "[deg] mean over lbd");
                
                sprintf (qc_name, "ESO QC VIS2_FT%s_P%d AVG", GRAVI_BASE_NAME[base], pol+1);
                cpl_propertylist_update_double (plist, qc_name, gravi_table_get_column_mean (oi_vis2_FT, "VIS2DATA", base, nbase));
                cpl_propertylist_set_comment (plist, qc_name, "mean over lbd");
                
                sprintf (qc_name, "ESO QC VIS2ERR_FT%s_P%d AVG", GRAVI_BASE_NAME[base], pol+1);
                cpl_propertylist_update_double (plist, qc_name, gravi_table_get_column_mean (oi_vis2_FT, "VIS2ERR", base, nbase));
                cpl_propertylist_set_comment (plist, qc_name, "mean over lbd");
                
                sprintf (qc_name, "ESO QC VISAMP_FT%s_P%d AVG", GRAVI_BASE_NAME[base], pol+1);
                cpl_propertylist_update_double (plist, qc_name, gravi_table_get_column_mean (oi_vis_FT, "VISAMP", base, nbase));
                cpl_propertylist_set_comment (plist, qc_name, "mean over lbd");
                
                sprintf (qc_name, "ESO QC VISAMPERR_FT%s_P%d AVG", GRAVI_BASE_NAME[base], pol+1);
                cpl_propertylist_update_double (plist, qc_name, gravi_table_get_column_mean (oi_vis_FT, "VISAMPERR", base, nbase));
                cpl_propertylist_set_comment (plist, qc_name, "mean over lbd");
                
                CPLCHECK_MSG("Cannot compute QC parameter for OI_VIS for FT");
            } /* End loop on base */

            /* 
             * Loop on triplet to compute OIT3 for FT
             */
            cpl_msg_info (cpl_func, "Compute QC OIT3 for FT");

            cpl_table * oi_T3_FT = gravi_data_get_oi_t3 (vis_data, GRAVI_FT, pol, npol_ft);
            
            for (int clo = 0; clo < nclo; clo++){
                
                sprintf (qc_name, "ESO QC T3PHI_FT%s_P%d AVG", GRAVI_CLO_NAME[clo], pol+1);
                cpl_propertylist_update_double (plist, qc_name, gravi_table_get_column_mean (oi_T3_FT, "T3PHI", clo, nclo));
                cpl_propertylist_set_comment (plist, qc_name, "[deg] mean over lbd");
                
                sprintf (qc_name, "ESO QC T3PHIERR_FT%s_P%d AVG", GRAVI_CLO_NAME[clo], pol+1);
                cpl_propertylist_update_double (plist, qc_name, gravi_table_get_column_mean (oi_T3_FT, "T3PHIERR", clo, nclo));
                cpl_propertylist_set_comment (plist, qc_name, "[deg] mean over lbd");
                
                CPLCHECK_MSG("Cannot compute QC parameter for OI_T3 for FT");
            } /* End loop on triplets */
            
            /* 
             * Loop on beams to compute OI_FLUX for FT
             */
            cpl_msg_info (cpl_func, "Compute QC OI_FLUX for FT");
		
            cpl_table * oi_flux_FT = gravi_data_get_oi_flux (vis_data, GRAVI_FT, pol, npol_ft);
            
            for (int tel = 0; tel < ntel; tel++){
                
                sprintf (qc_name, "ESO QC FLUX_FT%d_P%d AVG", tel+1, pol+1);
                cpl_propertylist_update_double (plist, qc_name, gravi_table_get_column_mean (oi_flux_FT, "FLUX", tel, ntel));
                cpl_propertylist_set_comment (plist, qc_name, "[e/total_int_time] mean over lbd");
                
                sprintf (qc_name, "ESO QC FLUXERR_FT%d_P%d AVG", tel+1, pol+1);
                cpl_propertylist_update_double (plist, qc_name, gravi_table_get_column_mean (oi_flux_FT, "FLUXERR", tel, ntel));
                cpl_propertylist_set_comment (plist, qc_name, "[e/total_int_time] mean over lbd");
                
                sprintf (qc_name, "ESO QC FLUXRATE_FT%d_P%d SUM", tel+1, pol+1);
                double flux_rate = cpl_array_get_mean (cpl_table_get_array (oi_flux_FT, "FLUX", tel)) *
                    cpl_array_get_size (cpl_table_get_array (oi_flux_FT, "FLUX", tel)) / cpl_table_get_double (oi_flux_FT, "INT_TIME", tel, &nv);
                cpl_propertylist_update_double (plist, qc_name, flux_rate);
                cpl_propertylist_set_comment (plist, qc_name, "[e/s] sum over lbd");
                
                CPLCHECK_MSG("Cannot compute QC parameter for OI_FLUX for FT");
            } /* End loop on beams */
            
        } /* end loop on pol */
    } /* End FT */


    
    /* 
     * Then with SC 
     */
    if (gravi_data_has_type (vis_data, "_SC") <= 0 ) {
		cpl_msg_info (cpl_func, "VIS data has no SC extensions");
    }
    else {
        
        /* Loop on polarisations */
        int npol_sc = gravi_pfits_get_pola_num (vis_header, GRAVI_SC);
        for (int pol = 0; pol < npol_sc; pol++) {

            /* 
             * Loop on bases to compute OIVIS2 and OIVIS for SC
             */
            cpl_msg_info (cpl_func, "Compute QC OIVIS2 and OIVIS for SC");

            cpl_table * oi_vis2_SC = gravi_data_get_oi_vis2 (vis_data, GRAVI_SC, pol, npol_sc);
            cpl_table * oi_vis_SC = gravi_data_get_oi_vis (vis_data, GRAVI_SC, pol, npol_sc);
            
            for (int base = 0; base < nbase; base++) {
                /* FIXME: repair these QC parameters, for instance by computing them in P2VMRED */
                
                // sprintf (qc_name, "ESO QC VFACTOR%s_P%d AVG", GRAVI_BASE_NAME[base], pol+1);
                // double vmean = gravi_table_get_column_mean (vis_SC, "V_FACTOR_WL", base, nbase);
                // cpl_propertylist_update_double (plist, qc_name, vmean);
                // cpl_propertylist_set_comment (plist, qc_name, "mean v-factor");
                // 
                // sprintf (qc_name, "ESO QC PFACTOR%s_P%d AVG", GRAVI_BASE_NAME[base], pol+1);
                // double pmean = gravi_table_get_column_mean (vis_SC, "P_FACTOR", base, nbase);
                // cpl_propertylist_update_double (plist, qc_name, pmean);
                // cpl_propertylist_set_comment (plist, qc_name, "mean p-factor");
                
                sprintf (qc_name, "ESO QC GD_SC%s_P%d AVG", GRAVI_BASE_NAME[base], pol+1);
                cpl_propertylist_update_double (plist, qc_name, gravi_table_get_column_mean (oi_vis_SC, "GDELAY", base, nbase));
                cpl_propertylist_set_comment (plist, qc_name, "[m] mean Group-Delay");
                
                sprintf (qc_name, "ESO QC VIS2_SC%s_P%d AVG", GRAVI_BASE_NAME[base], pol+1);
                cpl_propertylist_update_double (plist, qc_name, gravi_table_get_column_mean (oi_vis2_SC, "VIS2DATA", base, nbase));
                cpl_propertylist_set_comment (plist, qc_name, "mean over lbd");
                
                sprintf (qc_name, "ESO QC VIS2ERR_SC%s_P%d AVG", GRAVI_BASE_NAME[base], pol+1);
                cpl_propertylist_update_double (plist, qc_name, gravi_table_get_column_mean (oi_vis2_SC, "VIS2ERR", base, nbase));
                cpl_propertylist_set_comment (plist, qc_name, "mean over lbd");
                
                sprintf (qc_name, "ESO QC VISPHI_SC%s_P%d AVG", GRAVI_BASE_NAME[base], pol+1);
                cpl_propertylist_update_double (plist, qc_name, gravi_table_get_column_mean (oi_vis_SC, "VISPHI", base, nbase));
                cpl_propertylist_set_comment (plist, qc_name, "[deg] mean over lbd");
                
                sprintf (qc_name, "ESO QC VISPHIERR_SC%s_P%d AVG", GRAVI_BASE_NAME[base], pol+1);
                cpl_propertylist_update_double (plist, qc_name, gravi_table_get_column_mean (oi_vis_SC, "VISPHIERR", base, nbase));
                cpl_propertylist_set_comment (plist, qc_name, "[deg] mean over lbd");
                
                sprintf (qc_name, "ESO QC VISAMP_SC%s_P%d AVG", GRAVI_BASE_NAME[base], pol+1);
                cpl_propertylist_update_double (plist, qc_name, gravi_table_get_column_mean (oi_vis_SC, "VISAMP", base, nbase));
                cpl_propertylist_set_comment (plist, qc_name, "mean over lbd");
                
                sprintf (qc_name, "ESO QC VISAMPERR_SC%s_P%d AVG", GRAVI_BASE_NAME[base], pol+1);
                cpl_propertylist_update_double (plist, qc_name, gravi_table_get_column_mean (oi_vis_SC, "VISAMPERR", base, nbase));
                cpl_propertylist_set_comment (plist, qc_name, "mean over lbd");
                
                double coeff2 = gravi_table_get_value (oi_vis_SC, "PHASE_REF_COEFF", base, 2);
                sprintf (qc_name, "ESO QC PHASE_REF_COEFF2 SC%s_P%d", GRAVI_BASE_NAME[base], pol+1);
                cpl_propertylist_update_double (plist, qc_name, coeff2);
                cpl_propertylist_set_comment (plist, qc_name, "[rad] 2sd order of FT phase");
                
                CPLCHECK_MSG("Cannot set QC parameter for OI_VIS for SC");
            } /* End loop on base */
            
            /* 
             * Loop on triplet to compute OIT3 for SC
             */
            cpl_msg_info (cpl_func, "Compute QC OIT3 for SC");
            
            cpl_table * oi_T3_SC = gravi_data_get_oi_t3 (vis_data, GRAVI_SC, pol, npol_sc);
            
            for (int clo = 0; clo < nclo; clo++){
                
                sprintf (qc_name, "ESO QC T3PHI_SC%s_P%d AVG", GRAVI_CLO_NAME[clo], pol+1);
                cpl_propertylist_update_double (plist, qc_name, gravi_table_get_column_mean (oi_T3_SC, "T3PHI", clo, nclo));
                cpl_propertylist_set_comment (plist, qc_name, "[deg] mean over lbd");
                
                sprintf (qc_name, "ESO QC T3PHIERR_SC%s_P%d AVG", GRAVI_CLO_NAME[clo], pol+1);
                cpl_propertylist_update_double (plist, qc_name, gravi_table_get_column_mean (oi_T3_SC, "T3PHIERR", clo, nclo));
                cpl_propertylist_set_comment (plist, qc_name, "[deg] mean over lbd");
          
                CPLCHECK_MSG("Cannot set QC parameter for OI_T3 for SC");
            }/* End loop on triplets */
            
            /* 
             * Loop on beams to compute OI_FLUX for SC
             */
            cpl_msg_info (cpl_func, "Compute QC OI_FLUX for SC");
            
            cpl_table * oi_flux_SC = gravi_data_get_oi_flux (vis_data, GRAVI_SC, pol, npol_sc);
            
            for (int tel = 0; tel < ntel; tel++){
                
                sprintf (qc_name, "ESO QC FLUX_SC%d_P%d AVG", tel+1, pol+1);
                cpl_propertylist_update_double (plist, qc_name, gravi_table_get_column_mean (oi_flux_SC, "FLUX", tel, ntel));
                cpl_propertylist_set_comment (plist, qc_name, "[e/total_int_time] mean over lbd");
                
                sprintf (qc_name, "ESO QC FLUXERR_SC%d_P%d AVG", tel+1, pol+1);
                cpl_propertylist_update_double (plist, qc_name, gravi_table_get_column_mean (oi_flux_SC, "FLUXERR", tel, ntel));
                cpl_propertylist_set_comment (plist, qc_name, "[e/total_int_time] mean over lbd");
                
                sprintf (qc_name, "ESO QC FLUXRATE_SC%d_P%d SUM", tel+1, pol+1);
                double flux_rate = cpl_array_get_mean (cpl_table_get_array (oi_flux_SC, "FLUX", tel)) *
                    cpl_array_get_size(cpl_table_get_array (oi_flux_SC, "FLUX", tel)) / cpl_table_get_double (oi_flux_SC, "INT_TIME", tel, &nv);
                cpl_propertylist_update_double (plist, qc_name, flux_rate);
                cpl_propertylist_set_comment (plist, qc_name, "[e/s] sum over lbd");
                
                double ftpos_mean = cpl_table_get (oi_flux_SC, "FT_POS", tel, NULL);
                sprintf (qc_name, "ESO QC FT_POS SC%d_P%d", tel+1, pol+1);
                cpl_propertylist_update_double (plist, qc_name, ftpos_mean);
                cpl_propertylist_set_comment (plist, qc_name, "[V]");
                
                double oplair_mean = cpl_table_get (oi_flux_SC, "OPL_AIR", tel, NULL);
                sprintf (qc_name, "ESO QC OPL_AIR SC%d_P%d", tel+1, pol+1);
                cpl_propertylist_update_double (plist, qc_name, oplair_mean);
                cpl_propertylist_set_comment (plist, qc_name, "[m]");
                
                CPLCHECK_MSG("Cannot set QC parameter for OI_FLUX for SC");
            } /* End loop on beams */


        } /* end loop on pol */
    } /* End SC */
    
	gravi_msg_function_exit(1);
    return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Align the SC visibilities on the FT visibilities
 * 
 * @param vis_data: the VIS data to update in-place
 * 
 * Force the mean value of the VIS2DATA of the SC to match the mean
 * value of the VIS2DATA of the FT. Same for VISAMP.
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_normalize_sc_to_ft (gravi_data * vis_data)
{
  gravi_msg_function_start(1);
  cpl_ensure_code (vis_data, CPL_ERROR_NULL_INPUT);
  
  cpl_propertylist * hdr_data = gravi_data_get_header (vis_data);
  int pol, npol = gravi_pfits_get_pola_num (hdr_data, GRAVI_SC);
  double qFactor;

  if ( npol != gravi_pfits_get_pola_num (hdr_data, GRAVI_FT)) {
	return cpl_error_set_message(cpl_func,CPL_ERROR_ILLEGAL_INPUT,"polarisation of SC and FT shall be compatible");
  }
  
  /* Loop on polarisation. Assume the same polarisation
   * splitting for both SC and FT */
  for ( pol= 0 ; pol < npol ; pol++ ) {

	cpl_table * oi_wave_sc = gravi_data_get_oi_wave (vis_data, GRAVI_SC, pol, npol);
	cpl_table * oi_wave_ft = gravi_data_get_oi_wave (vis_data, GRAVI_FT, pol, npol);

	/* 
	 * Get data for VIS2 and VISAMP
	 */
	cpl_table * oi_vis2_sc = gravi_data_get_oi_vis2 (vis_data, GRAVI_SC, pol, npol);
	cpl_table * oi_vis2_ft = gravi_data_get_oi_vis2 (vis_data, GRAVI_FT, pol, npol);
	cpl_table * oi_vis_sc = gravi_data_get_oi_vis (vis_data, GRAVI_SC, pol, npol);
	cpl_table * oi_vis_ft = gravi_data_get_oi_vis (vis_data, GRAVI_FT, pol, npol);

	CPLCHECK_MSG("Cannot get data");
	
	/* Loop on baselines */
	for (cpl_size base = 0; base < cpl_table_get_nrow (oi_vis2_ft); base ++) {

	  /* Get the VIS2 of FT */
	  const cpl_array * vis2_ft = cpl_table_get_array (oi_vis2_ft, "VIS2DATA", base);
	  
	  /* Create the rebin VIS2 of SC */
	  cpl_array * vis2_lr = gravi_array_rebin (cpl_table_get_array (oi_vis2_sc, "VIS2DATA", base),
											   cpl_table_get_array (oi_vis2_sc, "VIS2ERR", base),
											   oi_wave_sc, oi_wave_ft);
	  
	  /* Compute the mean visibility loss over the band.
	   * FIXME: shall use  vis2_ft = (a.lbd+b) * vis2_lr, which could be done with cpl_matrix */
	  qFactor = cpl_array_get_mean (vis2_lr) / cpl_array_get_mean (vis2_ft);
	  cpl_msg_info (cpl_func, "vis2 %lli: qFactor = %f", base, qFactor);

	  /* Divide the SC data */
	  cpl_array_divide_scalar (cpl_table_get_data_array (oi_vis2_sc, "VIS2DATA")[base], qFactor);
	  cpl_array_divide_scalar (cpl_table_get_data_array (oi_vis2_sc, "VIS2ERR")[base], qFactor);
	  
	  /* Get the VISAMP of FT */
	  const cpl_array * vis_ft = cpl_table_get_array (oi_vis_ft, "VISAMP", base);
	  
	  /* Create the rebin VISAMP of SC */
	  cpl_array * vis_lr = gravi_array_rebin (cpl_table_get_array (oi_vis_sc, "VISAMP", base),
											  cpl_table_get_array (oi_vis_sc, "VISAMPERR", base),
											  oi_wave_sc, oi_wave_ft);
	  
	  /* Compute the mean visibility loss over the band. */
	  qFactor = cpl_array_get_mean (vis_lr) / cpl_array_get_mean (vis_ft);
	  cpl_msg_info (cpl_func, "visAmp %lli: qFactor = %f", base, qFactor);

	  /* Divide the SC data */
	  cpl_array_divide_scalar (cpl_table_get_data_array (oi_vis_sc, "VISAMP")[base], qFactor);
	  cpl_array_divide_scalar (cpl_table_get_data_array (oi_vis_sc, "VISAMPERR")[base], qFactor);
	  
	  FREE (cpl_array_delete, vis2_lr);
	  FREE (cpl_array_delete, vis_lr);
	} /* End loop on baselines */
	
  } /* End loop on pol */
  
  gravi_msg_function_exit(1);
  return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Recompute the TIME column of all OIFITS extension from the MJD
 *        column, following the OIFITS standard (number of second since
 *        the DATE-OBS at 00:00). 
 * 
 * @param vis_data: the gravi_data to update in-place
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_vis_mjd_to_time (gravi_data * vis_data)
{
  gravi_msg_function_start(1);
  cpl_ensure_code (vis_data, CPL_ERROR_NULL_INPUT);
  char date[90];

  /* Loop on extension */
  int next = gravi_data_get_size (vis_data);
  for (int ext = 0; ext < next; ext ++) {

	const char * extname = gravi_data_get_extname (vis_data,ext);
	if (!strcmp (extname, "OI_VIS") ||
		!strcmp (extname, "OI_VIS2") ||
		!strcmp (extname, "OI_T3") ||
		!strcmp (extname, "OI_FLUX")) {

	  /* Get the DATE-OBS in format YYYY-MM-DDT00:00:00.000 */
	  cpl_propertylist * plist = gravi_data_get_plist_x (vis_data, ext);
	  sprintf (date, "%.10sT00:00:00.000", cpl_propertylist_get_string (plist, "DATE-OBS"));

	  /* Get the MJD of this DATE-OBS */
	  double mjd0 = gravi_convert_to_mjd (date);
	  cpl_msg_debug (cpl_func, "DATE-OBS = %s  ->  mjd = %.3f", date, mjd0);

	  /* Compute TIME in [s] following the OIFITS standard */
	  cpl_table * oi_table = gravi_data_get_table_x (vis_data, ext);
	  cpl_size nrow = cpl_table_get_nrow (oi_table);
	  for (cpl_size row = 0; row < nrow; row++) {
		double mjd = cpl_table_get (oi_table, "MJD", row, NULL);
		cpl_table_set (oi_table, "TIME", row, (mjd-mjd0) * 24 * 3600);
	  }

	  /* Set units */
	  cpl_table_set_column_unit (oi_table, "TIME", "s");
	}
  } /* End loop on extensions */
  
  gravi_msg_function_exit(1);
  return CPL_ERROR_NONE;
}


/*----------------------------------------------------------------------------*/
/**
 * @brief Divide the OI_FLUX by OI_FLUX from the P2VM
 *        (no checks, no time distance...)
 * 
 * @param vis_data:   the VIS data whose OI_FLUX tables
 *                    will be updated in-place
 * @param p2vm_map:   the input P2VM calibration map to
 *                    load the reference OI_FLUX
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_flat_flux (gravi_data * vis_data, gravi_data * p2vm_map)
{
	gravi_msg_function_start(1);
    cpl_ensure_code (vis_data,  CPL_ERROR_NULL_INPUT);
    cpl_ensure_code (p2vm_map, CPL_ERROR_NULL_INPUT);

	int num_used_tf = 1;
	gravi_data **used_tf_data = &p2vm_map;

	/* Get the header */
	cpl_propertylist * hdr_data = gravi_data_get_header (vis_data);
	
	/* For each type of data SC / FT */
	int ntype_data = 2;
	for (int type_data = 0; type_data < ntype_data ; type_data ++) {

	  /* Loop on polarisation */
	  int npol = gravi_pfits_get_pola_num (hdr_data, type_data);
	  for (int pol= 0 ; pol < npol ; pol++ ) {

		/* Calibrate the FLUX as a real quantity */
		double delta_t = 10000.0;
		gravi_apply_tf_amp (vis_data, NULL, used_tf_data, num_used_tf,
							GRAVI_OI_FLUX_EXT,
							GRAVI_INSNAME(type_data, pol, npol),
							"FLUX", "FLUXERR", 4, delta_t);
		
		CPLCHECK_MSG("Cannot apply normalize flux");

	  }
	  /* End loop on polarisation */
	}
	/* End loop on data_type */

	gravi_msg_function_exit(1);
	return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Average amplitudes column of a multi-observation OIFITS table
 *        The averaged quantities are stored in the first nbase rows
 * 
 * @param oi_table:   the table to update (number of row will be reduced)
 * @param name:       name of column to average
 * @param err:        corresponding error column
 * @param nbase:      number of row per obs (ex: 6 for VIS and 4 for FLUX)
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_vis_average_amp (cpl_table *oi_table, const char *name,  const char *err, int nbase)
{
  cpl_ensure_code (oi_table, CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (name,     CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (err,      CPL_ERROR_NULL_INPUT);

  cpl_size nwave = cpl_table_get_column_depth (oi_table, name);
  cpl_array * weight = gravi_array_init_double (nwave, 0.0);
  cpl_array * value  = gravi_array_init_double (nwave, 0.0);

  /* Loop on base */
  int nrow = cpl_table_get_nrow (oi_table) / nbase;
  for (cpl_size base = 0; base < nbase ; base++) {
	cpl_array_fill_window (weight, 0, nwave, 0.0);
	cpl_array_fill_window (value,  0, nwave, 0.0);

	/* Loop on row and wave */
	for (cpl_size row = 0; row < nrow ; row++) {
	  const cpl_array * rval = cpl_table_get_array (oi_table, name, base + row*nbase);
	  const cpl_array * rerr = cpl_table_get_array (oi_table, err,  base + row*nbase);
	  const cpl_array * flag = cpl_table_get_array (oi_table, "FLAG", base + row*nbase);
	  for (cpl_size wave = 0; wave < nwave; wave++) {
		double w = pow (cpl_array_get (rerr, wave, NULL), -2);
		if (cpl_array_get (flag, wave, NULL)) w = 10e-20;
		double v = cpl_array_get (rval, wave, NULL);
		cpl_array_set (weight, wave, cpl_array_get (weight, wave, NULL) + w);
		cpl_array_set (value,  wave, cpl_array_get (value, wave, NULL) + v * w);
	  }
	}
	CPLCHECK_MSG("Cannot average amp");

	/* Set the mean */
	cpl_array_divide (value, weight);
	cpl_table_set_array (oi_table, name, base, value);
    
	/* Set the variance of the mean */
	cpl_array_power (weight, -0.5);
	cpl_table_set_array (oi_table, err,  base, weight);


	
  } /* End loop on base */

  FREE (cpl_array_delete, weight);
  FREE (cpl_array_delete, value);
  return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Average phases column of a multi-observation OIFITS table
 *        Phases are averaged with arg{<exp(i.phi)>}
 *        The averaged quantities are stored in the first nbase rows
 * 
 * @param oi_table:   the table to update (number of row will be reduced)
 * @param name:       name of column to average
 * @param err:        corresponding error column
 * @param nbase:      number of row per obs (ex: 6 for VIS and 4 for FLUX)
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_vis_average_phi (cpl_table *oi_table, const char *name,  const char *err, int nbase)
{
  cpl_ensure_code (oi_table, CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (name,     CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (err,      CPL_ERROR_NULL_INPUT);

  cpl_size nwave = cpl_table_get_column_depth (oi_table, name);
  cpl_array * weight = gravi_array_init_double (nwave, 0.0);
  cpl_array * value  = gravi_array_init_double_complex (nwave, 0.0 + I*0.0);
  
  /* Loop on base */
  int nrow = cpl_table_get_nrow (oi_table) / nbase;
  for (cpl_size base = 0; base < nbase ; base++) {
	cpl_array_fill_window (weight, 0, nwave, 0.0);
	cpl_array_fill_window_complex (value,  0, nwave, 0.0);
	
	/* Loop on row and wave */
	for (cpl_size row = 0; row < nrow ; row++) {
	  const cpl_array * rval = cpl_table_get_array (oi_table, name, base + row*nbase);
	  const cpl_array * rerr = cpl_table_get_array (oi_table, err,  base + row*nbase);
	  const cpl_array * flag = cpl_table_get_array (oi_table, "FLAG", base + row*nbase);
	  for (cpl_size wave = 0; wave < nwave; wave++) {
		double w = pow (cpl_array_get (rerr, wave, NULL), -2);
		if (cpl_array_get (flag, wave, NULL)) w = 10e-20;
		double complex v = cexp (1.*I * cpl_array_get (rval, wave, NULL) * CPL_MATH_RAD_DEG);
		cpl_array_set (weight, wave, cpl_array_get (weight, wave, NULL) + w);
		cpl_array_set_complex (value,  wave, cpl_array_get_complex (value, wave, NULL) + v * w);
	  }
	}
	CPLCHECK_MSG("Cannot average phi");

	/* Set the mean */
	gravi_table_set_array_phase (oi_table, name, base, value);
	/* Set the variance of the mean */
	cpl_array_power (weight, -0.5);
	cpl_table_set_array (oi_table, err,  base, weight);
	
  } /* End loop on base */

  FREE (cpl_array_delete, weight);
  FREE (cpl_array_delete, value);
  return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Average scalar column of a multi-observation OIFITS table.
 *        The averaged quantities are stored in the first nbase rows
 * 
 * @param oi_table:   the table to update (number of row will be reduced)
 * @param name:       name of column to average
 * @param err:        corresponding error column
 * @param nbase:      number of row per obs (ex: 6 for VIS and 4 for FLUX)
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_vis_average_value (cpl_table *oi_table, const char *name,  const char *err, int nbase)
{
  cpl_ensure_code (oi_table, CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (name,     CPL_ERROR_NULL_INPUT);
  int nv = 0;

  /* Loop on base */
  int nrow = cpl_table_get_nrow (oi_table) / nbase;
  for (cpl_size base = 0; base < nbase ; base++) {
	double weight = 0.0;
	double value  = 0.0;
	
	/* Loop on row */
	for (cpl_size row = 0; row < nrow ; row++) {
	  double w = (err!=NULL ? pow (cpl_table_get (oi_table, err, base + row*nbase, &nv), -2.0) : 1.0);
	  value  += cpl_table_get (oi_table, name, base + row*nbase, &nv) * w;
	  weight += w;
	}
	CPLCHECK_MSG("Cannot average value");

	/* Set the mean */
	cpl_table_set (oi_table, name, base, value / weight);
  } /* End loop on base */

  return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Force uncertainties
 * 
 * @param oi_data: the OIFITS map to co-add in-place
 * 
 * Read the parameter of the recipe and update the uncertainties on
 * the corresponding quantities.
 */
/*----------------------------------------------------------------------------*/
cpl_error_code gravi_force_uncertainties (gravi_data * oi_data,
                                          const cpl_parameterlist * parlist)
{
  gravi_msg_function_start(1);
  cpl_ensure_code (oi_data, CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (parlist, CPL_ERROR_NULL_INPUT);

  double err;

  /* Get header */
  cpl_propertylist * header = gravi_data_get_header (oi_data);
  
  /* VISAMPERR SC */
  err = gravi_param_get_double_default (parlist, "gravity.postprocess.visamperr-sc", -1.0);
  if ( err > 0) {
      cpl_msg_info (cpl_func,"Force VISAMPERR to %.e for all observations", err);

      int npol = gravi_pfits_get_pola_num (header, GRAVI_SC);
      for (int pol = 0; pol < npol; pol++) {
          cpl_table *  oi_table = gravi_data_get_oi_vis (oi_data, GRAVI_SC, pol, npol);
          cpl_array ** array = cpl_table_get_data_array (oi_table, "VISAMPERR");
          for (cpl_size row = 0; row<cpl_table_get_nrow (oi_table); row++) {
              cpl_array_fill_window (array[row], 0, CPL_SIZE_MAX, err);
          }
      }
  }

  /* VISPHIERR SC */
  err = gravi_param_get_double_default (parlist, "gravity.postprocess.visphierr-sc", -1.0);
  if ( err > 0) {
      cpl_msg_info (cpl_func,"Force VISPHIERR to %.e for all observations", err);

      int npol = gravi_pfits_get_pola_num (header, GRAVI_SC);
      for (int pol = 0; pol < npol; pol++) {
          cpl_table *  oi_table = gravi_data_get_oi_vis (oi_data, GRAVI_SC, pol, npol);
          cpl_array ** array = cpl_table_get_data_array (oi_table, "VISPHIERR");
          for (cpl_size row = 0; row<cpl_table_get_nrow (oi_table); row++) {
              cpl_array_fill_window (array[row], 0, CPL_SIZE_MAX, err);
          }
      }
  }

  /* FLUXERR SC */
  err = gravi_param_get_double_default (parlist, "gravity.postprocess.fluxerr-sc", -1.0);
  if ( err > 0) {
      cpl_msg_info (cpl_func,"Force FLUXERR to %.e for all observations", err);

      int npol = gravi_pfits_get_pola_num (header, GRAVI_SC);
      for (int pol = 0; pol < npol; pol++) {
          cpl_table *  oi_table = gravi_data_get_oi_flux (oi_data, GRAVI_SC, pol, npol);
          cpl_array ** array = cpl_table_get_data_array (oi_table, "FLUXERR");
          for (cpl_size row = 0; row<cpl_table_get_nrow (oi_table); row++) {
              cpl_array_fill_window (array[row], 0, CPL_SIZE_MAX, err);
          }
      }
  }

  /* VIS2ERR SC */
  err = gravi_param_get_double_default (parlist, "gravity.postprocess.vis2err-sc", -1.0);
  if ( err > 0) {
      cpl_msg_info (cpl_func,"Force VIS2ERR to %.e for all observations", err);

      int npol = gravi_pfits_get_pola_num (header, GRAVI_SC);
      for (int pol = 0; pol < npol; pol++) {
          cpl_table *  oi_table = gravi_data_get_oi_vis2 (oi_data, GRAVI_SC, pol, npol);
          cpl_array ** array = cpl_table_get_data_array (oi_table, "VIS2ERR");
          for (cpl_size row = 0; row<cpl_table_get_nrow (oi_table); row++) {
              cpl_array_fill_window (array[row], 0, CPL_SIZE_MAX, err);
          }
      }
  }
    
  gravi_msg_function_exit(1);
  return CPL_ERROR_NONE;
}


/*----------------------------------------------------------------------------*/
/**
 * @brief Coadd the observations together.
 * 
 * @param oi_data: the OIFITS map to co-add in-place
 * 
 * All observation are averaged together.  No checks of consistency !
 * The result always contains 6 VIS, 6 VIS2, 4 T3, and 4 FLUX,
 * whatever the initial content of the OIFITS. The MJD, UCOORD, VCOORD...
 * are averaged with equal weight (no weighting by the SNR).
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_average_vis (gravi_data * oi_data)
{
  gravi_msg_function_start(1);
  cpl_ensure_code (oi_data, CPL_ERROR_NULL_INPUT);

  gravi_msg_warning ("FIXME", "Average the different observation is EXPERIMENTAL");
  
  cpl_msg_warning (cpl_func, "FIXME: Weightening of UVCOORD and MJD is not done properly yet !");
  cpl_msg_warning (cpl_func, "FIXME: Integration of INT_TIME is not done properly yet !");

  int nbase = 6;
  cpl_table * oi_table;

  /* Create output data */
  cpl_propertylist * header = gravi_data_get_header (oi_data);

  /* Loop on oidata, type_data and polarisation */
  for (int type_data = 0; type_data < 2 ; type_data ++) {

    if (!gravi_data_has_type (oi_data, type_data == GRAVI_SC ? "_SC" : "_FT")) {
        cpl_msg_info (cpl_func, "OI_VIS has no %s, skip", GRAVI_TYPE(type_data));
        continue;
    }
    
	int npol = gravi_pfits_get_pola_num (header, type_data);
	for (int pol = 0 ; pol < npol ; pol++ ) {

      cpl_size nrow = cpl_table_get_nrow (gravi_data_get_oi_vis (oi_data, type_data, pol, npol))/nbase;
      if (nrow == 1) {
          cpl_msg_info (cpl_func, "OI_VIS %s has only one observation, skip", GRAVI_TYPE(type_data));
          continue;
      }

	  /* OI_VIS */
	  oi_table = gravi_data_get_oi_vis (oi_data, type_data, pol, npol);
	  gravi_vis_average_value (oi_table, "TIME", NULL, 6);
	  gravi_vis_average_value (oi_table, "MJD", NULL, 6);
	  gravi_vis_average_value (oi_table, "INT_TIME", NULL, 6);
	  gravi_vis_average_value (oi_table, "UCOORD", NULL, 6);
	  gravi_vis_average_value (oi_table, "VCOORD", NULL, 6);
	  gravi_vis_average_amp (oi_table, "VISAMP",  "VISAMPERR", 6);
	  gravi_vis_average_phi (oi_table, "VISPHI",  "VISPHIERR", 6);
	  gravi_vis_average_amp (oi_table, "RVIS",  "RVISERR", 6);
	  gravi_vis_average_amp (oi_table, "IVIS",  "IVISERR", 6);
	  // gravi_vis_average_amp (oi_table, "VISDATA",  "VISERR"); // FIXME: to be done !!
	  gravi_msg_warning ("FIXME", "VISDATA are not averaged !!!!");
	  cpl_table_erase_window (oi_table, 6, CPL_SIZE_MAX);
	  
	  CPLCHECK_MSG ("Cannot co-add OI_VIS");
	  
	  /* OI_VIS2 */
	  oi_table = gravi_data_get_oi_vis2 (oi_data, type_data, pol, npol);
	  gravi_vis_average_value (oi_table, "TIME", NULL, 6);
	  gravi_vis_average_value (oi_table, "MJD", NULL, 6);
	  gravi_vis_average_value (oi_table, "INT_TIME", NULL, 6);
	  gravi_vis_average_value (oi_table, "UCOORD", NULL, 6);
	  gravi_vis_average_value (oi_table, "VCOORD", NULL, 6);
	  gravi_vis_average_amp (oi_table, "VIS2DATA", "VIS2ERR", 6);
	  cpl_table_erase_window (oi_table, 6, CPL_SIZE_MAX);
	  
	  CPLCHECK_MSG ("Cannot co-add OI_VIS2");
	  
	  /* OI_FLUX */
	  oi_table = gravi_data_get_oi_flux (oi_data, type_data, pol, npol);
	  gravi_vis_average_value (oi_table, "TIME", NULL, 4);
	  gravi_vis_average_value (oi_table, "MJD", NULL, 4);
	  gravi_vis_average_value (oi_table, "INT_TIME", NULL, 4);
	  gravi_vis_average_amp (oi_table, "FLUX", "FLUXERR", 4);
	  cpl_table_erase_window (oi_table, 4, CPL_SIZE_MAX);
	  
	  CPLCHECK_MSG ("Cannot co-add OI_FLUX");
	  
	  /* OI_T3 */
	  oi_table = gravi_data_get_oi_t3 (oi_data, type_data, pol, npol);
	  gravi_vis_average_value (oi_table, "TIME", NULL, 4);
	  gravi_vis_average_value (oi_table, "MJD", NULL, 4);
	  gravi_vis_average_value (oi_table, "INT_TIME", NULL, 4);
	  gravi_vis_average_value (oi_table, "U1COORD", NULL, 4);
	  gravi_vis_average_value (oi_table, "V1COORD", NULL, 4);
	  gravi_vis_average_value (oi_table, "U2COORD", NULL, 4);
	  gravi_vis_average_value (oi_table, "V2COORD", NULL, 4);
	  gravi_vis_average_amp (oi_table, "T3AMP", "T3AMPERR", 4);
	  gravi_vis_average_phi (oi_table, "T3PHI", "T3PHIERR", 4);
	  cpl_table_erase_window (oi_table, 4, CPL_SIZE_MAX);
	  
	  CPLCHECK_MSG ("Cannot co-add OI_T3");
	  
	} /* End loop on polarisation */
  } /* End loop on type_data */

  gravi_msg_function_exit(1);
  return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Smooth amplitude column of OIFITS table.
 * 
 * @param oi_table:   the table to update (column depth will be modified)
 * @param name:       name of column to rebin
 * @param err:        corresponding error column
 * @param nsamp:      number of consecutive samples to bin into single bin
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_vis_smooth_amp (cpl_table * oi_table, const char * name, const char * err,
									 cpl_size nsamp)
{
  gravi_msg_function_start(1);
  cpl_ensure_code (oi_table, CPL_ERROR_NULL_INPUT);
  if (nsamp < 1) return CPL_ERROR_NONE;
  int nv;

  /* Get values */
  cpl_size nwave = cpl_table_get_column_depth (oi_table, name);
  cpl_size nrow = cpl_table_get_nrow (oi_table);
  cpl_ensure_code (nrow > 0, CPL_ERROR_ILLEGAL_INPUT);

  /* Get arrays */
  cpl_array ** v_array = cpl_table_get_data_array (oi_table, name);
  cpl_array ** e_array = cpl_table_get_data_array (oi_table, err);
  cpl_array ** f_array = cpl_table_get_data_array (oi_table, "FLAG");
  CPLCHECK_MSG ("Cannot get data");

  /* Allocate output */
  cpl_array * smo_array = cpl_array_duplicate (v_array[0]);
  cpl_array * err_array = cpl_array_duplicate (e_array[0]);
  
  /* Loop on rows */
  for (cpl_size row = 0 ; row < nrow ; row ++) {

    /* Median filter the uncertainties, to avoid
     * putting all on some sample */
    cpl_vector * i_vector = cpl_vector_new (nwave);
    for (cpl_size wave = 0; wave < nwave; wave++) 
        cpl_vector_set (i_vector, wave, cpl_array_get (e_array[row],wave,&nv));
    cpl_vector * o_vector;
    o_vector = cpl_vector_filter_median_create (i_vector, nsamp);

	/* Loop on  waves */
	for (cpl_size wave = 0 ; wave < nwave ; wave ++) {
        double sum = 0.0, weight = 0.0;

        /* Loop on samples to average */
		for (cpl_size samp = CPL_MAX(0,wave-nsamp) ; samp < CPL_MIN(nwave,wave+nsamp) ; samp ++) {
            if (cpl_array_get (f_array[row],samp,&nv)) {
                weight += 10e-20;
                sum += 0.0;
            } else {
                double w = pow (cpl_vector_get (o_vector,samp), -2);
                sum    += cpl_array_get (v_array[row],samp,&nv) * w;
                weight += w;
            }
		}

        cpl_array_set_double (smo_array, wave, sum / weight);
        cpl_array_set_double (err_array, wave, pow (weight, -0.5));
    }

    /* Set back */
    cpl_table_set_array (oi_table, name, row, smo_array);
    cpl_table_set_array (oi_table, err,  row, err_array);
    CPLCHECK_MSG ("Cannot smooth amp");
	
    FREE (cpl_vector_delete, i_vector);
    FREE (cpl_vector_delete, o_vector);
  } /* End loop on rows */

  FREE (cpl_array_delete, smo_array);
  FREE (cpl_array_delete, err_array);

  gravi_msg_function_exit(1);
  return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Smooth phase column of OIFITS table.
 * 
 * @param oi_table:   the table to update (column depth will be modified)
 * @param name:       name of column to rebin
 * @param err:        corresponding error column
 * @param nsamp:      number of consecutive samples to bin into single bin
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_vis_smooth_phi (cpl_table * oi_table, const char * name, const char * err,
									 cpl_size nsamp)
{
  gravi_msg_function_start(1);
  cpl_ensure_code (oi_table, CPL_ERROR_NULL_INPUT);
  if (nsamp < 1) return CPL_ERROR_NONE;

  /* Get values */
  cpl_size nwave = cpl_table_get_column_depth (oi_table, name);
  cpl_size nrow = cpl_table_get_nrow (oi_table);
  cpl_ensure_code (nrow > 0, CPL_ERROR_ILLEGAL_INPUT);
  int nv;

  /* Get arrays */
  cpl_array ** v_array = cpl_table_get_data_array (oi_table, name);
  cpl_array ** e_array = cpl_table_get_data_array (oi_table, err);
  cpl_array ** f_array = cpl_table_get_data_array (oi_table, "FLAG");
  CPLCHECK_MSG ("Cannot get data");
  
  /* Allocate output */
  cpl_array * smo_array = cpl_array_duplicate (v_array[0]);
  cpl_array * err_array = cpl_array_duplicate (e_array[0]);

  /* Loop on rows */
  for (cpl_size row = 0 ; row < nrow ; row ++) {

    /* Median filter the uncertainties, to avoid
     * putting all on some sample */
    cpl_vector * i_vector = cpl_vector_new (nwave);
    for (cpl_size wave = 0; wave < nwave; wave++) 
        cpl_vector_set (i_vector, wave, cpl_array_get (e_array[row],wave,&nv));
    cpl_vector * o_vector;
    o_vector = cpl_vector_filter_median_create (i_vector, nsamp);

	/* Loop on  waves */
	for (cpl_size wave = 0 ; wave < nwave ; wave ++) {
        double complex sum = 0.0 + I*0.0;
        double weight = 0.0;

        /* Loop on samples to average */
		for (cpl_size samp = CPL_MAX(0,wave-nsamp) ; samp < CPL_MIN(nwave,wave+nsamp) ; samp ++) {
            if (cpl_array_get (f_array[row],samp,&nv)) {
                weight += 10e-20;
                sum += 0.0;
            } else {
                double w = pow (cpl_vector_get (o_vector,samp), -2);
                sum    += cexp (1.*I* cpl_array_get (v_array[row],samp,&nv) * CPL_MATH_RAD_DEG) * w;
                weight += w;
            }
		}
        
        cpl_array_set_double (smo_array, wave, carg (sum) * CPL_MATH_DEG_RAD);
        cpl_array_set_double (err_array, wave, pow (weight, -0.5));
    }

    /* Set back */
    cpl_table_set_array (oi_table, name, row, smo_array);
    cpl_table_set_array (oi_table, err,  row, err_array);
    CPLCHECK_MSG ("Cannot smooth phi");
	
    FREE (cpl_vector_delete, i_vector);
    FREE (cpl_vector_delete, o_vector);
  } /* End loop on rows */

  FREE (cpl_array_delete, smo_array);
  FREE (cpl_array_delete, err_array);

  gravi_msg_function_exit(1);
  return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Smooth amp column of OIFITS table.
 * 
 * @param oi_table:   the table to update (column depth will be modified)
 * @param name:       name of column to rebin
 * @param err:        corresponding error column
 * @param nsamp:      number of consecutive samples to bin into single bin
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_vis_fit_amp (cpl_table * oi_table, const char * name,
                                  const char * err, cpl_size maxdeg)
{
  gravi_msg_function_start(1);
  cpl_ensure_code (oi_table, CPL_ERROR_NULL_INPUT);
  if (maxdeg < 0) return CPL_ERROR_NONE;

  /* Get values */
  cpl_size nwave = cpl_table_get_column_depth (oi_table, name);
  cpl_size nrow = cpl_table_get_nrow (oi_table);
  cpl_ensure_code (nrow > 0, CPL_ERROR_ILLEGAL_INPUT);
  int nv;

  /* Get arrays */
  cpl_array ** v_array = cpl_table_get_data_array (oi_table, name);
  cpl_array ** e_array = cpl_table_get_data_array (oi_table, err);
  cpl_array ** f_array = cpl_table_get_data_array (oi_table, "FLAG");
  CPLCHECK_MSG ("Cannot get data");
  
  /* Loop on rows */
  for (cpl_size row = 0 ; row < nrow ; row ++) {

      /* Create the vectors and matrix */
      cpl_matrix * coeff = cpl_matrix_new (nwave,maxdeg+1);
      cpl_matrix * rhs   = cpl_matrix_new (nwave,1);
      
      /* Fill */
      for (cpl_size wave = 0 ; wave < nwave ; wave ++) {
		  double weight = cpl_array_get (f_array[row],wave,&nv) ? 10e-20 : 
                          pow (cpl_array_get (e_array[row],wave,&nv), -2);
		  double value = cpl_array_get (f_array[row],wave,&nv) ? 0.0 : 
                         cpl_array_get (v_array[row],wave,&nv);

                  if (weight > 1e10)
                  {
                    cpl_msg_warning (cpl_func, "name = %s row = %lli wave = %lli "
                                     "has faulty uncertainty", name, row, wave);
                    weight = 0.0;
                    cpl_array_set (f_array[row],wave,1);
                  }
                  
          cpl_matrix_set (rhs, wave, 0, value * weight);
          for (cpl_size deg = 0; deg <= maxdeg; deg++) 
              cpl_matrix_set (coeff, wave, deg, pow ((double)wave,(double)deg) * weight);
          CPLCHECK_MSG ("Cannot fill");
      }

      //printf ("name = %s row = %i maxdeg = %i\n", name, row, maxdeg);
      //cpl_matrix_dump (rhs, NULL);
      //cpl_matrix_dump (coeff, NULL);

      /* Solve */
      cpl_errorstate prev_state = cpl_errorstate_get();
      cpl_matrix * solve = cpl_matrix_solve_normal (coeff, rhs);

      /* Dump errors */
      if ( !cpl_errorstate_is_equal (prev_state))
      {
          cpl_errorstate_dump (prev_state, 0, NULL);
          cpl_msg_error (cpl_func,"%s row=%lld",name,row);
      }
      CPLCHECK_MSG ("Cannot solve matrix");
      
      /* Evaluate */
      for (cpl_size wave = 0 ; wave < nwave ; wave ++) {
          double value = 0;
          for (cpl_size deg = 0; deg <= maxdeg; deg++)
              value += cpl_matrix_get (solve, deg, 0) * pow (wave, deg);
          cpl_array_set (v_array[row],wave,value);
          CPLCHECK_MSG ("Cannot evaluate");
      }
      
      FREE (cpl_matrix_delete, coeff);
      FREE (cpl_matrix_delete, rhs);
      FREE (cpl_matrix_delete, solve);
  }

  gravi_msg_function_exit(1);
  return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Smooth the SC table by nsamp consecutive spectral bins.
 * 
 * @param oi_data      VIS data to process, in-place
 * @param nsamp_vis    integer, the number of consecutive bin to smooth
 * @param nsamp_flx    integer, the number of consecutive bin to smooth
 * @param maxdeg       integer, fit order
 * @param mindeg       integer, fit order    
 * 
 * The OI_VIS, OI_VIS2, OI_FLUX, OI_T3
 * tables are updated accordingly. Note that this operation is not
 * flux conservative (FIXME: understand how to better averaged flux).
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_vis_smooth (gravi_data * oi_data,
                                 cpl_size nsamp_vis,
                                 cpl_size nsamp_flx,
                                 cpl_size maxdeg)
{
  gravi_msg_function_start(1);
  cpl_ensure_code (oi_data, CPL_ERROR_NULL_INPUT);
  
  cpl_table * oi_table;

  /* Create output data */
  cpl_propertylist * header = gravi_data_get_header (oi_data);

  int type_data = GRAVI_SC;
  int npol = gravi_pfits_get_pola_num (header, type_data);
  
  for (int pol = 0 ; pol < npol ; pol++ ) {

	/* OI_FLUX */
	oi_table = gravi_data_get_oi_flux (oi_data, type_data, pol, npol);
        gravi_vis_flag_nan (oi_table);
	gravi_vis_flag_lower (oi_table, "FLUXERR", "FLAG", 0.0);
	gravi_vis_smooth_amp (oi_table, "FLUX", "FLUXERR", nsamp_flx);
	gravi_vis_flag_relative_threshold (oi_table, "FLUXERR", "FLUX", "FLAG", 1.0);
	CPLCHECK_MSG ("Cannot resamp OI_FLUX");

	/* OI_VIS2 */
	oi_table = gravi_data_get_oi_vis2 (oi_data, type_data, pol, npol);
        gravi_vis_flag_nan (oi_table);
	gravi_vis_flag_lower (oi_table, "VIS2ERR", "FLAG", 0.);
        gravi_vis_flag_median (oi_table, "VIS2ERR", "FLAG", 5.0);    
	gravi_vis_smooth_amp (oi_table, "VIS2DATA", "VIS2ERR", nsamp_vis);
        gravi_vis_fit_amp (oi_table, "VIS2DATA", "VIS2ERR", maxdeg);
	gravi_vis_flag_threshold (oi_table, "VIS2ERR", "FLAG", 1.);
	CPLCHECK_MSG ("Cannot resamp OI_VIS2");

	/* OI_VIS */
	oi_table = gravi_data_get_oi_vis (oi_data, type_data, pol, npol);
        gravi_vis_flag_nan (oi_table);
	gravi_vis_flag_lower (oi_table, "VISAMPERR", "FLAG", 0.);
        gravi_vis_flag_median (oi_table, "VISPHIERR", "FLAG", 5.0);
	gravi_vis_smooth_amp (oi_table, "VISAMP", "VISAMPERR", nsamp_vis);
	gravi_vis_smooth_phi (oi_table, "VISPHI", "VISPHIERR", nsamp_vis);
        gravi_vis_fit_amp (oi_table, "VISAMP", "VISAMPERR", maxdeg);
        gravi_vis_fit_amp (oi_table, "VISPHI", "VISPHIERR", maxdeg);
	gravi_vis_smooth_amp (oi_table, "RVIS", "RVISERR", nsamp_flx);
	gravi_vis_smooth_amp (oi_table, "IVIS", "IVISERR", nsamp_flx);
	gravi_vis_flag_threshold (oi_table, "VISAMPERR", "FLAG", 1.);
    
        gravi_msg_warning ("FIXME", "VISDATA is not properly smooth !!");
	CPLCHECK_MSG ("Cannot resamp OI_VIS");
	
	/* OI_T3 */
	oi_table = gravi_data_get_oi_t3 (oi_data, type_data, pol, npol);
        gravi_vis_flag_nan (oi_table);
	gravi_vis_flag_lower (oi_table, "T3AMPERR", "FLAG", 0.0);
        gravi_vis_flag_median (oi_table, "T3PHIERR", "FLAG", 5.0);
	gravi_vis_smooth_amp (oi_table, "T3AMP", "T3AMPERR", nsamp_vis);
	gravi_vis_smooth_phi (oi_table, "T3PHI", "T3PHIERR", nsamp_vis);
	gravi_vis_fit_amp (oi_table, "T3AMP", "T3AMPERR", maxdeg);
	gravi_vis_fit_amp (oi_table, "T3PHI", "T3PHIERR", maxdeg);
	gravi_vis_flag_threshold (oi_table, "T3AMPERR", "FLAG", 1.0);
	CPLCHECK_MSG ("Cannot resamp OI_T3");
	  
  } /* End loop on polarisation */

  gravi_msg_function_exit(1);
  return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Rebin amplitude column of OIFITS table.
 * 
 * @param oi_table:   the table to update (column depth will be modified)
 * @param name:       name of column to rebin
 * @param err:        corresponding error column
 * @param nsamp:      number of consecutive samples to bin into single bin
 * @param nwave_new:  number of final bins to reach (FIXME: why needed ??)
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_vis_resamp_amp (cpl_table * oi_table, const char * name, const char * err,
									 cpl_size nsamp, cpl_size nwave_new)
{
  gravi_msg_function_start(1);
  cpl_ensure_code (oi_table, CPL_ERROR_NULL_INPUT);
  
  /* Loop on rows */
  cpl_size nrow = cpl_table_get_nrow (oi_table);
  cpl_ensure_code (nrow > 0, CPL_ERROR_ILLEGAL_INPUT);
  for (cpl_size row = 0 ; row < nrow ; row ++) {

	/* Loop on new waves */
	for (cpl_size wave = 0 ; wave < nwave_new ; wave ++) {
	  double sum = 0.0;
		double weight = 0.0;
		for (cpl_size samp = 0 ; samp < nsamp ; samp ++) {
		  double w = pow (gravi_table_get_value (oi_table,err,row,wave*nsamp+samp), -2.0);
		  if (gravi_table_get_value (oi_table,"FLAG",row,wave*nsamp+samp)) w = 10e-20;
		  sum    += gravi_table_get_value (oi_table,name,row,wave*nsamp+samp) * w;
		  weight += w;
		}
		gravi_table_set_value (oi_table,name,row,wave, sum / weight);
		gravi_table_set_value (oi_table,err,row,wave, pow (weight, -0.5));
	  }
	
  } /* End loop on rows */

  cpl_table_set_column_depth (oi_table, name, nwave_new);
  cpl_table_set_column_depth (oi_table, err,  nwave_new);

  gravi_msg_function_exit(1);
  return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Rebin phase column of OIFITS table (arg{<exp(i.phi)>})
 * 
 * @param oi_table:   the table to update (column depth will be modified)
 * @param name:       name of column to rebin
 * @param err:        corresponding error column
 * @param nsamp:      number of consecutive samples to bin into single bin
 * @param nwave_new:  number of final bins to reach (FIXME: why needed ??)
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_vis_resamp_phi (cpl_table * oi_table, const char * name, const char * err,
									 cpl_size nsamp, cpl_size nwave_new)
{
  gravi_msg_function_start(1);
  cpl_ensure_code (oi_table, CPL_ERROR_NULL_INPUT);
  
  /* Loop on rows */
  cpl_size nrow = cpl_table_get_nrow (oi_table);
  cpl_ensure_code (nrow > 0, CPL_ERROR_ILLEGAL_INPUT);
  for (cpl_size row = 0 ; row < nrow ; row ++) {

	/* Loop on new waves */
	for (cpl_size wave = 0 ; wave < nwave_new ; wave ++) {
	    double complex sum = 0.0;
		double weight = 0.0;
		for (cpl_size samp = 0 ; samp < nsamp ; samp ++) {
		  double w = pow (gravi_table_get_value (oi_table,err,row,wave*nsamp+samp), -2.0);
		  if (gravi_table_get_value (oi_table,"FLAG",row,wave*nsamp+samp)) w = 10e-20;
		  sum    += cexp (1.*I* gravi_table_get_value (oi_table,name,row,wave*nsamp+samp) * CPL_MATH_RAD_DEG) * w;
		  weight += w;
		}
		gravi_table_set_value (oi_table,name,row,wave, carg (sum) * CPL_MATH_DEG_RAD);
		gravi_table_set_value (oi_table,err,row,wave, pow (weight, -0.5));
	  }
	
  } /* End loop on rows */

  cpl_table_set_column_depth (oi_table, name, nwave_new);
  cpl_table_set_column_depth (oi_table, err,  nwave_new);

  gravi_msg_function_exit(1);
  return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Re-bin the SC table by nsamp consecutive spectral bins.
 * 
 * @param oi_data    VIS data to process, in-place
 * @param nsamp:       integer, the number of consecutive bin to sum
 * 
 * The OI_VIS, OI_VIS2, OI_FLUX, OI_T3 and OI_WAVELENGTHs
 * tables are updated accordingly. Note that this operation is not
 * flux conservative (FIXME: understand how to better averaged flux).
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_vis_resamp (gravi_data * oi_data, cpl_size nsamp)
{
  gravi_msg_function_start(1);
  cpl_ensure_code (oi_data, CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (nsamp>1,   CPL_ERROR_ILLEGAL_INPUT);
  
  cpl_table * oi_table;
  cpl_size nwave, nwave_new;
  int nv = 0;

  /* Create output data */
  cpl_propertylist * header = gravi_data_get_header (oi_data);

  int type_data = GRAVI_SC;
  int npol = gravi_pfits_get_pola_num (header, type_data);
  for (int pol = 0 ; pol < npol ; pol++ ) {

	/* OI_WAVELENGTH table */
	oi_table = gravi_data_get_oi_wave (oi_data, type_data, pol, npol);
	nwave = cpl_table_get_nrow (oi_table);

	/* New number of wave sample */
	nwave_new = nwave / nsamp;
    cpl_msg_info (cpl_func, "Resamp the SC data by %lld bins: %lld -> %lld",
                  nsamp, nwave, nwave_new);
	cpl_ensure_code (nwave_new > 1, CPL_ERROR_ILLEGAL_INPUT);
    
	/* Get the wavelength */
	for (cpl_size wave = 0 ; wave < nwave_new ; wave ++) {
	  /* Compute effective band first */
	  cpl_table_set (oi_table, "EFF_BAND", wave,
					 cpl_table_get (oi_table, "EFF_WAVE", wave*nsamp+nsamp-1, &nv) -
					 cpl_table_get (oi_table, "EFF_WAVE", wave*nsamp, &nv));
	  
	  /* Average them */
	  double mean = 0.0;
	  for (cpl_size samp = 0 ; samp < nsamp ; samp++)
		mean += cpl_table_get (oi_table, "EFF_WAVE", wave*nsamp+samp, &nv);
	  cpl_table_set (oi_table, "EFF_WAVE", wave, mean / nsamp);
	}

	/* Erase last entries */
	cpl_table_erase_window (oi_table, nwave_new, CPL_SIZE_MAX);
	CPLCHECK_MSG ("Cannot resamp OI_WAVELENGTH");

	/* OI_FLUX */
	oi_table = gravi_data_get_oi_flux (oi_data, type_data, pol, npol);
	gravi_vis_resamp_amp (oi_table, "FLUX", "FLUXERR", nsamp, nwave_new);
	cpl_table_set_column_depth (oi_table, "FLAG", nwave_new);
	gravi_vis_flag_relative_threshold (oi_table, "FLUXERR", "FLUX", "FLAG", 1.0);
	CPLCHECK_MSG ("Cannot resamp OI_FLUX");

	/* OI_VIS2 */
	oi_table = gravi_data_get_oi_vis2 (oi_data, type_data, pol, npol);
	gravi_vis_resamp_amp (oi_table, "VIS2DATA", "VIS2ERR", nsamp, nwave_new);
	cpl_table_set_column_depth (oi_table, "FLAG", nwave_new);
	gravi_vis_flag_threshold (oi_table, "VIS2ERR", "FLAG", 1.);
	CPLCHECK_MSG ("Cannot resamp OI_VIS2");

	/* OI_VIS */
	oi_table = gravi_data_get_oi_vis (oi_data, type_data, pol, npol);
	gravi_vis_resamp_amp (oi_table, "VISAMP", "VISAMPERR", nsamp, nwave_new);
	gravi_vis_resamp_amp (oi_table, "VISPHI", "VISPHIERR", nsamp, nwave_new);
	gravi_vis_resamp_amp (oi_table, "RVIS", "RVISERR", nsamp, nwave_new);
	gravi_vis_resamp_amp (oi_table, "IVIS", "IVISERR", nsamp, nwave_new);
	cpl_table_set_column_depth (oi_table, "FLAG", nwave_new);
	gravi_vis_flag_threshold (oi_table, "VISAMPERR", "FLAG", 1.);
    
    gravi_msg_warning ("FIXME", "VISDATA is not properly resampled !!");
	cpl_table_set_column_depth (oi_table, "VISDATA", nwave_new);
	cpl_table_set_column_depth (oi_table, "VISERR", nwave_new);
	CPLCHECK_MSG ("Cannot resamp OI_VIS");
	
	/* OI_T3 */
	oi_table = gravi_data_get_oi_t3 (oi_data, type_data, pol, npol);
	gravi_vis_resamp_amp (oi_table, "T3AMP", "T3AMPERR", nsamp, nwave_new);
	gravi_vis_resamp_amp (oi_table, "T3PHI", "T3PHIERR", nsamp, nwave_new);
	cpl_table_set_column_depth (oi_table, "FLAG", nwave_new);
	gravi_vis_flag_threshold (oi_table, "T3AMPERR", "FLAG", 1.0);
	CPLCHECK_MSG ("Cannot resamp OI_T3");
	  
  } /* End loop on polarisation */

  gravi_msg_function_exit(1);
  return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Duplicate the column FLUX into FLUXDATA, for OIFITS2 compliance
 * 
 * @param oi_data    VIS data to process, in-place
 * 
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_vis_copy_fluxdata (gravi_data * oi_data)
{
  gravi_msg_function_start(1);
  cpl_ensure_code (oi_data, CPL_ERROR_NULL_INPUT);
  
  /* header */
  cpl_propertylist * header = gravi_data_get_header (oi_data);

  /* Loop on oidata, type_data and polarisation */
  for (int type_data = 0; type_data < 2 ; type_data ++) {

    if (!gravi_data_has_type (oi_data, type_data == GRAVI_SC ? "_SC" : "_FT")) {
        cpl_msg_info (cpl_func, "OI_FLUX has no %s, skip", GRAVI_TYPE(type_data));
        continue;
    }
    
    int npol = gravi_pfits_get_pola_num (header, type_data);
    for (int pol = 0 ; pol < npol ; pol++ ) {

      /* OI_FLUX table */
      cpl_table * oi_flux;
      oi_flux = gravi_data_get_oi_flux (oi_data, type_data, pol, npol);

      /* Delete column if existing */
      if (cpl_table_has_column (oi_flux, "FLUXDATA") )
        cpl_table_erase_column (oi_flux, "FLUXDATA");

      /* Duplicate FLUX into FLUXDATA */
      cpl_table_duplicate_column (oi_flux, "FLUXDATA", oi_flux, "FLUX");
	  
    } /* End loop on polarisation */

  } /* End loop on SC/FT */
  
  gravi_msg_function_exit(1);
  return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Flag samples of OIFITS table which are NAN or NULL
 * 
 * @param oi_table:  the cpl_table to update, in-place
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_vis_flag_nan (cpl_table * oi_table)
{
    gravi_msg_function_start(1);
    cpl_ensure_code (oi_table, CPL_ERROR_NULL_INPUT);

    cpl_size nrow = cpl_table_get_nrow (oi_table);
    cpl_ensure_code (nrow > 0, CPL_ERROR_ILLEGAL_INPUT);
    
    int ncols = 10;
    const char * names[] = {"VIS2DATA","VIS2ERR","VISAMP","VISAMPERR",
                            "VISPHI","VISPHIERR","T3PHI","T3PHIERR",
                            "T3AMP","T3AMPERR"};

    /* Loop on columns */
    for (int c = 0; c < ncols; c++) {
        
        if (!cpl_table_has_column (oi_table,names[c])) continue;
        cpl_msg_info (cpl_func,"Check column %s",names[c]);

        /* Get data of this columns */
        cpl_size nwave = cpl_table_get_column_depth (oi_table, names[c]);
        cpl_array ** v_array = cpl_table_get_data_array (oi_table, names[c]);
        cpl_array ** f_array = cpl_table_get_data_array (oi_table, "FLAG");
        CPLCHECK_MSG ("Cannot get data");

        /* Loop on rows and waves */
        cpl_size ninvalid = 0;
        for (cpl_size row = 0; row < nrow ; row ++) {
            for (cpl_size wave = 0 ; wave < nwave ; wave ++) {

                /* Get value */
                int nv = 0;
                double value = cpl_array_get (v_array[row], wave, &nv);

                /* Check value */
                if (nv || isnan (value)) {
                    cpl_array_set (f_array[row], wave, 1.0);
                    cpl_array_set (v_array[row], wave, 0.0);
                    ninvalid ++;
                }
                CPLCHECK_MSG ("Cannot check data");
            }
        } /* End loop on rows and waves */

        /* Verbose */
        if (ninvalid > 0) {
            cpl_msg_warning (cpl_func, "Flag %lld invalid data (NAN or NULL) in %s", ninvalid, names[c]);
        }
        
    } /* End loop on columns */
     
    gravi_msg_function_exit(1);
    return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Flag samples of OIFITS table based on absolute threshold.
 * 
 * @param oi_table:  the cpl_table to update, in-place
 * @param name:      the column name to verify
 * @param flag:      the corresponding flag column (to update in-place)
 * @param value:     the threshold to compare
 *
 * For each sample in the column, if sample>value then the corresponding 
 * element in the flag column is set to 'T'.
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_vis_flag_threshold (cpl_table * oi_table, const char * data, const char *flag, double value)
{
  gravi_msg_function_start(0);
  cpl_ensure_code (oi_table, CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (data, CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (flag, CPL_ERROR_ILLEGAL_OUTPUT);
  
  /* Get pointer to speed up */
  int nv = 0;
  cpl_size nrow = cpl_table_get_nrow (oi_table);
  cpl_array ** pdata = cpl_table_get_data_array (oi_table, data);
  cpl_array ** pflag = cpl_table_get_data_array (oi_table, flag);
  
  CPLCHECK_MSG ("Cannot get data");
  
  cpl_size size = cpl_array_get_size (pdata[0]);
  
  /* Loop on row and index. Add to FLAG if data is above threshold */
  for ( cpl_size row = 0 ; row < nrow ; row ++ ) {
        if (pdata[row]==NULL) continue;
        
	for ( cpl_size indx = 0 ; indx < size ; indx ++ ) {
	  if ( cpl_array_get (pdata[row], indx, &nv) > value ) {
	       cpl_array_set (pflag[row], indx, cpl_array_get (pflag[row], indx, &nv) + 1 );
	  }
	}
  }

  gravi_msg_function_exit(0);
  return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Flag samples of OIFITS table based on absolute threshold.
 * 
 * @param oi_table:  the cpl_table to update, in-place
 * @param name:      the column name to verify
 * @param flag:      the corresponding flag column (to update in-place)
 * @param value:     the threshold to compare
 *
 * For each sample in the column, if sample<=value then the corresponding 
 * element in the flag column is set to 'T'.
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_vis_flag_lower (cpl_table * oi_table, const char * data, const char *flag, double value)
{
  gravi_msg_function_start(0);
  cpl_ensure_code (oi_table, CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (data, CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (flag, CPL_ERROR_ILLEGAL_OUTPUT);
  
  /* Get pointer to speed up */
  int nv = 0;
  cpl_size nrow = cpl_table_get_nrow (oi_table);
  cpl_array ** pdata = cpl_table_get_data_array (oi_table, data);
  cpl_array ** pflag = cpl_table_get_data_array (oi_table, flag);
  
  CPLCHECK_MSG ("Cannot get data");
  
  cpl_size size = cpl_array_get_size (pdata[0]);
  
  /* Loop on row and index. Add to FLAG if data is above threshold */
  for ( cpl_size row = 0 ; row < nrow ; row ++ ) {
        if (pdata[row]==NULL) continue;
        
	for ( cpl_size indx = 0 ; indx < size ; indx ++ ) {
	  if ( cpl_array_get (pdata[row], indx, &nv) <= value ) {
	       cpl_array_set (pflag[row], indx, cpl_array_get (pflag[row], indx, &nv) + 1 );
	  }
	}
  }

  gravi_msg_function_exit(0);
  return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Flag samples of OIFITS table based on runnning median.
 * 
 * @param oi_table:  the cpl_table to update, in-place
 * @param name:      the column name to verify
 * @param flag:      the corresponding flag column (to update in-place)
 * @param value:     the threshold to compare
 *
 * For each sample in the column, if sample>value*median then the corresponding 
 * element in the flag column is set to 'T'.
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_vis_flag_median (cpl_table * oi_table, const char * data, const char *flag, double value)
{
  gravi_msg_function_start(0);
  cpl_ensure_code (oi_table, CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (data, CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (flag, CPL_ERROR_ILLEGAL_OUTPUT);
  
  /* Get pointer to speed up */
  int nv = 0;
  cpl_size nrow = cpl_table_get_nrow (oi_table);
  cpl_array ** pdata = cpl_table_get_data_array (oi_table, data);
  cpl_array ** pflag = cpl_table_get_data_array (oi_table, flag);
  
  CPLCHECK_MSG ("Cannot get data");
  
  cpl_size size = cpl_array_get_size (pdata[0]);
  cpl_vector * i_vector = cpl_vector_new (size);
  
  /* Loop on row and index. Add to FLAG if data is above threshold */
  for ( cpl_size row = 0 ; row < nrow ; row ++ ) {
        if (pdata[row]==NULL || size<100) continue;

      /* Set */
      for (cpl_size indx = 0; indx < size; indx++) 
          cpl_vector_set (i_vector, indx, cpl_array_get (pdata[0],indx,&nv));

      /* Median filter over 8 pixels wide */
      cpl_vector * o_vector = cpl_vector_filter_median_create (i_vector, 4);

      /* Check whose pixel have a large values compare to this median */
      for ( cpl_size indx = 0 ; indx < size ; indx ++ ) {
          if ( cpl_array_get (pdata[row], indx, &nv) > value * cpl_vector_get (o_vector, indx)) {
              cpl_array_set (pflag[row], indx, cpl_array_get (pflag[row], indx, &nv) + 1 );
          }
      }

      FREE (cpl_vector_delete, o_vector);
  }

  FREE (cpl_vector_delete, i_vector);
  gravi_msg_function_exit(0);
  return CPL_ERROR_NONE;
}


/*----------------------------------------------------------------------------*/
/**
 * @brief Flag samples of OIFITS table based on relative threshold.
 * 
 * @param oi_table:  the cpl_table to update, in-place
 * @param name:      the column name to verify
 * @param name:      the corresponding error column
 * @param flag:      the corresponding flag column (to update in-place)
 * @param value:     the threshold to compare
 *
 * For each sample in the column, if err/sample>value then the corresponding 
 * element in the flag column is set to 'T'.
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_vis_flag_relative_threshold(cpl_table * oi_table, const char * err,
												 const char * data, const char *flag, double value)
{
  gravi_msg_function_start(0);
  cpl_ensure_code (oi_table, CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (data, CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (err,  CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (flag, CPL_ERROR_ILLEGAL_OUTPUT);
  
  /* Get pointer to speed up */
  int nv = 0;
  cpl_size row, nrow = cpl_table_get_nrow (oi_table);
  cpl_array ** perr  = cpl_table_get_data_array (oi_table, err);
  cpl_array ** pdata = cpl_table_get_data_array (oi_table, data);
  cpl_array ** pflag = cpl_table_get_data_array (oi_table, flag);

  CPLCHECK_MSG ("Cannot get data");
  
  cpl_size indx, size = cpl_array_get_size (pdata[0]);
  
  /* Loop on row and index. Add to FLAG if data is above threshold */
  for ( row = 0 ; row < nrow ; row ++ ) {
        if (perr[row]==NULL) continue;
        cpl_array * tmp = cpl_array_duplicate (perr[row]);
        cpl_array_divide (tmp, pdata[row]);
	for ( indx = 0 ; indx < size ; indx ++ ) {
	  if ( cpl_array_get (tmp, indx, &nv) > value) {
	       cpl_array_set (pflag[row], indx, cpl_array_get (pflag[row], indx, &nv) + 1 );
	  }
	}
	cpl_array_delete (tmp);
  }

  gravi_msg_function_exit(0);
  return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Erase observation from an OIFITS table
 *
 * @param oi_table     The oitable, modified in-place
 * @param flag_array   Array with nobs value
 * @param ntel         Number of row in oitable for each obs
 * 
 * The recipe first select all obs in oitable corresponding to flagged
 * observations (flag_array != 0). Then they are erased in-place.
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_vis_erase_obs (cpl_table * oi_table, cpl_array *flag_array, cpl_size ntel)
{
    gravi_msg_function_start(0);
    cpl_ensure_code (oi_table,   CPL_ERROR_NULL_INPUT);
    cpl_ensure_code (flag_array, CPL_ERROR_NULL_INPUT);

    /* Get nrow */
    cpl_size nrow = cpl_table_get_nrow (oi_table) / ntel;
    cpl_ensure_code (nrow == cpl_array_get_size (flag_array),
                     CPL_ERROR_ILLEGAL_INPUT);

    /* Loop and select */
    cpl_table_unselect_all (oi_table);
    for (cpl_size row = 0; row < nrow; row++) {
        if (cpl_array_get (flag_array, row, NULL) == 0) continue;
        cpl_table_or_selected_window (oi_table, row * ntel, ntel);
    }

    /* Delete selected */
    cpl_table_erase_selected (oi_table);
    CPLCHECK_MSG ("Cannot erase");
    
    gravi_msg_function_exit(0);
    return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Compute the mean of a column in OIFITS table, and save the
 *        result in the specified output table.
 *
 * @param out_table    Output table
 * @param in_table     Input table
 * @param name         Column name
 * @param ntel         Number of tel (or base, or triplet)
 *
 * The routine create a column in the output table with the same name and units
 * as the one in the input table. It is filled with an average of the signal
 * over all DITs, considering the column REJECTION_FLAG if present (non-zero
 * means frame is rejected). This averaging is performed independently for the
 * ntel tel (or base or triplet). The output table shall thus contain ntel
 * rows while the input table shall contain ntel*NDIT rows.
 *
 * Note that this routine is not optimized for performance and thus shall 
 * not be used on the FT tables.
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_vis_compute_column_mean (cpl_table * out_table,
                                              cpl_table * in_table,
                                              const char * name, int ntel)
{
    gravi_msg_function_start(0);
    cpl_ensure_code (out_table, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code (in_table,  CPL_ERROR_NULL_INPUT);
    cpl_ensure_code (name,      CPL_ERROR_NULL_INPUT);
    cpl_ensure_code (ntel == cpl_table_get_nrow (out_table),
                     CPL_ERROR_ILLEGAL_OUTPUT);

    /* Check if column exist */
    if ( !cpl_table_has_column (in_table, name)) {
        cpl_msg_info (cpl_func, "Cannot average column %s (not existing)", name);
        return CPL_ERROR_NONE;
    }

    /* Cast the type into an int, to avoid warnings */
	int type = cpl_table_get_column_type (in_table, name);
    cpl_size depth = cpl_table_get_column_depth (in_table, name);

    /* Get the number of rows */
    cpl_size nrow  = cpl_table_get_nrow (in_table) / ntel;
    cpl_ensure_code (nrow,  CPL_ERROR_ILLEGAL_INPUT);

    /* Get the column REJECTION_FLAG */
    int * flag = NULL;
    if (cpl_table_has_column (in_table, "REJECTION_FLAG"))
        flag = cpl_table_get_data_int (in_table, "REJECTION_FLAG");

    cpl_msg_info (cpl_func, "Average column: %s (%s REJECTION_FLAG)",
                  name, flag ? "with" : "without");

	switch (type) {
	case CPL_TYPE_DOUBLE:
	case CPL_TYPE_FLOAT:
	case CPL_TYPE_INT:
        /* Case scalar column */
        if (!cpl_table_has_column (out_table, name))
            cpl_table_new_column (out_table, name, CPL_TYPE_DOUBLE);
        for (int tel = 0; tel < ntel; tel++) {
            cpl_size nvalid = 0;
            double mean = 0.0;
            for (cpl_size row = 0; row < nrow; row++) {
                if (flag && flag[row*ntel+tel] != 0) continue;
                nvalid ++;
                mean += cpl_table_get (in_table, name, row*ntel+tel, NULL);
            }
            if (nvalid > 0) mean /= nvalid;
            cpl_table_set (out_table, name, tel, mean);
        }
        break;
        
	case CPL_TYPE_POINTER|CPL_TYPE_DOUBLE:
	case CPL_TYPE_POINTER|CPL_TYPE_FLOAT:
	case CPL_TYPE_POINTER|CPL_TYPE_INT:
        /* Case real array column */
        if (!cpl_table_has_column (out_table, name))
            cpl_table_new_column_array (out_table, name, CPL_TYPE_DOUBLE, depth);
        for (int tel = 0; tel < ntel; tel++) {
            cpl_size nvalid = 0;
            cpl_array * mean = gravi_array_init_double (depth, 0.0);
            for (cpl_size row = 0; row < nrow; row++) {
                if (flag && flag[row*ntel+tel] != 0) continue;
                nvalid ++;
                cpl_array_add (mean, cpl_table_get_array (in_table, name, row*ntel+tel));
                CPLCHECK_MSG ("Cannot add arrays...");
            }
            if (nvalid > 0) cpl_array_divide_scalar (mean, nvalid);
            cpl_table_set_array (out_table, name, tel, mean);
            FREE (cpl_array_delete, mean);
        }
        break;

	case CPL_TYPE_POINTER|CPL_TYPE_DOUBLE_COMPLEX:
	case CPL_TYPE_POINTER|CPL_TYPE_FLOAT_COMPLEX:
        /* Case complex array column */
        if (!cpl_table_has_column (out_table, name))
            cpl_table_new_column_array (out_table, name, CPL_TYPE_DOUBLE_COMPLEX, depth);
        for (int tel = 0; tel < ntel; tel++) {
            cpl_size nvalid = 0;
            cpl_array * mean = gravi_array_init_double_complex (depth, 0.0*I+0.0);
            for (cpl_size row = 0; row < nrow; row++) {
                if (flag && flag[row*ntel+tel] != 0) continue;
                nvalid ++;
                cpl_array_add (mean, cpl_table_get_array (in_table, name, row*ntel+tel));
            }
            if (nvalid > 0) cpl_array_divide_scalar (mean, nvalid);
            cpl_table_set_array (out_table, name, tel, mean);
            
            FREE (cpl_array_delete, mean);
        }
        break;
        
        cpl_msg_error (cpl_func, "Type column not yet supported...");
        cpl_error_set_message (cpl_func, CPL_ERROR_ILLEGAL_INPUT,"This type is not supported.");
        return CPL_ERROR_ILLEGAL_INPUT;
    }

    /* Copy units */
    cpl_table_set_column_unit (out_table, name, cpl_table_get_column_unit (in_table, name));

    gravi_msg_function_exit(0);
    return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Force all data in OI_TABLE to have the same TIME and MJD.
 * 
 * @param oi_data    VIS data to process, in-place
 * 
 * The TIME and MJD columns of all OI_VIS, OI_VIS2, OI_FLUX and
 * OI_T3 are averaged and replaced by these averages. Note that:
 * SC and FT are averaged into the same time/mjd, all observations are
 * averaged into the same time/mjd, all observables are averaged into
 * the same time/mjd.
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_vis_force_time (gravi_data * oi_data)
{
    gravi_msg_function_start(1);
    cpl_ensure_code (oi_data, CPL_ERROR_NULL_INPUT);
    
    /* Number of extensions */
    cpl_size nb_ext = gravi_data_get_size (oi_data);
    cpl_ensure_code (nb_ext>0,   CPL_ERROR_ILLEGAL_INPUT);

    /* Init averaging */
    double mean_mjd = 0.0;
    double mean_time = 0.0;
    cpl_size count = 0;

    /* Loop on extensions */
    for (int ext = 0; ext < nb_ext; ext++) {

        /* Check if data table */
        const char * extname = gravi_data_get_extname (oi_data, ext);
        if (!strcmp (extname, "OI_VIS") || !strcmp (extname, "OI_VIS2") ||
            !strcmp (extname, "OI_T3") || !strcmp (extname, "OI_FLUX")) {

            /* Average MJD and TIME */
            cpl_table * oi_table = gravi_data_get_table_x (oi_data, ext);
            mean_mjd  += cpl_table_get_column_mean (oi_table, "MJD");
            mean_time += cpl_table_get_column_mean (oi_table, "TIME");
            count ++;

            CPLCHECK_MSG ("Cannot get TIME or MJD...");
        }
    }

    /* Compute mean */
    cpl_ensure_code (count>0, CPL_ERROR_ILLEGAL_INPUT);
    mean_mjd /= count;
    mean_time /= count;

    /* Verbose */
    cpl_msg_info (cpl_func, "Mean MDJ = %g [mdj]", mean_mjd);
    cpl_msg_info (cpl_func, "Mean TIME = %g [s]", mean_time);

    /* Loop on extensions */
    for (int ext = 0; ext < nb_ext; ext++) {

        /* Check if data table */
        const char * extname = gravi_data_get_extname (oi_data, ext);
        if (!strcmp (extname, "OI_VIS") || !strcmp (extname, "OI_VIS2") ||
            !strcmp (extname, "OI_T3") || !strcmp (extname, "OI_FLUX")) {

            /* Set MJD and TIME */
            cpl_table * oi_table = gravi_data_get_table_x (oi_data, ext);
            cpl_table_fill_column_window (oi_table, "MJD",  0, CPL_SIZE_MAX, mean_mjd);
            cpl_table_fill_column_window (oi_table, "TIME", 0, CPL_SIZE_MAX, mean_time);
            
            CPLCHECK_MSG ("Cannot set average TIME or MJD...");
        }
    }
    
    gravi_msg_function_exit(1);
    return CPL_ERROR_NONE;
}  


/**@}*/
