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
 */
/**@{*/

/*-----------------------------------------------------------------------------
                                   Includes
 -----------------------------------------------------------------------------*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#define _XOPEN_SOURCE
#include <cpl.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <complex.h>

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
											int nseg, int nboot, int tel);
cpl_error_code gravi_t3_average_bootstrap(cpl_table * oi_t3_avg,
										  cpl_table * oi_vis,
										  cpl_table * oi_flux,
										  int nseg, int nboot,
										  int closure,
										  int use_vFactor,
                                          int use_pFactor);
cpl_error_code gravi_vis_average_bootstrap (cpl_table * oi_vis_avg,
											cpl_table * oi_vis2_avg,
											cpl_table * oi_vis,
											int nseg, int nboot,
											int base,
											int use_exp_phase,
											int use_vFactor,
                                            int use_pFactor,
											int use_debiasing);

cpl_error_code gravi_vis_average_amp (cpl_table *oi_table, const char *name,  const char *err, int nbase);
cpl_error_code gravi_vis_average_phi (cpl_table *oi_table, const char *name,  const char *err, int nbase);
cpl_error_code gravi_vis_average_value (cpl_table *oi_table, const char *name,  const char *err, int nbase);
cpl_error_code gravi_vis_resamp_amp (cpl_table * oi_table, const char * name, const char * err,
									 cpl_size nsamp, cpl_size nwave_new);
cpl_error_code gravi_vis_resamp_phi (cpl_table * oi_table, const char * name, const char * err,
									 cpl_size nsamp, cpl_size nwave_new);

cpl_error_code gravi_vis_compute_column_mean (cpl_table * out_table,
                                              cpl_table * in_table,
                                              const char * name, int ntel);

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
											int nseg,
											int nboot,
											int tel)
{
   gravi_msg_function_start(0);
   cpl_ensure_code (oi_flux_avg, CPL_ERROR_ILLEGAL_OUTPUT);
   cpl_ensure_code (oi_flux,     CPL_ERROR_NULL_INPUT);
   cpl_ensure_code (nseg>0,  CPL_ERROR_ILLEGAL_INPUT);
   cpl_ensure_code (nboot>0, CPL_ERROR_ILLEGAL_INPUT);
   cpl_ensure_code (tel>=0,  CPL_ERROR_ILLEGAL_INPUT);
  
  /* Tel for base and base for closure */
  int nv = 0, ntel = 4;
  cpl_size nvalid = 0;
  cpl_size nrow = cpl_table_get_nrow (oi_flux) / ntel;
  cpl_size nwave = cpl_table_get_column_depth (oi_flux, "FLUX");

  /* Pointer to columns, to speed-up */
  cpl_array ** pFLUX     = cpl_table_get_data_array (oi_flux, "FLUX");
  cpl_array ** pFLUXERR  = cpl_table_get_data_array (oi_flux, "FLUXERR");
  double * pINTTIME = cpl_table_get_data_double (oi_flux, "INT_TIME");
  double * pMJD     = cpl_table_get_data_double (oi_flux, "MJD");
  CPLCHECK_MSG ("Cannot get the data");

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
  cpl_size nrow_per_seg = CPL_MAX( nvalid / nseg, 1);
  nseg = nvalid / nrow_per_seg;

  /* Ensure there are at least 5 samples to bootstrap on,
   * if no add montecarlo samples */
  cpl_size nsamp = 5, nmontecarlo = CPL_MAX (nsamp - nseg, 0);

  cpl_msg_info ("Stat", "%6lld valid frames over %6lld (%5.1f%%), make %4d seg. of %5lld (miss %lld), add %lld MonteCarlo",
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

  /* 
   * (3) Save the results on the oi_t3_avg tables 
   */

  /* Save the cloture amplitude on the oi_T3 tables  */
  cpl_table_set_array (oi_flux_avg, "FLUX", tel, flux_res[2]);
  cpl_table_set_array (oi_flux_avg, "FLUXERR", tel, flux_res[3]);
  CPLCHECK_MSG("filling FLUX and FLUXERR");

  /* Flag the data with >100% error */
  gravi_vis_flag_relative_threshold (oi_flux_avg, "FLUXERR", "FLUX", "FLAG", 1.0);
  CPLCHECK_MSG("cannot flag baddata data");

  /* Compute the total integration time */
  cpl_msg_debug(cpl_func,"Total integration time = %.3f s", total_exptime);
  cpl_table_set_double (oi_flux_avg, "INT_TIME", tel, total_exptime);
  cpl_table_set_double (oi_flux_avg, "MJD", tel, mjd_avg);

  /* Set the TARGET_ID and STA_INDEX */
  cpl_table_set_int (oi_flux_avg, "TARGET_ID", tel, cpl_table_get_int (oi_flux, "TARGET_ID", tel, &nv));
  cpl_table_set_int (oi_flux_avg, "STA_INDEX", tel, cpl_table_get_int (oi_flux, "STA_INDEX", tel, &nv));

  FREELOOP (cpl_array_delete, flux_res, 4);
  cpl_free(flag);

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
										  int nseg,
										  int nboot,
										  int closure,
										  int use_vFactor,
                                          int use_pFactor)
{
  gravi_msg_function_start(0);
  cpl_ensure_code (oi_t3_avg, CPL_ERROR_ILLEGAL_OUTPUT);
  cpl_ensure_code (oi_vis,    CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (oi_flux,   CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (nseg>0,     CPL_ERROR_ILLEGAL_INPUT);
  cpl_ensure_code (nboot>0,    CPL_ERROR_ILLEGAL_INPUT);
  cpl_ensure_code (closure>=0, CPL_ERROR_ILLEGAL_INPUT);

  /* Tel for base and base for closure 
   * Assume all the observations are made with the same array geometry.
   * sta_index_t3[0] <-> tel1[clo[j][0]] = tel1[clo[j][2]]
   * sta_index_t3[1] <-> tel1[clo[j][1]] = tel2[clo[j][0]]
   * sta_index_t3[2] <-> tel2[clo[j][2]] = tel2[clo[j][1]]
   */
  int nv = 0, nbase = 6, ntel = 4;
  int base0 = GRAVI_CLO_BASE[closure][0];
  int base1 = GRAVI_CLO_BASE[closure][1];
  int base2 = GRAVI_CLO_BASE[closure][2];
  int ctel0 = GRAVI_CLO_TEL[closure][0];
  int ctel1 = GRAVI_CLO_TEL[closure][1];
  int ctel2 = GRAVI_CLO_TEL[closure][2];

  cpl_size nvalid = 0;
  cpl_size nrow = cpl_table_get_nrow (oi_vis) / nbase;
  cpl_size nwave = cpl_table_get_column_depth (oi_vis, "VISDATA");

  /* Pointer to column, to speed-up */
  cpl_array ** pVISDATA = cpl_table_get_data_array (oi_vis, "VISDATA");
  cpl_array ** pVISERR  = cpl_table_get_data_array (oi_vis, "VISERR");
  cpl_array ** pFLUX    = cpl_table_get_data_array (oi_flux, "FLUX");
  double * pINTTIME = cpl_table_get_data_double (oi_vis, "INT_TIME");
  double * pMJD     = cpl_table_get_data_double (oi_vis, "MJD");
  double * pUCOORD  = cpl_table_get_data_double (oi_vis, "UCOORD");
  double * pVCOORD  = cpl_table_get_data_double (oi_vis, "VCOORD");
  cpl_array ** pVFACTOR = use_vFactor?cpl_table_get_data_array (oi_vis, "V_FACTOR"):NULL;
  double * pPFACTOR = use_pFactor?cpl_table_get_data_double (oi_vis, "P_FACTOR"):NULL;
  CPLCHECK_MSG ("Cannot get the data");
  
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
  cpl_size nrow_per_seg = CPL_MAX (nvalid / nseg, 1);
  nseg = nvalid / nrow_per_seg;

  /* Ensure there are at least 5 samples to bootstrap on,
   * if no add montecarlo samples */
  cpl_size nsamp = 5, nmontecarlo = CPL_MAX (nsamp - nseg, 0);

  cpl_msg_info ("Stat", "%6lld valid frames over %6lld (%5.1f%%), make %4d seg. of %5lld (miss %lld), add %lld MonteCarlo",
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
          double VFACTOR0 = (use_vFactor?cpl_array_get (pVFACTOR[row * nbase + base0], w, NULL):1.0);
          double VFACTOR1 = (use_vFactor?cpl_array_get (pVFACTOR[row * nbase + base1], w, NULL):1.0);
          double VFACTOR2 = (use_vFactor?cpl_array_get (pVFACTOR[row * nbase + base2], w, NULL):1.0);

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
        gravi_array_threshold_min (f012_boot, 0.0, 1e-15);
        
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

  /* Flag the data with >100% error or >60deg error */
  gravi_vis_flag_threshold (oi_t3_avg, "T3PHIERR", "FLAG", 60.0);
  gravi_vis_flag_threshold (oi_t3_avg, "T3AMPERR", "FLAG", 1.0);
  CPLCHECK_MSG("cannot flag baddata data");
  
  /* Compute the total integration time and MJD */
  cpl_msg_debug(cpl_func,"Total integration time = %.3f s", total_exptime);
  cpl_table_set_double (oi_t3_avg, "INT_TIME", closure, total_exptime);
  cpl_table_set_double (oi_t3_avg, "MJD", closure, mjd_avg);
  cpl_table_set_double (oi_t3_avg, "U1COORD", closure, u1Coord);
  cpl_table_set_double (oi_t3_avg, "V1COORD", closure, v1Coord);
  cpl_table_set_double (oi_t3_avg, "U2COORD", closure, u2Coord);
  cpl_table_set_double (oi_t3_avg, "V2COORD", closure, v2Coord);
  
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
 * @param use_exp_phase:  use the PHASE_REF to rephase the DITs
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
											int nseg,
											int nboot,
											int base,
											int use_exp_phase,
											int use_vFactor,
                                            int use_pFactor,
											int use_debiasing)
{
  gravi_msg_function_start(0);
  cpl_ensure_code (oi_vis_avg,  CPL_ERROR_ILLEGAL_OUTPUT);
  cpl_ensure_code (oi_vis2_avg, CPL_ERROR_ILLEGAL_OUTPUT);
  cpl_ensure_code (oi_vis,      CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (nseg>0,  CPL_ERROR_ILLEGAL_INPUT);
  cpl_ensure_code (nboot>0, CPL_ERROR_ILLEGAL_INPUT);
  cpl_ensure_code (base>=0, CPL_ERROR_ILLEGAL_INPUT);
  
  /* Tel for base */
  int nv = 0, nbase = 6;
  
  cpl_size nvalid = 0;
  cpl_size nrow = cpl_table_get_nrow (oi_vis) / nbase;
  cpl_size nwave = cpl_table_get_column_depth (oi_vis, "VISDATA");

  /* Pointer to columns, to speed-up */
  cpl_array ** pVISDATA = cpl_table_get_data_array (oi_vis, "VISDATA");
  cpl_array ** pVISERR  = cpl_table_get_data_array (oi_vis, "VISERR");
  cpl_array ** pFNORM   = cpl_table_get_data_array (oi_vis, "F1F2");
  double * pINTTIME = cpl_table_get_data_double (oi_vis, "INT_TIME");
  double * pMJD     = cpl_table_get_data_double (oi_vis, "MJD");
  double * pUCOORD  = cpl_table_get_data_double (oi_vis, "UCOORD");
  double * pVCOORD  = cpl_table_get_data_double (oi_vis, "VCOORD");
  cpl_array ** pPHASEREF = use_exp_phase?cpl_table_get_data_array (oi_vis, "PHASE_REF"):NULL;
  cpl_array ** pVFACTOR = use_vFactor?cpl_table_get_data_array (oi_vis, "V_FACTOR"):NULL;
  double * pPFACTOR = use_pFactor?cpl_table_get_data_double (oi_vis, "P_FACTOR"):NULL;
  CPLCHECK_MSG ("Cannot get the data");

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
  cpl_size nrow_per_seg = CPL_MAX (nvalid / nseg, 1);
  nseg = nvalid / nrow_per_seg;

  /* Ensure there are at least 5 samples to bootstrap on,
   * if no add montecarlo samples */
  cpl_size nsamp = 5, nmontecarlo = CPL_MAX (nsamp - nseg, 0);

  cpl_msg_info ("Stat", "%6lld valid frames over %6lld (%5.1f%%), make %4d seg. of %5lld (miss %lld), add %lld MonteCarlo",
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
            CPLCHECK_MSG ("Cannot get data");

			/* Add noise if this is a Monte Carlo sample.
			 * APPROX: Noise is only added to the coherent flux */
			if ( seg > nseg-1 ) {
			  mR += 2 * eR * gravi_randn();
			  mI += 2 * eI * gravi_randn();
			}
    
    		/* Integrate <R> and <I> rephased by exp_phase if needed 
    		 * The phase is exp{R,I} * exp_phase */
            tR[w] += cos(PHASEREF) * mR - sin(PHASEREF) * mI;
            tI[w] += cos(PHASEREF) * mI + sin(PHASEREF) * mR;
    
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
        gravi_array_threshold_min (F1F2_boot, 0.0, 1e-15);
    			  
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
  gravi_vis_flag_threshold (oi_vis_avg,  "VISAMPERR", "FLAG", 1.0);
  CPLCHECK_MSG("cannot flag baddata data");
  
  /* Compute the total integration time */
  cpl_table_set_double (oi_vis_avg,  "INT_TIME", base, total_exptime);
  cpl_table_set_double (oi_vis2_avg, "INT_TIME", base, total_exptime);
  cpl_table_set_double (oi_vis_avg,  "MJD", base, mjd_avg);
  cpl_table_set_double (oi_vis2_avg, "MJD", base, mjd_avg);
  cpl_table_set_double (oi_vis_avg,  "UCOORD", base, uCoord);
  cpl_table_set_double (oi_vis_avg,  "VCOORD", base, vCoord);
  cpl_table_set_double (oi_vis2_avg, "UCOORD", base, uCoord);
  cpl_table_set_double (oi_vis2_avg, "VCOORD", base, vCoord);
  
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
  
  gravi_msg_function_exit(0);
  return CPL_ERROR_NONE;
}
 
/*----------------------------------------------------------------------------*/
/**
 * @brief The function average the individual frames of a P2VMREDUCED file
 *        into a final, single observation per base and per tel.
 * 
 * @param p2vmred_data  P2VMREDUCED data containing OI_FLUX and OI_VIS tables
 * 	  	                comming from the @c gravi_compute_p2vmred
 * @param parlist       Parameters of the recipe
 */
/*----------------------------------------------------------------------------*/

gravi_data * gravi_compute_vis (gravi_data * p2vmred_data,
                                const cpl_parameterlist * parlist)
{
	gravi_msg_function_start(1);
	cpl_ensure (p2vmred_data, CPL_ERROR_NULL_INPUT, NULL);
	cpl_ensure (parlist,      CPL_ERROR_NULL_INPUT, NULL);
	
	int nv, nbase = 6, ntel=4, nclo=4;
	char qc_name[100];

	/* 
	 * Prepare the output 
	 */
	
    cpl_msg_info(cpl_func, "Construction of the averaged output data");

    gravi_data * vis_data = gravi_data_new (0);
	cpl_propertylist * p2vmred_header = gravi_data_get_header (p2vmred_data);
    cpl_propertylist * vis_header = gravi_data_get_header (vis_data);
    cpl_propertylist_append (vis_header, p2vmred_header);
		
    /* Copy the oi tables needed in output data.
     * This will duplicate all OI_WAVELENGTH tables */
    gravi_data_copy_ext (vis_data, p2vmred_data, GRAVI_OI_ARRAY_EXT);
    gravi_data_copy_ext (vis_data, p2vmred_data, GRAVI_OI_TARGET_EXT);
    gravi_data_copy_ext (vis_data, p2vmred_data, GRAVI_OI_WAVELENGTH_EXT);
    
    CPLCHECK_NUL ("Cannot get tables for output data");

    /* Output header */
    cpl_propertylist * plist = gravi_data_get_header (vis_data);

    
    /* Start with FT */
    if (gravi_data_has_type (p2vmred_data, "_FT") <= 0 ) {
		cpl_msg_info (cpl_func, "P2VMRED data has no FT extensions");
    }
    else {
		/* Reduction parameters */
		int exp_phase_flag_ft = 1;
		int v_factor_flag_ft = 0;
		int p_factor_flag_ft = 0;
		int debiasing_flag_ft = gravi_param_get_bool (parlist, "gravity.vis.debias-ft");
        int nboot_ft = gravi_param_get_int (parlist, "gravity.vis.nboot");

        cpl_msg_info (cpl_func, "Bias subtraction of V2 for FT is %s",debiasing_flag_ft?"ENABLE":"DISABLE");
		cpl_msg_info (cpl_func, "reference_phase for FT is %s",exp_phase_flag_ft?"ENABLE":"DISABLE");

		CPLCHECK_NUL("Cannot get parameters");

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
            cpl_size nrow_ft  = cpl_table_get_nrow (vis_FT) / nbase;
            int nseg_ft = CPL_MIN (nrow_ft, 100);
            CPLCHECK_NUL ("Cannot get data");

            /* 
             * Loop on bases to compute OIVIS2 and OIVIS for FT
             */
            cpl_msg_info (cpl_func, "Compute OIVIS2 and OIVIS for FT");
            
            cpl_table * oi_vis2_FT = gravi_table_oi_create (nwave_ft, 1, GRAVI_OI_VIS2_EXT);
            cpl_table * oi_vis_FT = gravi_table_oi_create (nwave_ft, 1, GRAVI_OI_VIS_EXT);

            /* Get required data */
            int * reject_flag_ft = cpl_table_get_data_int (vis_FT, "REJECTION_FLAG");
            CPLCHECK_NUL ("Cannot get data");
            
            for (int base = 0; base < nbase; base++) {
                
                gravi_vis_average_bootstrap (oi_vis_FT, oi_vis2_FT, vis_FT,
                                             nseg_ft, nboot_ft, base,
                                             exp_phase_flag_ft,
                                             v_factor_flag_ft,
                                             p_factor_flag_ft,
                                             debiasing_flag_ft);
                CPLCHECK_NUL("Cannot average the FT frames");
                
                /* Add the QC parameters for FT */
                sprintf (qc_name, "ESO QC ACCEPTED_RATIO_FT%s_P%d", GRAVI_BASE_NAME[base], pol+1);
                double ratio = 0.0; for (cpl_size r=0; r<nrow_ft;r++) ratio += (reject_flag_ft[r*nbase+base]>0?0:1);
                cpl_propertylist_update_double (plist, qc_name, round(ratio / nrow_ft * 100.0 * 1e2) / 1e2);
                cpl_propertylist_set_comment (plist, qc_name, "[%] of accepted frames");
                
                sprintf (qc_name, "ESO QC VISPHIERR_FT%s_P%d AVG", GRAVI_BASE_NAME[base], pol+1);
                cpl_propertylist_update_double (plist, qc_name, cpl_array_get_mean(cpl_table_get_array (oi_vis_FT, "VISPHIERR", base)));
                cpl_propertylist_set_comment (plist, qc_name, "[deg] mean over lbd");
                
                sprintf (qc_name, "ESO QC VIS2_FT%s_P%d AVG", GRAVI_BASE_NAME[base], pol+1);
                cpl_propertylist_update_double (plist, qc_name, cpl_array_get_mean(cpl_table_get_array (oi_vis2_FT, "VIS2DATA", base)) );
                cpl_propertylist_set_comment (plist, qc_name, "mean over lbd");
                
                sprintf (qc_name, "ESO QC VIS2ERR_FT%s_P%d AVG", GRAVI_BASE_NAME[base], pol+1);
                cpl_propertylist_update_double (plist, qc_name, cpl_array_get_mean(cpl_table_get_array (oi_vis2_FT, "VIS2ERR", base)));
                cpl_propertylist_set_comment (plist, qc_name, "mean over lbd");
                
                sprintf (qc_name, "ESO QC VISAMP_FT%s_P%d AVG", GRAVI_BASE_NAME[base], pol+1);
                cpl_propertylist_update_double (plist, qc_name, cpl_array_get_mean(cpl_table_get_array (oi_vis_FT, "VISAMP", base)));
                cpl_propertylist_set_comment (plist, qc_name, "mean over lbd");
                
                sprintf (qc_name, "ESO QC VISAMPERR_FT%s_P%d AVG", GRAVI_BASE_NAME[base], pol+1);
                cpl_propertylist_update_double (plist, qc_name, cpl_array_get_mean(cpl_table_get_array (oi_vis_FT, "VISAMPERR", base)));
                cpl_propertylist_set_comment (plist, qc_name, "mean over lbd");
                
                CPLCHECK_NUL("Cannot compute QC parameter for OI_VIS for FT");
            } /* End loop on base */

            /* 
             * Loop on triplet to compute OIT3 for FT
             */
            cpl_msg_info (cpl_func, "Compute OIT3 for FT");

            cpl_table * oi_T3_FT = gravi_table_oi_create (nwave_ft, 1, GRAVI_OI_T3_EXT);
            
            for (int clo = 0; clo < nclo; clo++){
                
                gravi_t3_average_bootstrap (oi_T3_FT, vis_FT, flux_FT,
                                            nseg_ft, nboot_ft, clo,
                                            v_factor_flag_ft,
                                            p_factor_flag_ft);
                
                CPLCHECK_NUL("Cannot average t3 of FT");
                
                /* Add QC */
                sprintf (qc_name, "ESO QC T3PHI_FT%s_P%d AVG", GRAVI_CLO_NAME[clo], pol+1);
                cpl_propertylist_update_double (plist, qc_name, cpl_array_get_mean (cpl_table_get_array (oi_T3_FT, "T3PHI", clo)));
                cpl_propertylist_set_comment (plist, qc_name, "[deg] mean over lbd");
                
                sprintf (qc_name, "ESO QC T3PHIERR_FT%s_P%d AVG", GRAVI_CLO_NAME[clo], pol+1);
                cpl_propertylist_update_double (plist, qc_name, cpl_array_get_mean(cpl_table_get_array (oi_T3_FT, "T3PHIERR", clo)));
                cpl_propertylist_set_comment (plist, qc_name, "[deg] mean over lbd");
                
                CPLCHECK_NUL("Cannot compute QC parameter for OI_T3 for FT");
            } /* End loop on triplets */
            
            /* 
             * Loop on beams to compute OI_FLUX for FT
             */
            cpl_msg_info (cpl_func, "Compute OI_FLUX for FT");
		
            cpl_table * oi_flux_FT = gravi_table_oi_create (nwave_ft, 1, GRAVI_OI_FLUX_EXT);
            
            for (int tel = 0; tel < ntel; tel++){
                
                gravi_flux_average_bootstrap (oi_flux_FT, flux_FT,
                                              nseg_ft, nboot_ft, tel);
                
                CPLCHECK_NUL("Cannot average flux of FT");
                
                /* Add QC */
                sprintf (qc_name, "ESO QC FLUX_FT%d_P%d AVG", tel+1, pol+1);
                cpl_propertylist_update_double (plist, qc_name, cpl_array_get_mean (cpl_table_get_array (oi_flux_FT, "FLUX", tel)));
                cpl_propertylist_set_comment (plist, qc_name, "[e/total_int_time] mean over lbd");
                
                sprintf (qc_name, "ESO QC FLUXERR_FT%d_P%d AVG", tel+1, pol+1);
                cpl_propertylist_update_double (plist, qc_name, cpl_array_get_mean (cpl_table_get_array (oi_flux_FT, "FLUXERR", tel)));
                cpl_propertylist_set_comment (plist, qc_name, "[e/total_int_time] mean over lbd");
                
                sprintf (qc_name, "ESO QC FLUXRATE_FT%d_P%d SUM", tel+1, pol+1);
                double flux_rate = cpl_array_get_mean (cpl_table_get_array (oi_flux_FT, "FLUX", tel)) *
                    cpl_array_get_size (cpl_table_get_array (oi_flux_FT, "FLUX", tel)) / cpl_table_get_double (oi_flux_FT, "INT_TIME", tel, &nv);
                cpl_propertylist_update_double (plist, qc_name, flux_rate);
                cpl_propertylist_set_comment (plist, qc_name, "[e/s] sum over lbd");
                
                CPLCHECK_NUL("Cannot compute QC parameter for OI_FLUX for FT");
            } /* End loop on beams */
            
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
            
        } /* end loop on pol */
    } /* End FT */


    
    /* Then with SC */
    if (gravi_data_has_type (p2vmred_data, "_SC") <= 0 ) {
		cpl_msg_info (cpl_func, "P2VMRED data has no SC extensions");
    }
    else {
		/* Reduction parameters */
		int exp_phase_flag_sc = 1;
		int v_factor_flag_sc = strstr (gravi_param_get_string (parlist, "gravity.vis.vis-correction-sc"),"VFACTOR") ? 1 : 0;
		int p_factor_flag_sc = strstr (gravi_param_get_string (parlist, "gravity.vis.vis-correction-sc"),"PFACTOR") ? 1 : 0;
		int debiasing_flag_sc = gravi_param_get_bool (parlist, "gravity.vis.debias-sc");
		int nboot_sc = gravi_param_get_int (parlist, "gravity.vis.nboot");
		
		cpl_msg_info (cpl_func, "Bias subtraction of V2 for SC is %s",debiasing_flag_sc?"ENABLE":"DISABLE");
		cpl_msg_info (cpl_func, "reference_phase for SC is %s",exp_phase_flag_sc?"ENABLE":"DISABLE");
		cpl_msg_info (cpl_func, "vFactor correction for SC is %s",v_factor_flag_sc?"ENABLE":"DISABLE");
		cpl_msg_info (cpl_func, "pFactor correction for SC is %s",p_factor_flag_sc?"ENABLE":"DISABLE");

		CPLCHECK_NUL("Cannot get parameters");

        /* 
         * Loop on polarisations
         */
        int npol_sc = gravi_pfits_get_pola_num (p2vmred_header, GRAVI_SC);
        for (int pol = 0; pol < npol_sc; pol++) {
            cpl_msg_info (cpl_func, "Start SC polarisation %d over %d",pol+1, npol_sc);
            
            /* Get the input table of SC */
            cpl_table * vis_SC = gravi_data_get_oi_vis (p2vmred_data, GRAVI_SC, pol, npol_sc);
            cpl_table * flux_SC = gravi_data_get_oi_flux (p2vmred_data, GRAVI_SC, pol, npol_sc);
            cpl_table * oi_wavelengthsc = gravi_data_get_oi_wave (p2vmred_data, GRAVI_SC, pol, npol_sc);
            int nwave_sc = cpl_table_get_column_depth (vis_SC, "VISDATA");
            cpl_size nrow_sc  = cpl_table_get_nrow (vis_SC) / nbase;
            int nseg_sc = CPL_MIN (nrow_sc, 100);
            
            CPLCHECK_NUL ("Cannot get data");

            /* Get the wavelength array of FT and SC data in [m] 
             * Also compute the wavenumber for SC in [m^-1] */
            cpl_array * wavenumber_sc;
            wavenumber_sc = cpl_array_new (nwave_sc, CPL_TYPE_DOUBLE);
            for (cpl_size wave = 0; wave < nwave_sc; wave ++){
                cpl_array_set (wavenumber_sc, wave, 1./cpl_table_get (oi_wavelengthsc, "EFF_WAVE", wave, &nv));
            }

            CPLCHECK_NUL ("Cannot build the wave and wavenumber");

            /* 
             * Loop on bases to compute OIVIS2 and OIVIS for SC
             */
            cpl_msg_info (cpl_func, "Compute OIVIS2 and OIVIS for SC");

            cpl_table * oi_vis2_SC = gravi_table_oi_create (nwave_sc, 1, GRAVI_OI_VIS2_EXT);
            cpl_table * oi_vis_SC = gravi_table_oi_create (nwave_sc, 1, GRAVI_OI_VIS_EXT);
            
            /* Additional columns in final, averaged product */
            gravi_table_new_column (oi_vis_SC, "GDELAY", "m", CPL_TYPE_DOUBLE);
            gravi_table_new_column (oi_vis_SC, "PHASE", "rad", CPL_TYPE_DOUBLE);

            CPLCHECK_NUL("Cannot create columns in averaged OIFITS...");

            // Shall consider the case where OPD_DISP is NULL or absent
            // gravi_vis_compute_column_mean (oi_vis_SC, vis_SC, "OPD_DISP", 6);
            
            gravi_vis_compute_column_mean (oi_vis_SC, vis_SC, "OPD_MET_FC", 6);
            gravi_vis_compute_column_mean (oi_vis_SC, vis_SC, "PHASE_REF_COEFF", 6);
            gravi_vis_compute_column_mean (oi_vis_SC, vis_SC, "E_U", 6);
            gravi_vis_compute_column_mean (oi_vis_SC, vis_SC, "E_V", 6);
            gravi_vis_compute_column_mean (oi_vis_SC, vis_SC, "E_W", 6);
            gravi_vis_compute_column_mean (oi_vis_SC, vis_SC, "E_AZ", 6);
            gravi_vis_compute_column_mean (oi_vis_SC, vis_SC, "E_ZD", 6);
            
            CPLCHECK_NUL("Cannot create columns in averaged OIFITS...");
            
            for (int base = 0; base < nbase; base++) {
                
                gravi_vis_average_bootstrap (oi_vis_SC, oi_vis2_SC, vis_SC,
                                             nseg_sc, nboot_sc, base,
                                             exp_phase_flag_sc,
                                             v_factor_flag_sc,
                                             p_factor_flag_sc,
                                             debiasing_flag_sc);
                
                CPLCHECK_NUL("Cannot average the SC frames");
                
                /* 
                 * Re compute the astrometric phase from VISDATA to deal with absolute phase
                 * VISDATA as well as (R,I) remains unchanged.
                 */

                /* We duplicate VISDATA, to keep the value untouched in the table */
                cpl_array * visData_sc, * visErr_sc;
                visData_sc = cpl_array_cast (cpl_table_get_array (oi_vis_SC, "VISDATA", base), CPL_TYPE_DOUBLE_COMPLEX);
                
                /* Normalize the phasor to avoid bad pixels */
                visErr_sc = cpl_array_duplicate (visData_sc);
                cpl_array_abs (visErr_sc);
                cpl_array_divide (visData_sc, visErr_sc);
                
                /* Compute and remove the mean group delay in [m] */
                double mean_delay = 0.0;
                gravi_array_get_group_delay_loop (&visData_sc, wavenumber_sc, &mean_delay, 1, 2e-3, CPL_FALSE);
                gravi_array_multiply_phasor (visData_sc, - 2*I*CPL_MATH_PI * mean_delay, wavenumber_sc);
                cpl_msg_debug (cpl_func, "group-delay in SC is : %f [microns]", mean_delay * 1e6);
                
                /* Save this delay [m] */
                cpl_table_set (oi_vis_SC, "GDELAY", base, mean_delay);
                
                /* Compute and remove the mean phase in [rad] */
                double mean_phase = carg (cpl_array_get_mean_complex (visData_sc));
                cpl_array_multiply_scalar_complex (visData_sc, cexp(- I * mean_phase));
                cpl_msg_debug (cpl_func, "phase-delay in SC is : %f [deg]", mean_phase * 180 / CPL_MATH_PI);
                
                /* Save this phase [rad] */
                cpl_table_set (oi_vis_SC, "PHASE", base, mean_phase);
                
                /* Set back the phase in [deg] */
                cpl_array_arg (visData_sc);
                gravi_table_set_array_phase (oi_vis_SC, "VISPHI", base, visData_sc);
                cpl_array_delete (visData_sc);
                cpl_array_delete (visErr_sc);
                
                CPLCHECK_NUL("when computing the astrometric phase");

                /* Add QC parameters for SC */
                int * reject_flag_sc = cpl_table_get_data_int (vis_SC, "REJECTION_FLAG");
                
                CPLCHECK_NUL ("Cannot get data");
                
                sprintf (qc_name, "ESO QC VFACTOR%s_P%d AVG", GRAVI_BASE_NAME[base], pol+1);
                double vmean = gravi_table_get_column_mean (vis_SC, "V_FACTOR_WL", base, nbase);
                cpl_propertylist_update_double (plist, qc_name, vmean);
                cpl_propertylist_set_comment (plist, qc_name, "mean v-factor");
                
                sprintf (qc_name, "ESO QC PFACTOR%s_P%d AVG", GRAVI_BASE_NAME[base], pol+1);
                double pmean = gravi_table_get_column_mean (vis_SC, "P_FACTOR", base, nbase);
                cpl_propertylist_update_double (plist, qc_name, pmean);
                cpl_propertylist_set_comment (plist, qc_name, "mean p-factor");
                
                sprintf (qc_name, "ESO QC ACCEPTED_RATIO_SC%s_P%d", GRAVI_BASE_NAME[base], pol+1);
                double ratio = 0.0; for (cpl_size r=0; r<nrow_sc;r++) ratio += (reject_flag_sc[r*nbase+base]>0?0:1);
                cpl_propertylist_update_double (plist, qc_name, round(ratio / nrow_sc * 100.0 * 1e2) / 1e2);
                cpl_propertylist_set_comment (plist, qc_name, "[%] of accepted frames");
                
                sprintf (qc_name, "ESO QC GD_SC%s_P%d AVG", GRAVI_BASE_NAME[base], pol+1);
                cpl_propertylist_update_double (plist, qc_name, mean_delay);
                cpl_propertylist_set_comment (plist, qc_name, "[m] mean Group-Delay");
                
                sprintf (qc_name, "ESO QC VIS2_SC%s_P%d AVG", GRAVI_BASE_NAME[base], pol+1);
                cpl_propertylist_update_double (plist, qc_name, cpl_array_get_mean (cpl_table_get_array (oi_vis2_SC, "VIS2DATA", base)));
                cpl_propertylist_set_comment (plist, qc_name, "mean over lbd");
                
                sprintf (qc_name, "ESO QC VIS2ERR_SC%s_P%d AVG", GRAVI_BASE_NAME[base], pol+1);
                cpl_propertylist_update_double (plist, qc_name, cpl_array_get_mean (cpl_table_get_array (oi_vis2_SC, "VIS2ERR", base)));
                cpl_propertylist_set_comment (plist, qc_name, "mean over lbd");
                
                sprintf (qc_name, "ESO QC VISPHI_SC%s_P%d AVG", GRAVI_BASE_NAME[base], pol+1);
                cpl_propertylist_update_double (plist, qc_name, cpl_array_get_mean (cpl_table_get_array (oi_vis_SC, "VISPHI", base)));
                cpl_propertylist_set_comment (plist, qc_name, "[deg] mean over lbd");
                
                sprintf (qc_name, "ESO QC VISPHIERR_SC%s_P%d AVG", GRAVI_BASE_NAME[base], pol+1);
                cpl_propertylist_update_double (plist, qc_name, cpl_array_get_mean (cpl_table_get_array (oi_vis_SC, "VISPHIERR", base)));
                cpl_propertylist_set_comment (plist, qc_name, "[deg] mean over lbd");
                
                sprintf (qc_name, "ESO QC VISAMP_SC%s_P%d AVG", GRAVI_BASE_NAME[base], pol+1);
                cpl_propertylist_update_double (plist, qc_name, cpl_array_get_mean (cpl_table_get_array (oi_vis_SC, "VISAMP", base)));
                cpl_propertylist_set_comment (plist, qc_name, "mean over lbd");
                
                sprintf (qc_name, "ESO QC VISAMPERR_SC%s_P%d AVG", GRAVI_BASE_NAME[base], pol+1);
                cpl_propertylist_update_double (plist, qc_name, cpl_array_get_mean (cpl_table_get_array (oi_vis_SC, "VISAMPERR", base)));
                cpl_propertylist_set_comment (plist, qc_name, "mean over lbd");
                
                double coeff2 = gravi_table_get_value (oi_vis_SC, "PHASE_REF_COEFF", base, 2);
                sprintf (qc_name, "ESO QC PHASE_REF_COEFF2 SC%s_P%d", GRAVI_BASE_NAME[base], pol+1);
                cpl_propertylist_update_double (plist, qc_name, coeff2);
                cpl_propertylist_set_comment (plist, qc_name, "[rad] 2sd order of FT phase");
                
                CPLCHECK_NUL("Cannot set QC parameter for OI_VIS for SC");
            } /* End loop on base */
            
            /* 
             * Loop on triplet to compute OIT3 for SC
             */
            cpl_msg_info (cpl_func, "Compute OIT3 for SC");
            
            cpl_table * oi_T3_SC = gravi_table_oi_create (nwave_sc, 1, GRAVI_OI_T3_EXT);
            
            for (int clo = 0; clo < nclo; clo++){
                
                gravi_t3_average_bootstrap (oi_T3_SC, vis_SC, flux_SC,
                                            nseg_sc, nboot_sc, clo,
                                            v_factor_flag_sc,
                                            p_factor_flag_sc);
                
                CPLCHECK_NUL("Cannot average t3 of SC");
                
                /* Add QC */
                sprintf (qc_name, "ESO QC T3PHI_SC%s_P%d AVG", GRAVI_CLO_NAME[clo], pol+1);
                cpl_propertylist_update_double (plist, qc_name, cpl_array_get_mean(cpl_table_get_array (oi_T3_SC, "T3PHI", clo)));
                cpl_propertylist_set_comment (plist, qc_name, "[deg] mean over lbd");
                
                sprintf (qc_name, "ESO QC T3PHIERR_SC%s_P%d AVG", GRAVI_CLO_NAME[clo], pol+1);
                cpl_propertylist_update_double (plist, qc_name, cpl_array_get_mean(cpl_table_get_array (oi_T3_SC, "T3PHIERR", clo)));
                cpl_propertylist_set_comment (plist, qc_name, "[deg] mean over lbd");
          
                CPLCHECK_NUL("Cannot set QC parameter for OI_T3 for SC");
            }/* End loop on triplets */
            
            /* 
             * Loop on beams to compute OI_FLUX for SC
             */
            cpl_msg_info (cpl_func, "Compute OI_FLUX for SC");
            
            cpl_table * oi_flux_SC = gravi_table_oi_create (nwave_sc, 1, GRAVI_OI_FLUX_EXT);
            
            gravi_table_new_column (oi_flux_SC, "LKDT_MET_FC", "mjd", CPL_TYPE_DOUBLE);

            gravi_vis_compute_column_mean (oi_flux_SC, flux_SC, "OPD_MET_FC", 4);
            gravi_vis_compute_column_mean (oi_flux_SC, flux_SC, "FT_POS", 4);
            gravi_vis_compute_column_mean (oi_flux_SC, flux_SC, "SC_POS", 4);
            gravi_vis_compute_column_mean (oi_flux_SC, flux_SC, "OPL_AIR", 4);
            
            CPLCHECK_NUL ("Cannot create columns");
            
            for (int tel = 0; tel < ntel; tel++){
                
                gravi_flux_average_bootstrap (oi_flux_SC, flux_SC,
                                              nseg_sc, nboot_sc, tel);
                
                CPLCHECK_NUL("Cannot average flux of SC");

                /* Save the FC metrology lock date */
                double lockdate = gravi_pfits_get_metfc_lockmjd (p2vmred_header, tel);
                cpl_table_set (oi_flux_SC, "LKDT_MET_FC", tel, lockdate);
                
                /* Add QC */
                sprintf (qc_name, "ESO QC FLUX_SC%d_P%d AVG", tel+1, pol+1);
                cpl_propertylist_update_double (plist, qc_name, cpl_array_get_mean (cpl_table_get_array (oi_flux_SC, "FLUX", tel)));
                cpl_propertylist_set_comment (plist, qc_name, "[e/total_int_time] mean over lbd");
                
                sprintf (qc_name, "ESO QC FLUXERR_SC%d_P%d AVG", tel+1, pol+1);
                cpl_propertylist_update_double (plist, qc_name, cpl_array_get_mean (cpl_table_get_array (oi_flux_SC, "FLUXERR", tel)));
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
                
                CPLCHECK_NUL("Cannot set QC parameter for OI_FLUX for SC");
            } /* End loop on beams */

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
            cpl_propertylist_update_string (oivis_plist, "PHITYP","differential");
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

	/* Set the variance of the mean */
	if (err==NULL) weight = 1.0;
	cpl_table_set (oi_table, name, base, value / weight);
  } /* End loop on base */

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
  cpl_size row, nrow = cpl_table_get_nrow (oi_table);
  cpl_array ** pdata = cpl_table_get_data_array (oi_table, data);
  cpl_array ** pflag = cpl_table_get_data_array (oi_table, flag);
  
  CPLCHECK_MSG ("Cannot get data");
  
  cpl_size indx, size = cpl_array_get_size (pdata[0]);
  
  /* Loop on row and index. Add to FLAG if data is above threshold */
  for ( row = 0 ; row < nrow ; row ++ ) {
        if (pdata[row]==NULL) continue;
	for ( indx = 0 ; indx < size ; indx ++ ) {
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
 * Note that this routine is not optimized for performance and thus shall be
 * used on the FT tables.
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
    cpl_ensure_code (cpl_table_has_column (in_table, name),
                     CPL_ERROR_ILLEGAL_INPUT);
    cpl_ensure_code (ntel == cpl_table_get_nrow (out_table),
                     CPL_ERROR_ILLEGAL_OUTPUT);

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
            cpl_table_set (out_table, name, tel, mean / nvalid);
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
            cpl_array_divide_scalar (mean, nvalid);
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
            cpl_array_divide_scalar (mean, nvalid);
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


/**@}*/
