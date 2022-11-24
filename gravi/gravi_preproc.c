/* $Id: gravi_preproc.c,v 1.12 2011/04/31 06:10:40 nazouaoui Exp $
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
 * @defgroup gravi_preproc  Preprocessing functions
 *
 * This module implements the preprocessing of the data, which is mainly to convert raw
 * data into spectra. The main function @c gravi_extract_spectrum() is used at several
 * place in the pipeline. And the steps of the preprocessing are description in the sections :
 * - Algorithm/Spectrum Extraction
 * - Algorithm/Re-interpolation to a common wavelength
 *
 */
/**@{*/

/*
 * History 
 *    11/01/2019  Fix Warnings , unused parameter : profile_header
 *                                   brackets of if statements 
 */

/*-----------------------------------------------------------------------------
                                   Includes
 -----------------------------------------------------------------------------*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cpl.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <complex.h>

#include "gravi_data.h"
#include "gravi_calib.h"
#include "gravi_dfs.h"
#include "gravi_pfits.h"
#include "gravi_cpl.h"

#include "gravi_utils.h"
#include "gravi_preproc.h"

/*-----------------------------------------------------------------------------
                              Private prototypes
 -----------------------------------------------------------------------------*/
cpl_table * gravi_table_ft_format (cpl_table * table_ft, cpl_table * sky_table_std,
								   cpl_table * sky_table_avg, cpl_table * badft_table,
                                   int n_region, double gain, const cpl_parameterlist * parlist);

cpl_table * gravi_imglist_sc_collapse (cpl_table * profile_table,
                                       cpl_imagelist * raw_imglist,
                                       cpl_imagelist * rawVar_imglist,
                                       cpl_size startx);

cpl_table * gravi_imglist_sc_collapse_robust (cpl_table * profile_table,
                                              cpl_imagelist * raw_imglist,
                                              cpl_imagelist * rawVar_imglist,
                                              cpl_size startx);

cpl_error_code gravi_interpolate_spectrum_table (cpl_table * spectrum_table,
                                                 cpl_table * wave_table,
                                                 cpl_table ** oiwave_tables,
                                                 cpl_table * detector_table,
                                                 cpl_table * specflat_table,
                                                 int type_data);

/*-----------------------------------------------------------------------------
                              Function code
 -----------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/**
 * @brief Check if the pixel in the BADPIX map is a good pixel
 * 
 * @param bad_img   badpix as img
 * @param x         xpos
 * @param y         ypos
 * 
 * @return int 0/1
 */
/*---------------------------------------------------------------------------*/

int gravi_pixel_is_good (cpl_image * bad_img, int x, int y)
{
	int bad;
	int nv;

	bad = cpl_image_get (bad_img, x, y, &nv);
	if ((bad & BADPIX_DARK) || (bad & BADPIX_RMS) || (bad & BADPIX_FLAT))
		return 0;
	else
		return 1;
}


/*---------------------------------------------------------------------------*/
/**
 * @brief Remove the badpixel of the SC
 * 
 * @param imglist_sc      input data as imglist, remove inplace
 * @param bad_img         badpixel image
 *
 * The bad pixels of imglist_sc, identified by the map bad_img, are
 * re-interpolated from neighboring values with a special care of the 
 * spatial / spectral directions since the spectra are almost 1D.
 */
/*---------------------------------------------------------------------------*/

cpl_error_code gravi_remove_badpixel_sc (cpl_imagelist * imglist_sc, cpl_image * bad_img)
{
	gravi_msg_function_start(1);
	cpl_ensure_code (imglist_sc, CPL_ERROR_NULL_INPUT);
	cpl_ensure_code (bad_img,    CPL_ERROR_NULL_INPUT);
	
	cpl_image * img;
	int middel = 3, size = 7, nv, comp, badpix_comp;
	cpl_vector * x = NULL;

	cpl_size nrow = cpl_imagelist_get_size (imglist_sc);
	cpl_size nx = cpl_image_get_size_x (bad_img);
	cpl_size ny = cpl_image_get_size_y (bad_img);
	
	/* Check consistency of BADPIX and DATA */
	if ((nx != cpl_image_get_size_x (cpl_imagelist_get (imglist_sc, 0))) ||
		(ny != cpl_image_get_size_y (cpl_imagelist_get (imglist_sc, 0))) ){
	  return cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
								   "The image lists have not the same size");
	}

	/* Compute the amount of badpixels */
	cpl_size nbad = 0;
	for (cpl_size k = 0; k < ny; k++)
	  for (cpl_size i = 0; i < nx; i++)
	    if (cpl_image_get (bad_img, i+1, k+1, &nv) != 0) nbad++;

	/* Check the fraction of badpixels */
	if (nbad > 0.25 * nx*ny) {
	  return cpl_error_set_message (cpl_func, CPL_ERROR_ILLEGAL_INPUT,
					"Too many bad pixels (more than 25 percent)");
	}

	/* Loop on the pixels of the image */
	for (cpl_size k = 0; k < ny; k++){
	  for (cpl_size i = 0; i < nx; i++){
		
		/* This is a bad pixel */
		if (cpl_image_get (bad_img, i+1, k+1, &nv) != 0) {
		  
		  /* Loop on the images of the RAW imagelist */
		  for (cpl_size row = 0; row < nrow; row ++){
			
			img = cpl_imagelist_get (imglist_sc, row);
			x = cpl_vector_new (size-1);
			
			if ((i - middel) < 0 ){
			  comp = 0;
			  badpix_comp = 1;
			  for (cpl_size j = 0; j < size; j++){
				if (j != i) {
				  while (cpl_image_get (bad_img, i + badpix_comp, k + 1, &nv) != 0){
					badpix_comp ++;
				  }
				  cpl_vector_set(x, comp, cpl_image_get (img, i + badpix_comp, k + 1, &nv));
				  comp ++ ;
				}
			  }
			  CPLCHECK_MSG("Fail 1");
			}
			else if ((i + middel) >= nx){
			  comp = 0;
			  badpix_comp = 1;
			  for (cpl_size j = 0; (j < size); j++){
				if (j != (nx - i)){
				  while (cpl_image_get (bad_img, nx - (j + 1) + badpix_comp, k + 1, &nv) != 0){
					
					badpix_comp --;
				  }
				  cpl_vector_set(x, comp, cpl_image_get (img, nx - (j + 1) + badpix_comp, k + 1 , &nv));
				  comp ++;
				}
			  }
			  CPLCHECK_MSG("Fail 2");
			}
			else{
			  comp = 0;
			  badpix_comp = 1;
			  int test1 = 0, test2 = 0;
			  for (cpl_size j = 0; j < size; j++){
				if (j != middel){
				  if (i + j - middel + badpix_comp <= 1){
					while (cpl_image_get (bad_img, i + j - middel + badpix_comp, k+1, &nv) != 0){
					  badpix_comp ++;
					  test1 = 1;
					}
				  }
				  else if (i + j - middel + badpix_comp >= nx){
					while (cpl_image_get (bad_img, i + j - middel + badpix_comp, k+1, &nv) != 0){
					  badpix_comp --;
					  test2 = 1;
					}
				  }
				  else {
					if ((test1 == 0) && (test2 == 0)){
					  if (j > middel)
						badpix_comp ++;
					  else
						badpix_comp --;
					}
					else if (test1)
					  badpix_comp ++;
					else
					  badpix_comp --;
				  }
				  cpl_vector_set(x, comp, cpl_image_get (img, i + j - middel + badpix_comp, k+1, &nv));
				  comp ++;
				}
			  }
			  CPLCHECK_MSG("Fail 3");
			}
			
			double mean = cpl_vector_get_median (x);
			cpl_vector_delete (x);
			cpl_image_set (img, i + 1, k + 1, mean);
			
			CPLCHECK_MSG("Fail 4");
			
		  } /* End loop on row */
		} /* End case this is a bad pixel */
	  } /* End loop k */
	} /* End loop k */

	gravi_msg_function_exit(1);
	return CPL_ERROR_NONE;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Extract FT spectrum from PIX column
 * 
 * @param pix_table        Input raw data     [adu]
 * @param skystd_table     Input SKY std data [adu]
 * @param skyavg_table     Input SKY avg data [adu]
 * @param n_region         Number of regions to descramble
 * @param gain             Conversion gain in [ADU/e]
 * 
 * @return A cpl_table with TIME, DATA# and DATAERR#
 *
 * The TIME column is duplicated in output.
 * The pixels from the PIX column are descrambled into
 * DATA# = (RAW-SKY)/gain 
 * DATAERR# = sqrt (max(DATA#,0) + (SKYSTD/gain)^2)
 */
/*---------------------------------------------------------------------------*/

cpl_table * gravi_table_ft_format (cpl_table * pix_table,
								   cpl_table * skystd_table,
                                   cpl_table * skyavg_table,
                                   cpl_table * badft_table,
								   int n_region, double gain, const cpl_parameterlist * parlist)
{
  /* Verbose */
  gravi_msg_function_start(0);
  cpl_ensure (pix_table, CPL_ERROR_NULL_INPUT, NULL);
  cpl_ensure (skystd_table, CPL_ERROR_NULL_INPUT, NULL);
  cpl_ensure (skyavg_table, CPL_ERROR_NULL_INPUT, NULL);
  
  /* Get the number of frames */
  cpl_size nrow = cpl_table_get_nrow (pix_table);
  cpl_table * output_table = cpl_table_new (nrow);

  /* Get PIX of data */
  cpl_array ** arr_data = cpl_table_get_data_array (pix_table, "PIX");
  CPLCHECK_NUL ("Cannot get data");
  
  /* Get PIX of badft */
  cpl_array ** arr_badft = cpl_table_get_data_array (badft_table, "PIX");
  CPLCHECK_NUL ("Cannot get data");

  /* Get the dimensions for the descramling */
  cpl_size npix = cpl_table_get_column_depth (pix_table, "PIX");
  cpl_size sizex = cpl_table_get_column_dimension(pix_table, "PIX", 0);
  cpl_size sizey = cpl_table_get_column_dimension(pix_table, "PIX", 1);

  int npol = n_region > 24 ? 2 : 1;
  cpl_size ny = sizey / npol;

  /* Keep the ny_out to 5 for now to keep the same size.
   * could be ny if the pipeline accept 6 pixels spectra */
  cpl_size ny_out;
  if (gravi_param_get_bool(parlist, "gravity.preproc.extra-pixel-ft")) {
      ny_out =ny;
      cpl_msg_info (cpl_func, "Option extra-pixel-ft applied, using %lld pixels", ny_out);
      }
  else{
      ny_out = 5;
      cpl_msg_info (cpl_func, "Option extra-pixel-ft disabled, using %lld pixels", ny_out);
  }

  /* Number of RAW pixels */
  cpl_size n_output = 24;
  cpl_size nx = sizex / n_output;
  cpl_msg_info (cpl_func, "Descramble %lld x %lld as %lld outputs x %lld pixels x %lld channels x %i pol",
                sizex, sizey, n_output, nx, ny, npol);

  /* Copy the column TIME */
  if (cpl_table_has_column (pix_table, "TIME")) {
      cpl_table_duplicate_column (output_table, "TIME", pix_table, "TIME");
      CPLCHECK_NUL ("Cannot get TIME data");
  }

  /* Get pointer to the sky mean [e] out of the loop */
  double * pSky = cpl_calloc (sizex * sizey, sizeof(double));
  if (skyavg_table) {
      for (cpl_size pix = 0; pix < npix; pix++) {
          pSky[pix] = gravi_table_get_value (skyavg_table, "PIX", 0, pix) / gain;
      }
  }
  CPLCHECK_NUL ("Cannot get the sky data");

  /* Get pointer to the sky variance [e^2] out of the loop */
  double * pSkyVar = cpl_calloc (sizex * sizey, sizeof(double));
  if (skystd_table) {
      for (cpl_size pix = 0; pix < npix; pix++) {
          pSkyVar[pix] = gravi_table_get_value (skystd_table, "PIX", 0, pix);
          pSkyVar[pix] = pow (pSkyVar[pix] / gain, 2.0);
      }
  }
  CPLCHECK_NUL ("Cannot get sky variance");
  

  /* Loop on regions */
  for (cpl_size region = 0; region < n_output; region++) {
	
	/* Loop on polarisation */
	for (int pol = 0; pol < npol; pol ++) {
	  
	  /* Verbose every 6 regions */
	  if ( !(region+pol) || !((region*npol+pol+1)%6) )
		cpl_msg_info_overwritable (cpl_func,
								   "Extract region of FT %lld over %lld (fast-no-cpl)",
								   (region*npol+pol+1), n_output*npol);

	  /* Create DATA column */
	  const char * data = GRAVI_DATA[region*npol + pol];
	  cpl_table_new_column_array (output_table, data, CPL_TYPE_DOUBLE, ny_out);
	  cpl_array ** tData = cpl_table_get_data_array (output_table, data);
	  CPLCHECK_NUL ("Cannot create DATA column");

	  /* Create DATAERR column */
	  const char * dataerr = GRAVI_DATAERR[region*npol + pol];
	  cpl_table_new_column_array (output_table, dataerr, CPL_TYPE_DOUBLE, ny_out);
	  cpl_array ** tDataErr = cpl_table_get_data_array (output_table, dataerr);
	  CPLCHECK_NUL ("Cannot create DATAERR column");
	  	  
	  /* Compute the flux in in [e] by loop on row and spectral direction 
	   * - Desinterlace the FT frames 
       * - data     = data - mean_sky
       * - var_data = data - mean_sky + var_sky
       * If 2 pixels per element
       * - data     = data1 - mean_sky1 + data2 - mean_sky2
       * - var_data = data + var_sky1 + var_sky2
       * Only considere the pixel that are not flag as BAD
	   * Allocating the arrays wrap data takes most of
	   * the time of the extract (80%) */
	  
	  for (cpl_size row = 0; row < nrow; row ++) {
		double *pData    = cpl_malloc (ny_out * sizeof(double));
		double *pDataErr = cpl_malloc (ny_out * sizeof(double));
		for (cpl_size j = 0; j < ny_out; j++) {
		  long idx = sizex * (j + ny*pol) + region*nx;
		  double value = 0;
		  double var = 0;
          if (cpl_array_get (arr_badft[0], idx, NULL) == 0) { // if not bad
              value += cpl_array_get (arr_data[row], idx, NULL) / gain - pSky[idx];
              var += pSkyVar[idx];
          }
		  if (nx>1 && cpl_array_get (arr_badft[0], idx+1, NULL) == 0) { // if not bad and if there is 2 pixels
		      value += cpl_array_get (arr_data[row], idx+1, NULL) / gain - pSky[idx+1];
		      var += pSkyVar[idx+1];
		  }
		  pData[j]   = value;
		  pDataErr[j] = sqrt (CPL_MAX (value,0.0) + var);
		}
		tData[row]    = cpl_array_wrap_double (pData, ny_out);
		tDataErr[row] = cpl_array_wrap_double (pDataErr, ny_out);
		
		CPLCHECK_NUL ("Cannot extract and wrap data");
	  } /* End loop on rows */
	  
	  CPLCHECK_NUL ("Cannot compute region");
	}
	/* End loop on polarisation */
  }
  /* End loop on regions */

  FREE (cpl_free, pSky);
  FREE (cpl_free, pSkyVar);
  
  CPLCHECK_NUL ("Cannot format ft");
  
  gravi_msg_function_exit(0);
  return output_table;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Extract the SC spectrum with profile
 *
 * @param profile_table     Input table containing the profile images
 * @param raw_imagelist     Input imagelist of raw data
 * @param rawVar_imagelist  Input imagelist of variance
 * @param startx            Input left location of profile in image
 *
 * Extract all the spectrum using the spectrum extraction method
 * based on the optimal extraction algorithm (Horne, 1986). This
 * algorithm is based on the knowledge of a spatial profile. Note
 * the profile may not be optimal (ex: boxcard) but shall ensure 
 * flux-conservation from its normalization.
 *
 * The profile table shall contains DATA# columns, each containing
 * the image of the profile of this region. The return tables contains
 * DATA# column with flux, and DATAERR# columns with errors
 * (sqrt(variance)).
 *
 * startx is the left start column of the profile in images, in FITS
 * convention (1 for the first pixel).
 */
/*----------------------------------------------------------------------------*/

cpl_table * gravi_imglist_sc_collapse (cpl_table * profile_table,
                                       cpl_imagelist * raw_imglist,
                                       cpl_imagelist * rawVar_imglist,
                                       cpl_size startx)
{
    int nv;
    gravi_msg_function_start(1);

    cpl_ensure (profile_table,  CPL_ERROR_NULL_INPUT, NULL);
    cpl_ensure (raw_imglist,    CPL_ERROR_NULL_INPUT, NULL);
    cpl_ensure (rawVar_imglist, CPL_ERROR_NULL_INPUT, NULL);

    /* Get data */
    cpl_size n_region = gravi_spectrum_get_nregion (profile_table);
    cpl_size n_row = cpl_imagelist_get_size (raw_imglist);
    cpl_ensure (n_row > 0, CPL_ERROR_ILLEGAL_INPUT, NULL);
    cpl_ensure (n_region==24 || n_region==48, CPL_ERROR_ILLEGAL_INPUT, NULL);

    /* Ensure the input images are of type DOUBLE */
    cpl_type type0 = cpl_image_get_type (cpl_imagelist_get (raw_imglist, 0));
    cpl_type type1 = cpl_image_get_type (cpl_imagelist_get (rawVar_imglist, 0));
    cpl_ensure (type0 == CPL_TYPE_DOUBLE, CPL_ERROR_ILLEGAL_INPUT, NULL);
    cpl_ensure (type1 == CPL_TYPE_DOUBLE, CPL_ERROR_ILLEGAL_INPUT, NULL);
    
    /* Get output tables */
    cpl_table * spectrum_table = cpl_table_new (n_row);

	/* Get the profile dimension */
	cpl_size nx = cpl_table_get_column_dimension (profile_table,"DATA1",0);
	cpl_size ny = cpl_table_get_column_dimension (profile_table,"DATA1",1);
    CPLCHECK_NUL ("Cannot get profile dimension");
    
	/* Loop on regions (output of beam combiner) */
	for (cpl_size region = 0; region < n_region; region++){

	    /* Verbose every 6 regions */
	    if ( !region || !((region+1)%6) )
		  cpl_msg_info_overwritable(cpl_func, "Extract region of SC %lld over %lld",region+1,n_region);

		/* Construction of the DATA column */
		const char * regionNameData = GRAVI_DATA[region];
		cpl_table_new_column_array (spectrum_table, regionNameData, CPL_TYPE_DOUBLE, nx);
		
		/* Construction of the DATAERR column */
		const char * regionNameErr = GRAVI_DATAERR[region];
		cpl_table_new_column_array (spectrum_table, regionNameErr, CPL_TYPE_DOUBLE, nx);

		/* Get the profile of this region as a 2D image */
		cpl_imagelist * profile_imglist = gravi_imagelist_wrap_column (profile_table, GRAVI_DATA[region]);
		cpl_image * profile_img = cpl_imagelist_get (profile_imglist, 0);
		CPLCHECK_NUL ("Cannot get data");

		/* Get the bound of the profile of this region.
		 * x = spectral, y = spatial */
		int xmin = 1,  xmax = nx;
		int ymin = ny, ymax = 0;
		for (cpl_size jy = 1 ; jy <= ny ; jy++ ) {
            for (cpl_size ix = 1 ; ix <= nx ; ix++ ) {
                if ( !cpl_image_get (profile_img, ix, jy, &nv) ) continue;
                if ( jy<ymin ) ymin=jy;
                if ( jy>ymax ) ymax=jy;
            }
		}
		CPLCHECK_NUL ("Cannot get profile limits");

		/* Dump the found limits */
		cpl_msg_debug (cpl_func, "Found limits: x=[%4d, %4d] y=[%4d, %4d]  (nx=%lld,ny=%lld, FITS convention)",
                       xmin,xmax,ymin,ymax,nx,ny);

		/* Extract a cropped version of the profile */
		cpl_image * profile_crop = cpl_image_extract (profile_img, xmin, ymin, xmax, ymax);
        CPLCHECK_NUL ("Cannot extract profile");
		
		/* Get the pointer to rows, to avoid calling this into the loop */
		cpl_array ** tData    = cpl_table_get_data_array (spectrum_table, regionNameData);
		cpl_array ** tDataErr = cpl_table_get_data_array (spectrum_table, regionNameErr);

        /* Loop on frame */
		for (cpl_size row = 0; row < n_row; row ++) {
		  
		  /* Extracted flux in [e] for this frame 
		   * rawFlux = < image * profile > */
		  cpl_image * data = cpl_image_extract (cpl_imagelist_get (raw_imglist, row),
                                                            xmin+startx-1, ymin, xmax+startx-1, ymax);
		  cpl_image * rawFlux_profiled = cpl_image_multiply_create (data, profile_crop);
          
          
		  /* Extracted variance for this frame in [e] 
		   * rawVar =  < variance * profile^2 > */
		  cpl_image * variance = cpl_image_extract (cpl_imagelist_get (rawVar_imglist, row),
                                                           xmin+startx-1, ymin, xmax+startx-1, ymax);

          cpl_image_divide (rawFlux_profiled,variance);
                                            
		  cpl_image * rawVar_profiled = cpl_image_divide_create (profile_crop, variance);
		  cpl_image_multiply (rawVar_profiled, profile_crop);
		  cpl_image *rawErr = cpl_image_collapse_create (rawVar_profiled,0);
		  cpl_image_threshold (rawErr, 0.0, DBL_MAX, 0.0, DBL_MAX);
          CPLCHECK_NUL ("Cannot collapse variance");

		  cpl_image *rawFlux = cpl_image_collapse_create (rawFlux_profiled,0);
          cpl_image_divide (rawFlux,rawErr);
          CPLCHECK_NUL ("Cannot collapse flux");

		  cpl_image_power (rawErr, -0.5);

          cpl_image_delete (rawFlux_profiled);
		  cpl_image_delete (rawVar_profiled);

        

		  /* Fill the output table for the given region and frame : flux in [e] */
		  tData[row] = cpl_array_wrap_double (cpl_image_get_data_double (rawFlux), nx);

		  /* Fill the output table for the given region and frame : error in [e] */
		  tDataErr[row] = cpl_array_wrap_double (cpl_image_get_data_double (rawErr), nx);

		  /* Delete tmp images */
          cpl_image_delete (data);
          cpl_image_delete (variance);
		  FREE (cpl_image_unwrap, rawFlux);
		  FREE (cpl_image_unwrap, rawErr);
		}

 		/* Delete all arrays of the loop */
        FREE (gravi_imagelist_unwrap_images, profile_imglist);
		FREE (cpl_image_delete, profile_crop);
	}
	/* End loop on region */
    
    gravi_msg_function_exit(1);
    return spectrum_table;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Extract the SC spectrum with profile
 *
 * @param profile_table     Input table containing the profile images
 * @param raw_imagelist     Input imagelist of raw data
 * @param rawVar_imagelist  Input imagelist of variance
 * @param startx            Input left location of profile in image
 *
 * Extract all the spectrum using the spectrum extraction method
 * based on the optimal extraction algorithm (Horne, 1986). This
 * algorithm is based on the knowledge of a spatial profile. Note
 * the profile may not be optimal (ex: boxcard) but shall ensure 
 * flux-conservation from its normalization.
 *
 * The profile table shall contains DATA# columns, each containing
 * the image of the profile of this region. The return tables contains
 * DATA# column with flux, and DATAERR# columns with errors
 * (sqrt(variance)).
 *
 * startx is the left start column of the profile in images, in FITS
 * convention (1 for the first pixel).
 */
/*----------------------------------------------------------------------------*/

cpl_table * gravi_imglist_sc_collapse_robust (cpl_table * profile_table,
                                       cpl_imagelist * raw_imglist,
                                       cpl_imagelist * rawVar_imglist,
                                       cpl_size startx)
{
    int nv;
    gravi_msg_function_start(1);
    const int NB_ROBUST_ITERATION = 3;
    cpl_ensure (profile_table,  CPL_ERROR_NULL_INPUT, NULL);
    cpl_ensure (raw_imglist,    CPL_ERROR_NULL_INPUT, NULL);
    cpl_ensure (rawVar_imglist, CPL_ERROR_NULL_INPUT, NULL);

    /* Get data */
    cpl_size n_region = gravi_spectrum_get_nregion (profile_table);
    cpl_size n_frame = cpl_imagelist_get_size (raw_imglist);
    cpl_ensure (n_frame > 0, CPL_ERROR_ILLEGAL_INPUT, NULL);
    cpl_ensure (n_region==24 || n_region==48, CPL_ERROR_ILLEGAL_INPUT, NULL);

    /* Ensure the input images are of type DOUBLE */
    cpl_type type0 = cpl_image_get_type (cpl_imagelist_get (raw_imglist, 0));
    cpl_type type1 = cpl_image_get_type (cpl_imagelist_get (rawVar_imglist, 0));
    cpl_ensure (type0 == CPL_TYPE_DOUBLE, CPL_ERROR_ILLEGAL_INPUT, NULL);
    cpl_ensure (type1 == CPL_TYPE_DOUBLE, CPL_ERROR_ILLEGAL_INPUT, NULL);
    
    /* Get output tables */
    cpl_table * spectrum_table = cpl_table_new (n_frame);

	/* Get the profile dimension */
	cpl_size nx = cpl_table_get_column_dimension (profile_table,"DATA1",0);
	cpl_size ny = cpl_table_get_column_dimension (profile_table,"DATA1",1);
    CPLCHECK_NUL ("Cannot get profile dimension");
    
	/* Loop on regions (output of beam combiner) */
	for (cpl_size region = 0; region < n_region; region++){

	    /* Verbose every 6 regions */
	    if ( !region || !((region+1)%6) )
		  cpl_msg_info_overwritable(cpl_func, "Extract region of SC %lld over %lld",region+1,n_region);

		/* Construction of the DATA column */
		const char * regionNameData = GRAVI_DATA[region];
		cpl_table_new_column_array (spectrum_table, regionNameData, CPL_TYPE_DOUBLE, nx);
		
		/* Construction of the DATAERR column */
		const char * regionNameErr = GRAVI_DATAERR[region];
		cpl_table_new_column_array (spectrum_table, regionNameErr, CPL_TYPE_DOUBLE, nx);

		/* Get the profile of this region as a 2D image */
		cpl_imagelist * profile_imglist = gravi_imagelist_wrap_column (profile_table, GRAVI_DATA[region]);
		cpl_image * profile_img = cpl_imagelist_get (profile_imglist, 0);
		CPLCHECK_NUL ("Cannot get data");

		/* Get the bound of the profile of this region.
		 * x = spectral, y = spatial */
		int xmin = 1,  xmax = nx;
		int ymin = ny, ymax = 0;
		for (cpl_size jy = 1 ; jy <= ny ; jy++ ) {
            for (cpl_size ix = 1 ; ix <= nx ; ix++ ) {
                if ( !cpl_image_get (profile_img, ix, jy, &nv) ) continue;
                if ( jy<ymin ) ymin=jy;
                if ( jy>ymax ) ymax=jy;
            }
		}
		CPLCHECK_NUL ("Cannot get profile limits");

		/* Dump the found limits */
		cpl_msg_debug (cpl_func, "Found limits: x=[%4d, %4d] y=[%4d, %4d]  (nx=%lld,ny=%lld, FITS convention)",
                       xmin,xmax,ymin,ymax,nx,ny);

		/* Extract a cropped version of the profile */
		cpl_image * profile_crop = cpl_image_extract (profile_img, xmin, ymin, xmax, ymax);
        cpl_image * squared_profile_crop = cpl_image_power_create (profile_crop, 2);

        CPLCHECK_NUL ("Cannot extract profile");
		
		/* Get the pointer to frames, to avoid calling this into the loop */
		cpl_array ** tData    = cpl_table_get_data_array (spectrum_table, regionNameData);
		cpl_array ** tDataErr = cpl_table_get_data_array (spectrum_table, regionNameErr);

        /* Loop on frame */
		for (cpl_size frame = 0; frame < n_frame; frame ++) {
		    /* Extracted flux in [e] for this frame 
		     * rawFlux = < image * profile > */
		    cpl_image * data = cpl_image_extract (cpl_imagelist_get (raw_imglist, frame),
                                                              xmin+startx-1, ymin, xmax+startx-1, ymax);
		    cpl_image * rawFlux_profiled = cpl_image_multiply_create (data, profile_crop);
          
		    /* Extracted variance for this frame in [e] 
		     * rawVar =  < variance * profile^2 > */
		    cpl_image * variance = cpl_image_extract (cpl_imagelist_get (rawVar_imglist, frame),
                                                           xmin+startx-1, ymin, xmax+startx-1, ymax);

            cpl_image * weighted_flux_profile = cpl_image_divide_create (rawFlux_profiled,variance);
                                            
		    cpl_image * rawVar_profiled = cpl_image_divide_create (squared_profile_crop, variance);
		    cpl_image *rawErr = cpl_image_collapse_create (rawVar_profiled,0);
		    cpl_image_threshold (rawErr, 0.0, DBL_MAX, 0.0, DBL_MAX);
		    cpl_image_delete (rawVar_profiled);
            CPLCHECK_NUL ("Cannot collapse variance");

		    cpl_image *rawFlux = cpl_image_collapse_create (weighted_flux_profile,0);
		    cpl_image_delete (weighted_flux_profile);
            cpl_image_divide (rawFlux,rawErr);
            CPLCHECK_NUL ("Cannot collapse flux");

		    cpl_image_power (rawErr, -0.5);


        
            for (int n_it = 0; n_it < NB_ROBUST_ITERATION; n_it++){
                int numbad =0;
                double mad = 0;
                cpl_image * residuals = cpl_image_duplicate (data);
                double * residuals_values = cpl_image_get_data_double(residuals);

                int ncol = cpl_image_get_size_x(residuals);
                int nrow = cpl_image_get_size_y(residuals);
                for ( int col = 0; col < ncol; col++) {
                    double flux  = cpl_image_get ( rawFlux,col+1,1,&nv);
                    for ( int row = 0; row < nrow; row++) {
                        double model  = flux * cpl_image_get ( profile_crop, col+1, row+1, &nv);
                        double precision = cpl_image_get ( variance, col+1, row+1, &nv);
                        if (nv==1)
                            precision = 0.;
                        else
                            precision = precision <= 0 ? 0 : 1.0 / (precision);
                        residuals_values[col + row * ncol]  =   (residuals_values[col + row * ncol] - model)* sqrt(precision);
                    }
                }
                double med = cpl_image_get_mad(residuals, &mad);


                /* Verbose every 6 regions */
	            if ( !region || !((region+1)%12) ){
		            cpl_msg_info(cpl_func, "med %g",med);
		            cpl_msg_info(cpl_func, "mad %g",mad*CPL_MATH_STD_MAD);
                }

                for ( int col = 0; col < ncol; col++) {
                    for ( int row = 0; row < nrow; row++) {
                        if (mad==0){
                            residuals_values[col + row * ncol] = 1;
                        }else{
                        double precision = cpl_image_get ( variance, col+1, row+1, &nv);
                        if (nv==1)
                            precision = 0.;
                        else
                            precision = precision <= 0 ? 0 : 1.0 / (precision);
                        double r = residuals_values[col + row * ncol]/ (mad*CPL_MATH_STD_MAD) ; 
                        residuals_values[col + row * ncol] = precision  / ( 1.0 + pow( r  / (2.985) ,2.0 ) );  // weights
                        }
                    }
                }


                /* Verbose every 6 regions */
	            if ( !region || !((region+1)%12) ){
		            cpl_msg_info(cpl_func, "Weights mean %g",cpl_image_get_mean(residuals));
		            //cpl_msg_info(cpl_func, "mad %g",mad*CPL_MATH_STD_MAD);
                }

                cpl_image * weighted_flux_profile = cpl_image_multiply_create (rawFlux_profiled,residuals);
                CPLCHECK_NUL ("Cannot collapse variance 1");
                                           
		        cpl_image * rawVar_profiled = cpl_image_multiply_create (squared_profile_crop, residuals);
                CPLCHECK_NUL ("Cannot collapse variance 2");
    		    cpl_image * tmp_rawErr = cpl_image_collapse_create (rawVar_profiled,0);
                CPLCHECK_NUL ("Cannot collapse variance 3");
	    	    cpl_image_threshold (tmp_rawErr, 0.0, DBL_MAX, 0.0, DBL_MAX);
                CPLCHECK_NUL ("Cannot collapse variance 4");
		        cpl_image_delete (rawVar_profiled);
                CPLCHECK_NUL ("Cannot collapse variance 5");

    		    cpl_image *tmp_rawFlux = cpl_image_collapse_create (weighted_flux_profile,0);
	    	    cpl_image_delete (weighted_flux_profile);
                cpl_image_divide (tmp_rawFlux,tmp_rawErr);
                CPLCHECK_NUL ("Cannot collapse flux");

    		    cpl_image_power (tmp_rawErr, -0.5);

                cpl_image_copy(rawFlux,tmp_rawFlux,1,1);
                cpl_image_copy(rawErr,tmp_rawErr,1,1);
                
                cpl_image_delete (residuals);
                cpl_image_delete (tmp_rawFlux);
                cpl_image_delete (tmp_rawErr);
            }

		    /* Delete tmp images */
            cpl_image_delete (data);
            cpl_image_delete (variance);
            cpl_image_delete (rawFlux_profiled);

		    /* Fill the output table for the given region and frame : flux in [e] */
		    tData[frame] = cpl_array_wrap_double (cpl_image_get_data_double (rawFlux), nx);

		    /* Fill the output table for the given region and frame : error in [e] */
		    tDataErr[frame] = cpl_array_wrap_double (cpl_image_get_data_double (rawErr), nx);

		    FREE (cpl_image_unwrap, rawFlux);
		    FREE (cpl_image_unwrap, rawErr);
		}

 		/* Delete all arrays of the loop */
        FREE (gravi_imagelist_unwrap_images, profile_imglist);
		FREE (cpl_image_delete, squared_profile_crop);
		FREE (cpl_image_delete, profile_crop);
	}
	/* End loop on region */
    
    gravi_msg_function_exit(1);
    return spectrum_table;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Create the SPECTRUM gravi_data with extracted spectrum per region
 * 
 * @param raw_data         Input RAW gravi_data			
 * @param profile_data     FLAT calibration map, with profiles
 * @param dark_data        DARK calibration map
 * @param bad_map          BAD calibration map
 * @param sky_map          SKY calibration map
 * @param parlist          Input parameter list with :
 *                         - ditshift-sc : Shift the time of SC DITs by an integer
 *                         value to account for lost frames in exposure (issue on the
 *                         instrument side, report to instrument team). The
 *                         time of all DITs in exposure are increased by
 *                         ditshift x PERIOD. ditshift can be 0, positive
 *                         (system has lost one SC DIT), or negative
 *                         (SC desynchronized).
 * @param det_type         The detector to extract GRAVI_DET_SC, GRAVI_DET_FT or
 *                         GRAVI_DET_ALL
 *
 * @return The output gravi data contenning all the spectrums
 *
 * \exception CPL_ERROR_NULL_INPUT input data is missing
 *
 * It substrates the dark map. It converts data into photoelectrons using
 * the gain map and identify the bad pixels for a correction. Finally it extracts
 * the spectra with the calibrated spatial profile (for SC data) 
 */
/*----------------------------------------------------------------------------*/

gravi_data * gravi_extract_spectrum (gravi_data * raw_data,
									 gravi_data * profile_map,
									 gravi_data * dark_map,
									 gravi_data * bad_map,
									 gravi_data * sky_map,
                                     const cpl_parameterlist * parlist,
                                     enum gravi_detector_type det_type)
{
    gravi_msg_function_start(1);
	cpl_ensure (raw_data,     CPL_ERROR_NULL_INPUT, NULL);
	cpl_ensure (profile_map,  CPL_ERROR_NULL_INPUT, NULL);
	cpl_ensure (bad_map,      CPL_ERROR_NULL_INPUT, NULL);
	cpl_ensure (dark_map || sky_map, CPL_ERROR_NULL_INPUT, NULL);
    cpl_ensure (parlist,      CPL_ERROR_NULL_INPUT, NULL);

    /* Get header of input gravi_data */
    cpl_propertylist * raw_header = gravi_data_get_header (raw_data);
	/* cpl_propertylist * profile_header = gravi_data_get_header (profile_map); */

    /* Create output gravi_data */
	gravi_data * spectrum_data = gravi_data_new(0);
	cpl_propertylist * spectrum_header = gravi_data_get_header (spectrum_data);

    /* Dump all header from RAW product */
    cpl_propertylist_append (spectrum_header, raw_header);
    
	if ((det_type == GRAVI_DET_FT || det_type == GRAVI_DET_ALL) && 
	        gravi_data_has_detector (raw_data, GRAVI_FT)) {

		/*
		 *  Compute the flux and variance of each spectral element of FT
		 *  The spectrum are already extracted from spatial direction.
		 */

        /* Dupplicate necessary extension in the output spectrum_data */
        gravi_data_copy_ext (spectrum_data, raw_data, GRAVI_IMAGING_DETECTOR_FT_EXT);
        CPLCHECK_NUL ("Cannot duplicate necessary extensions");

        /* Get the data */
		cpl_table * imaging_data_ft = gravi_data_get_table (raw_data, GRAVI_IMAGING_DATA_FT_EXT);
		cpl_size n_region = cpl_table_get_nrow (gravi_data_get_table (raw_data,GRAVI_IMAGING_DETECTOR_FT_EXT));

		/* Get the background mean and std. In case we have both SKY and DARK we still use
		 * only the SKY because it is surely the closest in time (taken at night) and with the same setup
		 * FIXME: could be interesting to check the setup and time distance... */
        cpl_table * skyavg_table = NULL, * skystd_table = NULL;
		if (sky_map != NULL) {
		  cpl_msg_info (cpl_func, "Extract FT spectra with SKY as background and variance");
		  skyavg_table = gravi_data_get_table (sky_map, GRAVI_IMAGING_DATA_FT_EXT);
		  skystd_table = gravi_data_get_table (sky_map, GRAVI_IMAGING_ERR_FT_EXT);
		}
		else if (dark_map != NULL) {
		  if ( !gravi_data_is_internal(raw_data) )
			gravi_pfits_add_check (spectrum_header, "Extract FT spectra with DARK as background and variance");
		  else
			cpl_msg_info (cpl_func, "Extract FT spectra with DARK as background and variance");
		  skyavg_table = gravi_data_get_table (dark_map, GRAVI_IMAGING_DATA_FT_EXT);
		  skystd_table = gravi_data_get_table (dark_map, GRAVI_IMAGING_ERR_FT_EXT);
		} else {
		  cpl_error_set_message (cpl_func,CPL_ERROR_NULL_INPUT,
					 "DARK or SKY should be provided");
		  return NULL;
		}

		/* Get the FT gain in [ADU/e] */
		double gain_ft = gravi_pfits_get_ft_gain (raw_header);
		cpl_propertylist_update_double (spectrum_header, "ESO QC USEDGAIN FT", gain_ft);
		cpl_propertylist_set_comment (spectrum_header, "ESO QC USEDGAIN FT", "[adu/e-] value used for reduction");

		CPLCHECK_NUL("Problem while getting the tables");

        /* Get the badpixel image */
        cpl_table * badft_table = gravi_data_get_table(bad_map, GRAVI_IMAGING_DATA_FT_EXT);

		/* Convert PIX column to DATA# and DATAERR# */
		cpl_table * spectrum_ft;
        spectrum_ft = gravi_table_ft_format (imaging_data_ft, skystd_table, skyavg_table, badft_table, n_region, gain_ft, parlist);
        CPLCHECK_NUL ("Cannot format FT data");

        /* Set units */
        for (cpl_size reg=0; reg<n_region; reg++) {
            cpl_table_set_column_unit (spectrum_ft, GRAVI_DATA[reg], "e");
            cpl_table_set_column_unit (spectrum_ft, GRAVI_DATAERR[reg], "e");
        }
        
		/* Set the data in output array */
		cpl_propertylist * spectrum_plist = cpl_propertylist_new ();
        int nx = gravi_spectrum_get_nwave (spectrum_ft);
		cpl_propertylist_update_int (spectrum_plist, PROFILE_NX, nx);
		cpl_propertylist_update_int (spectrum_plist, PROFILE_STARTX, 1);
		cpl_propertylist_update_int (spectrum_plist, PROFILE_FULLSTARTX, 0);
        
		gravi_data_add_table (spectrum_data, spectrum_plist, GRAVI_SPECTRUM_DATA_FT_EXT, spectrum_ft);
	}

        
    if ((det_type == GRAVI_DET_SC || det_type == GRAVI_DET_ALL) && 
            gravi_data_has_detector (raw_data, GRAVI_SC)) {
        
        /*
         *  Compute the flux and variance of
         *  each spectral element of SC
         */
        
        /* Dupplicate necessary extension in the output spectrum_data */
        gravi_data_copy_ext (spectrum_data, raw_data, GRAVI_IMAGING_DETECTOR_SC_EXT);
        CPLCHECK_NUL ("Cannot duplicate necessary extensions");

        /* Get the profile data */
        cpl_table * profile_table = gravi_data_get_table (profile_map, GRAVI_PROFILE_DATA_EXT);
        cpl_ensure (profile_table, CPL_ERROR_ILLEGAL_INPUT, NULL);
        
        /* Get the data */
        cpl_imagelist * imaging_data = gravi_data_get_cube (raw_data, GRAVI_IMAGING_DATA_SC_EXT);
        cpl_propertylist * profile_plist = gravi_data_get_plist (profile_map, GRAVI_PROFILE_DATA_EXT);
        CPLCHECK_NUL ("Cannot get data");

        /* Get the SC gain in [ADU/e] */
        double gain_sc = gravi_pfits_get_sc_gain (raw_header);
        cpl_propertylist_update_double (spectrum_header, "ESO QC USEDGAIN SC", gain_sc);
        cpl_propertylist_set_comment (spectrum_header, "ESO QC USEDGAIN SC", "[adu/e-] value used for reduction");
        
        /* Get the badpixel image */
        cpl_image * badpix_img = gravi_data_get_img (bad_map, GRAVI_IMAGING_DATA_SC_EXT);
        
        /* Estimate the stellar flux in [e].
         * Give priority to SKY for estimating the background */
        cpl_image * skyavg_img;
        if (sky_map != NULL){
            cpl_propertylist * sky_header = gravi_data_get_header (sky_map);
            int sky_map_usefull=1;
            if (gravi_pfits_get_dit_sc (raw_header)!=gravi_pfits_get_dit_sc (sky_header))
            {
                cpl_msg_warning (cpl_func, "SC DIT is different from SKY SC DIT");
                sky_map_usefull=0;
            }
            if ( strcmp(gravi_pfits_get_spec_res (raw_header),gravi_pfits_get_spec_res (sky_header)) > 0)
            {
                cpl_msg_warning (cpl_func, "SC spectral resolution is different from SKY resolution :");
                cpl_msg_warning (cpl_func, "SC resolution is %s", gravi_pfits_get_spec_res (raw_header));
                cpl_msg_warning (cpl_func, "SKY resolution is %s", gravi_pfits_get_spec_res (sky_header));
                sky_map_usefull=0;
            }
            if ( strcmp(gravi_pfits_get_pola_mode (raw_header, GRAVI_DET_SC),gravi_pfits_get_pola_mode (sky_header, GRAVI_DET_SC)) > 0)
            {
                cpl_msg_warning (cpl_func, "SC POLA is different from SKY POLA :");
                cpl_msg_warning (cpl_func, "SC POLA is %s", gravi_pfits_get_pola_mode (raw_header, GRAVI_DET_SC));
                cpl_msg_warning (cpl_func, "SKY POLA is %s", gravi_pfits_get_pola_mode (sky_header, GRAVI_DET_SC));
                sky_map_usefull=0;
            }
            if ((sky_map_usefull==1)||(dark_map == NULL))
            {
                cpl_msg_info (cpl_func,"Extract SC spectra with SKY as background");
                skyavg_img = gravi_data_get_img (sky_map, GRAVI_IMAGING_DATA_SC_EXT);
            } else {
                if ( !gravi_data_is_internal(raw_data) )
                    gravi_pfits_add_check (spectrum_header, "Extract SC spectra with DARK as background");
                else
                    cpl_msg_info (cpl_func,"Extract SC spectra with DARK as background");
                skyavg_img = gravi_data_get_img (dark_map, GRAVI_IMAGING_DATA_SC_EXT);
            }
        } else if (dark_map != NULL) {
            if ( !gravi_data_is_internal(raw_data) )
                gravi_pfits_add_check (spectrum_header, "Extract SC spectra with DARK as background");
            else
                cpl_msg_info (cpl_func,"Extract SC spectra with DARK as background");
            skyavg_img = gravi_data_get_img (dark_map, GRAVI_IMAGING_DATA_SC_EXT);
        }
        
        /* Estimate the total variance in [e^2] composed of background variance 
         * and additional photonic variance. Give priority to DARK for estimating
         * the background variance */
        cpl_image * darkavg_img;
        cpl_image * darkstd_img;
        CPLCHECK_NUL ("Check for error to remove");
        if (dark_map != NULL) {
            cpl_msg_info (cpl_func,"Extract SC photonic variance with DARK as background and variance");
            darkavg_img = gravi_data_get_img (dark_map, GRAVI_IMAGING_DATA_SC_EXT);
            darkstd_img = gravi_data_get_img (dark_map, GRAVI_IMAGING_ERR_SC_EXT);
        }
        else if (sky_map != NULL) {
            gravi_pfits_add_check (spectrum_header, "Extract SC photonic variance with SKY as background and variance");
            darkavg_img = gravi_data_get_img (sky_map, GRAVI_IMAGING_DATA_SC_EXT);
            darkstd_img = gravi_data_get_img (sky_map, GRAVI_IMAGING_ERR_SC_EXT);
        }        

        /* Build the dimension for the raw_imglist, that shall match the profile dimension */
        cpl_size startx = gravi_pfits_get_startx (profile_plist);
        CPLCHECK_NUL ("The coordinate dimensions of the new window is missing");

        /* raw = (DATA - SKY) / gain   [e] */
        cpl_msg_info (cpl_func,"Compute flux image");
        cpl_imagelist * raw_imglist;
        raw_imglist = cpl_imagelist_duplicate (imaging_data);
        cpl_imagelist_subtract_image (raw_imglist, skyavg_img);
        cpl_imagelist_divide_scalar (raw_imglist, gain_sc);
        gravi_remove_badpixel_sc (raw_imglist, badpix_img);
        CPLCHECK_NUL ("Cannot extract the data");
        
        /* rawVar = (DATA - DARK) / gain    [e^2] */
        cpl_msg_info (cpl_func,"Compute variance image");
        cpl_imagelist * rawVar_imglist;
        rawVar_imglist = cpl_imagelist_duplicate (imaging_data);
        cpl_imagelist_subtract_image (rawVar_imglist, darkavg_img);
        cpl_imagelist_divide_scalar (rawVar_imglist, gain_sc);
        gravi_remove_badpixel_sc (rawVar_imglist, badpix_img);
        CPLCHECK_NUL ("Cannot extract the data");

        /* rawVar += (std(DARK) / gain) ** 2    [e^2] */
        cpl_image * darkvar_img = cpl_image_power_create (darkstd_img, 2.0);
        cpl_image_divide_scalar (darkvar_img, gain_sc * gain_sc);
        cpl_imagelist_add_image (rawVar_imglist, darkvar_img);
        FREE (cpl_image_delete, darkvar_img);
        CPLCHECK_NUL ("Cannot extract the data");

        /* Collapse images with spectrum and create
         * spectrum_table for SC */
        cpl_table * spectrum_table;
        spectrum_table = gravi_imglist_sc_collapse_robust (profile_table, raw_imglist,
                                                    rawVar_imglist, startx);
        CPLCHECK_NUL ("Cannot collapse the spectrum");

        cpl_size n_row    = cpl_table_get_nrow (spectrum_table);
        cpl_size n_region = cpl_table_get_ncol (profile_table);

        /* Set units */
        for (cpl_size reg=0; reg<n_region; reg++) {
            cpl_table_set_column_unit (spectrum_table, GRAVI_DATA[reg], "e");
            cpl_table_set_column_unit (spectrum_table, GRAVI_DATAERR[reg], "e");
        }

        /* Get the user-specified DIT shift if any */
        int ditshift = gravi_param_get_int_default (parlist, "gravity.preproc.ditshift-sc", 0);
        cpl_msg_info (cpl_func, "ditshift = %i", ditshift);

        if (ditshift !=0 )
            gravi_msg_warning ("CRITICAL","DITSHIFT is not set to zero... are you sure??");
        
        /* Compute the time of the middle of the SC DIT, in [us]
         * with respect to the RMN start PCR.ACQ.START */
        cpl_table_new_column (spectrum_table, "TIME", CPL_TYPE_INT);
        for (cpl_size row = 0; row < n_row; row ++) {
            cpl_table_set_int (spectrum_table, "TIME", row,
                               gravi_pfits_get_time_sc (raw_header, row + ditshift));
        }
	
        /* Check if the SC profil extraction is flux conservative */
        double full_flux_img = gravi_imagelist_get_flux (raw_imglist);
        double full_flux_reg = gravi_spectrum_get_flux (spectrum_table);
	
        cpl_msg_info (cpl_func, "Total flux in REGIONs: %.2e [e], in IMGs:%.2e [e]  (ratio=%.5f)",
                      full_flux_reg,full_flux_img,full_flux_reg/full_flux_img);
	
        cpl_propertylist_update_double (spectrum_header, "ESO QC TRANS PROFILE SC", full_flux_reg/full_flux_img);
        cpl_propertylist_set_comment (spectrum_header, "ESO QC TRANS PROFILE SC", "[e/e] at profile extraction");

        CPLCHECK_NUL ("Cannot verify if extraction was flux-conservative");

        /* Set the SPECTRUM_DATA_SC tables to the gravi_data */
        cpl_propertylist * spectrum_plist = cpl_propertylist_new ();
        cpl_propertylist_copy_property (spectrum_plist, profile_plist, PROFILE_FULLSTARTX);
        cpl_propertylist_copy_property (spectrum_plist, profile_plist, PROFILE_STARTX);
        cpl_propertylist_copy_property (spectrum_plist, profile_plist, PROFILE_NX);
        
        gravi_data_add_table (spectrum_data, spectrum_plist,
                              GRAVI_SPECTRUM_DATA_SC_EXT, spectrum_table);


        /* Apply the same extraction to the flat
         * FIXME: no units since it was not in e- */
        cpl_msg_info (cpl_func, "Extract the FLAT");
        
        cpl_table * spectrumflat_table;
        cpl_imagelist * flat_imglist = gravi_data_get_cube (profile_map, GRAVI_IMAGING_DATA_SC_EXT);
        spectrumflat_table = gravi_imglist_sc_collapse_robust (profile_table, flat_imglist,
                                                        flat_imglist, 1);
        CPLCHECK_NUL ("Cannot collapse the flat");

        /* Set the SPECTRUMFLAT_DATA_SC tables to the gravi_data */
        cpl_propertylist * spectrumflat_plist = cpl_propertylist_new ();
        gravi_data_add_table (spectrum_data, spectrumflat_plist,
                              GRAVI_SPECTRUMFLAT_DATA_SC_EXT, spectrumflat_table);

        /* Delete variable */
        FREE (cpl_imagelist_delete, raw_imglist);
        FREE (cpl_imagelist_delete, rawVar_imglist);
    } /* End SC case */

    gravi_msg_function_exit(1);
	return spectrum_data;
}



/*----------------------------------------------------------------------------*/
/**
 * @brief Substract metrology dark
 *
 * @param preproc_data     The table to manipulate
 * @param dark_data        DARK calibration map
 *
 * If Metrology dark is present, subtract DARK to METROLOGY data
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_subtract_met_dark (gravi_data * preproc_data,
                                                 gravi_data * dark_map)
{
    
    gravi_msg_function_start(1);
    cpl_ensure_code (dark_map , CPL_ERROR_NULL_INPUT);
    
    if (!gravi_data_has_extension (dark_map, GRAVI_METROLOGY_EXT)) {
        cpl_msg_warning (cpl_func,"The DARK data has no METROLOGY extension");
    }
    else
    {
        cpl_table * met_data_table = gravi_data_get_table (preproc_data, GRAVI_METROLOGY_EXT);
        cpl_table * met_dark_table = gravi_data_get_table (dark_map, GRAVI_METROLOGY_EXT);
        CPLCHECK_MSG("Problem while getting the tables");
        
        /* Get pointer to the data */
        const char * data_x="VOLT";
        cpl_size nrow  = cpl_table_get_nrow (met_data_table);
        cpl_array ** met_data_array = cpl_table_get_data_array (met_data_table, data_x);
        cpl_ensure_code (met_data_array, CPL_ERROR_ILLEGAL_INPUT);
        cpl_array ** met_dark_array = cpl_table_get_data_array (met_dark_table, data_x);
        cpl_ensure_code (met_dark_array, CPL_ERROR_ILLEGAL_INPUT);
        
        /* Subtract data */
        
        for (cpl_size j = 0; j < nrow ; j++)
        {
            cpl_array_subtract(met_data_array[j],met_dark_array[0]);
        }
        
        cpl_msg_info (cpl_func,"METROLOGY Volts have been subtracted from dark");
    }
 
    gravi_msg_function_exit(1);
    return CPL_ERROR_NONE;
    
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Re-interpolate in-place a spectrum table
 *
 * @param spectrum_table    The table to manipulate
 * @param wave_table        The input wave_table (one map per region)
 * @param oiwave_tables     The output OI_WAVEs table (same for all regions)
 * @param detector_rable    The input IMAGING DETECTOR table
 * @param specflat_table    The SPECTRUMFLAT, to allow unbiased effective wave
 * @param parlist           Input parameter list (no parameter used)
 *
 * All regions are re-interpolated into a common wavelength grid, defined
 * by the oiwave_table (OIFITS format). The input wave_table contains the 
 * the initial wavelength grid of each region.
 *
 * The wavelength bins in wave_table shall be monotonic incresing.
 * All wavelength bins of oiwave_table shall be within the min/max of 
 * wave_table (not extrapolation is allowed).
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_interpolate_spectrum_table (cpl_table * spectrum_table,
                                                 cpl_table * wave_table,
                                                 cpl_table ** oiwave_tables,
                                                 cpl_table * detector_table,
                                                 cpl_table * specflat_table,
                                                 int type_data)
{
    gravi_msg_function_start(1);

	/* Check needed data */
	cpl_ensure_code (spectrum_table,   CPL_ERROR_NULL_INPUT);
	cpl_ensure_code (wave_table,       CPL_ERROR_NULL_INPUT);
	cpl_ensure_code (oiwave_tables,    CPL_ERROR_NULL_INPUT);
	cpl_ensure_code (oiwave_tables[0], CPL_ERROR_NULL_INPUT);

    /* Sizes */
    cpl_size nb_region = gravi_spectrum_get_nregion (spectrum_table);
    cpl_size nb_pol    = gravi_spectrum_get_npol (spectrum_table);
    cpl_size nb_wave   = gravi_spectrum_get_nwave (spectrum_table);
    cpl_size nb_row    = cpl_table_get_nrow (spectrum_table);
    cpl_size nb_oiwave = cpl_table_get_nrow (oiwave_tables[0]);

    /* Ensure the two OI_WAVELENGTH have same size */
    if (nb_pol > 1) {
        cpl_ensure_code (oiwave_tables[1], CPL_ERROR_NULL_INPUT);
        cpl_ensure_code (cpl_table_get_nrow (oiwave_tables[0]) ==
                         cpl_table_get_nrow (oiwave_tables[1]),
                         CPL_ERROR_ILLEGAL_INPUT);
    }

    /* Only 24 and 48 regions are accepted */
	cpl_ensure_code (nb_region == 24 || nb_region == 48,
                     CPL_ERROR_ILLEGAL_INPUT);

    /* So far oversampling is not supported */
    cpl_msg_info (cpl_func,"nb_oiwave=%lld, nb_wave=%lld", nb_oiwave, nb_wave);
    cpl_ensure_code (nb_oiwave <= nb_wave, CPL_ERROR_ILLEGAL_INPUT);

    /* Verbose about the target effective map */
    if (specflat_table)
        cpl_msg_info (cpl_func, "Change the target wavelength because of FLAT");
    else
        cpl_msg_info (cpl_func, "Don't change the target wavelenght because of FLAT");
    

    /* Allocate memory for the temporary target wavelength (xout) 
     * Allocate memory for the temporary output (yout) */
    double * xout = cpl_malloc (nb_oiwave * sizeof(double));
    double * yout = cpl_malloc (nb_oiwave * sizeof(double));
    cpl_size * id     = cpl_malloc (nb_oiwave * sizeof(cpl_size));

    double * outputData = cpl_malloc (nb_oiwave * sizeof(double));
    double * outputErr = cpl_malloc (nb_oiwave * sizeof(double));

    /* Loop on polarisations */
    for (int pol = 0; pol < nb_pol; pol++) {

        /* Get max and min of the new wavelength table */
        double oiwave_min = cpl_table_get_column_min (oiwave_tables[pol], "EFF_WAVE");
        double oiwave_max = cpl_table_get_column_max (oiwave_tables[pol], "EFF_WAVE");
        
        cpl_msg_info(cpl_func,"Reinterpolate pol %i: oiwave_min=%g, oiwave_max=%g [um]",
                     pol, oiwave_min*1e6, oiwave_max*1e6);
        
        /* Loop on regions -- skip regions not in requested polarisation */
        for (cpl_size reg = 0; reg < nb_region; reg++){
            if (gravi_region_get_pol (detector_table, reg) != pol) continue;
            
            /* Verbose every regions */

            /* Name and dimension of this region */
            const char * data_x = GRAVI_DATA[reg];
            const char * data_errx = GRAVI_DATAERR[reg];

            /* Check size of SPECTRUM of this region */
            int nbd = cpl_table_get_column_depth (spectrum_table, data_x);
            cpl_ensure_code (nbd == nb_wave, CPL_ERROR_ILLEGAL_INPUT);
            
            /* Check size of WAVE of this region */
            int nbw = cpl_table_get_column_depth (wave_table, data_x);
            if (nbw != nb_wave)
                cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                  "Wavelength calibration is incompatible with the data. Make sure the \n"
                                  "calibration has been produced by the same version of the pipeline");
            cpl_ensure_code (nbw == nb_wave, CPL_ERROR_ILLEGAL_INPUT);

            /* Ensure we have no extrapolation and have increasing wavelengths */
            double wave_first = gravi_table_get_value (wave_table, data_x, 0, 0);
            double wave_last  = gravi_table_get_value (wave_table, data_x, 0, nb_wave-1);
            cpl_ensure_code (wave_first < wave_last,  CPL_ERROR_ILLEGAL_INPUT);
            //cpl_ensure_code (oiwave_min > wave_first, CPL_ERROR_ILLEGAL_INPUT);
            //cpl_ensure_code (oiwave_max < wave_last,  CPL_ERROR_ILLEGAL_INPUT);
            
//            /* Get directly the pointer to arrays to speed-up */
//            cpl_array ** inputData = cpl_table_get_data_array (spectrum_table, data_x);
//            cpl_array ** inputErr  = cpl_table_get_data_array (spectrum_table, data_errx);
//            CPLCHECK_MSG("Cannot get data");
            

            /* Get directly the pointer to arrays to speed-up */
            cpl_array ** inputData = cpl_table_get_data_array (spectrum_table, data_x);
            cpl_array ** inputErr  = cpl_table_get_data_array (spectrum_table, data_errx);
            const cpl_array * flat = NULL;
            if (specflat_table){
                flat = cpl_table_get_array (specflat_table, data_x, 0);
            }
            cpl_array * wave = cpl_table_get_data_array(wave_table, data_x)[0];
            CPLCHECK_MSG("Cannot get data");

            /* Verify the data are double since we will get pointer to them */
            cpl_ensure_code ((cpl_array_get_type (inputData[0]) == CPL_TYPE_DOUBLE) &&
                             (cpl_array_get_type (inputErr[0])  == CPL_TYPE_DOUBLE) &&
                             (cpl_array_get_type (wave)  == CPL_TYPE_DOUBLE),
                             CPL_ERROR_ILLEGAL_INPUT);

            
            /* Get the current wavelength as pointer */
            double * wave_ref = cpl_array_get_data_double (wave);
            CPLCHECK_MSG("Cannot get data");
            

            /* Check the Option
             * Go to 3 pixels interp if :
             *  1- if LOW ie. nwave < 50 (old: Option set)
             *  2- Flat is given in input
             *  3- Not FT data
             */
            //if ( (gravi_param_get_bool (parlist,"gravity.preproc.interp-3pixels") == TRUE)
            //        && (flat  != NULL) && ( nb_oiwave !=5 ))
            if ( ( nb_oiwave <= 50 ) && ( type_data == GRAVI_SC) && (flat  != NULL) ) {
                cpl_msg_info_overwritable (cpl_func, "Reinterpolate (3 pixels) region %lld "
                                                       "over %lld (%lld->%lld channels)",
                                                       reg+1,nb_region,nb_wave,nb_oiwave);
                /* Init the index and weights for the interpolation such that
                 * yout[iw] =   inputData[id[iw]-1] * weight_m[iw]
                 *            + inputData[id[iw]]   * weight_m[iw]
                 *            + inputData[id[iw]+1] * weight_p[iw]
                 * yerr[iw] =  sqrt( (inputErr[id[iw]-1] * weight_m[iw])^2
                 *            + (inputErr[id[iw]]   * weight_m[iw])^2
                 *            + (inputErr[id[iw]+1] * weight_p[iw])^2 )
                 * We here assume monotonicity and no extrapolation */
                double * weight_m  = cpl_malloc (nb_oiwave * sizeof(double));
                double * weight_i  = cpl_malloc (nb_oiwave * sizeof(double));
                double * weight_p  = cpl_malloc (nb_oiwave * sizeof(double));

                for (cpl_size iw = 0 ; iw < nb_oiwave ; iw ++)
                {
                	// xo is the wavenumber of the targeted wavelength (of indice iw)
                	double l0=cpl_table_get (oiwave_tables[pol], "EFF_WAVE", iw, NULL);

                	// id[iw] is the index of the closest wavenumber to xo
                	// id[iw] must be above 1, and below nb_wave-1
                	// wave is assumed to be in increasing order

                	if (iw == 0) id[iw] = 1;
                	else id[iw] = id[iw-1];
                	while ( (1./l0 < 0.5/wave_ref[id[iw]]+0.5/wave_ref[id[iw]+1]) && ( id[iw] < nb_wave-2 ) ){
                		id[iw]++;
                	}

                	// compute the wavenumbers (normalized to x0)
                	// FE: declaration should be further up
                	double xm=l0/wave_ref[id[iw]-1]-1.;
                	double xi=l0/wave_ref[id[iw]]-1.;
                	double xp=l0/wave_ref[id[iw]+1]-1.;
                	double xT=(xp-xm)/2;
                	double Fm = 1.;
                	double Fi = 1.;
                	double Fp = 1.;
                	double norm= 1.;

                	if (specflat_table)
                	{
                		Fm=cpl_array_get (flat, id[iw]-1, NULL);
                		Fi=cpl_array_get (flat, id[iw], NULL);
                		Fp=cpl_array_get (flat, id[iw]+1, NULL);
                	}

                	// Compute the weigths (initialized at zero)
                	weight_m[iw]=0;
                	weight_i[iw]=0;
                	weight_p[iw]=0;

                	// only compute the weights if the three flats are not null (bad pixel)
                	//  and if xo is below (xm+xi)/2
                	//  and if xo is above (xi+xp)/2
                	double fluxLimit = 40 ; // below this value the pixel is consider dead
                	if ( (Fi > fluxLimit) && (Fp > fluxLimit) && (Fm > fluxLimit) && (xm-xi>1e-4) && (xi-xp>1e-4) )
                	{

                		/* The three equations
                		 * 1-> fi = 1 =>  conservation
                		 * 2-> fm*xm + fi*xi + fp*xp = 0  => First moment
                		 * 3-> fm*xm^3 + fi*xi^3 + fp*xp^3 = 0 => Third moment
                         are now approximated by:
                		 * 1-> fi = 1 =>  conservation
                		 * 2-> fm*xm + fi*xi + fp*xp = 0  => First moment
                		 * 3-> fm*xm*(-xT-xi)^2 + fi*xi^3 + fp*xp*(xT-xi)^2 = 0 => Third moment
                		 * so the solution is simpler, and do not diverge: */

                		weight_p[iw]= (xT-2*xi)/(4*xp) ;
                		weight_i[iw]= 1. ;
                		weight_m[iw]= -(xT+2*xi)/(4*xm);

                		/* Acounting for flat
                		 * W/=Flat
                		 */
                		weight_p[iw]/=Fp;
                		weight_i[iw]/=Fi;
                		weight_m[iw]/=Fm;
                	}
                	else if ((Fi > fluxLimit) && (Fp > fluxLimit) && (xi>1e-4) && (xi-xp>1e-4) )
                	{ // left pixel is dead but the interpolation is between 2 good pixels
                		weight_p[iw]= -xi/xp ;
                		weight_i[iw]= 1. ;
                		weight_m[iw]= 0;

                		/* Acounting for flat
                		 * W/=Flat
                		 */
                		weight_p[iw]/=Fp;
                		weight_i[iw]/=Fi;
                		//cpl_msg_info(cpl_func,"Reg %d pix %d 2 right : i %g %g %g x %g %g %g", reg, iw, Fp, Fi, Fm, xp, xi, xm);
                	}
                    else if ( (Fm > fluxLimit) && (Fi > fluxLimit) && (xi<1e-4) && (xm-xi>1e-4))
                    { // right pixel is dead but the interpolation is between 2 good pixels

                		weight_p[iw]= 0;
                		weight_i[iw]= 1. ;
                		weight_m[iw]= -xi/xm ;

                		/* Acounting for flat
                		 * W/=Flat
                		 */
                		weight_m[iw]/=Fm;
                		weight_i[iw]/=Fi;
                		//cpl_msg_info(cpl_func,"Reg %d pix %d 2 left : i %g %g %g x %g %g %g", reg, iw, Fp, Fi, Fm, xp, xi, xm);
                    }
                    //else cpl_msg_info(cpl_func,"Reg %d pix %d no iterp : i %g %g %g x %g %g %g", reg, iw, Fp, Fi, Fm, xp, xi, xm);

                	// normalize the weight to conserve the variance
                	//norm = sqrt( pow(weight_m[iw],2) +  pow(weight_i[iw],2) + pow(weight_p[iw],2));
                	/* FE: the original 2 pixel interpolation applies normalization to preserve photon noise) */
                	norm = weight_m[iw] +  weight_i[iw] + weight_p[iw];
                	/* SL: changed to simple sum(weights) =1 */
                	if ( norm > 1e-7 )
                	{
                    	weight_p[iw]/=norm;
                    	weight_i[iw]/=norm;
                    	weight_m[iw]/=norm;
                	}

                	/* FE: for debugging puposes*/
                	// cpl_msg_info(cpl_func," iw: %i l0: %g id[iw]: %i wm: %g wi: %g wp: %g Fm: %g Fi: %g Fp: %g xm: %g xi: %g xp: %g norm: %g\n",
                	//	 iw, l0, id[iw], weight_m[iw], weight_i[iw], weight_p[iw], Fm, Fi, Fp, xm, xi, xp,norm);


                }// End loop on nb_oiwave

                /* Loop on frames */
                for (cpl_size j = 0; j < nb_row; j++){
                    double * yref;

                    /*
                     *  Compute the fluxes
                     */
                    yref = cpl_array_get_data_double (inputData[j]);
                    /* loop on wavelength */
                    for (cpl_size iw=0 ; iw < nb_oiwave ; iw ++) {
                        outputData[iw]=yref[id[iw]-1]*weight_m[iw] + yref[id[iw]]*weight_i[iw] + yref[id[iw]+1]*weight_p[iw];
                    }
                    /* Copy output array into input array */
                    for (cpl_size iw = 0 ; iw < nb_oiwave ; iw ++){
                        yref[iw]=outputData[iw];
                    }

                    /*
                     *  Compute the noise (sqrt of variance)
                     */
                    yref = cpl_array_get_data_double (inputErr[j]);
                    /* loop on wavelength */
                    for (cpl_size iw=0 ; iw < nb_oiwave ; iw ++) {
                        outputErr[iw]=sqrt(pow(yref[id[iw]-1]*weight_m[iw],2) + pow(yref[id[iw]]*weight_i[iw],2) + pow(yref[id[iw]+1]*weight_p[iw],2));
                    }
                    /* Copy output array into input array */
                    for (cpl_size iw = 0 ; iw < nb_oiwave ; iw ++){
                        yref[iw]=outputErr[iw];
                    }
                } /* end loop on frames */

                /* Free the weight for interpolation */
                FREE (cpl_free, weight_m);
                FREE (cpl_free, weight_i);
                FREE (cpl_free, weight_p);


            } /* End of interpolation on 3 pixels */
            else
            {
                /* starting interpolation 2 pixels */
                cpl_msg_info_overwritable (cpl_func, "Reinterpolate region %lld "
                                      "over %lld (%lld->%lld channels)",
                                      reg+1,nb_region,nb_wave,nb_oiwave);
                /* starting interpolation 2 pixels */
                /* Init the index and weights for the interpolation such that
                 * yout[iw] =   inputData[id[iw]-1] * weight_m[iw]
                 *            + inputData[id[iw]] * weight_p[iw]
                 * yerr[iw] =  sqrt( (inputErr[id[iw]-1] * weight_m[iw])^2
                 *            + (inputErr[id[iw]] * weight_p[iw])^2 )
                 * We here assume monotonicity and no extrapolation */
                double * weight_m  = cpl_malloc (nb_oiwave * sizeof(double));
                double * weight_p  = cpl_malloc (nb_oiwave * sizeof(double));
                
                for (cpl_size iw = 0 ; iw < nb_oiwave ; iw ++)
                {
                    // xo is the wavenumber of the targeted wavelength (of indice iw)
                    double l0=cpl_table_get (oiwave_tables[pol], "EFF_WAVE", iw, NULL);
                    
                    // id[iw] is the index of the closest wavenumber to xo
                    // id[iw] must be above 1, and below nb_wave-1
                    // wave is assumed to be in increasing order
                    
                    if (iw == 0) id[iw] = 1;
                    else id[iw] = id[iw-1];
                    while ( (1./l0 < 1./wave_ref[id[iw]]) && ( id[iw] < nb_wave-1 ) ){
                        id[iw]++;
                    }
                    
                    // compute the wavenumbers (normalized to x0)
                    // FE: declaration should be further up
                    double xm=l0/wave_ref[id[iw]-1]-1.;
                    double xp=l0/wave_ref[id[iw]]-1.;
                    double Fm = 1.;
                    double Fp = 1.;
                    double norm = 1.;
                    
                    
                    if (specflat_table)
                    {
                        Fm=cpl_array_get (flat, id[iw]-1, NULL);
                        Fp=cpl_array_get (flat, id[iw], NULL);
                    }
                    
                    // Compute the weigths (initialized at zero)
                    weight_m[iw]=0;
                    weight_p[iw]=0;
                    
                    // only compute the weights if the three flats are not null (bad pixel)
                    //  and if xo is below (xm+xi)/2
                    //  and if xo is above (xi+xp)/2
                    if ( (Fp>1e-4) && (Fm>1e-4) && (xm>=0.0) && (xp<=0.0) )
                    {
                        
                        /* The three equations
                         * 1-> fm + fp = xm - xp =>  conservation
                         * 2-> fm*xm + fp*xp = 0  => First moment */
                        
                        weight_p[iw]= xm ;
                        weight_m[iw]= -xp ;
                        
                        /* Acounting for flat
                         * W/=Flat
                         */
                        weight_p[iw]/=Fp;
                        weight_m[iw]/=Fm;
                        
                        // normalize the weight to conserve the variance
                        //norm = sqrt( pow(weight_m[iw],2) + pow(weight_p[iw],2));
                        norm = weight_m[iw]+ weight_p[iw];
                        /* FE: the original 2 pixel interpolation applies normalization to preserve photon noise) */
                        /* SL: changed to simple sum(weights) =1 */
                        weight_p[iw]/=norm;
                        weight_m[iw]/=norm;
                        
                    }
                    
                    /* FE: for debugging puposes
                    cpl_msg_info(cpl_func," reg: %i iw: %i l0: %g id[iw]: %i wm: %g wp: %g Fm: %g Fp: %g xm: %g  xp: %g \n",
                     iw, iw, l0, id[iw], weight_m[iw], weight_p[iw], Fm,  Fp, xm, xp);*/
                    
                    
                    
                }// End loop on nb_oiwave
                
                /* Loop on frames */
                for (cpl_size j = 0; j < nb_row; j++){
                    double * yref;
                    
                    /*
                     *  Compute the fluxes
                     */
                    yref = cpl_array_get_data_double (inputData[j]);
                    /* loop on wavelength */
                    for (cpl_size iw=0 ; iw < nb_oiwave ; iw ++) {
                        outputData[iw]=yref[id[iw]-1]*weight_m[iw] + yref[id[iw]]*weight_p[iw];
                    }
                    /* Copy output array into input array */
                    for (cpl_size iw = 0 ; iw < nb_oiwave ; iw ++){
                        yref[iw]=outputData[iw];
                    }
                    
                    /*
                     *  Compute the noise (sqrt of variance)
                     */
                    yref = cpl_array_get_data_double (inputErr[j]);
                    /* loop on wavelength */
                    for (cpl_size iw=0 ; iw < nb_oiwave ; iw ++) {
                        outputErr[iw]=sqrt(pow(yref[id[iw]-1]*weight_m[iw],2) + pow(yref[id[iw]]*weight_p[iw],2));
                    }
                    /* Copy output array into input array */
                    for (cpl_size iw = 0 ; iw < nb_oiwave ; iw ++){
                        yref[iw]=outputErr[iw];
                    }
                } /* end loop on frames */
                
                /* Free the weight for interpolation */
                FREE (cpl_free, weight_m);
                FREE (cpl_free, weight_p);
                
               /* end of interpolation on 2 pixels */
            }

            /* Modify the depth of the region (remove useless pixels) */
            if (nb_wave > nb_oiwave) {
                cpl_msg_debug (cpl_func,"Modify depth of column %s (%lld->%lld)", data_x, nb_wave, nb_oiwave);
                cpl_table_set_column_depth (spectrum_table, data_x, nb_oiwave);
                cpl_table_set_column_depth (spectrum_table, data_errx, nb_oiwave);
                CPLCHECK_MSG ("Cannot change column depth");
            }
            
        } /* End loop on regions i */
        
    } /* End loop on polarisation */
        
    /* Desallocate the modified wavelength array */
    FREE (cpl_free, yout);
    FREE (cpl_free, xout);

    /* Free the weight for interpolation */
    FREE (cpl_free, id);
    FREE (cpl_free, outputData);
    FREE (cpl_free, outputErr);


    gravi_msg_function_exit(1);
    return CPL_ERROR_NONE;
}


/*----------------------------------------------------------------------------*/
/**
 * @brief Regrid the regions into a common wavelength (in-place)
 * 
 * @param spectrum_data     The SPECTRUM data to regrid
 * @param wave_map          The WAVE calibration map (current grid)
 * @param p2vm_map or NULL  The P2VM calibration map (target grid)
 * @param det_type          The detector to align. If GRAVI_DET_SC, then only
 *                          the science detector extensions will be processed.
 *                          GRAVI_DET_FT will do the same for FT detector
 *                          and GRAVI_DET_ALL will do it for both. 
 * @param parlist           Input parameter list (no parameter used).
 * 
 * It re-samples
 * the spectral element according to the wavelength calibration.
 * The SPECTRUM are re-interpolated into the OI_WAVELENGTH of this
 * P2VM map. The target wavelength are the one of the OI_WAVELENGTH
 * from the P2VM map (assumed to be the same for both polarisation!!)
 *
 * This is done in-place to save enormous time, especially on FT.
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_align_spectrum (gravi_data * spectrum_data,
                                     gravi_data * wave_map,
                                     gravi_data * p2vm_map,
                                     enum gravi_detector_type det_type)
{
    gravi_msg_function_start(1);


	/* Check needed data */
	cpl_ensure_code (spectrum_data, CPL_ERROR_NULL_INPUT);
	cpl_ensure_code (wave_map,      CPL_ERROR_NULL_INPUT);
	cpl_ensure_code (p2vm_map,      CPL_ERROR_NULL_INPUT);

	/* Duplicate extensions in the spectrum data from p2vm_map 
    *  FIXME: shall copy only required from the type FT/SC */
	gravi_data_copy_ext (spectrum_data, p2vm_map, GRAVI_OI_WAVELENGTH_EXT);
    CPLCHECK_MSG ("Cannot copy OI_WAVELENGTH extention(s)");
	
	/* Loop on FT/SC */
    int init_type_data = 1;
    int end_type_data = 1;
    if(det_type == GRAVI_DET_SC || det_type == GRAVI_DET_ALL) {
        init_type_data = 0;
    }
    if(det_type == GRAVI_DET_FT || det_type == GRAVI_DET_ALL) {
        end_type_data = 1;
    }
	for (int type_data = init_type_data; type_data <= end_type_data; type_data++ ) {

		/* Check if SPECTRUM data exists */
        if (!gravi_data_has_spectrum (spectrum_data, type_data)) {
            cpl_msg_info (cpl_func,"No data for %s, skip it", GRAVI_TYPE(type_data));
            continue;
		}
		
		/* Get a timer */
	    cpl_msg_info (cpl_func, "Re-interpolate the %s",GRAVI_TYPE(type_data));

		/* Get all the input tables */
        cpl_table * spectrum_table, * wave_table, * detector_table;
        spectrum_table = gravi_data_get_spectrum_data (spectrum_data, type_data);
        wave_table = gravi_data_get_wave_data (wave_map, type_data);
        detector_table = gravi_data_get_imaging_detector (spectrum_data, type_data);

        cpl_table * specflat_table;
		if (type_data ==  GRAVI_SC) 
			specflat_table = gravi_data_get_table (spectrum_data, GRAVI_SPECTRUMFLAT_DATA_SC_EXT);
		else
            specflat_table = NULL;
        
		CPLCHECK_MSG ("Cannot load the extention data");

        /* Measure total flux before interpolation */
        double full_flux_reg0 = gravi_spectrum_get_flux (spectrum_table);

        /* Get the OI_WAVE table of each polarisation */
        int npol = gravi_spectrum_get_npol (wave_table);
        cpl_table ** oiwave_tables = gravi_data_get_oiwave_tables (p2vm_map, type_data, npol);
        
        /* Re-interpolate the spectrum into OI_WAVELENGTH
         * Run interpolation -- in-place to save enormous time */
        gravi_interpolate_spectrum_table (spectrum_table,
                                          wave_table,
                                          oiwave_tables,
                                          detector_table,
                                          specflat_table,
                                          type_data);
        
        FREE (cpl_free, oiwave_tables);
        CPLCHECK_MSG ("Cannot interpolate");

        /* Measure total flux after interpolation */
        double full_flux_reg1 = gravi_spectrum_get_flux (spectrum_table);

        /* Replace (FULLSTARTX,STARTX,NX) by NWAVE in the table plist */
        cpl_propertylist * spectrum_plist = gravi_data_get_spectrum_data_plist (spectrum_data, type_data);
        int nwave = gravi_spectrum_get_nwave (spectrum_table);
        cpl_propertylist_erase (spectrum_plist, PROFILE_FULLSTARTX);
        cpl_propertylist_erase (spectrum_plist, PROFILE_STARTX);
        cpl_propertylist_erase (spectrum_plist, PROFILE_NX);
        cpl_propertylist_update_int (spectrum_plist, "NWAVE", nwave);
        CPLCHECK_MSG ("Cannot replace STARTX by NWAVE");
            
        /* Check if the interpolation is flux conservative */
        cpl_msg_info (cpl_func, "Total flux in REGION1s: %.2e [e], in REGION0s:%.2e [e]  (ratio=%.5f)",
                      full_flux_reg1,full_flux_reg0,full_flux_reg1/full_flux_reg0);
        
        char qc_name[80];
        cpl_propertylist * spectrum_header = gravi_data_get_header (spectrum_data);
        sprintf (qc_name, "ESO QC TRANS INTERP %s",GRAVI_TYPE(type_data));
        cpl_propertylist_update_double (spectrum_header, qc_name, full_flux_reg1/full_flux_reg0);
        cpl_propertylist_set_comment (spectrum_header, qc_name, "[e/e] at interpolation");
        CPLCHECK_MSG ("Cannot set QC parameters");
        
    } /* End loop on SC / FT */
	
    gravi_msg_function_exit(1);
	return CPL_ERROR_NONE;
}

/**@}*/
