/* $Id: gravi_cpl.c,v 1.10 2014/11/12 15:10:40 nazouaoui Exp $
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
 * @defgroup gravi_cpl  Gravity CPL interface
 *
 * This module implements some interactions with cpl objects to address specific
 * gravity needs. The manipulated cpl object are @c cpl_matrix, @c cpl_array,
 * @c cpl_table and @c cpl_image
 *
 */
/**@{*/

/*-----------------------------------------------------------------------------
                                   Includes
 -----------------------------------------------------------------------------*/

#define _XOPEN_SOURCE // for M_PI symbol

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

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

#include "gravi_calib.h"
#include "gravi_vis.h"
#include "gravi_tf.h"
#include "gravi_eop.h"

/*-----------------------------------------------------------------------------
                              Private prototypes
 -----------------------------------------------------------------------------*/

double pythag(double , double );

cpl_matrix * svdcmp(cpl_matrix * , cpl_vector * , cpl_matrix * );

/*-----------------------------------------------------------------------------
                             Functions code
 -----------------------------------------------------------------------------*/


int gravi_array_threshold_min (cpl_array * array, double lo_cut,
                               double assign_lo_cut)
{
    cpl_ensure (array, CPL_ERROR_NULL_INPUT, -1);

    cpl_size size = cpl_array_get_size (array);
    cpl_size num = 0;
    
    for (cpl_size s = 0; s < size ; s++) {
        if (cpl_array_get (array, s, NULL) < lo_cut) {
            cpl_array_set (array, s, assign_lo_cut);
            num ++;
        }
    }
    
    return num;
}


/**
 * @brief Allocate a list of arrays, pre-filled with 0.0
 */

cpl_array ** gravi_array_new_list (int n, cpl_type type, int size)
{
  cpl_ensure (n>0,    CPL_ERROR_ILLEGAL_INPUT, NULL);
  cpl_ensure (size>0, CPL_ERROR_ILLEGAL_INPUT, NULL);
  
  cpl_array ** output;
  output = cpl_malloc (n * sizeof(cpl_array *));

  int out1;
  if ( type == CPL_TYPE_DOUBLE )
	for ( out1 = 0 ; out1 < n ; out1++ )  {
	  output[out1] = cpl_array_new (size, type);
	  cpl_array_fill_window_double (output[out1], 0, size, 0.0);
	}
  else if ( type == CPL_TYPE_DOUBLE_COMPLEX )
	for ( out1 = 0 ; out1 < n ; out1++ )  {
	  output[out1] = cpl_array_new (size, type);
	  cpl_array_fill_window_double_complex (output[out1], 0, size, 0.0 * I*0.0);
	}
  else {
	cpl_error_set_message (cpl_func,CPL_ERROR_ILLEGAL_INPUT,"This type is not supported.");
	FREE (cpl_free, output);
	return NULL;
  }
  
  return output;
}

cpl_error_code gravi_table_interpolate_column (cpl_table * to_table,
                                               const char * to_x,
                                               const char * to_y,
                                               const cpl_table * from_table,
                                               const char * from_x,
                                               const char * from_y)
{
    gravi_msg_function_start(1);
    cpl_ensure_code (to_table,   CPL_ERROR_NULL_INPUT);
    cpl_ensure_code (to_x,       CPL_ERROR_NULL_INPUT);
    cpl_ensure_code (to_y,       CPL_ERROR_NULL_INPUT);
    cpl_ensure_code (from_table, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code (from_x,     CPL_ERROR_NULL_INPUT);
    cpl_ensure_code (from_y,     CPL_ERROR_NULL_INPUT);

    /* Size */
    cpl_size nxref = cpl_table_get_nrow (from_table);
    cpl_size nxout = cpl_table_get_nrow (to_table);

    /* Xref vectors */
    cpl_vector * xref = cpl_vector_new (nxref);
    for (cpl_size row = 0; row < nxref; row++) {
        cpl_vector_set (xref, row, cpl_table_get (from_table, from_x, row, NULL));
        CPLCHECK_MSG ("Cannot get x data");
    }

    /* Xout vectors */
    cpl_vector * xout = cpl_vector_new (nxout);
    for (cpl_size row = 0; row < nxout; row++) {
        cpl_vector_set (xout, row, cpl_table_get (to_table, to_x, row, NULL));
        CPLCHECK_MSG ("Cannot get x data");
    }

    /* Allocate memory for the interpolation */
    cpl_vector * yref = cpl_vector_new (nxref);
    cpl_vector * yout = cpl_vector_new (nxout);
    cpl_bivector * fref = cpl_bivector_wrap_vectors (xref, yref);
    cpl_bivector * fout = cpl_bivector_wrap_vectors (xout, yout);

    /* Shall be a column of array */
    cpl_size depth = cpl_table_get_column_depth (from_table, from_y);
    if (depth == 0) 
    	cpl_error_set_message (cpl_func,CPL_ERROR_INVALID_TYPE ,
                               "Report this error to gravity.drs");

    /* Output is of type DOUBLE */
    const char * unit = cpl_table_get_column_unit (from_table, from_y);
    gravi_table_init_column_array (to_table, to_y, unit, CPL_TYPE_DOUBLE, depth);

    /* Loop on column depth */
    for (cpl_size size = 0; size < depth; size++) {
        
        /* Yref vector */
        for (cpl_size row = 0; row < nxref; row++) {
            double value = gravi_table_get_value (from_table, from_y, row, size);
            cpl_vector_set (yref, row, value);
            CPLCHECK_MSG ("Cannot get y data");
        }

        /* Interpolate linear */
        cpl_bivector_interpolate_linear (fout, fref);
        CPLCHECK_MSG ("Cannot interpolate");

        /* fill Yout */
        for (cpl_size row = 0; row < nxout; row++) {
            double value = cpl_vector_get (yout, row);
            gravi_table_set_value (to_table, to_y, row, size, value);
            CPLCHECK_MSG ("Cannot set y data");
        }
        
    } /* End loop on column depth */

    /* Free memory */
    FREE (cpl_bivector_delete, fref);
    FREE (cpl_bivector_delete, fout);

    gravi_msg_function_exit(1);
    return CPL_ERROR_NONE;
}

double gravi_table_get_column_mean (cpl_table * table, const char * name, int base, int nbase)
{
  cpl_ensure (table, CPL_ERROR_NULL_INPUT, 0.0);
  cpl_ensure (name,  CPL_ERROR_NULL_INPUT, 0.0);
  
  double mean = 0.0;
  cpl_size nrow = cpl_table_get_nrow (table) / nbase;
  cpl_ensure (nrow,  CPL_ERROR_ILLEGAL_INPUT, 0.0);
  
  cpl_type type  = cpl_table_get_column_type (table, name);
  cpl_size depth = cpl_table_get_column_depth (table, name);
  
  if (depth == 0 && type == CPL_TYPE_DOUBLE) {
	double * data = cpl_table_get_data_double (table, name);
	cpl_ensure (data, CPL_ERROR_ILLEGAL_INPUT, 0.0);
	for (cpl_size r=0; r<nrow;r++) mean += data[r*nbase+base];
  }
  else if (depth == 0 && type == CPL_TYPE_FLOAT) {
	float * data = cpl_table_get_data_float (table, name);
	cpl_ensure (data, CPL_ERROR_ILLEGAL_INPUT, 0.0);
	for (cpl_size r=0; r<nrow;r++) mean += data[r*nbase+base];
  }
  else if (depth == 0 && type == CPL_TYPE_INT) {
	int * data = cpl_table_get_data_int (table, name);
	cpl_ensure (data, CPL_ERROR_ILLEGAL_INPUT, 0.0);
	for (cpl_size r=0; r<nrow;r++) mean += data[r*nbase+base];
  }
  else if (depth > 0) {
      cpl_array ** arrays = cpl_table_get_data_array (table, name);
      cpl_ensure (arrays,  CPL_ERROR_ILLEGAL_INPUT, 0.0);
      cpl_array * output = cpl_array_duplicate (arrays[base]);
      cpl_ensure (output,  CPL_ERROR_ILLEGAL_INPUT, 0.0);
      for (cpl_size r=1; r<nrow;r++)
          cpl_array_add (output, arrays[r*nbase+base]);
      mean = cpl_array_get_mean (output) / nrow;
      FREE (cpl_array_delete, output);
  }
  else {
	cpl_error_set_message (cpl_func,CPL_ERROR_ILLEGAL_INPUT,"unknow type");
	return 0.0;
  }

  return mean / nrow;
}

double gravi_table_get_column_std (cpl_table * table, const char * name, int base, int nbase)
{
  cpl_ensure (table, CPL_ERROR_NULL_INPUT, 0.0);
  cpl_ensure (name,  CPL_ERROR_NULL_INPUT, 0.0);
  
  double mean  = 0.0;
  double mean2 = 0.0;
  cpl_size nrow = cpl_table_get_nrow (table) / nbase;
  cpl_ensure (nrow,  CPL_ERROR_ILLEGAL_INPUT, 0.0);
  
  cpl_type type  = cpl_table_get_column_type (table, name);
  cpl_size depth = cpl_table_get_column_depth (table, name);
  
  if (depth == 0 && type == CPL_TYPE_DOUBLE) {
	double * data = cpl_table_get_data_double (table, name);
	cpl_ensure (data, CPL_ERROR_ILLEGAL_INPUT, 0.0);
	for (cpl_size r=0; r<nrow;r++) mean  += data[r*nbase+base];
	for (cpl_size r=0; r<nrow;r++) mean2 += data[r*nbase+base] * data[r*nbase+base];
  }
  else if (depth == 0 && type == CPL_TYPE_FLOAT) {
	float * data = cpl_table_get_data_float (table, name);
	cpl_ensure (data, CPL_ERROR_ILLEGAL_INPUT, 0.0);
	for (cpl_size r=0; r<nrow;r++) mean += data[r*nbase+base];
	for (cpl_size r=0; r<nrow;r++) mean2 += data[r*nbase+base] * data[r*nbase+base];
  }
  else if (depth == 0 && type == CPL_TYPE_INT) {
	int * data = cpl_table_get_data_int (table, name);
	cpl_ensure (data, CPL_ERROR_ILLEGAL_INPUT, 0.0);
	for (cpl_size r=0; r<nrow;r++) mean += data[r*nbase+base];
	for (cpl_size r=0; r<nrow;r++) mean2 += data[r*nbase+base] * data[r*nbase+base];
  }
  else {
	cpl_error_set_message (cpl_func,CPL_ERROR_ILLEGAL_INPUT,"unknow type");
	return 0.0;
  }

  return sqrt (mean2 / nrow - mean*mean / nrow / nrow);
}


cpl_array * gravi_table_get_column_mean_array (cpl_table * table, const char * name, int base, int nbase)
{
  cpl_ensure (table, CPL_ERROR_NULL_INPUT, NULL);
  cpl_ensure (name,  CPL_ERROR_NULL_INPUT, NULL);
  
  cpl_size nrow = cpl_table_get_nrow (table) / nbase;
  cpl_ensure (nrow,  CPL_ERROR_ILLEGAL_INPUT, NULL);

  /* Get the pointer */
  cpl_array ** arrays = cpl_table_get_data_array (table, name);
  cpl_ensure (arrays,  CPL_ERROR_ILLEGAL_INPUT, NULL);

  /* Build the mean */
  cpl_array * output = cpl_array_duplicate (arrays[base]);
  for (cpl_size r=1; r<nrow;r++)
      cpl_array_add (output, arrays[r*nbase+base]);

  cpl_array_divide_scalar (output, nrow);
  return output;
}

  
double ** gravi_table_get_data_array_double(cpl_table * table, const char * name)
{
  cpl_ensure (table, CPL_ERROR_NULL_INPUT, NULL);
  cpl_ensure (name,  CPL_ERROR_NULL_INPUT, NULL);
  
  cpl_size nrow = cpl_table_get_nrow (table);
  cpl_array ** pdata = cpl_table_get_data_array (table, name);

  cpl_ensure (nrow,  CPL_ERROR_ILLEGAL_INPUT, NULL);
  cpl_ensure (pdata, CPL_ERROR_ILLEGAL_INPUT, NULL);

  /* Allocate memory, this will have to be desalocated */
  double ** data = cpl_malloc (sizeof(double*) * nrow);
  
  /* Get all the pointers */
  for (cpl_size row=0; row<nrow; row++) {
	data[row] = cpl_array_get_data_double (pdata[row]);
  }
  
  CPLCHECK_NUL("Cannot load the requested arrays");
  
  return data;
}

float ** gravi_table_get_data_array_float(cpl_table * table, const char * name)
{
  cpl_ensure (table, CPL_ERROR_NULL_INPUT, NULL);
  cpl_ensure (name,  CPL_ERROR_NULL_INPUT, NULL);
  
  cpl_size nrow = cpl_table_get_nrow (table);
  cpl_array ** pdata = cpl_table_get_data_array (table, name);

  cpl_ensure (nrow,  CPL_ERROR_ILLEGAL_INPUT, NULL);
  cpl_ensure (pdata, CPL_ERROR_ILLEGAL_INPUT, NULL);

  /* Allocate memory, this will have to be desalocated */
  float ** data = cpl_malloc (sizeof(float*) * nrow);
  
  /* Get all the pointers */
  for (cpl_size row=0; row<nrow; row++) {
	data[row] = cpl_array_get_data_float (pdata[row]);
  }
  
  CPLCHECK_NUL("Cannot load the requested arrays");
  
  return data;
}

float complex ** gravi_table_get_data_array_float_complex (cpl_table * table, const char * name)
{
  cpl_ensure (table, CPL_ERROR_NULL_INPUT, NULL);
  cpl_ensure (name,  CPL_ERROR_NULL_INPUT, NULL);
  
  cpl_size nrow = cpl_table_get_nrow (table);
  cpl_array ** pdata = cpl_table_get_data_array (table, name);

  cpl_ensure (nrow,  CPL_ERROR_ILLEGAL_INPUT, NULL);
  cpl_ensure (pdata, CPL_ERROR_ILLEGAL_INPUT, NULL);

  /* Allocate memory, this will have to be desalocated */
  float complex ** data = cpl_malloc (sizeof(float complex*) * nrow);
  
  /* Get all the pointers */
  for (cpl_size row=0; row<nrow; row++) {
	data[row] = cpl_array_get_data_float_complex (pdata[row]);
  }
  
  CPLCHECK_NUL ("Cannot load the requested arrays");
  
  return data;
}


double complex ** gravi_table_get_data_array_double_complex (cpl_table * table, const char * name)
{
  cpl_ensure (table, CPL_ERROR_NULL_INPUT, NULL);
  cpl_ensure (name,  CPL_ERROR_NULL_INPUT, NULL);
  
  cpl_size nrow = cpl_table_get_nrow (table);
  cpl_array ** pdata = cpl_table_get_data_array (table, name);

  cpl_ensure (nrow,  CPL_ERROR_ILLEGAL_INPUT, NULL);
  cpl_ensure (pdata, CPL_ERROR_ILLEGAL_INPUT, NULL);

  /* Allocate memory, this will have to be desalocated */
  double complex ** data = cpl_malloc (sizeof(double complex*) * nrow);
  
  /* Get all the pointers */
  for (cpl_size row=0; row<nrow; row++) {
	data[row] = cpl_array_get_data_double_complex (pdata[row]);
  }
  
  CPLCHECK_NUL ("Cannot load the requested arrays");
  
  return data;
}

int ** gravi_table_get_data_array_int(cpl_table * table, const char * name)
{
  cpl_ensure (table, CPL_ERROR_NULL_INPUT, NULL);
  cpl_ensure (name,  CPL_ERROR_NULL_INPUT, NULL);
  
  cpl_size nrow = cpl_table_get_nrow (table);
  cpl_array ** pdata = cpl_table_get_data_array (table, name);

  cpl_ensure (nrow,  CPL_ERROR_ILLEGAL_INPUT, NULL);
  cpl_ensure (pdata, CPL_ERROR_ILLEGAL_INPUT, NULL);

  /* Allocate memory, this will have to be desalocated */
  int ** data = cpl_malloc (sizeof(int*) * nrow);
  
  /* Get all the pointers */
  for (cpl_size row=0; row<nrow; row++) {
	data[row] = cpl_array_get_data_int (pdata[row]);
  }
  
  CPLCHECK_NUL("Cannot load the requested arrays");
  
  return data;
}

/* 
 * Define and init an array of DOUBLE 
 */
cpl_array * gravi_array_init_double (long n , double value)
{
  cpl_ensure (n>0, CPL_ERROR_ILLEGAL_INPUT, NULL);
  cpl_array * output = cpl_array_new (n, CPL_TYPE_DOUBLE);
  cpl_array_fill_window_double (output, 0, n, value);
  return output;
}

/* 
 * Define and init an array of INT
 */
cpl_array * gravi_array_init_int (long n, int value)
{
  cpl_ensure (n>0, CPL_ERROR_ILLEGAL_INPUT, NULL);
  cpl_array * output = cpl_array_new (n, CPL_TYPE_INT);
  cpl_array_fill_window_int (output, 0, n, value);
  return output;
}

/* 
 * Define and init an array of DOUBLE COMPLEX 
 */
cpl_array * gravi_array_init_double_complex (long n, double complex value)
{
  cpl_ensure (n>0, CPL_ERROR_ILLEGAL_INPUT, NULL);
  cpl_array * output = cpl_array_new (n, CPL_TYPE_DOUBLE_COMPLEX);
  cpl_array_fill_window_double_complex (output, 0, n, value);
  return output;
}			

/* 
 * Define and init an array of FLOAT COMPLEX 
 */
cpl_array * gravi_array_init_float_complex (long n, float complex value)
{
  cpl_ensure (n>0, CPL_ERROR_ILLEGAL_INPUT, NULL);
  cpl_array * output = cpl_array_new (n, CPL_TYPE_FLOAT_COMPLEX);
  cpl_array_fill_window_float_complex (output, 0, n, value);
  return output;
}			

/* 
 * Compute a new DOUBLE COMPLEX array filled with:  input_re + i * input_im 
 */
cpl_array * gravi_array_wrap_complex (cpl_array * input_re, cpl_array * input_im)
{
  cpl_ensure (input_re, CPL_ERROR_NULL_INPUT, NULL);
  cpl_ensure (input_im, CPL_ERROR_NULL_INPUT, NULL);
  
  cpl_size size_re = cpl_array_get_size (input_re);
  cpl_size size_im = cpl_array_get_size (input_im);

  cpl_ensure (size_re == size_im, CPL_ERROR_ILLEGAL_INPUT, NULL);
  
  cpl_array * output = cpl_array_new (size_re, CPL_TYPE_DOUBLE_COMPLEX);

  for (cpl_size  n = 0; n < size_re; n ++) {
	cpl_array_set_complex (output, n, 1.* I * cpl_array_get (input_im, n, NULL) +
                           cpl_array_get (input_re, n, NULL));
  }

  return output;
}

/* 
 * Compute a new FLOAT COMPLEX array filled with:  input_re + i * input_im 
 */
cpl_array * gravi_array_wrap_float_complex (cpl_array * input_re, cpl_array * input_im)
{
  cpl_ensure (input_re, CPL_ERROR_NULL_INPUT, NULL);
  cpl_ensure (input_im, CPL_ERROR_NULL_INPUT, NULL);
  
  cpl_size size_re = cpl_array_get_size (input_re);
  cpl_size size_im = cpl_array_get_size (input_im);

  cpl_ensure (size_re == size_im, CPL_ERROR_ILLEGAL_INPUT, NULL);
  
  cpl_array * output = cpl_array_new (size_re, CPL_TYPE_FLOAT_COMPLEX);
  
  int nv = 0.0;
  for (cpl_size n = 0; n < size_re; n ++) {
  	cpl_array_set_float_complex (output, n, 1.* I * cpl_array_get (input_im, n, &nv) +
  								  cpl_array_get (input_re, n, &nv));
  }

  return output;
}

/* 
 * Compute a new DOUBLE array filled with:  input_re^2 + input_im^2 
 */
cpl_array * gravi_array_compute_norm2 (cpl_array * input_re, cpl_array * input_im)
{
  cpl_ensure (input_re, CPL_ERROR_NULL_INPUT, NULL);
  cpl_ensure (input_im, CPL_ERROR_NULL_INPUT, NULL);

  cpl_size size_re = cpl_array_get_size (input_re);
  cpl_size size_im = cpl_array_get_size (input_im);

  cpl_ensure (size_re == size_im, CPL_ERROR_ILLEGAL_INPUT, NULL);
  
  cpl_array * output = cpl_array_new (size_re, CPL_TYPE_DOUBLE);
  cpl_array_fill_window_double (output, 0, size_re, 0.0);

  int nv = 0.0;
  for (cpl_size n = 0; n < size_re; n ++) {
	cpl_array_set_double (output, n, pow(cpl_array_get (input_im, n, &nv),2) +
						             pow(cpl_array_get (input_re, n, &nv),2) );
  }

  return output;
}

/* 
 * Set an pair of arrays as the real and imaginary part of
 * a FLOAT COMPLEX array into the table. Data are copied.
 */
cpl_error_code gravi_table_set_array_double_complex (cpl_table * table,
                                                     const char * name,
                                                     cpl_size row,
                                                     cpl_array * visR,
                                                     cpl_array * visI)
{
  cpl_ensure_code (table, CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (name,  CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (visR,  CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (visI,  CPL_ERROR_NULL_INPUT);
  
  cpl_array * tmp_cast = gravi_array_wrap_complex (visR, visI);
  cpl_table_set_array (table, name, row, tmp_cast);
  cpl_array_delete (tmp_cast);

  CPLCHECK_MSG("Cannot set float_complex array");
  
  return CPL_ERROR_NONE;
}

/* 
 * Set an array in radian into an table and convert it in DOUBLE
 * with units [deg]
 */
cpl_error_code gravi_table_set_array_phase (cpl_table * table, const char * name, cpl_size row, cpl_array * phase)
{
  cpl_ensure_code (table, CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (name,  CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (phase, CPL_ERROR_NULL_INPUT);

  cpl_array * tmp_cast;
  if (cpl_array_get_type (phase) == CPL_TYPE_FLOAT_COMPLEX ||
	  cpl_array_get_type (phase) == CPL_TYPE_DOUBLE_COMPLEX ) {
	tmp_cast = cpl_array_cast (phase, CPL_TYPE_DOUBLE_COMPLEX);
	cpl_array_arg (tmp_cast);
  } else 
	tmp_cast = cpl_array_cast (phase, CPL_TYPE_DOUBLE);
  
  cpl_array_multiply_scalar (tmp_cast, 180.0/ CPL_MATH_PI);
  cpl_table_set_array (table, name, row, tmp_cast);
  cpl_array_delete (tmp_cast);

  CPLCHECK_MSG("Cannot set phase array");
  
  return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Set string in table, ensuring fixed length (right space padding).
            see cpl_table_set_string.
  @param	len, the specified maximum length of the string.
  @return   CPL_ERRROR_NONE
 */
/*----------------------------------------------------------------------------*/
cpl_error_code gravi_table_set_string_fixlen (cpl_table *table, const char *name, int row, const char *value, int len)
{
  gravi_msg_function_start(0);
  cpl_ensure_code (table,  CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (name,   CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (row>-1, CPL_ERROR_ILLEGAL_INPUT);
  cpl_ensure_code (len>0,  CPL_ERROR_ILLEGAL_INPUT);

  char * str = cpl_sprintf ("%-*.*s", len, len, value);
  cpl_table_set_string (table, name, row, str);
  
  FREE (cpl_free, str);
  CPLCHECK_MSG("Cannot set string");
  
  gravi_msg_function_exit(0);
  return CPL_ERROR_NONE;
}

/* 
 * Rebin a cpl_array from the SC wavelength table to the FT wavelength table
 * For each FT bin, average all elements found in lbdEff+/-lbdBand/2
 * The output array is created, the input array is not modified.
 */
cpl_array * gravi_array_rebin (const cpl_array * input, const cpl_array * errs,
							   cpl_table * oi_wave_sc, cpl_table * oi_wave_ft)
{
  gravi_msg_function_start(0);
  cpl_ensure (input, CPL_ERROR_NULL_INPUT, NULL);
  cpl_ensure (errs,  CPL_ERROR_NULL_INPUT, NULL);
  cpl_ensure (oi_wave_sc, CPL_ERROR_NULL_INPUT, NULL);
  cpl_ensure (oi_wave_ft, CPL_ERROR_NULL_INPUT, NULL);
  
  /* Get the eff_wave and eff_band of the FT */
  float *effWave_sc = cpl_table_get_data_float (oi_wave_sc, "EFF_WAVE");
  float *effWave_ft = cpl_table_get_data_float (oi_wave_ft, "EFF_WAVE");
  float *effBand_ft = cpl_table_get_data_float (oi_wave_ft, "EFF_BAND");
  cpl_size nwave_ft = cpl_table_get_nrow (oi_wave_ft);
  cpl_size nwave_sc = cpl_table_get_nrow (oi_wave_sc);

  CPLCHECK_NUL ("Cannot get data");
  
  cpl_array * output = cpl_array_new (nwave_ft, CPL_TYPE_DOUBLE);
  cpl_array * weight = cpl_array_new (nwave_ft, CPL_TYPE_DOUBLE);

  /* Built the SC data at the spectral resolution of the FT  (vis2_sc_lr) */
  for (cpl_size wft = 0; wft < nwave_ft; wft++) {
	for (cpl_size wsc = 0; wsc < nwave_sc; wsc++)
	  if ( fabs (effWave_sc[wsc] - effWave_ft[wft]) < 0.5 * effBand_ft[wft] ) {
		
		double v = cpl_array_get (input, wsc, NULL);
		double w = cpl_array_get (errs, wsc, NULL);
		
		w = (w!=0?1./w:1.0);
		cpl_array_set (output, wft, cpl_array_get (output, wft, NULL) + v * w);
		cpl_array_set (weight, wft, cpl_array_get (weight, wft, NULL) + w);
		
		CPLCHECK_NUL("Cannot reduce the spectral resolution");
	  }
  }
  
  cpl_array_divide (output, weight);
  cpl_array_delete (weight);
  gravi_msg_function_exit(0);
  return output;
}


/* 
 * Add two columns from differents tables, without duplicating the data.
 * Use cpl_table_add_columns for columns in the same table.
 */
cpl_error_code gravi_table_add_columns (cpl_table * oi_vis1, const char *name1,
										cpl_table * oi_vis2, const char *name2)
{
  gravi_msg_function_start(0);
  cpl_ensure_code (oi_vis1, CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (oi_vis2, CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (name1,   CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (name2,   CPL_ERROR_NULL_INPUT);

  int row;
  cpl_msg_debug (cpl_func, "Colname (%s + %s) ",name1,name2);
  
  cpl_type type1 = cpl_table_get_column_type (oi_vis1, name1);
  cpl_type type2 = cpl_table_get_column_type (oi_vis2, name2);
  cpl_size nrow1 = cpl_table_get_nrow (oi_vis1);
  cpl_size nrow2 = cpl_table_get_nrow (oi_vis2);
  CPLCHECK_MSG("Cannot get type or nrow");

  if ( type1 != type2 || nrow1 != nrow2) {
	return cpl_error_set_message (cpl_func,CPL_ERROR_ILLEGAL_INPUT,
								  "input columns not conformables");
  }
  
  if ( type1 == CPL_TYPE_DOUBLE ) {
	double * data1 = cpl_table_get_data_double (oi_vis1, name1);
	double * data2 = cpl_table_get_data_double (oi_vis2, name2);
	CPLCHECK_MSG("Cannot load data");

	for (row=0 ; row<nrow1 ; row++) data1[row] += data2[row];
  }
  else if ( type1 == CPL_TYPE_FLOAT_COMPLEX ) {
	float complex * data1 = cpl_table_get_data_float_complex (oi_vis1, name1);
	float complex * data2 = cpl_table_get_data_float_complex (oi_vis2, name2);
	CPLCHECK_MSG("Cannot load data");

	for (row=0 ; row<nrow1; row++) data1[row] += data2[row];
  }
  else if ( type1 == CPL_TYPE_DOUBLE_COMPLEX ) {
	double complex * data1 = cpl_table_get_data_double_complex (oi_vis1, name1);
	double complex * data2 = cpl_table_get_data_double_complex (oi_vis2, name2);
	CPLCHECK_MSG("Cannot load data");

	for (row=0 ; row<nrow1; row++) data1[row] += data2[row];
  }
  else if ( type1 & CPL_TYPE_POINTER ) {
	cpl_array ** data1 = cpl_table_get_data_array (oi_vis1, name1);
	cpl_array ** data2 = cpl_table_get_data_array (oi_vis2, name2);
	CPLCHECK_MSG("Cannot load data");

	for (row=0 ; row<nrow1; row++) {
        cpl_array_add (data1[row], data2[row]);
        CPLCHECK_MSG("Cannot add data");
	}
  } else {
	return cpl_error_set_message (cpl_func,CPL_ERROR_ILLEGAL_INPUT,
                                  "unknow type -- report to DRS team");
  }
  /* End case  */

  gravi_msg_function_exit(0);
  return CPL_ERROR_NONE;
}



/* 
 * Smooth column input_name with a box card of length 2*nsmooth+1
 * Note that this is a running SUM, not a running MEAN.
 * It integrated 2 x nsmooth + 1 running samples
 * FIXME: check this running sum near the limits 
 */
cpl_error_code gravi_table_runint_column (cpl_table * oi_vis,
										  const char *input_name,
										  const char *output_name,
										  int nsmooth, int nbase)
{
  gravi_msg_function_start(0);
  cpl_ensure_code (oi_vis,      CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (input_name,  CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (output_name, CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (nsmooth>=0,   CPL_ERROR_ILLEGAL_INPUT);
  cpl_ensure_code (nbase>0,     CPL_ERROR_ILLEGAL_INPUT);

  int row, base;
  int row_add, row_sub;

  cpl_type type = cpl_table_get_column_type (oi_vis, input_name);
  cpl_size nrow = cpl_table_get_nrow (oi_vis) / nbase;

  if ( type == CPL_TYPE_DOUBLE ) {
	double * data = cpl_table_get_data_double (oi_vis, input_name);
	double * output = cpl_malloc (sizeof(double) * nrow * nbase);
	double buffer = 0.0;

	CPLCHECK_MSG("Cannot load data");

	/* Loop on base and rows */
	for ( base=0 ; base<nbase ; base++) {
	  buffer = 0.0;
	  
	  for ( row=-nsmooth ; row<nrow ; row++) {
		row_add = row+nsmooth;
		row_sub = row-nsmooth-1;
		if ( row_add < nrow ) buffer += data[row_add * nbase + base];
		if ( row_sub >= 0 )   buffer -= data[row_sub * nbase + base];
		if( row>=0 && row<nrow ) output[row * nbase + base] = (double)(buffer);
	  }
	  /* End loop on rows */
	}
	/* End loop on base */
	
	/* Wrap output column */
	if ( !strcmp (input_name, output_name)) {
	  for (row = 0; row < nrow * nbase; row++) data[row] = output[row];
	  cpl_free (output);
	} else {
        if (cpl_table_has_column (oi_vis, output_name))
            cpl_table_erase_column (oi_vis, output_name);
        cpl_table_wrap_double (oi_vis, output, output_name);
	}
	CPLCHECK_MSG("Cannot Wrap");
  }
  else if ( type == CPL_TYPE_DOUBLE_COMPLEX ) {
	double complex * data = cpl_table_get_data_double_complex (oi_vis, input_name);
	double complex * output = cpl_malloc (sizeof(double complex) * nrow * nbase);
	double complex buffer = 0.0;

	CPLCHECK_MSG("Cannot load data");

	/* Loop on rows */
	for ( base=0 ; base<nbase ; base++) {
	  buffer = 0.0;
	  
	  for ( row=-nsmooth ; row<nrow ; row++) {
		row_add = row+nsmooth;
		row_sub = row-nsmooth-1;
		if ( row_add < nrow ) buffer += data[row_add * nbase + base];
		if ( row_sub >= 0 )   buffer -= data[row_sub * nbase + base];
		if( row>=0 && row<nrow ) output[row * nbase + base] = (double complex)(buffer);
	  }
	}
	/* End loop on base */
	
	/* Wrap output column */
	if ( !strcmp (input_name, output_name)) {
	  for (row = 0; row < nrow * nbase; row++) data[row] = output[row];
	  cpl_free (output);
	} else {
        if (cpl_table_has_column (oi_vis, output_name))
            cpl_table_erase_column (oi_vis, output_name);
        cpl_table_wrap_double_complex (oi_vis, output, output_name);
	}
	CPLCHECK_MSG ("Cannot Wrap");
	
  } else {
	return cpl_error_set_message (cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                  "This type is not supported..."
                                  "report this error to DRS team !!");
  }

  gravi_msg_function_exit(0);
  return CPL_ERROR_NONE;
}

cpl_error_code gravi_table_smooth_column (cpl_table * oi_vis,
										  const char *input_name,
										  const char *output_name,
										  int nsmooth, int nbase)
{
  gravi_msg_function_start(0);
  cpl_ensure_code (oi_vis,      CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (input_name,  CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (output_name, CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (nsmooth>0,   CPL_ERROR_ILLEGAL_INPUT);
  cpl_ensure_code (nbase>0,     CPL_ERROR_ILLEGAL_INPUT);

  /* Performed the running integration and devide */
  gravi_table_runint_column (oi_vis, input_name, output_name, nsmooth, nbase);
  cpl_table_divide_scalar (oi_vis, output_name, 2.0*(double)nsmooth+1.0);
  
  gravi_msg_function_exit(0);
  return CPL_ERROR_NONE;
}

cpl_array * gravi_table_create_sigma_array (cpl_table * oi_wave)
{
  int nv = 0;
  cpl_size size = cpl_table_get_nrow (oi_wave);
  cpl_array * sigma = cpl_array_new (size, CPL_TYPE_DOUBLE);
  
  for (cpl_size w = 0; w<size; w++ )
	cpl_array_set (sigma, w, 1./cpl_table_get (oi_wave, "EFF_WAVE", w, &nv));
  
  CPLCHECK_NUL("Cannot compute the sigma array from EFF_WAVE");
  return sigma;
}

cpl_array * gravi_table_create_wave_array (cpl_table * oi_wave)
{
  cpl_ensure (oi_wave, CPL_ERROR_NULL_INPUT, NULL);
  
  int nv = 0;
  cpl_size size = cpl_table_get_nrow (oi_wave);
  cpl_array * wave = cpl_array_new (size, CPL_TYPE_FLOAT);
  
  for (cpl_size w = 0; w<size; w++ )
	cpl_array_set (wave, w, cpl_table_get (oi_wave, "EFF_WAVE", w, &nv));
  
  CPLCHECK_NUL("Cannot compute the wave array from EFF_WAVE");
  return wave;
}


cpl_error_code gravi_array_normalize_complex (cpl_array * input)
{
  cpl_ensure_code (input, CPL_ERROR_NULL_INPUT);
  
  int nv = 0;
  cpl_size size = cpl_array_get_size (input);

  double complex cpx = 0.0 * I + 0.0;
  for (cpl_size wave = 0; wave<size; wave++ ) {
	cpx = cpl_array_get_complex (input,wave,&nv);
	cpl_array_set_complex (input, wave, cpx / cabs (cpx));
	if (nv) {cpl_array_set_invalid (input, wave); nv=0;}
  }

  CPLCHECK_MSG("Cannot normalize complex");
  return CPL_ERROR_NONE;
}

cpl_array * gravi_array_create_inverse (cpl_array *input)
{
  cpl_ensure (input, CPL_ERROR_NULL_INPUT, NULL);
  
  cpl_size size = cpl_array_get_size (input);
  
  cpl_array * inv = cpl_array_new (size, CPL_TYPE_DOUBLE);
  cpl_array_fill_window (inv, 0, size, 1.0);
  cpl_array_divide (inv, input);

  CPLCHECK_NUL("Cannot compute the sigma array from wave");
  return inv;
}

/*
 * Unwrap an array of phase (in-place)
 */
cpl_error_code gravi_array_phase_unwrap (cpl_array * input)
{
    cpl_ensure_code (input, CPL_ERROR_NULL_INPUT);

	cpl_size size_x = cpl_array_get_size (input);
	double phi_i, phi_ii, d_phi;
	int k_wrap =0;

	for (cpl_size i_data = 1; i_data < size_x; i_data++){

		/* Evaluation of the phi(i_data) and
		 * phi(i_data-1) */

		phi_i = cpl_array_get(input, i_data, NULL);
		phi_ii = cpl_array_get(input, i_data-1, NULL);

		d_phi = phi_i + k_wrap * 2 * M_PI - phi_ii;

		if (d_phi > M_PI) k_wrap ++;
		if (d_phi < - M_PI) k_wrap --;

		cpl_array_set(input, i_data, cpl_array_get(input, i_data, NULL) + k_wrap * 2 * M_PI);
	}
	
	CPLCHECK_INT("Cannot unwrap the phase array");
    return CPL_ERROR_NONE;
}

/*
 * Fill in place the array with carg (cexp (1*I*array))
 */
cpl_error_code gravi_array_phase_wrap (cpl_array * input)
{
  cpl_ensure_code (input, CPL_ERROR_NULL_INPUT);
  
  cpl_size size = cpl_array_get_size (input);

  for (cpl_size wave = 0; wave<size; wave++ ) {
	cpl_array_set (input, wave, carg (cexp ( 1.0 * I * cpl_array_get (input, wave, NULL))));
  }
  
  CPLCHECK_MSG("Cannot wrap phase array");
  return CPL_ERROR_NONE;
}

/**
 * @brief Compute the complex exponention of an array:
 *        cexp (factor * input)
 */
cpl_array * gravi_array_cexp (double complex factor, const cpl_array * input)
{
  cpl_ensure (input, CPL_ERROR_NULL_INPUT, NULL);
  
  cpl_size size = cpl_array_get_size (input);
  
  cpl_array * output = cpl_array_new (size, CPL_TYPE_DOUBLE_COMPLEX);
  cpl_array_fill_window_complex (output, 0, size, (double complex)(0.0 + 0.0 * I));

  cpl_size n;
  int nv = 0.0;
  for (n = 0; n < size; n ++) {
	cpl_array_set_complex (output, n, cexp (factor * cpl_array_get (input, n, &nv) ) );
  }

  return output;
}

/**
 * @brief Multiply a REAL phase to a COMPLEX array, in-place:
 *        input = input * cexp (factor * phase)
 */
cpl_error_code gravi_array_multiply_phasor (cpl_array * input, double complex factor, cpl_array * phase)
{
  cpl_ensure_code (input, CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (phase, CPL_ERROR_NULL_INPUT);
  
  cpl_size size_in = cpl_array_get_size (input);
  cpl_size size_ph = cpl_array_get_size (phase);

  cpl_ensure_code (size_in == size_ph, CPL_ERROR_ILLEGAL_INPUT);
  
  cpl_size n;
  int nv = 0.0;
  for (n = 0; n < size_in; n ++) {
	cpl_array_set_complex( input, n,
						   cpl_array_get_complex( input, n, &nv) *
						   cexp ( factor * cpl_array_get( phase, n, &nv) ) );
  }

  CPLCHECK_MSG ("Cannot multiply phasor");
  return CPL_ERROR_NONE;
}

/**
 * @brief Add a REAL phase to a COMPLEX array, in-place:
 *        input = input + cexp (factor * phase)
 */
cpl_error_code gravi_array_add_phasor (cpl_array * input, double complex factor, cpl_array * phase)
{
  cpl_ensure_code (input, CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (phase, CPL_ERROR_NULL_INPUT);
  
  cpl_size size_in = cpl_array_get_size (input);
  cpl_size size_ph = cpl_array_get_size (phase);

  cpl_ensure_code (size_in == size_ph, CPL_ERROR_ILLEGAL_INPUT);
  
  cpl_size n;
  int nv = 0.0;
  for (n = 0; n < size_in; n ++) {
	cpl_array_set_complex( input, n,
						   cpl_array_get_complex( input, n, &nv) +
						   cexp ( factor * cpl_array_get( phase, n, &nv) ) );
    CPLCHECK_MSG ("Cannot add phasor");
  }

  return CPL_ERROR_NONE;
}

/**
 * @brief Add a pair of COMPLEX arrays to a COMPLEX array, in-place:
 *        input = input + add*conj(sub)
 */
cpl_error_code gravi_array_add_phasors (cpl_array * input, cpl_array * add, cpl_array * sub)
{
  cpl_ensure_code (input, CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (add,   CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (sub,   CPL_ERROR_NULL_INPUT);
  
  cpl_size size_in  = cpl_array_get_size (input);
  cpl_size size_add = cpl_array_get_size (add);
  cpl_size size_sub = cpl_array_get_size (sub);

  cpl_ensure_code (size_in == size_add, CPL_ERROR_ILLEGAL_INPUT);
  cpl_ensure_code (size_in == size_sub, CPL_ERROR_ILLEGAL_INPUT);
  
  cpl_size n;
  int nv = 0.0;
  for (n = 0; n < size_in; n ++) {
	cpl_array_set_complex (input, n,
						   cpl_array_get_complex (input, n, &nv) +
						   cpl_array_get_complex (add, n, &nv) *
                           conj (cpl_array_get_complex (sub, n, &nv)));
    CPLCHECK_MSG ("Cannot add phasors");
  }

  return CPL_ERROR_NONE;
}

/**
 * @brief Add a REAL phase to a REAL phase array, in-place:
 *        input = input + factor * phase
 */
cpl_error_code gravi_array_add_phase (cpl_array * input, double factor, cpl_array * phase)
{
  cpl_ensure_code (input, CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (phase, CPL_ERROR_NULL_INPUT);
  
  cpl_size size_in = cpl_array_get_size (input);
  cpl_size size_ph = cpl_array_get_size (phase);

  cpl_ensure_code (size_in == size_ph, CPL_ERROR_ILLEGAL_INPUT);

  cpl_array * tmp = cpl_array_duplicate (phase);
  
  cpl_array_multiply_scalar (tmp, factor);
  cpl_array_add (input, tmp);

  cpl_array_delete (tmp);
  
  CPLCHECK_MSG("Cannot add phase");
  return CPL_ERROR_NONE;
}

cpl_error_code gravi_array_multiply_conj (cpl_array * input1, cpl_array * input2)
{
  cpl_ensure_code (input1, CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (input2, CPL_ERROR_NULL_INPUT);

  cpl_size size_in1 = cpl_array_get_size (input1);
  cpl_size size_in2 = cpl_array_get_size (input2);

  cpl_ensure_code (size_in1 == size_in2, CPL_ERROR_ILLEGAL_INPUT);

  for (cpl_size w = 0 ; w < size_in1 ; w++)
	cpl_array_set_complex (input1, w, cpl_array_get_complex (input1, w, NULL) *
						   conj(cpl_array_get_complex (input2, w, NULL)));

  CPLCHECK_MSG ("Cannot multiply by conj");
  return CPL_ERROR_NONE;
}


/* 
 * Returned a smoothed version of an array of DOUBLE
 * Invalid data are not checked (all considered as valid)
 */
cpl_array * gravi_array_smooth (cpl_array * input, int nsmooth)
{
  cpl_ensure (nsmooth>0, CPL_ERROR_ILLEGAL_INPUT, NULL);
  cpl_ensure (input,     CPL_ERROR_NULL_INPUT, NULL);
  
  cpl_type type  = cpl_array_get_type (input);
  cpl_size row, n_row = cpl_array_get_size ( input );
  int row_add, row_sub;

  /* Allocate memory of the output -- no check of valid/invalid data */
  cpl_array * output = cpl_array_duplicate (input);

  if ( type == CPL_TYPE_DOUBLE ) {
	double * data_input  = cpl_array_get_data_double (input);
	double * data_output = cpl_array_get_data_double (output);
	double buffer = 0.0;
	int norm=0; // counter number of summed elements
	
	/* Running buffer -- FIXME: check arount limits */
	for ( row=-nsmooth ; row<n_row ; row++) {
	  row_add = row+nsmooth;
	  row_sub = row-nsmooth-1;
	  if ( row_add < n_row ) {
		  buffer += data_input[row_add];
		  norm+=1;
	  }
	  if ( row_sub >= 0 ) {
		  buffer -= data_input[row_sub];
		  norm-=1;
	  }
	  if( row>=0 && row<n_row ) data_output[row] = (double)(buffer)/norm;
	}
	/* End loop on rows */
	
  } else {
	cpl_error_set_message (cpl_func, CPL_ERROR_ILLEGAL_INPUT,
						   "This type is not supported... report this error to DRS team !!");
	return NULL;
  }

  return output;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Optimized computation of GDELAY for a list of arrays
 *
 * @param input        Pointer to arrays (complex coherent flux)
 * @param sigma        Array with the wavenumbers [m]
 * @param gd           Pointer to output GDELAYs  [m]
 * @param nrow         Size of input and gd
 * @param max_width    Maximum opd to explore     [m]
 * @param verbose      Flag to dump the computation time in log
 *
 * For each input array, the group-delay in [m] as the maximum of |Env(x)|
 * where Env(x) = < visdata(lbd) * exp(2i.pi * x / lbd) >
 * with <> sum over lbd
 * 
 * This is optimized to run on a column of arrays:
 * - guess the delay with the interspectra
 * - crude search over the max_width, resolution 2lbd
 * - fine search resolution lbd/10
 * - fine search resolution lbd/100
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_array_get_group_delay_loop (cpl_array ** input, cpl_array * sigma,
                                                 double * gd, cpl_size nrow,
                                                 double max_width,
                                                 int verbose)
{
  gravi_msg_function_start(verbose);
  cpl_ensure_code (input, CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (sigma, CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (gd, CPL_ERROR_ILLEGAL_OUTPUT);
  
  int nv = 0;
  cpl_size nsigma = cpl_array_get_size (sigma);
  double width1, step1, width2, step2, width3, step3;
  cpl_size w, s, nstep1, nstep2, nstep3;
  double x, gd0, gd1, gd2, gd3, current_max = -1.0;
  double lbd = 1.0 / cpl_array_get (sigma,nsigma/2,&nv);
  
  /* Width of a single spectral channel in [m] */
  double coherence = 0.5 * nsigma / fabs (cpl_array_get (sigma,0,&nv) - cpl_array_get (sigma,nsigma-1,&nv));

  /* We never explore more than max_width */
  width1 = CPL_MIN (coherence, max_width);
  step1  = 2.0 * lbd;
  nstep1 = (cpl_size)(width1/step1);

  /* Second pass */
  width2 = 3.0 * step1;
  step2  = 0.1 * lbd;
  nstep2 = (cpl_size)(width2/step2);

  /* Third pass */
  width3 = 3.0 * step2;
  step3  = 0.01 * lbd;
  nstep3 = (cpl_size)(width3/step3);

  /* Allocate memory for the grid search */
  double P = 0.0;
  double * sigdata  = cpl_malloc (sizeof(double complex) * nsigma);
  double complex * visdata  = cpl_malloc (sizeof(double complex) * nsigma);
  double complex * waveform1 = cpl_malloc (sizeof(double complex) * (nstep1+2) * nsigma);
  double complex * waveform2 = cpl_malloc (sizeof(double complex) * (nstep2+2) * nsigma);
  double complex * waveform3 = cpl_malloc (sizeof(double complex) * (nstep3+2) * nsigma);
  
  /* Copy data as double to secure their type */
  for (w=0; w<nsigma; w++) sigdata[w] = cpl_array_get (sigma, w, &nv);

  /* Sum of wavenumber differences */
  double ds = 0.0;
  for (w=1; w<nsigma; w++) ds += sigdata[w] - sigdata[w-1];
  ds /= (nsigma-1);

  /* Build waveform */
  cpl_msg_debug (cpl_func, "Build waveform for 3 pass -- %lli %lli %lli steps", nstep1, nstep2, nstep3);
  
  for (s=0, x = -width1/2.0; x < +width1/2.0; x+=step1)
	for (w=0; w<nsigma; w++) { waveform1[s] = cexp (-2.*I*CPL_MATH_PI * x * sigdata[w]); s++;}
  for (s=0, x = -width2/2.0; x < +width2/2.0; x+=step2)
	for (w=0; w<nsigma; w++) { waveform2[s] = cexp (-2.*I*CPL_MATH_PI * x * sigdata[w]); s++;}
  for (s=0, x = -width3/2.0; x < +width3/2.0; x+=step3) 
	for (w=0; w<nsigma; w++) { waveform3[s] = cexp (-2.*I*CPL_MATH_PI * x * sigdata[w]); s++;}
  
  cpl_msg_debug (cpl_func, "Loop on %lli rows to compute gdelay", nrow);
  
  /* Loop on rows */
  for (cpl_size row = 0; row<nrow; row++) {
    
    /* Copy data as double complex to secure their type and allow
     * in-place modification between the different grids */
    for (w=0; w<nsigma; w++) visdata[w] = cpl_array_get_complex (input[row], w, &nv);

    /* IOTA method in [m] -- as starting point, with 
     * equal weight to all channels to avoid badpixels */
    double complex is = 0.0 + I * 0.0;
    for (w=1; w<nsigma; w++) {
      is += visdata[w] * conj(visdata[w-1]) / CPL_MAX(cabs(visdata[w]) * cabs(visdata[w-1]), 1e-15);
    }
    gd0 = carg (is) / ds / CPL_MATH_2PI;

	/* Remove GD */
	for (w=0; w<nsigma; w++) visdata[w] *= cexp (-2.*I*CPL_MATH_PI*gd0*sigdata[w]);
    
	/* Loop on x to find the maximum of P(x) = |FT(input(sigma))| */
	for (current_max = -1.0, s = 0, x = -width1/2; x < +width1/2; x+=step1) {
	  double complex tmp = 0.0 * I + 0.0;
	  for (w=0; w<nsigma; w++) {tmp += visdata[w] * waveform1[s]; s++;}
	  P = cabs (tmp);
	  if ( P > current_max) { current_max = P; gd1 = x; }
	}

	/* Remove GD */
	for (w=0; w<nsigma; w++) visdata[w] *= cexp (-2.*I*CPL_MATH_PI*gd1*sigdata[w]);

	/* Loop on x to find the maximum of P(x) = |FT(input(sigma))| */
	for (current_max = -1.0, s = 0, x = -width2/2; x < +width2/2; x+=step2) {
	  double complex tmp = 0.0 * I + 0.0;
	  for (w=0; w<nsigma; w++) {tmp += visdata[w] * waveform2[s]; s++;}
	  P = cabs (tmp);
	  if ( P > current_max) { current_max = P; gd2 = x; }
	}

	/* Remove GD */
	for (w=0; w<nsigma; w++) visdata[w] *= cexp (-2.*I*CPL_MATH_PI*gd2*sigdata[w]);

	/* Loop on x to find the maximum of P(x) = |FT(input(sigma))| */
	for (current_max = -1.0, s = 0, x = -width3/2; x < +width3/2; x+=step3) {
	  double complex tmp = 0.0 * I + 0.0;
	  for (w=0; w<nsigma; w++) {tmp += visdata[w] * waveform3[s]; s++;}
	  P = cabs (tmp);
	  if ( P > current_max) { current_max = P; gd3 = x; }
	}
	
	gd[row] = gd0 + gd1 + gd2 + gd3;
	
	CPLCHECK_MSG("Cannot compute GD");
  } /* End loop on rows */

  /* Clean memory */
  FREE (cpl_free, visdata);
  FREE (cpl_free, sigdata);
  FREE (cpl_free, waveform1);
  FREE (cpl_free, waveform2);
  FREE (cpl_free, waveform3);

  gravi_msg_function_exit(verbose);
  return CPL_ERROR_NONE;
}

cpl_error_code gravi_table_compute_group_delay (cpl_table * table, const char *input,
												const char *output, cpl_table * oi_wave)
{
  gravi_msg_function_start(0);
  cpl_ensure_code (table,   CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (input,   CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (output,  CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (oi_wave, CPL_ERROR_NULL_INPUT);

  /* Create sigma */
  cpl_size nwave = cpl_table_get_nrow (oi_wave);
  cpl_array * sigma = cpl_array_new (nwave, CPL_TYPE_DOUBLE);
  for (cpl_size wave = 0; wave < nwave ; wave ++)
	  cpl_array_set (sigma, wave, 1./cpl_table_get (oi_wave, "EFF_WAVE", wave, NULL));

  /* If output columns doesn't exist, create it */
  gravi_table_new_column (table, output, "m", CPL_TYPE_DOUBLE);
  double * gdelay = cpl_table_get_data_double (table, output);

  /* Get data */
  cpl_array ** input_arrays = cpl_table_get_data_array (table, input);
  cpl_size nrow = cpl_table_get_nrow (table);

  CPLCHECK_MSG ("Cannot get data");

  /* Run */
  gravi_array_get_group_delay_loop (input_arrays, sigma, gdelay,
                                    nrow, 1.e-3, CPL_TRUE);
  FREE (cpl_array_delete, sigma);

  gravi_msg_function_exit(0);
  return CPL_ERROR_NONE;
}



/*----------------------------------------------------------------------------*/
/**
 * @brief    Check if two tables have the same content.
 * @param	first table to compare
 * @param    second table to compare
 * @return   1 / 0
 */
/*----------------------------------------------------------------------------*/
int gravi_table_are_equal (cpl_table * first, cpl_table * second)
{
  gravi_msg_function_start(0);
  cpl_ensure_code (first,  CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (second, CPL_ERROR_NULL_INPUT);
  
  cpl_size nrow = cpl_table_get_nrow (first);
  cpl_size ncol = cpl_table_get_ncol (first);
  int nv = 0;
  
  if (nrow != cpl_table_get_nrow (second)) {cpl_msg_info(cpl_func, "Different rows"); return 0;}
  if (ncol != cpl_table_get_ncol (second)) {cpl_msg_info(cpl_func, "Different cols"); return 0;}
  if (cpl_table_compare_structure (first, second)) {cpl_msg_info(cpl_func, "Different structure"); return 0;}

  cpl_array *names = cpl_table_get_column_names (first);

  /* Loop on colums */
  for (cpl_size c = 0; c<ncol; c++) {
	const char * name = cpl_array_get_string (names, c);
	cpl_msg_debug (cpl_func,"Now test %s", name);

    /* Cast the type into an int, to avoid warnings */
	int type = cpl_table_get_column_type (first, name);

	switch (type) {
	case CPL_TYPE_STRING:
	  for (cpl_size r = 0; r<nrow; r++)
		if (strcmp (cpl_table_get_string (first,  name, r),
					cpl_table_get_string (second, name, r) ) )
		  {cpl_msg_info (cpl_func,"Different values in column %s, row %lli", name,r); return 0;}
	  break;
	case CPL_TYPE_DOUBLE:
	case CPL_TYPE_FLOAT:
	case CPL_TYPE_INT:
	  for (cpl_size r = 0; r<nrow; r++)
		if (cpl_table_get (first,  name, r, &nv) !=
			cpl_table_get (second, name, r, &nv) )
		  {cpl_msg_info (cpl_func,"Different values in column %s, row %lli", name,r); return 0;}
	  break;
	case CPL_TYPE_POINTER|CPL_TYPE_STRING:
	  for (cpl_size r = 0; r<nrow; r++) 
		for (cpl_size p = 0; p<cpl_table_get_column_depth (first,name); p++) 
		  if ( strcmp (cpl_array_get_string (cpl_table_get_array (first,  name, r), p),
					   cpl_array_get_string (cpl_table_get_array (second, name, r), p)) )
			{cpl_msg_info (cpl_func,"Different values in column %s, row %lli", name,r); return 0;}
	  break;
	case CPL_TYPE_POINTER|CPL_TYPE_DOUBLE:
	case CPL_TYPE_POINTER|CPL_TYPE_FLOAT:
	case CPL_TYPE_POINTER|CPL_TYPE_INT:
	  for (cpl_size r = 0; r<nrow; r++) 
		for (cpl_size p = 0; p<cpl_table_get_column_depth (first,name); p++)  {
		  if (cpl_array_get (cpl_table_get_array (first,  name, r), p, NULL) !=
			  cpl_array_get (cpl_table_get_array (second, name, r), p, NULL) )
			{cpl_msg_info (cpl_func,"Different values in column %s, row %lli", name,r); return 0;}
		}
	  break;
	default:
	  cpl_error_set_message (cpl_func, CPL_ERROR_ILLEGAL_INPUT, "Cannot compare these tables (TBD, FIXME)");
	  return 0;
	}

  }/* End loop on columns */

  gravi_msg_function_exit(0);
  return 1;
}

cpl_error_code gravi_table_new_column (cpl_table * table, const char * name, const char * unit, cpl_type type)
{
  cpl_ensure_code (table, CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (name,  CPL_ERROR_NULL_INPUT);

  if ( cpl_table_has_column (table, name) &&
	   cpl_table_get_column_type (table, name) == type) {
	cpl_msg_info (cpl_func, "Column %s already exists", name);
  } else {
	cpl_table_new_column (table, name, type);
  }

  if (type == CPL_TYPE_DOUBLE_COMPLEX || type == CPL_TYPE_FLOAT_COMPLEX)
      cpl_table_fill_column_window_complex (table, name, 0, cpl_table_get_nrow (table), 0.0 + I*0.0);
  else
      cpl_table_fill_column_window (table, name, 0, cpl_table_get_nrow (table), 0.0);
  
  if (unit) cpl_table_set_column_unit (table, name, unit);
  
  return CPL_ERROR_NONE;
}

cpl_error_code gravi_table_new_column_array (cpl_table * table, const char * name, const char * unit, cpl_type type, cpl_size size)
{
  cpl_ensure_code (table, CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (name,  CPL_ERROR_NULL_INPUT);

  if ( cpl_table_has_column (table, name) )
      cpl_table_erase_column (table, name);
  
  cpl_table_new_column_array (table, name, type, size);
  if (unit) cpl_table_set_column_unit (table, name, unit);
  
  return CPL_ERROR_NONE;
}

cpl_error_code gravi_table_init_column_array (cpl_table * table, const char * name, const char * unit, cpl_type type, cpl_size size)
{
  cpl_ensure_code (table, CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (name,  CPL_ERROR_NULL_INPUT);

  if ( cpl_table_has_column (table, name) )
      cpl_table_erase_column (table, name);
  
  cpl_table_new_column_array (table, name, type, size);
  if (unit) cpl_table_set_column_unit (table, name, unit);

  cpl_array * array = cpl_array_new (size, type);
  cpl_array_fill_window (array, 0, size, 0.0);

  cpl_size nrow  = cpl_table_get_nrow (table);
  for (cpl_size row = 0; row < nrow; row++)
      cpl_table_set_array (table, name, row, array);
  
  FREE (cpl_array_delete, array);
  
  return CPL_ERROR_NONE;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Unwrap an imagelist an all its images
 *
 * @param imglist   The imagelist to unwrap
 *
 * All the images of imaglist are unwrapped, that is the data remain allocated
 * but the image structure is deleted (see cpl_image_unwrap).
 * Then the imagelist structure itself is deleted with cpl_imagelist_delete.
 */
/*---------------------------------------------------------------------------*/

cpl_error_code gravi_imagelist_unwrap_images (cpl_imagelist * imglist)
{
  cpl_ensure_code (imglist, CPL_ERROR_NULL_INPUT);
  
  cpl_size nrow = cpl_imagelist_get_size (imglist);
  
  for (cpl_size i = nrow-1; i>=0 ; i--) {
	cpl_image_unwrap (cpl_imagelist_unset (imglist, i));
	CPLCHECK_MSG("Cannot unset image");
  }

  cpl_imagelist_delete (imglist);
  CPLCHECK_MSG("Cannot delete imagelist");

  return CPL_ERROR_NONE;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Wrap a column array of a table into an imagelist
 *
 * @param table_data   The input table
 * @param data_x       The column name
 *
 * The array of each row are wrapped into image of their type, and append
 * into the returned imagelist. The returned pointer shall be deleted
 * with gravi_imagelist_unwrap_images
 */
/*---------------------------------------------------------------------------*/

cpl_imagelist * gravi_imagelist_wrap_column (cpl_table * table_data, const char * data_x)
{
	gravi_msg_function_start(0);
    cpl_ensure (table_data, CPL_ERROR_NULL_INPUT, NULL);
    cpl_ensure (data_x,     CPL_ERROR_NULL_INPUT, NULL);

	/* Get the type of the column 
     * Cast into an int to avoid warnings */
    int type_column = cpl_table_get_column_type (table_data, data_x);

    /* Get pointer to the data */
	cpl_size nrow  = cpl_table_get_nrow (table_data);
    cpl_array ** array = cpl_table_get_data_array (table_data, data_x);
    cpl_ensure (array, CPL_ERROR_ILLEGAL_INPUT, NULL);
    
    /* If the column has no-dimension, we fake them */
    cpl_size nx, ny;
    if (cpl_table_get_column_dimensions (table_data, data_x)<2) {
        nx = cpl_table_get_column_depth (table_data, data_x);
        ny = 1;
    } else {
        nx = cpl_table_get_column_dimension (table_data, data_x, 0);
        ny = cpl_table_get_column_dimension (table_data, data_x, 1);
    }
	CPLCHECK_NUL ("Cannot get dimension");

    /* Create output */
    cpl_image * img;
	cpl_imagelist * imglist = cpl_imagelist_new();
	
	/* compute the image list depending of the DATA type */
	switch (type_column)
	{
	    case CPL_TYPE_POINTER|CPL_TYPE_DOUBLE :
		  
	        for (cpl_size j = 0; j < nrow ; j++)
	        {
	            img = cpl_image_wrap_double (nx, ny, cpl_array_get_data_double(array[j]));
	            cpl_imagelist_set (imglist, img, j);
	        }

	    break;

	    case CPL_TYPE_POINTER|CPL_TYPE_INT :

	        for (cpl_size j = 0; j < nrow ; j++)
	        {
	            img = cpl_image_wrap_int (nx, ny, cpl_array_get_data_int(array[j]));
	            cpl_imagelist_set (imglist, img,j);
	        }

	    break;

	    case CPL_TYPE_POINTER|CPL_TYPE_FLOAT :
	        for (cpl_size j = 0; j < nrow ; j++)
	        {
	            img = cpl_image_wrap_float (nx, ny, cpl_array_get_data_float(array[j]));
	            cpl_imagelist_set (imglist, img, j);
	        }

	    break;

	    default:

	        cpl_error_set_message(cpl_func, CPL_ERROR_INVALID_TYPE,
	        		       "invalid type of image coming from %s", data_x);
	        cpl_imagelist_delete (imglist);
	        return NULL;
	    break;
	}

	gravi_msg_function_exit(0);
	return imglist;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Estimate the median noise in a window. This noise is estimated
 *        as the median{img**2}**0.5, so that it is a robust estimate.
 */
/*----------------------------------------------------------------------------*/

double gravi_image_get_noise_window (cpl_image *img,
                                     cpl_size llx, cpl_size lly,
                                     cpl_size urx, cpl_size ury)
{
    gravi_msg_function_start(0);
    cpl_ensure (img, CPL_ERROR_NULL_INPUT, -1);
    int nv;

    /* Extract values in vector */
    cpl_vector * flux = cpl_vector_new ((urx-llx+1)*(ury-lly+1));

    for (cpl_size v = 0, x = llx; x <= urx; x++) {
        for (cpl_size y = lly; y <= ury; y++) {
            cpl_vector_set (flux, v, cpl_image_get (img, x, y, &nv));
            v++;
            CPLCHECK_MSG ("Cannot fill vector");
        }
    }

    /* FIXME: remove the median before square */

    /* Compute typical error as the
     * median of spatial variation */
    cpl_vector_multiply (flux, flux);
    
    double RMS = sqrt (cpl_vector_get_median (flux));
    FREE (cpl_vector_delete, flux);
    
    gravi_msg_function_exit(0);
    return RMS;
}

/*----------------------------------------------------------------------------*/

cpl_error_code gravi_image_subtract_window (cpl_image * img1, const cpl_image * img2,
                                            cpl_size llx, cpl_size lly,
                                            cpl_size urx, cpl_size ury,
                                            cpl_size llx2, cpl_size lly2)
{
	gravi_msg_function_start(0);
    cpl_ensure_code (img1, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code (img2, CPL_ERROR_NULL_INPUT);

    /* Ensure size */
    cpl_size nx = cpl_image_get_size_x (img1);
    cpl_size ny = cpl_image_get_size_y (img1);
    urx = CPL_MIN (urx, nx);
    ury = CPL_MIN (ury, ny);
    llx2 -= llx;
    lly2 -= lly;

    int nv;
    for (cpl_size x=llx; x<=urx; x++) {
        for (cpl_size y=lly; y<=ury; y++) {
            cpl_image_set (img1, x, y,
                           cpl_image_get (img1,x,y,&nv) -
                           cpl_image_get (img2,x+llx2,y+lly2,&nv));
        }
    }
    
	gravi_msg_function_exit(0);
	return CPL_ERROR_NONE;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Create an image from a column array in table
 *
 * @param table_data   The input table
 * @param data_x       The column name
 * @param row          The row to extract
 *
 * The array is copied into image of same type. The returned pointer
 * shall be deleted with cpl_image_delete
 *
 * FIXME: to be optimized
 */
/*---------------------------------------------------------------------------*/

cpl_image * gravi_image_from_column (cpl_table * table_data, const char * data_x, cpl_size row)
{
	gravi_msg_function_start(0);
    cpl_ensure (table_data, CPL_ERROR_NULL_INPUT, NULL);
    cpl_ensure (data_x,     CPL_ERROR_NULL_INPUT, NULL);
    cpl_ensure (row < cpl_table_get_nrow (table_data), CPL_ERROR_ILLEGAL_INPUT, NULL);
    
    cpl_imagelist * wrap_imglist = gravi_imagelist_wrap_column (table_data, data_x);
    cpl_ensure (wrap_imglist, CPL_ERROR_ILLEGAL_INPUT, NULL);

    cpl_image * out_img = cpl_image_duplicate (cpl_imagelist_get (wrap_imglist, row));
    gravi_imagelist_unwrap_images (wrap_imglist);
    
	gravi_msg_function_exit(0);
	return out_img;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Create an imagelist from a column array in table
 *
 * @param table_data   The input table
 * @param data_x       The column name
 *
 * The array of each row are copied into image of their type, and append
 * into the returned imagelist. The returned pointer shall be deleted
 * with cpl_imagelist_delete
 *
 * FIXME: to be optimized
 */
/*---------------------------------------------------------------------------*/

cpl_imagelist * gravi_imagelist_from_column (cpl_table * table_data, const char * data_x)
{
	gravi_msg_function_start(0);
    cpl_ensure (table_data, CPL_ERROR_NULL_INPUT, NULL);
    cpl_ensure (data_x,     CPL_ERROR_NULL_INPUT, NULL);
    
    cpl_imagelist * wrap_imglist = gravi_imagelist_wrap_column (table_data, data_x);
    cpl_ensure (wrap_imglist, CPL_ERROR_ILLEGAL_INPUT, NULL);
    
    cpl_imagelist * out_imglist = cpl_imagelist_duplicate (wrap_imglist);
    gravi_imagelist_unwrap_images (wrap_imglist);
    
	gravi_msg_function_exit(0);
	return out_imglist;
}


/*---------------------------------------------------------------------------*/
/**
 * @brief Return an array ready for cpl_table_set_column_dimension
 *
 * @param table        The input table
 * @param name         The column name
 *
 * The return array (CPL_TYPE_INT) is of size the number returned
 * by cpl_table_get_column_dimensions, and each values corresponds to
 * the dimension returned by cpl_table_get_column_dimension
 */
/*---------------------------------------------------------------------------*/

cpl_array * gravi_table_get_column_dimension (const cpl_table * table, const char * name)
{
    cpl_ensure (table, CPL_ERROR_NULL_INPUT, NULL);
    cpl_ensure (name,  CPL_ERROR_NULL_INPUT, NULL);
    
    cpl_size ndim = cpl_table_get_column_dimensions (table, name);
    cpl_array * dimension = cpl_array_new (ndim, CPL_TYPE_INT);
    for (cpl_size dim = 0; dim < ndim; dim++) {
        int value = cpl_table_get_column_dimension (table, name, dim);
        cpl_array_set (dimension, dim, value);
    }

    return dimension;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Wrap the data of na image into an array
 *
 * @param img     The image to wrap
 *
 * @return An array with the same type as the image @c img.
 * @note the returned array must unwrap using the function @c cpl_array_unwrap().
 * @see  cpl_array_unwrap().
 */
/*---------------------------------------------------------------------------*/

cpl_array * gravi_array_wrap_image (cpl_image * img)
{
	cpl_ensure (img, CPL_ERROR_NULL_INPUT, NULL);

	cpl_type type_img = cpl_image_get_type (img);
	int x = cpl_image_get_size_x (img);
	int y = cpl_image_get_size_y (img);

	CPLCHECK_NUL ("Cannot get data");

	cpl_array * array = NULL;
	switch (type_img){
	case CPL_TYPE_FLOAT :
		array = cpl_array_wrap_float (cpl_image_get_data_float(img), x*y);
		break;
	case CPL_TYPE_DOUBLE :
		array = cpl_array_wrap_double (cpl_image_get_data_double(img), x*y);
		break;
	case CPL_TYPE_INT :
		array = cpl_array_wrap_int (cpl_image_get_data_int(img), x*y);
		break;
	default :
		cpl_error_set_message(cpl_func, CPL_ERROR_INVALID_TYPE,
			                               "invalid type of image");
		return NULL;
	}

	return array;
}


/*---------------------------------------------------------------------------*/
/**
 * @brief  Copy the content of two vector into a newly allocated 2D matrix
 * @param  vector1
 * @param  vector2
 * 
 * @return  The newly allocated matrix
 * 
 * @note The returned matrix must be deallocated using the function
 * @see cpl_matrix_delete()
 */
/*---------------------------------------------------------------------------*/

cpl_matrix * get_matrix_from_vector(cpl_vector * vector1,
									cpl_vector * vector2)
{
	cpl_ensure (vector1, CPL_ERROR_NULL_INPUT, NULL);
	
	cpl_matrix * matrix, * matrix_wrap;
	int size = cpl_vector_get_size(vector1);
	double * data1, * data2;

	if( vector2 != NULL ){
	    int size2 = cpl_vector_get_size(vector2);
	    cpl_ensure (size == size2, CPL_ERROR_ILLEGAL_INPUT, NULL);
		
		data1 = cpl_malloc(2 * size * sizeof(double));
		memcpy(data1, cpl_vector_get_data(vector1), size * sizeof(double));

		data2 = cpl_vector_get_data(vector2);
		memcpy(data1 + size, data2, size * sizeof(double));

		matrix_wrap = cpl_matrix_wrap(2, size, data1);
		matrix = cpl_matrix_transpose_create(matrix_wrap);

		cpl_matrix_unwrap(matrix_wrap);
	}
	else{
		data1 = cpl_malloc(size * sizeof(double));
		memcpy(data1, cpl_vector_get_data(vector1), size * sizeof(double));
		
		matrix = cpl_matrix_wrap(size, 1, data1);
	}
	cpl_free(data1);

	return matrix;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Create a vector from the row @c index of the column @c regname
 *
 * @param  spectrum_table The table loaded from the SPECTRUM_DATA field
 * @param  regname        The column's name on the table
 * @param  index          The index on the row (0 for a scalar column)
 * 
 * @return An allocated vector with values extracted from the table
 *
 * Create a vector from the row @c index of the column @c regname
 * The returned vector must be deallocated.
 */
/*---------------------------------------------------------------------------*/

cpl_vector * gravi_table_get_vector (cpl_table * spectrum_data,
                                     cpl_size index,
		                             const char * regname)
{
	/* Check the inputs */
	cpl_ensure (spectrum_data, CPL_ERROR_NULL_INPUT, NULL);
	cpl_ensure (regname,       CPL_ERROR_NULL_INPUT, NULL);
	cpl_ensure (index>=0,   CPL_ERROR_ILLEGAL_INPUT, NULL);

	/* Get the array size */
    cpl_size size = cpl_table_get_column_depth (spectrum_data, regname);
	cpl_ensure (index<size, CPL_ERROR_ILLEGAL_INPUT, NULL);

	/* allocate the vector */
	cpl_size nrow = cpl_table_get_nrow (spectrum_data);
	cpl_vector * data_value = cpl_vector_new (nrow);

    if (size > 0) {

        /* Extract the data from the column region */
        cpl_array ** column_array = cpl_table_get_data_array (spectrum_data, regname);
        cpl_ensure (column_array, CPL_ERROR_ILLEGAL_INPUT, NULL);

        /* Loop on row */
        double value;
        for (cpl_size row = 0; row < nrow; row ++){
            value = cpl_array_get (column_array[row], index, NULL);
            cpl_vector_set (data_value, row, value);
        }
    } else if (cpl_table_get_column_type (spectrum_data, regname)
               == CPL_TYPE_DOUBLE) {

        double * data = cpl_table_get_data_double (spectrum_data, regname);
        for (cpl_size row = 0; row < nrow; row++) {
            cpl_vector_set (data_value, row, data[row]);
        }

    } else if (cpl_table_get_column_type (spectrum_data, regname)
               == CPL_TYPE_INT) {

        int * data = cpl_table_get_data_int (spectrum_data, regname);
        for (cpl_size row = 0; row < nrow; row++) {
            cpl_vector_set (data_value, row, data[row]);
        }

    } else {
        cpl_error_set_message (cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                               "This type is not supported"
                               "(report to DRS team).");
        FREE (cpl_vector_delete, data_value);
        return NULL;
    }

    return (data_value);
}

cpl_vector * gravi_table_get_vector_scalar (cpl_table * table,
                                            const char * name,
                                            cpl_size base,
                                            cpl_size nbase)
{
	/* Check the inputs */
	cpl_ensure (table,   CPL_ERROR_NULL_INPUT, NULL);
	cpl_ensure (name,    CPL_ERROR_NULL_INPUT, NULL);
	cpl_ensure (base>=0, CPL_ERROR_ILLEGAL_INPUT, NULL);
	cpl_ensure (nbase>0, CPL_ERROR_ILLEGAL_INPUT, NULL);

	/* Ensure it is scalar */
    cpl_size size = cpl_table_get_column_depth (table, name);
	cpl_ensure (size==0, CPL_ERROR_ILLEGAL_INPUT, NULL);

	/* Allocate the vector */
	cpl_size nrow = cpl_table_get_nrow (table) / nbase;
	cpl_vector * vector = cpl_vector_new (nrow);

    /* Get the type */
    cpl_type type = cpl_table_get_column_type (table, name);

    if (type == CPL_TYPE_DOUBLE) {
        double * data = cpl_table_get_data_double (table, name);
        for (cpl_size row = 0; row < nrow; row++)
            cpl_vector_set (vector, row, data[row*nbase+base]);
    }
    else if (type == CPL_TYPE_FLOAT) {
        float * data = cpl_table_get_data_float (table, name);
        for (cpl_size row = 0; row < nrow; row++)
            cpl_vector_set (vector, row, data[row*nbase+base]);
    }
    else if (type == CPL_TYPE_INT) {
        int * data = cpl_table_get_data_int (table, name);
        for (cpl_size row = 0; row < nrow; row++)
            cpl_vector_set (vector, row, data[row*nbase+base]);
    }
    else {
        cpl_error_set_message (cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                               "This type is not supported"
                               "(report to DRS team).");
        FREE (cpl_vector_delete, vector);
        return NULL;
    }

    return vector;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Create a vector from the row @c index of the column @c regname
 *
 * @param  spectrum_table The table loaded from the SPECTRUM_DATA field
 * @param  index          The index on the row
 * @param  regname1       The column's name on the table
 * @param  regname2       The column's name on the table
 * 
 * @return An allocated vector with values extracted from the table
 *         regname1[index] - regname2[index]
 *
 * Create a vector from the row @c index of the column @c regname
 * The returned vector must be deallocated.
 */
/*---------------------------------------------------------------------------*/

cpl_vector * gravi_table_get_vector_diff (cpl_table * spectrum_data, int index,
                                          const char * regname1,
                                          const char * regname2)
{
	/* Check the inputs */
	cpl_ensure (spectrum_data, CPL_ERROR_NULL_INPUT, NULL);
	cpl_ensure (regname1,      CPL_ERROR_NULL_INPUT, NULL);
	cpl_ensure (regname2,      CPL_ERROR_NULL_INPUT, NULL);
	cpl_ensure (index>=0,      CPL_ERROR_ILLEGAL_INPUT, NULL);

	/* Get the array size */
    cpl_size size1 = cpl_table_get_column_depth (spectrum_data, regname1);
    cpl_size size2 = cpl_table_get_column_depth (spectrum_data, regname2);
	cpl_ensure (index<size1, CPL_ERROR_ILLEGAL_INPUT, NULL);
	cpl_ensure (index<size2, CPL_ERROR_ILLEGAL_INPUT, NULL);

	/* Extract the data from the column region */
	cpl_array ** column_array1 = cpl_table_get_data_array (spectrum_data, regname1);
	cpl_array ** column_array2 = cpl_table_get_data_array (spectrum_data, regname2);
	cpl_ensure (column_array1, CPL_ERROR_ILLEGAL_INPUT, NULL);
	cpl_ensure (column_array2, CPL_ERROR_ILLEGAL_INPUT, NULL);

	/* allocate the vector */
	cpl_size nrow = cpl_table_get_nrow (spectrum_data);
	cpl_vector * data_value = cpl_vector_new (nrow);

    /* Loop on row */
	double value;
	for (cpl_size row = 0; row < nrow; row++){
		value = cpl_array_get (column_array1[row], index, NULL) - cpl_array_get (column_array2[row], index, NULL);
		cpl_vector_set (data_value, row, value);
	}

	return (data_value);
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Fill entire image with value
 */
/*---------------------------------------------------------------------------*/

cpl_error_code gravi_image_fill (cpl_image * img, double value)
{
    cpl_ensure_code (img, CPL_ERROR_NULL_INPUT);

    cpl_size nx = cpl_image_get_size_x (img);
    cpl_size ny = cpl_image_get_size_y (img);
    cpl_image_fill_window (img, 1,1,nx,ny, value);

    return CPL_ERROR_NONE;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Wrap matrix into an image (data not duplicated)
 */
/*---------------------------------------------------------------------------*/

cpl_image * gravi_image_wrap_matrix (cpl_matrix * matrix)
{
    cpl_ensure (matrix, CPL_ERROR_NULL_INPUT, NULL);
    cpl_size nx = cpl_matrix_get_ncol (matrix);
    cpl_size ny = cpl_matrix_get_nrow (matrix);
    return cpl_image_wrap_double (nx, ny, cpl_matrix_get_data (matrix));
}

cpl_image * gravi_image_from_matrix (cpl_matrix * matrix)
{
    cpl_ensure (matrix, CPL_ERROR_NULL_INPUT, NULL);
    cpl_size nx = cpl_matrix_get_ncol (matrix);
    cpl_size ny = cpl_matrix_get_nrow (matrix);

    cpl_image * image = cpl_image_new (nx,ny,CPL_TYPE_DOUBLE);
    for (cpl_size i = 0; i < nx; i++)
        for (cpl_size j = 0; j < ny; j++)
            cpl_image_set (image, i+1, j+1, cpl_matrix_get (matrix,j,i));
    
    return image;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Wrap vector into an image (data not duplicated)
 */
/*---------------------------------------------------------------------------*/

cpl_image * gravi_image_wrap_vector (cpl_vector * vector)
{
    cpl_ensure (vector, CPL_ERROR_NULL_INPUT, NULL);
    cpl_size nx = cpl_vector_get_size (vector);
    cpl_size ny = 1;
    return cpl_image_wrap_double (nx, ny, cpl_vector_get_data (vector));
}

cpl_image * gravi_image_from_vector (cpl_vector * vector)
{
    cpl_ensure (vector, CPL_ERROR_NULL_INPUT, NULL);
    cpl_size nx = cpl_vector_get_size (vector);
    cpl_size ny = 1;
    
    cpl_image * image = cpl_image_new (nx,ny,CPL_TYPE_DOUBLE);
    for (cpl_size i = 0; i < nx; i++)
            cpl_image_set (image, i+1, 1, cpl_vector_get (vector,i));
    
    return image;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Compute the quantile of an image
 *
 * @param img     input image
 * @param thr     quantile fraction (0 < thr < 1)
 * 
 * Compute the value of the image corresponding to the quantile fraction
 * 'thr', that is such that their is nx*ny*thr pixels with values lower
 * the returned value.
 */
/*---------------------------------------------------------------------------*/

double gravi_image_get_quantile (const cpl_image * img, double thr)
{
  cpl_ensure (img, CPL_ERROR_NULL_INPUT, -1);
  cpl_ensure (thr>0 && thr<1, CPL_ERROR_ILLEGAL_INPUT, -1);

  cpl_size nx = cpl_image_get_size_x (img);
  cpl_size ny = cpl_image_get_size_y (img);
  cpl_size nq = (cpl_size)(thr * nx * ny);

  cpl_ensure (nq>=0 && nq<=nx*ny, CPL_ERROR_ILLEGAL_INPUT, -1);

  /* Create vector and fill */
  int nv;
  cpl_vector * vect = cpl_vector_new (nx*ny);
  for (cpl_size ix = 0; ix < nx ; ix++) 
	for (cpl_size iy = 0; iy < ny ; iy++) 
	  cpl_vector_set (vect, ix * ny + iy, cpl_image_get (img, ix+1, iy+1, &nv));

  /* Sort and get quartil */
  cpl_vector_sort (vect, CPL_SORT_ASCENDING);
  double value = cpl_vector_get (vect, nq);
  cpl_vector_delete (vect);
  
  return value;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief    Compute the value of the vector corresponding
 *           to the quantile 'thr' (0 < thr < 1)
 * @return   quantile value
 */
/*---------------------------------------------------------------------------*/

double gravi_array_get_quantile (cpl_array * arr, double thr)
{
  cpl_ensure (arr, CPL_ERROR_NULL_INPUT, -1);
  cpl_ensure (thr>0 && thr<1, CPL_ERROR_ILLEGAL_INPUT, -1);

  cpl_size nx = cpl_array_get_size (arr);
  cpl_size nq = (cpl_size)(thr * nx);

  cpl_ensure (nq>=0 && nq<=nx, CPL_ERROR_ILLEGAL_INPUT, -1);

  /* Create vector and fill */
  int nv;
  cpl_vector * vect = cpl_vector_new (nx);
  for (cpl_size ix = 0; ix < nx ; ix++) 
	cpl_vector_set (vect, ix, cpl_array_get (arr, ix, &nv));

  /* Sort and get quartil */
  cpl_vector_sort (vect, CPL_SORT_ASCENDING);
  double value = cpl_vector_get (vect, nq);
  cpl_vector_delete (vect);

  return value;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Return the index of maximum in a vector. If several indexes
 *        exists with the maximum value, the smallest index is returned.
 * @param vector     Input vector (size)
 * @return The index of maximum or NULL on error
 */
/*---------------------------------------------------------------------------*/

cpl_size gravi_vector_get_maxpos (cpl_vector * vector)
{
    cpl_ensure (vector, CPL_ERROR_NULL_INPUT, -1);
    cpl_size size = cpl_vector_get_size (vector);

    cpl_size pos = 0;
    double value = cpl_vector_get (vector, 0);
    for (int s = 1; s < size; s++) {
        if (cpl_vector_get (vector, s) > value) {
            pos = s;
            value = cpl_vector_get (vector, s);
        }
    }

    return pos;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Return the mean of a vector after extrema and RMS clipping
 * 
 * @param vector_in    input vector
 * @param percent      0-1, fraction of extremes values rejected
 * @param nsigma       sigma clipping
 * 
 * @return The mean of accepted values after
 */
/*---------------------------------------------------------------------------*/
double gravi_vector_get_mean_clip (cpl_vector * vector_in,
                                   double percent,
                                   double nsigma)
{
    cpl_ensure (vector_in,     CPL_ERROR_NULL_INPUT, 0.0);
    cpl_ensure (percent > 0,   CPL_ERROR_ILLEGAL_INPUT, 0.0);
    cpl_ensure (percent < 0.5, CPL_ERROR_ILLEGAL_INPUT, 0.0);
    cpl_ensure (nsigma > 0,    CPL_ERROR_ILLEGAL_INPUT, 0.0);

    /* Sort */
    cpl_vector * sort_vector = cpl_vector_duplicate (vector_in);
    cpl_vector_sort (sort_vector, CPL_SORT_ASCENDING);

    /* Clip extrems values */
    cpl_size size = cpl_vector_get_size (vector_in);
    cpl_size sizeout = size*(1-percent*2);
    cpl_size start = (size-sizeout)/2;

    cpl_vector * vector = cpl_vector_new (sizeout);
    for (cpl_size i = 0 ; i < sizeout ; i++)
        cpl_vector_set (vector, i, cpl_vector_get (sort_vector, i+start));

    /* Clip above several sigmas */
    cpl_vector * vector_med = cpl_vector_new (sizeout);
    double med = cpl_vector_get_median (vector);
    double rms = nsigma * cpl_vector_get_stdev (vector);
    
    cpl_size size_med = 0;
    for (cpl_size i = 0 ; i < cpl_vector_get_size (vector) ; i++)
        if ( (cpl_vector_get (vector, i) > med-rms) &&
             (cpl_vector_get (vector, i) < med+rms) ) {
            cpl_vector_set (vector_med, size_med, cpl_vector_get (vector, i));
            size_med++;
        }
    cpl_vector_set_size (vector_med, size_med);
        
    /* Compute mean of accepted values */
    double output = cpl_vector_get_mean (vector_med);

    FREE (cpl_vector_delete, vector_med);
    FREE (cpl_vector_delete, sort_vector);
    return output;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Extract part of a vector
 * @param vector     Input vector (size)
 * @param start      Starting index (0...size-1)
 * @param step       Step index, so the returned values are every step
 * @return A new vector with size size/step
 */
/*---------------------------------------------------------------------------*/

cpl_vector * gravi_vector_extract (const cpl_vector * vector, int start, int step)
{
    cpl_size size = cpl_vector_get_size (vector);
    cpl_size newsize = size / step;
    cpl_vector * out = cpl_vector_new (newsize);

    for (int s = 0; s < newsize; s++)
        cpl_vector_set (out, s, cpl_vector_get (vector, s*step+start));

    return out;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief    Unwrap a phase vector following a guess vector. The difference
 *           is actually unwrap and shall thus be slowly evolving.
 * 			 The operration if performed in-place.
 * @param    vector to unwrap [rad]
 * @param    ref is the guess 
 * @param    ref_to_phase is the coefficient to convert ref into [rad]
 */
/*---------------------------------------------------------------------------*/

cpl_error_code gravi_vector_unwrap_with_guess (cpl_vector * vector, cpl_vector * ref, double ref_to_phase)
{
  cpl_ensure_code (vector, CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (ref,    CPL_ERROR_NULL_INPUT);

  cpl_size nrow = cpl_vector_get_size (vector);

  double referenced, referenced_prev = 0.0;
  cpl_size wrap = 0;

  for (cpl_size row = 0 ; row < nrow; row ++) {
    double phase_ref = ref_to_phase * cpl_vector_get (ref, row);

    /* Referenced phase in radian in [0,2pi] */
    referenced = (cpl_vector_get (vector, row) - phase_ref);
    referenced = fmod (referenced, CPL_MATH_2PI);
    if (referenced < 0) referenced += CPL_MATH_2PI;
    
    /* Check if referenced_phase is wrapp */
    if ( referenced - referenced_prev >  CPL_MATH_PI ) wrap --;
    if ( referenced - referenced_prev < -CPL_MATH_PI ) wrap ++;
    referenced_prev = referenced;

    /* Set back in-place */
    cpl_vector_set (vector, row, phase_ref + referenced + wrap * CPL_MATH_2PI);
  }

  CPLCHECK_MSG ("Cannot unwrap with guess");

  return CPL_ERROR_NONE;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Multiply scalar or array column by scalar
 *
 * @param table       The table to modify in-place
 * @param name        The name of the column to modify
 * @param base        The starting row
 * @param nbase       The row increment
 * @param value       Value to multiply
 * 
 * Each row in base::nbase (pythonic notation) of the column name 
 * is being multiplied by value.
 */
/*---------------------------------------------------------------------------*/

cpl_error_code gravi_table_multiply_scalar (cpl_table * table, const char * name,
                                            int base, int nbase, double value)
{
  cpl_ensure_code (table, CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (name,  CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (nbase==1 || nbase==4 || nbase==6, CPL_ERROR_ILLEGAL_INPUT);
  cpl_ensure_code (base>=0 && base <nbase,           CPL_ERROR_ILLEGAL_INPUT);

  cpl_size nrow = cpl_table_get_nrow (table) / nbase;

  if (cpl_table_get_column_depth (table, name) > 0) {
      cpl_array ** array = cpl_table_get_data_array (table, name);  
      for (cpl_size row = 0 ; row < nrow ; row ++) {
          cpl_array_multiply_scalar (array[row*nbase+base], value);
          CPLCHECK_MSG ("Cannot multiply (array may not be numerical)");
      }
  } else if(cpl_table_get_column_type (table, name) == CPL_TYPE_DOUBLE) {
      double * array = cpl_table_get_data_double (table, name);
      for (cpl_size row = 0 ; row < nrow ; row ++) {
          array[row*nbase+base] *= value;
      }
  } else if(cpl_table_get_column_type (table, name) == CPL_TYPE_INT) {
      int * array = cpl_table_get_data_int (table, name);
      for (cpl_size row = 0 ; row < nrow ; row ++) {
          array[row*nbase+base] *= value;
      }
  } else {
      return cpl_error_set_message (cpl_func, CPL_ERROR_INVALID_TYPE,
                                    "Column type is not supported");
  }

  CPLCHECK_MSG ("Cannot multiply");
  return CPL_ERROR_NONE;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Multply scalar or array column by scalar
 *
 * @param table       The table to modify in-place
 * @param name        The name of the column to modify
 * @param base        The starting row
 * @param nbase       The row increment
 * @param value       Value to add
 * 
 * Each row in base::nbase (pythonic notation) of the column name 
 * is being added by value.
 */
/*---------------------------------------------------------------------------*/

cpl_error_code gravi_table_add_scalar (cpl_table * table, const char * name,
                                       int base, int nbase, double value)
{
  cpl_ensure_code (table, CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (name,  CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (nbase==1 || nbase==4 || nbase==6, CPL_ERROR_ILLEGAL_INPUT);
  cpl_ensure_code (base>=0 && base <nbase,           CPL_ERROR_ILLEGAL_INPUT);

  cpl_size nrow = cpl_table_get_nrow (table) / nbase;

  if (cpl_table_get_column_depth (table, name) > 0) {
      cpl_array ** array = cpl_table_get_data_array (table, name);  
      for (cpl_size row = 0 ; row < nrow ; row ++) {
          cpl_array_add_scalar (array[row*nbase+base], value);
          CPLCHECK_MSG ("Cannot multiply (array may not be numerical)");
      }
  } else if(cpl_table_get_column_type (table, name) == CPL_TYPE_DOUBLE) {
      double * array = cpl_table_get_data_double (table, name);
      for (cpl_size row = 0 ; row < nrow ; row ++) {
          array[row*nbase+base] += value;
      }
  } else if(cpl_table_get_column_type (table, name) == CPL_TYPE_INT) {
      int * array = cpl_table_get_data_int (table, name);
      for (cpl_size row = 0 ; row < nrow ; row ++) {
          array[row*nbase+base] += value;
      }
  } else {
      return cpl_error_set_message (cpl_func, CPL_ERROR_INVALID_TYPE,
                                    "Column type is not supported");
  }

  CPLCHECK_MSG ("Cannot multiply");
  return CPL_ERROR_NONE;
}


/*---------------------------------------------------------------------------*/
/**
 * @brief Linear interpolation of matrix column
 *
 * @param matrix      The input matrix to interpolate (size n,m)
 * @param xref        The vector of current x axis (size m)
 * @param xout        The vector of target x axis (size k)
 * 
 * @return A newly allocated matrix of size (n,k)
 *
 * The routine interpolate each column of the matrix using
 * cpl_bivector_interpolate_linear
 */
/*---------------------------------------------------------------------------*/

cpl_matrix * gravi_matrix_interpolate_col (cpl_matrix * matrix,
                                           cpl_vector * xref,
                                           cpl_vector * xout)
{
    gravi_msg_function_start(0);
    
    cpl_ensure (matrix, CPL_ERROR_NULL_INPUT, NULL);
    cpl_ensure (xref,   CPL_ERROR_NULL_INPUT, NULL);
    cpl_ensure (xout,   CPL_ERROR_NULL_INPUT, NULL);

    cpl_size nrow = cpl_matrix_get_nrow (matrix);
    cpl_size ncol = cpl_matrix_get_ncol (matrix);
    cpl_size nxref = cpl_vector_get_size (xref);
    cpl_size nxout = cpl_vector_get_size (xout);

    cpl_ensure (ncol == nxref, CPL_ERROR_ILLEGAL_INPUT, NULL);
    cpl_ensure (nxout > 0,     CPL_ERROR_ILLEGAL_INPUT, NULL);

    /* Allocate memory */
    cpl_matrix * outmatrix = cpl_matrix_new (nrow, nxout);
    cpl_vector * yref = cpl_vector_new (nxref);
    cpl_vector * yout = cpl_vector_new (nxout);
    cpl_bivector * fref = cpl_bivector_wrap_vectors (xref, yref);
    cpl_bivector * fout = cpl_bivector_wrap_vectors (xout, yout);

    /* Loop on rows */
    for (cpl_size row = 0; row < nrow; row++) {
        for (cpl_size x = 0; x < nxref; x++)
            cpl_vector_set (yref, x, cpl_matrix_get (matrix, row, x));
        cpl_bivector_interpolate_linear (fout, fref);
        for (cpl_size x = 0; x < nxout; x++)
            cpl_matrix_set (outmatrix, row, x, cpl_vector_get (yout, x));
        CPLCHECK_NUL ("Cannot interpolate matrix");
    }

    FREE (cpl_bivector_unwrap_vectors, fref);
    FREE (cpl_bivector_unwrap_vectors, fout);
    FREE (cpl_vector_delete, yref);
    FREE (cpl_vector_delete, yout);
    
    gravi_msg_function_exit(0);
    return outmatrix;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief  Invers a matrix with singular value decomposition
 *
 * @param  a_in    The input matrix to invert
 * 
 * @return Inverse matrix. In case of error a NULL is returned.
 */
/*---------------------------------------------------------------------------*/

cpl_matrix * gravi_matrix_invertSV_create (cpl_matrix *a_in)
{
    gravi_msg_function_start(0);
    
	cpl_ensure (a_in, CPL_ERROR_NULL_INPUT, NULL);
	
	cpl_size m = cpl_matrix_get_nrow (a_in);
	cpl_size n = cpl_matrix_get_ncol (a_in);

	cpl_vector * w = cpl_vector_new (n);
	cpl_matrix * v = cpl_matrix_new (n, n);

	/*
	 * fill nr matrix with data
	 * w are the singular values of the matrix
	 */
	cpl_matrix * a_out = svdcmp (a_in, w, v);
	CPLCHECK_NUL ("Error in inverse SV");

	/* Check the singular values */
	for (cpl_size ii = 0; ii < n; ii++) {
	  if (cpl_vector_get (w, ii) < 0.1)
		cpl_msg_warning (cpl_func, "Singular Value %lld = %e",
                         ii, cpl_vector_get (w, ii));
	}

	/*
	 * Compute inverse (vt*1/w*at)
	 */
    cpl_matrix * a_inv = cpl_matrix_new (n,m);
    double * a_inv_data = cpl_matrix_get_data (a_inv);
    
	for (cpl_size j = 0; j < m; j++) {
		for (cpl_size i = 0; i < n; i++){
			double wv_at = 0;

			for (cpl_size ii = 0; ii < n; ii++)
				if (cpl_vector_get (w, ii) > 1e-15)
                    wv_at += cpl_matrix_get (v, i, ii) /
                        cpl_vector_get (w, ii) *
                        cpl_matrix_get (a_out, j, ii);

			a_inv_data[j + i * m] = wv_at;
		}
	}

    /* Delete intermediate matrix */
	cpl_vector_delete (w);
	cpl_matrix_delete (v);
	cpl_matrix_delete (a_out);

    gravi_msg_function_exit(0);
	return a_inv;
}

/*
 * pythag fonction from numerical recipe in C,
 * necessary for gravi_matrix_invertSV_create
 */
double pythag(double a, double b)
{
	double absa=fabs(a);
	double absb=fabs(b);
	if (absa > absb)
		return absa*sqrt(1.0+pow(absb/absa, 2));
	else
		return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+pow(absa/absb, 2)));
}

/*
 * svdcmp fonction from numerical recipe in C
 * w and v are pre-allocated results,
 * necessary for gravi_matrix_invertSV_create
 */
cpl_matrix * svdcmp (cpl_matrix * a_in, cpl_vector * w, cpl_matrix * v)
{
    gravi_msg_function_start(0);
	cpl_ensure (a_in, CPL_ERROR_NULL_INPUT, NULL);
	
	int flag, i, its, j, jj, k, l, nm, n, m;
	double anorm, c, f, g, h, s, scale, x, y, z;
	cpl_vector * rv1;
	cpl_matrix * a;

	a = cpl_matrix_duplicate(a_in);
	m = cpl_matrix_get_nrow (a_in);
	n = cpl_matrix_get_ncol (a_in);
	rv1 = cpl_vector_new(n);
    
	g = scale = anorm = 0.0;
	for (i = 0; i < n; i++) {
		l = i + 1;
		cpl_vector_set(rv1, i, scale*g);
		g = s = scale = 0.0;
		if (i < m) {
			for (k = i; k < m; k++)
				scale += fabs(cpl_matrix_get(a, k, i));
			if (scale) {

				for (k = i; k < m; k++) {

					cpl_matrix_set (a, k, i, cpl_matrix_get(a, k, i)/scale);
					s += cpl_matrix_get(a, k, i) * cpl_matrix_get(a, k, i);
				}
				f = cpl_matrix_get(a, i, i);
				g = -SIGN(sqrt(s),f);
				h = f * g - s;
				cpl_matrix_set (a, i, i, f - g);
				for (j = l; j < n; j++) {
					for (s = 0.0, k = i; k < m; k++)
						s += cpl_matrix_get(a, k, i) * cpl_matrix_get(a, k, j);

					f = s / h;
					for (k = i; k < m; k++)
						cpl_matrix_set (a, k, j, cpl_matrix_get(a, k, j) +
												 f * cpl_matrix_get(a, k, i));

				}
				for (k = i; k < m; k++)
					cpl_matrix_set (a, k, i, cpl_matrix_get(a, k, i) * scale);

			}
		}

		cpl_vector_set(w, i, scale *g);

		g = s = scale = 0.0;
		if (i < m && i != (n - 1)) {
			for (k = l; k < n; k++)
				scale += fabs(cpl_matrix_get(a, i, k));;
			if (scale) {
				for (k = l; k < n; k++) {
					cpl_matrix_set (a, i, k, cpl_matrix_get(a, i, k) / scale);
					s += cpl_matrix_get(a, i, k) * cpl_matrix_get(a, i, k);
				}
				f = cpl_matrix_get(a, i, l);
				g = -SIGN(sqrt(s),f);
				h = f * g - s;
				cpl_matrix_set (a, i, l, f - g);
				for (k = l; k < n; k++)
					cpl_vector_set(rv1, k, cpl_matrix_get(a, i, k)/h);
				for (j = l; j < m; j++) {
					for (s = 0.0, k = l; k < n; k++)
						s += cpl_matrix_get(a, j, k) * cpl_matrix_get(a, i, k);
					for (k = l; k < n; k++)
						cpl_matrix_set (a, j, k, cpl_matrix_get(a, j, k) + s * cpl_vector_get(rv1, k));
				}
				for (k = l; k < n; k++)
					cpl_matrix_set (a, i, k, cpl_matrix_get(a, i, k) * scale);
			}
		}

		anorm = fmax(anorm,(fabs(cpl_vector_get(w, i)) +
										fabs(cpl_vector_get(rv1, i))));

	}

	for (i = (n - 1); i >= 0; i--) {
		if (i < n) {
			if (g) {
				for (j = l; j < n; j++)
					cpl_matrix_set(v, j, i, (cpl_matrix_get(a, i, j) /
												cpl_matrix_get(a, i, l)) / g);

				for (j = l;j < n; j++) {
					for (s = 0.0, k = l; k < n; k++)
						s += cpl_matrix_get(a, i, k) * cpl_matrix_get(v, k, j);

					for (k = l; k < n; k++)
						cpl_matrix_set(v, k, j, cpl_matrix_get(v, k, j) + s *
													  cpl_matrix_get(v, k, i));

				}
			}
			for (j = l; j < n; j++) {
				cpl_matrix_set(v, i, j, 0.0);
				cpl_matrix_set(v, j, i, 0.0);
			}
		}
		cpl_matrix_set(v, i, i, 1.0);

		g = cpl_vector_get(rv1, i);
		l = i;
	}

	for (i = (IMIN(m,n) - 1); i >= 0; i--) {
		l = i + 1;
		g = cpl_vector_get(w, i);
		for (j = l; j < n; j++)
			cpl_matrix_set(a, i, j, 0.0);

		if (g) {
			g = 1.0/g;
			for (j = l; j < n; j++) {
				for (s = 0.0, k = l; k < m; k++)
					s += cpl_matrix_get(a, k, i) * cpl_matrix_get(a, k, j);

				f = (s / cpl_matrix_get(a, i, i)) * g;

				for (k = i; k < m; k++)
					cpl_matrix_set(a, k, j, cpl_matrix_get(a, k, j) + f *
														cpl_matrix_get(a, k, i));

			}
			for (j = i; j < m; j++)
				cpl_matrix_set(a, j, i, cpl_matrix_get(a, j, i) * g);

		}
		else
			for (j = i; j < m; j++)
				cpl_matrix_set(a, j, i, 0.0);

		cpl_matrix_set(a, i, i, cpl_matrix_get(a, i, i) + 1);

	}

	for (k = (n - 1); k >= 0; k--) {
		for (its = 1; its <= 60; its ++) {
			flag = 1;

			for (l = k; l >= 0; l--) {
				nm = l;
				if ((fabs(cpl_vector_get(rv1, l)) + anorm) == anorm) {
					flag = 0;
					break;
				}
				if ((fabs(cpl_vector_get(w, nm)) + anorm) == anorm)
					break;
			}

			if (flag) {
				c = 0.0;
				s = 1.0;
				for (i = l; i <= k; i++) {
					f = s * cpl_vector_get(rv1, i);
					cpl_vector_set(rv1, i, c * cpl_vector_get(rv1, i));
					if ((fabs(f) + anorm) == anorm)
						break;
					g = cpl_vector_get(w, i);
					h = pythag(f,g);
					cpl_vector_set(w, i, h);
					h = 1.0 / h;
					c = g * h;
					s = -f * h;
					for (j = 0; j < m; j++) {
						y = cpl_matrix_get(a, j, nm);
						z = cpl_matrix_get(a, j, i);
						cpl_matrix_set(a, j, nm, y*c+z*s);
						cpl_matrix_set(a, j, i, z*c-y*s);
					}
				}
			}
			z = cpl_vector_get(w, k);
			if (l == k) {
				if (z < 0.0) {
					cpl_vector_set(w, k, -z);
					for (j = 0; j < n; j++)
						cpl_matrix_set(v, j, k, -cpl_matrix_get(v, j, k));

				}
				break;
			}

			if (its == 120) {
				cpl_error_set_message(cpl_func,
						CPL_ERROR_ILLEGAL_INPUT,
						"no convergence in 120 svdcmp iterations");
				cpl_vector_delete(rv1);
				cpl_matrix_delete(a);
				return NULL;
			}
			x = cpl_vector_get(w, l);
			nm = k ;//- 1;
			y = cpl_vector_get(w, nm); ;
			g = cpl_vector_get(rv1, nm);
			h = cpl_vector_get(rv1, k);
			f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
			g = pythag(f,1.0);
			f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;
			c = s = 1.0;
			for (j = l; j < nm; j++) {
				i = j+1;
				g = cpl_vector_get(rv1, i);
				y = cpl_vector_get(w, i);
				h = s * g;
				g = c * g;
				z = pythag(f, h);
				cpl_vector_set(rv1, j, z);
				c = f / z;
				s = h / z;
				f = x * c + g * s;
				g = g * c - x * s;
				h = y * s;
				y *= c;
				for (jj = 0; jj < n; jj++) {
					x = cpl_matrix_get(v, jj, j);
					z=cpl_matrix_get(v, jj, i);
					cpl_matrix_set(v, jj, j, x*c+z*s);
					cpl_matrix_set(v, jj, i, z*c-x*s);
				}
				z = pythag(f, h);
				cpl_vector_set(w, j, z);

				if (z) {
					z = 1.0 / z;
					c = f * z;
					s = h * z;
				}
				f = c * g + s * y;
				x = c * y - s * g;
				for (jj = 0; jj < m; jj++) {
					y = cpl_matrix_get(a, jj, j);
					z = cpl_matrix_get(a, jj, i);
					cpl_matrix_set(a, jj, j, y*c+z*s);
					cpl_matrix_set(a, jj, i, z*c-y*s);
				}
			}
			cpl_vector_set(rv1, l, 0.0);
			cpl_vector_set(rv1, k, f);
			cpl_vector_set(w, k, x);
		}

	}
	cpl_vector_delete(rv1);
	
    gravi_msg_function_exit(0);
	return a;
}


/*----------------------------------------------------------------------------*/
/**
 * @brief   Extract rows from table based on the TIME column.
 * 
 * @param	table        input table
 * @param   start, end   max and min accepted time.
 *
 * @return A new table, with rows matching TIME>=start && TIME<end
 */
/*----------------------------------------------------------------------------*/
cpl_table * gravi_table_extract_time_interval (cpl_table *table, double start, double end)
{
    gravi_msg_function_start(1);
	cpl_ensure (table, CPL_ERROR_NULL_INPUT, NULL);
    
    /* Select only interested */
    cpl_table_select_all (table);
    cpl_table_and_selected_double (table, "TIME", CPL_NOT_LESS_THAN, start);
    cpl_table_and_selected_double (table, "TIME", CPL_LESS_THAN, end);
    cpl_table * out = cpl_table_extract_selected (table);
    
    gravi_msg_function_exit(1);
    return out;
}

/**@}*/
