/* $Id: gravi_cpl.h,v 1.12 2014/11/12 06:10:40 nazouaoui Exp $
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

#ifndef GRAVI_CPL_H_
#define GRAVI_CPL_H_

/*-----------------------------------------------------------------------------
                                   Includes
 -----------------------------------------------------------------------------*/

#include <cpl.h>

/*-----------------------------------------------------------------------------
                                       Prototypes
 -----------------------------------------------------------------------------*/

cpl_array ** gravi_array_new_list (int n, cpl_type type, int size);

double ** gravi_table_get_data_array_double (cpl_table * table, const char * name);
float ** gravi_table_get_data_array_float (cpl_table * table, const char * name);
float complex ** gravi_table_get_data_array_float_complex (cpl_table * table, const char * name);
double complex ** gravi_table_get_data_array_double_complex (cpl_table * table, const char * name);
int ** gravi_table_get_data_array_int (cpl_table * table, const char * name);

double gravi_table_get_column_mean (cpl_table * table, const char * name, int base, int nbase);
double gravi_table_get_column_std (cpl_table * table, const char * name, int base, int nbase);
cpl_array * gravi_table_get_column_mean_array (cpl_table * table, const char * name, int base, int nbase);

#define gravi_table_get_value(table,name,row,value) cpl_array_get (cpl_table_get_array (table,name,row), value, NULL)
#define gravi_table_set_value(table,name,row,value,val) cpl_array_set (cpl_table_get_data_array (table,name)[row], value, val);

int gravi_array_threshold_min (cpl_array * array, double lo_cut,
                               double assign_lo_cut);

cpl_array * gravi_array_init_double (long n , double value);
cpl_array * gravi_array_init_int (long n, int value);
cpl_array * gravi_array_init_double_complex (long n, double complex value);
cpl_array * gravi_array_init_float_complex (long n, float complex value);

cpl_array * gravi_array_wrap_complex (cpl_array * input_re, cpl_array * input_im);
cpl_array * gravi_array_wrap_float_complex (cpl_array * input_re, cpl_array * input_im);
cpl_array * gravi_array_compute_norm2 (cpl_array * input_re, cpl_array * input_im);
cpl_array * gravi_array_smooth (cpl_array * input, int nsmooth);
cpl_array * gravi_array_rebin (const cpl_array * input, const cpl_array * errs,
							   cpl_table * oi_wave_sc, cpl_table * oi_wave_ft);
cpl_array * gravi_array_create_inverse (cpl_array *input);
cpl_array * gravi_array_cexp (double complex factor, const cpl_array * input);

cpl_error_code gravi_array_normalize_complex (cpl_array * input);
cpl_error_code gravi_array_phase_unwrap (cpl_array * input);
cpl_error_code gravi_array_phase_wrap (cpl_array * input);
cpl_error_code gravi_array_multiply_phasor (cpl_array * input, double complex factor, cpl_array * phase);
cpl_error_code gravi_array_add_phasor (cpl_array * input, double complex factor, cpl_array * phase);
cpl_error_code gravi_array_add_phase (cpl_array * input, double factor, cpl_array * phase);
cpl_error_code gravi_array_multiply_conj (cpl_array * input1, cpl_array * input2);
cpl_error_code gravi_array_add_phasors (cpl_array * input, cpl_array * add, cpl_array * sub);

cpl_error_code gravi_table_set_array_double_complex (cpl_table * table, const char * name,
                                                     cpl_size row, cpl_array * visR, cpl_array * visI);
cpl_error_code gravi_table_set_array_phase (cpl_table * table, const char * name,
											cpl_size row, cpl_array * phase);
cpl_error_code gravi_table_set_string_fixlen (cpl_table *table, const char *name,
											  int row, const char *value, int len);
int gravi_table_are_equal (cpl_table * first, cpl_table * second);

cpl_error_code gravi_table_new_column (cpl_table * table, const char * name, const char * unit, cpl_type type);
cpl_error_code gravi_table_new_column_array (cpl_table * table, const char * name, const char * unit, cpl_type type, cpl_size size);
cpl_error_code gravi_table_init_column_array (cpl_table * table, const char * name, const char * unit, cpl_type type, cpl_size size);

cpl_error_code gravi_table_interpolate_column (cpl_table * to_table,
                                               const char * to_x,
                                               const char * to_y,
                                               const cpl_table * from_table,
                                               const char * from_x,
                                               const char * from_y);

cpl_error_code gravi_table_add_columns (cpl_table * oi_vis1, const char *name1,
										cpl_table * oi_vis2, const char *name2);
cpl_error_code gravi_table_runint_column (cpl_table * oi_vis, const char *input_name,
										  const char *output_name, int nsmooth, int nbase);
cpl_error_code gravi_table_smooth_column (cpl_table * oi_vis, const char *input_name,
										  const char *output_name, int nsmooth, int nbase);

cpl_error_code gravi_table_multiply_scalar (cpl_table * table, const char * name, int base, int nbase, double value);
cpl_error_code gravi_table_add_scalar (cpl_table * table, const char * name, int base, int nbase, double value);

cpl_array * gravi_table_create_sigma_array (cpl_table * oi_wave);
cpl_array * gravi_table_create_wave_array (cpl_table * oi_wave);

cpl_error_code gravi_array_get_group_delay_loop (cpl_array ** input, cpl_array * sigma,
						 double * gd, cpl_size nrow,
						 double max_width, int verbose);

cpl_error_code gravi_table_compute_group_delay (cpl_table * table, const char *input,
						const char *output, cpl_table * oi_wave);

cpl_matrix * get_matrix_from_vector (cpl_vector * , cpl_vector * );

cpl_imagelist * gravi_imagelist_from_column (cpl_table * table_data, const char * data_x);
cpl_imagelist * gravi_imagelist_wrap_column (cpl_table * table_data, const char * data_x);
cpl_error_code gravi_imagelist_unwrap_images (cpl_imagelist * imglist);
    
cpl_vector * gravi_table_get_vector (cpl_table * , cpl_size , const char * );
cpl_vector * gravi_table_get_vector_scalar (cpl_table * table,
                                            const char * name,
                                            cpl_size base,
                                            cpl_size nbase);

cpl_vector * gravi_table_get_vector_diff (cpl_table * spectrum_data, int index,
                                          const char * regname1,
                                          const char * regname2);

cpl_array * gravi_array_wrap_image (cpl_image *);

cpl_error_code gravi_image_subtract_window (cpl_image * img1, const cpl_image * img2,
                                            cpl_size llx, cpl_size lly,
                                            cpl_size urx, cpl_size ury,
                                            cpl_size llx2, cpl_size lly2);
cpl_error_code gravi_image_fill (cpl_image * img, double value);
cpl_image * gravi_image_wrap_matrix (cpl_matrix * matrix);
cpl_image * gravi_image_from_matrix (cpl_matrix * matrix);
cpl_image * gravi_image_wrap_vector (cpl_vector * vector);
cpl_image * gravi_image_from_vector (cpl_vector * vector);
cpl_image * gravi_image_from_column (cpl_table * table_data, const char * data_x, cpl_size row);

double gravi_image_get_quantile (const cpl_image * img, double value);
double gravi_array_get_quantile (cpl_array * arr, double value);

cpl_size gravi_vector_get_maxpos (cpl_vector * vector);
cpl_vector * gravi_vector_extract (const cpl_vector * vector, int start, int step);
cpl_error_code gravi_vector_unwrap_with_guess (cpl_vector * vector, cpl_vector * ref, double ref_to_phase);


cpl_array * gravi_table_get_column_dimension (const cpl_table * table, const char * name);

cpl_matrix * gravi_matrix_interpolate_col (cpl_matrix * matrix,
                                           cpl_vector * xref,
                                           cpl_vector * xout);

cpl_matrix * gravi_matrix_invertSV_create (cpl_matrix * a_in);


#endif /* GRAVI_CPL_H_ */


