/* $Id: gravi_data.h, 2011/04/07 11:22:28 nazouaoui $
 *
 * This file is part of the ESO Common Pipeline Library
 * Copyright (C) 2001-2008 European Southern Observatory
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

#ifndef GRAVI_DATA_H_
#define GRAVI_DATA_H_

CPL_BEGIN_DECLS

/*-----------------------------------------------------------------------------
                                   Includes
 ----------------------------------------------------------------------------*/

#include <cpl.h>
#include <stdio.h>
#include <string.h>

/*-----------------------------------------------------------------------------
                                   New types
 ----------------------------------------------------------------------------*/

typedef struct _gravi_data_ gravi_data;

/*-----------------------------------------------------------------------------
                                   Defines
 ----------------------------------------------------------------------------*/

#define gravi_data_get_oi_wave(data, type, pol, npol) gravi_data_get_oi_table (data, GRAVI_OI_WAVELENGTH_EXT, GRAVI_INSNAME(type,pol,npol))
#define gravi_data_get_oi_vis(data, type, pol, npol) gravi_data_get_oi_table (data, GRAVI_OI_VIS_EXT, GRAVI_INSNAME(type,pol,npol))
#define gravi_data_get_oi_vis2(data, type, pol, npol) gravi_data_get_oi_table (data, GRAVI_OI_VIS2_EXT, GRAVI_INSNAME(type,pol,npol))
#define gravi_data_get_oi_t3(data, type, pol, npol) gravi_data_get_oi_table (data, GRAVI_OI_T3_EXT, GRAVI_INSNAME(type,pol,npol))
#define gravi_data_get_oi_flux(data, type, pol, npol) gravi_data_get_oi_table (data, GRAVI_OI_FLUX_EXT, GRAVI_INSNAME(type,pol,npol))

#define gravi_data_get_wave_fibre(data, type) gravi_data_get_table (data, type==GRAVI_SC?GRAVI_WAVE_FIBRE_SC_EXT:GRAVI_WAVE_FIBRE_FT_EXT)
#define gravi_data_get_wave_fibre_plist(data, type) gravi_data_get_plist (data, type==GRAVI_SC?GRAVI_WAVE_FIBRE_SC_EXT:GRAVI_WAVE_FIBRE_FT_EXT)
#define gravi_data_get_wave_data(data, type) gravi_data_get_table (data, type==GRAVI_SC?GRAVI_WAVE_DATA_SC_EXT:GRAVI_WAVE_DATA_FT_EXT)
#define gravi_data_get_wave_data_plist(data, type) gravi_data_get_plist (data, type==GRAVI_SC?GRAVI_WAVE_DATA_SC_EXT:GRAVI_WAVE_DATA_FT_EXT)
#define gravi_data_has_wave(data, type) gravi_data_has_extension (data, type==GRAVI_SC?GRAVI_WAVE_DATA_SC_EXT:GRAVI_WAVE_DATA_FT_EXT)

#define gravi_data_get_p2vm_data(data, type) gravi_data_get_table (data, type==GRAVI_SC?GRAVI_P2VM_DATA_SC_EXT:GRAVI_P2VM_DATA_FT_EXT)
#define gravi_data_has_p2vm(data, type) gravi_data_has_extension (data, type==GRAVI_SC?GRAVI_P2VM_DATA_SC_EXT:GRAVI_P2VM_DATA_FT_EXT)

#define gravi_data_get_imaging_detector(data, type) gravi_data_get_table (data, type==GRAVI_SC?GRAVI_IMAGING_DETECTOR_SC_EXT:GRAVI_IMAGING_DETECTOR_FT_EXT)
#define gravi_data_has_detector(data, type) gravi_data_has_extension (data, type==GRAVI_SC?GRAVI_IMAGING_DETECTOR_SC_EXT:GRAVI_IMAGING_DETECTOR_FT_EXT)

#define gravi_data_get_spectrum_data(data, type) gravi_data_get_table (data, type==GRAVI_SC?GRAVI_SPECTRUM_DATA_SC_EXT:GRAVI_SPECTRUM_DATA_FT_EXT)
#define gravi_data_get_spectrum_data_plist(data, type) gravi_data_get_plist (data, type==GRAVI_SC?GRAVI_SPECTRUM_DATA_SC_EXT:GRAVI_SPECTRUM_DATA_FT_EXT)
#define gravi_data_has_spectrum(data, type) gravi_data_has_extension (data, type==GRAVI_SC?GRAVI_SPECTRUM_DATA_SC_EXT:GRAVI_SPECTRUM_DATA_FT_EXT)

#define gravi_data_get_profile_plist(data) gravi_data_get_plist (data, GRAVI_PROFILE_DATA_EXT)

#define gravi_data_get_oi_wave_plist(data, type, pol, npol) gravi_data_get_oi_plist (data, GRAVI_OI_WAVELENGTH_EXT, GRAVI_INSNAME(type,pol,npol))
#define gravi_data_get_oi_vis_plist(data, type, pol, npol) gravi_data_get_oi_plist (data, GRAVI_OI_VIS_EXT, GRAVI_INSNAME(type,pol,npol))
#define gravi_data_get_oi_vis2_plist(data, type, pol, npol) gravi_data_get_oi_plist (data, GRAVI_OI_VIS2_EXT, GRAVI_INSNAME(type,pol,npol))
#define gravi_data_get_oi_t3_plist(data, type, pol, npol) gravi_data_get_oi_plist (data, GRAVI_OI_T3_EXT, GRAVI_INSNAME(type,pol,npol))
#define gravi_data_get_oi_flux_plist(data, type, pol, npol) gravi_data_get_oi_plist (data, GRAVI_OI_FLUX_EXT, GRAVI_INSNAME(type,pol,npol))

#define gravi_data_get_header(data) gravi_data_get_plist (data, GRAVI_PRIMARY_HDR_EXT)
#define gravi_data_get_extname(data,ext) gravi_pfits_get_extname (gravi_data_get_plist_x(data, ext))
#define gravi_data_get_qc(data) gravi_plist_get_qc (gravi_data_get_plist(data, GRAVI_PRIMARY_HDR_EXT) )
#define gravi_data_is_internal(data) gravi_pfits_is_calib (gravi_data_get_plist(data, GRAVI_PRIMARY_HDR_EXT))
#define gravi_data_get_img(data,ext) cpl_imagelist_get (gravi_data_get_cube (data,ext), 0)
#define gravi_data_get_spec_res(data) gravi_pfits_get_spec_res (gravi_data_get_plist (data, GRAVI_PRIMARY_HDR_EXT))

/*----------------------------------------------------------------------------
                            Public prototypes
 ----------------------------------------------------------------------------*/

/* 
 * gravi data constructors and destructor
 */

gravi_data * gravi_data_new(int);
gravi_data * gravi_data_duplicate(const gravi_data *);
cpl_error_code gravi_data_append (gravi_data * first, const gravi_data * second, int force);
void gravi_data_delete(gravi_data *);

/* 
 * Load and save gravi data 
 */

gravi_data * gravi_data_load(const char * filename);
gravi_data * gravi_data_load_ext(const char * filename, 
                                 const char * extensions_regexp);
gravi_data * gravi_data_load_frame (cpl_frame * frame, cpl_frameset * used_frameset);
int gravi_data_patch (gravi_data * file_to_patch, cpl_frameset * patch_frameset);
gravi_data * gravi_data_load_rawframe (cpl_frame * frame, cpl_frameset * used_frameset);
gravi_data * gravi_data_load_rawframe_ext (cpl_frame * frame,
                                           cpl_frameset * used_frameset,
                                           const char * extensions_regexp);

cpl_error_code gravi_data_save_new (gravi_data 		  * self,
									cpl_frameset 	  * allframes,
									const char 		  * filename,
                                    const char        * suffix,
									const cpl_parameterlist * parlist,
									cpl_frameset	  * usedframes,
									cpl_frame * frame,
									const char 		  * recipe,
									cpl_propertylist  * applist,
									const char        * proCatg);

cpl_error_code gravi_data_save_data(gravi_data * ,
		                            const char * ,
		                            unsigned );

/*
 * Element access
 */

int gravi_data_has_extension (gravi_data * , const char * );
int gravi_data_has_type (gravi_data * self, const char * type);
int gravi_data_get_size (const gravi_data *);
int gravi_data_get_size_table (const gravi_data *);

cpl_propertylist * gravi_data_get_plist(gravi_data *, const char *);
cpl_propertylist * gravi_data_get_plist_x(gravi_data* , int );

cpl_table * gravi_data_get_table_x(gravi_data* , int );
cpl_imagelist * gravi_data_get_cube_x(gravi_data* , int );

cpl_table * gravi_data_get_table(gravi_data*, const char *);
cpl_imagelist * gravi_data_get_cube(gravi_data* , const char * );
cpl_propertylist * gravi_data_get_oi_plist(gravi_data * ,
		                     const char * , const char * );
cpl_table * gravi_data_get_oi_table(gravi_data * ,
		                     const char * , const char * );
cpl_table ** gravi_data_get_oiwave_tables (gravi_data * data, int type_data, int npol);

/*
 * Inserting and removing elements
 */

cpl_error_code gravi_data_erase_x (gravi_data *, int);
cpl_error_code gravi_data_erase (gravi_data *, const char *);
cpl_error_code gravi_data_erase_type (gravi_data * self, const char * type);

cpl_error_code gravi_data_add_table (gravi_data * self,
                                     cpl_propertylist * plist,
                                     const char * extname,
                                     cpl_table * table);

cpl_error_code gravi_data_add_cube (gravi_data * self,
                                    cpl_propertylist * plist,
                                    const char * extname,
                                    cpl_imagelist * imglist);
    
cpl_error_code gravi_data_add_img (gravi_data * self,
                                   cpl_propertylist * plist,
                                   const char * extname,
                                   cpl_image * image);

cpl_error_code gravi_data_copy_ext (gravi_data * output,
										 gravi_data * input,
										 const char * name);

cpl_error_code gravi_data_move_ext (gravi_data * output,
                                    gravi_data * input,
                                    const char * name);

/*
 * Work on gravi_data
 */

cpl_error_code gravi_data_clean_for_astro (gravi_data * data);
cpl_error_code gravi_data_check_consistency (gravi_data * data);
cpl_error_code gravi_data_detector_cleanup (gravi_data * data,
                                            const cpl_parameterlist * parlist);
cpl_error_code gravi_data_dump(gravi_data *self);
cpl_error_code gravi_data_dump_mode (gravi_data * data);


#endif

