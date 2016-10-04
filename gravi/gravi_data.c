/* $Id: gravi_data_3.c, 2011/05/18 15:21:53 nazouaoui Exp $
 *
 * This file is part of the ESO Common Pipeline Library
 * Copyright (C) 2001-2005 European Southern Observatory
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
 * @defgroup gravi_data Properties
 *
 * This module implements the gravi data type. The type @c gravi_data is
 * basically a structure of variables container which consists of type
 * cpl_propertylist and cpl_table. The contents of
 * the @c gravi_data structure is loaded from the FITS file so the property
 * lists will be stocked in the first variable @c propertylist, and the tables
 * associate at each property list will be stocked in the variable @c table.
 * The fields are similar to ordinary cpl_propertylist and cpl_table variables
 * so it can use all the features of these variables.
 */
/**@{*/

/*-----------------------------------------------------------------------------
                                   Includes
 -----------------------------------------------------------------------------*/

#include <cpl.h>
#include <stdio.h>
#include <string.h>
#include <config.h>
#include <ctype.h>
#include <time.h>
#include <math.h>

#include "gravi_dfs.h"
#include "gravi_data.h"
#include "gravi_pfits.h"
#include "gravi_utils.h"
#include "gravi_cpl.h"

/*-----------------------------------------------------------------------------
                                   Defines
 -----------------------------------------------------------------------------*/

#define GRAVI_DATA_SIZE 30

struct _gravi_data_{
	cpl_propertylist * primary_hdr;
	int nb_ext;
	cpl_propertylist ** exts_hdrs;
	cpl_imagelist ** exts_imgl;
	cpl_table ** exts_tbs;
};

/* ----------------------------------------------------------------------------
                            Private prototypes
   ---------------------------------------------------------------------------- */

int gravi_data_is_oi_ext (cpl_propertylist * hdr);
cpl_error_code gravi_data_check_savetypes(cpl_propertylist * hdr, cpl_table * oi_table);

/*-----------------------------------------------------------------------------
                             Functions code
 -----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
 * @brief Create an empty gravi_data.
 *
 * @return The newly created gravi_data type.
 *
 * The function allocates memory for a gravi_data.
 * The returned gravi_data must be deleted using the  destructor
 * @b gravi_data_delete(). The allocated structure can contain
 * up to GRAVI_DATA_SIZE extension, which with its own header
 * (propertylist). The extensions can be IMAGE or BINTABLE.
 */
/*----------------------------------------------------------------------------*/

gravi_data * gravi_data_new (int nb_ext)
{
    gravi_msg_function_start(0);
	cpl_ensure (nb_ext==0, CPL_ERROR_ILLEGAL_INPUT, NULL);
	
    gravi_data *self = cpl_malloc (sizeof (gravi_data));
    self->primary_hdr = cpl_propertylist_new ();
    self->nb_ext = 0;

    self->exts_hdrs = cpl_malloc (GRAVI_DATA_SIZE * sizeof(cpl_propertylist*));
    self->exts_tbs = cpl_malloc (GRAVI_DATA_SIZE * sizeof(cpl_table*));
    self->exts_imgl = cpl_malloc (GRAVI_DATA_SIZE * sizeof(cpl_imagelist*));

    for(int i = 0; i < GRAVI_DATA_SIZE; i++){
    	self->exts_tbs[i] = NULL;
    	self->exts_hdrs[i] = NULL;
    	self->exts_imgl[i] = NULL;
    }

    gravi_msg_function_exit(0);
    return self;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Delete a gravi data.

 * @param self  The gravi data to delete.
 *
 * The function destroys a gravi data of any kind. All the members
 * including their values are properly deallocated. If the gravi data @em self
 * is @c NULL, nothing is done and no error is set.
 */
/*----------------------------------------------------------------------------*/

void gravi_data_delete (gravi_data *self)
{
    gravi_msg_function_start(0);
	
    if (self) {
	  /* Delete main header */
	  FREE (cpl_propertylist_delete, self->primary_hdr);
	  /* Delete data */
	  FREELOOP (cpl_propertylist_delete, self->exts_hdrs, GRAVI_DATA_SIZE);
	  FREELOOP (cpl_table_delete, self->exts_tbs, GRAVI_DATA_SIZE);
	  FREELOOP (cpl_imagelist_delete, self->exts_imgl, GRAVI_DATA_SIZE);
	  /* Delete structure */
	  FREE (cpl_free, self);
	}
	
    gravi_msg_function_exit(0);
    return;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Check if EXTNAME starts with 'OI_'  (OIFITS extension)
 */
/*----------------------------------------------------------------------------*/

int gravi_data_is_oi_ext (cpl_propertylist * hdr)
{
  if (hdr==NULL) { return 0; }
  if (!cpl_propertylist_has (hdr, "EXTNAME")) { return 0; }
  if (strncmp (gravi_pfits_get_extname (hdr), "OI_", 3)) { return 0; }
  return 1;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Set the savetypes of the OIFITS table
 *
 * @param hdr:  property list of the table, to check it is an OIFITS table
 * @param oi_table: the OIFITS table to verify
 *
 * Ensure the FLAG columns are saved as BOOL, the STA_INDEX are saved as 
 * SHORT, the TARGET_ID are saved as SHORT and the MNT_STA are saved as SHORT.
 * This is to follow properly the OIFITS standard.
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_data_check_savetypes (cpl_propertylist * hdr, cpl_table * oi_table)
{
  cpl_ensure_code (hdr,      CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (oi_table, CPL_ERROR_NULL_INPUT);
  
  /* If not an OI_FITS table, nothing to do */
  if (!gravi_data_is_oi_ext (hdr)) return CPL_ERROR_NONE;

  if ( cpl_table_has_column (oi_table, "FLAG") ) {
	cpl_table_set_column_savetype(oi_table, "FLAG", CPL_TYPE_BOOL);
  }

  if ( cpl_table_has_column (oi_table, "STA_INDEX") ) {
	cpl_table_set_column_savetype(oi_table, "STA_INDEX", CPL_TYPE_SHORT);
  }
  
  if ( cpl_table_has_column (oi_table, "TARGET_ID") ) {
	cpl_table_set_column_savetype(oi_table, "TARGET_ID", CPL_TYPE_SHORT);
  }
  
  if ( cpl_table_has_column (oi_table, "MNT_STA") ) {
	cpl_table_set_column_savetype(oi_table, "MNT_STA", CPL_TYPE_SHORT);
  }

  /* Check errors */
  if ( cpl_error_get_code() != CPL_ERROR_NONE ) {
	cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
						  "Cannot reinstall the savetypes");
	return cpl_error_get_code();
  }

  return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Create a copy of the gravi data.
 *
 * @param self     The gravi data to duplicate.
 * @return The copy of self.
 *
 * The function returns a copy of the gravi data @em self. The copy is a
 * deep copy, i.e. all fields members are copied.  I return
 * @c NULL in case of an error. In the latter case an appropriate
 * error code is set.
 */
/*----------------------------------------------------------------------------*/

gravi_data * gravi_data_duplicate (const gravi_data *self)
{
    gravi_msg_function_start(1);
	cpl_ensure (self, CPL_ERROR_NULL_INPUT, NULL);
	
	gravi_data * copy = gravi_data_new(0);
    copy->nb_ext = self->nb_ext;
	
    if (self->primary_hdr) {
        cpl_propertylist_delete (copy->primary_hdr);
    	copy->primary_hdr = cpl_propertylist_duplicate(self->primary_hdr);
    } else {
        cpl_msg_warning (cpl_func,"a data without header shall not exist");
    }
	
	CPLCHECK_NUL ("Cannot duplicate header");

    for (int i = 0; i < copy->nb_ext ; i++){
    	if (self->exts_hdrs[i])
    		copy->exts_hdrs[i] = cpl_propertylist_duplicate(self->exts_hdrs[i]);
    	if (self->exts_tbs[i]) {
    		copy->exts_tbs[i] = cpl_table_duplicate(self->exts_tbs[i]);
			gravi_data_check_savetypes (copy->exts_hdrs[i], copy->exts_tbs[i]);
		}
    	if (self->exts_imgl[i])
    		copy->exts_imgl[i] = cpl_imagelist_duplicate(self->exts_imgl[i]);
		
		CPLCHECK_NUL ("Cannot duplicate extension");
    }

    gravi_msg_function_exit(1);
    return copy;
}


/*----------------------------------------------------------------------------*/
/**
 * @brief Append a gravi_data into another existing one.
 *
 * @param first: the main gravi_data (will grow)
 * @param second: the gravi_data to append to first
 *
 * The function only works for OIFITS so far.
 * The function skips some table that shall not be merged (OI_ARRAY,
 * OI_WAVELENGTH, OI_TARGET...). It appends the table extension and
 * the image extention by duplicating the data (@em second may be
 * properly desallocated).
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_data_append (gravi_data * first, const gravi_data * second, int force)
{
  gravi_msg_function_start(1);
  cpl_ensure_code (first,  CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (second, CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (first->nb_ext == second->nb_ext, CPL_ERROR_INCOMPATIBLE_INPUT);

  cpl_msg_warning (cpl_func, "Append data: keep only the HEADER of the first data");
  cpl_msg_warning (cpl_func, "Append data: FIXME: only works for OIFITS data");

  for (int i = 0; i < first->nb_ext ; i++) {

	/* Check EXTNAME */
	const char * name1 = gravi_pfits_get_extname (first->exts_hdrs[i]);
	const char * name2 = gravi_pfits_get_extname (second->exts_hdrs[i]);
	cpl_ensure_code (!strcmp (name1, name2), CPL_ERROR_INCOMPATIBLE_INPUT);
	
	/* These tables shall not be merged, but shall be equal */
	if ( !strcmp (name1,"OI_WAVELENGTH") ||
		 !strcmp (name1, "OI_TARGET") ||
		 !strcmp (name1, "OI_ARRAY") ||
		 !strcmp (name1, "IMAGING_DETECTOR_SC") ||
		 !strcmp (name1, "IMAGING_DETECTOR_FT") ||
		 !strcmp (name1, "ARRAY_DESCRIPTION") ||
		 !strcmp (name1, "ARRAY_GEOMETRY") ||
		 !strcmp (name1, "OPTICAL_TRAIN") ) {
	  if (force) {
		cpl_msg_info (cpl_func,"Don't check table %s", name1);
		continue;
	  } 
	  
	  cpl_msg_info (cpl_func,"Check table %s", name1);
	  if ( !gravi_table_are_equal (first->exts_tbs[i], second->exts_tbs[i]) ) {
		cpl_error_set_message (cpl_func, CPL_ERROR_ILLEGAL_INPUT, "Tables %s are different", name1);
		return 0;
	  }
	  continue;
	}

	/* Append as table */
	if ( first->exts_tbs[i] != NULL ) {
	  cpl_msg_info (cpl_func,"Merge table %s", name1);
	  cpl_table_insert (first->exts_tbs[i], second->exts_tbs[i], CPL_SIZE_MAX);

	  CPLCHECK_MSG ("Cannot append table");
	  continue;
	}

	/* Append as imagelist (duplicate each image) */
	if ( first->exts_imgl[i] != NULL ) {
	  cpl_msg_info (cpl_func,"Merge imglist %s", name1);
	  cpl_size n_first  = cpl_imagelist_get_size (first->exts_imgl[i]);
	  cpl_size n_second = cpl_imagelist_get_size (second->exts_imgl[i]);
	  CPLCHECK_MSG ("Cannot get the size of imglist");

	  for (cpl_size iimg = 0; iimg < n_second; iimg++) {
		cpl_image * img = cpl_image_duplicate (cpl_imagelist_get (second->exts_imgl[i], iimg));
		cpl_imagelist_set (first->exts_imgl[i], img, n_first + iimg);
	  }

	  CPLCHECK_MSG ("Cannot append imglist");
	  continue;
	}	
	
	CPLCHECK_MSG ("Cannot append data");
  }
  
  gravi_msg_function_exit(1);
  return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Low-level function to load FITS file.
 *
 * @param filename  Name of the input file.
 * @return The newly created gravi data or NULL if an error occurred.
 *
 * The function returns a gravi data created by reading the FITS file.
 * Currently only the FITS file format is supported. All the gravi data
 * memory members are allocated.
 */
/*----------------------------------------------------------------------------*/

gravi_data * gravi_data_load (const char * filename)
{
    gravi_msg_function_start(0);
	cpl_ensure (filename, CPL_ERROR_NULL_INPUT, NULL);

	cpl_msg_debug (cpl_func, "Load file : %s", filename);

	/* Find a number of extension on the FITS file */
	int nb_ext = cpl_fits_count_extensions (filename);
	if (nb_ext == -1){
		cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT, "no extension in this file");
		return NULL;
	}

	/* Create data */
	gravi_data * self = gravi_data_new (0);

	/* Load primary header */
    cpl_propertylist * header = cpl_propertylist_load (filename, 0);
	cpl_propertylist_append (self->primary_hdr, header);
    FREE (cpl_propertylist_delete, header);

	/* Loop on extensions */
	for (int i = 0; i < nb_ext ; i++ ){

	  /* Load header of this extension */
	  self->exts_hdrs[i] = cpl_propertylist_load (filename, i+1);
	  CPLCHECK_NUL ("Cannot load header");

	  /* Load extension as table */
	  if (gravi_pfits_get_extension_type (self->exts_hdrs[i]) == 2) {
		self->exts_tbs[i] = cpl_table_load (filename, i+1, 0);
		gravi_data_check_savetypes (self->exts_hdrs[i], self->exts_tbs[i]);
		CPLCHECK_NUL ("Cannot load bintable");
	  }
	  /* Load extension as imagelist */
	  else if (gravi_pfits_get_extension_type (self->exts_hdrs[i]) == 3) {
		self->exts_imgl[i] = cpl_imagelist_load (filename,CPL_TYPE_DOUBLE,i+1);
		CPLCHECK_NUL ("Cannot load imagelist");
	  }
	  /* error */
	  else {
		cpl_error_set_message (cpl_func, CPL_ERROR_ILLEGAL_INPUT, "The dimension of the extension is wrong");
		gravi_data_delete (self);
		return NULL;
	  }
	  
	  self->nb_ext ++;
	}

    gravi_msg_function_exit(0);
	return self;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Dump some information about data in messagin
 *
 * @param data       The gravi data to dump
 *
 * If the information is missing, the function does nothing
 * and return successfuly.
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_data_dump_mode (gravi_data * data)
{
    cpl_ensure_code (data, CPL_ERROR_NULL_INPUT);
	cpl_errorstate prestate = cpl_errorstate_get();
    
    /* Dump minimum info */
    cpl_propertylist * header = gravi_data_get_header (data);
    const char * res   = gravi_pfits_get_spec_res (header);
    const char * scpol = gravi_pfits_get_pola_mode (header, GRAVI_SC);
    const char * ftpol = gravi_pfits_get_pola_mode (header, GRAVI_FT);

    if (cpl_errorstate_is_equal(prestate)) {
        cpl_msg_info (cpl_func, "(insmode: %s %s %s)", res, scpol, ftpol);
    } else {
        cpl_errorstate_set (prestate);
    }
    
    return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Load a FITS file and create a gravi_data
 *
 * @param filename        Name of the input file.
 * @param used_frameset   If not NULL, this frameset is append with frame
 * @return gravi_data
 *
 * Load and create a gravi data type from input frame.
 */
/*----------------------------------------------------------------------------*/

gravi_data * gravi_data_load_frame (cpl_frame * frame,
									cpl_frameset * used_frameset)
{
  cpl_ensure (frame, CPL_ERROR_NULL_INPUT, NULL);
  
  const char * filename = cpl_frame_get_filename (frame);

  cpl_msg_info (cpl_func, "Load file %s (%s)", FILESHORT(filename), cpl_frame_get_tag (frame));
  
  if (used_frameset) cpl_frameset_insert (used_frameset, cpl_frame_duplicate (frame));

  /* Load data */
  gravi_data * data = gravi_data_load (filename);
  CPLCHECK_NUL ("Cannot load data");

  /* Dump minimum info */
  gravi_data_dump_mode (data);
  
  return data;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Load a RAW FITS file and create a gravi_data
 *
 * @param filename        Name of the input file.
 * @param used_frameset   If not NULL, this frameset is append with frame
 * @return gravi_data
 *
 * Load and create a gravi data type from input RAW frame.
 * Data integrity is verified.
 */
/*----------------------------------------------------------------------------*/

gravi_data * gravi_data_load_rawframe (cpl_frame * frame,
									   cpl_frameset * used_frameset)
{
  cpl_ensure (frame, CPL_ERROR_NULL_INPUT, NULL);
  
  const char * filename = cpl_frame_get_filename (frame);
  cpl_msg_info (cpl_func, "Load RAW file %s (%s)", FILESHORT(filename), cpl_frame_get_tag (frame));
  
  if (used_frameset) cpl_frameset_insert (used_frameset, cpl_frame_duplicate (frame));

  /* Load data */
  gravi_data * data = gravi_data_load (filename);
  CPLCHECK_NUL ("Cannot load data");

  /* Dump minimum info */
  gravi_data_dump_mode (data);
  
  /* Check consistency */
  gravi_data_check_consistency (data);
  CPLCHECK_NUL ("Cannot check data consistency");

  /* Delete unused data so far */
  gravi_data_erase (data, GRAVI_IMAGING_DATA_ACQ_EXT);
  gravi_data_erase (data, "ACQ_ABS_REF_POSITION");
  CPLCHECK_NUL ("Cannot erase useless data");

  return data;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Dump the overall structure of a gravi_data in stdout.
 * @param self : gravi_data to dump
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_data_dump (gravi_data *self)
{
  cpl_ensure_code (self, CPL_ERROR_NULL_INPUT);
  
  cpl_msg_info (cpl_func,"-------------------------");
  cpl_msg_info (cpl_func,"POINTER: %p", self);
  cpl_msg_info (cpl_func,"HEADER:  %p", self->primary_hdr);
  cpl_msg_info (cpl_func,"nb_ext = %i",self->nb_ext);
  
  for (int i=0; i<self->nb_ext; i++){
    cpl_msg_info (cpl_func,"%i: %s - %p %p %p", i,
		  gravi_pfits_get_extname (self->exts_hdrs[i]),
		  self->exts_hdrs[i],
		  self->exts_tbs[i],
		  self->exts_imgl[i]);
  }
  cpl_msg_info (cpl_func,"-------------------------");
	
  return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Get the number of extension in a gravi_data
 *
 * @param self     The input gravi_data.
 * @return The gravi data's current size, or 0 if the structure is empty.
 */
/*----------------------------------------------------------------------------*/

int gravi_data_get_size(const gravi_data *self)
{
    cpl_ensure (self, CPL_ERROR_ILLEGAL_INPUT, 0L);
    return (int)self->nb_ext;
}


/*----------------------------------------------------------------------------*/
/**
 * @brief Save a gravi data in a FITS file
 *
 * @param self      The gravi data to save or NULL if empty
 * @param filename  Name of the file to write
 * @param mode      The desired output options (combined with bitwise or)
 *
 * This function saves a gravi data to a FITS file, using cfitsio and save
 * all the tables members with its properties list.
 * Supported output modes are CPL_IO_CREATE (create a new file) and
 * CPL_IO_EXTEND  (append to an existing file).
 *
 * This is a low-level routine, wich does not provide a CPL compliant
 * product (see gravi_data_save_new).
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_data_save_data(gravi_data * self,
		                            const char * filename,
		                            unsigned mode)
{
    cpl_ensure_code (self,     CPL_ERROR_NULL_INPUT);
    cpl_ensure_code (filename, CPL_ERROR_NULL_INPUT);

	/* Create the file and save the first property list and table field with
	 * primary header entries. */
	if (self->exts_tbs[0] != NULL)
		cpl_table_save(self->exts_tbs[0], self->primary_hdr,
					   self->exts_hdrs[0], filename, mode);
	else if (self->exts_imgl[0] != NULL){
		cpl_propertylist_save (self->primary_hdr, filename, mode);
		cpl_imagelist_save (self->exts_imgl[0], filename,
				   CPL_TYPE_DOUBLE, self->exts_hdrs[0], CPL_IO_EXTEND);
	}
	else {
		cpl_error_set_message(cpl_func, CPL_ERROR_NULL_INPUT,
				                       "one of the inputs at least is NULL");
		return CPL_ERROR_NULL_INPUT;
	}


	/* Save the remainging extension */
	for (int i = 1; i < self->nb_ext; i++){
        
		if (gravi_pfits_get_extension_type (self->exts_hdrs[i]) == 2)
			cpl_table_save(self->exts_tbs[i], NULL,
		                      self->exts_hdrs[i], filename, CPL_IO_EXTEND);
		else if (gravi_pfits_get_extension_type (self->exts_hdrs[i]) == 3)
			cpl_imagelist_save (self->exts_imgl[i], filename,
					   CPL_TYPE_DOUBLE, self->exts_hdrs[i], CPL_IO_EXTEND);
		else {
			cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
					                "The dimension of the extension is wrong");
			gravi_data_delete(self);
			return CPL_ERROR_ILLEGAL_INPUT;
		}
	}

	return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
  * @brief Save a gravi data in a CPL-complian FITS file
  * 
  * @param self      The gravi data to save or NULL if empty
  * @param allframes  The list of input frames for the recipe
  * @param filename   Filename, or build from frame and proCatg if NULL
  * @param parlist    The list of input parameters
  * @param usedframes The list of raw/calibration frames used for this product
  * @param frame      The reference frame to build the header
  * @param recipe     The recipe name
  * @param applist    Optional propertylist to append to primary header or NULL
  * @param proCatg    Optional string coding the PRO.CATG
  */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_data_save_new (gravi_data 		  * self,
									cpl_frameset 	  * allframes,
									const char 		  * filename,
									const cpl_parameterlist * parlist,
									cpl_frameset	  * usedframes,
									cpl_frame * frame,
									const char 		  * recipe,
									cpl_propertylist  * applist,
									const char        * proCatg)
{
	gravi_msg_function_start(0);
	cpl_ensure_code (filename || frame, CPL_ERROR_NULL_INPUT);
	cpl_ensure_code (self,              CPL_ERROR_NULL_INPUT);
	cpl_ensure_code (proCatg,           CPL_ERROR_NULL_INPUT);
	
	cpl_frameset * frameset;
	int ext, i = 0, j = 0;

	/* If the optional propertylist is not given, we simply
	 * extract the QC parameters from the saved product. */
	if (applist == NULL) {
        applist = gravi_data_get_qc (self);
        
	} else {
        applist = cpl_propertylist_duplicate (applist);
    }

	/* Add the product CATG to the header */
	cpl_propertylist_append_string (applist, CPL_DFS_PRO_CATG, proCatg);

	/* Copy the DATE-OBS and NIGHT-OBS parameters if present */
	cpl_propertylist * hdr = gravi_data_get_header (self);
	if ( cpl_propertylist_has (hdr, GRAVI_NIGHT_OBS) ) {
	  cpl_propertylist_copy_property (applist, hdr, GRAVI_NIGHT_OBS);
	}
	if ( cpl_propertylist_has (hdr, "DATE-OBS") ) {
	  cpl_propertylist_copy_property (applist, hdr, "DATE-OBS");
	}

    /* Create keywords for OIFITS comliancy if VIS product */
    if (strstr (proCatg,"VIS") || strstr (proCatg,"TF")) {
        cpl_propertylist * tmp = gravi_plist_get_oifits_keywords (hdr);
        cpl_propertylist_append (applist, tmp);
        FREE (cpl_propertylist_delete, tmp);
    }

	/* Select the name extension depending on this catg */
	char catg_ext[800];
	for(i = 0; proCatg[i]; i++) if (proCatg[i]!='_') { catg_ext[j] = tolower(proCatg[i]); j++; }
	catg_ext[j] = '\0';


	/* If filename is NULL */
	if (filename == NULL && frame != NULL) {
	  filename = cpl_frame_get_filename (frame);
	} else if (filename == NULL && frame == NULL) {
	  cpl_error_set_message(cpl_func, CPL_ERROR_NULL_INPUT, "Need a frame if filename is void");
	  return CPL_ERROR_NULL_INPUT;
	}

	/* Use either the input filename or a name based on recipe */
	char * product_name = NULL;
	
	if ( cpl_parameterlist_find_const (parlist,"gravity.dfs.static-name") &&
		 gravi_param_get_bool (parlist,"gravity.dfs.static-name")) {
	  
	  /* Use the recipe name and add the selected catg_ext */
	  product_name = cpl_sprintf ("%s_%s.fits", recipe, catg_ext);
	}
	else {
	  
	  /* Remove the extension (last '.') and add the selected catg_ext */
	  char * filenoext = cpl_strdup (FILESHORT(filename));
	  char * lastdot = strrchr (filenoext, '.');
	  if (lastdot != NULL) *lastdot = '\0';
	  product_name = cpl_sprintf ("%s_%s.fits", filenoext, catg_ext);
	  FREE (cpl_free, filenoext);
	}
	
	cpl_msg_info (cpl_func, "Save file to %s", product_name);

	if (usedframes && frame==NULL) {
	  /* frameset is only the usedframes */
	  frameset = cpl_frameset_duplicate (usedframes);
	}
	else if (usedframes && frame) {
	  /* frameset is the usedframes, and
	   * frame is ensured to be the first */
	  frameset = cpl_frameset_new ();
	  cpl_frameset_insert (frameset, cpl_frame_duplicate (frame));
	  for (int f = 0; f < cpl_frameset_get_size (usedframes); f++) 
		if (strcmp (cpl_frame_get_filename (cpl_frameset_get_position (usedframes, f)), cpl_frame_get_filename (frame)))
		  cpl_frameset_insert (frameset, cpl_frame_duplicate (cpl_frameset_get_position (usedframes, f)));
	}
	else if (usedframes==NULL && frame) {
	  /* frameset is only the frame */
	  frameset = cpl_frameset_new ();
	  cpl_frameset_insert (frameset, cpl_frame_duplicate (frame));
	}
	else {
	  /* This shall not happen */
	  cpl_error_set_message(cpl_func, CPL_ERROR_NULL_INPUT, "Bug: need a frame or a usedframes");
	  return CPL_ERROR_NULL_INPUT;
	}

	/* Save the primary header */
	cpl_dfs_save_propertylist (allframes, NULL, parlist,
							   frameset, frame, recipe, applist,
							   NULL, PACKAGE_STRING, product_name);
	CPLCHECK_MSG("Cannot save the first extension primary header");

	/* Save the extensions */
	for (ext = 0; ext < self->nb_ext; ext ++)
	{
		if (self->exts_tbs[ext] != NULL)
			cpl_table_save (self->exts_tbs[ext], NULL,
							self->exts_hdrs[ext], product_name, CPL_IO_EXTEND);
		else if (self->exts_imgl[ext] != NULL)
			cpl_imagelist_save (self->exts_imgl[ext], product_name,
								cpl_image_get_type (cpl_imagelist_get (self->exts_imgl[ext], 0)),
								self->exts_hdrs[ext], CPL_IO_EXTEND);
		CPLCHECK_MSG("Cannot save the extension");
	}
	
	FREE (cpl_frameset_delete, frameset);
	FREE (cpl_free, product_name);
    FREE (cpl_propertylist_delete, applist);

	gravi_msg_function_exit(0);
	return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Perform self-bias correction to the SC raw data
 * 
 * @param data     gravi_data to be processed (in-place).
 * 
 * This function applies the self bias-pixel correction to the IMAGING_DATA_SC
 * images, in-place. The bias is estimated independently for each image in the
 * imagelist (DIT). It is the median of the bias-pixel, whose location depend
 * on the spectral setup (HIGH has bias-column, while MED and LOW have
 * bias-lines). The functions add QC paramter to recover the mean subtracted
 * bias.
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_data_detector_cleanup (gravi_data * data,
                                            const cpl_parameterlist * parlist)
{
  gravi_msg_function_start(1);
  cpl_ensure_code (data, CPL_ERROR_NULL_INPUT);

  int nv = 0;

  /* Read header and number of regions */
  cpl_propertylist * header = gravi_data_get_header (data);
  cpl_size nreg = cpl_table_get_nrow (gravi_data_get_table (data, GRAVI_IMAGING_DETECTOR_SC_EXT));
  const char * resolution = gravi_pfits_get_resolution (header);

  CPLCHECK_MSG ("Cannot get data");
  
  /* Get the data as imagelist */
  cpl_imagelist * imglist = gravi_data_get_cube (data, GRAVI_IMAGING_DATA_SC_EXT);
  cpl_size nframe = cpl_imagelist_get_size (imglist);
  cpl_size nx = cpl_image_get_size_x (cpl_imagelist_get (imglist, 0));
  cpl_size ny = cpl_image_get_size_y (cpl_imagelist_get (imglist, 0));

  CPLCHECK_MSG ("Cannot get data");
  
  /* To save the list of bias correction */
  cpl_vector * bias_list = cpl_vector_new (nframe);

  /* Some hardcoded information */
  cpl_size nx_bias_hr = 4;
  cpl_size ny_reg_mr = (ny - 1) / nreg;

  // cpl_msg_info (cpl_func, "nreg=%lli, nx=%lli, ny=%lli, ny_reg_mr=%lli", nreg, nx, ny, ny_reg_mr);

  if ( !strcmp(resolution, "HIGH") ) {
   /* High Resolution - the first 4 columns
    * are for bias. FIXME: make sure this is also true in
    *  in HIGH-COMBINED */
	
	cpl_vector * bias = cpl_vector_new (nx_bias_hr * ny);
	
	/* Loop on frames */
	for (cpl_size f = 0; f < nframe; f++) {
	  cpl_image * frame = cpl_imagelist_get (imglist, f);

	  /* Get the bias pixels */
	  for (cpl_size x = 0; x < nx_bias_hr; x++)
		for (cpl_size y = 0; y < ny; y++) {
		  cpl_vector_set (bias, x * ny + y, cpl_image_get (frame, x+1, y+1, &nv));
		  CPLCHECK_MSG ("Cannot get the bias pixels");
		}

	  /* Remove the median of bias pixels to image */
	  double bias_med = cpl_vector_get_median (bias);
	  cpl_vector_set (bias_list, f, bias_med);
	  cpl_image_subtract_scalar (frame, bias_med);
	  CPLCHECK_MSG ("Cannot subtract bias");
	}
    FREE (cpl_vector_delete, bias);
  }
  else if ( !strcmp (gravi_param_get_string_default (parlist,
                     "gravity.preproc.bias-method","MEDIAN"),
                     "MEDIAN_PER_COLUMN")) {
    /* Low and Medium - the first line of
     * each region is bias. Use the median of the bias-pixel
     * per column !!*/

    gravi_msg_fixme ("Remove bias per column -- experimental");
		
	cpl_vector * bias = cpl_vector_new (nx);
	cpl_vector * bias_column = cpl_vector_new (nreg);

	/* Loop on frames */
	for (cpl_size f = 0; f < nframe; f++) {
	  cpl_image * frame = cpl_imagelist_get (imglist, f);

      /* Loop on columns */
	  for (cpl_size x = 0; x < nx; x++)  {
          
          /* Get the bias pixels of this column */
          for (cpl_size y = 0; y < nreg; y++) {
              double value = cpl_image_get (frame, x+1, (y+1)*ny_reg_mr+1, &nv);
              cpl_vector_set (bias_column, y, value);
              CPLCHECK_MSG ("Cannot get the bias pixels");
          }

          /* Compute the bias of this column, and save it */
          double bias_med = cpl_vector_get_median (bias_column);
          cpl_vector_set (bias, x, bias_med);

          /* Remove the bias from this column */
          for (cpl_size y = 0; y < ny; y++) {
              double value = cpl_image_get (frame, x+1, y+1, &nv);
              cpl_image_set (frame, x+1, y+1, value - bias_med);
              CPLCHECK_MSG ("Cannot set the bias corrected pixels");
          }
      } /* End loop on columns */
	  
	  /* Remove the median of bias pixels to image */
	  double bias_mean = cpl_vector_get_mean (bias);
	  cpl_vector_set (bias_list, f, bias_mean);
	}
    FREE (cpl_vector_delete, bias);
    FREE (cpl_vector_delete, bias_column);
  }
  else {
    /* Low and Medium - the first line of
     * each region is bias */
		
	cpl_vector * bias = cpl_vector_new (nx * nreg);

	/* Loop on frames */
	for (cpl_size f = 0; f < nframe; f++) {
	  cpl_image * frame = cpl_imagelist_get (imglist, f);

	  /* Get the bias pixels */
	  for (cpl_size x = 0; x < nx; x++) 
	  for (cpl_size y = 0; y < nreg; y++) {
		  cpl_vector_set (bias, x * nreg + y, cpl_image_get (frame, x+1, (y+1)*ny_reg_mr+1, &nv));
		  CPLCHECK_MSG ("Cannot get the bias pixels");
	  } 
	  
	  /* Remove the median of bias pixels to image */
	  double bias_med = cpl_vector_get_median (bias);
	  cpl_vector_set (bias_list, f, bias_med);
	  cpl_image_subtract_scalar (frame, bias_med);
	  CPLCHECK_MSG ("Cannot subtract bias");
	}
    FREE (cpl_vector_delete, bias);
  }

  /* Add some QC */
  cpl_propertylist_update_double (header, "ESO QC PIXBIAS AVG", cpl_vector_get_mean (bias_list));
  cpl_propertylist_set_comment (header, "ESO QC PIXBIAS AVG", "[adu] avg over the bias pixel");
  
  cpl_propertylist_update_double (header, "ESO QC PIXBIAS RMS", cpl_vector_get_stdev (bias_list));
  cpl_propertylist_set_comment (header, "ESO QC PIXBIAS RMS", "[adu] rms over frames");
  cpl_vector_delete (bias_list);
  
  gravi_msg_function_exit(1);
  return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/

/* Return the first extension exactly matching EXTNAME=name */
inline static int _gravi_data_find (const gravi_data *self, const char *name)
{
    int ext = 0;
    while (ext < self->nb_ext){
           cpl_propertylist *p = self->exts_hdrs[ext];
           const char *key = cpl_propertylist_get_string (p, "EXTNAME");
           if (strcmp(key, name) == 0) break;
           ext ++;
    }
    return ext;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Copy extensions from one data to another.
 *
 * @param output:  gravi_data to insert the extension
 * @param input:   gravi_data to read the extension
 * @param name:    EXTNAME of the extension to copy
 * 
 * Deep copy of the extension NAME from input to output
 * Copy all if several extension with same EXTNAME are found.
 * Silent if the copy cannot be done.
 */
/*---------------------------------------------------------------------------*/

cpl_error_code gravi_data_copy_ext (gravi_data * output,
                                    gravi_data * input,
                                    const char * name)
{
  cpl_ensure_code (output, CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (input,  CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (name,   CPL_ERROR_NULL_INPUT);

  int n = gravi_data_get_size (input);
  for (int j = 0; j < n; j++){
	
	cpl_propertylist * plist = gravi_data_get_plist_x (input, j);
	const char * plist_name = gravi_pfits_get_extname (plist);

	CPLCHECK_MSG ("Cannot get input data");
	
	if (plist_name == NULL) continue;
	
	if (!(strcmp (plist_name, name))) {
	  int type_data = gravi_pfits_get_extension_type (plist);
	  
	  if (type_data == 2) {
          gravi_data_add_table (output, cpl_propertylist_duplicate (plist), NULL, 
                                cpl_table_duplicate (gravi_data_get_table_x (input, j)));
	  }
	  else if (type_data == 3) {
          gravi_data_add_cube (output, cpl_propertylist_duplicate (plist), NULL,
                               cpl_imagelist_duplicate (gravi_data_get_cube_x (input, j)));
      }

	  CPLCHECK_MSG ("Cannot copy extension");
	}
  }

  return CPL_ERROR_NONE;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Move extensions from one data to another
 *
 * @param output:  gravi_data to insert the extension
 * @param input:   gravi_data to remove the extension
 * @param name:    EXTNAME of the extension to move
 * 
 * Move the extension NAME from input to output
 * Move all if several extension with same EXTNAME are found.
 * Silent if no move could be done.
 */
/*---------------------------------------------------------------------------*/

cpl_error_code gravi_data_move_ext (gravi_data * output,
                                    gravi_data * input,
                                    const char * name)
{
    gravi_msg_function_start(0);
    cpl_ensure_code (output, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code (input,  CPL_ERROR_NULL_INPUT);
    cpl_ensure_code (name,   CPL_ERROR_NULL_INPUT);
    
    int n = gravi_data_get_size (input);
    
    int jnew = 0;
    for (int j = 0; j < n; j++){
        
        cpl_propertylist * plist = gravi_data_get_plist_x (input, j);
        const char * plist_name = gravi_pfits_get_extname (plist);
        CPLCHECK_MSG ("Cannot get input");
        
        /* Case we need to move this extension */
        if ( (plist_name!=NULL) &&
             !(strcmp (plist_name, name))) {

            cpl_msg_info (cpl_func, "Move extension %s", plist_name);
            int type_data = gravi_pfits_get_extension_type (plist);
            if (type_data == 2) {
                gravi_data_add_table (output, plist, NULL, 
                                      gravi_data_get_table_x (input, j));
            }
            else if (type_data == 3) {
                gravi_data_add_cube (output, plist, NULL,
                                     gravi_data_get_cube_x (input, j));
            }
            CPLCHECK_MSG ("Cannot move extension");
            
        } else {
            
            /* Move pointer in the input */
            input->exts_hdrs[jnew] = input->exts_hdrs[j];
            input->exts_imgl[jnew] = input->exts_imgl[j];
            input->exts_tbs[jnew]  = input->exts_tbs[j];
            jnew ++;
        }
        
    } /* End loop on extension */
    
    /* Cleanup input further away */
    input->nb_ext = jnew;
    for (int j = input->nb_ext; j < GRAVI_DATA_SIZE; j++) {
        input->exts_hdrs[j] = NULL;
        input->exts_imgl[j] = NULL;
        input->exts_tbs[j]  = NULL;
    }

    gravi_msg_function_exit(0);
    return CPL_ERROR_NONE;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Check if data has extension with given EXTNAME
 *
 * @param raw_calib:  gravi_data to search in
 * @param ext_name:   name of the extension to check
 * @return 1 if extension exists, 0 otherwise.
 */
/*---------------------------------------------------------------------------*/

int gravi_data_has_extension (gravi_data * raw_calib, const char * ext_name)
{
    cpl_ensure (raw_calib, CPL_ERROR_NULL_INPUT, 0);
    cpl_ensure (ext_name,  CPL_ERROR_NULL_INPUT, 0);
	
	int test = 1;
	int ext = _gravi_data_find (raw_calib, ext_name);

	if (ext == raw_calib->nb_ext) {
		test = 0;
	}

	return test;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Return the number of ext whose EXTNAME and INSNAME match 'type'
 *
 * @param self         gravi_data to search in
 * @param type         name to search (ex: "_SC", "_FT", "METROLOGY")
 * @return The number of extension matching
 */
/*---------------------------------------------------------------------------*/

int gravi_data_has_type (gravi_data * self, const char * type)
{
    cpl_ensure (self, CPL_ERROR_NULL_INPUT, -1);
    cpl_ensure (type, CPL_ERROR_NULL_INPUT, -1);
    cpl_ensure (strlen (type)>1, CPL_ERROR_ILLEGAL_INPUT, -1);

    int counter = 0;

    /* Loop on extension */
    for (int ext = 0; ext < gravi_data_get_size (self) ; ext ++) {
	
        cpl_propertylist * plist = gravi_data_get_plist_x (self, ext);

        /* Check if EXTNAME or INSNAME contains the 'type' 
         * Warning that the number of extension is changed in the loop */
        if ( (cpl_propertylist_has (plist, "INSNAME") &&
              strstr (gravi_pfits_get_insname (plist), type) ) ||
             (cpl_propertylist_has (plist, "EXTNAME") &&
              strstr (gravi_pfits_get_extname (plist), type) ) )
        {
            cpl_msg_debug (cpl_func,"Find '%s' ", gravi_pfits_get_extname (plist));
            counter ++;
        }
    }
    return counter;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Get the propertylist of an extension by position.
 *
 * @param self      The gravi data to search in.
 * @param i         The selected extension
 * 
 * @return The requested property list, or @c NULL on error.
 *
 * The function returns a pointer to the property list corresponding to the 
 * extension i.
 * The returned propertylist is still owned by @em self, i.e. it
 * must not be deleted through the returned pointer.
 */
/*---------------------------------------------------------------------------*/

cpl_propertylist * gravi_data_get_plist_x (gravi_data* self, int i)
{
    cpl_ensure (self,           CPL_ERROR_NULL_INPUT, NULL);
    cpl_ensure (i>=0,           CPL_ERROR_ILLEGAL_INPUT, NULL);
    cpl_ensure (i<self->nb_ext, CPL_ERROR_ILLEGAL_INPUT, NULL);

	return self->exts_hdrs[i];
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Get the table of an extension by position.
 *
 * @param self      The gravi data to search in.
 * @param i         The selected extension
 * 
 * @return The requested table, or @c NULL on error.
 *
 * The function returns a pointer to the table corresponding to the 
 * extension i, or NULL on error or if the extension is not
 * a table. The returned table is still owned by the gravi_data set,
 * i.e. it must not be deleted through the returned pointer.
 */
/*---------------------------------------------------------------------------*/

cpl_table * gravi_data_get_table_x (gravi_data* self, int i)
{
    cpl_ensure (self,           CPL_ERROR_NULL_INPUT, NULL);
    cpl_ensure (i>=0,           CPL_ERROR_ILLEGAL_INPUT, NULL);
    cpl_ensure (i<self->nb_ext, CPL_ERROR_ILLEGAL_INPUT, NULL);
	
	return self->exts_tbs[i];
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Get the cube of an extension by position.
 *
 * @param self      The gravi data to search in.
 * @param i         The selected extension
 * 
 * @return The requested imglist, or @c NULL on error.
 *
 * The function returns a pointer to the imglist corresponding to the 
 * extension i, or NULL on error or if the extension is not
 * a imaglist. The returned imglist is still owned by the gravi_data,
 * i.e. it must not be deleted through the returned pointer.
 */
/*---------------------------------------------------------------------------*/

cpl_imagelist * gravi_data_get_cube_x (gravi_data* self, int i)
{
    cpl_ensure (self,           CPL_ERROR_NULL_INPUT, NULL);
    cpl_ensure (i>=0,           CPL_ERROR_ILLEGAL_INPUT, NULL);
    cpl_ensure (i<self->nb_ext, CPL_ERROR_ILLEGAL_INPUT, NULL);
	
	return self->exts_imgl[i];
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Get an OI_FITS table from EXTNAME and INSNAME
 *
 * @param self      The gravi data to search in.
 * @param extname   The requested EXTNAME
 * @param insname   The requested INSNAME
 * 
 * @return The requested table, or @c NULL on error.
 *
 * The function returns a pointer to the table corresponding to the 
 * extension with the requested INSNAME and EXTNAME.
 * The returned table is still owned by the gravi_data,
 * i.e. it must not be deleted through the returned pointer.
 */
/*---------------------------------------------------------------------------*/

cpl_table * gravi_data_get_oi_table (gravi_data * self,
                                     const char * extname,
                                     const char * insname)
{
    cpl_ensure (self,    CPL_ERROR_NULL_INPUT, NULL);
    cpl_ensure (extname, CPL_ERROR_NULL_INPUT, NULL);
    cpl_ensure (insname, CPL_ERROR_NULL_INPUT, NULL);
	
	int nb_ext = self->nb_ext;
	int pos = -1;
	for (int i = 0; i < nb_ext; i++){

		if (cpl_propertylist_has (self->exts_hdrs[i], "INSNAME")){
		if ((!strcmp (gravi_pfits_get_insname (self->exts_hdrs[i]), insname)) &&
				(!strcmp (gravi_pfits_get_extname (self->exts_hdrs[i]), extname))){
			pos = i;
			break;
		}
		}
	}

    if (pos < 0){
       	cpl_msg_warning(cpl_func, "The extention %s doesn't exist "
       			"with the insname %s",
       			extname, insname);
      	return NULL;
    }

    return self->exts_tbs[pos];
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Get the propertylist from EXTNAME and INSNAME
 *
 * @param self      The gravi data to search in.
 * @param extname   The requested EXTNAME
 * @param insname   The requested INSNAME
 * 
 * @return The requested propertylist, or @c NULL on error.
 *
 * The function returns a pointer to the first propertylist
 * corresponding to the  extension with the requested
 * INSNAME and EXTNAME.
 * The returned pointer is still owned by the gravi_data,
 * i.e. it must not be deleted.
 */
/*---------------------------------------------------------------------------*/

cpl_propertylist * gravi_data_get_oi_plist (gravi_data * self,
                                            const char * extname,
                                            const char * insname)
{
    cpl_ensure (self,    CPL_ERROR_NULL_INPUT, NULL);
    cpl_ensure (extname, CPL_ERROR_NULL_INPUT, NULL);
    cpl_ensure (insname, CPL_ERROR_NULL_INPUT, NULL);

	int nb_ext = self->nb_ext;
	int pos = -1;
	for (int i = 0; i < nb_ext; i++){
	  if (cpl_propertylist_has (self->exts_hdrs[i], "INSNAME")){
		if ((!strcmp (gravi_pfits_get_insname (self->exts_hdrs[i]), insname)) &&
			(!strcmp (gravi_pfits_get_extname (self->exts_hdrs[i]), extname))){
		  pos = i;
		  break;
		}
	  }
	}

    if (pos < -1){
       	cpl_msg_warning(cpl_func, "The extention %s doesn't exist "
       			"with the insname %s",
       			extname, insname);
      	return NULL;
    }

    return self->exts_hdrs[pos];
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Get the propertylist from EXTNAME
 *
 * @param self      The gravi data to search in.
 * @param extname   The requested EXTNAME
 * 
 * @return The requested propertylist, or @c NULL on error.
 *
 * The function returns a pointer to the first propertylist
 * corresponding to the extension with the requested
 * EXTNAME. The propertylist of the main HEADER can be requested
 * with extname = GRAVI_PRIMARY_HDR_EXT
 * The returned pointer is still owned by the gravi_data,
 * i.e. it must not be deleted.
 */
/*---------------------------------------------------------------------------*/

cpl_propertylist * gravi_data_get_plist (gravi_data * self,
                                                const char * extname)
{
    cpl_ensure (self,    CPL_ERROR_NULL_INPUT, NULL);
    cpl_ensure (extname, CPL_ERROR_NULL_INPUT, NULL);

	if (!strcmp(extname, GRAVI_PRIMARY_HDR_EXT))
		return self->primary_hdr;

	int pos = _gravi_data_find(self, extname);

	cpl_ensure (pos<self->nb_ext, CPL_ERROR_ILLEGAL_INPUT, NULL);
	
    return self->exts_hdrs[pos];
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Return a pointer on a table extension by its EXTNAME.
 *
 * @param self     The gravi data to search in.
 * @param name     The EXTNAME to search
 *
 * @return A pointer to the table, or NULL.
 *
 * The function returns a pointer to the table with name
 * @em name, or @c NULL if it does not exist or is not a bintable.
 * @note The returned table is still owned by the gravi data set, i.e. the
 * obtained table must not be deleted through the returned handle.
 */
/*---------------------------------------------------------------------------*/

cpl_table * gravi_data_get_table (gravi_data* self, const char * extname)
{
    cpl_ensure (self,    CPL_ERROR_NULL_INPUT, NULL);
    cpl_ensure (extname, CPL_ERROR_NULL_INPUT, NULL);

	int pos = _gravi_data_find(self, extname);

	cpl_ensure (pos<self->nb_ext, CPL_ERROR_ILLEGAL_INPUT, NULL);
	
    if ((self->exts_tbs[pos] == NULL) ||
    		(gravi_pfits_get_extension_type (self->exts_hdrs[pos]) != 2)){
       	cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
       					"No BINTABLE at the extension %s", extname);
      	return NULL;
    }

    return self->exts_tbs[pos];
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Return a pointer on an IMAGE extension by its EXTNAME.
 *
 * @param self     The gravi data to search in.
 * @param name     The EXTNAME to search
 *
 * @return A pointer to the imagelist, or NULL.
 *
 * The function returns a pointer to the imagelist with name
 * @em name, or @c NULL if it does not exist or is not a IMAGE.
 * @note The returned imagelist is still owned by the gravi data set,
 * i.e. the obtained table must not be deleted through the returned pointer.
 */
/*---------------------------------------------------------------------------*/

cpl_imagelist * gravi_data_get_cube (gravi_data* self, const char * extname)
{
    cpl_ensure (self,    CPL_ERROR_NULL_INPUT, NULL);
    cpl_ensure (extname, CPL_ERROR_NULL_INPUT, NULL);

	int pos = _gravi_data_find(self, extname);

	cpl_ensure (pos<self->nb_ext, CPL_ERROR_ILLEGAL_INPUT, NULL);
	
    if ((self->exts_imgl[pos] == NULL) ||
    		(gravi_pfits_get_extension_type (self->exts_hdrs[pos]) != 3)){
       	cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
       					"No IMAGE at the extension %s", extname);
      	return NULL;
    }

    return self->exts_imgl[pos];
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Erase an extension by its position
 * 
 * @param self       The gravi_data to manipulate in-place
 * @param pos        The extension to delete
 * 
 * The function erases the extension in position @c pos.
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_data_erase_x (gravi_data * self, int pos)
{
    cpl_ensure_code (self,   CPL_ERROR_NULL_INPUT);
    cpl_ensure_code (pos>-1, CPL_ERROR_ILLEGAL_INPUT);

	/* Delete data or move pointers */
	int comp = 0;
	for (int i = 0; i < self->nb_ext; i++) {
		if (i != pos){
			self->exts_hdrs[comp] = self->exts_hdrs[i];
			self->exts_tbs[comp] = self->exts_tbs[i];
			self->exts_imgl[comp] = self->exts_imgl[i];
			comp++;
		}
		else{
			if ((gravi_pfits_get_extension_type (self->exts_hdrs[pos]) == 2))
			  FREE (cpl_table_delete, self->exts_tbs[i]);
			else if (gravi_pfits_get_extension_type (self->exts_hdrs[pos]) == 3)
			  FREE (cpl_imagelist_delete, self->exts_imgl[i]);
			CPLCHECK_MSG("Cannot delete data");

			FREE (cpl_propertylist_delete, self->exts_hdrs[i]);
			CPLCHECK_MSG("Cannot delete header");
		}
	}

	/* Cleanup further away */
	for (int i = comp ; i<GRAVI_DATA_SIZE; i++) {
			self->exts_hdrs[i] = NULL;
			self->exts_tbs[i] = NULL;
			self->exts_imgl[i] = NULL;
	}

	/* Change number of extension */
	self->nb_ext = comp;

	return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Erase an extension by its EXTNAME
 * 
 * @param self       The gravi_data to manipulate in-place
 * @param name       The string to search
 * 
 * The function erases the first (and only the first) extension
 * whose EXTNAME exactly matches @c name.
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_data_erase (gravi_data * self, const char * extname)
{
    cpl_ensure_code (self,    CPL_ERROR_NULL_INPUT);
    cpl_ensure_code (extname, CPL_ERROR_NULL_INPUT);

	/* Find the position of the property list */
	int pos = _gravi_data_find(self, extname);

	/* Check pos */
	if (pos == self->nb_ext) {
   	    cpl_msg_info (cpl_func,"Cannot delete '%s' (not found)",extname);
		return CPL_ERROR_NONE;
	} else {
	    cpl_msg_info (cpl_func,"Delete '%s' found in ext[%i] (over %i) ",extname, pos, self->nb_ext);
	}
	
	/* Erase */
	gravi_data_erase_x (self, pos);
	
	CPLCHECK_MSG ("Cannot erase this extension name");
	return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Erase all extension related to an instrument (SC, FT, MET...)
 * 
 * @param self       The gravi_data to manipulate in-place
 * @param type       The string to search ("SC", "FT", "MET"...)
 * 
 * The function erases all extension whose INSNAME or EXTNAME
 * parameters could be matched to @c type.
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_data_erase_type (gravi_data * self, const char * type)
{
  cpl_ensure_code (self, CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (type, CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (strlen (type)>1, CPL_ERROR_ILLEGAL_INPUT);

  /* Loop on extension */
  for (int ext = 0; ext < gravi_data_get_size (self) ; ext ++) {
	
	cpl_propertylist * plist = gravi_data_get_plist_x (self, ext);

	/* Check if EXTNAME or INSNAME contains the 'type' 
	 * Warning that the number of extension is changed in the loop */
	if ( (cpl_propertylist_has (plist, "INSNAME") &&
		  strstr (gravi_pfits_get_insname (plist), type) ) ||
		 (cpl_propertylist_has (plist, "EXTNAME") &&
		  strstr (gravi_pfits_get_extname (plist), type) ) )
    {
          cpl_msg_info (cpl_func,"Delete '%s' ", gravi_pfits_get_extname (plist));
          gravi_data_erase_x (self, ext);
          ext--;

          CPLCHECK_MSG ("Cannot erase this type");
    }
  }
  return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Add a BINTABLE extension in gravi_data
 * 
 * @param self       The gravi_data to insert in
 * @param plist      The (optional) plist to associate, not duplicated
 * @param extname    The (optional) EXTNAME
 * @param imglist    The table to insert, not duplicated
 * 
 * The table and its associated plist are append into the gravi_data, without
 * duplication (table and plist shall *not* be deleted).
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_data_add_table (gravi_data * self,
                                     cpl_propertylist * plist,
                                     const char * extname,
                                     cpl_table * table)
{
    cpl_ensure_code (self,  CPL_ERROR_NULL_INPUT);
    cpl_ensure_code (table, CPL_ERROR_NULL_INPUT);

	/* Set the header */
    if (plist == NULL) plist = cpl_propertylist_new ();
    if (extname != NULL) cpl_propertylist_update_string (plist, "EXTNAME", extname);
	
	/* Add header */
	cpl_propertylist_update_string (plist, "XTENSION", "BINTABLE");
    self->exts_hdrs[self->nb_ext] = plist;
	
	/* Add an EXTNAME to UNKNOWN */
	if (!cpl_propertylist_has(plist, "EXTNAME") ) {
        cpl_msg_error (cpl_func,"FIXME: set a table without EXTNAME !!");
        cpl_propertylist_update_string (plist, "EXTNAME", "UNKNOWN");
    }
	
	/* OIFITS: add OI_REVN to OI_* tables */
	const char * plist_name = 0;
	if (cpl_propertylist_has (plist, "EXTNAME"))
	  plist_name = gravi_pfits_get_extname (plist);
	
	if (plist_name && !(strcmp(plist_name,  GRAVI_OI_ARRAY_EXT) &&
						strcmp (plist_name, GRAVI_OI_TARGET_EXT) &&
						strcmp (plist_name, GRAVI_OI_WAVELENGTH_EXT) &&
						strcmp (plist_name, GRAVI_OI_T3_EXT) &&
						strcmp (plist_name, GRAVI_OI_VIS2_EXT) &&
						strcmp (plist_name, GRAVI_OI_VIS_EXT))) {
		cpl_propertylist_update_int (plist, "OI_REVN", 1);
	}
	
	/* OIFITS: check is saved as short */
	gravi_data_check_savetypes (plist, table);
	
	self->exts_tbs[self->nb_ext] = table;
	self->exts_imgl[self->nb_ext] = NULL;
	
	self->nb_ext ++;
	
	CPLCHECK_MSG ("Cannot add the table");
	return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Add an IMAGE (imagelist) extension in gravi_data
 * 
 * @param self       The gravi_data to insert in
 * @param plist      The (optional) plist to associate, not duplicated
 * @param extname    The (optional) EXTNAME
 * @param imglist    The imagelist to insert, not duplicated
 * 
 * The image and its associated plist are append into the gravi_data, without
 * duplication (image and plist shall *not* be deleted).
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_data_add_cube (gravi_data * self,
                                    cpl_propertylist * plist,
                                    const char * extname,
                                    cpl_imagelist * imglist)
{
  
    cpl_ensure_code (self,  CPL_ERROR_NULL_INPUT);
    cpl_ensure_code (imglist, CPL_ERROR_NULL_INPUT);

    if (plist == NULL) plist = cpl_propertylist_new ();
    if (extname != NULL) cpl_propertylist_update_string (plist, "EXTNAME", extname);

	/* Add header */
	self->exts_hdrs[self->nb_ext] = plist;
	cpl_propertylist_update_string (self->exts_hdrs[self->nb_ext],"XTENSION", "IMAGE");
	
	/* Add data */
	self->exts_imgl[self->nb_ext] = (imglist);
	self->exts_tbs[self->nb_ext] = NULL;
	
	self->nb_ext ++;

	CPLCHECK_MSG ("Cannot add the cube");
	return CPL_ERROR_NONE;

}

/*----------------------------------------------------------------------------*/
/**
 * @brief Add an IMAGE (single image) extension in gravi_data
 * 
 * @param self       The gravi_data to insert in
 * @param plist      The (optional) plist to associate, not duplicated
 * @param extname    The (optional) EXTNAME
 * @param image      The image to insert, not duplicated
 * 
 * The image and its associated plist are append into the gravi_data, without
 * duplication (image and plist shall *not* be deleted).
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_data_add_img (gravi_data * self,
                                   cpl_propertylist * plist,
                                   const char * extname,
                                   cpl_image * image) {
  
    cpl_ensure_code (self,  CPL_ERROR_NULL_INPUT);
    cpl_ensure_code (image, CPL_ERROR_NULL_INPUT);

    if (plist == NULL) plist = cpl_propertylist_new ();
    if (extname != NULL) cpl_propertylist_update_string (plist, "EXTNAME", extname);

	/* Add header */
	self->exts_hdrs[self->nb_ext] = plist;
	cpl_propertylist_update_string (self->exts_hdrs[self->nb_ext],"XTENSION", "IMAGE");

    /* Convert to imagelist */
    cpl_imagelist * imglist = cpl_imagelist_new ();
    cpl_imagelist_set (imglist, image, 0);

	/* Add data */
	self->exts_imgl[self->nb_ext] = (imglist);
	self->exts_tbs[self->nb_ext] = NULL;
	
	self->nb_ext ++;

	CPLCHECK_MSG ("Cannot add the image");
	return CPL_ERROR_NONE;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Clean the data to keep only OIFITS extensions related to SC.
 *
 * @param self   The gravi data to be cleaned in-place
 *
 * Remove all tables not related to SC. That is keep: OI_ARRAY, OI_TARGET,
 * and keep OI_WAVELENGTH, OI_VIS, OI_VIS2, OI_T3, OI_FLUX whose INSNAME
 * is the SC.
 */
/*---------------------------------------------------------------------------*/

cpl_error_code gravi_data_clean_for_astro (gravi_data * data)
{
	gravi_msg_function_start(1);
    cpl_ensure_code (data, CPL_ERROR_NULL_INPUT);

	/* Loop on extension in file ?? */
	for (int i = 0; i < gravi_data_get_size (data); i++) {

	    /* Keep OI_ARRAY and OI_TARGET */
		cpl_propertylist * plist = gravi_data_get_plist_x (data, i);
		const char * plist_name = gravi_pfits_get_extname (plist);
		if (!(strcmp (plist_name, GRAVI_OI_ARRAY_EXT)) ||
			!(strcmp (plist_name, GRAVI_OI_TARGET_EXT)) ) {
		  cpl_msg_debug (cpl_func,"NAME: %s kept", plist_name);
		  continue;
		}

		/* Keep all INSNAME_SC */
		if (cpl_propertylist_has (plist, "INSNAME") &&
			( (!strcmp (gravi_pfits_get_insname (plist), INSNAME_SC_P1)) ||
			  (!strcmp (gravi_pfits_get_insname (plist), INSNAME_SC_P2)) ||
			  (!strcmp (gravi_pfits_get_insname (plist), INSNAME_SC)) ) ) {
		  cpl_msg_debug (cpl_func,"NAME: %s kept", plist_name);
		  continue;
		}

		/* Delete */
		cpl_msg_debug (cpl_func,"NAME: %s deleted", plist_name);
		FREE (cpl_propertylist_delete, data->exts_hdrs[i]);
		FREE (cpl_table_delete, data->exts_tbs[i]);
		FREE (cpl_imagelist_delete, data->exts_imgl[i]);
	}

	/* Loop on extension in file to move the
	 * extension consecutively */
	int j = 0;
	for (int i = 0; i < gravi_data_get_size (data); i++) {
	  /* This one exist, copy in last place */
	  if (data->exts_hdrs[i] || data->exts_tbs[i] || data->exts_imgl[i]) {
		if ( i!=j ) {
		  data->exts_hdrs[j] = data->exts_hdrs[i]; data->exts_hdrs[i] = NULL;
		  data->exts_tbs[j] = data->exts_tbs[i]; data->exts_tbs[i] = NULL;
		  data->exts_imgl[j] = data->exts_imgl[i]; data->exts_imgl[i] = NULL;
		}
		j++;
	  }
	}

	data->nb_ext = j;

	gravi_msg_function_exit(1);
	return CPL_ERROR_NONE;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Verify the integrity of RAW data
 *
 * @param data: RAW data to verify
 * 
 * Check the consistency of table in the data, by looking at the
 * temporal coverage of the real-time signals METROLOGY, OPDC, FT...
 * (shall encompass the science), and the steps between consecutive
 * frames in these real-time signals (shall be constant).
 */
/*---------------------------------------------------------------------------*/

cpl_error_code gravi_data_check_consistency (gravi_data * data)
{
  gravi_msg_function_start(0);
  cpl_ensure_code (data, CPL_ERROR_NULL_INPUT);
  
  char qc_msg[80];

  /* Add the QC parameter, start with 0 warning */
  cpl_propertylist * header = gravi_data_get_header (data);
  cpl_propertylist_append_int (header, "ESO QC CHECK FLAGS", 0);
  
  /* Check number of extensions */
  int nb_ext = 12;
  const char *extname[] = {"ARRAY_DESCRIPTION", "ARRAY_GEOMETRY", "OPTICAL_TRAIN",
						   "IMAGING_DATA_ACQ", "ACQ_ABS_REF_POSITION",
						   "OPDC", "FDDL", "METROLOGY",
						   "IMAGING_DATA_SC","IMAGING_DETECTOR_SC",
						   "IMAGING_DATA_FT","IMAGING_DETECTOR_FT"};
  
  for (int ext = 0; ext< nb_ext; ext++)  {
	if ( !gravi_data_has_extension (data, extname[ext]) ) {
	  sprintf (qc_msg, "%s is missing", extname[ext]);
	  gravi_pfits_add_check (header, qc_msg);
	}
  }

  CPLCHECK_MSG ("Cannot check consistency of tables...");
   
  /* Build the time of the first and last frame of SC in
   * [us] with respect to PRC.ACQ.START. We need RMN data over
   * the entire first and last DITs of SC */
  int    ndit     = cpl_propertylist_get_int (header, "ESO DET2 NDIT");
  double scdit    = gravi_pfits_get_sc_dit (header);
  double start_sc = gravi_pfits_get_sc_time (header, 0) - scdit/2. * 1e6;
  double end_sc   = gravi_pfits_get_sc_time (header, ndit-1) + scdit/2. * 1e6;
  
  cpl_msg_debug (cpl_func,"start = %8.0f and end = %10.0f (%s [us] from PCR.ACQ.START)",
                 start_sc, end_sc, "SC");
  
  /* Check the TIME content of RMN tables */
  int nb_ext_rmn = 4;
  const char *ext_rmn[] = {"OPDC", "FDDL", "METROLOGY", "IMAGING_DATA_FT"};
  
  for (int ext = 0; ext< nb_ext_rmn; ext++) {

	/* Get this RMN table */
	if (!gravi_data_has_extension (data, ext_rmn[ext])) continue;
	cpl_table * table = gravi_data_get_table (data, ext_rmn[ext]);

	/* Get the first and last time sample, in [us] */
	double min_table = cpl_table_get_column_min (table, "TIME");
	double max_table = cpl_table_get_column_max (table, "TIME");

	cpl_msg_debug (cpl_func,"start = %8.0f and end = %10.0f (%s [us] from PCR.ACQ.START)",
                   min_table, max_table, ext_rmn[ext]);

	/* Check if it covers the SC data */
	if ( min_table > start_sc) {
	  sprintf (qc_msg, "%s starts *after* the science exposure", ext_rmn[ext]);
	  gravi_pfits_add_check (header, qc_msg);
	}
	if ( max_table < end_sc) {
	  sprintf (qc_msg, "%s finish *before* the science exposure", ext_rmn[ext]);
	  gravi_pfits_add_check (header, qc_msg);
	}

	/* Get data */
	cpl_size nrow = cpl_table_get_nrow (table);
	int * time = cpl_table_get_data_int (table, "TIME");
	CPLCHECK_MSG ("Cannot get data");

	/* Compute the median step, over the first 10000 samples */
	cpl_size nsamp = CPL_MIN (10000, nrow-1);
	cpl_array * delta = cpl_array_new (nsamp, CPL_TYPE_INT);
	for (cpl_size row = 0; row < nsamp ; row++) {
	  cpl_array_set_int (delta, row, time[row+1] - time[row]);
	}
	double median_delta = cpl_array_get_median (delta);
	FREE (cpl_array_delete, delta);

	/* Check missing samples (delta > 1.5 * median) */
	int nwrong = 0;
	for (cpl_size row = 0; row < nrow-1 ; row++) {
	 double current_delta = time[row+1] - time[row];
	 if (current_delta > 1.5 * median_delta) nwrong ++;
	}
	
	if (nwrong) {
	  sprintf (qc_msg, "%s has %i wrong steps", ext_rmn[ext], nwrong);
	  gravi_pfits_add_check (header, qc_msg);
	}
	
  } /* End loop on RMN tables */

  CPLCHECK_MSG ("Cannot check consistency of RMN tables...");
	
  gravi_msg_function_exit(0);
  return CPL_ERROR_NONE;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Get pointer to the OI_WAVELENGTH tables of both polarisations
 *
 * @param data:         Input gravi_data
 * @param type_data:    GRAVI_SC or GRAVI_FT
 * @param npol:         Number of polarisation in data (1 or 2)
 * @return A pointer to the one or two OI_WAVELENGHT table
 *
 * The routine allocate memory for 1 or 2 pointers on the tables.
 * Thus the return value shall be desallocated with cpl_free
 */
/*---------------------------------------------------------------------------*/

cpl_table ** gravi_data_get_oiwave_tables (gravi_data * data, int type_data, int npol)
{
    cpl_ensure (data, CPL_ERROR_NULL_INPUT, NULL);    
    cpl_ensure (npol == 1 || npol == 2, CPL_ERROR_ILLEGAL_INPUT, NULL);
    cpl_ensure (type_data == GRAVI_SC || type_data == GRAVI_FT, CPL_ERROR_ILLEGAL_INPUT, NULL);

    cpl_table ** oiwave_tables = cpl_calloc (npol, sizeof(cpl_table *));
    for (int pol = 0; pol < npol; pol++) 
        oiwave_tables[pol] = gravi_data_get_oi_wave (data, type_data, pol, npol);

    return oiwave_tables;
}



/**@}*/
