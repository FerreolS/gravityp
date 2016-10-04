/* $Id: gravi_data-test.c,v 1.59 2011/08/16 17:43:49 nazouaoui Exp $
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
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */


/*
 * gravi-test.c
 *
 *  Created on: 17 août 2011
 *      Author: nabih
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <assert.h>

#include "gravi-test.h"

#define COMPUTE_FILES	0


int gravi_array_compare(cpl_array * array_1, cpl_array * array_2){

	if ((array_1 == NULL) || (array_2 == NULL))
		return 1;

//	else if ((array_1 == NULL) && (array_2 != NULL)){
//		return 1;
//	}
//	else if ((array_1 != NULL) && (array_2 == NULL)){
//		return 1;
//	}
	cpl_type type_array_1 = cpl_array_get_type(array_1);
	cpl_type type_array_2 = cpl_array_get_type(array_2);

	if (strcmp(cpl_type_get_name (type_array_1), cpl_type_get_name (type_array_2))){
		cpl_msg_info(cpl_func,"The arrays have not the same type type1 = %s and type2 = %s\n",
				cpl_type_get_name (type_array_1), cpl_type_get_name (type_array_2));
		return 0;
	}

	if (cpl_array_get_size(array_1) != cpl_array_get_size(array_2)){
		cpl_msg_info(cpl_func,"the arrays don't have the same size\n");
		return 0;
	}

	int size = cpl_array_get_size(array_1);
	int j;

	/* compute the image list depending of the DATA type */
	double * array_double_1;
	double * array_double_2;
	double complex * array_complex_1;
	double complex * array_complex_2;
	int * array_int_1;
	int * array_int_2;
	float * array_float_1;
	float * array_float_2;
	char ** array_string_1;
	char ** array_string_2;
	switch (type_array_1)
	{
	    case CPL_TYPE_DOUBLE :

	    	array_double_1 = cpl_array_get_data_double(array_1);
	    	array_double_2 = cpl_array_get_data_double(array_2);
	        for (j = 0; j < size ; j++){
	        	if (array_double_1[j]  > 1e-10){
					if (fabs((array_double_1[j] - array_double_2[j]) / array_double_1[j]) > 0.1){
						cpl_msg_info(cpl_func,"The arrays are not the same on the %d nd index:"
								" %e and %e\n", j, array_double_1[j], array_double_2[j]);
						return 0;
					}
	        	}
	        	else {
	        		if (array_double_2[j] > 1e-10){
						cpl_msg_info(cpl_func,"The arrays are not equal at zero on the %d nd index:"
								" %e and %e\n", j, array_double_1[j], array_double_2[j]);
						return 0;
	        		}
	        	}


	        }

	    break;

	    case CPL_TYPE_DOUBLE_COMPLEX :

	    	array_complex_1 = cpl_array_get_data_double_complex(array_1);
	    	array_complex_2 = cpl_array_get_data_double_complex(array_2);
	        for (j = 0; j < size ; j++){
	        	if (cabs(array_complex_1[j]) > 1e-10){
					if ((carg(array_complex_1[j]) - carg(array_complex_2[j])) > 0.1){
						cpl_msg_info(cpl_func,"The argument of complex arrays are not the same on the %d nd index:"
								" %e and %e\n", j,carg( array_complex_1[j]), carg(array_complex_2[j]));
						return 0;
					}
					if (fabs((cabs(array_complex_1[j]) - cabs(array_complex_2[j]))) /
							cabs(array_complex_1[j]) > 0.1){
						cpl_msg_info(cpl_func,"The module of complex arrays are not the same on the %d nd index:"
								" %e and %e\n", j, cabs(array_complex_1[j]), cabs(array_complex_2[j]));
						return 0;
					}
	        	}
	        	else {
	        		if (cabs(array_complex_2[j]) > 1e-10) {
						cpl_msg_info(cpl_func,"The module of complex arrays are not equal at zero on the %d nd index:"
								" %e and %e\n", j, cabs(array_complex_1[j]), cabs(array_complex_2[j]));
						return 0;
	        		}
	        	}
	        }

	    break;

	    case CPL_TYPE_INT :
	    	array_int_1 = cpl_array_get_data_int(array_1);
			array_int_2 = cpl_array_get_data_int(array_2);
			for (j = 0; j < size ; j++){
				if (array_int_1[j] != 0) {
					if ((abs(array_int_1[j] - array_int_2[j]) / array_int_1[j]) > 1){
						cpl_msg_info(cpl_func,"The arrays are not the same on the %d nd index: "
								"%d and %d\n", j, array_int_1[j], array_int_2[j]);
						return 0;
					}
				}
				else {
					if (array_int_2[j] > 1e-10){
						cpl_msg_info(cpl_func,"The arrays are not equal at zero on the %d nd index: "
								"%d and %d\n", j, array_int_1[j], array_int_2[j]);
						return 0;
					}
				}


			}


	    break;

	    case CPL_TYPE_FLOAT :

	    	array_float_1 = cpl_array_get_data_float(array_1);
			array_float_2 = cpl_array_get_data_float(array_2);
			for (j = 0; j < size ; j++){
				if (array_float_1[j] > 1e-10){
					if (fabs((array_float_1[j] - array_float_2[j]) / array_float_1[j])> 0.1){
						cpl_msg_info(cpl_func,"The arrays are not the same on the %d nd index: "
								"%e and %e\n", j, array_float_1[j], array_float_2[j]);
						return 0;
					}
				}
				else {
					if (array_float_2[j] > 1e-10){
						cpl_msg_info(cpl_func,"The arrays are not equal at zero on the %d nd index: "
								"%e and %e\n", j, array_float_1[j], array_float_2[j]);
						return 0;
					}
				}

			}

	    break;

	    case CPL_TYPE_STRING :
	    	array_string_1 = cpl_array_get_data_string(array_1);
	    	array_string_2 = cpl_array_get_data_string(array_2);
	    	for (j = 0; j < size; j++){
	    		if(strcmp(*(array_string_1 + j), *(array_string_2 + j))){
	    			cpl_msg_info(cpl_func,"The arrays are not the same on the %d nd index: "
	    					"%s and %s\n", j, array_string_1[j], array_string_2[j]);
	    			return 0;
	    		}
	    	}
	    	break;

	    default:

	        cpl_msg_info(cpl_func,"invalid type of array %s\n", cpl_type_get_name (type_array_2));
	        return -1;
	    break;
	}

	return 1;

}


int gravi_table_compare(cpl_table * table1, cpl_table * table2 ){

	if (cpl_table_compare_structure(table1, table2) == 0){
		int ncol = cpl_table_get_ncol(table1);
		int nrow = cpl_table_get_nrow(table1);
		cpl_array * table_names = cpl_table_get_column_names (table1);
		int i, j;
		for (i = 0; i < ncol; i++){
			const char * col_name = cpl_array_get_string(table_names, i);

            /* Cast the type on int to avoid warning */
            int type = cpl_table_get_column_type(table1, col_name);
			cpl_array ** array_table1;
			cpl_array ** array_table2;
			double * double_table1;
			double * double_table2;
			int * int_table1;
			int * int_table2;
			float * float_table1;
			float * float_table2;
			char ** string_table1;
			char ** string_table2;

			switch(type){
				case CPL_TYPE_DOUBLE | CPL_TYPE_POINTER :
				case CPL_TYPE_POINTER | CPL_TYPE_INT:
				case CPL_TYPE_FLOAT | CPL_TYPE_POINTER :
				case CPL_TYPE_DOUBLE_COMPLEX | CPL_TYPE_POINTER :
				case CPL_TYPE_FLOAT_COMPLEX | CPL_TYPE_POINTER :
					array_table1 = cpl_table_get_data_array(table1, col_name);
					array_table2 = cpl_table_get_data_array(table2, col_name);
					for (j = 0; j < nrow; j++){
						if (gravi_array_compare(array_table1[j]
						                             , array_table2[j]) != 1){
							cpl_msg_info(cpl_func,"Case array : on the column %s and in the row %d"
									" the arrays are different", col_name, j);
							return 0;
						}
					}
				break;
				case CPL_TYPE_DOUBLE :
					double_table1 = cpl_table_get_data_double(table1, col_name);
					double_table2 = cpl_table_get_data_double(table2, col_name);
					for (j = 0; j < nrow; j++){
						if (double_table1[j] > 1e-10){
							if (fabs((double_table1[j] - double_table2[j]) / double_table1[j]) > 0.1){
								cpl_msg_info(cpl_func,"Case double : on the column %s and in the row %d"
										" the double are different", col_name, j);
								return 0;
							}
						}
						else {
							if (double_table2[j] > 1e-10){
								cpl_msg_info(cpl_func,"Case double : on the column %s and in the row %d"
										" the double are different", col_name, j);
								return 0;
							}
						}

					}
				break;
				case CPL_TYPE_INT :
					int_table1 = cpl_table_get_data_int(table1, col_name);
					int_table2 = cpl_table_get_data_int(table2, col_name);
					for (j = 0; j < nrow; j++){
						if (int_table1[j] != 0){
							if (abs(int_table1[j] - int_table2[j]) / int_table1[j]  > 1){
								cpl_msg_info(cpl_func,"Case int : on the column %s and in the row %d"
										" the integers are different\n", col_name, j);
								return 0;
							}
						}
						else {
							if (int_table2[j] != 0){
								cpl_msg_info(cpl_func,"Case int : on the column %s and in the row %d"
										" the integers are different\n", col_name, j);
								return 0;
							}
						}

					}
				break;

				case CPL_TYPE_FLOAT :
					float_table1 = cpl_table_get_data_float(table1, col_name);
					float_table2 = cpl_table_get_data_float(table2, col_name);
					for (j = 0; j < nrow; j++){
						if (float_table1[j] > 1e-10){
							if (fabs((float_table1[j] - float_table2[j]) / float_table1[j]) > 0.1){
								cpl_msg_info(cpl_func,"Case int : on the column %s and in the row %d"
										" the integers are different\n", col_name, j);
								return 0;
							}
						}
						else {
							if (float_table2[j] > 1e-10){
								cpl_msg_info(cpl_func,"Case int : on the column %s and in the row %d"
										" the integers are different\n", col_name, j);
								return 0;
							}
						}
					}
				break;

				case CPL_TYPE_STRING :
					string_table1 = cpl_table_get_data_string(table1, col_name);
					string_table2 = cpl_table_get_data_string(table2, col_name);
					for (j = 0; j < nrow; j++){
						if (strcmp(string_table1[j], string_table2[j]) != 0){
							cpl_msg_info(cpl_func,"Case string : on the column %s and in the row %d"
									" the strings are different\n", col_name, j);
							return 0;
						}
					}
				break;

				default:

					cpl_msg_info(cpl_func,"invalid type %s\n", cpl_type_get_name (type));
					return -1;
				break;
			}
		}
		cpl_array_delete(table_names);
		return 1;

	}
	else{
		cpl_msg_info(cpl_func,"The tables don't have the same number of columns, "
				"with the same names or the same types\n");
		return 0;
	}


}

int gravi_propertylist_compare(cpl_propertylist * plist1, cpl_propertylist * plist2){

	int size = cpl_propertylist_get_size (plist1);
	if ( size != cpl_propertylist_get_size(plist2)){
		cpl_msg_info(cpl_func,"the property lists don't have the same size\n");
		return 0;
	}
	int i;
	cpl_property * p;

	const char * name;
	for (i = 0; i < size; i++){
		p = cpl_propertylist_get(plist1, i);
		name = cpl_property_get_name(p);
		if(!cpl_propertylist_has (plist2, name)){
			cpl_msg_info(cpl_func,"The property %s is not present on the second "
					"property list \n", name);
			return 0;
		}

		cpl_type type1 = cpl_propertylist_get_type (plist1, name);
		cpl_type type2 = cpl_propertylist_get_type (plist2, name);
		if (strcmp(cpl_type_get_name(type1), cpl_type_get_name(type2))){
			cpl_msg_info(cpl_func,"The property %s in the first property list "
					"have not the same type as in the second property list \n"
					                                                    , name);
			return 0;
		}

		int p_int;
		double p_double;
		float p_float;
		const char * p_string;
		long p_long;
		char p_char;
		switch(type1) {
		case CPL_TYPE_DOUBLE :
			p_double = cpl_propertylist_get_double(plist1, name);
			if (p_double != cpl_propertylist_get_double(plist2, name)){
				cpl_msg_info(cpl_func,"case Double : the property %s have not the same value\n", name);
				return 0;
			}
			break;
		case CPL_TYPE_INT :
			p_int = cpl_propertylist_get_int(plist1, name);
			if (p_int != cpl_propertylist_get_int(plist2, name)){
				cpl_msg_info(cpl_func,"case interger : the property %s have not the same value\n", name);
				return 0;
			}
			break;
		case CPL_TYPE_FLOAT :
			p_float = cpl_propertylist_get_float(plist1, name);
			if (p_float != cpl_propertylist_get_float(plist2, name)){
				cpl_msg_info(cpl_func,"case float : the property %s have not the same value\n", name);
				return 0;
			}
			break;
		case CPL_TYPE_STRING :
			p_string = cpl_propertylist_get_string(plist1, name);
			if (strcmp(p_string, cpl_propertylist_get_string(plist2, name))){
				cpl_msg_info(cpl_func,"case string : the property %s have not the same value\n", name);
				return 0;
			}
			break;
		case CPL_TYPE_LONG :
			p_long = cpl_propertylist_get_long(plist1, name);
			if (p_long != cpl_propertylist_get_long(plist2, name)){
				cpl_msg_info(cpl_func,"case long : the property %s have not the same value\n", name);
				return 0;
			}
			break;
		case CPL_TYPE_CHAR :
			p_char = cpl_propertylist_get_char(plist1, name);
			if (p_char != cpl_propertylist_get_char(plist2, name)){
				cpl_msg_info(cpl_func,"case char : the property %s have not the same value\n", name);
				return 0;
			}
			break;
		default :
			cpl_msg_info(cpl_func,"invalid type of property %s\n", cpl_type_get_name (type1));
			return -1;
			break;

		}

	}

	return 1;
}

int gravi_data_compare(gravi_data * data1, gravi_data * data2){


	int size = gravi_data_get_size(data1);
	if (size != gravi_data_get_size(data1)){
		cpl_msg_info(cpl_func,"The elements don't have the same size\n");
		return 0;
	}
	cpl_propertylist * p1;
	p1 = gravi_data_get_plist(data1, GRAVI_PRIMARY_HDR_EXT);
	cpl_propertylist * p2;
	p2 = gravi_data_get_plist(data2, GRAVI_PRIMARY_HDR_EXT);
	if (gravi_propertylist_compare(p1, p2) != 1){
		cpl_msg_info(cpl_func,"The %s property list is diffetent "
				     "between the two elements\n", GRAVI_PRIMARY_HDR_EXT);
		return 0;
	}

	int i, ext;
	const char * name;
	for(i = 0; i < size; i++){
		p1 = gravi_data_get_plist_x(data1, i);//data1->exts_hdrs[i];
		name = cpl_propertylist_get_string(p1, "EXTNAME");
		ext = 0;
		while( ext != size ){
			p2 = gravi_data_get_plist_x(data1, ext);//data2->exts_hdrs[ext];
			if (!strcmp(name, cpl_propertylist_get_string(p2, "EXTNAME"))){
				break;
			}
			ext++;
		}
		if (ext == size){
			cpl_msg_info(cpl_func,"The seconde element don't contain the %s "
					       "property list\n", name);
			return 0;
		}
//		if (gravi_propertylist_compare(data1->exts_hdrs[i],
//				                          data2->exts_hdrs[ext]) != 1){
		if (gravi_propertylist_compare(p1, p2) != 1){
			cpl_msg_info(cpl_func,"The %s property list is not the same for the "
					"twe elements\n", name);
			return 0;
		}
//		if (gravi_table_compare(data1->exts_tbs[i], data2->exts_tbs[ext]) != 1){
		if (gravi_table_compare(gravi_data_get_table_x(data1, i),
				gravi_data_get_table_x(data2, ext)) != 1){
			cpl_msg_info(cpl_func,"The table associated to %s property list is not"
					" the same for the two elements\n", name);
			return 0;
		}

	}
	return 1;

}


