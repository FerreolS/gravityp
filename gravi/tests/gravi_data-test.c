/* $Id: gravi_data-test.c,v 1.59 2011/08/04 17:06:49 nazouaoui Exp $
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
 * $Author: nazouaoui $
 * $Date: 2011/08/04 17:06:49 $
 * $Revision: 1.0 $
 * $Name: gravi $
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <assert.h>

#include <cpl_test.h>
#include <cpl.h>
#include "gravi-test.c"
#include "gravi_data.h"
#include "gravi_utils.h"

/* *+

The VERBOSE macro is no longer used, instead
export CPL_MSG_LEVEL=info
can be used at run-time.
(At compile time, the default message level can be
controlled via cpl_test_init()).

#define VERBOSE
+* */

void printResetError(char* file, int line);
int gravi_data_test(void);


void printResetError(char* file, int line)
{
	if (cpl_error_get_code())
	{
		printf("Error in file %s (%d): message \"%s\" in function %s\n",
				file, line, cpl_error_get_message(), cpl_error_get_function());
		cpl_error_reset();
	}
}

/*
 * gravi_data-test.C
 *
 *  Created on: 4 août 2011
 *      Author: nabih
 */


int gravi_data_test(){
	gravi_data * data1, * data2;
	cpl_propertylist *plist1, *plist2;
	cpl_table * table1, *table2, * table3, *table_test;
    const char *keys[] = {
        "a", "b", "c", "d", "e", "f", "g"
    };
    int flag=EXIT_SUCCESS;

    const char *comments[] = {
        "A character value",
        "A boolean value",
        "A integer value",
        "A long integer value",
        "A floating point number",
        "A double precision number",
        "A string value"
    };
//	cpl_test_init(PACKAGE_BUGREPORT, CPL_MSG_WARNING);

	/*
	 *  Testing a new data with empty elements.
	 */

	test_data(data1, gravi_data_new(0),
			                        "Creating a data with empty elements... ", flag);
	printResetError(__FILE__, __LINE__);
	test_ivalue(0, gravi_data_get_size(data1), "Check zero data element... ", flag);
	printResetError(__FILE__, __LINE__);

	/*
	 *  Testing a data accessors.
	 */
	plist1 = cpl_propertylist_new();

    cpl_propertylist_append_string(plist1, "EXTNAME", "Plist1");
    cpl_propertylist_append_string(plist1, "XTENSION", "TABLE");
	cpl_propertylist_append_char(plist1, keys[0], 'a');
    cpl_propertylist_set_comment(plist1, keys[0], comments[0]);

//    cpl_propertylist_append_bool(plist1, keys[1], 0);
//    cpl_propertylist_set_comment(plist1, keys[1], comments[1]);

    cpl_propertylist_append_int(plist1, keys[2], -1);
    cpl_propertylist_set_comment(plist1, keys[2], comments[2]);

    plist2 = cpl_propertylist_new();
    cpl_propertylist_append_string(plist2, "EXTNAME", "Plist2");
    cpl_propertylist_append_string(plist2, "XTENSION", "TABLE");
    cpl_propertylist_append_float(plist2, keys[4], -1.2345);
    cpl_propertylist_set_comment(plist2, keys[4], comments[4]);

    cpl_propertylist_append_double(plist2, keys[5], 2.456987);
    cpl_propertylist_set_comment(plist2, keys[5], comments[5]);

    cpl_propertylist_append_string(plist2, keys[6], comments[6]);
    cpl_propertylist_set_comment(plist2, keys[6], comments[6]);

    table1 = cpl_table_new(10);
    table2 = cpl_table_new(10);
    table3 = cpl_table_new(10);
    cpl_table_new_column(table1, "Logical", CPL_TYPE_INT);
    cpl_table_set_int(table1, "Logical", 0, 1);
    cpl_table_set_int(table1, "Logical", 1, 0);
    cpl_table_set_int(table1, "Logical", 2, -1);
    cpl_table_set_int(table1, "Logical", 3, 0);
    cpl_table_set_int(table1, "Logical", 4, 0);
    cpl_table_set_int(table1, "Logical", 5, 0);
    cpl_table_set_int(table1, "Logical", 6, 0);
    cpl_table_set_int(table1, "Logical", 7, 0);
    cpl_table_set_int(table1, "Logical", 8, 0);
    cpl_table_set_int(table1, "Logical", 9, 1);
    cpl_table_set_column_savetype(table1, "Logical", CPL_TYPE_BOOL);
    cpl_table_new_column(table2, "Byte", CPL_TYPE_INT);
    cpl_table_set_int(table2, "Byte", 0, 1);
    cpl_table_set_int(table2, "Byte", 1, 30);
    cpl_table_set_int(table2, "Byte", 2, -50);
    cpl_table_set_int(table2, "Byte", 3, 0);
    cpl_table_set_int(table2, "Byte", 4, 130);
    cpl_table_set_int(table2, "Byte", 5, -130);
    cpl_table_set_int(table2, "Byte", 6, 0);
    cpl_table_set_int(table2, "Byte", 7, 0);
    cpl_table_set_int(table2, "Byte", 8, 0);
    cpl_table_set_int(table2, "Byte", 9, 1);
    cpl_table_set_column_savetype(table2, "Byte", CPL_TYPE_CHAR);
    cpl_table_new_column(table3, "Ubyte", CPL_TYPE_INT);
    cpl_table_set_int(table3, "Ubyte", 0, 1);
    cpl_table_set_int(table3, "Ubyte", 1, 30);
    cpl_table_set_int(table3, "Ubyte", 2, -50);
    cpl_table_set_int(table3, "Ubyte", 3, 0);
    cpl_table_set_int(table3, "Ubyte", 4, 130);
    cpl_table_set_int(table3, "Ubyte", 5, -130);
    cpl_table_set_int(table3, "Ubyte", 6, 0);
    cpl_table_set_int(table3, "Ubyte", 7, 258);
    cpl_table_set_int(table3, "Ubyte", 8, 0);
    cpl_table_set_int(table3, "Ubyte", 9, 1);
    cpl_table_set_column_savetype(table3, "Ubyte", CPL_TYPE_UCHAR);

    test(gravi_data_add_table(data1, plist1, NULL, table1),
    		                          "Add the first element on the data... ", flag);
    printResetError(__FILE__, __LINE__);
    test_ivalue(1, gravi_data_get_size(data1),
    		                        "Check there are one element on data... ", flag);
    printResetError(__FILE__, __LINE__);

    test(gravi_data_add_table(data1, plist2, NULL, table2),
    		                          "Add the second element on the data... ", flag);
    printResetError(__FILE__, __LINE__);

    test_ivalue(2, gravi_data_get_size(data1),
    		                        "Check there are two elements on data... ", flag);
    printResetError(__FILE__, __LINE__);

    test_data(table_test, gravi_data_get_table(data1, "Plist1"),
    		                                    "Get the table from data... ", flag);
    printResetError(__FILE__, __LINE__);


    test_ivalue(1, gravi_table_compare(table_test, table1), "Check that the two "
			                                  "tables are the same... " , flag);
    printResetError(__FILE__, __LINE__);

    test(gravi_data_add_table(data1, NULL, "Plist3", table3),
    		                          "Replace the first table by another...", flag);

    test_data(table_test, gravi_data_get_table(data1, "Plist3"),
    		                               "Get the table from data... ", flag);
    test_ivalue(1, gravi_table_compare(table_test, table3), "Check that the two "
                                                "tables are the same... " , flag);

    /*
	 *  Testing functions using data.
	 */

    test(gravi_data_erase(data1, "Plist2"), "Erase an element from data... ", flag);

    test_ivalue(2, gravi_data_get_size(data1), "Check the new size of data...", flag);

    test_data(data2, gravi_data_duplicate(data1), "Duplicate the data... ", flag);
    test_ivalue(1, gravi_data_compare(data1, data2), "Check no difference between the structures... ", flag);

    gravi_data_delete(data2);

    /*
	 * Load and save data.
	 */

    test_data(data2, gravi_data_load(DATADIR "Dark.fits"),
    		"Load the data from a FITS file... ", flag);

    /* Construction of the frame set */

	cpl_frame * _frame;
	cpl_frameset * frameset;

	/* Add frame to the frame set */

	frameset = cpl_frameset_new();
	_frame = cpl_frame_new();

	cpl_frame_set_filename(_frame, DATADIR "Dark.fits");
	cpl_frame_set_tag(_frame, "MASTER_BIAS");
	cpl_frame_set_type(_frame, CPL_FRAME_TYPE_IMAGE);
	cpl_frame_set_group(_frame, CPL_FRAME_GROUP_RAW);
	cpl_frame_set_level(_frame, CPL_FRAME_LEVEL_FINAL);

	cpl_frameset_insert(frameset, _frame);

    /* Construction of the parameter list */
	cpl_parameter *p[2];

    cpl_parameterlist *list = NULL;

    /* Create a parameter list */

    list = cpl_parameterlist_new();

    /* Append parameters to the list */

    p[0] = cpl_parameter_new_value("a", CPL_TYPE_STRING, "parameter1",
                                   "None", "value1");
    cpl_parameterlist_append(list, p[0]);

    p[1] = cpl_parameter_new_value("b", CPL_TYPE_STRING, "parameter2",
                                   "None", "value2");
    cpl_parameterlist_append(list, p[1]);

    /* Construction of the property list to add on the FITS file */

    cpl_propertylist * applist = cpl_propertylist_new();

    cpl_propertylist_append_string(applist, "OBJECTF", "Testing the data "
    		                                   "functions load and save");
    cpl_propertylist_append_string(applist, CPL_DFS_PRO_CATG, "DARK");
    test(gravi_data_save_data(data2, "test_data_save_data.fits", CPL_IO_CREATE),
    		"Save the data... ", flag);

    test(gravi_data_save_new(data2, frameset, "test_data_save.fits", NULL, list,
        frameset, _frame, "test_save", applist, "DARK"), "Save the data with all options... ", flag);

//    cpl_end();
//    exit(EXIT_SUCCESS);

	cpl_propertylist_delete(applist);
    cpl_parameterlist_delete(list);
	cpl_frameset_delete(frameset);
	gravi_data_delete(data1);
	gravi_data_delete(data2);
    unlink("test_data_save_data.fits");
    unlink("test_data_save_dark.fits");

    return 0;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Unit tests of gravi_data module
 */
/*----------------------------------------------------------------------------*/

int main(void)
{
    int flag;

    cpl_test_init(PACKAGE_BUGREPORT, CPL_MSG_INFO);

    flag=gravi_data_test();

    if (flag == EXIT_FAILURE)
    {
    	cpl_test_end(0);
    	exit(EXIT_FAILURE);
    }
    cpl_test_end(0);
    exit(EXIT_SUCCESS);
}
