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

#ifndef GRAVI_TEST_C
#define GRAVI_TEST_C


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
#include <gravi_data.h>

/* *+

The VERBOSE macro is no longer used, instead
export CPL_MSG_LEVEL=info
can be used at run-time.
(At compile time, the default message level can be
controlled via cpl_test_init()).

#define VERBOSE
+* */

/*
 * Test for functions returning a generic pointer to data.
 *
 * r = variable to store returned pointer to data - expected non-NULL
 * f = function call
 * m = message
 */



#define test_data(r,f,m, flag)                        \
    do {                                        \
        cpl_msg_info("test_data", "%s", m);     \
        r = f;                                  \
        cpl_assert(r != NULL);                  \
        if (cpl_error_get_code() != CPL_ERROR_NONE) flag=EXIT_FAILURE; \
        cpl_test_error(CPL_ERROR_NONE);         \
    } while (0)                                 \


/*
 * Test for functions returning 0 on success.
 *
 * f = function call
 * m = message
 */

#define test(f,m, flag)                               \
     do {                                       \
         cpl_msg_info("test", "%s", m);         \
         cpl_test_zero(f);                      \
         if (cpl_error_get_code() != CPL_ERROR_NONE) flag=EXIT_FAILURE; \
         cpl_test_error(CPL_ERROR_NONE);        \
     } while (0)

/*
 * Test for expected failure in functions returning 0 on success.
 *
 * e = expected error code
 * f = function call
 * m = message
 */

#define test_failure(e,f,m, flag)                     \
     do {                                       \
         cpl_msg_info("test_failure", "%s", m); \
         if (cpl_error_get_code() != e) flag=EXIT_FAILURE; \
		 cpl_test_eq_error(f, e);               \
     } while (0)

/*
 * Test for expected failure in functions returning a pointer on success.
 *
 * e = expected error code
 * f = function call
 * m = message
 */

#define test_pfailure(e,f,m,flag)                     \
     do {                                       \
         cpl_msg_info("test_failure", "%s", m); \
         f; \
         cpl_error_code err = cpl_error_get_code (); \
         cpl_test_eq_error(err, e);               \
         if (err != e) flag=EXIT_FAILURE;			\
     } while (0)


/*
 * Test for functions returning an expected integer value.
 *
 * e = expected value
 * f = function call
 * m = message
 */

#define test_ivalue(e,f,m, flag)                      \
     do {                                       \
         cpl_msg_info("test_ivalue", "%s", m);  \
         cpl_test_eq(f, e);                     \
         if (f != e) flag = EXIT_FAILURE; 		\
        cpl_test_error(CPL_ERROR_NONE);        \
     } while (0)

/*
 * Test for functions returning an expected pointer value.
 *
 * e = expected value
 * f = function call
 * m = message
 */

#define test_pvalue(e,f,m, flag)                      \
     do {                                       \
         cpl_msg_info("test_pvalue", "%s", m);  \
         cpl_test_eq_ptr(f, e);                 \
         if (cpl_error_get_code() != CPL_ERROR_NONE) flag = EXIT_FAILURE; 	\
         cpl_test_error(CPL_ERROR_NONE);        \
     } while (0)

/*
 * Test for functions returning an expected floating point value.
 *
 * e = expected value
 * t = tolerance on expected value
 * f = function call
 * m = message
 */

#define test_fvalue(e,t,f,m, flag)                    \
     do {                                       \
         cpl_msg_info("test_fvalue", "%s", m);  \
         cpl_test_abs(f, e, t);                 \
         if (cpl_error_get_code() != CPL_ERROR_NONE) flag = EXIT_FAILURE; 	\
         cpl_test_error(CPL_ERROR_NONE);        \
     } while (0)

/*
 * Test for functions returning an expected complex value.
 *
 * e = expected value
 * t = tolerance on expected value
 * f = function call
 * m = message
 */

#define test_cvalue(e,t,f,m, flag)                    \
     do {                                       \
         cpl_msg_info("test_cvalue", "%s", m);  \
         cpl_test_abs_complex(f, e, t);         \
         if (cpl_error_get_code() != CPL_ERROR_NONE) flag = EXIT_FAILURE; 	\
         cpl_test_error(CPL_ERROR_NONE);        \
     } while (0)

/*
 * Test for functions returning an expected character string.
 *
 * e = expected value
 * f = function call
 * m = message
 */

#define test_svalue(e,f,m, flag)                      \
     do {                                       \
         cpl_msg_info("test_svalue", "%s", m);  \
         cpl_test_eq_string(f, e);              \
         if (cpl_error_get_code() != CPL_ERROR_NONE) flag = EXIT_FAILURE; 	\
         cpl_test_error(CPL_ERROR_NONE);        \
     } while (0)

#define GRAVI_INS_TRAIN_EXT_NAME "INS_TRAIN"
#define GRAVI_INS_DESC_EXT_NAME "INS_DESCRIPTION"
#define GRAVI_PRIMARY_HDR_EXT "PRIMARY_HDR"
#define GRAVI_IMAGING_DATA_NAME_EXT "IMAGING_DATA"
#define GRAVI_IMAGING_DETECTOR_NAME_EXT "IMAGING_DETECTOR"

int gravi_array_compare(cpl_array * , cpl_array * );
int gravi_table_compare(cpl_table *, cpl_table *);
int gravi_propertylist_compare(cpl_propertylist *, cpl_propertylist *);
int gravi_data_compare(gravi_data *, gravi_data *);


#endif

