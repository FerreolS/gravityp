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
 * gravi_utils-test.c
 *
 *  Created on: 16 août 2011
 *      Author: nabih
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <cpl_test.h>
#include <cpl.h>
#include "gravi_strehl.h"

//Prototypes
int gravi_strehl_tests(void);
int gravi_strehl_test_1();
int gravi_strehl_test_2();


int gravi_strehl_test_1()
{
    gravi_strehl(1000.);
   
    return EXIT_SUCCESS;
}

int gravi_strehl_test_2()
{
    if (gravi_strehl(-1.) == CPL_ERROR_NONE)
        return EXIT_FAILURE;
   
    return EXIT_SUCCESS;
}

int gravi_strehl_tests(void)
{
    
    gravi_strehl_test_1();

    gravi_strehl_test_2();

    return EXIT_SUCCESS;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Unit tests of gravi_utils module
 */
/*----------------------------------------------------------------------------*/

int main(void)
{
    int flag;
    cpl_test_init(PACKAGE_BUGREPORT, CPL_MSG_INFO);

    flag=gravi_strehl_tests();

    if (flag == EXIT_FAILURE)
    {
    	cpl_test_end(0);
    	exit(EXIT_FAILURE);
    }
    cpl_test_end(0);
    exit(EXIT_SUCCESS);
}
