/* $Id: gravi_dfs-test.c,v 1.3 2007/07/30 07:08:14 llundin Exp $
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

/*
 * $Author: llundin $
 * $Date: 2007/07/30 07:08:14 $
 * $Revision: 1.3 $
 * $Name: HEAD $
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/*-----------------------------------------------------------------------------
                                Includes
 -----------------------------------------------------------------------------*/

#include <stdlib.h>
#include <string.h>

#include <cpl.h>

#include <gravi_image.h>
#include <gravi_dfs.h>

/*----------------------------------------------------------------------------*/
/**
 * @defgroup gravi_image_test  Unit test of gravi_image
 *
 */
/*----------------------------------------------------------------------------*/

/**@{*/

/*----------------------------------------------------------------------------*/
/**
  @brief    Textual representation of CPL frame group
  @param    group     to convert
  @return   textual representation
 */
/*----------------------------------------------------------------------------*/
static const char *
frame_group_to_string(cpl_frame_group group)
{
    switch(group) {
    case CPL_FRAME_GROUP_RAW:
        return CPL_FRAME_GROUP_RAW_ID;
        break;
    case CPL_FRAME_GROUP_NONE:
        return "NONE";
        break;
    case CPL_FRAME_GROUP_CALIB:
        return CPL_FRAME_GROUP_CALIB_ID;
        break;
    case CPL_FRAME_GROUP_PRODUCT:
        return CPL_FRAME_GROUP_PRODUCT_ID;
        break;
    default:
        return "???";
        break;
    }
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Unit test of gravi_image
 */
/*----------------------------------------------------------------------------*/
static void test_image(void)
{
#ifdef YORICK_BIN
    const char *const fctid = "test_image";
    const char *const test_subject = "gravi_image";
    cpl_errorstate prestate = cpl_errorstate_get();
    cpl_image * image;

    /* Test with invalid input */
   image=gravi_image(NULL, NULL, "");

    if (cpl_errorstate_is_equal(prestate)) {
    	cpl_msg_error(fctid, "Function %s did not fail on NULL input",
    	                      test_subject);
		cpl_end();
		exit(EXIT_FAILURE);
    }

    cpl_errorstate_set(prestate);

    /* Test with valid input */
    {
        cpl_frame *frame = cpl_frame_new();
        cpl_parameter *p;
        cpl_parameterlist *parlist;

        cpl_frame_set_filename(frame, DATADIR "data1.oifits");
        cpl_frame_set_tag(frame, GRAVI_MIRA_INPUT_PROCATG);

        parlist=cpl_parameterlist_new();
        gravi_parameter_add_image (parlist);

        /* Call the function */
        image=gravi_image(frame, parlist, "FKV0497");
        if (!cpl_errorstate_is_equal(prestate)) {
        	cpl_msg_error(fctid, "Function %s failed",
        	                      test_subject);
    		cpl_end();
    		exit(EXIT_FAILURE);
        }

        /* Verify results */
        if (image == NULL) {
                cpl_msg_error(fctid, "Return null image");
                cpl_errorstate_dump(prestate, CPL_FALSE, NULL);
                cpl_end();
                exit(EXIT_FAILURE);       
            }

        cpl_image_delete(image);
        cpl_parameterlist_delete(parlist);
        cpl_frame_delete(frame);
        if (!cpl_errorstate_is_equal(prestate)) 
        {
            cpl_msg_error(fctid, "Function %s failed on valid input",
                          test_subject);
            cpl_end();
            exit(EXIT_FAILURE);
        }

    }
    
#endif
    return;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Unit tests of gravi_dfs module
 */
/*----------------------------------------------------------------------------*/

int main(void)
{
#if defined CPL_VERSION_CODE && CPL_VERSION_CODE >= CPL_VERSION(4, 0, 0)
    cpl_test_init(PACKAGE_BUGREPORT, CPL_MSG_INFO);
#else
    cpl_init();
#endif

    test_image();

    cpl_test_end(0);
    exit(EXIT_SUCCESS);
}

/**@}*/
