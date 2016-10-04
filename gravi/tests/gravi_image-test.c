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
    const char *const fctid = "test_image";
    const char *const test_subject = "gravi_image";
    cpl_errorstate prestate = cpl_errorstate_get();
    cpl_image * image;

    /* Test with invalid input */
   image=gravi_image(NULL, NULL);

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
        /* --isotropic */
        p = cpl_parameter_new_value("gravi.gravi_image_reconstruct.isotropic_option",
                CPL_TYPE_BOOL, "a flag", "gravi.gravi_image_reconstruct", FALSE);
        cpl_parameterlist_append(parlist, p);

        /* --pixelsize */
        p = cpl_parameter_new_value("gravi.gravi_image_reconstruct.pixelsize",
                CPL_TYPE_DOUBLE, "size of the pixel (milliarcseconds)",
                "gravi.gravi_image_reconstruct", 0.2);
        cpl_parameterlist_append(parlist, p);

        /* --dim */
        p = cpl_parameter_new_value("gravi.gravi_image_reconstruct.dim",
                CPL_TYPE_INT, "number of pixels per side of the image",
                "gravi.gravi_image_reconstruct", 100);
        cpl_parameterlist_append(parlist, p);

        /* --regul */
        p = cpl_parameter_new_value("gravi.gravi_image_reconstruct.regul",
                CPL_TYPE_STRING, "name of regularization method",
                "gravi.gravi_image_reconstruct", "totvar");
        cpl_parameterlist_append(parlist, p);

        /* --regul_mu */
        p = cpl_parameter_new_value("gravi.gravi_image_reconstruct.regul_mu",
                CPL_TYPE_DOUBLE, "global regularization weight",
                "gravi.gravi_image_reconstruct", 1E4);
        cpl_parameterlist_append(parlist, p);

        /* --maxeval */
        p = cpl_parameter_new_value("gravi.gravi_image_reconstruct.maxeval",
                CPL_TYPE_INT, "maximum number of evaluations of the objective function",
                "gravi.gravi_image_reconstruct", 2000);
        cpl_parameterlist_append(parlist, p);

	    /* --timeout */
	    p = cpl_parameter_new_value("gravi.gravi_image_reconstruct.timeout",
            	CPL_TYPE_DOUBLE, "Maximum execution time of Mira process (s)",
            	"gravi.gravi_image_reconstruct", 60.);
	    cpl_parameterlist_append(parlist, p);


        /* Call the function */
        image=gravi_image(frame, parlist);
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
        
    }
    
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
    cpl_test_init(PACKAGE_BUGREPORT, CPL_MSG_WARNING);
#else
    cpl_init();
#endif

    test_image();

    cpl_test_end(0);
    exit(EXIT_SUCCESS);
}

/**@}*/
