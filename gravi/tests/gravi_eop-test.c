/* 
 * This file is part of the GRAVITY Pipeline
 * Copyright (C) 2016 European Southern Observatory
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 */


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/*-----------------------------------------------------------------------------
                                   Includes
 -----------------------------------------------------------------------------*/

#define _XOPEN_SOURCE 500
#include <stdlib.h>
#include <time.h>
#include "gravi_eop.h"

/*-----------------------------------------------------------------------------
                                   Static functions
 -----------------------------------------------------------------------------*/
static void gravi_eop_all_test(void);
static void gravi_eop_retrieve_eop_test(void);


/*-----------------------------------------------------------------------------
                                  Main
 -----------------------------------------------------------------------------*/
int main (void)
{

    cpl_test_init(PACKAGE_BUGREPORT, CPL_MSG_WARNING);

    gravi_eop_all_test();

    return cpl_test_end(0);
}

static void gravi_eop_all_test(void)
{

    //Tests for the online retrieval of the EOP data
    gravi_eop_retrieve_eop_test();
    
    return;
}

static void gravi_eop_retrieve_eop_test(void)
{

    //Initialize random number generator
    srandom(time(NULL));

    //Length of the raw retrieved data
    int data_length;

    //The FTP retrieval is tried a few times, in case there is a
    //transient network problem (ftp overloaded, for instance).
    int max_tries = 5;
    int itry = 1;
    //Get how many tests have failed so far
    int test_failed_before = cpl_test_get_failed();
    char * raw_text = NULL;
    while(itry < max_tries && raw_text == NULL)
    {
        cpl_msg_debug(__func__, "Trying EOP data retrieval. Trial %d",itry);
        //Test retrieving the EOP data from the server
        raw_text = gravity_eop_download_finals2000A(
                "ftp.eso.org",
                "/pub/dfs/pipelines/gravity/finals2000A.data",
                &data_length);
        if(raw_text == NULL)
        {
            itry++;
            //Reset the error status before trying again (with some random delay)
            if(itry < max_tries - 1)
            {
                int sleep_seconds = (int)(200.*random()/RAND_MAX);
                cpl_msg_debug(__func__, "Sleeping %d seconds before retrying",
                        sleep_seconds);
                sleep(sleep_seconds);
                errno = 0;
                cpl_error_reset();
            }
            else
                cpl_msg_debug(__func__, "Maximum number of retries (%d) "
                              "reached", max_tries);
        }
    }

    //A copy needs to be done since further cpl_test functions will set it to zero
    int errno_after = errno;

    //Check no cpl error is set
    cpl_test_error(CPL_ERROR_NONE);

    //Check there if the buffer is not empty
    cpl_test(data_length > 0);

    //Check pointer is not null
    cpl_test_nonnull(raw_text);

    //Check that no system calls have set errno
    cpl_test_zero(errno_after);

    //If those tests passed (no more failed tests than before),
    //then continue with further tests on the retrieved data
    if(cpl_test_get_failed() == test_failed_before)
    {

        //the first characters of raw_text is a number
        char * endptr;
        int num = strtol(raw_text, &endptr, 10);
        cpl_test(num > 0);

        //Test conversion to a table
        cpl_table * eop_table =
                gravity_eop_data_totable (raw_text, data_length);

        //Test that the table has more than one row
        cpl_test(cpl_table_get_nrow(eop_table) > 0);

        //Test that PMX, PMY, DUT don't have nonsense values (corrections are small)
        cpl_test(cpl_table_get_column_max(eop_table, "PMX") < 10);
        cpl_test(cpl_table_get_column_min(eop_table, "PMX") > -10);
        cpl_test(cpl_table_get_column_max(eop_table, "PMY") < 10);
        cpl_test(cpl_table_get_column_min(eop_table, "PMY") > -10);
        cpl_test(cpl_table_get_column_max(eop_table, "DUT") < 10);
        cpl_test(cpl_table_get_column_min(eop_table, "DUT") > -10);

        //Test that MJD increases monotonically
        for(cpl_size i_row = 1; i_row < cpl_table_get_nrow(eop_table); i_row++)
        {
            int null;
            cpl_test(cpl_table_get_double(eop_table, "MJD", i_row, &null) >
            cpl_table_get_double(eop_table, "MJD", i_row - 1, &null));
        }
        cpl_table_delete(eop_table);

        //Test conversion with NULL pointer
        eop_table = gravity_eop_data_totable (NULL, data_length);
        cpl_test_error(CPL_ERROR_NULL_INPUT);
        cpl_test_null(eop_table);

        //Test with wrong data length
        eop_table = gravity_eop_data_totable (raw_text, data_length - 1);
        cpl_test_error(CPL_ERROR_NULL_INPUT);
        cpl_test_null(eop_table);

        //Test with wrong data length and NULL pointer
        eop_table = gravity_eop_data_totable (NULL, data_length - 1);
        cpl_test_error(CPL_ERROR_NULL_INPUT);
        cpl_test_null(eop_table);

        free(raw_text);
    }

    //Test with wrong HOST name
    gravity_eop_download_finals2000A(
            "invalid_host.nowhere",
            "/products/eop/rapid/standard/finals2000A.data",
            &data_length);

    cpl_test_error(CPL_ERROR_DATA_NOT_FOUND);

    //Test with wrong URL
    gravity_eop_download_finals2000A(
            "ftp.eso.org",
            "/invalid/path",
            &data_length);

    cpl_test_error(CPL_ERROR_DATA_NOT_FOUND);

    //Test with wrong path
    gravity_eop_download_finals2000A(
            "ftp.eso.org",
            "WRONG FORMATTED PATH",
            &data_length);

    cpl_test_error(CPL_ERROR_DATA_NOT_FOUND);

    //Test with wrong URL and PATH
    gravity_eop_download_finals2000A(
            "invalid_host.nowhere",
            "/invalid/path",
            &data_length);

    cpl_test_error(CPL_ERROR_DATA_NOT_FOUND);

    //Test with NULL pointers
    gravity_eop_download_finals2000A(
            "ftp.eso.org",
            "/products/eop/rapid/standard/finals2000A.data",
            NULL);

    cpl_test_error(CPL_ERROR_NULL_INPUT);

    gravity_eop_download_finals2000A(
            NULL,
            "/products/eop/rapid/standard/finals2000A.data",
            &data_length);

    cpl_test_error(CPL_ERROR_NULL_INPUT);

    gravity_eop_download_finals2000A(
            "ftp.eso.org",
            NULL,
            &data_length);

    cpl_test_error(CPL_ERROR_NULL_INPUT);

    return;
}

