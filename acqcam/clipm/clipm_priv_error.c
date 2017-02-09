/*********************************************************************
 * E.S.O. - VLT project
 *
 * "@(#) $Id: clipm_priv_error.c 169734 2008-07-02 16:20:41Z hlorch $"
 *
 * CLIPM messaging
 *
 * who       when        what
 * --------  ----------  ----------------------------------------------
 * hlorch    2008-02-29  created
 */

/*-----------------------------------------------------------------------------
    Includes
 -----------------------------------------------------------------------------*/

#include "clipm_priv_error.h"

#include <string.h>

/*-----------------------------------------------------------------------------
    Global Definitions
 -----------------------------------------------------------------------------*/

/* Don't document these values with doxygen. */
const char  _CLIPM_MSG_ERR_UNEXPECTED[] =
                "unexpected error, aborting."
                " Please report to the CLIP team.",
            _CLIPM_MSG_ERR_HANDLING[] =
                "error handling bug found, aborting."
                " Please report to the CLIP team.",
            CLIPM_MSG_ERR_2ROWXY[] =
                "must contain 2 rows (x,y)",
            CLIPM_MSG_ERR_DIFFSIZES[] =
                "different sizes",
            CLIPM_MSG_ERR_DIFFTYPES[] =
                "different types";

/*-----------------------------------------------------------------------------
    Implementation
 -----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
 * @ingroup  clipm_priv_error
 * @internal
 * @brief   Print two messages into one string. If both messages are not empty,
 *          then they are separated by a colon.
 * @param   outstr  String to print into,
 *                  can be NULL
 * @param   msg1    Input string 1
 * @param   msg2    Input string 2
 * @param   maxlen  Maximum no. of characters to print, must fulfill
 *                  @a maxlen <= sizeof(@a outstr) - 1
 * @return  Nothing
 */
/*----------------------------------------------------------------------------*/
void        _clipm_priv_error_sprint_messages(
                                            char        *outstr,
                                            const char  *msg1,
                                            const char  *msg2,
                                            int         maxlen)
{
    if (outstr == NULL)
        return;
    
    outstr[0] = '\0';
    if (msg1 == NULL || msg1[0] == '\0')
    {
        if (msg2 != NULL)
        {
            strncpy(outstr, msg2, maxlen);
            outstr[maxlen] = '\0';
        }
    }
    else
    {
        strncpy(outstr, msg1, maxlen);
        outstr[maxlen] = '\0';
        
        if (msg2 != NULL && msg2[0] != '\0')
        {
            int l;
            l = strlen(outstr);
            strncat(outstr, ": ", maxlen-l);
            l += 2;
            l = l > maxlen ? maxlen : 1;
            strncat(outstr, msg2, maxlen - l);
        }
    }
    
    return;
}



