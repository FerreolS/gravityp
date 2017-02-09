/*******************************************************************************
* E.S.O. - VLT project
*
* tims.h
*
* "@(#) $Id: timsCCS.h 50 1995-02-02 11:37:06Z therlin $"
*
* who        when      what
* ---------  --------  ----------------------------------------------
* B.GILLI    29/04/93  Preliminary version
* B.GILLI    27/11/93  Elaborated version for ccs design document.
* B.GILLI    08/12/93  modified timsSleep in timerSleep, and timerSendMessage
*                      in timerSendCommand and timerSendReply
* B.GILLI    23/09/94  moved #include "msg.h" into timer.h file.
* MCOMIN     21/11/94  Add definitions for error mnemonics
* B.GUSTAFSSON 28/11/94 Created CCS specific include file
*/


/************************************************************************
 * CCS/TIME SYSTEM - Interface File
 */


/*
 * Define error mnemonics 
 */

#include "timsErrors.h"

/* 
 * constants 
 */
#define TIMS_POINT_NAME  "<alias>xntpd"
#define TIMS_MODE_ATTR   "timsmode"
#define TIMS_OFFSET_ATTR "timsoffset"


/* 
 * data types
 */

/* 
 * -------- end of data type definitions 
 */
/* 
 * macros 
 */

/* 
 * functions prototypes 
 */

ccsCOMPL_STAT timsUTCToST(       /* Convert UTC into ST  */
   ccsTIMEVAL         *UTC,      /* UTC to be converted, if NULL get it   */
   vltDOUBLE          *ST,       /* returned ST, other format ??? */
   ccsERROR           *error
   );

/* 
 * message definitions
 */


/* 
 * error definitions
 */

