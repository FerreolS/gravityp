#ifndef TIMS_H
#define TIMS_H
/*******************************************************************************
* E.S.O. - VLT project
*
* tims.h
*
* "@(#) $Id: tims.h 30977 1999-11-23 08:28:45Z bgilli $"
*
* who        when      what
* ---------  --------  ----------------------------------------------
* B.GILLI    29/04/93  Preliminary version
* B.GILLI    27/11/93  Elaborated version for ccs design document.
* B.GILLI    08/12/93  modified timsSleep in timerSleep, and timerSendMessage
*                      in timerSendCommand and timerSendReply
* B.GILLI    23/09/94  moved #include "msg.h" into timer.h file.
* MCOMIN     21/11/94  Add definitions for error mnemonics
* B.GUSTAFSSON 29/11/94 Including specific CCS/LCC include files
* B.GUSTAFSSON 01/12/94 Added C++ support
* B.GILLI    03/07/95  added prototypes of JD conversion subs.
* M.GAJARDO  15/11/96  timsTimeToIsoString/timsIsoStringtoTime prototypes added
* M.GAJARDO  03/12/96  timsMILLISECOND and timsMICROSECOND defined
* M.GAJARDO  15/12/97  Added timsIsoStringToDateTime, timsDateTimeToUTC
* T.EBERT    09/12/98  Added check for NOT_Y2K
* B.GILLI    15/06/99  fully removed non ISO compliant functions.
* B.GILLI    19/11/99  changed time parm name to avoid shadowing pb.
*/


/************************************************************************
 * CCS/TIME SYSTEM - Interface File
 */

#ifdef __cplusplus
extern "C" {
#endif

/* 
 * header files
 */
#include "ccs.h"
#include "time.h"


/*
 * Message constants
 */

#define timsMODE_UTC_STRING	   "UTC"
#define timsMODE_LOCAL_STRING	   "Local"
#define timsMODE_SIMULATION_STRING "Simulation"


/* 
 * data types
 */
typedef enum {
    timsDATE = 1,              /* yy-mm-dd  The date at midnight at the    */ 
                               /*           beginning of that day.         */
    timsTIME_OF_DAY,           /* hh:mm:ss.uuuuuu  The time since midnight */
                               /*                  the same day.           */
    timsABS_TIME               /* yy-mm-dd hh:mm:ss.uuuuuu  The combination*/
                               /*                  of the two previous one.*/
} timsTIME_FORMAT;


typedef enum {                  /* Working mode for the time system    */
    timsUTC = 1,                /* WS is synchronized with time Server */
    timsLOCAL,                  /* WS is not synchronized              */
    timsSIMULATION              /* For compatibility with LCU          */
} timsMODE;
/* 
 * -------- end of data type definitions 
 */
/* 
 * macros 
 */
#define timsMILLISECOND  3
#define timsMICROSECOND  6 

/* 
 * functions prototypes 
 */

ccsCOMPL_STAT timsGetUTC(        /*returns UTC in internal representation.  */
   ccsTIMEVAL         *UTC,      /* returned value       */
   timsMODE           *mode,     /* Current mode         */
   ccsERROR           *error     /* error structure      */
   );

ccsCOMPL_STAT timsSetUTC(        /*Set UTC  Only if needed for simulation*/
   const ccsTIMEVAL   *UTC,      /* time to be set        */
   ccsERROR           *error     /* error structure      */
   );

ccsCOMPL_STAT timsSetMode(        /*Set operating mode for time system  */
   timsMODE           mode,       /* mode to be set        */
   ccsERROR           *error     /* error structure      */
   );

ccsCOMPL_STAT timsGetMode(        /*Get operating mode for time system  */
   timsMODE           *mode,      /* Current mode        */
   ccsERROR           *error     /* error structure      */
   );
ccsCOMPL_STAT timsAddTime(       /* Add two time values in internal format.  */
   const ccsTIMEVAL         *time1,    
   const ccsTIMEVAL         *time2,
   ccsTIMEVAL               *time3,     /* time = time1 + time2    */
   ccsERROR           *error
   );
ccsCOMPL_STAT timsSubTime(       /* Subtract time values in internal format.  */
   const ccsTIMEVAL         *time1,    
   const ccsTIMEVAL         *time2,
   ccsTIMEVAL               *time3,     /* time = time1 - time2    */
   ccsERROR           *error
   );

ccsCOMPL_STAT timsTimeToIsoString(  /* Convert an internal representation into
                                    a printable one  */
   const ccsTIMEVAL   *time1,     /* Time in internal format */
   timsTIME_FORMAT    format,    /* Format to be used as output */
   vltUINT8           precision, /* presision on fraction of seconds [0-6]*/
   char               *string,   /* String to put result */
   ccsERROR           *error
   );
ccsCOMPL_STAT timsIsoStringToTime(  /* Convert an external representation into
                                    an internal one  */
   const char         *string,   /* String containing ASCII time */
   timsTIME_FORMAT    format,    /* Format of that string   */
   ccsTIMEVAL         *time1,     /* Time in internal format */
   ccsERROR           *error
   );

ccsCOMPL_STAT timsIsoStringToDateTime(/* Split IsoSting into Time and Date */
    const char       *string,       /* String containing ASCII time     */
    timsTIME_FORMAT  format,        /* Format of that string            */
    char             *date,          /* String containing date */
    char             *timestr,       /* String containing date */
    ccsERROR         *error         /* Returned error structure         */
    );
ccsCOMPL_STAT timsDateTimeToUTC(     /* Concatenate date+time ISO format */
                                    /* into internal representation     */
    const char       *date,         /* String containing date     */
    const char       *timestr,         /* String containing time     */
    ccsTIMEVAL       *timeUTC,         /* Time in internal format          */
    ccsERROR         *error         /* Returned error structure         */
    );
 


/* Convert Universal Time into Julian Date */
vltDOUBLE timsUTCToJD(ccsTIMEVAL *ut);    /* time in UT                      */
void timsJDToUTC( vltDOUBLE jd, ccsTIMEVAL *utc ); /*time in JD to UTC mode  */

/* Convert Universal Time into Modified Julian Date */
vltDOUBLE timsUTCToMJD(ccsTIMEVAL *ut);   /* time in UT                      */
void timsMJDToUTC( vltDOUBLE mjd, ccsTIMEVAL *utc ); /*time in MJD to UTC    */

/* 
 * message definitions
 */

/* 
 * error definitions
 */


#ifdef MAKE_VXWORKS
#include "timsLCC.h"
#else
#include "timsCCS.h"
#endif

#ifdef __cplusplus
}
#endif

#endif /*!TIMS_H*/
