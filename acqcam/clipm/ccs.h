/*************************************************************************
* E.S.O. - VLT project
#
# "@(#) $Id: ccs.h 178833 2009-01-16 14:36:11Z mcomin $" 
*
* ccs.h  -  CCS/Common Definitions - Interface File
*
* who        when      what
* ---------  --------  ----------------------------------------------
* B.GILLI    29/04/93  Preliminary version
* B.GILLI    13/05/93  Added ccsTIMEVAL type and TBD type (just for compiler)
* M.COMIN    27/05/93  Added ccsENV_TYPE data type. 
* M.COMIN    01/07/93  Complete the list of data types.
*                      All data types have prefix 'vlt'
*                      Correct definition for enum. ccsENV_TYPE
*
* M.COMIN    05/07/93  Changed type of tv_sec in ccsTIMEVAL, from 'long'
*                      to 'unsigned long'.
*
* M.COMIN    14/07/93  Typedef of DB_XREF moved to db.h.
*                      
* A.WALLANDER19/07/93  Correct prototypes for ccsInit and ccsExit.
*
* M.COMIN    20/07/93  Changed error structure : ccsERROR.
*                      Moved typedef logLOC_ID from log.h to here and 
*                      renamed it as ccsLOG_ID.
*                      Added define ccsLOCID_LEN.
*                      Add Task Option data structures for UNIX and VxWorks.
*
* M.COMIN    10/08/93  Add routine ccsGetMyProcId()
*                      Changed declaration of taskname in ccsVXM_TASK_OPTIONS
*
* M.COMIN    19/08/93  Add definition for TRUE/FALSE.
*                      Changed error structure : ccsERROR.
*
* A.WALLANDER20/08/93  Changed declaration of run option.
*
* A.WALLANDER26/08/93  VLT-SBI-MAIL-026. run option union,
*                                        obituary constants
*                                        new type for stack id
*
* M.COMIN     1/10/93  Add file of error definitions.
*
* M.COMIN    07/10/93  Add definition for modules : UNIX, VXWORKS, RTAP
*
*
*---------------------
*
* M.COMIN    25/01/94  Modifications according to the design 1.0
*                      - Change CCS module names according to PS
*                      - Add define ccsLOCAL_ENV
*
* M.COMIN    16/02/94  Add ccsDebugPrintf() and ccsGetInstallPath()
* M.COMIN    17/02/94  Add ccsOpenFile() 
* M.COMIN    09/05/94  Complete list of CCS module names and error offsets
* M.COMIN    29/09/94  Add definitions for 'evt' module
*
* M.COMIN    02/11/94  Add function ccsGetProcName 
* G.CHIOZZI  07/11/94  Added prototype for ccsGetTimeOfDay()
* M.COMIN    21/11/94  Error offsets and definitions moved to the module
*                      they belong to. 
* B.GUSTAFSSON 29/11/94 Including specific CCS/LCC include files
* B.GUSTAFSSON 01/12/94 Added C++ support
* ---------------------
* M.COMIN      26/03/96  Modification of ccsERROR structure
*              29/04/96  Modification of ccsERROR structure to allows multiple
*                        error stacks.
* bgustafs     03/05/96  modified error stack
* bgilli       13/05/96  added ccsCCS_VERSION for ext. compilation control
* bgustafs     14/05/96  modified error stack size on LCU
* bgilli       14/05/96  redefined ccsCCS_VERSION as 9606
* mcomin       15/05/96  Use different  stack sizes for CCS and CCS_LITE
* mcomin       04/02/97  Add define ccsSAMP
* bgilli       17/03/97  redefined ccsCCS_VERSION as 9705
* mcomin       18/03/97  Add CCS standard states
*              05/02/98  Definition of standard states moved into ccsStates.h
*              27/02/98  Add define for vltSCALAR_MAX_SIZE
*              24/06/98  Merging for CCS_LITE
* tebert       30/06/98  Changed size of error stack to default for CCS LITE
* mcomin       28/07/98  Changed typedef of vltINT8 for PPC porting
*              30/07/98  Changed indentation for #define (cc68k does not compile)
* tebert       10/08/98  Added new type of environment for CCS LITE (WSL)
* bgilli       18/08/98  redefined ccsCCS_VERSION as 199810     
* bgilli       27/08/98  added stuff for ccsWaitForProcStat.
* tebert       31/05/99  redefined ccsCCS_VERSION to  199910
* mcomin       13/12/99  redefined ccsCCS_VERSION to  200002
* bgilli       25/04/00  redefined ccsCCS_VERSION to  200011, removed book
* tebert       18/07/00  added ccsHIS
* tebert       06/09/01  ccsSEND_OBITUARY/ccsRECV_OBITUARY modified to cope
*                        with CCS LITE
* bgilli       02/12/03  added macro __FILE_LINE__ taken from eccs
* bgilli       21/09/04  moved general interest macros to vltPortGeneral.h
****************************************************************************/

#ifndef CCS_H
#define CCS_H

#ifdef __cplusplus
extern "C" {
#endif

#include "vltPortGeneral.h"
#include "ccsStates.h"

/************************************************************************
 *                           CCS  Constants                             *
 ************************************************************************/

#define ccsCCS_VERSION      200011
 
#define ccsDBITEMNAME_LEN   255   /* NOT USED : Only for compatibility  */
                                  /* with LCC                           */

#define ccsENVNAME_LEN        7   /* max. length of an environment name */ 
#define ccsPROCNAME_LEN      19   /* max. length of a process name      */
#define ccsMODULEID_LEN       7   /* max. length of a module name       */
                                  /* 6 characters + 1 byte alignement   */

#define ccsCMD_LEN            7   /* max. length of a command            */
#define ccsLOCID_LEN         80   /* max. length of a location specification */


#define ccsFALSE              0   /* Definitions according to VxWorks rules */
#define ccsTRUE               1           

#define ccsSEND_OBITUARY    0x04  /* send obituary on termination         */
#define ccsRECV_OBITUARY    0x01  /* receive obituary from other proc/env */


#define ccsMAX_DE_TO_STR_CHARS ((256 * 4) + 1)
#define vltSCALAR_MAX_SIZE  256   /* Defines the maximum size for a data element */

/*
 *   Define module names for RTAP, UNIX and VxWORKS
 */ 

#define  modRTAP        "RTAP"
#define  modUNIX        "UNIX"
#define modVXWORKS	"VXWORKS"

/*
 *   Define CCS module names
 */ 

#define   ccsGEN    "ccs"
#define   ccsDB     "db"
#define   ccsMSG    "msg" 
#define   ccsERR    "err"
  /* #define   ccsBOOK   "book" */
#define   ccsALRM   "alrm"
#define   ccsSCAN   "scan"
#define   ccsTIMS   "tims"
#define   ccsUIF    "uif"
#define   ccsLOG    "log"
#define   ccsEVT    "evt"
#define   ccsCMD    "cmd"
#define   ccsSAMP   "samp"
#define   ccsHIS    "his"

/*  Define Environment : Empty string means local environment */

#define ccsLOCAL_ENV  "" 

/* The macro __FILE_LINE__ produces a string in the format: "FileName:Line" */
/* with the name and line of the file where it has been called.             */
/* It is usefull in error messages as locId                                         */
/* (the other macros defined here are just for support, since I have not    */
/* found a cleaner way to define __FILE_LINE__ macro                        */
/* Copied from eccs.h                                                       */

#define _ccs_tostr(a) #a
#define _ccs_tostr_pass2(a) _ccs_tostr(a) 

#ifndef __FILE_LINE__ 
     #define __FILE_LINE__   __FILE__ ":" _ccs_tostr_pass2(__LINE__)
#endif






/*  CCS time structure */

typedef vltTIMEVAL ccsTIMEVAL;



/*
 * Direct address specification
 */

typedef struct {                    /* Direct point specification        */
     vltUINT16         plin;        /* Point internal number (rtPlin)    */
     vltUINT8          ain;         /* Attribute internal number (rtAin) */
     vltUINT8          reserved;
     vltUINT16         startRec;    /* First record to be accessed       */
     vltUINT16         endRec;      /* Last record to be accessed        */
     vltUINT8          startField;  /* First field to be accessed        */
     vltUINT8          endField;    /* Last field to be accessed         */
} ccsPOINT_ID;


typedef struct {
    vltUINT16         nameType;
    ccsPOINT_ID       pointId;
} vltDBXREF;

typedef char ccsENVNAME[ccsENVNAME_LEN+1];        /* Environment(Node) name */
typedef char ccsPROCNAME[ccsPROCNAME_LEN+1];      /* Process name           */
typedef unsigned char ccsPROCNUM;                 /* Process number         */
typedef char ccsMODULEID[ccsMODULEID_LEN+1];      /* Software module name   */
typedef char ccsLOC_ID[ccsLOCID_LEN+1];           /* Subroutine identifier  */

typedef char ccsDBITEMNAME[ccsDBITEMNAME_LEN+1];  /* Name of database item  */
                                                     /* Includes:           */
                                                     /* - environment name  */ 
                                                     /* - point name        */
                                                     /* - attribute name    */
                                                  /* NOT used by WS         */
/*
 *    Error Structure Definition.  It is formed by two parts
 *    -  Information on the stack location    
 *    -  Routine dependent error parameters 
 */

typedef union {
    struct {
        vltUINT16   hostId;         /* host identification: 16 lower bits   */
                                    /* of internet address                  */
        vltUINT16   localNumber;    /* local stack sequence number          */
    } local;
    vltUINT32   id;                 /* Complete stack id                    */
} ccsSTACK_ID;


typedef struct {                    

    vltSTRING32   timeStamp;      /* Time when the error occurred         */
    ccsENVNAME   envName;         /* Where to log errors for this stack   */
    ccsSTACK_ID  stackId;         /* on which stack                       */
    vltUINT16    sequenceNumber;  /* Sequence number in the stack         */
 
    ccsPROCNUM   procNum;         /* Process number                       */
    ccsPROCNAME  procName;        /* Name of the process                  */
    ccsLOC_ID    location;        /* Location where the error occurred :  */
				  /* - Subroutine, Function, etc ..       */ 
    ccsMODULEID  moduleId;        /* Addresses the error messages file    */
    vltINT16     errorNumber;     /* Addresses the right message          */
    vltINT8      severity;
    vltSTRING256  runTimePar;     /* Run Time Information about the error */
    vltLOGICAL   errStackEmpty;   /* Private flag used by LCC : for the   */
                                  /* normal user the stack is empty  only */
                                  /* if the sequence number = 0           */
} ccsSTACK_ELEM;                 

#ifdef MAKE_VXWORKS
#   define errSTACK_SIZE 10
#else
#   define errSTACK_SIZE 18
#endif



typedef struct {                            /* Error Stack Structure         */
    ccsSTACK_ELEM  errStack[errSTACK_SIZE]; /* Time when the error occurred  */
    vltINT16       errStackSize;            /* Current Stack Size            */
    vltLOGICAL     errStackOverflow;   /* Flags if a stack overflow occurred */
    vltLOGICAL     errReplyFlag;       /* Identifies whether the current     */
                                       /* stack comes from an error reply    */
} ccsERROR_STACK;                         


typedef struct {                    /* Returned error structure             */

    vltSTRING32      timeStamp;     /* Time when the error occurred         */
    ccsENVNAME      envName;        /* Where to log errors for this stack   */
    ccsSTACK_ID     stackId;        /* on which stack                       */
    vltUINT16       sequenceNumber; /* Sequence number in the stack         */
    ccsPROCNUM      procNum;        /* Process number                       */
    ccsPROCNAME     procName;       /* Name of the process                  */
    ccsLOC_ID       location;       /* Location where the error occurred :  */
				    /* - Subroutine, Function, etc ..       */ 
    ccsMODULEID     moduleId;       /* Addresses the error messages file    */
    vltINT16        errorNumber;    /* Addresses the right message          */
    vltINT8         severity;
    vltSTRING256     runTimePar;    /* Run Time Information about the error */
    vltLOGICAL      errStackEmpty;  /* Indicates that the stack is empty */
    ccsERROR_STACK  stack;          /* Stack associated to current error    */
} ccsERROR;                         


/*
 *   Definition of the routine completion status
 */
typedef enum {
    FAILURE = 1,
    SUCCESS 
} ccsCOMPL_STAT;       /* Completion status returned by subroutines */


typedef enum {
    ENV_UNKNOWN = 1,   /* Environment specification not recognized */
    WS,                /* Environment type = Workstation           */ 
    LCU,              /* Environment type = LCU                   */
    WSL               /* Environment type = Workstation (CCS-LITE) */
} ccsENV_TYPE;         /* List of Environment Type */


/*
 *   Data structure run options to schedule a process.   
 *   Depending if destination environment is unix or VxWorks a 
 *   union is used.
 */
typedef struct {                  /*  VxWorks Task options   */
    vltUINT32  priority;          /*  Task priority          */
    vltUINT32  options;           /*  Task options           */
    vltUINT32  stacksize;         /*  Stack size             */
    char       *taskname;         /*  Pointer to task name   */
} ccsVXWORKS_OPTIONS;

typedef struct {                  /* Unix Process Options    */
    vltINT32   flags;
} ccsUNIX_OPTIONS;  

typedef union {                   /* Combined definition     */  
    ccsVXWORKS_OPTIONS  vxto;
    ccsUNIX_OPTIONS     unixto;
} ccsRUNOPT;

typedef enum {
    ccsPROC_EXIT = 1,   /* Process not existing */
    ccsPROC_REACHABLE,  /* PING can be dropped into inpu queue */ 
    ccsPROC_ANSWERING   /* PING with reply successfull */
} ccsPROCSTAT;          /* Possible states for a process */


/************************************************************************
 *                     CCS Functions prototypes                         *
 ************************************************************************/

/*  Application must call the following two routines at the beginning   */
/*  and at the end.                                                     */

ccsCOMPL_STAT ccsInit (            /* Initialize the CCS environment        */
   const ccsPROCNAME procName,     /* process name under which to register  */
   vltUINT32         obituary,     /* flag for obituary handling (future)   */
   void      (*breakHandler)(int), /* Break handler, NULL = default handler */
   void      (*killHandler)(int),  /* Kill handler, NULL = default handler  */
   ccsERROR          *error 
    );

ccsCOMPL_STAT ccsExit (        /* Release CCS resources used by application */
    ccsERROR    *error
    );


ccsCOMPL_STAT ccsGetMyProcId(  /* Get information about 'my process' */
    ccsENVNAME  envName,       /* Environment name     */
    ccsPROCNUM  *procNum,      /* Process Id number    */
    ccsPROCNAME procName,      /* Process name         */
    ccsERROR    *error         /* Returned data structure */
    );

ccsCOMPL_STAT ccsWaitForProcStat( /* wait for process in given state */
    ccsENVNAME         envName,
    ccsPROCNAME        procName,
    ccsPROCSTAT        procStat,   /* state to wait for */
    vltINT32           howLong,    /* How long to wait in seconds */
    ccsERROR          *error
    );
#ifdef MAKE_VXWORKS
#include "ccsLCC.h"
#else
#include "ccsCCS.h"
#ifdef CCS_LIGHT
#include "ccsLite.h"
#else
#define ENVNAME "RTAPENV"
#endif
#endif

#ifdef __cplusplus
}
#endif

#endif /*!CCS_H*/
