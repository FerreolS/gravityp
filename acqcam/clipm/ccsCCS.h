/*******************************************************************************
* E.S.O. - VLT project
#
# "@(#) $Id: ccsCCS.h 266574 2015-03-19 15:08:02Z tebert $" 
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
* M.COMIN    02/11/94  Add function ccsGetMyName 
* G.CHIOZZI  07/11/94  Added prototype for ccsGetTimeOfDay()
* M.COMIN    21/11/94  Error offsets and definitions moved to the module
*                      they belong to. 
* B.GUSTAFSSON 29/11/94 Created CCS specific include file
*
*---------------------
* M.COMIN    02/01/95   Add prototypes for Shared Memory Area
* M.COMIN    12/01/95   Add prototypes for ccsGetProcName()
* T.EBERT    13/06/95   Add prototype for ccsGetMyEnvName()
* M.COMIN    17/06/95   Add prototype for ccsFindFile()
* T.EBERT    20/10/95   Added definition of logErrOffSet
* J.I.Leon   06/11/95   Added definition of ccsGetProcNum
* B. Gilli   25/04/00   Removed bookErrOffset
*                  
****************************************************************************/

#include <stdio.h>

#include "ccsSHMemory.h"

/************************************************************************
 *                           CCS  Errors                                *
 ************************************************************************/

/*
 *    Define error offset for each module : for documentation only
 *    These offsets are defined in the error definition file of the 
 *    specific module 
 *
#define  ccsErrOffset     0
#define  msgErrOffset   100
#define  dbErrOffset    200
#define  timsErrOffset  300
#define  scanErrOffset  400
#define  alrmErrOffset  500
#define  evtErrOffset   700
#define  cmdErrOffset   800
#define  logErrOffset   900
*/

#include "ccsErrors.h"       /* General CCS Errors    */


/************************************************************************
 *                           CCS  Constants                             *
 ************************************************************************/


/************************************************************************
 *                          CCS   Data  Types                           *
 ************************************************************************/

typedef int TBD;               /* To avoid compiler to complaints !!! */


/************************************************************************
 *                     CCS Functions prototypes                         *
 ************************************************************************/


vltLOGICAL ccsIsInitialized(void);


ccsCOMPL_STAT ccsGetProcName (   /* Returns the of a process for a given env.*/
    const ccsENVNAME   envName,  /* Env. Name where the process runs */
          ccsPROCNUM   procNum,  /* Process Identifier Number */ 
          ccsPROCNAME  procName, /* Returned Process Name     */
          ccsERROR    *error     /* Returned error structure  */    
    );


ccsCOMPL_STAT ccsGetProcNum (    /* Returns the number of a process for a given env.*/
    const ccsENVNAME   envName,  /* Env. Name where the process runs */
    const ccsPROCNAME  procName, /* Process Name */
          ccsPROCNUM  *procNum,  /* Returned Process Number   */
          ccsERROR    *error     /* Returned error structure  */
    );

ccsCOMPL_STAT ccsGetMyEnvName (   /* Returns the name of the local environment */
          ccsENVNAME   envName,  /* Env. Name where the process runs */
          ccsERROR    *error     /* Returned error structure  */    
    );


ccsCOMPL_STAT ccsStoreMyName (   /* Store in a static variable the process */
    char  *procName              /* name : typically argv[0]               */
    );


ccsCOMPL_STAT ccsGetMyName (    /* Returns argv[0] as the name of the */
    char *procName              /* current process                    */
    );

ccsCOMPL_STAT ccsFindFile (     /* Find a file in the VLT hierarchy */
   const   char *fileName,      /* Name of the file to be found     */
   char         *filePath,      /* Returned full path of the file   */
   vltINT8      *mode,          /*   */
   vltUINT32    *size,          /*   */
   vltINT8      *dirFlag        /*   */
   );


ccsCOMPL_STAT ccsOpenFile ( /* Opens a file in the INTROOT or VLTROOT area  */
    const char *fileName,   /* filename, including sub path if necessary    */
    const char *mode,       /* Open mode : "r", "w", "a", "r+", "w+", "a+"  */
    FILE       **fileHandle,  /* Returned file handle     */
    ccsERROR   *error         /* Returned error structure */
    );


void ccsDebugPrintf(           /* Execute a printf when the program is   */
                               /* is compiled with the '-D DEBUG' option */
                               
    char *format,              /* printf format parameter        */
    ...                        /* list of runt time parameters   */
    );


ccsCOMPL_STAT ccsGetShmPtr (/* Returns the pointer to the CCS Shared Memory */
    char      **shmPtr      /* Returned Pointer to Shared Memory */
    );


/*  The following routines are used by other CCS modules  */

int ccsGetTimeOfDay(ccsTIMEVAL *tp);


ccsCOMPL_STAT ccsGetInstallPath (   /* Get the installation paths */
    char **vltPath,                 /* VLTROOT path  */
    char **intPath,                 /* INTROOT path  */
    char **vltData                  /* VLTDATA path  */    
    );


ccsCOMPL_STAT ccsGetShmId (   /* Returns identifier to the CCS Shared Memory */
    int       *shmId          /* Returned Shared Memory Identifier */
    );
    
ccsCOMPL_STAT ccsFileExist ( /* returns SUCCESS if the file exists */
   char *file                /* Full file path name */
   );
