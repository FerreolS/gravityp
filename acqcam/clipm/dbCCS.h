/**************************************************************************
* E.S.O. - VLT project
#
# "@(#) $Id: dbCCS.h 77119 2003-01-27 22:06:08Z tebert $" 
*
* <db.h>  -  CCS/ON-LINE DATABASE Interface File
*
* who        when      what
* ---------  --------  ----------------------------------------------
* B.GILLI    29/04/93  Preliminary version
* M.COMIN    27/05/93  Second Preliminary version
* B.GILLI    27/05/93  for PDR
* M.COMIN    28/05/93  Removed 'Global Data'.
* M.COMIN    07/06/93  Introduced typedef    dbHANDLE
*                      Introduced define     LOCAL_ENV
*                      Introduced typedef    dbVIEW
*                      Introduced typedef    dbFIELD
*                      Typedef dbITEMNAME shortened into dbITEM
*                      New routine dbGetHandle().
*                      New routine dbExit().
*                      Final version for routines dbAddVector(), dbAddTable(),
*                      dbAddScalar() and dbAddPoint().
*                      Parameter is added to routine dbCopyAttr() 
*
* M.COMIN    08/06/93  New routine dbNameToType() 
*
* M.COMIN    30/06/93  Introduced typedef    dbPOINT
*                      Changed definition of typedef dbDIRADDRESS
*                      Redefined param buffer in dbWriteSymbolic(), 
*                      dbReadSymbolic(), dbWriteDirect(), dbReadDirect().
*                      New routine dbMemCpy().
*
* M.COMIN    01/07/93  All parameter definition like  *names[dbMAX_ATTR_CNT] 
*                      are changed to  *names.
*                      Changed all data types according to ccs.h
*                      Parameter 'attrName' added to dbGetDirAddr()
*
* M.COMIN    05/07/93  Add define for direct addressing levels.
*                      Add define EMPTY_STRING.
*
* M.COMIN    07/07/93  New parameter added to routine dbGetFieldNames() 
* 
* M.COMIN    14/07/93  Typedef of DB_XREF moved here from  ccs.h.   
*                      Functions related to the scan system removed ... they
*                      will be included in a separate file. 
*
* M.COMIN    19/07/93  New parameter added to routine dbCopyAttr()      
*
* M.COMIN    20/07/93  Add routine dbFillBuf().
* M.COMIN    26/07/93  Add parameter error to dbExit()
* M.COMIN    12/08/93  Changed name of 1-st parameter in dbWriteDirect()
*                      and dbReadDirect()
*
* M.COMIN    14/08/93  All string parameters of type dbSYMADDRESS,dbITEM,
*                      ccsENVNAME, etc .. are passed without '*'.
*
*                      All string input data structures have the qualifier 
*                      'const' as prefix.
* 
*                      Changed tupe of param lockTime in dbXferLock() 
*                      to unsigned long.
* 
* M.COMIN   19/08/93   Add new typedef dbITEM. 
*                      Add  rtap/database.h as include file. 
*  
*                      Changed type of attrName to dbITEM in routines 
*                      dbReadSymbolic(),dbWriteSymbolic() and dbGetDirAddr().
*
* M.COMIN   15/09/93   Changed type of member of dbDIRADDRESS according to
*                      VLT data types (rtDbDirectType replaced with vltUINT16)
*
* --------------------
*
* M.COMIN   20/01/94   Changed according to the CCS design
* M.COMIN   28/03/94   Add define dbALL_ENV. 
* M.COMIN   19/04/94   Add define dbEMPTY
* M.COMIN   15/06/94   Add error to dbExit().
*                      dbITEM changed to (char *)
*
* M.COMIN   12/07/94   Add function prototypes for event handling
*                      dbParseEvent() dbAttachEvent() dbEnableEvent(0
*                      dbDisableEvent()
*                      Add data structures for event handling.
*
* M.COMIN   13/07/94   Complete the list of query routines.
*                      Add dbStrToDe() and dbDeToStr()
*
* M.COMIN   15/07/94   Add routines to handle lists of attributes
*
* M.COMIN   18/07/94   Add routines dbMultiRead(0 and dbMultiWrite(3)
*                      Add routines dbGetFamilyNames() and dbPlinToName(3)
*
* G.CHIOZZI 06/09/95   Modified dbLISTID structure to fulfill implementation
*                      requirements 
*
* T.HERLIN  20/09/94   Changed dbPlinToName() to dbDirAddrToName().
*                      Renamed dbNameToAlias() to dbAliasToName()
*                      Changed dbGetFamilyNames() 'childCnt' value to ptr.
*                      Removed database event related information. This is
*                      information is now residing in the evt.h of the
*                      event module.
* T.HERLIN  23/09/94   Changed interfaces of:
*                      dbGetFamilyNames(),dbGetFamily(),dbGetAttrNames()
*                      to return pointer to array allocated within the
*                      library.
* T.HERLIN  29/09/94   Removed myplin from dbGetParent() due to LCC compliance
* G.CHIOZZI 19/10/94   Added defines for TRUE and FALSE
* M.COMIN   21/11/94   Add definitions for error menmonics
* B.GUSTAFSSON 29/11/94 Created CCS specific include file
* J.I.LEON  31/05/95   Add conditional compilation for CCS LITE
* T.EBERT   30/06/95   include Rtap header files by ccsRtap.h
* T.EBERT   09/02/96   Added prototype of dbReopen()
* T.EBERT   16/04/96   merged db.h/dbCCS.h/dbLCC.h (multi-read/write)
* T.EBERT   25/10/96   Added defines of dbATTR_XXXX (CCS-LITE)
*
* gfilippi  20/04/97   FUNCTIONS USED by INS not defined for CCSlite
*
* bgilli    05/05/98   Mauro seems to have used this file to move dbPrivate
*                      declarations half public!!!
*                      Added dbParseClassAttr from dbPrivate.h
* T.EBERT   31/08/98   Added dbHandlToEnvname()
* B.Gilli   16/05/00   Added define for branch separator.

************************************************************************
*/

/*
 *  Header Files
 */
#include "ccsRtap.h"
#include "dbErrors.h"

/* Define specific delimiters used to specify an attribute */

#define ATTR_CHAR      '.'
#define ATTR_STR       "."    
#define ENV_DELIMITER  '@'
#define BRANCH_CHAR      ':'
#define BRANCH_STR       ":"
   

/************************************
 *  Define mnemonics for database   *
 *  access commands                 *
 ************************************/
#define  LCU_READ          "DBREAD"
#define  LCU_WRITE         "DBWRIT"
#define  LCU_MREAD         "DBMREAD"
#define  LCU_MWRITE        "DBMWRIT"
#define  LCU_ATTR_INFO     "DBGAINF"
#define  LCU_DIR_ADDR      "DBGDIRA"
#define  LCU_ALIAS         "DBGALS"
#define  LCU_ALTONAM       "DBATON"
#define  LCU_GATTRNAM      "DBGANAM"
#define  LCU_GCWP          "DBGCWP"
#define  LCU_FAMILY        "DBGFAM"
#define  LCU_FAM_NAMES     "DBGFAMN"
#define  LCU_FIELD_NAM     "DBGFNAM"
#define  LCU_GET_PARENT    "DBGPAR"
#define  LCU_LOCK_POINT    "DBLCKP"
#define  LCU_UNLOCK_POINT  "DBULCKP"
#define  LCU_SCWP          "DBSCWP"
#define  LCU_DIRTONAME     "DBDATON"

#define  CCS_GETCE         "DBGECE"
#define  CCS_SETCE         "DBSECE"

/*
 ************************************************ 
 *  Database types and constant definitions     *
 ************************************************
 */
/*
 * Field descriptor data type (for tables).
 */
typedef rtTableField  dbFIELD;

/*
 * element types.
 */
typedef rtDbXref  DB_XREF;

/* 
 * Define quality of an attribute
 */
#define dbATTR_OK        rtATTR_OK        /* no known errors in data */
#define dbATTR_SUSPECT   rtATTR_SUSPECT   /* depends on a suspect,   */
	                                  /* error or disabled value */
#define dbATTR_ERROR     rtATTR_ERROR     /* calc. engine got math error */
#define dbATTR_DISABLED  rtATTR_DISABLED  /* scan point disabled or ce   */
				          /* operation indicator disabled */
/*
 * this variable is a pointer to an array containing the sizes of 
 * all database types:
 * vltLOGICAL: size = dbDataElemSize[dbLOGICAL]
 * This way has to be used, because the array rtDateElemSize[] provided
 * by RTAP is declared as extern and the array dbDataElemSize[] has to
 * be initialized by the compiler!
 */
extern vltUINT16 *dbDataElemSize;

/*
 ***********************************************************************
 * MISCELLANEOUS ROUTINES : these routines provide the application     *
 * programmer general services.                                        *
 ***********************************************************************
 */
ccsCOMPL_STAT dbExit (        /* Close a specific database connection or all */
  const ccsENVNAME  envName,  /* Name of the environment */
        ccsERROR   *error     /* Returned error structure. */
 ); 
       
ccsCOMPL_STAT dbOpen(          /* Establishes a connection to a database  */
   const ccsENVNAME  envName,  /* Name of the environment  */
   ccsERROR         *error     /* Returned error structure. */
   );

ccsCOMPL_STAT dbReopen(        /* Re-opens a connection to a database */
   const ccsENVNAME   envName, /* Name of the environment             */
         ccsERROR    *error    /* error data structure                */
   );


/*
 ***********************************************************************
 * INTERNAL ROUTINES : these routines are internally used by other     *
 * CCS modules                                                         *
 ***********************************************************************
 */

ccsCOMPL_STAT dbParseEnv (          /* Extract the environment specified    */
                                    /* in the symbolic address              */
   const dbSYMADDRESS  fullName,    /* Full point Name with env. specifier  */
         dbSYMADDRESS  **pointName, /* Point address without env. specifier */
         ccsENVNAME    envName,     /* Environment  Name                    */
         ccsERROR      *error       /* Returned error structure             */ 
);


ccsCOMPL_STAT dbGetHandle(      /* Get the handle to an open database */
  const ccsENVNAME  envName,    /* Environment name : NULL if local   */ 
        dbHANDLE    *dbHandle,  /* Returned handle to the database    */
        ccsENV_TYPE *envType,   /* Returned environment type          */
        ccsERROR    *error      /* Returned error structure           */
);


ccsCOMPL_STAT dbDirAddrToNameOnEnv(  /* returns sybolic name of direct       */
        dbHANDLE       dbHandle,     /* address in remote database           */
        dbDIRADDRESS   *dirAddr, 
	dbVIEW         addrType,
	dbSYMADDRESS   pointName,       
	ccsERROR       *error
);

/*
 ***********************************************************************
 * FUNCTIONS USED by TCS                                               *
 ***********************************************************************
 */

void dbParseStatus(          /* Enable/disable parsing       */
   vltLOGICAL enable
);

ccsCOMPL_STAT dbXrefToAlias( 
    const ccsENVNAME env,
    const DB_XREF    *addr,
    char       *string,
    ccsERROR   *error
);

ccsCOMPL_STAT dbGetDirAddrOnEnv( 
    const dbSYMADDRESS  symName,
    const dbHANDLE      dbHandle,
    const ccsENV_TYPE   envType,
    ccsENVNAME          envName,
    dbDIRADDRESS       *dirAddr,
    ccsERROR           *error
);


/*
 ***********************************************************************
 *                                            *
 ***********************************************************************
 */

#ifndef CCS_LIGHT

    ccsCOMPL_STAT dbAddTable(            /* Add TABLE attribute to a point   */
      const dbSYMADDRESS pointName,      /* Reference to the database point  */
      const dbATTRIBUTE  attrName,       /* Name of the attribute            */
            vltUINT16    recordCnt,      /* Number of records (table rows)   */
            vltUINT8     fieldCnt,       /* Number of fields  (table columns)*/
            dbFIELD      *fieldType,     /* field descriptions               */
            ccsERROR     *error          /* returned error structure         */
        );

    ccsCOMPL_STAT dbLockConfig(          /* Lock a database                 */
      const ccsENVNAME   envName,        /* Name of the db(=env) to lock    */
            ccsERROR     *error          /* Returned error structure.       */
        );     

    ccsCOMPL_STAT dbUnlockConfig(        /* UnLock a database               */
      const ccsENVNAME   envName,        /* Name of the db(=env) to unlock  */
            ccsERROR     *error          /* Returned error structure.       */
        );

#endif

ccsCOMPL_STAT dbParseClassAttr(      /* Parses and attribute name taking    */
        dbSYMADDRESS  point,         /* into account that it can be a static*/
        ccsERROR      *error         /* attribute. Returns an alias SYMADDR */
    );

ccsCOMPL_STAT dbHandleToEnvname (
    dbHANDLE     dbHandle,
    ccsENVNAME   envName,
    ccsENV_TYPE *envType,
    ccsERROR    *error
);
