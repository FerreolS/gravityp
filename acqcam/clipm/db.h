/**************************************************************************
* E.S.O. - VLT project
#
# "@(#) $Id: db.h 179781 2009-02-03 14:37:49Z mcomin $" 
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
* B.GUSTAFSSON 29/11/94 Including specific CCS/LCC include files
* B.GUSTAFSSON 01/12/94 Added C++ support
* B.GUSTAFSSON 02/12/94 Made different definition of dbATTRIBUTE for LCC and 
*                       CCS
* gfilippi  05/03/95  do not include dbCCS.h if RTAP is not available.
* T.EBERT   16/04/96   merged db.h/dbCCS.h/dbLCC.h (multi-read/write)
* T.EBERT   22/04/96   added dbMAX_ELEMENT_SIZE
* T.EBERT   10/07/96   completed list of db-type-strings
* B.GILLI   19/02/03   increased dbMAX_LIST to 20(SPR 20030054)
* B.GILLI   20/07/07   increased dbMAX_LIST to 100(SPR 20070128) allways
*
************************************************************************
*/
#ifndef DB_H
#define DB_H

#ifdef __cplusplus
extern "C" {
#endif

/*
 *  Header Files
 */

#include "ccs.h"
#include "dbDefine.h"

/*
 * Database  view Classes 
 */
 
typedef enum {
    VIEW_UNKNOWN = 1,      /* Database view not recognized     */
    ALIAS,                 /* Point address is an alias        */
    RELATIVE,              /* Point address is relative to ..  */
    ABSOLUTE               /* Point address is absolute        */
} dbVIEW;                  /* List of RTAP defined views       */

/*
 * Residence flag for a Point's values
 */

typedef vltUINT8 dbRESIDENCE;

/*
 * types of attributes
 */

typedef    vltUINT8    dbATTRTYPE;

/*
 * element types.
 */

typedef vltUINT8    dbTYPE;


/* String containing the hierarchical symbolic address,  including record   */
/* and field information according to the  naming rules.                    */
/* Addressing with alias is included in the symbolic addressing by prefixing*/
/* the symbolic address with the view specifier <alias>                     */

typedef  char   dbSYMADDRESS[dbITEM_LEN+1];

/*
 * Different length of definition of dbATTRIBUTE for LCC and CCS
 */

#ifdef MAKE_VXWORKS
typedef  char   dbATTRIBUTE[dbATTRIBUTE_LEN+1]; /* Name of the attribute */
                                                /* including range       */
#else
typedef  char   dbATTRIBUTE[dbITEM_LEN+1];      /* Name of the attribute  */
                                                /* including range        */
#endif

typedef     char    *dbITEM;	    /* Used to reference an attribute   */
                                    /* including a range specification  */
                                   
typedef     char    dbFIELDNAME[dbFIELDNAME_LEN+1];
				    /* Name of a table field 		*/

typedef ccsPOINT_ID	 dbPOINT;   /* Direct point specification       */

typedef struct {                    /* Direct address      */
     vltUINT16         addrLevel;   /* Addressing Level    */
     vltUINT16         reserved;
     dbPOINT           pointId;     /* Point Specification */
}    dbDIRADDRESS;


/*
 * Message formats
 */
typedef struct {                    /* Parameters for DBREAD command	*/
    dbSYMADDRESS  name;             /* Name of entry in database        */
    vltUINT16     maxLen;           /* Maximum number of bytes to return*/
} dbDBREAD_PARAM;

typedef struct {                    /* Reply for RDBREAD command	*/
    dbATTRTYPE    attrType;         /* Type of attribute		*/
    vltUINT16     bufLen;           /* Number of bytes returned		*/
} dbDBREAD_REPLY;

typedef struct {                    /* Parameters of DBWRIT command	*/
    dbSYMADDRESS  name;             /* Name of entry in database        */
    dbATTRTYPE    attrType;         /* Type of attribute		*/
    vltUINT16     recCnt;           /* Number of records written	*/
    vltUINT16     bufLen;           /* Number of bytes in buffer	*/
} dbDBWRIT_PARAM;


typedef union                       /* Value structure for database   */
    {                               /* element                        */
    vltLOGICAL     tLogical;
    vltINT8        tInt8;
    vltUINT8       tUint8;
    vltINT16       tInt16;
    vltUINT16      tUint16;
    vltINT32       tInt32;
    vltUINT32      tUint32;
    vltFLOAT       tFloat;
    vltDOUBLE      tDouble;
    ccsTIMEVAL     tTime;
    unsigned char  *strValue;   /* For all string types */
    vltPOLAR       *polarVal;   /* vltPOLAR: dynamically allocated */
    vltRECTANGULAR *rectVal;    /* vltRECTANGULAR : dynamically allocated */
    vltDBXREF      *xrefVal;    /* vltDBXREF : dynamically allocated */
    } dbVALUE;

/*
 ***********************************************************************
 * MULTIPLE READ/WRITE STRUCTURES                                      *
 ***********************************************************************
 */
#define  dbREAD    1
#define  dbWRITE   2

#define  dbBLOCK   1
#define  dbNORMAL  2

#define  dbMAX_LIST 100

/*
 * Database handle data type
 */
typedef vltINT32 dbHANDLE;

typedef struct dbLIST {          /* Direct point specification            */
     vltUINT8       index;       /* List Identifier                       */
     ccsENVNAME     envName;     /* Name of the RTAP environement         */
     ccsENV_TYPE    envType;     /* Type  "  "   "       "                */
     vltUINT8       writeMode;   /* Read or write mode                    */
     vltUINT32      blockCnt;    /* Numb. of attr. with SUCCESS operation */
     dbHANDLE       dbHandle;    /* Database handle                       */
     vltUINT32      count;       /* number of elements in the array       */   
     vltUINT32      allocCount;  /* number of allocated elements          */   
     void          *list;        /* List of attributes to be read/written */
                                 /* it is an opaque structure, known only */
                                 /* inside the functions implementing     */
                                 /* multiple read and write               */
     dbSYMADDRESS  *names;       /* List of attributes names. It is an    */
                                 /* array with the same number of elements*/
                                 /* The two arrays are kept parallel and  */
                                 /* are not put in a unique atructure to  */
                                 /* Better handle the list structure      */
    } *dbLISTID;

typedef struct {               /* Attribute description                */
     vltUINT32     quality;
     vltINT32      bufSize;    /* Buffer size                          */
     void         *buffer;     /* Pointer to the data buffer   	       */
                               /* with the current values              */
     vltUINT16     recordCnt;  /* Number of record                     */
     dbTYPE       *dataType;   /* Pointer to the array with data types */
     dbATTRTYPE    attrType;   /* Type of attribute                    */
     vltINT32      actual;     /* Actual number of read/written bytes  */
     ccsCOMPL_STAT status;     /* Completion status of the operation   */
} dbATTRINFO;

/*
 ***********************************************************************
 *  ONLINE DATABASE   : Application  Programmer  Routines              *
 ***********************************************************************
 */

/*
 * The following routines allow the user to read/write RTAP database
 * elements by means of the 'Symbolic' names.  
 *                                                                      
 * These routines support :
 *
 * -  Local and remote access.
 * -  Symbolic addressing using the full path name or a one relative
 *    to a current working point.
 * -  Addressing with Aliases.
 */

ccsCOMPL_STAT dbReadSymbolic(     /* Read a single attribute by name        */
  const dbSYMADDRESS pointName,   /* Name of the database point             */
  const char         *attrName,   /* Attribute Name                         */
        dbTYPE       *dataType,   /* Return the data types of attribute     */
        char         *buffer,     /* Where to put data                      */
        vltINT32     size,        /* How big in bytes is buffer             */
        vltINT32     *actual,     /* How many bytes have been put in buffer */
        vltUINT16    *recordCnt,  /* How many records have been read        */
        dbATTRTYPE   *attrType,   /* Scalar, vector or array (returned)     */ 
        ccsERROR     *error       /* returned error structure.              */
);                                /* If the address refers to a remote node */
                                  /* on the LCU, the request is translated  */
                                  /* into the appropriate command .         */

ccsCOMPL_STAT dbWriteSymbolic(    /* Write a single attribute by name       */
  const dbSYMADDRESS  pointName,  /* Name of the database point             */
  const char          *attrName,  /* Attribute Name                         */
  const dbTYPE        *dataType,  /* Data types of attribute  (input)       */
  const char          *buffer,    /* data to be written                     */
        vltINT32      size,       /* Buffer lenght in bytes                 */
        vltINT32      *actual,    /* How many bytes have been written       */
        vltUINT16     recordCnt,  /* How many records to write (input)      */
        dbATTRTYPE    attrType,   /* Scalar, vector or array (input)        */ 
        ccsERROR      *error      /* returned error structure.              */
);                                /* If the address refers to a remote node */
                                  /* on the LCU, the request is translated  */
                                  /* into the appropriate command .         */

/*
 * The following routines allow the user to read/write RTAP database
 * elements by means of internal indexes.
 *
 * These routines support :
 *
 * -  Local access.
 * -  Addressing with internal indexes.
 */                                                                       

ccsCOMPL_STAT dbReadDirect(      /* read an attribute by internal address  */
  const dbDIRADDRESS *dirAddr,   /* Internal address of the attribute      */
        dbTYPE       *dataType,  /* Return the data types of attribute     */
        char         *buffer,    /* Where to put data                      */
        vltINT32     size,       /* How big in bytes is buffer             */
        vltINT32     *actual,    /* How many bytes have been read          */
        vltUINT16    *recordCnt, /* How many records have been read        */
        dbATTRTYPE   *attrType,  /* Scalar, vector or array (returned)     */
        ccsERROR     *error      /* returned error structure.              */
);

ccsCOMPL_STAT dbWriteDirect(    /* write a single attribute by internal add.*/
  const dbDIRADDRESS *dirAddr,  /* Full address of the entry to be accessed */
  const dbTYPE       *dataType, /* Data types of attributes (input)         */
  const char         *buffer,   /* data to be written                       */
        vltINT32     size,      /* How big in bytes is buffer               */
        vltINT32     *actual,   /* How many bytes have been written         */
        vltUINT16    recordCnt, /* How many records to write (input)        */
        dbATTRTYPE   attrType,  /* Scalar, vector or array (input)          */
        ccsERROR     *error     /* returned error structure.                */
    );
 
/*
 ***********************************************************************
 * MISCELLANEOUS ROUTINES : these routines provide the application     *
 * programmer general services.                                        *
 ***********************************************************************
 */
dbTYPE dbNameToType (      /* Returns the CCS data type from the data  */
   char  *name             /* type name such as : dbvltINT8,dbFLOAT .. */
   );

void dbMemCpy (             /* Copy single value from buf2 to buf1  */
                            /* according to data type               */    
   char    *buf1,           /* Pointer to the destination parameter */
   char    *buf2,           /* Pointer to the value to be copied    */
   dbTYPE  dataType         /* CCS data type                        */
   );

ccsCOMPL_STAT dbFillBuf (   /* Fill Up the user buffer according  */
                            /* to data type                       */
   char       **userBuf,    /* User buffer                        */
   char       *value,       /* Pointer to value to be copied      */
   vltINT32   *bufSize,     /* Returned current buffer size       */      
   dbTYPE     dataType      /* CCS data type                      */
   );

ccsCOMPL_STAT dbStrToDe(    /* Convert a string into the equivalent */
                            /* data element value.                  */
   dbTYPE   dataType,       /* Data Type the string is to be converted into */                            
   char    *string,         /* String to be converted               */
   void    *buffer          /* Pointer to the data element value    */
   );

ccsCOMPL_STAT dbDeToStr(    /* Formats a data element value into a string */
   dbTYPE     dataType,     /* Type of data to be formatted */
   void      *buffer,       /* Pointer tothe data element   */
   char      *string,       /* ASCII formatted value        */
   vltINT8    fieldWidth,   /* Specifies the formatting     */
   vltUINT16  precision     /* Maximum number of significant digits. */
   );

/*
 ***********************************************************************
 * QUERY OPERATIONS to get information about the DB structure ,like :  *
 * size, type, attributes...                                           *
 ***********************************************************************
 */
ccsCOMPL_STAT dbGetAlias(           /* Get the point's alias name   */
  const dbSYMADDRESS  pointName,    /* Point name                */
        dbSYMADDRESS  alias,        /* Returned alias            */
        ccsERROR      *error        /* Returned error structure. */
  );

ccsCOMPL_STAT dbGetAttrNames(       /* Get the list of attributes of a point */
  const dbSYMADDRESS pointName,     /* Point name                            */
        vltUINT8     *attrCount,    /* number of attributes                  */
        vltUINT8     **ains,        /* Attributes Identification Numbers     */
        dbATTRIBUTE  **names,       /* Attributes  names                     */
        ccsERROR     *error         /* returned error structure.             */
  );

ccsCOMPL_STAT dbGetAttrInfo(        /* Get info about one attribute         */
  const dbSYMADDRESS pointName,     /* Point name                           */
  const dbATTRIBUTE  attrName,      /* Attribute name   (optional)          */
        dbATTRTYPE   *type,         /* attribute type :scalar,vector,table  */
        vltUINT8     *fieldCnt,     /* number of fields                     */
        vltUINT16    *recCnt,       /* number of records                    */
        vltUINT32    *recSize,      /* size in bytes of a record or element */
        vltUINT16    *recsUsed,     /* number of records used : specifies   */
                                    /* the index of the last used record.   */
        dbTYPE       *dataType,     /* returned data type.                  */
        ccsERROR     *error         /* returned error structure.            */
  );

ccsCOMPL_STAT dbGetCwp(             /* Get the current working point    */
  const ccsENVNAME    envName,      /* Which database ??                */
        dbSYMADDRESS  pointName,    /* Hierarchical name returned       */
        ccsERROR      *error        /* returned error structure.        */
  );

ccsCOMPL_STAT dbGetDirAddr(         /* Get the direct address for a point or*/
                                    /* or an attribute                      */
  const dbSYMADDRESS  pointName,    /* Point name                           */
  const char          *attrName,    /* Attribute name                       */
        dbDIRADDRESS  *dirAddr,     /* Full direct address of the point     */
        ccsERROR      *error        /* returned error structure.            */
  );

ccsCOMPL_STAT dbGetFieldNames(      /* Get field name of a table attr.      */
  const dbSYMADDRESS  pointName,    /* Point name                           */
  const dbATTRIBUTE   attrName,     /* attribute name                       */
        vltUINT8      *fieldCnt,    /* number of fields                     */
        dbFIELDNAME   **fieldName,  /* array of names                       */
        dbTYPE        **dataType,   /* describes the record format          */
        ccsERROR      *error        /* returned error structure.            */
  );

ccsCOMPL_STAT dbGetFamily(       /* Get direct address of myself, my parent */
                                 /* and the children points                 */
  const dbSYMADDRESS  pointName,    /* Point name                           */
        vltUINT16     *myPlin,      /* My direct address                    */
        vltUINT16     *parentPlin,  /* Parent direct address                */
        vltUINT16     *childCnt,    /* number of child  points              */
        vltUINT16     **childPlin,  /* Array of child direct addresses      */
        ccsERROR      *error        /* returned error structure.            */
  );

ccsCOMPL_STAT dbGetFamilyNames(      /* Get symbolic address of the parent */
                                     /* and the children points            */
   const dbSYMADDRESS  pointName,    /* Point name                         */
   dbVIEW              addrType,     /* Address Type : ALIAS, ABSOLUTE,    */
                                     /*                RELATIVE            */
   dbSYMADDRESS        parentName,   /* Parent Name                        */
   vltUINT16          *childCnt,     /* number of child  points            */
   dbSYMADDRESS       **childName,   /* Array with of children point names */
   ccsERROR           *error         /* returned error structure.          */
   );  

ccsCOMPL_STAT dbDirAddrToName(       /* Get symbolic address from PLIN     */
   dbDIRADDRESS      *dirAddr,       /* Pointer to direct address          */
   dbVIEW            addrType,       /* Address Type : ALIAS, ABSOLUTE,    */
                                     /*                RELATIVE            */
   dbSYMADDRESS      pointName,      /* Symbolic Point Name                */
   ccsERROR          *error          /* returned error structure.          */
   );  

ccsCOMPL_STAT dbGetParent(          /* Get my PLIN and my parent's  PLIN   */
  const dbSYMADDRESS  pointName,    /* Point name                          */
        vltUINT16     *parentPlin,  /* Parent direct address               */
        ccsERROR      *error        /* returned error structure.           */
  );

ccsCOMPL_STAT dbAliasToName(        /* get the symbolic address from an alias */
  const dbSYMADDRESS   alias,       /* Point's alias name        */
  dbSYMADDRESS         pointName,   /* Point's symbolic address  */
  ccsERROR             *error       /* returned error structure. */
  );

ccsCOMPL_STAT dbGetClassID(
   const dbSYMADDRESS  pointName,
   dbDIRADDRESS       *dirAddr,
   vltUINT16          *classId, 
   ccsERROR           *error
   );

/*
 ***********************************************************************
 * Control operations to lock/unlock points or branches ...            *
 ***********************************************************************
 */

ccsCOMPL_STAT dbLockPoint(           /* Lock a point in the database    */
  const dbSYMADDRESS  pointName,     /* Reference to the database point */
  const vltUINT32     seconds,       /* Timeout for retry               */
        ccsERROR      *error         /* Returned error structure.       */
    );      

ccsCOMPL_STAT dbUnlockPoint(         /* UnLock a point in the database  */
  const dbSYMADDRESS  pointName,     /* Reference to the database point */
        ccsERROR      *error         /* Returned error structure.       */
    ); 

ccsCOMPL_STAT dbSetCwp(              /* Set the Current Working Point   */
  const dbSYMADDRESS pointName,      /* Reference to the database point */
        ccsERROR     *error          /* Returned error structure.       */
    );

/*
 ***********************************************************************
 * MULTIPLE READ/WRITE FUNCTIONS                                       *
 ***********************************************************************
 */
ccsCOMPL_STAT dbListCreate(    /* Create an empty a list of attributes */
   const char       *listName,      /* Name of the list              */
   const ccsENVNAME  envName,       /* Name of the RTAP environement */
   vltINT8      useFlag,          
   dbLISTID    *listId,        /* Returned list identifier      */
   ccsERROR    *error          /* Returned error structure      */    
   );

ccsCOMPL_STAT dbListGetId(     /* Get the list identifier from the name */
   const char  *listName,      /* Name of the list              */
   dbLISTID    *listId,        /* Returned list identifier      */
   ccsERROR    *error          /* Returned error structure      */    
   );

ccsCOMPL_STAT dbListDestroy(   /* Destroy a list of attributes */
   dbLISTID   listId,          /* List Identifier          */
   ccsERROR  *error            /* Returned error structure */    
   );

ccsCOMPL_STAT dbListChangeEnv( /* Replace the environment name with another */
   dbLISTID     listId,        /* List Identifier          */
   ccsENVNAME   newEnvName,    /* New environment name     */
   ccsERROR    *error          /* Returned error structure */ 
   );

ccsCOMPL_STAT dbListAdd(       /* Insert a new element in the list */
   dbLISTID        listId,     /* List Identifier                  */
   dbSYMADDRESS    attrName,   /* Symbolic address of the database attribute */  
   vltINT32        bufSize,    /* Buffer size                      */ 
   void           *buffer,     /* Pointer to the buffer with current values  */
   ccsERROR       *error       /* Returned error structure */
   );
   
ccsCOMPL_STAT dbListModify(    /* Modify the configuration parameters of */
                               /* an element of the list                 */
   dbLISTID        listId,     /* List Identifier          */
   dbSYMADDRESS    attrName,   /* Symbolic address of the database attribute */
   vltINT32        bufSize,    /* Buffer size              */ 
   void           *buffer,     /* Pointer to the buffer with current values  */
   ccsERROR       *error       /* Returned error structure */
   );
   
ccsCOMPL_STAT dbListRemove(    /* Remove an element from the list */
   dbLISTID      listId,       /* List Identifier          */
   dbSYMADDRESS  attrName,     /* Symbolic address of the database attribute */
   ccsERROR     *error         /* Returned error structure */ 
   );

ccsCOMPL_STAT dbListExtract(   /* Get the information returned after */
                               /* a read/write operation             */
   dbLISTID        listId,     /* List Identifier                    */      
   dbSYMADDRESS    attrName,   /* Symbolic address of the database attribute */
   dbATTRINFO     *attrInfo,
   ccsERROR       *error       /* Returned error structure */  
   );

ccsCOMPL_STAT dbMultiWrite(    /* Write list of attributes */
    dbLISTID   listId,         /* List Identifier          */
    vltUINT8   writeMode,      /* Write mode : dbBLOCK, dbNORMAL */
    ccsERROR   *error          /* Returned error structure */ 
    );

ccsCOMPL_STAT dbMultiRead(     /* Read list of attributes */
    dbLISTID    listId,        /* List Identifier          */
    ccsERROR   *error          /* Returned error structure */ 
    );

ccsCOMPL_STAT dbGetCe(const dbSYMADDRESS  pointName,
                      const dbATTRIBUTE   attribute,
                            char         *definition,
		            vltINT32      size,
                            ccsERROR     *error);

ccsCOMPL_STAT dbSetCe(const dbSYMADDRESS  pointName,
                      const dbATTRIBUTE   attribute,
                            char         *definition,
                            ccsERROR     *error);

#ifdef MAKE_VXWORKS
#include "dbLCC.h"
#else

#ifdef CCS_LIGHT
#include "dbLite.h"
#endif

/*
 * normally the following is included, but in a non-RTAP installation (LCU only),
 * the CCS part would produce errors, because it uses RTAP include files
 * Executable that have to run on both cases, like lcc<xx> shall 
 * define the flag RTAP_EXCLUDE
 */
#ifndef RTAP_EXCLUDE
#include "dbCCS.h"
#endif

#endif /* MAKE_VXWORKS */

#ifdef __cplusplus
}
#endif

#endif /*!DB_H*/
