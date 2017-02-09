/**************************************************************************
* E.S.O. - VLT project
#
# "@(#) $Id: dbLite.h 120020 2005-07-12 13:46:06Z tebert $" 
*
* <dbLite.h>  -  CCS_Lite/ON-LINE DATABASE Interface File
*
* who        when      what
* ---------  --------  ----------------------------------------------
* B.GILLI    18/07/00  Added standard header with RCS
* B.GILLI    18/09/01  dbSwapBytesInBuffer returns a completion code.
* B.GILLI    04/06/02  suppressed ifdef for class/poinClass, always use poinClass.
************************************************************************
*/
#ifndef CCS_LITE_PUB_H
#define CCS_LITE_PUB_H

#include <string.h>
#include "ccsLite.h"

#define rtMAX_PROC_NAME_LEN     20
/*
 * minimum values of numeric data element types
 */
#define rtMIN_LOGICAL           0  
#define rtMIN_INT8             -128 
#define rtMIN_UINT8             0 
#define rtMIN_INT16            -32768 
#define rtMIN_UINT16            0   
#define rtMIN_INT32             INT_MIN 
#define rtMIN_UINT32            0  
#define rtMIN_FLOAT            -FLT_MAX  
#define rtMIN_DOUBLE           -DBL_MAX 

/*
 * maximum values of numeric data element types
 */
#define rtMAX_LOGICAL           1     
#define rtMAX_INT8              127  
#define rtMAX_UINT8             255    
#define rtMAX_INT16             32767  
#define rtMAX_UINT16            65535    
#define rtMAX_INT32             INT_MAX
/* 
 * Unsigned longs are limited to 2G for the calculation engine.
 * This makes all calculation in the CE able to be done as longs 
 * and allows for more error checking.
 */
#define rtMAX_UINT32            INT_MAX  
#define rtMAX_FLOAT             FLT_MAX 
#define rtMAX_DOUBLE            DBL_MAX 


#define rtMAX_POINT_CNT         65535   
#define rtMAX_ATTR_CNT          255  
#define rtMAX_ELEM_CNT          65535
#define rtMAX_FIELD_CNT         255
#define rtMAX_RECORD_CNT        65535
#define rtMAX_FIELD_CNT         255
#define rtMAX_CHILDREN          512 

#define rtSCALAR  0
#define rtVECTOR  1
#define rtTABLE   2

#define rtATTR_OK        0 
#define rtATTR_SUSPECT   1
#define rtATTR_ERROR     2 
#define rtATTR_DISABLED  3

#define rtDIR_INVALID    0
#define rtDIR_PLIN       1
#define rtDIR_AIN        2
#define rtDIR_REC_EL     3
#define rtDIR_FIELD      4
#define rtUNDEFINED             ((rtDeType) 0)
 
#define rtLOGICAL               ((rtDeType) 1)
#define rtINT8                  ((rtDeType) 2)
#define rtUINT8                 ((rtDeType) 3)
#define rtINT16                 ((rtDeType) 4)
#define rtUINT16                ((rtDeType) 5)
#define rtINT32                 ((rtDeType) 6)
#define rtUINT32                ((rtDeType) 7)
#define rtFLOAT                 ((rtDeType) 8)
#define rtDOUBLE                ((rtDeType) 9)
#define rtPOLAR                 ((rtDeType) 10)
#define rtRECTANGULAR           ((rtDeType) 11)
#define rtBYTES4                ((rtDeType) 16)
#define rtBYTES8                ((rtDeType) 17)
#define rtBYTES12               ((rtDeType) 18)
#define rtBYTES16               ((rtDeType) 19)
#define rtBYTES20               ((rtDeType) 20)
#define rtBYTES32               ((rtDeType) 21)
#define rtBYTES48               ((rtDeType) 22)
#define rtBYTES64               ((rtDeType) 23)
#define rtBYTES80               ((rtDeType) 24)
#define rtBYTES128              ((rtDeType) 25)
#define rtBYTES256              ((rtDeType) 26) 
#define rtDB_XREF               ((rtDeType) 28)
#define rtDATE                  ((rtDeType) 29)
#define rtTIME_OF_DAY           ((rtDeType) 30)
#define rtABS_TIME              ((rtDeType) 31)

#define rtMAX_DE_TYPE           31

#define rtQUERY_CONN_INFO       0       /* get connection information */
#define rtQUERY_CATEG_NAMES     1       /* get category names (all 32) */
#define rtQUERY_GROUP_NAMES     2       /* get group names (all 4) */
#define rtQUERY_SPACE           3       /* get free space info */
#define rtQUERY_CFI             4       /* get Config. interlock data */
 
/* Point level */
#define rtQUERY_ALIAS           10      /* get Point Alias */
#define rtQUERY_ATTR_CNT        11      /* get Attribute count of Point */
#define rtQUERY_ATTR_NAMES      12      /* get Attribute names for Point */
#define rtQUERY_ATTR_ORDER      13      /* get Attribute order list */
#define rtQUERY_CATEGORIES      14      /* get Point categories */
#define rtQUERY_CE_OPER         15      /* get CE operation indicator */
#define rtQUERY_EXPR_ORDER      16      /* get expression order type */
#define rtQUERY_FIRST_CHILD     17      /* get first child PLIN */
#define rtQUERY_LRL             18      /* get logical reference list */
#define rtQUERY_NEXT_SIBLING    19      /* get next sibling PLIN */
#define rtQUERY_PARENT          20      /* get parent PLIN */
#define rtQUERY_PT_NAME         21      /* get rtPointName of Point */
#define rtQUERY_RESIDENCE       22      /* get Point residence */
#define rtQUERY_USAGE           23      /* get Point usage flags */
#define rtQUERY_ATTR_IDS        24      /* get Attribute ID numbers (AINs) */
#define rtQUERY_ALPHA_ATTRS     25      /* get sorted attribute names */
#define rtQUERY_PT_CLASS        26      /* get point class */
#define rtQUERY_PTS_IN_CLASS    27      /* get points matching a class */
#define rtQUERY_PTS_IN_CATEG    28      /* get points matching categories */
#define rtQUERY_CE_DEP_REF      29      /* get pts with CE ref's to a pt */
#define rtQUERY_CE_DEP_UPD      30      /* get pts CE would update for a pt */
 
/* Attribute level */
#define rtQUERY_ATTRIBUTE       40      /* get Attribute specs (type, size) */
#define rtQUERY_ATTR_ACCESS     41      /* get read write access to Attribute */
#define rtQUERY_ATTR_NAME       42      /* get name of Attribute */
#define rtQUERY_DE_TYPE         43      /* get data-element type(s) */
#define rtQUERY_EVENT           44      /* get event triggers */
#define rtQUERY_FIELD_NAMES     45      /* get field names */
#define rtQUERY_GROUPS          46      /* get Attribute access groups */
#define rtQUERY_DEFINITION      47      /* get unparsed definition */
/* any address level */
#define rtQUERY_DIRECT          60      /* get direct equiv of sym addr */
#define rtQUERY_SYM_ABS         61      /* get abs-path equiv of given addr */
#define rtQUERY_SYM_ALIAS       62      /* get alias equiv of given addr */
#define rtQUERY_SYM_REL         63      /* get rel-path equiv of given addr */
#define rtQUERY_DIRECT_ATTR     64      /* get fully-specified direct addr */

#define rtCONTROL_CE_OPER       0       /* set CE operation indicator */
#define rtCONTROL_LOCK_PT       1       /* lock Points */
#define rtCONTROL_SET_CWP       2       /* set current working position */
#define rtCONTROL_SNAPSHOT      3       /* force snapshot */
#define rtCONTROL_UNLOCK_PT     4       /* unlock Points */
#define rtCONTROL_XFER_LOCK     5       /* transfer locks to another process */
#define rtCONTROL_RUN_CE        6       /* run calculation engine on pts */
#define rtCONTROL_ADD_USAGE     7       /* add usage flags */
#define rtCONTROL_DEL_USAGE     8       /* delete usage flags */
#define rtCONTROL_SET_CFI       9       /* set configuration interlock */
#define rtCONTROL_REL_CFI       10      /* release configuration interlock */
#define rtCONTROL_SET_USAGE     11      /* set usage flags array */
#define rtCONTROL_ENABLE_SNAPS  12      /* enable snapshotting */
#define rtCONTROL_DISABLE_SNAPS 13      /* disable snapshotting */
#define rtCONTROL_REUSE_PLINS   14      /* allow reuse of PLINs */
#define rtCONTROL_RDR_CFI       15      /* obtain a multi-reader CFI */
#define rtCONTROL_CLR_CFI       16      /* clear a dead process's zombie CFI */

#define rtCONFIG_CATEGORIES     1       /* set point categories */
#define rtCONFIG_EXPR_ORDER     2       /* change expression order indic */
#define rtCONFIG_DEFINITION     3       /* change expression */
#define rtCONFIG_GROUPS         4       /* set attr access groups */
#define rtCONFIG_PT_NAME        5       /* change the name of a point */
#define rtCONFIG_MOVE_POINT     6       /* move a point (branch) */
#define rtCONFIG_ALIAS          7       /* change the alias of a point */
#define rtCONFIG_ADD_NULL_PT    8       /* add a null point */
#define rtCONFIG_RESIDENCE      9       /* change point residence */
#define rtCONFIG_DEL_BRANCH     10      /* delete branch */
#define rtCONFIG_COPY_BRANCH    11      /* copy branch */
#define rtCONFIG_COPY_POINT     12      /* copy point */
#define rtCONFIG_PT_CLASS       13      /* change point class */
#define rtCONFIG_DEL_ATTR       14      /* delete attribute */
#define rtCONFIG_COPY_ATTR      15      /* copy attribute */
#define rtCONFIG_ATTR_NAME      16      /* change attribute name */
#define rtCONFIG_ADD_SCALAR     20      /* add scalar attribute */
#define rtCONFIG_ADD_VECTOR     21      /* add vector attribute */
#define rtCONFIG_ADD_TABLE      22      /* add table attribute */
#define rtCONFIG_DEL_BR_CHK     23      /* check if branch deletable */

#define rtCE_DISABLED            0
#define rtCE_ENABLED             1
#define rtCE_DISABLED_NORUN      2
#define rtCE_ENABLED_OPT_A       3

#define rtNATURAL                0
#define rtUSER_DEFINED           1

#define rtMAX_DEFINITION_SIZE 2048         
#define rtMAX_MATCH_CNT       1024

typedef struct 
    {
    rtUInt16 defnSize;
    rtChar   defn[rtMAX_DEFINITION_SIZE];
    } rtQueryDefinition;

typedef rtUInt16  rtCeOperation;
typedef rtUInt16  rtExprOrder;
typedef vltINT32  rtDbConnection;
typedef char      rtPointName[20];
typedef char      rtAlias[20];
typedef char      rtAttributeName[20];
typedef char      rtFieldName[20];
typedef char      rtDbViewName[20];
typedef vltUINT8  rtAttributeType;
typedef vltUINT32 rtDbQuality;
typedef vltUINT16 rtPlin;                   
typedef vltUINT8  rtAin;
typedef vltUINT16 rtDbDirectType;
typedef vltUINT32 rtDbCqAction;
typedef vltUINT8  rtDeType;

typedef struct
    {
    rtUInt16        pointClass;
    rtPlin          plinCnt;
    rtPlin          plins[rtMAX_MATCH_CNT];
    } rtQueryPtsInClass;

typedef enum
    {
    rtSYM_HIERARCHICAL,
    rtSYM_ALIAS,
    rtSYM_UNKNOWN
    } rtDbSymbolicType;

typedef enum
    {
    rtDB_DIRECT,
    rtDB_SYMBOLIC
    } rtDbAddressType;

typedef struct
    {
    rtDbSymbolicType   type;
    char              *name;
    } rtDbSymbolic;

typedef struct
    {
    rtPlin          plin;
    rtAin           ain;
    vltUINT8        reserved;
    vltUINT16       startRecEl;
    vltUINT16       endRecEl; 
    vltUINT8        startField; 
    vltUINT8        endField;
    } rtPointId;


typedef struct
    {
    rtDbDirectType  type;
    vltUINT16       reserved;
    rtPointId      *pointId;
    } rtDbDirect;

typedef struct
    {
    rtDbAddressType         type;
    union
        {
        rtDbDirect      direct; 
        rtDbSymbolic    symbolic;
        } ds;
    } rtDbAddress;

typedef struct
    {
    rtDbAddress     addr; 
    rtDbQuality     quality;
    vltUINT16       recordCnt;
    rtAttributeType attrType;
    vltUINT8        reserved;
    rtDeType       *deTypes;
    vltUINT8       *buffer;
    vltINT32        size;
    vltINT32        actual;
    } rtDbReadWrite;

typedef struct
    {
    rtDbDirectType  nameType;
    rtPointId       pointId;
    } rtDbXref;

typedef struct
    {
    rtDbCqAction    action; 
    vltUINT16       addrCnt;
    vltUINT8        reserved[2];
    rtDbAddress    *addr;
    vltUINT8       *buffer; 
    vltUINT32       size;
    vltUINT32       actual;
    } rtDbCq;

typedef struct
    {
    rtFieldName     name;
    rtDeType        deType;
    vltUINT8        reserved[3];
    } rtTableField;

typedef struct
    {
    rtPlin          currPoint;
    vltUINT32       groups:4;       /* rtUInt8 */
    vltUINT32       reserved:4;     /* rtUInt8 */
    vltUINT8        security;
    vltUINT32       ptReadCats;
    vltUINT32       ptWriteCats;
    ccsENVNAME      procEnv;
    } rtQueryConnInfo;

typedef struct
    {
    rtAttributeType type;
    vltUINT8        fieldCnt;
    vltUINT16       recElCnt;
    vltUINT32       recordSize;
    vltUINT16       recsUsed;
    vltUINT8        reserved[2];
    } rtQueryAttribute;

typedef struct
    {
    vltUINT8      deCnt;
    rtDeType      deTypes[rtMAX_FIELD_CNT];
    } rtQueryDeType;

typedef struct
    {
    rtAin             ain;
    vltUINT8          reserved;
    rtAttributeName   name;
    } rtAinAndName;

typedef struct
    {
    rtAin           nameCnt;
    vltUINT8        reserved;
    rtAinAndName    attrs[rtMAX_ATTR_CNT];
    } rtQueryAlphaAttrs;

typedef struct
    {
    rtPlin          thisPoint;  
    rtPlin          parent;
    rtPlin          childCnt;
    rtPlin          children[rtMAX_CHILDREN];
    } rtQueryLrl;

typedef struct
    {
    vltUINT8         nameCnt;
    rtFieldName      names[rtMAX_FIELD_CNT];
    } rtQueryFieldNames;

typedef rtDbXref        rtQueryDbDirect;

typedef struct {
    enum {
        READER,
        CONFIG 
    }                   ilkType;  
    rtUInt32            count;
    rtUInt32            pid;
    rtChar              procName[rtMAX_PROC_NAME_LEN];
    rtProcessNum        procNum;
    rtLogical           terminated; 
    rtUInt8             reserved[2];
    rtTime              time;
} rtQueryCfi;

typedef rtQueryCfi      rtControlClrCfi;

typedef union
    {
    rtLogical       deLogical;
    rtInt8          deInt8;
    rtUInt8         deUInt8;
    rtInt16         deInt16;
    rtUInt16        deUInt16;
    rtInt32         deInt32;
    rtUInt32        deUInt32;
    rtFloat         deFloat;
    rtDouble        deDouble;
    rtRectangular   deRectangular;
    rtPolar         dePolar;
    rtBytes4        deBytes4;
    rtBytes8        deBytes8;
    rtBytes12       deBytes12;
    rtBytes16       deBytes16;
    rtBytes20       deBytes20;
    rtBytes32       deBytes32;
    rtBytes48       deBytes48;
    rtBytes64       deBytes64;
    rtBytes80       deBytes80;
    rtBytes128      deBytes128;
    rtBytes256      deBytes256;
    rtDbXref        deDbXref;
    rtTime          deDate;
    rtTime          deTimeOfDay;
    rtTime          deAbsTime;
    } rtAllDeTypes;

#include "evtLite.h"

typedef struct 
    {
    rtProcessId  trigger;
    rtInt        triggerUId;
    vltUINT16    pIndex;
    vltUINT16    iIndex;
    vltLOGICAL   used;
    union
       {
        rtEventScalarInfo       scalar;
        rtEventVectorInfo       vector;
        rtEventTableInfo        table;
        } svt;
    } dbEVT_QUEUE_ITEM;

typedef struct
    {
    rtPlin     plinT;
    rtAin      ainT;
    rtPlin     plin;
    rtAin      ain;
    rtUInt32   sequence;
    rtTime     chgTime;
    rtUInt16   oldIndex;
    rtUInt16   index;
    vltLOGICAL used;
    } dbALRM_QUEUE_ITEM;

#ifdef CCS_PTHREADS
#define  dbMAX_CONNECTIONS  50
#define  dbLINK_DOWN      0
#define  dbLINK_UP        1
#define  dbLINK_NEW       2
#define  dbLINK_CLOSED    3

typedef struct {
    ccsENVNAME  envName;          /* Environment name                */
    vltUINT8    status;           /* Status of the link : UP/DOWN    */
    ccsENV_TYPE type;             /* Type of environment : WS or LCU */
    dbHANDLE    linkId;           /* Handle to access the database   */
    } dbLinkStatus;
#endif

rtDbConnection rtOpenDatabase    (char *env);
vltINT32       rtCloseDatabase   (rtDbConnection dbConnection);
vltINT32       rtReadDatabase    (rtDbConnection dbConnection, rtDbReadWrite *rw);
vltINT32       rtDeToStr         (rtDeType type, void *value, char *string, 
                                  vltINT32 fieldWidth, vltUINT16 precision);
vltINT32       rtQueryDatabase   (rtDbConnection dbConnection, rtDbCq *dbCq);
vltINT32       rtControlDatabase (rtDbConnection dbConnection, rtDbCq *dbCq);
vltINT32       rtConfigDatabase  (rtDbConnection dbConnection, rtDbCq *dbCq);
vltINT32       rtStrToType       (char *string, rtDeType deType, void *buffer);
vltINT32       rtWriteDatabase   (rtDbConnection dbConnection, rtDbReadWrite *dbRw);
vltINT32       rtMultiReadDb     (rtDbConnection dbConnection, rtDbReadWrite *rw, 
                                  vltUINT32 blkCnt);
vltINT32       rtMultiWriteDb    (rtDbConnection dbConnection, rtDbReadWrite *rw, 
                                  vltUINT32 blkCnt);
vltINT32       rtUnitWriteDb     (rtDbConnection dbConnection, rtDbReadWrite *rw, 
                                  vltUINT32 blkCnt);
vltINT32       rtCoerceDataElems (const void *oldData, rtDeType oldType, 
                                  void *newData, rtDeType newType);

extern vltUINT16  rtDataElemSize[];
;

/* These functions are for internal use ONLY and MUST NOT BE CALLED BY A PROCESS DIRECTLY !!!! */
ccsCOMPL_STAT dbAttachEvent  (rtDbConnection dbConnection, rtDbAddress *addr, rtHangerStatus status, 
                              vltUINT8 enable, vltUINT16 pIndex, vltUINT16 iIndex, 
                              dbDIRADDRESS *address, vltINT32 *index, dbATTRTYPE *type);
ccsCOMPL_STAT dbEnableEvent  (rtDbConnection dbConnection, dbDIRADDRESS *dbAddress, vltINT32 index);
ccsCOMPL_STAT dbDisableEvent (rtDbConnection dbConnection, dbDIRADDRESS *dbAddress, vltINT32 index);
ccsCOMPL_STAT dbDetachEvent  (rtDbConnection dbConnection, dbDIRADDRESS *dbAddress, vltINT32 index);
vltLOGICAL    dbIsOpen();
void         *dbGetShmPtr();
ccsCOMPL_STAT dbGetCwpPlin(vltUINT16 *plin);
ccsCOMPL_STAT dbTriggerCe (vltUINT16 *plins, vltUINT16 plinCnt);
#if defined(i386) || defined(__i386)
ccsCOMPL_STAT dbSwapBytesInBuffer ( dbTYPE *dataType, char *buffer, int size);
#endif /* i386 stuff */
#endif
