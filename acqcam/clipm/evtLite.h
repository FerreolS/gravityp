#ifndef EVT_LITE_H
#define EVT_LITE_H
/*******************************************************************************
* E.S.O. - VLT project
*
* "@(#) $Id: evtLite.h 35831 2000-06-13 12:21:57Z bgilli $"
*
* who       when      what
* --------  --------  ----------------------------------------------
* vlt       18/06/98  created
* bgilli    06/06/00  moved evtMSG_STRUCT & evt CMDs to evtPrivate.h
*/

/************************************************************************
 *
 *----------------------------------------------------------------------
 */


typedef enum
    {
    rtSYSTEM_EVENT=0,
    rtAPPLICATION_EVENT,
    rtDATABASE_EVENT=3    /* supported only */
    } rtEventType;

#define rtHANGER_INACTIVE               0
#define rtHANGER_ACTIVE                 1
typedef rtUInt16 rtHangerStatus;

#define rtDB_LESS               1
#define rtDB_EQUAL              2
#define rtDB_GREATER            4
#define rtDB_NO_EVENT           0
#define rtDB_GREATER_OR_EQUAL   (rtDB_GREATER+rtDB_EQUAL)
#define rtDB_NOT_EQUAL          (rtDB_GREATER+rtDB_LESS)
#define rtDB_LESS_OR_EQUAL      (rtDB_LESS+rtDB_EQUAL)
#define rtDB_ANY_WRITE          (rtDB_GREATER+rtDB_LESS+rtDB_EQUAL)
#define rtDB_VECTOR             7
#define rtDB_TABLE              7

typedef struct
    {
    rtProcessEnv    env;
    rtInt32         event;
    } rtEventId;

typedef struct
    {
    rtEventId       id;
    rtEventType     type;
    rtUInt8         enable;
    rtUInt8         trigger;
    rtUInt8         reserved[2];
    rtInt           triggerUId;
    rtProcessId     triggerSource;
    rtProcessId     hangerSource;
    rtUInt16        eventSize;      /* of rtEventInfo */
    rtUInt16        hangerSize;
    } rtEventHeader;

typedef struct
    {
    rtDeType        deType;
    rtBitfld        oldQuality:2;   /* rtDbQuality */
    rtBitfld        newQuality:2;   /* rtDbQuality */
    rtBitfld        reserved:4;     /* rtUInt8 */
    rtUInt8         reserved1[2];
    rtUInt8         values [512];
    } rtEventScalarInfo;
 
typedef struct
    {
    rtUInt16        startElement;
    rtUInt16        endElement;
    } rtEventVectorInfo;
 
typedef struct
    {
    rtUInt16        startRecord;
    rtUInt16        endRecord;
    rtUInt8         startField;
    rtUInt8         endField;
    rtUInt8         reserved[6];
    } rtEventTableInfo;

typedef struct
    {
    rtPlin                  plin;
    rtAin                   ain;
    rtAttributeType         type;
    union
       {
        rtEventScalarInfo       scalar;
        rtEventVectorInfo       vector;
        rtEventTableInfo        table;
        } svt;
    }       rtEventDbInfo;

typedef union
    {
    rtEventDbInfo   db;
    rtUInt8         applic; 
    } rtEventInfo;

typedef rtUInt8 rtHangerMsg;


/* Functions */


#endif /*!EVT_LITE_H*/
