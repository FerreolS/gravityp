#ifndef CCS_LITE_H
#define CCS_LITE_H
/*******************************************************************************
* E.S.O. - VLT project
*
* "@(#) $Id: ccsLite.h 266574 2015-03-19 15:08:02Z tebert $"
*
* who       when      what
* --------  --------  ----------------------------------------------
* tebert    29/04/98  created
* tebert    21/05/99  removed unused parameters from ccsDB_DATA
* tebert    21/05/99  added cwp to structure ccsPROCESS_TABLE
* tebert    22/05/99  added functions ccsSetCwp() and ccsGetCwp()
* bgilli    01/07/99  added rtE_CE_TYPE error mnemonic
* bgilli    05/07/99  Merged new errors for T. Ebert.
* tebert    14/07/99  added rtGetPid()
* tebert    02/10/99  added support of LINUX
* bgilli    06/10/99  added rtE_DB_NOT_LOCKED  
* bgilli    23/11/99  cleaned warnings revealed by new gcc compiler. 
* bgilli    18/01/00  Some more cleaning: (void) functions without params.
* tebert    02/03/00  Added new error codes for CE
* bgilli    15/03/00  Added byte swapping macros for Linux endianness.
* tebert    18/07/00  Added CCSHISMANAGER
* tebert    14/08/00  Added restart flag to environment table
* tebert    05/09/00  Added ccsBaseName()
* tebert    05/09/00  Added ccsQueueMon
* tebert    05/03/01  Added management structures for alarm system
* bgilli    07/02/02  Added CCS_SWAP_FLOAT & CCS_SWAP_DOUBLE; courtesy ALO
* bgilli    05/07/02  Removed special Linux time related stuff.
* tebert    04/11/02  ccsSINGLE_ERROR structure size of progname corrected
* tebert    03/05/11  added internal function ccsForceRead()
* tebert    15/01/13  merged with multithreading code VLTSW2011
* tebert    30/01/13  added ccsShmAttachReadOnly() (VLTSW-9587)
* tebert    18/04/13  VLTSW-9605
*/

/************************************************************************
 *
 *----------------------------------------------------------------------
 */

#define ccsMAX_PROCS     255
#define ENVNAME          "RTAPENV"
#define ENVTBLNAME       "CcsEnvTable"
#define CCSSCHEDULER     "ccsScheduler"
#define CCSQSERVER       "qsemu"
#define CCSPIPEDIR       "/tmp/"
#define CCSEVENTMANAGER  "evtEventConfg"
#define CCSDBMANAGER     "dbMQDBM"
#define CCSSCANMANAGER   "scanMngr"
#define CCSHISMANAGER    "hisDHMngr"
#define ENVTBLLIST       "CcsEnvList"
#define CCSJANITOR       "ccsJanitor"
#define CCSQUEUEMON      "ccsQueueMon"

#define ccsMSG_QUERY_TYPE      1
#define ccsMSG_QUERY_PROCNAME  12000
#define ccsMSG_QUERY_PROCNUM   13000
#define ccsMSG_QUERY_PID       14000
#define ccsMSG_ACK_TYPE        10000
#define ccsMSG_REGISTER_TYPE   11000
#define ccsMSG_CMD_TYPE        12000
#define ccsMSG_SCHEDULE_CMD    15000

#define ccsQSERVER_PROCNUM     3
#define ccsTIMEOUT_ACK         9000
#define ccsTIMEOUT             5000

#define ccsSCHED_CMD_ENTER      "ENTER"
#define ccsSCHED_CMD_NEXT_STEP  "NEXT" 
#define ccsSCHED_CMD_SHUTDOWN   "SDOWN"
#define ccsSCHED_CMD_EXIT       "EXIT"
#define ccsSCHED_CMD_RESTART    "SREST"
#define ccsSCHED_CMD_SCHEDULE   "SCHED"

#define ccsOBI_CARE_ON     0x01
#define ccsOBI_CARE_OFF    0x02
#define ccsOBI_CLEANUP_ON  0x04
#define ccsOBI_CLEANUP_OFF 0x08

#define rtHPCompatibility  ccsCompatibility
#define rtRegisterForDebug ccsRegister
#define rtSetMyName        ccsSetMyName
#define rtGetMyName        ccsGetMyNameInt
#define rtGetMyProcid      ccsGetMyProcid
#define rtGetProcName      ccsGetProcNameInt
#define rtGetProcNum       ccsGetProcNumInt
#define rtCompareProcids   ccsCompareProcids
#define rtGetPidByName     ccsGetPidByName
#define rtMsgSend          ccsMsgSend
#define rtMsgRecv          ccsMsgRecv
#define rtMonitorQueue     ccsMonitorQueue
#define rtScheduleProcess  ccsScheduleProcess
#define rtPrintError       ccsPrintOldError
#define rtSingleError      ccsSINGLE_ERROR
#define rtErrno            ccsErrno
#define rtError            ccsOLD_ERROR
#define rtStartNextPhase   ccsStartNextPhase
#define rtProcessNum       ccsPROCNUM
#define rtGetUid           ccsGetUid 

#define rtRESPOND       (vltUINT8) 2       /* acknowledgement desired */
#define rtNO_RESPOND    (vltUINT8) 1       /* don't acknowledge */
#define rtRESPONSE      (vltUINT8) 0       /* acknowledgement message */

#define rtACKNOWLEDGE           64
#define rtTK_MESSAGE            65
#define rtEVENT_MESSAGE         66
#define rtPS_OBITUARY           67
#define rtASCII_MESSAGE         71
#define rtUSER_MESSAGE          128

#define rtQ_SERVER_PING_REQUEST   37
#define rtQ_SERVER_PING_RESPONSE  38
#define rtQ_SERVER_ACKNOWLEDGE    39

#define rtNO_WAIT       (vltINT32) -1      /* don't wait if no message */
#define rtWAIT_SIG      (vltINT32) -2      /* wait unless signal recv'ed*/
#define rtNO_TIMEOUT    (vltINT32) -3      /* wait forever for message */

#define rtMINPRIO       (vltINT32)100      /* lowest priority allowed */

#define rtBLOCKING      (vltUINT8) 0        /* block on full/empty queue */
#define rtNO_BLOCK      (vltUINT8) 1       /* fail on full/empty queue */

#define rtNO            (rtInt)'N'
#define rtPROGRAM       (rtInt)'P'
#define rtENVIRONMENT   (rtInt)'E'

#define rtSUCCESS       0
#define rtFAILED       -1

#define rtREMOTE_QUERY_PROCNUM      11
#define rtREMOTE_QUERY_PROCNAME     12
#define rtREMOTE_PROCNUM            13
#define rtREMOTE_PROCNAME           14
#define rtREMOTE_QUERY_PID          24
#define rtREMOTE_PID                25

#define rtMAX_DE_TO_STR_CHARS ((256 * 4) + 1)
#define rtMSGSIZE             8192

#define rtMERGE_RUNSTRING   0x04

#define rtRTAP_QUEUE_SERVER "qsemu"
#define ccsDEAD_KEY         "VLTSW_9605"

#if defined(i386) || defined(__i386)
/*
 * byte swapping macros needed for LINUX on PCs
 */

/* most significant byte of 2-byte integer */
#define CCS_MSB_INT16(x)        (((x) >>  8) & 0xff)
/* least significant byte of 2-byte integer */
#define CCS_LSB_INT16(x)         ((x)        & 0xff)
/* most significant word of 2-word integer */
#define CCS_MSW_INT32(x) (((x) >> 16) & 0xffff)
/* least signifcant byte of 2-word integer */
#define CCS_LSW_INT32(x)  ((x)        & 0xffff)

/* swap the MSB with the LSB of a 16 bit integer */
#define CCS_SWAP_INT16(x) (CCS_MSB_INT16(x) | (CCS_LSB_INT16(x) << 8))

/* swap all bytes of a 32 bit integer */
#define CCS_SWAP_INT32(x) (CCS_SWAP_INT16(CCS_LSW_INT32(x)) << 16) | \
                             CCS_SWAP_INT16(CCS_MSW_INT32(x))
/*void ccs_swap_2int32(int32_t* i1, int32_t* i2)
  { int32_t i=*i1; *i1=*i2; *i2=i; };*/
                          
/* swap all bytes of two consecutive 32 bit integers */
#define CCS_SWAP_2INT32(x_ptr,y_ptr) {int32_t i=*x_ptr;\
                                     *x_ptr=*y_ptr;\
                                     *y_ptr=i;}\
                           *x_ptr = CCS_SWAP_INT32(*x_ptr);         \
                           *y_ptr = CCS_SWAP_INT32(*y_ptr) 
/* courtesy alongino */			   
#define CCS_SWAP_FLOAT(v)  {vltINT32 *i1; i1=(vltINT32 *)&(v); \
                            *i1 = CCS_SWAP_INT32(*i1);} 
#define CCS_SWAP_DOUBLE(v) {vltINT32 *i1,*i2; i1=(vltINT32 *)&(v); \
                            i2=i1+1;CCS_SWAP_2INT32(i1,i2);} 
#endif

#define CCS_ALIGN(n) (((long)(n+7)) & ~7)

#ifndef rtTRUE
#define rtTRUE  (1)
#endif

#ifndef rtFALSE
#define rtFALSE (0)
#endif

typedef vltUINT16      rtUInt16;
typedef vltINT16       rtInt16;
typedef vltINT8        rtInt8;
typedef vltUINT8       rtUInt8;
typedef char           rtChar;
typedef vltINT32       rtInt;
typedef vltUINT32      rtUInt32;
typedef vltINT32       rtInt32;
typedef ccsENVNAME     rtProcessEnv;
typedef unsigned int   rtBitfld;
typedef vltDOUBLE      rtDouble;
typedef vltFLOAT       rtFloat;
typedef vltLOGICAL     rtLogical;
typedef vltBYTES4      rtBytes4;
typedef vltBYTES8      rtBytes8;
typedef vltBYTES12     rtBytes12;
typedef vltBYTES16     rtBytes16;
typedef vltBYTES20     rtBytes20;
typedef vltBYTES32     rtBytes32;
typedef vltBYTES48     rtBytes48;
typedef vltBYTES64     rtBytes64;
typedef vltBYTES80     rtBytes80;
typedef vltBYTES128    rtBytes128;
typedef vltBYTES256    rtBytes256;
typedef vltPOLAR       rtPolar;
typedef vltRECTANGULAR rtRectangular; 
typedef ccsTIMEVAL     rtTime;

typedef struct
    {
    ccsENVNAME    procEnv;        /* environment name */
    ccsPROCNUM    procNum;        /* local process number */
    vltUINT8      reserved;
    } rtProcessId;

#define rtANYBODY        (rtProcessId *)0       /* receive msg from anyone */
#define rtANY_RESPONSE   (-ccsMAX_PROCS)         /* receive any response */
#define rtCHRON_ORDER    (vltINT32)0
#define rtPRIOR_ORDER    (vltINT32)-rtMINPRIO

typedef struct
    {
    rtUInt16 currentBytes;
    rtUInt16 queueSize;
    rtUInt16 lastSndPid;
    rtUInt32 lastSndTime;
    rtUInt32 lastRcvTime;    
    rtUInt32 lastChangeTime;
    } rtMsgQstats;

typedef struct ccsNode		/* Node of a linked list. */
    {
    struct ccsNode *next;	/* Points at the next node in the list */
    struct ccsNode *previous;	/* Points at the previous node in the list */
    } ccsNODE;

typedef struct			/* Header for a linked list. */
    {
    ccsNODE node;	       	/* Header list node */
    int     count;		/* Number of nodes in list */
    } ccsLIST;

typedef struct
    {
    ccsNODE  node;
    int      size;
    void    *address;
    } ccsMEM_BLK;

typedef struct
    {
    int        shmId;
    int        semId;
    vltLOGICAL status;
    } ccsSCAN_DATA;

typedef struct
    {
    int        shmId;
    int        dbConfSem;
    int        dbWriteSem;    
    void      *shmAddr;
    void      *mapRootPtr;
    void      *mapTable;
    void      *dbShmVars;
    ccsLIST    usedBlocks;
    ccsLIST    freeBlocks;
    vltLOGICAL status;
    } ccsDB_DATA;

typedef struct
    {
    int        cntIn;
    int        cntOut;
    int        queueSize;
    void      *shmAddr;
    int        evtSemId;
    int        triggerSem;
    vltLOGICAL status;
    } ccsEVT_DATA;

typedef struct
    {
    int        cntIn;
    int        cntOut;
    int        queueSize;
    void      *shmAddr;
    int        triggerSem;
    vltLOGICAL status;
    } ccsALRM_DATA;

typedef struct 
    {
    ccsPROCNAME procName;
    vltLOGICAL  fix;
    vltLOGICAL  status;
    vltLOGICAL  obi;
    vltLOGICAL  care;
    vltUINT8    restart;
    ccsTIMEVAL  timeStamp;
    int         msgId;
    int         pid;
    int         uid;
    int         gid;
    int         monSem;
    vltINT32    pipeFd[2];
    vltINT32    monPid;
    int         janitorSem;
    void       *cwp;
    ccsPROCNAME oldName;
    } ccsPROCESS_TABLE;

typedef struct
    {
    ccsPROCESS_TABLE proc[ccsMAX_PROCS];
    ccsDB_DATA       db;
    ccsEVT_DATA      evt;
    ccsSCAN_DATA     scan;
    ccsALRM_DATA     alrm;
    int              ccsMonShmId;
    } ccsGLOBAL_SHM;

typedef struct
    {
    vltUINT16   phase;
    vltUINT16   index;
    ccsPROCNUM  procNum;
    char        command[512];
    } ccsPROC_AUTO;

typedef struct
    {
    ccsPROCNAME procName;
    int         msgId;
    int         pid;
    int         uid;
    int         gid;
    vltUINT32   obi;
    } ccsENTER;

typedef struct  
    {
    vltINT16      sourceId;
    ccsPROCNUM    procNum;
    vltUINT8      reserved;
    vltINT32      errorNumber;
    char          sourceName[20];
    ccsENVNAME    envName;
    char          progName[ccsPROCNAME_LEN+1];
    char          functionName[80];
    char          errorString[128];
    char          messageString[256];
    char          location[80];
    } ccsSINGLE_ERROR;

typedef struct 
    {
    vltUINT16        errorCount;
    vltUINT16        ignoreError;
    vltUINT32        errorTime;
    ccsSINGLE_ERROR  errorList[10];
    } ccsOLD_ERROR;

typedef struct 
    {
    ccsOLD_ERROR     error;
    } ccsACKNOWLEDGE;

typedef struct 
    {
    ccsENVNAME  env;
    ccsPROCNAME proc;
    } ccsENVPROC_NAME;

typedef struct
    {
    vltINT32     priority;              /* Internal use only                 */
    rtProcessId  source;                /* Sender of the message             */
    rtProcessId  dest;                  /* Recipient                         */
    vltINT32     sourceUId;             /* User Id                           */
    vltUINT16    msgId;                 /* To distinguish identical commands */
    vltUINT8     msgType;               /* Defined by message system         */
    vltUINT8     responseFlag;          /* Command/Reply                     */
    vltINT32     typeHdrSize;    
    vltINT32     msgBodySize;
    } rtMsgHeader;

typedef struct  {
    vltINT32       exitStatus;          /* how it terminated                 */
    rtProcessId    deceased;            /* process that died                 */
                                        /* deceased proc_num : positive if   */
                                        /* process obituary, 0 if environment*/
                                        /* obituary                          */
    ccsPROCNAME    procName;            /* process Name of the dead process  */
    vltINT8        restart;             /* will it be restarted              */
    vltINT8        debug;               /* did it register for debug         */
    vltUINT8       reserved[4];
    } rtObituary;

typedef struct
    {
    rtMsgHeader msgHeader;
    union
        {
        char            str[rtMSGSIZE - sizeof (rtMsgHeader)];
        ccsACKNOWLEDGE  ack;
        rtObituary      obi;
        int             param[1];
        } body;
    } ccsMSG_ACK;

typedef struct
    {
    ccsTIMEVAL      time;           /* makes ID unique */
    rtProcessId     procid;         /* who made request */
    vltUINT8        reserved[2];
    } rtTimerId;

typedef struct
    {
    rtTimerId       id;             /* request id */
    ccsTIMEVAL      time;           /* when alarm expired */
    vltUINT8        reserved[4];
    } rtTimerMessage;

#ifndef CCS_PTHREADS
extern ccsOLD_ERROR  ccsErrno;
#else  
extern __thread ccsOLD_ERROR  ccsErrno;
#endif

int            ccsGetSemaphore      (int num);
int            ccsGetSemaphoreEmpty (int num);
ccsCOMPL_STAT  ccsSemTake           (int semid, int num, int waitx);
ccsCOMPL_STAT  ccsSemTakeNu         (int semid, int num, int waitx);
ccsCOMPL_STAT  ccsSemGive           (int semid, int num);
ccsCOMPL_STAT  ccsSemGiveNu         (int semid, int num);
ccsCOMPL_STAT  ccsSemDelete         (int semid);
int            ccsGetShm            (int size);
int            ccsGetCcsShmId       (void);
ccsCOMPL_STAT  ccsShmAttach         (int shmId, void **addr);
ccsCOMPL_STAT  ccsShmAttachReadOnly (int shmId, void **addr);
ccsCOMPL_STAT  ccsShmDetach         (void *addr);
ccsCOMPL_STAT  ccsShmDelete         (int shmId);
int            ccsGetMsgId          (int key);
void           ccsDelMsgId          (int key);
ccsCOMPL_STAT  ccsReadNamedPipe     (char *msg);
void           ccsUnregister        (ccsERROR *error);
void           ccsForceRead         (void);

void           ccsCompatibility  (int   argc, char **argv);
rtProcessId   *ccsRegister       (char *env,  char  *name);
vltINT32       ccsSetMyName      (char *procName);
rtProcessId   *ccsGetMyProcid    (void);
char          *ccsGetMyNameInt   (void);
char          *ccsGetProcNameInt (rtProcessId *processId);
ccsPROCNUM     ccsGetProcNumInt  (char *env, char *name);
int            ccsCompareProcids (rtProcessId *id1, rtProcessId *id2);
int            ccsGetPidByName   (char *procName);
vltINT32       ccsMsgSend        (vltUINT8 blockCntrl, void *message, vltUINT32 msgSize);
rtMsgHeader   *ccsMsgRecv        (vltINT32 timeout, rtProcessId *source, vltINT32 priority);
int            ccsMonitorQueue   (int state);
ccsCOMPL_STAT  ccsShutdown       (ccsENVNAME envName);
rtProcessId   *ccsScheduleProcess(char *env, char *runString, int flags);
void           ccsPrintOldError  (FILE *fp);
vltINT32       ccsStartNextPhase (void);
int            ccsGetUid         (void);
ccsCOMPL_STAT  ccsSetCwp         (void *node);
void          *ccsGetCwp         (void);
void           ccsSetObiStatus   (vltUINT32 obi);
void           ccsSetCoreProc    (void);
char          *ccsGetDeadName    (void);
rtInt          rtGetQueueStats   (rtProcessId *procId, rtMsgQstats *mqStats);
rtInt          rtGetPid          (rtChar *env, rtChar *procName);
void           rtSetError        (rtChar *funcName, rtInt errorNum, rtChar *message);
void           rtErrorContext    (rtChar *funcName, rtInt errorNum, rtChar *message);
void           rtLogError        (rtInt  printError, rtChar *message);
rtInt          rtSetRestart      (rtInt  type);

char         *ccsInternalFindPath(char *fileName);
ccsCOMPL_STAT ccsInternalTestPath(char *fileName);
char         *ccsBaseName        (char *name);

int ccsMemInit(void *addr, int size, ccsLIST *usedBlocks, ccsLIST *freeBlocks);
void *ccsMemCalloc(ccsLIST *usedBlocks, ccsLIST *freeBlocks, int elements, int elsize);
void *ccsMemRealloc(ccsLIST *usedBlocks, ccsLIST *freeBlocks, void *oldPtr, int size);
void ccsMemFree(ccsLIST *usedBlocks, ccsLIST *freeBlocks, void *address);
void ccsMemStatistics (ccsLIST *usedBlocks, ccsLIST *freeBlocks,vltLOGICAL print);

void             lstInit         (ccsLIST *pList);
void             lstAdd          (ccsLIST *pList, ccsNODE *pNode);
int              lstCount        (ccsLIST *pList);
void             lstDelete       (ccsLIST *pList, ccsNODE *pNode);
ccsNODE         *lstFirst        (ccsLIST *pList);
ccsNODE         *lstLast         (ccsLIST *pList);
ccsNODE         *lstGet          (ccsLIST *pList);
void             lstInsert       (ccsLIST *pList, ccsNODE *pPrev, ccsNODE *pNode);
ccsNODE         *lstNext         (ccsNODE *pNode);
ccsNODE         *lstPrevious     (ccsNODE *pNode);
void             lstFree         (ccsLIST *pList);


#define rtE_NAME_TOO_LONG           0x00300007
#define rtE_ENV_VAR                 0x00000012
#define rtE_NOENV                   0x0030000A
#define rtE_NO_RTAPDIR              0x0030000C
#define rtE_INVLD_RTAPDIR           0x0030000D
#define rtE_DB_INV_PARAMETER        0x00520040
#define rtE_NOT_IMPL                0x0000000F
#define rtE_DB_INV_BUFF             0x00520009
#define rtE_BUFF_SIZE               0x00000014
#define rtE_DB_SMALL_BUFF           0x0052000D
#define rtE_DATA                    0x00000017
#define rtE_UTL_DE_COERCE           0x00200009
#define rtE_DB_MATCH_NULL_PTR       0x00520050
#define rtE_NOLOC                   0x00300008
#define rtE_NORMT                   0x00300009
#define rtE_REMOTE_DOWN             0x00100006
#define rtE_REMOTE_HOST             0x00100007
#define rtE_REMOTE_ENV              0x00100008
#define rtE_NO_QUEUE_SERVER         0x00100009
#define rtE_ENV_NOT_ACTIVE          0x0010000A
#define rtE_SYSERR                  0x00000011
#define rtE_SS_SSMA_NOT_OPEN        0x00900000
#define rtE_DB_NOT_OPEN             0x00520000
#define rtE_DB_INV_CONN             0x00520001
#define rtE_DB_INV_NAME             0x00520034
#define rtE_DB_SYM_SYNTAX           0x0052000F
#define rtE_DB_ADDR_NSF             0x00520010
#define rtE_DB_ADDR_2MUCH           0x00520011
#define rtE_DB_DUP_HNAME            0x00520018
#define rtE_DB_DUP_ALIAS            0x00520019
#define rtE_DB_INV_ALIAS            0x00520036
#define rtE_DB_INV_POINT            0x00520005
#define rtE_DB_DUP_ATTR_NAME        0x00520046
#define rtE_DB_INV_ATTR             0x00520006
#define rtE_DB_INV_ATYPE            0x00520002
#define rtE_DB_INV_ADDR             0x00520003
#define rtE_DB_BUFF_ADDR            0x00520020
#define rtE_DB_WR_ATTY              0x00520028
#define rtE_DB_WR_DETY              0x00520029
#define rtE_DB_INV_DE_TYPE          0x00520047
#define rtE_DB_INV_RECEL            0x00520007
#define rtE_DB_WR_RECEL             0x0052002A
#define rtE_DB_INV_FIELD            0x00520008
#define rtE_DB_INV_ACTION           0x00520012
#define rtE_DB_ATTR_TYPE            0x00520021
#define rtE_DB_INV_UNLOCK           0x00520024
#define rtE_DB_NOT_LOCKED           0x00520025
#define rtE_DB_GROUP                0x00520027
#define rtE_DB_CATEGORY             0x0052000B
#define rtE_DB_NO_PLINS             0x00520014
#define rtE_DB_REUSE_PLIN           0x00520015
#define rtE_DB_NO_ALIASES           0x00520017
#define rtE_DB_SHM_RAM_V            0x0052001A
#define rtE_DB_SHM_RAM_H            0x0052001B
#define rtE_DB_DUP_FLD_NAME         0x00520049
#define rtE_DB_CNFG_ROOT            0x00520016
#define rtE_DB_DEL_CHILD            0x00520037
#define rtE_DB_DEL_DEPEND           0x00520038
#define rtE_DB_DEL_DEFN             0x00520039
#define rtE_DB_DEL_EVENT            0x0052003A
#define rtE_DB_INV_RESID            0x0052001C
#define rtE_DB_BAD_BLKCNT           0x0052004B
#define rtE_SEND_MSG                0x00000009
#define rtE_RECV_MSG                0x0000000A
#define rtE_MSG_TYPE                0x0000000B
#define rtE_MSGOVFL                 0x00100000
#define rtE_MSGSIZE                 0x00100002
#define rtE_MSGTO                   0x00100003
#define rtE_NOQUEUE                 0x00100001
#define rtE_MSG_Q_FULL              0x0010000B
#define rtE_MSG_Q_EMPTY             0x0010000C
#define rtE_MSG_DEST                0x0010000E
#define rtE_SEC_NO_FILE             0x00B00000
#define rtE_SEC_NO_DIR              0x00B00001
#define rtE_EVENT_BAD_ADDR          0x00800001
#define rtE_EVENT_NOEM              0x0080000A
#define rtE_EVENT_TOOBIG            0x00800007
#define rtE_EVENT_SEND_FAIL         0x00800010
#define rtE_EVENT_NOT_CONFIG        0x00800005
#define rtE_EVENT_NOMEM             0x0080000C
#define rtE_EVENT_NONE              0x00800000
#define rtE_ARG                     0x00000016
#define rtE_TOOBIG                  0x00400004
#define rtE_BADTIME                 0x00400005 
#define rtE_NOID                    0x00400006
#define rtE_EXISTS                  0x00400007
#define rtE_UTL_BADTIME             0x00200008
#define rtE_UTL_BAD_VAL             0x00200004
#define rtE_CE_TYPE                 0x0070001F
#define rtE_PARSE_EXPR              0x00600006
#define rtE_PARSE_FUNCTION          0x0060000C
#define rtE_DB_CFI_SET              0x00520043
#define rtE_DB_CFI_NOTSET           0x00520044
#define rtE_DB_CFI_FAILURE          0x00520045
#define rtE_DB_ADDR_CNT             0x0052001D
#define rtE_PARSE_SYNTAX            0x00600002
#define rtE_PARSE_INTERNAL          0x00600004
#define rtE_PARSE_EXPR              0x00600006
#define rtE_PARSE_TOO_MANY          0x0060000A
#define rtE_PARSE_TOO_FEW           0x00600009
#define rtE_PARSE_INPUT             0x00600003
#define rtE_PARSE_BRACKET           0x00600008
#define rtE_PARSE_OPERATION         0x00600001
#define rtE_PARSE_STATE             0x0060000B
#define rtE_PARSE_FUNCTION          0x0060000C
#define rtE_PARSE_NO_TYPE           0x0060000D
#define rtE_PARSE_TABLE             0x0060000E
#define rtE_PARSE_STRING            0x0060000F
#define rtE_PARSE_SCAN              0x00600010
#define rtE_PARSE_ELEMENT           0x00600011
#define rtE_PARSE_EMPTY             0x00600012
#define rtE_PARSE_LAST              0x00600014
#define rtE_SS_INVLD_SCAN_LINK      0x00900017
#define rtE_SS_NO_SCAN_CNFG_PT      0x00900012
#define rtE_PARSE_NO_TYPE           0x0060000D
#define rtE_PARSE_SYNTAX            0x00600002
#define rtE_PARSE_ELEMENT           0x00600011
#define rtE_DB_INV_CE_OPER          0x0052003B
#define rtE_SS_NO_SCAN_MNGR         0x0090000A
#define rtE_SS_INVLD_PARAMETER      0x00900004
#define rtE_SS_INVLD_SSMA_OFFSET    0x00900001 
#define rtE_SS_INVLD_SSMA_QUALITY   0x00900003
#define rtE_SS_INVLD_BUFFER_SIZE    0x00900010
#define rtE_SS_IMPROPER_SM_RESP     0x00900015
#define rtE_SS_INVLD_CONNECTION     0x0090000C
#define rtE_SS_INVLD_ACTION         0x0090000F
#define rtE_SS_INVLD_ADDR_COUNT     0x0090000E
#define rtE_SS_INVLD_ADDRESS        0x00900011
#define rtE_SS_REQUEST_FAILED       0x00900013
#define rtE_SS_NO_ST_RESPONSE       0x00900022
#define rtE_DB_EXPR_ORDER           0x0052003C
#define rtE_RESTRT                  0x00300004
#define rtE_NOCHG                   0x00300005
#define rtE_MAXPROC                 0x00300002
#define rtE_RUNNING                 0x00300000 

#endif /*!CCS_LITE_H*/
