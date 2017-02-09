#ifndef gvacqDEFINES_H
#define gvacqDEFINES_H
/*******************************************************************************
* E.S.O. - VLT project
*
* "@(#) $Id: gvacqDefines.h 293324 2017-02-06 14:33:13Z tott $"
*
*/
/** \internal
* who       when        what
* --------  ----------  ----------------------------------------------
* amorim    2017-02-06  added SETSTAR2
* amorim    2017-01-26  added SETSKY
* narsi     2016-03-25  added gvacqNUM_REC_ABS_SLOPES
* narsi     2016-03-19  added gvacqNUM_REC_ABS_SPOT_NUM, gvacqNUM_REC_PT_SPOT_NUM and gvacqNUM_REC_ABS_SLOPES 
* narsi     2016-03-19  added gvacqNUM_REC_DEFOCUS and gvacqNUM_REC_REF_DEFOCUS
* ekw       2016-03-11  add SAVEIMG
* ekw       2016-03-01  add TEST2
* ekw       2014-01-20  clean up
* jott      2014-01-20  added support to register at rtdcoreServer
* jott      2013-07-01  added definitions for the gvoacqProcessImageAC database
* jott      2013-05-29  created.
* ekw       2013-07-17  adding modes/commands
*/

/*#ifndef __cplusplus
#error This is a C++ include file and cannot be used from plain C
#endif
*/

/*
 * Basic application headers
 */
#include "evhStates.h"


/** Module name. */
#define gvacqMOD   "gvacq"

#define gvacqCONCAT_SVN_STRING "SVN Version " MOD_VERSION " -> compiled at " __DATE__ " "  __TIME__
#define gvacqSVN_ID()                      \
   static const char *svnId __attribute__((unused)) = gvacqCONCAT_SVN_STRING;


/** Database root point for the gvacqControl process required for GRAVITY observations */
#define gvacqDB_CONTROL_ROOT        "<alias>gvo:gvacqControl"

/** Database root point for the telescope environment */
/**define gvacqDB_TELESCOPE_ROOT      "@wsimiss:Appl_data:VLTI" */
#define gvacqDB_TELESCOPE_ROOT      ":Appl_data:VLTI"

/** Database root point for the telescope environment */
#define gvacqDB_GRAVITY_ROOT      "<alias>GRAVITY"

/** Database root point for the gvacqControl process required for GRAVITY observations */
#define gvacqDB_CONTROL_ROOT        "<alias>gvo:gvacqControl"



/** Attribute name for gvacqControl process sub-state. */
#define gvacqDB_SUBSTATE            "subState"

/** Attribute name for gvacqControl process operational mode */
#define gvacqDB_OPMODE              "opMode"

/** Database point name for the data of the gvoacqProcessImageAC library function */
#define gvacqDB_POINT_PROCIMAGEAC   "corrections"


/** Sub-state names and identifiers */
#define gvacqSUBSTATE_STR_UNK            evhSTATE_STR_UNK
#define gvacqSUBSTATE_STR_ERROR          evhSTATE_STR_ERROR
#define gvacqSUBSTATE_STR_IDLE           evhSTATE_STR_IDLE
#define gvacqSUBSTATE_STR_BUSY           evhSTATE_STR_BUSY
#define gvacqSUBSTATE_STR_INITIALIZE     "INITIALIZING"
#define gvacqSUBSTATE_STR_STOPPED        "STOPPED"

#define gvacqSUBSTATE_ID_UNK             evhSTATE_UNK     // = 0
#define gvacqSUBSTATE_ID_ERROR           evhSTATE_ERROR   // = 6
#define gvacqSUBSTATE_ID_IDLE            evhSTATE_IDLE    // = 7
#define gvacqSUBSTATE_ID_BUSY            evhSTATE_BUSY    // = 9
#define gvacqSUBSTATE_ID_INITIALIZE      20
#define gvacqSUBSTATE_ID_STOPPED         21


/**Definitions for commands supported by the gvacqControl process.*/
#define gvacqCMD_DUMMY               "DUMMY"  
#define gvacqCMD_SETMODE             "SETMODE"
#define gvacqCMD_CLRMODE             "CLRMODE"
#define gvacqCMD_GETMODE             "GETMODE"
#define gvacqCMD_ENABLE              "ENABLE"
#define gvacqCMD_DISABLE             "DISABLE"
#define gvacqCMD_ATTACH              "ATTACH"
#define gvacqCMD_DETACH              "DETACH"
#define gvacqCMD_STORSKY             "STORSKY"
#define gvacqCMD_BLINK               "BLINK"
#define gvacqCMD_ADRCORR             "ADRCORR"
#define gvacqCMD_SETHK               "SETHK"
#define gvacqCMD_SETFIW              "SETFIW"
#define gvacqCMD_RESTORE             "RESTORE"
#define gvacqCMD_TEST2               "TEST2"
#define gvacqCMD_SAVEIMG             "SAVEIMG"
#define gvacqCMD_SETSKY              "SETSKY"
#define gvacqCMD_SETSTAR2            "SETSTAR2"

/** Definitions for operational modes supported by the gvacqControl process.*/
#define gvacqOPMODE_STR_NOTSET         "NOT_SET"
#define gvacqOPMODE_STR_PUPILTRACK     "PUPILTRACK"
#define gvacqOPMODE_STR_ABERRATION     "ABERRATION"
#define gvacqOPMODE_STR_ALL            "ALL"
#define gvacqOPMODE_STR_FIELDTRACK     "FIELDTRACK"
/*ErW temporary change to start other modes as well */
#define gvacqOPMODE_STR_FIBERCOUPLERREFPOS	"FIBERCOUPLERREFPOS"
#define gvacqOPMODE_STR_CALABERRATIONSENSOR	"CALABERRATIONSENSOR"
#define gvacqOPMODE_STR_CALBARYCENTER		"CALBARYCENTER"

#define gvacqOPMODE_ID_NOTSET         0
#define gvacqOPMODE_ID_PUPILTRACK     1
#define gvacqOPMODE_ID_ABERRATION     2
#define gvacqOPMODE_ID_ALL            3
#define gvacqOPMODE_ID_FIELDTRACK     4
/*ErW temporary change to start other modes as well */
#define gvacqOPMODE_ID_FIBERCOUPLERREFPOS	5
#define gvacqOPMODE_ID_CALABERRATIONSENSOR	6
#define gvacqOPMODE_ID_CALBARYCENTER		7

#define gvacqOPMODE_PARA_NOTSET         "NOT_SET"
#define gvacqOPMODE_PARA_PUPILTRACK     "PUPILTRACK"
#define gvacqOPMODE_PARA_ABERRATION     "ABERRATION"
#define gvacqOPMODE_PARA_ALL            "ALL"
#define gvacqOPMODE_PARA_FIELDTRACK     "FIELDTRACK"
/*ErW temporary change to start other modes as well */
#define gvacqOPMODE_PARA_FIBERCOUPLERREFPOS	"FIBERCOUPLERREFPOS"
#define gvacqOPMODE_PARA_CALABERRATIONSENSOR	"CALABERRATIONSENSOR"
#define gvacqOPMODE_PARA_CALBARYCENTER		"CALBARYCENTER"

/** Definitions for event names and related database source attributes,
 *  which can be enabled or disabled via command. 
 */  
#define gvacqEVENT_DET1_IMAGE_STR    "ACQ_IMAGE"
#define gvacqEVENT_DET1_IMAGE_SRC    "<alias>ngcircon_NGCIR1:exposure.newDataFileName"
#define gvacqEVENT_DET1_IMAGE_FLAG   "evtFlagDet1Image"


/** Definitions for cameras to be attached via rtdcoreServer. */
#define gvacqCAMERA_DUMMYACQ         "DummyACQ"
#define gvacqCAMERA_ACQ              "NGCIR1"

#define ACQ_FITS_EXT "IMAGING_DATA_ACQ"

/** Definitions of the number of records for the gvoacqProcessImageAC tables */ 
#define gvacqNUM_REC_ABS_REF                154
#define gvacqNUM_REC_ABS_REF_ERROR          gvacqNUM_REC_ABS_REF
#define gvacqNUM_REC_ZERNIKE                28
#define gvacqNUM_REC_ZERNIKE_ERROR          gvacqNUM_REC_ZERNIKE
#define gvacqNUM_REC_PT_REF_POSITION        32
#define gvacqNUM_REC_PT_REF_POSITION_ERROR  gvacqNUM_REC_PT_REF_POSITION
#define gvacqNUM_REC_PT_POSITION            3
#define gvacqNUM_REC_PT_POSITION_ERROR      gvacqNUM_REC_PT_POSITION
#define gvacqNUM_REC_PT_REF_SPOT            8
#define gvacqNUM_REC_PT_REF_SPOT_ERROR      gvacqNUM_REC_PT_REF_SPOT
#define gvacqNUM_REC_PT_BARY_CENTER         8
#define gvacqNUM_REC_FI                     2
#define gvacqNUM_REC_FI_ERROR               gvacqNUM_REC_FI 
#define gvacqNUM_REC_TELESCOPES             4
#define gvacqNUM_REC_XYPOS                  2
#define gvacqNUM_REC_DEFOCUS                4
#define gvacqNUM_REC_REF_DEFOCUS            4
#define gvacqNUM_REC_ABS_SPOT_NUM           77
#define gvacqNUM_REC_PT_SPOT_NUM            16
#define gvacqNUM_REC_ABS_SLOPES            154 /* maximum 136 slopes for 9x9 lenslet array */
/** VLTSW 2011 does not use latest CPL 6.3.1. */ 
#define gvacqCPL_631  0

#endif /* gvacqDEFINES_H */
