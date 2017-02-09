#ifndef evhStates_H
#define evhStates_H
/*******************************************************************************
* E.S.O. - VLT project
*
* "@(#) $Id: evhStates.h 8941 1997-04-17 11:40:16Z gchiozzi $"
*
* who       when      what
* --------  --------  ----------------------------------------------
* gchiozzi  17/04/97  Patched to be used without including ccs.h
* gchiozzi  22/05/96  Removed c++ style comments
* gchiozzi  23/05/95  created
*/

/************************************************************************
 *
 * This include file contains standard definitions for state codes and 
 * names.
 * It is used both by C++ code and by dbLoader database definitions
 *----------------------------------------------------------------------
 */

/************************************
 * Defines for commonly used states *
 ************************************/

/* This should be replaced by */
/* #include "ccsStates.h"     */

#ifndef ccsSTATE_UNK
#define    ccsSTATE_UNK       0
#define    ccsSTATE_OFF       1
#define    ccsSTATE_LOADED    2
#define    ccsSTATE_STANDBY   3
#define    ccsSTATE_ONLINE    4
#define    ccsSTATE_STANDALONE 11
#endif

#define    evhSTATE_UNK        ccsSTATE_UNK
#define    evhSTATE_OFF        ccsSTATE_OFF
#define    evhSTATE_LOADED     ccsSTATE_LOADED
#define    evhSTATE_STANDBY    ccsSTATE_STANDBY
#define    evhSTATE_ONLINE     ccsSTATE_ONLINE
#define    evhSTATE_READY      5
#define    evhSTATE_ERROR      6
#define    evhSTATE_IDLE       7
#define    evhSTATE_WAITING    8
#define    evhSTATE_BUSY       9
#define    evhSTATE_Q_CHECK    10
#define    evhSTATE_STANDALONE ccsSTATE_STANDALONE

/******************************************
 * Defines for commonly used states names *
 ******************************************/

#define    evhSTATE_STR_UNK       "UNKNOWN"
#define    evhSTATE_STR_OFF       "OFF"
#define    evhSTATE_STR_LOADED    "LOADED"
#define    evhSTATE_STR_STANDBY   "STANDBY"
#define    evhSTATE_STR_ONLINE    "ONLINE"
#define    evhSTATE_STR_READY     "READY"
#define    evhSTATE_STR_ERROR     "ERROR"
#define    evhSTATE_STR_IDLE      "IDLE"
#define    evhSTATE_STR_WAITING   "WAITING"
#define    evhSTATE_STR_BUSY      "BUSY"
#define    evhSTATE_STR_Q_CHECK   "QUEUE_CHECK"

/*
 * evhDB_TASK.C contains the fndNAME_AND_INDEX table for name-index conversion
 * The table must be updated whenever these definitions change
 */

#endif /*!evhStates_H*/
