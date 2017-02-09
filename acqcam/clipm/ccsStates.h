/*************************************************************************
* E.S.O. - VLT project
*
* "@(#) $Id: ccsStates.h 15267 1998-03-19 08:29:42Z mcomin $" 
*
* ccsStates.h  -  
*
* who        when      what
* ---------  --------  ----------------------------------------------
* mcomin     05/02/98  Creation
*            18/03/98  Changed definition of  #ifndef CCS_STATES
*                  
****************************************************************************/

#ifndef CCS_STATES
#define CCS_STATES

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Defines for commonly used states
 */
#define    ccsSTATE_UNK       0
#define    ccsSTATE_OFF       1
#define    ccsSTATE_LOADED    2
#define    ccsSTATE_STANDBY   3
#define    ccsSTATE_ONLINE    4
#define    ccsSTATE_STANDALONE 11

/*
 * Defines for commonly used states names
 */

#define    ccsSTATE_STR_UNK           "UNKNOWN"
#define    ccsSTATE_STR_OFF           "OFF"
#define    ccsSTATE_STR_LOADED        "LOADED"
#define    ccsSTATE_STR_STANDBY       "STANDBY"
#define    ccsSTATE_STR_ONLINE        "ONLINE"
#define    ccsSTATE_STR_STANDALONE    "STANDALONE"

#ifdef __cplusplus
}
#endif

#endif /*!CCS_STATES*/
