/*******************************************************************************
* E.S.O. - VLT project
#
# "@(#) $Id: ccsSHMemory.h 266574 2015-03-19 15:08:02Z tebert $" 
*
* ccsSHMemory.h  -  Definition of the CCS Shared memory Area
*
* who        when      what
* ---------  --------  ----------------------------------------------
* M.COMIN    12/12/94  Created
* tebert     03/08/14  VLTSW-10807
*                   
****************************************************************************/
#ifndef CCSSHM_H
#define CCSSHM_H

#define  ccsSHMTest "Test SUCCESSFUL"
#define  ccsSHM_DEFAULT_SIZE 10000

/*************************************
 *  Definition of the single parts   *
 *************************************/
typedef struct  {                   
    int        semId;
    vltUINT32  stackId;    /* Complete stack id   */ 
} ccsSHMError; 

/*************************************
 *  Definition of Shared Memory Area *
 *************************************/
typedef struct  {           /* Structure of the SHM Area  */
    int          shmId;
    vltUINT32    size;
    ccsSHMError  error;
    char         test[80];  /* Used by the test procedure */    
} ccsSHMArea;  

#endif  /* !CCSSHM_H */
