/**************************************************************************
* E.S.O. - VLT project
#
# "@(#) $Id: dbDefine.h 213210 2011-04-19 08:37:01Z tebert $" 
*
* <dbDefine.h>  -  CCS/ON-LINE DATABASE Interface File
*
* who        when      what
* ---------  --------  ----------------------------------------------
* mcomin     15/09/2006 created 
*
************************************************************************
*/
#ifndef DB_DEFINE
#define DB_DEFINE

/*
 ************************************************ 
 *  Database types and constant definitions     *
 ************************************************
 */

#ifndef TRUE
#define TRUE  1
#define FALSE 0
#endif

#define dbEMPTY     ""              /* Dummy value for optional args    */
                                    /* not provided                     */
#define dbALL_ENV   "allEnvs"

/*
 * Residence flag for a Point's values
 */
#define  dbRAM     0
#define  dbDISK    1

/*
 * Direct address specification level.   The first value is not a valid direct
 * address type, but can be used to mark an DB_XREF as invalid or unused.
 */
#define  dbDIR_INVALID   0
#define  dbDIR_PLIN      1
#define  dbDIR_AIN       2
#define  dbDIR_REC       3
#define  dbDIR_FIELD     4

/*
 * types of attributes
 */
#define    dbSCALAR   0
#define    dbVECTOR   1
#define    dbTABLE    2

/*
 * element types.
 */

#define dbUNDEFINED             ((dbTYPE) 0)

#define dbLOGICAL               ((dbTYPE) 1)
#define dbINT8                  ((dbTYPE) 2)
#define dbUINT8                 ((dbTYPE) 3)
#define dbINT16                 ((dbTYPE) 4)
#define dbUINT16                ((dbTYPE) 5)
#define dbINT32                 ((dbTYPE) 6)
#define dbUINT32                ((dbTYPE) 7)
#define dbFLOAT                 ((dbTYPE) 8)
#define dbDOUBLE                ((dbTYPE) 9)

#define dbPOLAR                 ((dbTYPE) 10)
#define dbRECTANGULAR           ((dbTYPE) 11)

#define dbBYTES4                ((dbTYPE) 16)
#define dbBYTES8                ((dbTYPE) 17)
#define dbBYTES12               ((dbTYPE) 18)
#define dbBYTES16               ((dbTYPE) 19)
#define dbBYTES20               ((dbTYPE) 20)
#define dbBYTES32               ((dbTYPE) 21)
#define dbBYTES48               ((dbTYPE) 22)
#define dbBYTES64               ((dbTYPE) 23)
#define dbBYTES80               ((dbTYPE) 24)
#define dbBYTES128              ((dbTYPE) 25)
#define dbBYTES256              ((dbTYPE) 26)

#define dbXREF                  ((dbTYPE) 28)

#define dbDATE                  ((dbTYPE) 29)
#define dbTIME_OF_DAY           ((dbTYPE) 30)
#define dbABS_TIME              ((dbTYPE) 31)

#define  dbITEM_LEN         255  /* Max len for an attribute specification */
#define  dbMAX_ATTR_CNT     255  /* max number of attributes per point     */  
#define  dbMAX_FIELD_CNT    255  /* Max number of fields for a table       */
#define  dbFIELDNAME_LEN    19   /* Max len for a field name               */
#define  dbATTRIBUTE_LEN    19   /* Max length of an attribute name        */
#define  dbMAX_CHILDREN     255  /* from rtMAX_CHILDREN */

#define  dbMAX_ELEMENT_SIZE 256  /* biggest element size */

/*
 * Message constants
 */
#define dbFROM_FILE_TRUE_STRING  "True"
#define dbFROM_FILE_FALSE_STRING "False"

#define dbFILE_SOURCE_STRING     "File"
#define dbBUFFER_SOURCE_STRING   "Buffer"

#define dbMODE_BRANCH_STRING     "Branch"
#define dbMODE_POINT_STRING      "Point"

#define dbATTR_SCALAR_STRING     "Scalar"
#define dbATTR_VECTOR_STRING     "Vector"
#define dbATTR_TABLE_STRING      "Table"

#define dbADDRESS_PLIN_STRING 	 "plin"
#define dbADDRESS_AIN_STRING 	 "ain"
#define dbADDRESS_RECORD_STRING	 "record"
#define dbADDRESS_FIELD_STRING	 "field"

#define dbABSOLUTE_STRING	 "Absolute"
#define dbRELATIVE_STRING	 "Relative"
#define dbALIAS_STRING		 "Alias"

#define dbTYPE_LOGICAL_STRING	  "dbLOGICAL"
#define dbTYPE_INT8_STRING	  "dbINT8"
#define dbTYPE_UINT8_STRING	  "dbUINT8"
#define dbTYPE_INT16_STRING	  "dbINT16"
#define dbTYPE_UINT16_STRING	  "dbUINT16"
#define dbTYPE_INT32_STRING	  "dbINT32"
#define dbTYPE_UINT32_STRING	  "dbUINT32"
#define dbTYPE_FLOAT_STRING	  "dbFLOAT"
#define dbTYPE_DOUBLE_STRING	  "dbDOUBLE"
#define dbTYPE_BYTES4_STRING	  "dbBYTES4"
#define dbTYPE_BYTES8_STRING	  "dbBYTES8"
#define dbTYPE_BYTES12_STRING	  "dbBYTES12"
#define dbTYPE_BYTES16_STRING	  "dbBYTES16"
#define dbTYPE_BYTES20_STRING	  "dbBYTES20"
#define dbTYPE_BYTES32_STRING	  "dbBYTES32"
#define dbTYPE_BYTES48_STRING	  "dbBYTES48"
#define dbTYPE_BYTES64_STRING	  "dbBYTES64"
#define dbTYPE_BYTES80_STRING	  "dbBYTES80"
#define dbTYPE_BYTES128_STRING	  "dbBYTES128"
#define dbTYPE_BYTES256_STRING	  "dbBYTES256"
#define dbTYPE_XREF_STRING        "dbXREF"
#define dbTYPE_DATE_STRING        "dbDATE"
#define dbTYPE_TIME_OF_DAY_STRING "dbTIME_OF_DAY"
#define dbTYPE_ABS_TIME_STRING    "dbABS_TIME"
#define dbTYPE_POLAR_STRING       "dbPOLAR"
#define dbTYPE_RECTANGULAR_STRING "dbRECTANGULAR"

#endif /*!DB_DEFINE*/

