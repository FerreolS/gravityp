#ifndef gvacqTYPES_H
#define gvacqTYPES_H
/*******************************************************************************
* E.S.O. - VLT project
*
* "@(#) $Id: gvacqTypes.h 285696 2016-07-13 13:48:43Z ewiezorr $"
*
*/
/** \internal
* who       when        what
* --------  ----------  ----------------------------------------------
* narsi     2016-03-25  added acAbsSlopes, acAbsDefocusGaussFit and acAbsTargetPosDefocus
* narsi     2016-03-23  added acAbsTargetPos, acAbsTargetSpotsFlux and acPtSpotsFlux
* narsi     2016-03-19  added acAbsDefocus and acAbsRefDefocus
* amorim    2015-11-12  added acFiFWHM and acFiFlux
* jott      2013-07-01  created.
*/

/* #ifndef __cplusplus
#error This is a C++ include file and cannot be used from plain C
#endif
*/

#include <vltPortGeneral.h>
#include "gvacqDefines.h"


/**
 * Structure which can be used to read a single record
 * from the database tables used for the gvoacqProcessImageAC 
 * library function
 */
typedef struct
    {
    vltDOUBLE    arm1;
    vltDOUBLE    arm2;
    vltDOUBLE    arm3;
    vltDOUBLE    arm4;
    } gvacqAC_IMAGE_DATA_REC;


/**
 * Structures for reading the gvoacqProcessImageAC tables in one go.
 */  
typedef struct
    {
    gvacqAC_IMAGE_DATA_REC  value[gvacqNUM_REC_ABS_REF];
    } gvacqAC_ABS_REF;

typedef struct
    {
    gvacqAC_IMAGE_DATA_REC  value[gvacqNUM_REC_ABS_REF_ERROR];
    } gvacqAC_ABS_REF_ERR;

typedef struct
    {
    gvacqAC_IMAGE_DATA_REC  value[gvacqNUM_REC_ABS_REF];
    } gvacqAC_ABS_TARGET_POS;

typedef struct
    {
    gvacqAC_IMAGE_DATA_REC  value[gvacqNUM_REC_ABS_SPOT_NUM];
    } gvacqAC_ABS_SPOTS_FLUX;

typedef struct
    {
    gvacqAC_IMAGE_DATA_REC  value[gvacqNUM_REC_ABS_SLOPES];
    } gvacqAC_ABS_SPOTS_SLOPES;

typedef struct
    {
    gvacqAC_IMAGE_DATA_REC  value[gvacqNUM_REC_DEFOCUS];
    } gvacqAC_ABS_Defocus_GaussFit;

typedef struct
    {
    gvacqAC_IMAGE_DATA_REC  value[gvacqNUM_REC_ZERNIKE];
    } gvacqAC_ABS_ZERNIKE;

typedef struct
    {
    gvacqAC_IMAGE_DATA_REC  value[gvacqNUM_REC_ZERNIKE_ERROR];
    } gvacqAC_ABS_ZERNIKE_ERR;


typedef struct
    {
    gvacqAC_IMAGE_DATA_REC  value[gvacqNUM_REC_PT_REF_POSITION];
    } gvacqAC_PT_REF;

typedef struct
    {
    gvacqAC_IMAGE_DATA_REC  value[gvacqNUM_REC_PT_SPOT_NUM];
    } gvacqAC_PT_SPOTS_FLUX;


typedef struct
    {
    gvacqAC_IMAGE_DATA_REC  value[gvacqNUM_REC_PT_REF_POSITION_ERROR];
    } gvacqAC_PT_REF_ERR;


typedef struct
    {
    gvacqAC_IMAGE_DATA_REC  value[gvacqNUM_REC_PT_POSITION];
    } gvacqAC_PT_POS;


typedef struct
    {
    gvacqAC_IMAGE_DATA_REC  value[gvacqNUM_REC_PT_POSITION_ERROR];
    } gvacqAC_PT_POS_ERR;


typedef struct
    {
    gvacqAC_IMAGE_DATA_REC  value[gvacqNUM_REC_PT_REF_SPOT];
    } gvacqAC_PT_REF_SPOT;

typedef struct
    {
    gvacqAC_IMAGE_DATA_REC  value[gvacqNUM_REC_PT_REF_SPOT_ERROR];
    } gvacqAC_PT_REF_SPOT_ERR;

typedef struct
    {
	gvacqAC_IMAGE_DATA_REC  value[gvacqNUM_REC_PT_BARY_CENTER];
    } gvacqAC_PT_PT_BARY_CEN;

typedef struct
    {
	gvacqAC_IMAGE_DATA_REC  value[gvacqNUM_REC_REF_DEFOCUS];
    } gvacqAC_Ref_Defocus;

typedef struct
    {
	gvacqAC_IMAGE_DATA_REC  value[gvacqNUM_REC_REF_DEFOCUS];
    } gvacqAC_TargetPos_Defocus;

typedef struct
    {
    gvacqAC_IMAGE_DATA_REC  value[gvacqNUM_REC_FI];
    } gvacqAC_FI;

typedef struct
    {
    gvacqAC_IMAGE_DATA_REC  value[gvacqNUM_REC_FI_ERROR];
    } gvacqAC_FI_ERR;

typedef vltINT32 gvacqAC_FI_TR_STS[gvacqNUM_REC_TELESCOPES];

typedef vltINT32 gvacqAC_PU_TR_STS[gvacqNUM_REC_TELESCOPES];

typedef vltINT32 gvacqAC_AB_TR_STS[gvacqNUM_REC_TELESCOPES];

typedef vltINT32 gvacqAC_PT_CAL_STS[gvacqNUM_REC_TELESCOPES];

typedef vltDOUBLE gvacqAC_FI_ADR_CORR[gvacqNUM_REC_XYPOS];

typedef vltDOUBLE gvacqAC_FI_FWHM[gvacqNUM_REC_TELESCOPES];

typedef vltDOUBLE gvacqAC_FI_Flux[gvacqNUM_REC_TELESCOPES];

typedef vltDOUBLE gvacqAC_Defocus[gvacqNUM_REC_DEFOCUS];

typedef vltINT32 gvacqAC_Ftopt_TR_STS[gvacqNUM_REC_TELESCOPES];

typedef vltDOUBLE gvacqAC_double;

typedef vltINT32 gvacqAC_int;
typedef vltINT32 gvacqAC_ABS_LENS_SZ[gvacqNUM_REC_TELESCOPES];

typedef struct {
	gvacqAC_ABS_REF			  *acAbsRef;
	gvacqAC_ABS_REF_ERR		  *acAbsRefErr;
	gvacqAC_ABS_ZERNIKE		  *acAbsZernike;
	gvacqAC_ABS_ZERNIKE		  *acAbsZernikeErr;
        gvacqAC_int				  *acAbsLensletSize;
        gvacqAC_ABS_LENS_SZ		  *acAbsLensletSizeArray;
	gvacqAC_PT_REF			  *acPtRef;
	gvacqAC_PT_REF_ERR		  *acPtRefErr;
	gvacqAC_PT_POS			  *acPtPos;
	gvacqAC_PT_POS_ERR		  *acPtPosErr;
	gvacqAC_PT_REF_SPOT		  *acPtRefSpot;
	gvacqAC_PT_REF_SPOT_ERR	  *acPtRefSpotErr;
	gvacqAC_FI			      *acFiPos;
	gvacqAC_FI			      *acFiPos2;
	gvacqAC_FI_ERR			  *acFiPosErr;
	gvacqAC_FI			      *acFiObj;
	gvacqAC_FI			      *acFiObj2;
	gvacqAC_FI_ERR			  *acFiObjErr;
	gvacqAC_FI_TR_STS         *acFiTrSts;
	gvacqAC_PU_TR_STS         *acPuTrSts;
	gvacqAC_AB_TR_STS         *acAbTrSts;
	gvacqAC_PT_CAL_STS        *acPtCalSts;
	gvacqAC_double            *acFiAdrAngle;
	gvacqAC_FI_ADR_CORR       *acFiAdrCorr;
	gvacqAC_FI_FWHM           *acFiFWHM;
	gvacqAC_FI_Flux           *acFiFlux;
	gvacqAC_double            *acFiAdrHK;
	gvacqAC_double            *acFiAdrLambdaACQ;
	gvacqAC_double            *acFiAdrLambdaBC;
	gvacqAC_PT_PT_BARY_CEN    *acPtBaryCenter;
    gvacqAC_PT_SPOTS_FLUX     *acPtSpotsFlux;
        gvacqAC_Defocus           *acAbsDefocus;
        gvacqAC_Ref_Defocus       *acAbsRefDefocus;
    gvacqAC_TargetPos_Defocus   *acAbsTargetPosDefocus;
    gvacqAC_ABS_TARGET_POS    *acAbsTargetPos;
    gvacqAC_ABS_SPOTS_FLUX    *acAbsTargetSpotsFlux;
    gvacqAC_ABS_SPOTS_SLOPES  *acAbsSlopes;
    gvacqAC_ABS_Defocus_GaussFit *acAbsDefocusGaussFit;
} gvacqALL_DATA;

#endif  /* !gvacqTYPES_H */





