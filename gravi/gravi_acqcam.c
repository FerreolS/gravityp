/* $Id: gravi_vis.c,v 1.10 2014/11/12 15:10:40 nazouaoui Exp $
 *
 * This file is part of the GRAVI Pipeline
 * Copyright (C) 2002,2003 European Southern Observatory
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/**
 * @defgroup gravity_acqcam  TBD
 */
/**@{*/

/*-----------------------------------------------------------------------------
                                   Includes
 -----------------------------------------------------------------------------*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#define _XOPEN_SOURCE
#include <cpl.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "gravi_data.h"
#include "gravi_pfits.h"

#include "gravi_cpl.h"
#include "gravi_utils.h"

#include "gravi_acqcam.h"

/* ACQ CAM instrument workstation code */
#include "vltPortGeneral.h"
#include "gvacqTypes.h"
#include "gvoProcessImageAC.h"

/*-----------------------------------------------------------------------------
                               Private prototypes
 -----------------------------------------------------------------------------*/

cpl_table * gravi_acqcam_table_create (cpl_imagelist * acqcam_imglist, cpl_propertylist * header);

int AcqCamImConvert(cpl_image **DetPointer, cpl_propertylist * header);

/*-----------------------------------------------------------------------------
                             Functions code
 -----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
 * @brief Reduce the ACQ camera images
 *  
 * @param data        The gravi_data input/output
 * 
 * The routine read the images from the ACQ camera, reduce these images
 * with the XXXXX algorithm. The resulting tables is append into
 * data.
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_reduce_acqcam (gravi_data * data, gravi_data * acqcam_map, char * acqcam_map_fits)
{
    gravi_msg_function_start(1);
    cpl_ensure_code (data, CPL_ERROR_NULL_INPUT);

    /* Check if extension exist */
    if (!gravi_data_has_extension (data, GRAVI_IMAGING_DATA_ACQ_EXT)) {
        cpl_msg_warning (cpl_func, "Cannot reduce the ACQCAM, not data");
        return CPL_ERROR_NONE;
    }

    /* Get the data and header */
    cpl_imagelist * acqcam_imglist;
    acqcam_imglist = gravi_data_get_cube (data, GRAVI_IMAGING_DATA_ACQ_EXT);

    cpl_propertylist * data_header;
    data_header = gravi_data_get_header (data);
    CPLCHECK_MSG ("Cannot get data or header");

    double FIfieldWindowSizeAndSigma[2] = { 110, 3.0 };
    setFIfitWindowSizeAndSigma(FIfieldWindowSizeAndSigma);

    int coordinates[8] = { 287, 803, 1257, 1699, 228, 222, 204, 214 };
    setCoordinates(coordinates);

    double FIminFWHM     = 1.5;
    setFIminFWHM(FIminFWHM);

    int    FIwindowGauss = 16;
    setFIwindowGauss(FIwindowGauss);

    cpl_imagelist * acqcam_map_imglist;
    acqcam_map_imglist = gravi_data_get_cube (acqcam_map, GRAVI_MAP_DATA_ACQ_EXT);

    /*Deal with Flats, Sky and Dead */
    enableFromImagelist(acqcam_map_imglist);

    int    ABSLensletSizeIN       = 10;
    setABSLensletSize(ABSLensletSizeIN);

    FILE *dbcfg;

    int intABScoordinates[154 * 4];
    dbcfg = fopen("/home/aamorim/gvacq/dbl/ImageABSCoordinates.dbcfg", "r");
    int i,j;
    char *line;
    char buffer[100];
    for (i = 0; i < 2; i++) {
        line = fgets(buffer, 100, dbcfg);
    }
    int nvar;
    for (i = 0; i <  4; i++) {
       for (j = 0; j < 154; j++) {
           nvar = fscanf(dbcfg, "%d\n", intABScoordinates + i*154 + j);
       }
    }
    fclose(dbcfg);
    double ABScoordinates[154 * 4];
    for (i = 0; i < 154 * 4; i++) {
        ABScoordinates[i] = intABScoordinates[i];
    }

    int    ABSwindowSize          = 176;
    cpl_table *ABSWindowPos_table = gravi_data_get_table (acqcam_map, "ABSwindowPos");
    int* ABSwindowPos =	cpl_table_get_data_int (ABSWindowPos_table, "ABSwindowPos");

//AAMORIM
//char AB_inv_Z2S_68_136Fits[256]      = "../../Fits/gvacq_ABS_inv_Z2S_68_136.fits";
printf("XXXXXXXXXXXXXXXXXXXXXXx %s XXXXXXXXXXXXXXXXXXXXXx\n",acqcam_map_fits);

    double DummyDefocusFactory[32] = {
        16, -16, 0, 0, 0, 0, 16, -16,
        16, -16, 0, 0, 0, 0, 16, -16,
        16, -16, 0, 0, 0, 0, 16, -16,
        16, -16, 0, 0, 0, 0, 16, -16
    };
    setABSrefCoordinates(ABSwindowPos,ABScoordinates,ABSwindowSize,acqcam_map_fits,DummyDefocusFactory);

    double PT_spots_scan_sigma[4] = { 10.0, 9.0, 8.0, 7.0 };
    setPTspotsScanSigma(PT_spots_scan_sigma);

    /* Setting window coordinates to main function */
    int    CheckIteration_in = 50;
    setCheckIteration(CheckIteration_in);

    /* Create the output table */
    cpl_table * acqcam_table;
    acqcam_table = gravi_acqcam_table_create (acqcam_imglist, data_header);
	CPLCHECK_MSG ("Cannot compute acqcam_table");

    /* Add this output table in the gravi_data */
	gravi_data_add_table (data, NULL, GRAVI_OI_VIS_ACQ_EXT, acqcam_table);
	CPLCHECK_MSG ("Cannot add acqcam_table in data");
    
    gravi_msg_function_exit(1);
    return CPL_ERROR_NONE;   
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Reduce the ACQ camera images
 *  
 * @param acqcam_imglist     The imagelist with the ACQ images
 * @param header             The header corresponding to these images
 * 
 * The routine creates a cpl_table with NDIT x 4 rows, ordered
 * T1,T2,T3,T4,T1,T2... and with the following columns:
 * 
 * TIME [us]
 * FIELD_X, FIELD_Y  FIELD2_X FIELD2_Y
 * PUPIL_X, PUPIL_Y, PUPIL_Z 
 * ABERATION, array of 27 zernick coefficients [???]
 */
/*----------------------------------------------------------------------------*/

cpl_table * gravi_acqcam_table_create (cpl_imagelist * acqcam_imglist,
                                       cpl_propertylist * header)
{
    gravi_msg_function_start(1);
    
    cpl_ensure (acqcam_imglist, CPL_ERROR_NULL_INPUT, NULL);
    cpl_ensure (header,         CPL_ERROR_NULL_INPUT, NULL);

    cpl_size ntel = 4, nzernick = 27;
    cpl_size nrow = cpl_imagelist_get_size (acqcam_imglist);

    /* Here we shall build the table */
    cpl_table * output_table;
    output_table = cpl_table_new (nrow * ntel);

    /* Time column shall contain the time from PCR.ACQ.START in [us] */
    cpl_table_new_column (output_table, "TIME", CPL_TYPE_INT);
    cpl_table_set_column_unit (output_table, "TIME", "us");

    /* Field positions (or we use array of 2) */
    cpl_table_new_column (output_table, "FIELD_X", CPL_TYPE_DOUBLE);
    cpl_table_new_column (output_table, "FIELD_Y", CPL_TYPE_DOUBLE);
    cpl_table_new_column (output_table, "FIELD2_X", CPL_TYPE_DOUBLE);
    cpl_table_new_column (output_table, "FIELD2_Y", CPL_TYPE_DOUBLE);

    /* Pupil positions (or we use array of 3)  */
    cpl_table_new_column (output_table, "PUPIL_X", CPL_TYPE_DOUBLE);
    cpl_table_new_column (output_table, "PUPIL_Y", CPL_TYPE_DOUBLE);
    cpl_table_new_column (output_table, "PUPIL_Z", CPL_TYPE_DOUBLE);

    /* Aberation as arrays of 27 coefficients */
    cpl_table_new_column_array (output_table, "ABERATION", CPL_TYPE_DOUBLE, nzernick);

    int    PTwindowGaussUT        = 8;   /* Gauss fitting window size PT*/
    int    PTwindowGaussAT        = 12;  /* Gauss fitting window size PT*/
    int    PTwindowGauss;
    PTwindowGauss = PTwindowGaussAT;
    double lFIPixelScale = 17.78e-3;
    if (!strcmp(gravi_pfits_get_telescope(header),"ESO-VLTI-U1234")) {
              PTwindowGauss = PTwindowGaussUT;
              lFIPixelScale = 4.44*17.78e-3;
    }
    setPTwindowGauss(PTwindowGauss);

    /* Preparing variable for compatibility with online computation */
    int blinkOFF=1;
    cpl_image * DetPointer;
    double fluxON, fluxOFF;
    int    llx, lly, urx, ury;
    cpl_error_code error = 0;


    double PT5thSpotcoordinates[32] =
    {
        264.6,   264.1,  319.8,  320.8,
        731.1,   731.0,  787.4,  787.5,
        1201.6, 1201.1, 1257.5, 1257.9,
        1670.2, 1670.2, 1726.3, 1726.3,
        1347.7, 1403.4, 1405.5, 1348.3,
        1343.3, 1398.8, 1399.0, 1342.9,
        1346.6, 1403.8, 1403.5, 1346.7,
        1354.8, 1411.6, 1411.3, 1354.7
    };
    /* read values from fits header ESO ACQ PTFC REFPOS1...32 */
    int i,j;
    for (i=0; i< 4; i++ ) 
      for (j=0; j<4 ; j++ ) {
         PT5thSpotcoordinates[j+4*i] = gravi_pfits_get_ptfc_acqcam (header, i+4*j+1);
         PT5thSpotcoordinates[j+4*i+16] = gravi_pfits_get_ptfc_acqcam (header, i+4*j+1+16);
    }
    setPTrefCoordinates(PT5thSpotcoordinates);
    

    double lAltitude         = 20.0;                                               /* degrees */
    double lTemperature      = 11.5;                                               /* Centrigrade */
    double lPressure         = 743.0;                                              /*milli bars */
    double lHumidity         = 14.5;                                               /* percentage */
    double lParallacticAngle = 30;                                                 /* degrees */

    double         lPupilRotationAngles[4]={0,0,0,0};
    double         lImageRotationAngles[4]={0,0,0,0};
    double         lKmirrorRotation[4]={0,0,0,0};
    double         lkmirrorOffset[4]={0,0,0,0};
    double         lrotationOffset[4]={0,0,0,0};

    /* the ESO ARD CORX  and CORY are also defined in the fits header. We could pick them from there */
    setTelescopeEnvironment(lTemperature,lPressure,lHumidity,lAltitude,lParallacticAngle,lFIPixelScale,lImageRotationAngles,lPupilRotationAngles,lKmirrorRotation,lkmirrorOffset,lrotationOffset);

    double HKmagnitude_in    = 5.0;
    setHKmagnitude(HKmagnitude_in);

    double PolyACQ[11]       = { 1.37822, 0.77782, 0.730068, 0.928234, 0, 0, 0, 0, 0, 0, 0 }; /* array moved to database */
    double PolyBC[11]        = { 2.2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    setHKpolynomial(PolyACQ, PolyBC);

    /* end of defining image processing parameters */


    /* defining allDataVars structure */
    /** Member attributes used to read and write database lists. */
    gvacqAC_ABS_REF              mAcAbsRef;
    gvacqAC_ABS_REF_ERR          mAcAbsRefErr;
    gvacqAC_ABS_ZERNIKE          mAcAbsZernike;
    gvacqAC_ABS_ZERNIKE          mAcAbsZernikeErr;
    gvacqAC_PT_REF               mAcPtRef;
    gvacqAC_PT_REF_ERR           mAcPtRefErr;
    gvacqAC_PT_POS               mAcPtPos;
    gvacqAC_PT_POS_ERR           mAcPtPosErr;
    gvacqAC_PT_REF_SPOT          mAcPtRefSpot;
    gvacqAC_PT_REF_SPOT_ERR      mAcPtRefSpotErr;
    gvacqAC_PT_PT_BARY_CEN       mAcPtBaryCenter;
    gvacqAC_PT_SPOTS_FLUX        mAcPtSpotsFlux;
    gvacqAC_FI                   mAcFiPos;
    gvacqAC_FI                   mAcFiPos2;
    gvacqAC_FI_ERR               mAcFiPosErr;
    gvacqAC_FI                   mAcFiObj;
    gvacqAC_FI                   mAcFiObj2;
    gvacqAC_FI_ERR               mAcFiObjErr;
    gvacqAC_FI_TR_STS            mAcFiTrSts;
    gvacqAC_PU_TR_STS            mAcPuTrSts;
    gvacqAC_AB_TR_STS            mAcAbTrSts;
    gvacqAC_PT_CAL_STS           mAcPtCalSts;
    gvacqAC_double               mAcFiAdrAngle;
    gvacqAC_FI_ADR_CORR          mAcFiAdrCorr;
    gvacqAC_FI_FWHM              mAcFiFWHM;
    gvacqAC_FI_Flux              mAcFiFlux;
    gvacqAC_double               mAcFiAdrHK;
    gvacqAC_double               mAcFiAdrLambdaACQ;
    gvacqAC_double               mAcFiAdrLambdaBC;
    gvacqAC_int                  mAcAbsLensletSize;
    gvacqAC_ABS_LENS_SZ          mAcAbsLensletSizeArray;
    gvacqAC_Defocus              mAcAbsDefocus;
    gvacqAC_Ref_Defocus          mAcAbsRefDefocus;
    gvacqAC_TargetPos_Defocus    mAcAbsTargetPosDefocus;
    gvacqAC_ABS_TARGET_POS       mAcAbsTargetPos;
    gvacqAC_ABS_SPOTS_FLUX       mAcAbsTargetSpotsFlux;
    gvacqAC_ABS_SPOTS_SLOPES     mAcAbsSlopes;
    gvacqAC_ABS_Defocus_GaussFit mAcAbsDefocusGaussFit;

    gvacqALL_DATA allDataVars = {
        &mAcAbsRef,              &mAcAbsRefErr,      &mAcAbsZernike,
        &mAcAbsZernikeErr,       &mAcAbsLensletSize, &mAcAbsLensletSizeArray,
        &mAcPtRef,               &mAcPtRefErr,       &mAcPtPos,              &mAcPtPosErr,
        &mAcPtRefSpot,           &mAcPtRefSpotErr,
        &mAcFiPos,               &mAcFiPos2,         &mAcFiPosErr,           &mAcFiObj,         &mAcFiObj2,   &mAcFiObjErr,
        &mAcFiTrSts,             &mAcPuTrSts,        &mAcAbTrSts,            &mAcPtCalSts,
        &mAcFiAdrAngle,          &mAcFiAdrCorr,      &mAcFiFWHM,             &mAcFiFlux,        &mAcFiAdrHK,
        &mAcFiAdrLambdaACQ,      &mAcFiAdrLambdaBC,
        &mAcPtBaryCenter,        &mAcPtSpotsFlux,    &mAcAbsDefocus,         &mAcAbsRefDefocus,
        &mAcAbsTargetPosDefocus, &mAcAbsTargetPos,
        &mAcAbsTargetSpotsFlux,  &mAcAbsSlopes,      &mAcAbsDefocusGaussFit
    };

    for (i = 0; i < 4; i++) {
        allDataVars.acPtRefSpot->value[i].arm1 = PT5thSpotcoordinates[i];
        allDataVars.acPtRefSpot->value[i].arm2 = PT5thSpotcoordinates[i + 4];
        allDataVars.acPtRefSpot->value[i].arm3 = PT5thSpotcoordinates[i + 8];
        allDataVars.acPtRefSpot->value[i].arm4 = PT5thSpotcoordinates[i + 12];

        allDataVars.acPtRefSpot->value[i + 4].arm1 = PT5thSpotcoordinates[i + 16];
        allDataVars.acPtRefSpot->value[i + 4].arm2 = PT5thSpotcoordinates[i + 4 + 16];
        allDataVars.acPtRefSpot->value[i + 4].arm3 = PT5thSpotcoordinates[i + 8 + 16];
        allDataVars.acPtRefSpot->value[i + 4].arm4 = PT5thSpotcoordinates[i + 12 + 16];
    }

    int    RefWindowPosPT[8]      = { 292, 759, 1230, 1698, 1376, 1371, 1375, 1383 };
    setRefWindowPosPT(RefWindowPosPT);

    int    PTImageWindowSize      = 220; /* pupil tracker image window */
    setPTImageWindowSize(PTImageWindowSize);

    llx=  RefWindowPosPT[0] - PTImageWindowSize / 2;
    urx = RefWindowPosPT[3] + PTImageWindowSize / 2;
    lly = RefWindowPosPT[4] - PTImageWindowSize / 2;
    ury = RefWindowPosPT[7] + PTImageWindowSize / 2;

    int converted=0;
    /* tell if the first image is blink on or off */
    if (nrow>1) {
          DetPointer = cpl_imagelist_get(acqcam_imglist,0);
          converted=AcqCamImConvert(&DetPointer,header);
          fluxON     = cpl_image_get_flux_window(DetPointer, llx, lly, urx, ury);
          if (converted) cpl_image_delete(DetPointer);
          DetPointer=cpl_imagelist_get(acqcam_imglist,1);
          converted=AcqCamImConvert(&DetPointer,header);
          fluxOFF    = cpl_image_get_flux_window(DetPointer, llx, lly, urx, ury);
          if (converted) cpl_image_delete(DetPointer);
          if (fluxON > fluxOFF) blinkOFF = 0;
    }

    if (nrow == 0)
        gravi_msg_warning ("INFO", "No ACQ_CAM images in the observation");

    /* Loop on DIT in cube */
    for (cpl_size row = 0; row < nrow; row++) {

        /* Fill the TIME column (same value for all beams) */
        double time = gravi_pfits_get_time_acqcam (header, row);

        for (int tel = 0; tel < ntel; tel ++)
            cpl_table_set (output_table, "TIME", row*ntel+tel, time);
        DetPointer=cpl_imagelist_get(acqcam_imglist,row);
        converted=AcqCamImConvert(&DetPointer,header);
        error = gvoacqProcessImageAC(DetPointer, blinkOFF, ACMODE_FIELD, &allDataVars);
        error = gvoacqProcessImageAC(DetPointer, blinkOFF, ACMODE_PUPIL, &allDataVars);
        error = gvoacqProcessImageAC(DetPointer, blinkOFF, ACMODE_ABERR, &allDataVars);
        storeSky(DetPointer, STORSKY_PUPIL);
    
        int tel=0;
        cpl_table_set (output_table, "FIELD_X", row*ntel+tel, (allDataVars.acFiPos)->value[0].arm1);
        cpl_table_set (output_table, "FIELD_Y", row*ntel+tel, (allDataVars.acFiPos)->value[1].arm1);
        cpl_table_set (output_table, "FIELD2_X", row*ntel+tel, (allDataVars.acFiPos2)->value[0].arm1);
        cpl_table_set (output_table, "FIELD2_Y", row*ntel+tel, (allDataVars.acFiPos2)->value[1].arm1);
        tel=1;
        cpl_table_set (output_table, "FIELD_X", row*ntel+tel, (allDataVars.acFiPos)->value[0].arm2);
        cpl_table_set (output_table, "FIELD_Y", row*ntel+tel, (allDataVars.acFiPos)->value[1].arm2);
        cpl_table_set (output_table, "FIELD2_X", row*ntel+tel, (allDataVars.acFiPos2)->value[0].arm2);
        cpl_table_set (output_table, "FIELD2_Y", row*ntel+tel, (allDataVars.acFiPos2)->value[1].arm2);
        tel=2;
        cpl_table_set (output_table, "FIELD_X", row*ntel+tel, (allDataVars.acFiPos)->value[0].arm3);
        cpl_table_set (output_table, "FIELD_Y", row*ntel+tel, (allDataVars.acFiPos)->value[1].arm3);
        cpl_table_set (output_table, "FIELD2_X", row*ntel+tel, (allDataVars.acFiPos2)->value[0].arm3);
        cpl_table_set (output_table, "FIELD2_Y", row*ntel+tel, (allDataVars.acFiPos2)->value[1].arm3);
        tel=3;
        cpl_table_set (output_table, "FIELD_X", row*ntel+tel, (allDataVars.acFiPos)->value[0].arm4);
        cpl_table_set (output_table, "FIELD_Y", row*ntel+tel, (allDataVars.acFiPos)->value[1].arm4);
        cpl_table_set (output_table, "FIELD2_X", row*ntel+tel, (allDataVars.acFiPos2)->value[0].arm4);
        cpl_table_set (output_table, "FIELD2_Y", row*ntel+tel, (allDataVars.acFiPos2)->value[1].arm4);
    
        tel=0;
        cpl_table_set (output_table, "PUPIL_X", row*ntel+tel, (allDataVars.acPtPos)->value[0].arm1);
        cpl_table_set (output_table, "PUPIL_Y", row*ntel+tel, (allDataVars.acPtPos)->value[1].arm1);
        cpl_table_set (output_table, "PUPIL_Z", row*ntel+tel, (allDataVars.acPtPos)->value[2].arm1);
        tel=1;
        cpl_table_set (output_table, "PUPIL_X", row*ntel+tel, (allDataVars.acPtPos)->value[0].arm2);
        cpl_table_set (output_table, "PUPIL_Y", row*ntel+tel, (allDataVars.acPtPos)->value[1].arm2);
        cpl_table_set (output_table, "PUPIL_Z", row*ntel+tel, (allDataVars.acPtPos)->value[2].arm2);
        tel=2;
        cpl_table_set (output_table, "PUPIL_X", row*ntel+tel, (allDataVars.acPtPos)->value[0].arm3);
        cpl_table_set (output_table, "PUPIL_Y", row*ntel+tel, (allDataVars.acPtPos)->value[1].arm3);
        cpl_table_set (output_table, "PUPIL_Z", row*ntel+tel, (allDataVars.acPtPos)->value[2].arm3);
        tel=3;
        cpl_table_set (output_table, "PUPIL_X", row*ntel+tel, (allDataVars.acPtPos)->value[0].arm4);
        cpl_table_set (output_table, "PUPIL_Y", row*ntel+tel, (allDataVars.acPtPos)->value[1].arm4);
        cpl_table_set (output_table, "PUPIL_Z", row*ntel+tel, (allDataVars.acPtPos)->value[2].arm4);
    
        
        if (converted) cpl_image_delete(DetPointer);
        if (blinkOFF) blinkOFF=0;
        else blinkOFF=1;

    } /* End loop on DIT in cube */

    return output_table;
}

int AcqCamImConvert(cpl_image **DetPointer, cpl_propertylist * header){
    int i,j,nx,ny,frame,framex,framey,strx,stry,converted;
    cpl_type type;
    cpl_image *DetPointerNew;
    int sizex,sizey,x0,y0i,pis_rejected;
    char name[256];
    float content;

    nx=cpl_image_get_size_x(*DetPointer);
    ny=cpl_image_get_size_y(*DetPointer);
    converted=0;

    /* Check if conversion is needed */
    if ( nx != 2048 || ny != 1536 ){
       /* if yes create a new image with the full size */
       converted=1;
       type=cpl_image_get_type(*DetPointer);
       DetPointerNew=cpl_image_new(2048,1536,type); 
       sizex = cpl_propertylist_get_int(header, "ESO DET1 FRAMES NX");
       sizey = cpl_propertylist_get_int(header, "ESO DET1 FRAMES NY");

       for (frame=0; frame<16; frame++) {
           framex=frame-(frame/4)*4;
           framey=frame/4;
           sprintf(name,"ESO DET1 FRAM%d STRX",frame+1);
           strx = cpl_propertylist_get_int(header, name);
           sprintf(name,"ESO DET1 FRAM%d STRY",frame+1);
           stry = cpl_propertylist_get_int(header, name);
           for (i=1; i<=sizex; i++)
               for (j=1; j<=sizey; j++) {
                   content=cpl_image_get(*DetPointer,i+framex*sizex,j+framey*sizey,&pis_rejected);
                   cpl_image_set(DetPointerNew,i+strx,j+stry,content);
               } 
       }      

       *DetPointer=DetPointerNew;
    } 
    return converted;    
}

/**@}*/
