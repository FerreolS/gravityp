/************************************************************************
 *   This file is part of the E.S.O. - VLTI GRAVITY project
 *
 *  "@(#) $Id: gvoacqPreprocess.c 274087 2015-09-26 22:08:53Z tott $
 *   
 *   who       when        what
 *   --------  ----------  ----------------------------------------------
 *   narsi    2013-05-27  created
 *   narsi    2013-06-11  Added preprocessing stubs
 *   ekw      2014-02-11  Adopt to shared memory and Fits handling
 */

/****************************************************************************
 *   NAME 
 *   gvoacqPreprocess - procedere to compute pre processing of the detector image
 *  
 *   SYNOPSIS
 *   #include <gvoProcessImageAC.h>
 * 
 *   DESCRIPTION
 *   Given a pointer to the acquisition camera image, will preprocess the image like dead pixel correction,
 *   sky correction and flat fielding.
 *
 *   FILES
 *
 *   ENVIRONMENT
 * 
 *   RETURN VALUES
 *
 *   CAUTIONS (optional)
 *
 *   EXAMPLES (optional)
 *
 *   SEE ALSO (optional)
 *
 *   BUGS (optional)
 *
 *------------------------------------------------------------------------
 */


/*-----------------------------------------------------------------------------
                               Includes
 ----------------------------------------------------------------------------*/

#include <gvoProcessImageAC.h> /* */



/* ----------------------------------------------------------------------*/
/**
 * @internal
 * @brief  Loads the master flat file
 * @param  MasterFlatFits filename to the master flat image.
 * @param  MasterFlat (output) master flat data
 * @return error
 *  --TOBE DONE--
 * Define the procederes to interact with database 
 * @see    gvoacqLoadDeadPixel, gvoacqLoadSky
 */

cpl_error_code gvoacqLoadFlat(char * MasterFlatFits, cpl_image * MasterFlat)
{
    cpl_error_code   error =CPL_ERROR_NONE;
    
    MasterFlat = cpl_image_load(MasterFlatFits, CPL_TYPE_DOUBLE, 0, 0);
    cpl_image_delete(MasterFlat);
    return error;
}
/* End of gvoacqLoadFlat */


/* ----------------------------------------------------------------------*/
/**
 * @internal
 * @brief  Loads the sky  files and return sky buffer
 * @param  MasterSkyFits  filename to the sky frames
 * @param  MasterSky (output) master sky data
 * @return error
 *  --TOBE DONE--
 * 1. Define the procederes to interact with database 
 * 2. Read a signle sky fits image and return sky buffer
 * @see    gvoacqLoadDeadPixel, gvoacqLoadFlat
 */

cpl_error_code gvoacqLoadSky(char * MasterSkyFits,  cpl_image * MasterSkytmp)
{
    cpl_error_code   error =CPL_ERROR_NONE;
    
    MasterSkytmp = cpl_image_load(MasterSkyFits, CPL_TYPE_DOUBLE, 0, 0);
    cpl_image_delete(MasterSkytmp);
   
    return error;
}
/* End of gvoacqLoadSky */


/* ----------------------------------------------------------------------*/
/**
 * @internal
 * @brief  Loads the sky  files and makes the master sky file and updates the database.
 * @param  MasterSkyFits  filename to the sky frames
 * @param  MasterSky (output) master sky data
 * @return error
 *  --TOBE DONE--
 * 1. Define the procederes to interact with database 
 * 2. Write the master sky file to the database
 * @see    gvoacqLoadDeadPixel, gvoacqLoadFlat
 */

cpl_error_code gvoacqMakeSky(char * MasterSkyFits, 
			     int numberOfSkyFrames,  
			     cpl_image * MasterSky)
{
    cpl_error_code   error =CPL_ERROR_NONE;
    int i;

    cpl_image * MasterSkyTemp = cpl_image_new(2048,2048,CPL_TYPE_DOUBLE);

    MasterSky = cpl_image_new(2048,2048,CPL_TYPE_DOUBLE);

    for(i=0; i<numberOfSkyFrames; i++)
        {
        error = gvoacqLoadSky(MasterSkyFits, MasterSkyTemp);
        error = cpl_image_add(MasterSky, MasterSkyTemp);
        }
    error = cpl_image_multiply_scalar(MasterSky, 1/numberOfSkyFrames);
 
/* TODO
       write the MasterSky buffer to the database
       For now it can save to the hard disk
*/
    /* cpl_image_save(MasterSky, "MASTER_Sky_DB.fits", CPL_BPP_IEEE_DOUBLE, NULL, CPL_IO_CREATE); */
    
/* FE Test 20150925 */
/*cpl_image_save(MasterSky, "Sky.fits", CPL_BPP_IEEE_DOUBLE, NULL, CPL_IO_CREATE);*/

    cpl_image_delete(MasterSky);
    cpl_image_delete(MasterSkyTemp);
    
    return error;
}
/* End of gvoacqLoadSky */




/* ----------------------------------------------------------------------*/
/**
 * @internal
 * @brief  Loads the dead pixel map file and updates the database
 * @param  DeadPixelMapFits makes the master sky file
 * @param  DeadPixelMap (output) dead pixel map  data
 * @return error
 * @see    gvacqLoadSky, gvacqLoadFlat
 *
 */

cpl_error_code gvoacqLoadDeadPixel(char * DeadPixelMapFits, cpl_mask * DeadPixelMap)
{
    cpl_error_code   error = CPL_ERROR_NONE;
   
    cpl_image * DeadPixelMapImage = cpl_image_load(DeadPixelMapFits, CPL_TYPE_DOUBLE, 0, 0);
    
    if(DeadPixelMapImage == NULL)
	{
	cpl_image_delete(DeadPixelMapImage);  
	return -1;
	}
    
    int Msizey =cpl_image_get_size_y (DeadPixelMapImage);
    int Msizex=	cpl_image_get_size_x (DeadPixelMapImage);
    
    /* 
    double Max, Min;
    Max=cpl_image_get_max(DeadPixelMapImage); 
    Min=cpl_image_get_mean(DeadPixelMapImage);
    DeadPixelMap=cpl_mask_threshold_image_create(DeadPixelMapImage, Min, Max);
    */
    
    
    if(Msizex == 2048 && Msizey == 1536){
    error = cpl_mask_threshold_image(DeadPixelMap, DeadPixelMapImage, 0, 2, CPL_BINARY_1);
    }else printf("Reading bad DeadPixelMap.fits at %d %d %s:%s:%d \n", Msizey, Msizex, __FILE__,__FUNCTION__, __LINE__);
    
    cpl_image_delete(DeadPixelMapImage);    
    return error;
}
/* End of gvacqDeadPixel */




/* ----------------------------------------------------------------------*/
/**
 * @internal
 * @brief  Does the preprocessing operations like flat-fielding etc.. on raw detector image 
 * @param  DetPointer (input and ouput) pointer to the  detector image.
 * @return error
 * @see    gvoacqFieldImager, gvoacqPupilTracker, gvoacqAberrationSensor
 *
 * @par Error Handling:
 * CPL_ERROR_NULL_INPUT if(one of ) the input pointer (s) is NULL
 * CPL_ERROR_NULL
 */

cpl_error_code gvoacqPreprocess(cpl_image *DetPointer, 
				cpl_image *MasterFlat,
				cpl_image *MasterSky,
				cpl_mask  *DeadPixelMap)
{
    /* test entries */

/* FE Test 20150925 */
/*cpl_image_save(MasterSky, "MASTER_Sky_DB.fits", CPL_BPP_IEEE_DOUBLE, NULL, CPL_IO_CREATE);*/
    
    if (MasterFlat == NULL || MasterSky == NULL || DeadPixelMap == NULL ) {
    printf("ERROR:: -> Please make sure the database reading Flat or Sky or DeadPixelMap correcly. Error occured at %s:%s:%d\n", __FILE__, __FUNCTION__, __LINE__);
    
    return -1;
    }
     
    cpl_error_code   error = 0;

    int Isizey =cpl_image_get_size_y (DetPointer);
    int Isizex=	cpl_image_get_size_x (DetPointer);
    int Fsizey =cpl_image_get_size_y (MasterFlat);
    int Fsizex=	cpl_image_get_size_x (MasterFlat);
    int Ssizey =cpl_image_get_size_y (MasterSky );
    int Ssizex=	cpl_image_get_size_x (MasterSky );
/*    int Dsizey =cpl_mask_get_size_y (DeadPixelMap); */
/*    int Dsizex=	cpl_mask_get_size_x (DeadPixelMap); */

    if (Ssizey != Isizey || Ssizex != Isizex) {
      printf("Sky image size doesn't fit detector readout\n");
      error = 1;
    }
    if (Fsizey != Isizey || Fsizex != Isizex) {
      printf("Flat image size doesn't fit detector readout\n");
      error = 1;
    }
    /* if (Dsizey != Isizey || Dsizex != Isizex) { */
    /*   printf("Dead pixel image size doesn't fit detector readout\n"); */
    /*   error = 1; */
    /* } */

    if (error != 0) {
      return error;
    } else {
      error = CPL_ERROR_NONE;
      error = cpl_image_subtract(DetPointer, MasterSky);
      if (error != CPL_ERROR_NONE) {
	return error;
      }
/*
      error =cpl_image_divide(DetPointer, MasterFlat);
      if (error != CPL_ERROR_NONE) {
	return error;
      }
*/
      /* cpl_image* Im = cpl_image_new_from_mask(DeadPixelMap); */
      /* if(cpl_image_get_max(Im) > 0 ){ */
      /* 	error = cpl_image_reject_from_mask (DetPointer, DeadPixelMap); */
      /* 	error = cpl_detector_interpolate_rejected(DetPointer); */
      /* } */
      /* cpl_image_delete(Im); */
    }

/* FE Test 20150925 */
/*cpl_image_save(DetPointer, "OnOff.fits", CPL_BPP_IEEE_DOUBLE, NULL, CPL_IO_CREATE);*/
   
    return error;
}
/**@}*/







