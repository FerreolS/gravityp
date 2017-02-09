/************************************************************************
 *   This file is part of the E.S.O. - VLTI GRAVITY project
 *
 *  "@(#) $Id: gvoacqWrite.c 263764 2015-01-14 15:58:37Z ewieprec $
 *   
 *   who       when        what
 *   --------  ----------  ----------------------------------------------
 *   narsi    2012-05-15  created

 */

/****************************************************************************
 *   NAME 
 *   gvoacqWrite - procedere to write the data into ASCII file and to instrument database
 *  
 *   SYNOPSIS
 *   #include "gvoProcessImageAC.h"
 * 
 *   DESCRIPTION
 *   Write the data to instrument database with relevant quantities.
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


/**
 * @defgroup process gravity observation software for acquisition camera
 * @ingroup image_analysis
 *
 * This module provides functions for acquision camera data reduction
 *
 * @par Synopsis:
 * @code
 #include "gvoProcessImageAC.h"
 * @endcode
 */



/*-----------------------------------------------------------------------------
                               Includes
 ----------------------------------------------------------------------------*/

#include <gvoProcessImageAC.h> /* */

/* ----------------------------------------------------------------------*/
/**
 * @internal
 * @brief  Write the data to a ASCII file
 * @param  acqStruct structure to the data. 
 * @return error
 */

cpl_error_code gvacqWrite2File(struct acqStruct *acq, char * option)
{
    
    cpl_error_code error=CPL_ERROR_NONE;    

    char Query = 'N';
    int ZernikeNum=28;
    int i;
 		int ABSSpotsNum=77;
		int PT_centroid_num=32;
    int PT_refernce_spot_num=8;

    
    /* Processed results writing to text file for now */
    FILE       * ACQdata;  /* pointer to a text file */
    ACQdata = fopen( "gvoacqData.txt", "a+" );
    

    time_t now = time(0);
    char Tbuffer[100];
    strftime(Tbuffer, 100, "# Date: %Y:%m:%d, Time: %H:%M:%S.00", localtime(&now));
    /* fprintf(ACQdata, "%10s \n", Tbuffer); */
   

    if(strcmp(option, "start") == 0) Query= 'A'; 
    else if(strcmp(option, "PTCalib")==0) Query = 'C';
    else if(strcmp(option, "ABSCalib")==0) Query = 'H';
    else if(strcmp(option, "FI")==0)Query = 'I';
    else if(strcmp(option, "PT")==0) Query = 'P';
    else if(strcmp(option, "ABS")==0) Query = 'Z';
    else if(strcmp(option, "stop")==0) Query = 'N';
    else if(strcmp(option, "DeadPixel")==0) Query = 'D'; 
    else if(strcmp(option, "Sky")==0) Query = 'S'; 
    else if(strcmp(option, "Flat")==0) Query = 'F';
    
    
    switch(Query)
	{
	case 'A':
	    fprintf(ACQdata, "\n"
		    "# ****************************  begin pupil tracker  data  ****************************\n"
		    "#The pupil tracker reference spots in pixel (x, y centriods in total 34 rows and  errors)\n"
		    "#S.No     UT1 (1 & 2 columns)       UT2 (3 & 4 columns)          UT3 (5 & 6 columns)           UT4(7 & 8 columns)\n");
	    
	    for(i=0; i<PT_centroid_num; i++)
		{
		fprintf(ACQdata,"%02d %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e  %14.7e \n", 
			i+1, 
			(acq->PT_REF)[i] , 
			(acq->PT_REF_ERROR)[i],
			(acq->PT_REF)[i+PT_centroid_num],
			(acq->PT_REF_ERROR)[i+PT_centroid_num],
			(acq->PT_REF)[i+3*PT_centroid_num],
			(acq->PT_REF_ERROR)[i+2*PT_centroid_num],
			(acq->PT_REF)[i+3*PT_centroid_num],
			(acq->PT_REF_ERROR)[i+3*PT_centroid_num]);
		}
	    break;
	    
	    
	    fprintf(ACQdata," \n"
		    "#Lateral  and longitudinal pupil positions x, y & z are 1, 2 & 3 rows respectively  (defocus example; Units: meter)\n"
		    "#S.No      UT1           UT2          UT3              UT4\n");
	    
	    for(i=0; i<3; i++){
	    fprintf(ACQdata, "%02d %+14.7e %+14.7e %+14.7e %+14.7e \n", i+1,
		    (acq->PT_POSITION)[i], 
		    (acq->PT_POSITION)[i+3], 
		    (acq->PT_POSITION)[i+6], 
		    (acq->PT_POSITION)[i+9] );
	    }

	    fprintf(ACQdata," \n"
		    "#Lateral  and longitudinal pupil position errors x, y & z are 1, 2 & 3 rows respectively  (defocus example; Units: meter)\n"
		    "#S.No      UT1           UT2          UT3              UT4\n");
	    
	    for(i=0; i<3; i++){
	    fprintf(ACQdata, "%02d %+14.7e %+14.7e %+14.7e %+14.7e \n", i+1,
		    (acq->PT_POSITION_ERROR)[i], 
		    (acq->PT_POSITION_ERROR)[i+3], 
		    (acq->PT_POSITION_ERROR)[i+6], 
		    (acq->PT_POSITION_ERROR)[i+9] );
	    }
	    
	    fprintf(ACQdata, "\n"
		    "#Reference spot position (5th spot) for each UT (x & y centriods) (in pixels)\n"
		    "#S.No          UT1           UT2            UT3           UT4   \n");
	    
	    for(i=0; i<PT_refernce_spot_num; i++)
		{
		fprintf(ACQdata,"%02d  %+14.7e %+14.7e %+14.7e %+14.7e \n", 
			i+1, 
			(acq->PT_REF_SPOT)[i],
			(acq->PT_REF_SPOT)[i+PT_refernce_spot_num],
			(acq->PT_REF_SPOT)[i+2*PT_refernce_spot_num],
			(acq->PT_REF_SPOT)[i+3*PT_refernce_spot_num]);
		}
	    fprintf(ACQdata, "\n"
		    "#Reference spot position error (5th spot) for each UT (x & y centriods) (in pixels)\n"
		    "#S.No          UT1           UT2            UT3           UT4   \n");
	    
	    for(i=0; i<PT_refernce_spot_num; i++)
		{
		fprintf(ACQdata,"%02d  %+14.7e %+14.7e %+14.7e %+14.7e \n", 
			i+1, 
			(acq->PT_REF_SPOT_ERROR)[i],
			(acq->PT_REF_SPOT_ERROR)[i+PT_refernce_spot_num],
			(acq->PT_REF_SPOT_ERROR)[i+2*PT_refernce_spot_num],
			(acq->PT_REF_SPOT_ERROR)[i+3*PT_refernce_spot_num]);
		}
	    
	    fprintf(ACQdata, "\n" 
		    "#********************************* begin aberration sensor data ***********************************\n"
		    "#The aberration sensor reference case centriod data for  77 spots (154 slopes for each UT; pixels)  \n");
	    
	    fprintf(ACQdata,"#S.No          UT1           UT2            UT3           UT4   \n");
	    for(i=0; i<ABSSpotsNum*2; i++)
		{
		fprintf(ACQdata,"%03d %+14.7e %+14.7e %+14.7e %+14.7e \n", 
			i+1, 
			(acq->ABS_REF)[i] , 
			(acq->ABS_REF)[i+ABSSpotsNum*2],
			(acq->ABS_REF)[i+2*ABSSpotsNum*2],
			(acq->ABS_REF)[i+3*ABSSpotsNum*2]);
	    }
	    
	    
	    fprintf(ACQdata, "\n"
		    "#The aberration sensor fitted Zernike Coefficients (77 measured here for each UT; defocus aberration case here; Units: lambda)  \n"
		    "#S.No          UT1           UT2            UT3           UT4   \n");
	    for(i=0; i<ZernikeNum; i++)
		{
		fprintf(ACQdata,"%02d %+14.7e %+14.7e %+14.7e %+14.7e \n", 
			i+1, 
			(acq->ABS_ZERNIKE)[i] , 
			(acq->ABS_ZERNIKE)[i+ZernikeNum],
			(acq->ABS_ZERNIKE)[i+2*ZernikeNum],
			(acq->ABS_ZERNIKE)[i+3*ZernikeNum]);
		}
	    
	    
	    fprintf(ACQdata, "\n"
		    "#The aberration sensor fitted Zernike Coefficient error (68 measured here for each UT; defocus aberration case here; Units: lambda)  \n"
		    "#S.No          UT1           UT2            UT3           UT4   \n");
	    for(i=0; i<ZernikeNum; i++)
		{
		fprintf(ACQdata,"%02d %+14.7e %+14.7e %+14.7e %+14.7e \n", 
			i+1, 
			(acq->ABS_ZERNIKE_ERROR)[i] , 
			(acq->ABS_ZERNIKE_ERROR)[i+ZernikeNum],
			(acq->ABS_ZERNIKE_ERROR)[i+2*ZernikeNum],
			(acq->ABS_ZERNIKE_ERROR)[i+3*ZernikeNum]);
		}
	    
	    
	    fprintf(ACQdata, "\n" "*********************** begin field imager data ********************************\n"
		    "#The field imager: brightest object position and its error for each UT (in pixels) \n"
		    "T.No Postion(X)   Position error(X)   Position(Y)    Position error(Y)\n");
	    
	    for(i=0; i<4;i++)
		{
		fprintf(ACQdata,"%d %14.7e %14.7e %14.7e  %14.7e\n", i+1, 
			(acq->FI)[i], 
			(acq->FI_ERROR)[i], 
			(acq->FI)[i+4],
			(acq->FI_ERROR)[i+4] );
		}

	    /* pupil tracker */
	case 'C':
	    fprintf(ACQdata, "\n"
		    "# ****************************  begin pupil tracker  data  ****************************\n"
		    "#The pupil tracker reference spots in pixel (x, y centriods in total 34 rows and  errors)\n"
		    "#S.No     UT1 (1 & 2 columns)       UT2 (3 & 4 columns)          UT3 (5 & 6 columns)           UT4(7 & 8 columns)\n");
	    
	    for(i=0; i<8; i++)
		{
		fprintf(ACQdata,"%02d %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e  %14.7e \n", 
			i+1, 
			(acq->PT_REF)[i] , 
			(acq->PT_REF_ERROR)[i],
			(acq->PT_REF)[i+PT_refernce_spot_num],
			(acq->PT_REF_ERROR)[i+PT_refernce_spot_num],
			(acq->PT_REF)[i+2*PT_refernce_spot_num],
			(acq->PT_REF_ERROR)[i+2*PT_refernce_spot_num],
			(acq->PT_REF)[i+3*PT_refernce_spot_num],
			(acq->PT_REF_ERROR)[i+3*PT_refernce_spot_num]);
		}
	    break;
	    

	case 'P':
	    fprintf(ACQdata," \n"
		    "#Lateral  and longitudinal pupil positions x, y & z are 1, 2 & 3 rows respectively \n"
		    "units: meter; +ve for diverging beam, -ve for converging beam \n"
		    "#S.No      UT1           UT2          UT3              UT4\n");
	    
	    for(i=0; i<3; i++){
	    fprintf(ACQdata, "%02d %+14.7e %+14.7e %+14.7e %+14.7e \n", i+1,
		    (acq->PT_POSITION)[i], 
		    (acq->PT_POSITION)[i+3], 
		    (acq->PT_POSITION)[i+2*3], 
		    (acq->PT_POSITION)[i+3*3] );
	    }
	    
	    fprintf(ACQdata," \n"
		    "#Lateral  and longitudinal pupil position errors x, y & z are 1, 2 & 3 rows respectively \n"
		    "#S.No      UT1           UT2          UT3              UT4\n");
	    
	    for(i=0; i<3; i++)
		{
		fprintf(ACQdata, "%02d %+14.7e %+14.7e %+14.7e %+14.7e \n", i+1,
			(acq->PT_POSITION_ERROR)[i], 
			(acq->PT_POSITION_ERROR)[i+3], 
			(acq->PT_POSITION_ERROR)[i+2*3], 
			(acq->PT_POSITION_ERROR)[i+3*3] );
		}
	    
	    fprintf(ACQdata, "\n"
		    "#Reference spot position (5th spot) for each UT (x & y centriods) (in pixels)\n"
		    "#S.No          UT1           UT2            UT3           UT4   \n");
	    
	    for(i=0; i<PT_refernce_spot_num; i++)
		{
		fprintf(ACQdata,"%02d  %+14.7e %+14.7e %+14.7e %+14.7e \n", 
			i+1, 
			(acq->PT_REF_SPOT)[i],
			(acq->PT_REF_SPOT)[i+PT_refernce_spot_num],
			(acq->PT_REF_SPOT)[i+2*PT_refernce_spot_num],
			(acq->PT_REF_SPOT)[i+3*PT_refernce_spot_num]);
		}
	    fprintf(ACQdata, "\n"
		    "#Reference spot position error (5th spot) for each UT (x & y centriods) (in pixels)\n"
		    "#S.No          UT1           UT2            UT3           UT4   \n");
	    
	    for(i=0; i<PT_refernce_spot_num; i++)
		{
		fprintf(ACQdata,"%02d  %+14.7e %+14.7e %+14.7e %+14.7e \n", 
			i+1, 
			(acq->PT_REF_SPOT_ERROR)[i],
			(acq->PT_REF_SPOT_ERROR)[i+PT_refernce_spot_num],
			(acq->PT_REF_SPOT_ERROR)[i+2*PT_refernce_spot_num],
			(acq->PT_REF_SPOT_ERROR)[i+3*PT_refernce_spot_num]);
		}
	    break;

	    /* aberration sensor */
	case 'H':
	    fprintf(ACQdata, "\n" 
		    "#********************************* begin aberration sensor data ***********************************\n"
		    "#The aberration sensor reference case centriod data for  68 spots (136 slopes for each UT; pixels)  \n");
	    
	    fprintf(ACQdata,"#S.No          UT1           UT2            UT3           UT4   \n");
	    for(i=0; i<ABSSpotsNum*2; i++)
		{
		fprintf(ACQdata,"%03d %+14.7e %+14.7e %+14.7e %+14.7e \n", 
			i+1, 
			(acq->ABS_REF)[i] , 
			(acq->ABS_REF)[i+ABSSpotsNum*2],
			(acq->ABS_REF)[i+2*ABSSpotsNum*2],
			(acq->ABS_REF)[i+3*ABSSpotsNum*2]);
	    }
	    
	    break;
	    
	case 'Z':
	    fprintf(ACQdata, "\n"
		    "#The aberration sensor fitted Zernike Coefficients (68 measured here for each UT; defocus aberration case here; Units: lambda)  \n"
		    "#S.No          UT1           UT2            UT3           UT4   \n");
	    for(i=0; i<ZernikeNum; i++)
		{
		fprintf(ACQdata,"%02d %+14.7e %+14.7e %+14.7e %+14.7e \n", 
			i+1, 
			(acq->ABS_ZERNIKE)[i] , 
			(acq->ABS_ZERNIKE)[i+ZernikeNum],
			(acq->ABS_ZERNIKE)[i+2*ZernikeNum],
			(acq->ABS_ZERNIKE)[i+3*ZernikeNum]);
		}
	    

	    fprintf(ACQdata, "\n"
		    "#The aberration sensor fitted Zernike Coefficient error (68 measured here for each UT; defocus aberration case here; Units: lambda)  \n"
		    "#S.No          UT1           UT2            UT3           UT4   \n");
	    for(i=0; i<ZernikeNum; i++)
		{
		fprintf(ACQdata,"%02d %+14.7e %+14.7e %+14.7e %+14.7e \n", 
			i+1, 
			(acq->ABS_ZERNIKE_ERROR)[i] , 
			(acq->ABS_ZERNIKE_ERROR)[i+ZernikeNum],
			(acq->ABS_ZERNIKE_ERROR)[i+2*ZernikeNum],
			(acq->ABS_ZERNIKE_ERROR)[i+3*ZernikeNum]);
		}
	    break;

/* field imager */
	case 'I':
	    fprintf(ACQdata, "\n" "*********************** begin field imager data ********************************\n"
		    "#The field imager: brightest object position and its error for each UT (in pixels) \n"
		    "T.No Postion(X)   Position error(X)   Position(Y)    Position error(Y)\n");
	    
	    for(i=0; i<4;i++)
		{
		fprintf(ACQdata,"%d %14.7e %14.7e %14.7e  %14.7e\n", i+1, 
			(acq->FI)[i], 
			(acq->FI_ERROR)[i], 
			(acq->FI)[i+4],
			(acq->FI_ERROR)[i+4] );
		}
	    break;

/* nothing */
	case 'N':
	    break;

/* sky */
	case 'S':
        
	    break;
/* Dead pixel */

	case 'D':
	    break;
	    
/* Flat */
	case 'F':
	    
	    
	    break;
	    
	}

    fclose(ACQdata);

    return error;
}
/* ----------------------------------------------------------------------*/


/* ----------------------------------------------------------------------*/
/**
 * @internal
 * @brief  Write/updates the database with the input attribute
 * @param  attribute name of the attribute. Example: DeadPixelMap
 * @param  pointer to the value of the attribute. Example: pointer to the DeadPixelMap.fits
 * @return error
 */

cpl_error_code gvacqWrite2Db(char * attribute, double * value)
{
    cpl_error_code   error=CPL_ERROR_NONE;

    return error;
}
/* End of gvoacqWrite2Db*/



