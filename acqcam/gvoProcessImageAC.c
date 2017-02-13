 /************************************************************************
 *   This file is part of the E.S.O. - VLTI GRAVITY project
 *
 *  "@(#) $Id: gvoProcessImageAC.c 288965 2016-10-09 12:58:58Z aramirez $
 *   
 *   who       when        what
 *   --------  ----------  ----------------------------------------------
 *   feisenha 2016-10-08  added correct signs for kmirror and rotator
 *   amorin   2016-10-08  updated N and Zenith calculation
 *   narsi    2016-05-13  The conversion of ADR angle to gvrtdACQ pixel coordinates done using ZenithAngle
 *   narsi    2016-05-06  Computing proper Zenith and North angles; taken from gvrtdACQ.tcl
 *   ekw      2016-04-28  Save image separate depending on mode
 *   narsi    2016-03-26  Added new database points and its filling functionality
 *   narsi    2016-03-23  Edited  fillPT_REF, created fillABS_TARGET_POS_AND_FLUX
 *                         for database filling of flux of spots of PT and ABS 
 *   narsi    2016-03-19  Added database filling, set and get functions for defocus estimation
 *   rdembet  2016-03-19  Added OLDB writing at the very end to trigger volac
 *   narsi    2016-03-18  Added database setup functions for defocus estimation
 *   narsi    2016-03-18  Added defocus estimation for calib and Aberration Sensor
 *   feisenha 2016-03-18  Added INS ID and Changed DPR CATG to TEST for saveimg
 *                        Added dbWrite of filename to trigger volac
 *   ekw      2016-03-14  Added all CB and functions for SAVEIMG command
 *                        Merge with nasis changes 
 *   feisenha 2016-02-25  Added MasterSkyABS
 *   feisenha 2015-10-31  Adjusted sign of angles used in ADR calculation
 *   narsi    2015-08-28 Changed order of FI and Pupil to reduce delay for FI
 *   narsi    2014-07-30 VLT2014 adapted
 *   narsi    2015-01-20 ADR code is included
 *   narsi    2014-07-30 VLT2014 adapted
 *   narsi    2014-03-05  add set functions for database coordinates reading
 *   narsi    2014-03-05  add correlation centrioding library
 *   narsi    2014-01-05  update field imager and clean old funcs
 *   ekw      2014-01-20  Master clean up of all interfaces
 *   ekw      2014-02-05  Adopt to Shared Memory access
 *   narsi    2013-12-17  add set coordinates for ABS, PT
 *   narsi    2013-12-17  PT ref spots are measuring from 5th spot data
 *   narsi    2013-12-17  abs upgrade to Galatic Center/extended images
 *   ekw      2013-11-15  Error handling FieldImager and format code acording to standard
 *   ekw      2013-11-12  add enable architecture
 *   narsi    2013-11-04  upgrade the field imager and error handling
 *   ekw      2013-10-02  Add setCoordinates
 *   narsi    2013-09-19  data base filling corrections 
 *   narsi    2013-09-16  K-Mirror rotation for pupil tracker
 *   ekw      2013-06-28  data base filling
 *   narsi    2013-06-12  Zernike approach for pupil tracker
 *   narsi    2013-04-26  field imager for Galactic field
 *   narsi    2013-03-15  Modified i/p and o/p arguments of funcs 
 *                         and defined new stubs for DeadPixelMap, Sky, Flat
 *   narsi    2013-01-29  gvacqAberrationSensor
 *   narsi    2012-11-20  gvacqPupilTracker
 *   narsi    2012-10-27  created
 */

/****************************************************************************
 *   NAME 
 *   gvoacqProcessImageAC - procedere to compute acquistion camera functions
 *  
 *   SYNOPSIS
 *   #include <gvoProcessImageAC.h>
 * 
 *   DESCRIPTION
 *   Given a pointer to the acquisition camera image, will update the 
 *   instrument database with relevant quantities.
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

/*
 * System Headers
 */

/**
 * @defgroup process gravity observation software for acquisition camera
 * @ingroup image_analysis
 *
 * This module provides functions for acquision camera data reduction
 *
 * @par Synopsis:
 * @code
 * @endcode
 */

/*-----------------------------------------------------------------------------
 Includes
 ----------------------------------------------------------------------------*/

#include "gvoProcessImageAC.h" /* */
#include "cpl_apertures.h"
#include <math.h>
#include "tims.h"
#include <string.h>    /* memcpy ... */
#include <ccs.h>
#include <db.h>

/*-----------------------------------------------------------------------------
 Defines
 ----------------------------------------------------------------------------*/

#define POSIX_SOURCE 1
/* #define VLT2011 */
#define VLT2014  


#ifdef VLT2011
struct _cpl_apertures_ {
    /* Number of apertures */
    int         naperts;

    /* All positions are in FITS conventions : (1,1) is the lower left) */
    /* Aperture center: x=sum(x*1)/sum(1) */
    double  *   x;
    double  *   y;
    /* Aperture weighted center : xcentroid=sum(x*f(x)) / sum(f(x)) */
    double  *   xcentroid;
    double  *   ycentroid;

    int     *   npix;
    int     *   left_x;
    int     *   left_y;
    int     *   right_x;
    int     *   right_y;
    int     *   top_x;
    int     *   top_y;
    int     *   minpos_x;
    int     *   minpos_y;
    int     *   maxpos_x;
    int     *   maxpos_y;
    int     *   bottom_x;
    int     *   bottom_y;
    double  *   max_val;
    double  *   min_val;
    double  *   mean;
    double  *   varsum;
    double  *   median;
    double  *   stdev;
    double  *   flux;
};
typedef long long  cpl_size;
#endif


#define ZernikeNum 28   
#define ABSNUMSPOTS 77


static int ADRenabled = 0;
static int enabled=0; 
static int IMAGE_X=2048;/* pixes */
static int IMAGE_Y=1536;/* pixes */
static double ABSiStart[ABSNUMSPOTS*2*4];
static double ABSRefSpotsCenters[ABSNUMSPOTS*2*4];
static double ABSDefocusFactory[4*4];
static double RefDefocusPositionsMeasured[4*4];
static int ABSwindow[8]={218,681,1155,1612,978,976,981,969}; /* pixes coordinates */
static int ABSwindowSize=176; /* pixes */
static int ABS9x9Slopes=ABSNUMSPOTS*2;
static int PRINTTEXT=0;
static int PRINTMSG=0;


short store_sky = 0;
short save_image_field = 0;
short save_image_pupil = 0;
short save_image_aberr = 0;
char file_info [256];

int Extended_source = 0;
cpl_image * MasterFlat = NULL;
cpl_image * MasterSkyPupil = NULL;
cpl_image * MasterSkyField = NULL;
cpl_image * MasterSkyABS = NULL;
cpl_mask * DeadPixelMap = NULL;

cpl_error_code eval_pos_with_rotationAND5thspot(double *, double, double *);
cpl_error_code extract_PT_spot_positions(double, double*, double*);
cpl_error_code coordinate_rotate(double*, double*, double, double *);

double PT5thSpotcoordinates[32];
/*DB*/
static double Spot5RefData[32]; /* 4 Tel. 4 (FC spots) *2(x,y) */
static double FIfwhm_minthreshold;
static int Check_iteration=50;
static double ZenithDistance; /* degrees, taken as 90-altitude  */
double ZenithAngle[4]=   {0.,0.,0.,0.0};    
static double Temperature; /* degree centigrade  */
static double Pressure; /* units Millibars  */
static double Humidity; /* in % */
static double ParAngleIn[4]={0,0,0,0};/* degrees  */
static double FIPixelScale;/* arcsec  */
static double PosAngleIn[4] = {0., 0, 0,0}; /* degrees */
static double KmirrorRotation[4]; /* Steps  */
static double HKmagnitude;
static int ABSLensletSize=6;

char AB_inv_Z2S_68_136Fits[256];

int FIErrorStatus[4];
int PTErrorStatus[4];
int ABSErrorStatus[4];
double PTFitStateErrorThreshold=10;
int PTCalibErrorStatus[4];
int ABSCalibErrorStatus[4];
double PolyACQ[11]={1.52385,0.0348105,-0.00569223,0.000729958,3.02968e-05,-2.37115e-05,
		    1.28736e-06,2.88224e-07,-4.86121e-08,6.68202e-09,-5.14515e-10}; 
double PolyBC[11]={2.2, 0, 0, 0, 0 ,0 ,0 ,0 ,0 ,0,0};



/*-----------------------------------------------------------------------------
 Private function prototypes
 ----------------------------------------------------------------------------*/
void copyAbsRefPosition(gvacqALL_DATA *dbCopy);



/* ----------------------------------------------------------------------*/
/**
 * @internal 
 * @brief  Estimates the pupil, field, wavefront correction parameters
 * @param  invertFlag used for Blink mode, in the process of removing background pupil tracker window
 * @param  acMode acquisition camera guiding mode
 * @param  allDataVars parameter to write the all measured values to the database
 * @return status error 
 * @see    gvacqFieldImager, gvacqPupilTracker, gvoacqAberrationSensor
 * This function takes the query string 
 * start i.e. Does preprocess, pupil tracker, aberration sensor, field imager
 * 
 * It can be called with different modes, those are set by func setACmode
 * 1. "ACMODE_FIELD" = Field tracking
 * 2. "ACMODE_PUPIL" = Pupil tracking
 * 3. "ACMODE_ABERR"= Shack-Hartmann wavefront coefficients in Zernike poly
 * 4. "ACMODE_FIBREF" = computation of pupil tracker reference positions 
 * 5. "ACMODE_ABBCAL" = computation of Shack-Hartmann reference positions 
 * 6. "ADRenabled" = computes star drift caused by atmospheric differential refraction 
 * Please note above are global variables defined in header.
 *
 * After computation of tracking parameters, these are updated to the database and
 * also it can be written to the text file by selecting PRINTEXT=1 above to this function
 * 
 * If theree modes ACMODE_FIELD, ACMODE_PUPIL, ACMODE_ABERR are called it executes
 * field tracking and its dabase update, pupil tracking and its database update and aberration measurement
 * and database update respectively.
 *
 * When there is lot of background in the pupil tracker window, Blink mode is executed by using invertFlag=1
 * For this usage the pupil tracker windows are already subtracted by background, have a look at func storeSky
 *
 * Also note that, MasterFlat, MasterSkyPupil, MasterSkyField, MasterSkyABS, 
 * DeadPixelMap are passed by global variables 
 *
 * For database reading from setup, it is recommnded to see individual functions  
 * gvacqFieldImager, gvacqPupilTracker and 
 * gvoacqAberrationSensor.
 *
 * @par Error Handling:
 * Possible #_cpl_error_code_ se in this function:
 * CPL_ERROR_NONE initialise error parameter as none
 */

extern int gvoacqProcessImageAC(cpl_image * DetPointer, int invertFlag, int acMode,
			gvacqALL_DATA* allDataVars) {

	/*Test input parameter set */
	if (DetPointer == NULL ) {
        printf("ERROR:: -> Detector buffer submitted is NULL at gvoacqProcessImageAC. \
Error occured at %s:%s:%d\n", 
               __FILE__, __FUNCTION__, __LINE__);	
        return cpl_error_get_code();
	}
	
	if(cpl_image_get_type(DetPointer) != CPL_TYPE_DOUBLE){
	printf("ERROR: The software works only with type double image. Detected %u at %s:%s:%d \n", 
	       cpl_image_get_type(DetPointer), __FILE__, __FUNCTION__, __LINE__);
	
	return cpl_error_get_code();
	}
    
	int size_x=cpl_image_get_size_x(DetPointer);
	int size_y=cpl_image_get_size_y(DetPointer);
	
	if(size_x != IMAGE_X || size_y != IMAGE_Y) {
        printf("ERROR::-> Image size detected is %d x %d, but expected 2048x1536. Detected at %s:%s:%d \n", 
               size_x, size_y, __FILE__, __FUNCTION__, __LINE__);
        
        return cpl_error_get_code();
	}
	
        
	
	/* Intialising the memory for parameters
	 */
	double ObjectPosition[4 * 2]; /* 4 Tel.  and 1 object centriods (x &y) */
	double ObjectPosition2[4 * 2]; /* 4 Tel.  and 1 object centriods (x &y) for star 2 */
	double ObjectPositionError[4 * 2]; /* 4 Tel.  and 1 object centriod errors (x &y) */
	double FiFWHM[4]; /* 4 Tel FI FWHM */
	double FiFlux[4]; /* 4 Tel FI Flux */
	
	double PTRefPosition[4 * 32]; /* 4 Tel.  and 16 spots (32 centriods x &y) */
	double PTRefPositionError[4 * 32]; /* 4 Tel.  and 16 spots error (32 centriods x &y) */
	double PupilPosition[4 * 3]; /* 4 Tel.  lateral (x & y) and longitudinal pupil positions */
	double PupilPositionError[4 * 3]; /* 4 Tel.  lateral (x & y) and longitudinal pupil position erros */
	double barycentre[4 * 8]; /* 4 Tel.  and 4 barycentre  (8 centriods x &y) */
	double barycentre_error[4 * 8]; /* 4 Tel.  and 4 barycentre error (8 centriods x &y) */
	double PTSpotsFlux[4*32]; /* 4 Tel.  and 16 spots (32 centriods x &y) */

	double ABSRefPosition[4*ABSNUMSPOTS*2]; /* 4 Tel.  and 77 spots (154 slopes x &y) */
	double ABSRefPositionError[4 * ABSNUMSPOTS*2]; /* 4 Tel.  and 68 spots (136 slopes error x &y) */
	double ZernikeCoeff[4 * ZernikeNum]; /* 4 Tel. 68 Zernike coefficients */
	double ZernikeCoeffError[4 * ZernikeNum]; /* 4 Tel. 68 Zernike coefficients error */
	double AutoCorrelationDefocus[4]; /* For Target wavefront or image 4 Tel. */
	double TargetDefocusPositions[16]; /* target auto-correlation spot positions */
	double Spot5RefData_error[4 * 8]; /* FC coordinates error */
	int ABSLensletSizeArray[4]; /* used/calculated ABSLensletSizes */
	double ABSSpotsCurrentPosition[4*ABSNUMSPOTS*2]; /* 4 Tel.  and 77 spots (154 slopes x &y) */
	double ABSSpotsFlux[4*ABSNUMSPOTS]; /* 4 Tel.  and 16 spots (32 centriods x &y) */
	double ABSSlopes[ABS9x9Slopes*4];

	int i,j; 
	cpl_image * DetPointerPupil = NULL;

	cpl_error_code error = CPL_ERROR_NONE;
	cpl_error_code ADRerror = CPL_ERROR_NONE;
	struct acqStruct ACQData;

	/*----------------------------------------------------------------------*/

	/* Treat STORSKY command*/
	/* storage to disk anymore */

	if (store_sky) { 
        storeSky(DetPointer, store_sky);
        store_sky = 0;
	}
    
	if (acMode & ( ACMODE_PUPIL | ACMODE_FIBREF )) {
		DetPointerPupil = cpl_image_duplicate(DetPointer); /* must be deleted before return */
	}
    
	/* preform Preprocessing if ENABLED */
	if (enabled) 
	{
		if (acMode & ACMODE_FIELD)
		{
			error = gvoacqPreprocess(DetPointer, MasterFlat, MasterSkyField, DeadPixelMap);
			if (error != 0)
			{
				cpl_msg_error(cpl_func, "Preprocessing failed, error code: %d\n",error);
			}

			if (save_image_field) {
			  SaveImage(DetPointer,"FIELD");
			  save_image_field = 0;
			}

		}            
		if (acMode & ( ACMODE_PUPIL | ACMODE_FIBREF ))
		{
			error = gvoacqPreprocess(DetPointerPupil, MasterFlat, MasterSkyPupil, DeadPixelMap);
			if (PRINTMSG) printf("Background substracted for pupil"); 
			if (error != 0)
			{
				printf("Preprocessing failed, error code: %d\n",error);
			}
			if (invertFlag) {
				error = cpl_image_multiply_scalar(DetPointerPupil,-1.0);
				if (error != 0)
				{
					printf("Inverting sky image failed, error code: %d\n",error);
				}
				/*		  cpl_image_save(DetPointerPupil, 
						  "PreProcess_Output_Pup_I.fits", 
						  CPL_BPP_IEEE_DOUBLE, NULL, CPL_IO_CREATE);*/
			} else {
				/*		  cpl_image_save(DetPointerPupil, 
						  "PreProcess_Output_Pup_N.fits", 
						  CPL_BPP_IEEE_DOUBLE, NULL, CPL_IO_CREATE);*/
			}
			

			if (save_image_pupil) {
			  SaveImage(DetPointerPupil,"PUPIL");
			  save_image_pupil = 0;
			}

		}
		if (acMode & ( ACMODE_ABERR | ACMODE_ABBCAL ))
		{
			error = gvoacqPreprocess(DetPointer, MasterFlat, MasterSkyABS, DeadPixelMap);
			if (error != 0)
			{
				cpl_msg_error(cpl_func, "Preprocessing failed, error code: %d\n",error);
			}

			if (save_image_aberr) {
			  SaveImage(DetPointer, "ABERR");
			  save_image_aberr = 0;
			}

		}            


	}/*End of enabled */

	/*----------------------------------------------------------------------*/
	/*
	 *  MAIN EXECUTION BLOCK depending on MODE
	 *  old 'P' , new ACMODE_PUPIL : PUPILTRACK
	 *  old 'Z' , new ACMODE_FIELD : ZERNICKE PARAMETER
	 *  old 'I' , new ACMODE_ABERR : FIELD TRACKING mode
	 */
    
	for (i=0; i<gvacqNUM_REC_TELESCOPES; i++) {
        FIErrorStatus[i] = 999;
        PTErrorStatus[i] = 999;
        ABSErrorStatus[i] = 999;
        PTCalibErrorStatus[i] = 999;
        ABSCalibErrorStatus[i] = 999;
	}

	/** 
	 * ----------------------------------------------------------------------
	 *  Field Tracking mode
	 * ------------------------------------------------------------------------
	**/
	if (acMode & ACMODE_FIELD) {
        if(PRINTMSG) cpl_msg_info(cpl_func, "set T, P, H, Z, ParalacticAng PosAng FIPixelScale %f %f %f %f %f %f %f\n",  
				  Temperature, Pressure, Humidity, ZenithDistance, 
				  ParAngleIn[0], PosAngleIn[0],FIPixelScale );

	error = gvacqFieldImager(DetPointer, FIfwhm_minthreshold, ObjectPosition, ObjectPosition2,
				 ObjectPositionError, FiFWHM, FiFlux, FIErrorStatus);
        
            if(PRINTMSG)printf("Tel4 FiFWHM %f FiFlux %f\n",  FiFWHM[3],FiFlux[3]);

	    gvacq_error_print(error, "FieldImager", __LINE__ - 1, __FILE__);
	    if(PRINTMSG){
            printf(">> FI positions: %f %f %f %f %f %f %f %f \n", ObjectPosition[0], ObjectPosition[1],  
                   ObjectPosition[2],  ObjectPosition[3],  ObjectPosition[4], ObjectPosition[5], 
                   ObjectPosition[6],  ObjectPosition[7]);
            printf(">> FI positions(2): %f %f %f %f %f %f %f %f \n", ObjectPosition2[0], ObjectPosition2[1],  
                   ObjectPosition2[2],  ObjectPosition2[3],  ObjectPosition2[4], ObjectPosition2[5], 
                   ObjectPosition2[6],  ObjectPosition2[7]);
	    }
        
	   
        
	    ACQData.FI = ObjectPosition;
	    ACQData.FI2 = ObjectPosition2;
	    ACQData.FI_ERROR = ObjectPositionError;
	    error = fillFI(allDataVars, ObjectPosition, ObjectPosition2, ObjectPositionError, FiFWHM, FiFlux);
/* ADR code */	

        double FIXYoffsetPX[8]; 
        double FIXYoffsetArcSec;
        double EffectiveLambdaa[2]={0.0, 0.0};
	double deg2rad = 0.0174532925199;
        ADRerror = ADRCorrPixelsXY(HKmagnitude, PolyACQ, PolyBC, Pressure, Temperature, 
                                   Humidity, ZenithDistance, PosAngleIn, ParAngleIn, 
                                   FIPixelScale, EffectiveLambdaa, &FIXYoffsetArcSec, FIXYoffsetPX);
	

	/* FIXYoffsetPX re-calculating from gvrtdACQ.tcl formula using the ADR in arcsec */
	for(i=0; i<4; i++){
	 if(FIPixelScale !=0){
	   FIXYoffsetPX[i]= cos(ZenithAngle[i]*deg2rad)*FIXYoffsetArcSec/FIPixelScale;
/* aiginst to intutive acq rtd; -1 sign is to match the gvrtdACQ.tcl directions*/
	   FIXYoffsetPX[4+i]= -sin(ZenithAngle[i]*deg2rad)*FIXYoffsetArcSec/FIPixelScale; 
	 }else{
	   FIXYoffsetPX[i]=0;
	   FIXYoffsetPX[i+4]=0;
	 }
	}

	    if(error ==0) {
            ACQData.ADR_pixels = FIXYoffsetPX;	
            ACQData.ADR_arcsec = &FIXYoffsetArcSec;
            ACQData.EffectiveLambda=EffectiveLambdaa;
	    } else {
            printf("ERROR::-> in ADR calculation \n");
	    }
	    if(PRINTMSG) printf(">> H-K=%f,  ADR lambda[0,1]=%f, %f, angle and corrections[0,1]  %f %f %f \n",
				HKmagnitude, EffectiveLambdaa[0], EffectiveLambdaa[1], 
				FIXYoffsetArcSec, FIXYoffsetPX[0],  FIXYoffsetPX[1]);
        
	    if (error != CPL_ERROR_NONE) {
            cpl_msg_error(cpl_func, "Detected at %s:%s:%d\n", __FILE__, __FUNCTION__,
                          __LINE__);
	    }
	    if (PRINTTEXT) error = gvacqWrite2File(&ACQData, "FI");
        

            if (PRINTMSG) { /* AAMORIM debug */
	        printf("ADR Test T1:%f %f %f %f \n", ObjectPosition[0],ObjectPosition[4],FIXYoffsetPX[0],FIXYoffsetPX[4]); 
    	        printf("ADR Test T2:%f %f \n", ObjectPosition[1],ObjectPosition[5]); 
	        printf("ADR Test T3:%f %f \n", ObjectPosition[2],ObjectPosition[6]); 
	        printf("ADR Test T4:%f %f \n", ObjectPosition[3],ObjectPosition[7]); 
	        printf("ADR Test T1:%f %f \n", ObjectPosition2[0],ObjectPosition2[4]); 
	        printf("ADR Test T2:%f %f \n", ObjectPosition2[1],ObjectPosition2[5]); 
	        printf("ADR Test T3:%f %f \n", ObjectPosition2[2],ObjectPosition2[6]); 
	        printf("ADR Test T4:%f %f \n", ObjectPosition2[3],ObjectPosition2[7]); 
            }

	    if(ADRenabled){
            for(i=0; i<4; i++){
	       ObjectPosition[i]=ObjectPosition[i]-FIXYoffsetPX[i];
	       ObjectPosition[i+4]=ObjectPosition[i+4]-FIXYoffsetPX[4+i];
	       ObjectPosition2[i]=ObjectPosition2[i]-FIXYoffsetPX[i];
	       ObjectPosition2[i+4]=ObjectPosition2[i+4]-FIXYoffsetPX[4+i];
            }
	    }
/* 
   ADR offset pixels are stored for only telescope in database
   In the gvrtdtcl, currenly the adr plot is plotted with objectposition=FIpostion-adr
 */
	    double FIXYoffsetPXNew[2];
	    FIXYoffsetPXNew[0]=FIXYoffsetPX[0];
	    FIXYoffsetPXNew[1]=FIXYoffsetPX[4];
	    /* NEW OBJECT POSITION ADR CORRECTED */
	    error = fillFIobj(allDataVars, ObjectPosition, ObjectPosition2, ObjectPositionError, EffectiveLambdaa, 
			      &FIXYoffsetArcSec, FIXYoffsetPXNew, HKmagnitude);
	    if(PRINTMSG){
                   printf(">> FI positions + ADR: %f %f %f %f %f %f %f %f \n", 
		   ObjectPosition[0], ObjectPosition[1],  
                   ObjectPosition[2],  ObjectPosition[3],  
		   ObjectPosition[4], ObjectPosition[5], 
                   ObjectPosition[6],  ObjectPosition[7]);
	    }
	    if(PRINTMSG){
                   printf(">> FI positions2 + ADR: %f %f %f %f %f %f %f %f \n", 
		   ObjectPosition2[0], ObjectPosition2[1],  
                   ObjectPosition2[2],  ObjectPosition2[3],  
		   ObjectPosition2[4], ObjectPosition2[5], 
                   ObjectPosition2[6],  ObjectPosition2[7]);
	    }

	    
	    /* Dubugging messages */
	if(PRINTMSG){
           if(ADRenabled !=0) for(i=0; i<4; i++) FIErrorStatus[i] += ADRerror;
           cpl_msg_info(cpl_func, "FIErrorStatus= %d %d %d %d  \n",
                     FIErrorStatus[0], FIErrorStatus[1], FIErrorStatus[2], FIErrorStatus[3]);	  
	}
	
	for (i=0; i<gvacqNUM_REC_TELESCOPES; i++) {
          if(ADRenabled) { 
	      (*allDataVars->acFiTrSts)[i] = FIErrorStatus[i]+ADRerror;
          }else {
	      (*allDataVars->acFiTrSts)[i] = FIErrorStatus[i];		
          }
	}

	}/* End of acMode & ACMODE_FIELD */
  


	/** 
	 * ----------------------------------------------------------------------
	 *  Pupil Tracking mode
	 * ------------------------------------------------------------------------
	**/
	if (acMode & ACMODE_PUPIL) {
        /* Initialize array with zeros */	
	for(i=0; i<32*4; i++) PTRefPosition[i]=0.0;
	for(i=0; i<16*4; i++) PTSpotsFlux[i]=0.0;
	for(i=0; i<8*4; i++) barycentre[i]=0.0;

        error = gvacqPupilTracker(DetPointerPupil,
				  Spot5RefData,
				  FIfwhm_minthreshold,
				  PTFitStateErrorThreshold,
                                  PupilPosition,
                                  PupilPositionError, 
                                  barycentre, 
                                  barycentre_error, 
                                  PTRefPosition,
                                  PTRefPositionError,
				  PTSpotsFlux,
				  PTErrorStatus);
	/*
	for(i=0; i<12; i++){
	if(fabs(PupilPosition[i]) > 100){
	if(PTErrorStatus[0]==0 && i <3)              PTErrorStatus[0]=-3;
	if(PTErrorStatus[1]==0 && (i >=3 || i<6) )   PTErrorStatus[1]=-3;
	if(PTErrorStatus[2]==0 && (i >=6 || i<9) )   PTErrorStatus[2]=-3;
	if(PTErrorStatus[3]==0 && (i >=9 || i<12) )  PTErrorStatus[3]=-3;
	}
	}
	*/
        

	if(PRINTMSG) printf("%f %f %f %f %f %f %f %f %f %f %f %f \n", 
			    PupilPosition[0],  PupilPosition[1],  PupilPosition[2],  
			    PupilPosition[3],  PupilPosition[4], PupilPosition[5], 
			    PupilPosition[6],  PupilPosition[7], PupilPosition[8], 
			    PupilPosition[9],  PupilPosition[10], PupilPosition[11]);

       
        
        /* Write to a txt file or database */
        ACQData.PT_POSITION = PupilPosition;
        ACQData.PT_POSITION_ERROR = PupilPositionError;
        /* next two lines only needed for gvacqWrite2File "PT" mode */
        ACQData.PT_REF_SPOT = PTRefPosition;
        ACQData.PT_REF_SPOT_ERROR = PTRefPositionError;
        /*Write pupil shifts measured*/
        error = fillPT_POSITION(allDataVars, PupilPosition, PupilPositionError);
        /*Write pupil tracker spots positions (16(x,y) spots for each Tel.) */
        error = fillPT_REF(allDataVars, PTRefPosition,PTRefPositionError,PTSpotsFlux);
        /*Write barycenter positions for each Tel.) */
       error = fillPT_BARYCENTER(allDataVars, barycentre, barycentre_error);
        
        if (error != CPL_ERROR_NONE) {
            cpl_msg_error(cpl_func, "Detected at %s:%s:%d\n", __FILE__, __FUNCTION__,
                          __LINE__);
	  }
        if (PRINTTEXT) error = gvacqWrite2File(&ACQData, "PT");

	if(PRINTMSG){
        cpl_msg_info(cpl_func, "PTErrorStatus=%d %d %d %d   \n",
                     PTErrorStatus[0], PTErrorStatus[1], PTErrorStatus[2], PTErrorStatus[3] );	  
	}
	
	for (i=0; i<gvacqNUM_REC_TELESCOPES; i++) {
        (*allDataVars->acPuTrSts)[i] = PTErrorStatus[i];	
	}

	}/*End of  acMode & ACMODE_PUPIL */
    

           /** 
	    * ----------------------------------------------------------------------
	    *  Aberration Tracking mode
	    * ------------------------------------------------------------------------
	    **/

	if (acMode & ACMODE_ABERR) {
        
        copyAbsRefPosition(allDataVars);

        for(i=0; i<4*ZernikeNum; i++) ZernikeCoeff[i]=0.0;
	for(i=0; i<4*ABSNUMSPOTS*2; i++) ABSSpotsCurrentPosition[i]=ABSRefSpotsCenters[i];
	for(i=0; i<4*ABSNUMSPOTS; i++) ABSSpotsFlux[i]=0.0;

        error = gvoacqAberrationSensor(DetPointer, ABSwindow, ABSwindowSize, 
				       ABSRefSpotsCenters, AB_inv_Z2S_68_136Fits, 
				       Extended_source, RefDefocusPositionsMeasured, 
				       &ABSLensletSize, ZernikeCoeff, ZernikeCoeffError, 
				       AutoCorrelationDefocus, TargetDefocusPositions, 
				       ABSLensletSizeArray, ABSSpotsCurrentPosition, ABSSlopes, 
				       ABSSpotsFlux, ABSErrorStatus);

 

       

        gvacq_error_print(error, "AberrationSensor", __LINE__ - 1, __FILE__);
        if(PRINTMSG){
            double WavefrontRMS[4]={0.0,0.0,0.0,0.0};
            if (error == CPL_ERROR_NONE) {
                for(j=0; j<4; j++) {
                    for(i=1; i<ZernikeNum; i++){
                        WavefrontRMS[j] += sqrt(ZernikeCoeff[i+ j*ZernikeNum]*ZernikeCoeff[i+j*ZernikeNum]);
                    }
                    /* printf(" Tel. %d, wavefront rms error = %14e lambda \n", j+1, WavefrontRMS[j]/1.65); */
                }    
            }
        }
        /* Write to a txt file or database */
        ACQData.ABS_ZERNIKE = ZernikeCoeff;
        ACQData.ABS_ZERNIKE_ERROR = ZernikeCoeffError;
        *allDataVars->acAbsLensletSize = ABSLensletSize;
        error = fillZERNIKE(allDataVars, ZernikeCoeff, ZernikeCoeffError, 
			    ABSLensletSizeArray, AutoCorrelationDefocus);
	
        if (PRINTTEXT) error = gvacqWrite2File(&ACQData, "ABS"); 
        if (error != CPL_ERROR_NONE) {
            cpl_msg_error(cpl_func, "Detected at %s:%s:%d\n", __FILE__, __FUNCTION__,
                          __LINE__);
        }
	error=fillABS_TARGET_POS_AND_FLUX(allDataVars, ABSSpotsCurrentPosition, ABSSlopes, ABSSpotsFlux);
	
	if (error != CPL_ERROR_NONE) {
	cpl_msg_error(cpl_func, "Detected at %s:%s:%d\n", __FILE__, __FUNCTION__,
		      __LINE__);
        }

	if(PRINTMSG){
        cpl_msg_info(cpl_func, "ABSErrorStatus= %d %d %d %d\n",
                     ABSErrorStatus[0], ABSErrorStatus[1], ABSErrorStatus[2], ABSErrorStatus[3]);	  
	}

	for (i=0; i<gvacqNUM_REC_TELESCOPES; i++) {
        (*allDataVars->acAbTrSts)[i] = ABSErrorStatus[i];
	}

	}/* acMode & ACMODE_ABERR */
    

	/** 
	 * ----------------------------------------------------------------------
	 *  reference fiber spots calinration for pupil tracking
	 * ------------------------------------------------------------------------
	 **/
	
	if (acMode & ACMODE_FIBREF) {
        acMode = acMode & (~ ACMODE_FIBREF); /* single shot mode => clear mode bit */
        cpl_msg_info(cpl_func, "----Compute the fiber coupler reference positions -----");
        
        /*
         * The fiber coupler reference spot corrdinates are the reference for the
         * pupil trcaker to measure the lateral and longitudinal shifts.
         *
         * gvacqPupilTrackerFCspots functions reads the existing coordinates from the datebase
         * and apply gauss fits on the calibration exposure image taken in the lab.
         *
         */
	    
	    
	/* if(PRINTMSG) for(i=0; i<16; i++) printf("INITIAL %f %f\n", Spot5RefData[i], Spot5RefData[16+i]); */
     
     cpl_error_code error = gvacqPupilTrackerFCspots(DetPointerPupil, Spot5RefData, /* Xc[0-15], Yc[15-32] */
                                                     Spot5RefData, /* Xcm[0-15], Ycm[15-32] */
                                                     Spot5RefData_error, PTCalibErrorStatus);
     
     gvacq_error_print(error, "FiberRefPos", __LINE__ - 1, __FILE__);
     if(PRINTMSG) for(i=0; i<16; i++) printf("%f %f\n", Spot5RefData[i], Spot5RefData[16+i]); 
	 
     ACQData.PT_REF = Spot5RefData;  
     ACQData.PT_REF_ERROR = Spot5RefData_error;  
     if (PRINTTEXT) error = gvacqWrite2File(&ACQData, "PTCalib");
     if (error != CPL_ERROR_NONE) {
         cpl_msg_error(cpl_func, "Detected at %s:%s:%d\n", __FILE__, __FUNCTION__,
                       __LINE__);
     }
     error = fillREF_SPOT(allDataVars, Spot5RefData, Spot5RefData_error);
	 
     if (error != CPL_ERROR_NONE)
         cpl_msg_error(cpl_func, "Detected at %s:%s:%d\n", __FILE__, __FUNCTION__,__LINE__);



     if(PRINTMSG){
     cpl_msg_info(cpl_func, "PTCalibErrorStatus= %d %d %d %d \n",
		  PTCalibErrorStatus[0], PTCalibErrorStatus[1], 
		  PTCalibErrorStatus[2], PTCalibErrorStatus[3]);
	}
     
     for (i=0; i<gvacqNUM_REC_TELESCOPES; i++) {
     (*allDataVars->acPtCalSts)[i] = PTCalibErrorStatus[i];
     }
     
	}/* End of acMode & ACMODE_FIBREF */
    

	/** 
	 * ----------------------------------------------------------------------
	 *  Shack-Hartmann reference spots calibration for aberration tracking
	 * ------------------------------------------------------------------------
	 **/
	
	if (acMode & ACMODE_ABBCAL) {
        cpl_msg_info(cpl_func, "---------Computing the aberration sensor reference positions -----");
        
        error = gvoacqAberrationSensorCalib(DetPointer, ABSwindow, ABSwindowSize, ABSiStart,  
					    Extended_source, ABSLensletSize, ABSRefPosition,
                                            ABSRefPositionError, ABSDefocusFactory, 
					    RefDefocusPositionsMeasured);
        /* gvacq_error_print(error, "AberrationCal", __LINE__ - 1, __FILE__); */
        if(PRINTMSG)gvacq_error_print(error, "AberrationCal", __LINE__ - 1, __FILE__);
        
        ACQData.ABS_REF = ABSRefPosition;
        ACQData.ABS_REF_ERROR = ABSRefPositionError;
	    
        if (PRINTTEXT) error = gvacqWrite2File(&ACQData, "ABSCalib");
        error = fillABS_REF(allDataVars, ABSRefPosition, 
			    ABSRefPositionError, RefDefocusPositionsMeasured);

        /*----------------------------------------------------------------------*/	    
        if (error != CPL_ERROR_NONE)
            cpl_msg_error(cpl_func, "Detected at %s:%s:%d\n", __FILE__, __FUNCTION__,__LINE__);
	}/* End of  acMode & ACMODE_ABBCAL */

	
	if (error != CPL_ERROR_NONE) {
		cpl_msg_error(cpl_func, "Detected at %s:%s:%d\n", __FILE__, __FUNCTION__,
                      __LINE__);
	} else{
	        if(PRINTMSG) cpl_msg_info(cpl_func, "Successfull exit from func gvoProcessImage ..");
    }
    
	if (DetPointerPupil != NULL) {
        cpl_image_delete(DetPointerPupil);
	}
	
	
/* 	return cpl_error_get_code(); */
	return CPL_ERROR_NONE;
    
}



/**
 * @brief Passed sky background buffers to local variables for the usage of 
 * field tracking and pupil tracking 
 * @param type STORSKY_FIELD or STORSKY_PUPIL or STORSKY_ABS
 *
 * In a result variables MasterSkyField, MasterSkyPupil are filled
 * This function is called in toplevel function gvacqMAIN_TASK.C:ImageReadyCB
 */
void storeSky (cpl_image *DetPointer, int type) {
    if (type & STORSKY_FIELD) { 
		cpl_msg_info(cpl_func, "Recorded new sky image for field tracking\n");
		
		if(DetPointer == NULL) {
            printf("Detector buffer is null \n");
		}
		
        if (MasterSkyField != NULL ) {
            cpl_image_delete(MasterSkyField);
            MasterSkyField = NULL;
		}
        
		MasterSkyField = cpl_image_duplicate(DetPointer);
    }
    
    if (type & STORSKY_PUPIL) { 
		cpl_msg_info(cpl_func, "Recorded new sky image for pupil tracking\n");
		
        if (MasterSkyPupil != NULL ) {
            cpl_image_delete(MasterSkyPupil);
            MasterSkyPupil = NULL;
		}
		MasterSkyPupil = cpl_image_duplicate(DetPointer);
    }

    if (type & STORSKY_ABS) { 
		cpl_msg_info(cpl_func, "Recorded new sky image for aberration sensor\n");
		
        if (MasterSkyABS != NULL ) {
            cpl_image_delete(MasterSkyABS);
            MasterSkyABS = NULL;
		}
		MasterSkyABS = cpl_image_duplicate(DetPointer);
    }

}/* End of stroreSky */

/**
 * @brief Save ACQ Image
 *
 * In a result variables MasterSkyField, MasterSkyPupil are filled
 * This function is called in toplevel function gvacqMAIN_TASK.C:ImageReadyCB
 */
void saveImage(cpl_image *DetPointer, char* fileinfo) {

  printf(" IN saveImage wit fileinfo : %s", fileinfo);
/*
  if (type & STORSKY_FIELD) {
    cpl_msg_info(cpl_func, "Recorded new sky image for field tracking\n");

    if (DetPointer == NULL) {
      printf("Detector buffer is null \n");
    }

  if (type & STORSKY_PUPIL) {
    cpl_msg_info(cpl_func, "Recorded new sky image for pupil tracking\n");

    if (MasterSkyPupil != NULL) {
      cpl_image_delete(MasterSkyPupil);
      MasterSkyPupil = NULL;
    }
    MasterSkyPupil = cpl_image_duplicate(DetPointer);
  }
*/
}/* End of saveImage */

/**
 * @brief  Quickly scans if there exits any objects in the input image
 * @param  fwhm_threshold Thershold fwhm is used to reject hot pixels and cosmic events
 * @return error
 * 
 * Used for field imager and pupil tracker
 */
cpl_error_code Check_objects_found(cpl_image *FI_image, double fwhm_threshold){
    
    
    cpl_error_code error = CPL_ERROR_NONE;
    int  j, k;
    
#ifdef VLT2014
    cpl_size  px1, py1; 
#else  	 
    int px1, py1;
#endif
    
    
    int im_size_x, im_size_y,  status;
    double fwhm_xx, fwhm_yy;
    
    /*test entries */
    if (FI_image ==NULL ) return -1;
    
    im_size_x = cpl_image_get_size_x(FI_image);
    im_size_y = cpl_image_get_size_y(FI_image);
    status =0;
    
    
    do {
        error = cpl_image_get_maxpos(FI_image, &px1, &py1); /* Get X and Y pixel positions*/
        
        if (error != CPL_ERROR_NONE)
            cpl_msg_error(cpl_func, "Detected at %s:%s:%d\n", __FILE__,
                          __FUNCTION__, __LINE__);
        
        
        /*To avoid the objects/pixels which are edge of the window */
        if (im_size_x - px1 < 3 || px1 < 3 || im_size_y - py1 < 3 || py1 < 3) {
            for (j = 1; j < im_size_x; j++) {
                error = cpl_image_set(FI_image, 1, j, 0.0);
                if (error != CPL_ERROR_NONE)
                    cpl_msg_error(cpl_func, "Detected at %s:%s:%d\n", __FILE__,
                                  __FUNCTION__, __LINE__);
                error = cpl_image_set(FI_image, 2, j, 0.0);
                if (error != CPL_ERROR_NONE)
                    cpl_msg_error(cpl_func, "Detected at %s:%s:%d\n", __FILE__,
                                  __FUNCTION__, __LINE__);
                error = cpl_image_set(FI_image, 3, j, 0.0);
                if (error != CPL_ERROR_NONE)
                    cpl_msg_error(cpl_func, "Detected at %s:%s:%d\n", __FILE__,
                                  __FUNCTION__, __LINE__);
                error = cpl_image_set(FI_image, j, 1, 0.0);
                if (error != CPL_ERROR_NONE)
                    cpl_msg_error(cpl_func, "Detected at %s:%s:%d\n", __FILE__,
                                  __FUNCTION__, __LINE__);
                error = cpl_image_set(FI_image, j, 2, 0.0);
                if (error != CPL_ERROR_NONE)
                    cpl_msg_error(cpl_func, "Detected at %s:%s:%d\n", __FILE__,
                                  __FUNCTION__, __LINE__);
                error = cpl_image_set(FI_image, j, 3, 0.0);
                if (error != CPL_ERROR_NONE)
                    cpl_msg_error(cpl_func, "Detected at %s:%s:%d\n", __FILE__,
                                  __FUNCTION__, __LINE__);
                error = cpl_image_set(FI_image, j, im_size_y - 1, 0.0);
                if (error != CPL_ERROR_NONE)
                    cpl_msg_error(cpl_func, "Detected at %s:%s:%d\n", __FILE__,
                                  __FUNCTION__, __LINE__);
                error = cpl_image_set(FI_image, j, im_size_y - 2, 0.0);
                if (error != CPL_ERROR_NONE)
                    cpl_msg_error(cpl_func, "Detected at %s:%s:%d\n", __FILE__,
                                  __FUNCTION__, __LINE__);
                error = cpl_image_set(FI_image, j, im_size_y - 3, 0.0);
                if (error != CPL_ERROR_NONE)
                    cpl_msg_error(cpl_func, "Detected at %s:%s:%d\n", __FILE__,
                                  __FUNCTION__, __LINE__);
                error = cpl_image_set(FI_image, im_size_x - 1, j, 0.0);
                if (error != CPL_ERROR_NONE)
                    cpl_msg_error(cpl_func, "Detected at %s:%s:%d\n", __FILE__,
                                  __FUNCTION__, __LINE__);
                error = cpl_image_set(FI_image, im_size_x - 2, j, 0.0);
                if (error != CPL_ERROR_NONE)
                    cpl_msg_error(cpl_func, "Detected at %s:%s:%d\n", __FILE__,
                                  __FUNCTION__, __LINE__);
                error = cpl_image_set(FI_image, im_size_x - 3, j, 0.0);
                if (error != CPL_ERROR_NONE)
						cpl_msg_error(cpl_func, "Detected at %s:%s:%d\n", __FILE__,
                                      __FUNCTION__, __LINE__);
            }
            error = cpl_image_get_maxpos(FI_image, &px1, &py1);
            if (error != CPL_ERROR_NONE)
                cpl_msg_error(cpl_func, "Detected at %s:%s:%d\n", __FILE__,
							__FUNCTION__, __LINE__);
        }
        
        error = cpl_image_get_fwhm(FI_image, px1, py1, &fwhm_xx, &fwhm_yy);
        if (error != CPL_ERROR_NONE)
            cpl_msg_error(cpl_func, "Detected at %s:%s:%d\n", __FILE__,
                          __FUNCTION__, __LINE__);
        
        if (fwhm_xx <= FIfwhm_minthreshold || fwhm_yy <= FIfwhm_minthreshold) {
            for (j = 1; j < 2; j++) {
                for (k = 1; k < 2; k++) {
                    error = cpl_image_set(FI_image, px1 - 1 + j, py1 - 1 + k, 0.0);
                    if (error != CPL_ERROR_NONE)
                        cpl_msg_error(cpl_func, "Detected at %s:%s:%d\n", __FILE__,
                                      __FUNCTION__, __LINE__);
                }
            }
            
        }
        
        error = cpl_image_get_maxpos(FI_image, &px1, &py1); /* Get X and Y pixel positions*/
        if (error != CPL_ERROR_NONE)
            cpl_msg_error(cpl_func, "Detected at %s:%s:%d\n", __FILE__,
                          __FUNCTION__, __LINE__);
        error = cpl_image_get_fwhm(FI_image, px1, py1, &fwhm_xx, &fwhm_yy);
        if (error != CPL_ERROR_NONE)
            cpl_msg_error(cpl_func, "Detected at %s:%s:%d\n", __FILE__,
                          __FUNCTION__, __LINE__);
        
        status += 1;
        if (status >= Check_iteration) {
            
            return -2;			       
			
        }
    } while (fwhm_xx <= fwhm_threshold || fwhm_yy <= fwhm_threshold); /*end of do while  */
    
    
    return error;
}

/* 
   Disable the preprocessing mode
*/
void disable() {
	/* change flag for pre-processing */
	enabled = 0;
}

/*
  Which source you are operating for aberration sensor
  Extended source = 1 : uses correlation algorithm to compute slopes
  Point source = 0: Gaussfit is used to compute slopes
 */
void setExtendedSourceForABS(int Extended_source_in) {
    Extended_source = Extended_source_in;
    if(PRINTMSG) {cpl_msg_info(cpl_func,"Extended source status  %d \n",  Extended_source );}
}


/* 
   print error message
 */
void gvacq_error_print(cpl_error_code error, char * message, int line,
                       char * file) {
	if (error != CPL_ERROR_NONE) {
		cpl_msg_error(cpl_func, "Detected %s  at %s:%d \n", message, file, line);
	}

}

/**
 * @internal
 * @brief Returns Master flat image buffer using fits file name
 */
void readMasterFlat(char * MasterFlatFits, cpl_image * MasterFlat) {
	cpl_error_code error = gvoacqLoadFlat(MasterFlatFits, MasterFlat);
	if(error !=0) printf("Error in loading MasterFlat:%s:%s:%d\n", 
			     __FILE__,__FUNCTION__,__LINE__);
}

/**
 * @internal
 * @brief Returns DeadPixel map buffer using fits file name
 */
void readDeadPixel(char * DeadPixelMapFits, cpl_mask * DeadPixelMap) {
	cpl_error_code error = gvoacqLoadDeadPixel(DeadPixelMapFits, DeadPixelMap);
	if(error !=0) printf("Error in loading DeadPixelMap:%s:%s:%d\n", 
			     __FILE__,__FUNCTION__,__LINE__);
    int Msizex =cpl_mask_get_size_x (DeadPixelMap);
	int Msizey=	cpl_mask_get_size_y (DeadPixelMap);
    
	if (Msizex != IMAGE_X || Msizey != IMAGE_Y){
        printf("Error: size of readDeadPixel %d %d \n", Msizex, Msizey); 
	}
}

/* acqcam_map (MasterSkyField=0, MasterSkyPupil=1, MasterSkyABS=2, MasterFlat=3, DeadPixel=4) */
void enableFromImagelist(cpl_imagelist * acqcam_map){

	  if (MasterSkyField != NULL) {
	    cpl_image_delete(MasterSkyField);
	    MasterSkyField = NULL;
	  }
	  if (MasterSkyPupil != NULL) {
	    cpl_image_delete(MasterSkyPupil);
	    MasterSkyPupil = NULL;
	  }
	  if (MasterSkyABS != NULL) {
	    cpl_image_delete(MasterSkyABS);
	    MasterSkyABS = NULL;
	  }
	  if (acqcam_map == NULL || cpl_imagelist_get_size (acqcam_map)< 4){
            printf("Error reading the calibration imagelist. At %s:%s:%d\n", __FILE__, __FUNCTION__, __LINE__);
	  } else {
          MasterSkyField = cpl_imagelist_get(acqcam_map,0);
          int Sky_F_size_x=cpl_image_get_size_x (MasterSkyField);
          int Sky_F_size_y=cpl_image_get_size_y (MasterSkyField);
          if(Sky_F_size_x != IMAGE_X && Sky_F_size_y != IMAGE_Y) 
            printf("Error in reading MasterSkyField array size [%d %d]\n", Sky_F_size_x, Sky_F_size_y );
          MasterSkyPupil = cpl_imagelist_get(acqcam_map,1);
          int Pupil_F_size_x=cpl_image_get_size_x (MasterSkyField);
          int Pupil_F_size_y=cpl_image_get_size_y (MasterSkyField);
          if(Pupil_F_size_x != IMAGE_X && Pupil_F_size_y != IMAGE_Y) 
            printf("Error read MasterSkyFitsPupil array size [%d %d]\n", Pupil_F_size_x, Pupil_F_size_y );
          MasterSkyABS = cpl_imagelist_get(acqcam_map,2);
          int ABS_F_size_x=cpl_image_get_size_x (MasterSkyABS);
          int ABS_F_size_y=cpl_image_get_size_y (MasterSkyABS);
          if(ABS_F_size_x != IMAGE_X && ABS_F_size_y != IMAGE_Y) 
            printf("Error read MasterSkyFitsABS array size [%d %d]\n", ABS_F_size_x, ABS_F_size_y );
	  if (MasterFlat != NULL) {
	    cpl_image_delete(MasterFlat);
	    MasterFlat = NULL;
	  }
          MasterFlat = cpl_imagelist_get(acqcam_map,3);
          int Flat_F_size_x=cpl_image_get_size_x (MasterSkyField);
          int Flat_F_size_y=cpl_image_get_size_y (MasterSkyField);
          if(Flat_F_size_x != IMAGE_X && Flat_F_size_y != IMAGE_Y) 
            printf("Error in reading MasterFlat array size [%d %d]\n",Flat_F_size_x, Flat_F_size_y );
	  if (DeadPixelMap != NULL) {
	    cpl_mask_delete(DeadPixelMap);
	    DeadPixelMap = NULL;
	  }
	  DeadPixelMap = cpl_mask_new(IMAGE_X,IMAGE_Y);
          cpl_image * DeadPixelMapImage =  cpl_imagelist_get(acqcam_map,4);
          cpl_error_code error = cpl_mask_threshold_image(DeadPixelMap, DeadPixelMapImage, 0, 2, CPL_BINARY_1);
	  if(error !=0) printf("Error in loading DeadPixelMap:%s:%s:%d\n", __FILE__,__FUNCTION__,__LINE__);
          cpl_image_delete(DeadPixelMapImage);

          int Msizex =cpl_mask_get_size_x (DeadPixelMap);
          int Msizey= cpl_mask_get_size_y (DeadPixelMap);
          if(Msizex !=IMAGE_X && Msizex !=IMAGE_Y){
            printf("Error in reading DeadPixelMap array size [%d %d] exptected %d %d\n",
                   Msizex, Msizey, IMAGE_X, IMAGE_Y);
	  }
          enabled=1;
	  cpl_msg_info(cpl_func, "gvacqControl  ENABLED ...\n");
        } 
}

/**
 * @internal
 * @brief enabling the image buffers for preprocessing of the detector image,
 * before applying data reduction algorithms.
 * @param MasterFlatFits file and location of Master Flat image data
 * @param MasterSkyFitsField file and location of Master sky image data for field tracker
 * @param MasterSkyFitsPupil file and location of Master sky image data for pupil tracker 
 * @param DeadPixelMapFits file and location of Dead Pixel Map image data
 * called in gvacqMAIN_TASK.C
*/
void enable(    char* MasterFlatFits, char* MasterSkyFitsField, 
		char* MasterSkyFitsPupil, char* MasterSkyFitsABS,
		char* DeadPixelMapFits) {
    
	/* TODO reffits, PT.., AB.. still to be loaded HERE */
	/* load calibration files in memory */
    printf("enable = %d \n", 	enabled);
    
	/* MASTERSKY - Read the Fits file (MasterSkyFits) */
	cpl_msg_info(cpl_func, "MasterSkyField file as  %s\n", MasterSkyFitsField);
	cpl_msg_info(cpl_func, "MasterSkyPupil file as  %s\n", MasterSkyFitsPupil);
	cpl_msg_info(cpl_func, "MasterSkyABS file as  %s\n", MasterSkyFitsABS);

	if (MasterSkyField != NULL) {
	    cpl_image_delete(MasterSkyField);
	    MasterSkyField = NULL;
	}
	if (MasterSkyPupil != NULL) {
	    cpl_image_delete(MasterSkyPupil);
	    MasterSkyPupil = NULL;
	}
	if (MasterSkyABS != NULL) {
	    cpl_image_delete(MasterSkyABS);
	    MasterSkyABS = NULL;
	}
	if(MasterSkyFitsField == NULL){
        printf("Error in reading MasterSkyFitsField file. At %s:%s:%d\n", __FILE__, __FUNCTION__, __LINE__);
	} else {
        MasterSkyField = cpl_image_load(MasterSkyFitsField, CPL_TYPE_DOUBLE, 0, 0);
        int Sky_F_size_x=cpl_image_get_size_x (MasterSkyField);
        int Sky_F_size_y=cpl_image_get_size_y (MasterSkyField);
        if(Sky_F_size_x != IMAGE_X && Sky_F_size_y != IMAGE_Y) 
            printf("Error in reading MasterSkyField array size [%d %d]\n",
                   Sky_F_size_x, Sky_F_size_y );
	}
	if(MasterSkyFitsPupil == NULL){
        printf("Error in reading MasterSkyFitsPupil file. At %s:%s:%d\n", __FILE__, __FUNCTION__, __LINE__);
	} else {
        MasterSkyPupil = cpl_image_load(MasterSkyFitsPupil, CPL_TYPE_DOUBLE, 0, 0);
        int Pupil_F_size_x=cpl_image_get_size_x (MasterSkyField);
        int Pupil_F_size_y=cpl_image_get_size_y (MasterSkyField);
        if(Pupil_F_size_x != IMAGE_X && Pupil_F_size_y != IMAGE_Y) 
            printf("Error in reading MasterSkyFitsPupil array size [%d %d]\n",
                   Pupil_F_size_x, Pupil_F_size_y );
	}
	if(MasterSkyFitsABS == NULL){
        printf("Error in reading MasterSkyFitsABS file. At %s:%s:%d\n", __FILE__, __FUNCTION__, __LINE__);
	} else {
        MasterSkyABS = cpl_image_load(MasterSkyFitsABS, CPL_TYPE_DOUBLE, 0, 0);
        int ABS_F_size_x=cpl_image_get_size_x (MasterSkyABS);
        int ABS_F_size_y=cpl_image_get_size_y (MasterSkyABS);
        if(ABS_F_size_x != IMAGE_X && ABS_F_size_y != IMAGE_Y) 
            printf("Error in reading MasterSkyFitsABS array size [%d %d]\n",
                   ABS_F_size_x, ABS_F_size_y );
	}
	
	/* MASTERFLAT - Fill the MasterFlat image reading the fits file MasterFlatFits */
	cpl_msg_info(cpl_func, "MasterFlat file as  %s\n", MasterFlatFits);
    
	if (MasterFlat != NULL) {
	    cpl_image_delete(MasterFlat);
	    MasterFlat = NULL;
	}
	if(MasterFlatFits == NULL){
        printf("Error in reading MasterFlatFits file. At %s:%s:%d\n", __FILE__, __FUNCTION__, __LINE__);
	}else {
        MasterFlat = cpl_image_load(MasterFlatFits, CPL_TYPE_DOUBLE, 0, 0);
        int Flat_F_size_x=cpl_image_get_size_x (MasterSkyField);
        int Flat_F_size_y=cpl_image_get_size_y (MasterSkyField);
        if(Flat_F_size_x != IMAGE_X && Flat_F_size_y != IMAGE_Y) 
            printf("Error in reading MasterFlat array size [%d %d]\n",
                   Flat_F_size_x, Flat_F_size_y );
	}
	if(MasterFlat == NULL) 
	    printf("Error in reading MasterFlat file. At %s:%s:%d\n", __FILE__, __FUNCTION__, __LINE__); 
    
	/* DEADPIXEL - Make an empty temporary image and extract bitmap- Then fill this Bitmap */
    
	cpl_msg_info(cpl_func, "DeadPixelMap file as  %s\n", DeadPixelMapFits);
	
      
	if (DeadPixelMap != NULL) {
	    cpl_mask_delete(DeadPixelMap);
	    DeadPixelMap = NULL;
	}
	DeadPixelMap = cpl_mask_new(IMAGE_X,IMAGE_Y);
    
	if(DeadPixelMapFits == NULL){
        printf("Error in reading DeadPixelMapFits file. At %s:%s:%d\n", __FILE__, __FUNCTION__, __LINE__);
	}else {
        readDeadPixel(DeadPixelMapFits, DeadPixelMap);
        int Msizex =cpl_mask_get_size_x (DeadPixelMap);
        int Msizey= cpl_mask_get_size_y (DeadPixelMap);
        if(Msizex !=IMAGE_X && Msizex !=IMAGE_Y){
            printf("Error in reading DeadPixelMap array size [%d %d] exptected %d %d\n",
                   Msizex, Msizey, IMAGE_X, IMAGE_Y);
        }
	}
	if(DeadPixelMap == NULL) 
	    printf("Error in reading DeadPixelMap file. At %s:%s:%d\n", __FILE__, __FUNCTION__, __LINE__); 
    
	enabled=1;
	cpl_msg_info(cpl_func, "gvacqControl  ENABLED ...\n");
    
}

/**
 * @internal
 * @brief store sky for background substraction.
*/
void storsky(char *type) {
    if (strcmp(type, "ALL") == 0) {
        store_sky = STORSKY_FIELD | STORSKY_PUPIL | STORSKY_ABS ;
    } else if (strcmp(type, "FIELD") == 0) {
        store_sky = STORSKY_FIELD;
    } else if (strcmp(type, "PUPIL") == 0) {
        store_sky = STORSKY_PUPIL;
    } else if (strcmp(type, "ABS") == 0) {
        store_sky = STORSKY_ABS;
    } else {
        store_sky = 0;
        cpl_msg_error(cpl_func, "gvacqControl: unknows STORSKY type >%s< ...\n", type);
    }
    cpl_msg_info(cpl_func, "gvacqControl  STORSKY %s ...\n", type);
}


/**
 * @internal
 * @brief  Database entry set to golabl variable Check_iteration.
 * @param CheckIteration_in The number of times checked for objects by rejecting hot pixels and cosmic events 
 * Called in gvacqMAIN_TASK.C
*/
void setCheckIteration(int CheckIteration_in) {
    Check_iteration  = CheckIteration_in;
    if(PRINTMSG) {cpl_msg_info(cpl_func,"Check_iteration  %d \n", Check_iteration   );}
}

/**
 * @internal
 * @brief Setting up database entries for measuring the effective wavelengths using H-K color magnitude
 * @param PolyACQ_in Acquisition camera effective wavelngth modelling polynomial 
 * in function of H-K color magintude
 * @param PolyBC_in Beam combiner effective wavelngth modelling polynomial in function of H-K color magintude
 * Called in gvacqMAIN_TASK.C
*/
void setHKpolynomial(double * PolyACQ_in, double *PolyBC_in) {
    int i;
    for(i=0; i<11; i++){
        PolyACQ[i]= PolyACQ_in[i];
        PolyBC[i]= PolyBC_in[i];
    }
}

/**
 * @internal
 * @brief Setting up H-K color magnitude database value to local variable
 * @param HKmagnitude_in H-K color magnitude
 * Called in gvacqMAIN_TASK.C
*/
void setHKmagnitude(double HKmagnitude_in){
    HKmagnitude =HKmagnitude_in; /*a double number */
}

/**
 * @internal
 * @brief Setting up ADR status of calling
 * @param adrCorrState 1 called, 0 off
 * Called in gvacqMAIN_TASK.C
*/
void setAdrCorrState(int adrCorrState) {
    ADRenabled = adrCorrState;
}

/**
 * @internal
 * @brief Setting up observatory weather and geometry parameters to local variables 
 * @param lTemperature temeparture (k)
 * @param lPressure    Pressure millibars
 * @param lHumidity Humidity percentage 
 * @param lAltitude Altitude. Zenith Distance = 90-Altitude 
 * @param lParallacticAngle Angle between Zenith to North axis in sky
 * @param lImageRotationAngles Image rotation angles from VLTI side, which are useful to meaure position angle 
 * (angle between North axis in sky to the detector Y-axi of the Field imager)
 * @param lKmirrorRotation K-mirror totation of GRAVITY in degrees 
 * Called in gvacqMAIN_TASK.C
*/
void setTelescopeEnvironment(double lTemperature,
                             double lPressure,
                             double lHumidity,
                             double lAltitude,
                             double lParallacticAngle,
                             double lFIPixelScale,
                             double *lImageRotationAngles,
			     double *lPupilRotationAngles,
                             double *lKmirrorRotation,
                             double *lkmirrorOffset,
                             double *lrotationOffset)
{
    int i;
    double imgSign     = -1.0;
    double pupSign     = -1.0;
    double frImgSign   = -1.0;
    double frPupSign   = -1.0;
    double kmirSign    = -1.0;
    double rotSign     = -1.0;

    double imgOffset  =180.0;
    double pupOffset  =180.0;
    double NorthAng[4]={0.,0.,0.,0.0};
                   

    /* ZenithAngle = 90-lAltitude; *//* degrees  */
    Temperature=lTemperature; /* Kelvin ?  */
    Pressure=lPressure; /* Millibars ? */
    Humidity=lHumidity; /* percentage */
    
    FIPixelScale=lFIPixelScale;

/*    PosAngleIn = lImageRotationAngles[0]-lKmirrorRotation[0];*/

    
    for(i=0; i<4; i++){
        KmirrorRotation[i]=lKmirrorRotation[i]; /* degrees */
	NorthAng[i]=imgSign * lImageRotationAngles[i] + frImgSign * lKmirrorRotation[i] +  kmirSign * lkmirrorOffset[i] + rotSign * lrotationOffset[i] + imgOffset;
	ZenithAngle[i]=pupSign * lPupilRotationAngles[i] + frPupSign * lKmirrorRotation[i] + kmirSign * lkmirrorOffset[i] + rotSign * lrotationOffset[i] +  pupOffset;
    }
   

    for(i=0; i<4; i++){
    PosAngleIn[i] = NorthAng[i]; /* Posangle is angle between North axis and detector Y-axis */
    ParAngleIn[i] = NorthAng[i]-ZenithAngle[i];/* ParAng is angle between North axis and Zenith axis */
    }
    
    
    ZenithDistance = 90-lAltitude; /* Zenith distance: angle between source to the vertical */  
    if(ZenithDistance>75) {
    /* printf("Zenith distance is larger than 75, now set to 75 degrees\n"); */
    ZenithDistance=75;
    }

}



/**
 * @internal
 * @brief  Set the pupil tracker reference spot (5th spot corrdinates)
 * @param  PT5thSpotcoordinates reference spot (5th spot corrdinates) [4Tel.*4spots*2(x,y)]
 * @see    setABSrefCoordinates, setFIrefCoordinates
 * 
 * Reads the coordinates from PT5thSpotcoordinates and set to the Spot5Refdata variable
 * @par Error Handling:
 * CPL_ERROR_NULL_INPUT if(one of ) the input pointer (s) is NULL
 */
void setPTrefCoordinates(double *PT5thSpotcoordinates) {
	int i;
	for (i = 0; i < 32; i++) {
		Spot5RefData[i] = PT5thSpotcoordinates[i];
	}
}

/*
  The minimum FWHM threshold used to reject hot pixels and cosmic events 
*/
void setFIminFWHM(double FIminFWHM_in) {
    FIfwhm_minthreshold = FIminFWHM_in;
    /*FIfwhm_minthreshold=1.5;*/
    if(PRINTMSG) {cpl_msg_info(cpl_func,"FIfwhm_minthreshold  =%f \n", FIfwhm_minthreshold);}
}

/*
  The minimum FWHM threshold used to reject hot pixels and cosmic events 
*/
void setPTFitStateErrorThreshold(double value_in) {
    PTFitStateErrorThreshold = value_in;
    if(PRINTMSG) {cpl_msg_info(cpl_func,"FitStateErrorThreshold  =%f \n", PTFitStateErrorThreshold);}
}


/*----------------------------------------Database values filling/writing--------------------------------*/
/**
 * @internal
 * @brief  Pupil tracker spots [16 x 2(x,y)] are written to the database.
 * @param  allDataVars database structure
 * @param  PTRefPosition Pupil tracker spots positions [4 Tel. *16 spots *2 (x,y)]
 * @param  PTRefPosition Pupil tracker spots positions error [4 Tel. *16 spots *2 (x,y)]
*/
int fillPT_REF(gvacqALL_DATA *allDataVars, double PTRefPosition[],
               double PTRefPositionError[], double PTSpotsFlux[]) {
	int i;
	for (i = 0; i < 32; i++) {
		((gvacqALL_DATA) *allDataVars).acPtRef->value[i].arm1 = PTRefPosition[i];
		((gvacqALL_DATA) *allDataVars).acPtRef->value[i].arm2 = PTRefPosition[i
                                                                              + 32];
		((gvacqALL_DATA) *allDataVars).acPtRef->value[i].arm3 = PTRefPosition[i + 2
                                                                              * 32];
		((gvacqALL_DATA) *allDataVars).acPtRef->value[i].arm4 = PTRefPosition[i + 3
                                                                              * 32];
        
		((gvacqALL_DATA) *allDataVars).acPtRefErr->value[i].arm1
            = PTRefPositionError[i];
		((gvacqALL_DATA) *allDataVars).acPtRefErr->value[i].arm2
            = PTRefPositionError[i + 32];
		((gvacqALL_DATA) *allDataVars).acPtRefErr->value[i].arm3
            = PTRefPositionError[i + 2 * 32];
		((gvacqALL_DATA) *allDataVars).acPtRefErr->value[i].arm4
            = PTRefPositionError[i + 3 * 32];
	}
	
	for (i = 0; i < 16; i++) {
	((gvacqALL_DATA) *allDataVars).acPtSpotsFlux->value[i].arm1 = PTSpotsFlux[i];
	((gvacqALL_DATA) *allDataVars).acPtSpotsFlux->value[i].arm2 = PTSpotsFlux[i + 16];
	((gvacqALL_DATA) *allDataVars).acPtSpotsFlux->value[i].arm3 = PTSpotsFlux[i + 2* 16];
	((gvacqALL_DATA) *allDataVars).acPtSpotsFlux->value[i].arm4 = PTSpotsFlux[i + 3* 16];
        
	}

	return 0;
}


/**
 * @internal
 * @brief  Pupil tracker measured pupil positions are written to the database
 * @param  allDataVars database structure
 * @param  PupilPosition Pupil lateral and longitudinal pupil shifts  [4 Tel. *  3(x,y,z)]
 * @param  PupilPositionError Pupil lateral and longitudinal pupil shifts errors  [4 Tel. *  (x,y,z)]
*/
int fillPT_POSITION(gvacqALL_DATA *allDataVars, double PupilPosition[],
                    double PupilPositionError[]) {
    
	int i;
	if (PRINTMSG) printf("-------------------------------------------pupil tracker data \n");
    
    
	for (i = 0; i < 3; i++) {
		((gvacqALL_DATA) *allDataVars).acPtPos->value[i].arm1 = PupilPosition[i];
		((gvacqALL_DATA) *allDataVars).acPtPos->value[i].arm2
            = PupilPosition[i + 3];
		((gvacqALL_DATA) *allDataVars).acPtPos->value[i].arm3
            = PupilPosition[i + 6];
		((gvacqALL_DATA) *allDataVars).acPtPos->value[i].arm4
            = PupilPosition[i + 9];
        
		((gvacqALL_DATA) *allDataVars).acPtPosErr->value[i].arm1
            = PupilPositionError[i];
		((gvacqALL_DATA) *allDataVars).acPtPosErr->value[i].arm2
            = PupilPositionError[i + 3];
		((gvacqALL_DATA) *allDataVars).acPtPosErr->value[i].arm3
            = PupilPositionError[i + 6];
		((gvacqALL_DATA) *allDataVars).acPtPosErr->value[i].arm4
            = PupilPositionError[i + 9];
	}
	return 0;
}

/**
 * @internal
 * @brief  Pupil tracker reference spots (5th spot) are written to the database 
 * @param  allDataVars database structure
 * @param  Reference_spots Reference spots of pupil tracker (so called 5th spot) [4 Tel. *4 spots *2 (x,y)]
 * @param  Reference_spots_error Reference spots error of pupil tracker (so called 5th spot) 
 * [4 Tel. *4 spots *2 (x,y)]
*/

int fillREF_SPOT(gvacqALL_DATA *allDataVars, double Reference_spots[],
                 double Reference_spots_error[]) {
    
    cpl_msg_info(cpl_func,"======================================== > fillREF_SPOT <=============="); 
	int i;
	for (i = 0; i < 4; i++) { /*spots*/
        
		((gvacqALL_DATA) *allDataVars).acPtRefSpot->value[i].arm1
            = Reference_spots[i];
		((gvacqALL_DATA) *allDataVars).acPtRefSpot->value[i].arm2
            = Reference_spots[i + 4];
		((gvacqALL_DATA) *allDataVars).acPtRefSpot->value[i].arm3
            = Reference_spots[i + 8];
		((gvacqALL_DATA) *allDataVars).acPtRefSpot->value[i].arm4
            = Reference_spots[i + 12];
        
		((gvacqALL_DATA) *allDataVars).acPtRefSpot->value[i+4].arm1
            = Reference_spots[i + 16];
        ((gvacqALL_DATA) *allDataVars).acPtRefSpot->value[i+4].arm2
            = Reference_spots[i + 4 + 16];
        ((gvacqALL_DATA) *allDataVars).acPtRefSpot->value[i+4].arm3
            = Reference_spots[i + 8 + 16];
        ((gvacqALL_DATA) *allDataVars).acPtRefSpot->value[i+4].arm4
            = Reference_spots[i + 12 + 16];
        
    ((gvacqALL_DATA) *allDataVars).acPtRefErr->value[i].arm1
        = Reference_spots_error[i];
    ((gvacqALL_DATA) *allDataVars).acPtRefErr->value[i].arm2
        = Reference_spots_error[i + 4];
    ((gvacqALL_DATA) *allDataVars).acPtRefErr->value[i].arm3
        = Reference_spots_error[i + 8];
    ((gvacqALL_DATA) *allDataVars).acPtRefErr->value[i].arm4
        = Reference_spots_error[i + 12];
    
    ((gvacqALL_DATA) *allDataVars).acPtRefErr->value[i+4].arm1
        = Reference_spots_error[i + 16];
    ((gvacqALL_DATA) *allDataVars).acPtRefErr->value[i+4].arm2
        = Reference_spots_error[i + 4 + 16];
    ((gvacqALL_DATA) *allDataVars).acPtRefErr->value[i+4].arm3
        = Reference_spots_error[i + 8 + 16];
    ((gvacqALL_DATA) *allDataVars).acPtRefErr->value[i+4].arm4
        = Reference_spots_error[i + 12 + 16];
	}
    
	return 0;
}

/**
 * @internal
 * @brief  Aberration tracker computed Zernike coefficients are written to the database
 * @param  allDataVars database structure
 * @param  ZernikeCoeff Zernike coefficients array [4*28]
 * @param  ZernikeCoeff Zernike coefficients array error [4*28]
*/
int fillZERNIKE(gvacqALL_DATA * allDataVars, double ZernikeCoeff[],
		double ZernikeCoeffError[], int AbsLensletSizeArray[4], double AutoCorrelationDefocus[4]) {
	int i;
	for (i = 0; i < ZernikeNum; i++) {
		((gvacqALL_DATA) *allDataVars).acAbsZernike->value[i].arm1
            = ZernikeCoeff[i];
		((gvacqALL_DATA) *allDataVars).acAbsZernike->value[i].arm2 = ZernikeCoeff[i
                                                                                  + ZernikeNum];
		((gvacqALL_DATA) *allDataVars).acAbsZernike->value[i].arm3 = ZernikeCoeff[i
                                                                                  + 2 * ZernikeNum];
		((gvacqALL_DATA) *allDataVars).acAbsZernike->value[i].arm4 = ZernikeCoeff[i
                                                                                  + 3 * ZernikeNum];
        
		((gvacqALL_DATA) *allDataVars).acAbsZernikeErr->value[i].arm1
            = ZernikeCoeffError[i];
		((gvacqALL_DATA) *allDataVars).acAbsZernikeErr->value[i].arm2
            = ZernikeCoeffError[i + ZernikeNum];
		((gvacqALL_DATA) *allDataVars).acAbsZernikeErr->value[i].arm3
            = ZernikeCoeffError[i + 2 * ZernikeNum];
		((gvacqALL_DATA) *allDataVars).acAbsZernikeErr->value[i].arm4
            = ZernikeCoeffError[i + 3 * ZernikeNum];
	}
	for (i=0; i<4; i++) {
		(*allDataVars->acAbsLensletSizeArray)[i] = AbsLensletSizeArray[i];
		(*allDataVars->acAbsDefocus)[i] = AutoCorrelationDefocus[i];
	}
	
	
	return 0;
}

/**
 * @internal
 * @brief  Field tracker measured values are written to the database entries when ADR not activated
 * @param  allDataVars database structure
 * @param  ObjectPosition Star position 
 * @param  ObjectPosition2 Star2 position 
 * @param  ObjectPositionError Star position error
*/
int fillFI(gvacqALL_DATA *allDataVars, double ObjectPosition[], double ObjectPosition2[],
           double ObjectPositionError[], double FiFWHM[], double FiFlux[]) {
	int i, i4;
	for (i = 0; i < 2; i++) {
		i4 = i * 4;
		((gvacqALL_DATA) *allDataVars).acFiPos->value[i].arm1 = ObjectPosition[i4];
		((gvacqALL_DATA) *allDataVars).acFiPos->value[i].arm2 = ObjectPosition[i4 + 1];
		((gvacqALL_DATA) *allDataVars).acFiPos->value[i].arm3 = ObjectPosition[i4 + 2];
		((gvacqALL_DATA) *allDataVars).acFiPos->value[i].arm4 = ObjectPosition[i4 + 3];
        
		((gvacqALL_DATA) *allDataVars).acFiPos2->value[i].arm1 = ObjectPosition2[i4];
		((gvacqALL_DATA) *allDataVars).acFiPos2->value[i].arm2 = ObjectPosition2[i4 + 1];
		((gvacqALL_DATA) *allDataVars).acFiPos2->value[i].arm3 = ObjectPosition2[i4 + 2];
		((gvacqALL_DATA) *allDataVars).acFiPos2->value[i].arm4 = ObjectPosition2[i4 + 3];

		((gvacqALL_DATA) *allDataVars).acFiPosErr->value[i].arm1
            = ObjectPositionError[i4];
		((gvacqALL_DATA) *allDataVars).acFiPosErr->value[i].arm2
            = ObjectPositionError[i4 + 1];
		((gvacqALL_DATA) *allDataVars).acFiPosErr->value[i].arm3
            = ObjectPositionError[i4 + 2];
		((gvacqALL_DATA) *allDataVars).acFiPosErr->value[i].arm4
            = ObjectPositionError[i4 + 3];
	}
	for (i = 0; i < 4; i++) {
	        (*allDataVars->acFiFWHM)[i] = FiFWHM[i];
	        (*allDataVars->acFiFlux)[i] = FiFlux[i];
        }
    
	return 0;
}

/**
 * @internal
 * @brief  Field tracker measured values are written to the database entries when ADR activated
 * @param  allDataVars database structure
 * @param  ObjectPosition Star position + ADR shift 
 * @param  ObjectPosition2 Star2 position + ADR shift 
 * @param  ObjectPositionError Star position + ADR shift corresponding error
 * @param  effLambda Effective lambda [ACQ, Beam combiner]
 * @param  ADRangle ADR angle in arc-sec
 * @param  ADRxyCorr Shifts in pixels with ADR effects
 * @param  HKmagnitude H-K color magnitude (input only)
 */
int fillFIobj (gvacqALL_DATA *allDataVars, double ObjectPosition[],  double ObjectPosition2[], double ObjectPositionError[],
               double *effLambda, double *ADRangle, double *ADRxyCorr, double HKmagnitude) {
	int i, i4;
	for (i = 0; i < 2; i++) {
		i4 = i * 4;
		((gvacqALL_DATA) *allDataVars).acFiObj->value[i].arm1 = ObjectPosition[i4];
		((gvacqALL_DATA) *allDataVars).acFiObj->value[i].arm2 = ObjectPosition[i4 + 1];
		((gvacqALL_DATA) *allDataVars).acFiObj->value[i].arm3 = ObjectPosition[i4 + 2];
		((gvacqALL_DATA) *allDataVars).acFiObj->value[i].arm4 = ObjectPosition[i4 + 3];

		((gvacqALL_DATA) *allDataVars).acFiObj2->value[i].arm1 = ObjectPosition2[i4];
		((gvacqALL_DATA) *allDataVars).acFiObj2->value[i].arm2 = ObjectPosition2[i4 + 1];
		((gvacqALL_DATA) *allDataVars).acFiObj2->value[i].arm3 = ObjectPosition2[i4 + 2];
		((gvacqALL_DATA) *allDataVars).acFiObj2->value[i].arm4 = ObjectPosition2[i4 + 3];
        
		((gvacqALL_DATA) *allDataVars).acFiObjErr->value[i].arm1
				= ObjectPositionError[i4];
		((gvacqALL_DATA) *allDataVars).acFiObjErr->value[i].arm2
				= ObjectPositionError[i4 + 1];
		((gvacqALL_DATA) *allDataVars).acFiObjErr->value[i].arm3
            = ObjectPositionError[i4 + 2];
		((gvacqALL_DATA) *allDataVars).acFiObjErr->value[i].arm4
            = ObjectPositionError[i4 + 3];
	}
	*allDataVars->acFiAdrLambdaACQ = effLambda[0];
	*allDataVars->acFiAdrLambdaBC = effLambda[1];
	*allDataVars->acFiAdrAngle = *ADRangle;
	*allDataVars->acFiAdrHK = HKmagnitude;
	(*allDataVars->acFiAdrCorr)[0] = ADRxyCorr[0];
	(*allDataVars->acFiAdrCorr)[1] = ADRxyCorr[1];
	if (PRINTMSG) printf("fillFIobj: %f %f \n", ADRxyCorr[0], ADRxyCorr[1]);

	return 0;
}

/**
 * @internal
 * @brief  AberrationSensorCalib measured spots reference positions writing to the database entry.
 * @param  allDataVars database structure
 * @param  ABSRefPosition Shack-Hartmann spots positions 
 * @param  ABSRefPositionError Shack-Hartmann spots positions error
 */
int fillABS_REF(gvacqALL_DATA *allDataVars, double ABSRefPosition[],
		double ABSRefPositionError[], double RefDefocusPositionsMeasured[16]) {
	int i;
	for (i = 0; i < ABSNUMSPOTS*2; i++) {
		((gvacqALL_DATA) *allDataVars).acAbsRef->value[i].arm1 = ABSRefPosition[i];
		((gvacqALL_DATA) *allDataVars).acAbsRef->value[i].arm2 = ABSRefPosition[i + 1 * ABSNUMSPOTS*2];
		((gvacqALL_DATA) *allDataVars).acAbsRef->value[i].arm3 = ABSRefPosition[i + 2 * ABSNUMSPOTS*2];
		((gvacqALL_DATA) *allDataVars).acAbsRef->value[i].arm4 = ABSRefPosition[i + 3 * ABSNUMSPOTS*2];
        
		((gvacqALL_DATA) *allDataVars).acAbsRefErr->value[i].arm1
            = ABSRefPositionError[i];
		((gvacqALL_DATA) *allDataVars).acAbsRefErr->value[i].arm2
            = ABSRefPositionError[i + 1 * ABSNUMSPOTS*2];
		((gvacqALL_DATA) *allDataVars).acAbsRefErr->value[i].arm3
				= ABSRefPositionError[i + 2 * ABSNUMSPOTS*2];
		((gvacqALL_DATA) *allDataVars).acAbsRefErr->value[i].arm4
            = ABSRefPositionError[i + 3 * ABSNUMSPOTS*2];
	}

	for (i = 0; i < 4; i++) {
	((gvacqALL_DATA) *allDataVars).acAbsRefDefocus->value[i].arm1
            = RefDefocusPositionsMeasured[i];
	
	((gvacqALL_DATA) *allDataVars).acAbsRefDefocus->value[i].arm2
            = RefDefocusPositionsMeasured[i+4];
	
	((gvacqALL_DATA) *allDataVars).acAbsRefDefocus->value[i].arm3
            = RefDefocusPositionsMeasured[i+4*2];
	
	((gvacqALL_DATA) *allDataVars).acAbsRefDefocus->value[i].arm4
            = RefDefocusPositionsMeasured[i+4*3];
	}
	return 0;
}


/**
 * @internal
 * @brief  AberrationSensorCalib measured spots reference positions writing to the database entry.
 * @param  allDataVars database structure
 * @param  ABSRefPosition Shack-Hartmann spots positions 
 */
int fillABS_TARGET_POS_AND_FLUX(gvacqALL_DATA *allDataVars, double ABSTargetPosition[], 
				double ABSSlopes[], double ABSSpotsFlux[]) {
	int i;
	for (i = 0; i < ABSNUMSPOTS*2; i++) {
	((gvacqALL_DATA) *allDataVars).acAbsTargetPos->value[i].arm1 = ABSTargetPosition[i];
	((gvacqALL_DATA) *allDataVars).acAbsTargetPos->value[i].arm2 = ABSTargetPosition[i + 1 * ABSNUMSPOTS*2];
	((gvacqALL_DATA) *allDataVars).acAbsTargetPos->value[i].arm3 = ABSTargetPosition[i + 2 * ABSNUMSPOTS*2];
	((gvacqALL_DATA) *allDataVars).acAbsTargetPos->value[i].arm4 = ABSTargetPosition[i + 3 * ABSNUMSPOTS*2];
	}
	
	for (i = 0; i < ABSNUMSPOTS; i++) {
	((gvacqALL_DATA) *allDataVars).acAbsTargetSpotsFlux->value[i].arm1 = ABSSpotsFlux[i];
	((gvacqALL_DATA) *allDataVars).acAbsTargetSpotsFlux->value[i].arm2 = ABSSpotsFlux[i + 1 * ABSNUMSPOTS];
	((gvacqALL_DATA) *allDataVars).acAbsTargetSpotsFlux->value[i].arm3 = ABSSpotsFlux[i + 2 * ABSNUMSPOTS];
	((gvacqALL_DATA) *allDataVars).acAbsTargetSpotsFlux->value[i].arm4 = ABSSpotsFlux[i + 3 * ABSNUMSPOTS];
	}
	
	for (i = 0; i < ABS9x9Slopes; i++) {
	((gvacqALL_DATA) *allDataVars).acAbsSlopes->value[i].arm1 = ABSSlopes[i];
	((gvacqALL_DATA) *allDataVars).acAbsSlopes->value[i].arm2 = ABSSlopes[i + 1 * ABSNUMSPOTS*2];
	((gvacqALL_DATA) *allDataVars).acAbsSlopes->value[i].arm3 = ABSSlopes[i + 2 * ABSNUMSPOTS*2];
	((gvacqALL_DATA) *allDataVars).acAbsSlopes->value[i].arm4 = ABSSlopes[i + 3 * ABSNUMSPOTS*2];
	}
	return 0;
}




/**
 * @internal
 * @brief  Puptiltracker bary center positions writing to the database entry.
 * @param  allDataVars database structure
 * @param  bary center positions
 * @param  bary center positions error
 */
int fillPT_BARYCENTER (gvacqALL_DATA *allDataVars, double baryCenter[], double baryCenterErr[]) {
	int i;
	for (i = 0; i < 4; i++) { /*spots*/

		((gvacqALL_DATA) *allDataVars).acPtBaryCenter->value[i].arm1
            = baryCenter[i];
		((gvacqALL_DATA) *allDataVars).acPtBaryCenter->value[i].arm2
            = baryCenter[i + 4];
		((gvacqALL_DATA) *allDataVars).acPtBaryCenter->value[i].arm3
            = baryCenter[i + 8];
		((gvacqALL_DATA) *allDataVars).acPtBaryCenter->value[i].arm4
            = baryCenter[i + 12];

		((gvacqALL_DATA) *allDataVars).acPtBaryCenter->value[i+4].arm1
            = baryCenter[i + 16];
        ((gvacqALL_DATA) *allDataVars).acPtBaryCenter->value[i+4].arm2
            = baryCenter[i + 4 + 16];
        ((gvacqALL_DATA) *allDataVars).acPtBaryCenter->value[i+4].arm3
            = baryCenter[i + 8 + 16];
        ((gvacqALL_DATA) *allDataVars).acPtBaryCenter->value[i+4].arm4
            = baryCenter[i + 12 + 16];
	}
	return 0;
}

/* ----------------------------------------------------------------------*/
/**
 * @internal
 * @brief  Set the aberration sensor reference spot (5th spot corrdinates)
 * @param  ABSwindowPos  aberration sensor image positions for first telescope [llx, lly, urx, ury]
 * @param  ABScoordinates  the spot window start locations [Xstart(77), Ystart(77)]
 * @see    setPTrefCoordinates, setFIrefCoordinates
 *
 * Reads the coordinates from ABSwindowPos and ABScoordinates and set to the
 * ABSwindow, ABSiStart, ABSjStart variables respectively
 * 
 * @par Error Handling:
 */
void setABSrefCoordinates(int * ABSwindowPos, double *ABScoordinates, 
			  int ABSwindowSize_in, char  AB_inv_Z2S_68_136Fits_in[], 
			  double *ABSDefocusFactroryCorrd) {
	int i;
	
	if(PRINTMSG) {
        cpl_msg_info(cpl_func,"Set ABS window Pos X (T1-T4):  %d %d %d %d\n", 
                     ABSwindowPos[0], ABSwindowPos[1], ABSwindowPos[2], ABSwindowPos[3]);	
        cpl_msg_info(cpl_func,"Set ABS window Pos Y (T1-T4):  %d %d %d  %d  \n", 
                     ABSwindowPos[4], ABSwindowPos[5], ABSwindowPos[6], ABSwindowPos[7]);	
        cpl_msg_info(cpl_func,"ABSwindowSize_in = %d \n", ABSwindowSize_in);
        
	}
	for (i = 0; i < ABSNUMSPOTS*2*4; i++) {
        ABSiStart[i] = ABScoordinates[i];
        ABSRefSpotsCenters[i]=ABScoordinates[i]; /* Hack */
    }
	
	for (i = 0; i < 8; i++)
		ABSwindow[i] = ABSwindowPos[i];
	ABSwindowSize = ABSwindowSize_in;
	if(AB_inv_Z2S_68_136Fits_in == NULL) {
        printf("Error in reading the AB_inv_Z2S_68_136Fits matrix \n");
	}else{
        if(PRINTMSG) cpl_msg_info(cpl_func,"ABS matrix file read as %s \n",  
				  AB_inv_Z2S_68_136Fits_in);
	}
	for(i=0; i<255; i++) AB_inv_Z2S_68_136Fits[i] = AB_inv_Z2S_68_136Fits_in[i];
	for(i=0; i<16; i++) ABSDefocusFactory[i]= ABSDefocusFactroryCorrd[i];
	for(i=0; i<16; i++) printf(" %f", ABSDefocusFactory[i]);
}

/**
 * @internal
 * @brief The AberrationSensor function uses the latest AberrationSensorCalib computed reference positions.
 *    For this, the reference positions are copied to local variable ABSRefSpotsCenters and passed to 
 *   AberrationSensor function.
 * @param dbCopy database array
 */
void copyAbsRefPosition(gvacqALL_DATA *dbCopy) {
	int i=0;	

	if (((gvacqALL_DATA) *dbCopy).acAbsRef->value[0].arm1 != 0.0) { /* no calibration done yet */ 
	    for (i = 0; i < ABSNUMSPOTS*2; i++) {
            ABSRefSpotsCenters[i]= ((gvacqALL_DATA) *dbCopy).acAbsRef->value[i].arm1;
            ABSRefSpotsCenters[i+ABSNUMSPOTS*2]=((gvacqALL_DATA) *dbCopy).acAbsRef->value[i].arm2;
            ABSRefSpotsCenters[i+ABSNUMSPOTS*2*2]=((gvacqALL_DATA) *dbCopy).acAbsRef->value[i].arm3;
            ABSRefSpotsCenters[i+ABSNUMSPOTS*2*3]=((gvacqALL_DATA) *dbCopy).acAbsRef->value[i].arm4;
	    }
	    
	    for (i = 0; i < 4; i++) {
	    RefDefocusPositionsMeasured[i]=((gvacqALL_DATA) *dbCopy).acAbsRefDefocus->value[i].arm1;
	    RefDefocusPositionsMeasured[i+4]=((gvacqALL_DATA) *dbCopy).acAbsRefDefocus->value[i].arm2;
	    RefDefocusPositionsMeasured[i+4*2]=((gvacqALL_DATA) *dbCopy).acAbsRefDefocus->value[i].arm3;
	    RefDefocusPositionsMeasured[i+4*3]=((gvacqALL_DATA) *dbCopy).acAbsRefDefocus->value[i].arm4;
	    }
	}

}

void setABSLensletSize(int ABSLensletSizeIN){
    ABSLensletSize = ABSLensletSizeIN;
    printf("LensletSize chosen = %d \n",  ABSLensletSize);
}

/*
 * Trigger , that the next incoming/ pre-processed image will be saved
 */
void triggerSaveImage(char* fileinfo, short saveField, short savePupil, short saveAberr){

   save_image_field = saveField;
   save_image_pupil = savePupil;
   save_image_aberr = saveAberr;

   sprintf(file_info, "%s", fileinfo);
}

/*
 * Actual saving an Image
 * TODO Error handling
 */
void SaveImage(cpl_image *DetPointer, char* type){
  char filename[256];
  ccsTIMEVAL currentTime;
  timsMODE mode;
  vltSTRING80 times;

  ccsERROR myError;

  /* Make filename with time */
  timsGetUTC(&currentTime, &mode, &myError);
  /* memcpy(times, NULL, 80); */
  timsTimeToIsoString(&currentTime, timsABS_TIME, 6, times, &myError);

  sprintf(filename,"%s/SYSTEM/DETDATA/GRAVI_%s.%s.fits", getenv("INS_ROOT"),type,times);
  /* sprintf(filename,"/data/GRAVITY/INS_ROOT/SYSTEM/DETDATA/GRAVI_%s.%s.fits", type,times); */

  /* Fill FITS Header */
  cpl_propertylist *pl = cpl_propertylist_new();
  cpl_propertylist_append_string(pl, "FILEINFO", file_info);
  cpl_propertylist_append_string(pl, "ESO DPR CATG", "TEST");
  cpl_propertylist_append_string(pl, "ESO INS ID", "GRAVITY/176642T");
  double mjd = timsUTCToMJD(&currentTime);
  cpl_propertylist_append_double(pl, "MJD-OBS", mjd); 

  /* Save image with Fits header entries */
  cpl_image_save(DetPointer, filename, CPL_BPP_IEEE_DOUBLE, pl, CPL_IO_DEFAULT);
  cpl_propertylist_delete(pl);

  /* Put filename to database to trigger volac */
  /*  dbWrite(filename, "@wgv:Appl_data:GRAVITY:VOLAC", "newDataFileName"); */
  /* rdembet - 2016/03/19 : OLDB writing with dbWriteSymbolic */
  	dbATTRTYPE    attrType;
	dbTYPE        dataType[2];
  	vltINT32      actual;
  	vltUINT16     recCnt;
  	
  	dataType[0] = dbBYTES256;
  	dataType[1] = 0;

  	recCnt   = 1;
  	attrType = dbSCALAR;

  	if (dbWriteSymbolic("@wgv:Appl_data:GRAVITY:VOLAC","newDataFileName",
			    dataType,filename,strlen(filename),
			    &actual,recCnt,attrType,&myError) == FAILURE) {
      	   printf("ERROR DURING FILENAME OLDB WRITING !!\n");
      	}

/* command: msgSend "" gvacqControl SAVEIMG "TEST"*/
}

/*___oOo___*/
