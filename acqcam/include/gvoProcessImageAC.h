/*********************************************************************
 * E.S.O. - VLT project
 *
 * "@(#) $Id: gvoProcessImageAC.h 293322 2017-02-06 14:22:13Z pgarcia $"
 *
 * Estimates the wavefront and telescope parametres for guiding 
 *
 * who       when        what
 * --------  ----------  ----------------------------------------------
 * narsi     2016-05-06  Reading extra VLTI pupilRotIP database point using setTelescopeEnvironment
 * ekw       2016-04-28  Save image separate depending on mode
 * narsi     2016-03-27  edited fillABS_REF to write ABSSlopes to database 
 * narsi     2016-03-22  edited fillPT_REF to write PTSpotsFlux to database 
 * narsi     2016-03-21  edited fillABS_REF to write ABSSpotsFlux to database 
 * narsi     2016-03-19  edited ABS function definations to suit database filling  
 * ekw       2016-03-11  added triggerSaveImage, Save Image
 *                             Modify for Narsi gvoacqAberrationSensor
 * feisenha  2016-02-25  added STORSKY_ABS
 * jott      2014-01-27  added function gvacqProcessShmImageAC()
 * narsi     2013-05-27  defins added
 * narsi     2012-10-27  created
 */

#ifndef PROCESS_IMAGE_AC_H
#define PROCESS_IMAGE_AC_H 
//#define PACKAGE_BUGREPORT "narsireddy.anugu@fe.up.pt"

/************************************************************************
 * gviProcessImageAC.h - brief description
 *----------------------------------------------------------------------
 */

#define _POSIX_SOURCE 1

#include <clipm.h>
#include <cpl_test.h>
#include <time.h>
#include <string.h>
#include "gvacqDefines.h"
#include "gvacqTypes.h"
#define gvoacqWRITE2FILE
#define gvoacqWRITE2DB

/* STORESKY type must be powers of 2 */
enum storsky_tpes {
	STORSKY_FIELD = 1,
	STORSKY_PUPIL = 2,
	STORSKY_ABS   = 4
};

struct acqStruct
{
    double *ABS_REF;
    double *ABS_REF_ERROR;
    double *ABS_ZERNIKE;
    double *ABS_ZERNIKE_ERROR;
    double *PT_REF;
    double *PT_REF_ERROR;
    double *PT_POSITION;
    double *PT_POSITION_ERROR;
    double *PT_REF_SPOT;
    double *PT_REF_SPOT_ERROR;
    double *FI;
    double *FI2;
    double *FI_ERROR;
    char   *FILENAME;
    double *ADR_pixels;
    double *ADR_arcsec;
    double *EffectiveLambda;
};

enum acmode_bits {
  ACMODE_NONE = 0,
  ACMODE_PUPIL = 1,
  ACMODE_FIELD = 2,
  ACMODE_ABERR = 4,
  ACMODE_FIBREF = 16,
  ACMODE_ABBCAL = 32,
  ACMODE_BARCAL = 64
};

void usage(void);
void setCoordinates(int* coordinates);
void setPTrefCoordinates(double *PTcoordinates);
void setABSrefCoordinates(int *ABSwindowPos, double *ABScoordinates,
			  int ABSwindowSize, char ABS_inv_Z2S_68_136Fits[], 
double *ABSDefocusFactroryCorrd);
int setABS_REF(gvacqALL_DATA *allDataVars);
void setFIfitWindowSizeAndSigma(double * FIfieldWindowSizeAndSigma);
void setFIStar2fromStar1(int tel, int dx, int dy);
void setRefWindowPosPT(vltINT32 * RefWindowPosPT);	  
void setPTImageWindowSize(vltINT32 PTImageWindowSize); 
void setPTwindowGauss(vltINT32 PTwindowGauss);	  
void setFIwindowGauss(vltINT32 FIwindowGauss);	  
void setPTspotsScanSigma(double * PTspotsScanSigma);
void setFIminFWHM(double FIminFWHM);   
void setPTFitStateErrorThreshold(double value_in);
void setCheckIteration(int CheckIteration_in);
void setExtendedSourceForABS (int ExtendedSource);
void setTelescopeEnvironment(double lTemperature,
			       double lPressure,
			       double lHumidity,
			       double lAltitude,
			       double lParallacticAngle,
			       double lFIPixelScale,
			       double *lImageRotationAngles,
			     double *lPupiltionAngles,
			     double *lKmirrorRotation,
                             double *lkmirrorOffset,
                             double *lrotationOffset);


void storeSky(cpl_image *image, int type);

void setHKmagnitude(double HKmagnitude_in);
void setHKpolynomial(double *PolyACQ_in, double *PolyBC_in);
void setAdrCorrState(int adrCorrState);
void setABSLensletSize(int ABSLensletSizeIN);

cpl_error_code Check_objects_found(cpl_image *FI_image, double fwhm_threshold);
cpl_error_code ADRCorrPixelsXY(double H_K, double *ACQPolyCoeff, double *BCPolyCoeff, 
			       double P, double t, 
			       double humidity, double Zenith, double* PosAng, 
			       double *ParAng, 
			       double FIPixelScale, double * EffectiveLambda, 
			       double *XYoffsetArcSec, double *XYoffsetPX);

void enable(char* MasterFlatFitsI,
	    char* MasterSkyFitsField,
	    char* MasterSkyFitsPupil,
	    char* MasterSkyFitsABS,
	    char* DeadPixelMapFitsI);
void disable();

void gvacq_error_print(cpl_error_code error, char * message, int line,
                        char * file);

cpl_error_code gvoacqPreprocess(cpl_image *, cpl_image *, cpl_image *,  cpl_mask *);

extern int gvoacqProcessImageAC(cpl_image*      DetPointer,
				int invertFlag,
				int acMode,
				gvacqALL_DATA *allDataVars);

extern void storsky(char *type);
extern void triggerSaveImage(char *filename, short saveField, short save_Pupil, 
			     short saveAberr);
void SaveImage(cpl_image*      DetPointer, char * type);

extern int gvacqProcessShmImageAC(cpl_image*      detPointer,
                                  unsigned char   newFileName[256], 
                                  int*            KmirrorRotation, 
                                  double*         ABSRefPositionRead, 
                                  double*         ABSRefPositionErrorRead, 
                                  char            QueryString,
                                  gvacqALL_DATA*  allDataVars);



cpl_error_code gvacqFieldImager(cpl_image *, double, double *, double *, double*,
				double *, double *, int *);

cpl_error_code gvacqPupilTracker(cpl_image * DetPointer,
                                 double * Spot5RefData,
                                 double   FIfwhm_minthreshold,
				 double PTFitStateErrorThreshold,
                                 double * PupilPosition,
                                 double * PupilPositionError, 
                                 double * barycentre,
                                 double * barycentre_error, 
                                 double * PTRefPosition,
                                 double * PTRefPositionError,
				 double * SpotsFlux,
				 int *PTErrorStatus );

cpl_error_code PupilWindowFlux(cpl_image * PupilDetector,
                                       int *Window_coordinates,
                                       int Window_size,
                                       double *FluxONStatus);

cpl_error_code gvacqPupilTrackerCalib(cpl_image *, int *,  double *, double *);

cpl_error_code gvacqPupilTrackerReferenceSpot(cpl_image *,
					       double    *,
					       double    *, 
					       double    *);

cpl_error_code gvacqPupilTrackerFCspots (cpl_image *image, 
					  double    *GuessCentre, /* Xc[0-3], Yc[4-7] */
					  double    *FCspots,  /* Xcm[0-3], Ycm[4-7] */
					 double    *FCspots_error, 
					 int * PTCalibErrorStatus);



cpl_error_code gvoacqAberrationSensor(cpl_image *detector_image, 
				      int * ABSwindow,
				      int ABSwindowSize,
				      double *ABSRefSpotsCenters,
				      char *AB_inv_Z2S_68_136Fits,
				      int Extended_source,
				      double * RefDefocusPositionsMeasured,
				      int *ABSLensletSize,
				      double *ZernikeCoeff,
				      double *ZernikeCoeffError,
				      double * AutoCorrelationDefocus,
				      double *TargetDefocusPositions,
				      int ABSLensletSizeArray[4],
				      double *ABSSpotsPos,
				      double *ABSSlopes,
				      double *ABSSpotsFlux,
				      int *ABSErrorStatus);

cpl_error_code gvoacqAberrationSensorCalib(cpl_image *, int * ABSwindow, int ABSwindowSize, 
                                           double *, int Extended_source, int, double *, 
					   double *, double *, double *);
cpl_error_code splineResample(cpl_image *image, long dimx, double is, double ie, 
               double js, double je, int sampling, double *res);


cpl_error_code gvoacqLoadFlat(char *, cpl_image *);
cpl_error_code gvoacqLoadSky(char *, cpl_image *);
cpl_error_code gvoacqMakeSky(char *, int, cpl_image *);
cpl_error_code gvoacqLoadDeadPixel(char *,  cpl_mask * );


void readMasterSky(char *, cpl_image *);
void readMasterFlat(char *, cpl_image *);
void readDeadPixel(char *, cpl_mask *);


int fillPT_REF      (gvacqALL_DATA *allDataVars, double PTRefPosition[], 
		     double PTRefPositionError[], double PTSpotsFlux[]);
int fillABS_TARGET_POS_AND_FLUX(gvacqALL_DATA *allDataVars, double ABSRefPosition[],  
				double ABSSlopes[], double ABSSpotsFlux[]);

int fillPT_POSITION (gvacqALL_DATA *allDataVars, double PupilPosition[]  , 
		     double PupilPositionError[]);

int fillREF_SPOT    (gvacqALL_DATA *allDataVars, double Reference_spots[], 
		     double reference_Spots_error[]);

int fillZERNIKE     (gvacqALL_DATA *allDataVars, double ZernikeCoeff[]   , 
		     double ZernikeCoeffError[], int *AbsLensletSizeArray, 
		     double AutoCorrelationDefocus[]);

int fillFI          (gvacqALL_DATA *allDataVars, double ObjectPosition[] ,  double ObjectPosition2[] ,
		     double ObjectPositionError[] , double FiFWHM[], double FiFlux[]);

int fillFIobj (gvacqALL_DATA *allDataVars, double ObjectPosition[], double ObjectPosition2[],
	       double ObjectPositionError[],
	       double *effLambda, double *ADRangle, double *ADRxyCorr, double HKmagnitude);
int fillABS_REF     (gvacqALL_DATA *allDataVars, double ABSRefPosition[] , 
		     double ABSRefPositionError[], double DefocusReference[16]);
int fillPT_BARYCENTER (gvacqALL_DATA *allDataVars, double baryCenter[], 
		       double baryCenterErr[]);
double *ReadTextData(char *filename);

unsigned char *getrefPTAlgo();

#ifdef gvoacqWRITE2FILE
cpl_error_code gvoacqWrite2Db(char *, double * );
#endif /* !gvoacqWRITE2FILE */

#ifdef gvoacqWRITE2DB
cpl_error_code gvacqWrite2File(struct acqStruct *, char*);
#endif /* !gvoacqWRITE2DB */

#endif /*!PROCESS_IMAGE_AC_H*/
/*___oOo___*/
