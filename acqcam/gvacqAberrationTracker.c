 /************************************************************************
 *   This file is part of the E.S.O. - VLTI GRAVITY project
 *
 *   
 *   who       when        what
 *   --------  ----------  ----------------------------------------------
 *   amorim   2016-05-15  Fixed auto determination of window size (was crashing
 *   narsi    2016-05-13  Measured independent focus is written to the Defocus database point
 *   narsi    2016-05-12  Measuring defocus only using four slopes to Zernikes matrix (not 28)
 *   narsi    2016-05-08  Documentation and cleaning
 *   narsi    2016-05-08  Debugged focus measurement using auto-correlation
 *   narsi    2016-04-13  special badpixel handling, slopes ouliers correction, SpecialFocus
 *   narsi    2016-03-27  Measuring defocus from 4 spots
 *   narsi    2016-03-22  Writing target Shack-Hartmann spots positions and their flux to the database
 *   narsi    2016-03-20  cpl_image_fft taking large time for image size of 256 px, now doing with 128 px
 *   narsi    2016-03-19  Added defocus Calib and estimation functions 
 *   narsi    2016-03-19  Improved Zernikes to Slopes matrix and modified gvacqAberrationSensorZ2SLoad
                          accordingly
 *   narsi    2016-03-18  Fixed coredump caused when FlagDefocusWithAutoCorrelation=0
 *   narsi    2016-03-13  Defocus estimation with auto-correlation
 *   narsi    2016-02-20  Implemented variable number of spots in WFS
 *   narsi    2015-09-20  Splitted the gvacqProcessImageAC.c into gvacqAberrationTracker.c
 *   narsi    2015-09-20  created
 */

/****************************************************************************
 *   NAME 
 *   gvacqAberrationTracker.c - Measures Zernike coefficients which corresponding to telescope input beams abberations
 *  
 *   SYNOPSIS
 *   #include <gvoProcessImageAC.h>
 * 
 *   DESCRIPTION
 *   It provieds the wavefront aberrations of input telescope beams. 
 *   The aberrations are represented in 28 Zernike coefficients. 
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
 * @defgroup process GRAVITY observation software for acquisition camera
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


/*-----------------------------------------------------------------------------
 Defines
 ----------------------------------------------------------------------------*/
#define VLT2014
#define PI 3.14159265359
#define ZernikeNum 28       /* Number of Zernikes to measure */
#define Nslopes_ 136
#define f_ABSlenslet 4.97e-3  /* aberration sensor focal length in meter */
#define ABSNUMSPOTS 77 /* number of spots avilable in Shack-Hartmann  */
static int RON=12;
#define PRINTMSG 0

/*-----------------------------------------------------------------------------
 Private function prototypes
 ----------------------------------------------------------------------------*/
cpl_error_code gvacqAberrationSensorSpots2Slopes(cpl_image *fimage,
						 int *imistart2, int *imjstart2, 
						 int nsubs, int binxy2, int yoff,
						 int *validsubs, double *mesvec, 
						 double *mesvec_error);

cpl_error_code gvacqAberrationSensorSpots2CPL_CENTROID(cpl_image *fimage,
		int *imistart2, int *imjstart2, int nsubs, int binxy2, int yoff,
		int *validsubs, double *mesvec, double *mesvec_error);

cpl_error_code gvacqAberrationSensorSpotsMultiGauss(cpl_image *fimage,
						    int *imistart2, int *imjstart2, int * validsubs,
						    double *mesvec, double *mesvec_error, double *Flux);
cpl_error_code CorrCen(cpl_image *image, cpl_image * template,
		double InterpSampling, int key, double *XYcen);
double corr_coeff_with_interp(cpl_image *Image, cpl_image *template, int i,
				      int j, int win, double step, int key);

cpl_error_code gvacqAberrationSensorSlopes2Zernike(double *, double *,
						   char *, int, int , double *, double *);

cpl_error_code gvacqcorr_template_image(cpl_image *fimage, int *imistart2,
		int *imjstart2, int nsubs, int *validsubs, cpl_image *template);
cpl_error_code gvacqcorr_template_image_ext(cpl_image *fimage, int *imistart2,
					    int *imjstart2, int nsubs, int *validsubs, cpl_image *template);
cpl_error_code gvacqCorr_center_error(cpl_image *, double, double, double*);
cpl_error_code ParabolaFit(cpl_image *, double *);
cpl_error_code FFTCorrCen(cpl_image *template, cpl_image *image, double *XYcen);
cpl_error_code image_roll_double(double * image, int im_x);
cpl_error_code roll2d(double * a, long ioff, long istd, long n, long n2,
		double * ws);
int gvoacqAberrationSensorSpots2Slopes(double *, int *, int *, int, int, int,
		int, int, double *, double *, int *, double *, double *);

cpl_error_code gvoacqAberrationSensorSlopes2Zernike(double *, double *,
						    const char *,
						    double *, double *);
cpl_error_code gvacqAberrationSensorSlopes2ZernikeFocusOnly(double *slopes,
							  double *slopes_error, 
							  char * AB_inv_Z2S_68_136Fits,
							  int Slopes_size,
							  int ABSlensletSize,
							  double *Zernike_coefficient_array, 
							  double *Zernike_coefficient_array_error);
cpl_error_code FFTCorrCenInterp(cpl_image *, cpl_image *, cpl_image *, double, double *); 
int roundd(double value);

cpl_error_code gvacqAberrationSensorZ2SLoad(int ArraySize,
					    int *llx, int *lly, int *urx, int *ury);

int GetValidSubsWithABSLensletSize(int ABSLensletSize, int *validsubs, int *Nslopes);
cpl_error_code gvacqAberrationSensorSpotsFlux(cpl_image *fimage,
					      int *imistart2, int *imjstart2, int *validsubs,
					       double *Flux);
cpl_error_code gvacqAberrationSensorEstimateABSLensletSize(cpl_image *detector_image, 
							   int * ABSwindow,
							   int  ABSwindowSize,
							   double *ABSRefSpotsCenters,
							   int* ABSLensletSizeArray, 
							   int * ABSErrorStatus);
cpl_error_code FocusWithAutoCorrelation(cpl_image *detector_image, 
					int * ABSwindow,
					int  ABSwindowSize,
					double * RefDefocusPositionsMeasured,
					double* FocusArray,
					double *TargetPositions,
					int * ABSErrorStatus);
cpl_image * FFTCorr(cpl_image *ref, cpl_image *im);
cpl_error_code MultipleGaussFit(cpl_image * CorrPeak, double *locX, double *locY, 
				int Npeaks, int window, double *FitPositions, 
				double *FitPositions_error);
cpl_error_code FocusWithAutoCorrelationCalib(cpl_image *detector_image,
					     int * ABSwindow,
					     int  ABSwindowSize, double* FactoryRefPositions,
					     double* RefPositionsMeasured,
					     double* RefPositionsMeasuredError,
					     int * ABSErrorStatus);
int VerifySlopes(double *copy_double, int arraySize, double * AvgAndRMS);
cpl_error_code SpecialFocus(cpl_image *detector_image, 
					int * ABSwindow,
					int  ABSwindowSize,
					double *RefDefocusPositionsMeasured,
					double* FocusArray,
					int * ABSErrorStatus);
int BadPixDetection(double *ima, int *bpm, double nsigma, int stride, int npix);


/*
 *******************************************************************************
 *                                                                             *
 * Aberration Sensor                                                           *
 *                                                                             *
 *******************************************************************************
 */

/**
 * @internal
 * @brief  Measures Zernike coefficients to the four telescopes aberrated wavefront
 * @param  detector_image   pointer to the  detector image
 * @param  ABSwindow aberation sensor window centers
 * @param  Aberration sensor window size
 * @param  ABSRefSpotsCenters Reference centers
 * @param  AB_inv_Z2S_68_136Fits slopes2zernike matrix
 * @param  Extended_source source of operatiom, if point source = 0, extended = 1
 * @param  ZernikeCoeff (output) returns Zernike coefficients [68] for the 4 telescopes
 * @param  ZernikeCoeffError (output) returns Zernike coefficients [68] errors for the 4 telescopes
 * @return status error
 * @see    gvacqAberrationSensorCalib, gvacqPupilTracker
 * 
 * Functions pre-required:
 * setABSrefCoordinates
 *
 * @par Error Handling:
 * CPL_ERROR_NULL_INPUT if(one of ) the input pointer (s) is NULL
 * CPL_ERROR_NULL
 */
cpl_error_code gvoacqAberrationSensor(cpl_image *detector_image, 
				      int * ABSwindow,
				      int  ABSwindowSize,
				      double * ABSRefSpotsCenters,
				      char * AB_inv_Z2S_68_136Fits,
				      int Extended_source,
				      double *RefDefocusPositionsMeasured,
				      int*ABSLensletSizeIn,
				      double * ZernikeCoeff,
				      double * ZernikeCoeffError,
				      double * AutoCorrelationDefocus,
				      double *TargetDefocusPositions,
				      int ABSLensletSizeArray[4],
				      double  *ABSSpotsPos,
				      double  *ABSSlopes,
				      double *ABSSpotsFlux,
				      int * ABSErrorStatus) {


    /* test entries */
    /* cpl_ensure_code(detector_image || AB_inv_Z2S_68_136Fits, CPL_ERROR_NULL_INPUT); */
    int i, ii, jj;
    if (detector_image == NULL){
        printf("ERROR::-> Error in reading detector buffer to ABS. Error occured at %s:%s:%d\n", 
               __FILE__, __FUNCTION__, __LINE__);
        for (i=0; i<4; i++) ABSErrorStatus[i]=-1;
        return -1;
    }
    
    if (AB_inv_Z2S_68_136Fits == NULL){
        printf("ERROR::-> Error in reading AB_inv_Z2S_68_136Fits matrix. Error occured at %s:%s:%d\n", 
	   __FILE__, __FUNCTION__, __LINE__);
        for (i=0; i<4; i++) ABSErrorStatus[i]=-1;
        return -1;
    }
    
    if(ABSwindow[0] <219-50 || ABSwindow[0] > 219+50 ||
       ABSwindow[1] <680-50 || ABSwindow[1] > 680+50 ||
       ABSwindow[2] <1156-50 || ABSwindow[2] >1156+50 ||
       ABSwindow[3] <1615-50 || ABSwindow[3] > 1615+50 ||
       ABSwindow[4] <978-50 || ABSwindow[4] > 978+50 ||
       ABSwindow[5] <977-50 || ABSwindow[5] > 977+50 ||
       ABSwindow[6] <984-50 || ABSwindow[6] > 984+50 ||
       ABSwindow[6] <976-50 || ABSwindow[6] > 976+50 ){
        for (i=0; i<4; i++) ABSErrorStatus[i]=-1;
        printf("ERROR: The aberration sensor window positions read from database are wrong. \
Please check it!!! \n");
        return -1;
    }
    
    
    
    for(i=0; i<154*4; i++){
        if(ABSRefSpotsCenters[i] <=0 || ABSRefSpotsCenters[i] >= 200) {
            printf("ERROR: The aberration sensor ref spots positions read from database are wrong. \
Please check it!!!! \n");
            for (i=0; i<4; i++) ABSErrorStatus[i]=-1;
            return -1;
        }
    } 
    
    if(ABSwindowSize != 176){
        printf("ERROR: The aberration window size expected is 176, but given from database is %d \n", 
	       ABSwindowSize);
        for (i=0; i<4; i++) ABSErrorStatus[i]=-1;
        return -1;
    }

   
    /*initialise variables */
    cpl_error_code error = CPL_ERROR_NONE;
    int window = 16; /* no of pixels per subaperture */
    int nsub = ABSNUMSPOTS; /* no of spots avilable with 9x9 SH */
   
    /* number of Zernikes we want to output */
    int no_Zernikes_aimed = ZernikeNum;
    double Z_coeff[no_Zernikes_aimed];
    double Z_coeff_error[no_Zernikes_aimed];
    /* double * SHimage; */
    
    int validsubs[nsub];
    int yoffset = 0;
    int  ABSiLoc[ABSNUMSPOTS];
    int  ABSjLoc[ABSNUMSPOTS];
    int vsc;
    int Nslopes;
    int nvalidsubs = 0;
    int ABSLensletSizeInternal;
    int ABSLensletSize;

    for(i=0; i<4; i++) ABSErrorStatus[i] =0;
    for(i=0; i<16; i++) TargetDefocusPositions[i]=0;
    
    /*
 *******************************************************************************
 Measurement 1: defocus estimation from auto-correlation method
 *******************************************************************************
 */

    /*
     error=FocusWithAutoCorrelation(detector_image, 
				    ABSwindow,
				    ABSwindowSize,
				    RefDefocusPositionsMeasured,
				    AutoCorrelationDefocus,
				    TargetDefocusPositions,
				    ABSErrorStatus);   
    */
     /*
*******************************************************************************
Measurement 2: defocus estimation from four spots centroid measurement
*******************************************************************************
*/
     /*
       error=SpecialFocus(detector_image, 
				   ABSwindow,
				   ABSwindowSize,
				   ABSRefSpotsCenters,
				   AutoCorrelationDefocus,
				   ABSErrorStatus);   
     */
     /*
*******************************************************************************
Measurement 3: the usual estimation of Zernike coefficients method 
*******************************************************************************
*/

     for (i = 0; i < no_Zernikes_aimed*4; i++) {
	ZernikeCoeff[i]=0.0;
	ZernikeCoeffError[i]=0.0;
    }
   
    ABSLensletSize = *ABSLensletSizeIn;
    
    if(ABSLensletSize <0) {
    cpl_msg_error(cpl_func, "No good choice, ABSLensletSize should be postive.... %s:%s:%d\n", 
		  __FILE__, __FUNCTION__, __LINE__);
    return -1;
    }
    if(ABSLensletSize>9) ABSLensletSize=9;

    /* If ABSLensletSize =0, estimate automatically the ABSLensletSize */
    if(ABSLensletSize >= 0 && ABSLensletSize < 3){
    	error = gvacqAberrationSensorEstimateABSLensletSize(detector_image,
							    ABSwindow,
							    ABSwindowSize,
							    ABSRefSpotsCenters,
							    ABSLensletSizeArray,
							    ABSErrorStatus);
    	if(ABSErrorStatus[0] !=0 && ABSErrorStatus[1] !=0 && 
	   ABSErrorStatus[2] !=0 && ABSErrorStatus[3] !=0 ){
    		cpl_msg_error(cpl_func, "Aberration Sensor Error: No good sub-apertures %s:%s:%d\n",
    				__FILE__, __FUNCTION__, __LINE__);
    		*ABSLensletSizeIn=ABSLensletSizeArray[0];
    		return -4;
    	}
    } else {
    	int iii;
    	for (iii=0; iii<4 ; iii++) {
    		ABSLensletSizeArray[iii] = ABSLensletSize;
    	}
    }
    
    
    /*
      The no of pixels separation between two telescope pupil images is 472 pixels.
    */
    
    for (ii = 1; ii <= 4; ii++) {/* For each telescope */
        /* printf("%d %d %d %d \n", ABSwindow[ii-1], ABSwindow[ii-1]+ABSwindowSize-1, 
           ABSwindow[ii+3], ABSwindow[ii+3]+ABSwindowSize-1); */
        
    for (i = 0; i < nsub; i++) validsubs[i]=0;
    
    /* If ABSLensletSize =0, estimate automatically the ABSLensletSize */
    ABSLensletSizeInternal=ABSLensletSizeArray[ii-1];
/*
    if(ABSLensletSize >= 0 && ABSLensletSize < 3){
    ABSLensletSizeInternal=ABSLensletSizeArray[ii-1];
    *ABSLensletSizeIn=ABSLensletSizeInternal;
    }else ABSLensletSizeInternal =ABSLensletSize;
*/
   
   

    if(ABSLensletSizeInternal >=3 ){
    for(i=0;i<nsub; i++) validsubs[i]=0;
    error=GetValidSubsWithABSLensletSize(ABSLensletSizeInternal, validsubs, &Nslopes);
   

    int no_of_slopes = Nslopes;
    double centroid_array[Nslopes];
    double centroid_array_FocusOnly[Nslopes];
    double centroid_array_error[Nslopes];
    double Reference_slopes[Nslopes];
    double AvgAndRMS[4];
    int nsigma=5;
    cpl_mask * Mask = NULL;

    nvalidsubs=0;
    for (i = 0; i < nsub; i++) nvalidsubs += validsubs[i];
    double SingleTelSpotsFlux[nvalidsubs];
    /*if(PRINTMSG) printf("ABSLensletSize nvalidsubs Nslopes %d %d %d \n", 
      ABSLensletSizeInternal, nvalidsubs, Nslopes);*/

    /*
     * initialise arrays to measure slopes
     */
   
    for (i = 0; i < no_of_slopes; i++) {
        centroid_array[i] = 0.0;
        centroid_array_error[i] = 0.0;
    }

        cpl_image *window_abs = cpl_image_extract(detector_image, 
                                                  ABSwindow[ii-1], 
                                                  ABSwindow[ii+3], 
                                                  ABSwindow[ii-1]+ABSwindowSize-1, 
                                                  ABSwindow[ii+3]+ABSwindowSize-1);
        

 
	/*********************
	  Hot pixel estimation: START
	***********************/

	double *image_data=cpl_image_get_data_double(window_abs);
	int npix_size_x = (int)cpl_image_get_size_y(window_abs);
	int npix=(int)cpl_image_get_size_x(window_abs)*(int)cpl_image_get_size_y(window_abs);
	
	cpl_image *bpm=cpl_image_new(npix_size_x, npix_size_x, CPL_TYPE_INT);
	int *bpm_data=cpl_image_get_data_int(bpm);

	int nbpm= BadPixDetection(image_data, bpm_data, nsigma, npix_size_x, npix);

	cpl_image *bpm_image= cpl_image_wrap_int (npix_size_x, npix_size_x, bpm_data);
	Mask=cpl_mask_threshold_image_create(bpm_image, 0, 2);

	if(nbpm < 3000) error = cpl_image_reject_from_mask(window_abs, Mask);
	/*error=  cpl_detector_interpolate_rejected(window_abs);*/
         

	cpl_mask_delete(Mask);
	cpl_image_unwrap(bpm_image);
	cpl_image_delete(bpm);
	
      	/*********************
	  Hot pixel estimation: END
	***********************/
	
      
        /* For ABSRefSpotsCenters visit setABS_REF. Reference spots centers/slopes  */
        for(i=0; i<ABSNUMSPOTS; i++) {
            ABSiLoc[i]= roundd(ABSRefSpotsCenters[(ii-1)*ABSNUMSPOTS*2 + i]); /* ABS i spots locations */
            ABSjLoc[i]= roundd(ABSRefSpotsCenters[(ii-1)*ABSNUMSPOTS*2 + ABSNUMSPOTS+ i]); 
/* ABS j spots locations */
            
        }
        
        vsc=0;
        for (i = 0; i <ABSNUMSPOTS; i++) {
            if (validsubs[i]) {
                vsc++;
                Reference_slopes[vsc-1]= ABSRefSpotsCenters[(ii-1)*ABSNUMSPOTS*2 + i];
                Reference_slopes[nvalidsubs+vsc-1]= ABSRefSpotsCenters[(ii-1)*ABSNUMSPOTS*2 + ABSNUMSPOTS+i];
            }
        }
       
        /* cpl_ensure_code(window_abs  , CPL_ERROR_NULL_INPUT); */
        if(window_abs != NULL){
            if(Extended_source)
                {
                    error = gvacqAberrationSensorSpots2Slopes(window_abs, ABSiLoc, ABSjLoc, 
                                                              nsub, window, yoffset, validsubs, 
                                                              centroid_array, centroid_array_error);
                    
                }else
                { /*  point source Gauss fit: time taken compute 0.15 sec*/
                    /* SHimage=cpl_image_get_data_double(window_abs);*/	
                    error = gvacqAberrationSensorSpotsMultiGauss(window_abs,  ABSiLoc, ABSjLoc, 
								 validsubs, centroid_array, 
								 centroid_array_error, SingleTelSpotsFlux);
		   
                }
        }else printf("window_abs is NULL at %s:%s:%d\n", __FILE__, __FUNCTION__, __LINE__);
        
	if(error !=0) ABSErrorStatus[ii-1]+=error;
        cpl_image_delete(window_abs); /* free up the memory */
        
	vsc=0;
        for (i = 0; i <ABSNUMSPOTS; i++) {
            if (validsubs[i]) {
                vsc++;
                ABSSpotsPos[(ii-1)*ABSNUMSPOTS*2 + i]=centroid_array[vsc-1];
		ABSSpotsPos[(ii-1)*ABSNUMSPOTS*2 + ABSNUMSPOTS+i]=centroid_array[nvalidsubs+vsc-1];
		ABSSpotsFlux[(ii-1)*ABSNUMSPOTS + i]=SingleTelSpotsFlux[vsc-1];
            }
        }
	
        for(i=0; i<Nslopes; i++) {
	/*printf("%d %f %f %f %f\n", i, centroid_array[i], Reference_slopes[i], 
	  centroid_array[i]- Reference_slopes[i], centroid_array_error[i]);*/
	centroid_array[i] = (centroid_array[i]- Reference_slopes[i]);
	ABSSlopes[i+(ii-1)*154]=centroid_array[i];
	centroid_array_error[i]=centroid_array_error[i];
        }

	error= VerifySlopes(centroid_array, Nslopes/2, AvgAndRMS);
	double NSigma=15;
	int nn=0;
	for(i=0; i<Nslopes/2; i++){
	if(fabs(centroid_array[i]-AvgAndRMS[0]) > AvgAndRMS[2]*NSigma) 
	    {centroid_array[i]=AvgAndRMS[0]; nn++;}
	if(fabs(centroid_array[i+Nslopes/2]-AvgAndRMS[1]) > AvgAndRMS[3]*NSigma) 
	    {centroid_array[i+Nslopes/2]=AvgAndRMS[1]; nn++;}
	}

        for (i = 0; i < no_Zernikes_aimed; i++) {
	Z_coeff[i] = 0.0;
        Z_coeff_error[i] = 0.0;
	}
	
	double ZOne2Four[4]={0,0,0,0};
        /* Slopes 2 Zernike */
	for(i=0; i<Nslopes; i++) centroid_array_FocusOnly[i]=centroid_array[i];


	/* Measuring full 28 Zernike coefficients */
	error = gvacqAberrationSensorSlopes2Zernike(centroid_array,
                                                    centroid_array_error, 
                                                    AB_inv_Z2S_68_136Fits, /* First array used*/
						    Nslopes,
						    ABSLensletSizeInternal,
                                                    Z_coeff, Z_coeff_error);

	/* Measuring tip-tilt and defocus only */
	error = gvacqAberrationSensorSlopes2ZernikeFocusOnly(centroid_array_FocusOnly,
                                                    centroid_array_error, 
                                                    AB_inv_Z2S_68_136Fits, /* Second array used*/
						    Nslopes,
						    ABSLensletSizeInternal,
                                                    ZOne2Four, Z_coeff_error);

        
	AutoCorrelationDefocus[ii-1] = ZOne2Four[3];

        if (PRINTMSG) printf("Focus-->Only: %0.3f; Total: %0.3f \n",ZOne2Four[3], Z_coeff[3]);

        for (jj = 0; jj < no_Zernikes_aimed; jj++) {
	
            ZernikeCoeff[jj + (ii - 1) * no_Zernikes_aimed] = Z_coeff[jj];
            ZernikeCoeffError[jj + (ii - 1) * no_Zernikes_aimed] = Z_coeff_error[jj];
        }
	

    }else {
    ABSErrorStatus[ii-1]=-4;
    error=-1;
    }
    

    if (error  != 0 ) {
    ABSErrorStatus[ii-1]+=error; 
            if(error==-1) 
		cpl_msg_error(cpl_func, "Aberration Sensor Error: File reading error, at %s:%s:%d\n", 
                                        __FILE__, __FUNCTION__, __LINE__);
            if(error==-2) 
		cpl_msg_error(cpl_func, "Aberration Sensor Error:INCOMPATABLE input, at %s:%s:%d\n", 
                                        __FILE__, __FUNCTION__, __LINE__);
            if(error==-3) 
		cpl_msg_error(cpl_func, "Aberration Sensor Error: NO object  data found, at %s:%s:%d\n", 
                                        __FILE__, __FUNCTION__, __LINE__);
	    if(error==-4) 
		cpl_msg_error(cpl_func, "Aberration Sensor Error: No good sub-apertures %s:%s:%d\n", 
			    __FILE__, __FUNCTION__, __LINE__);
    }/* End of error handling */
        
    }/* for ii */
    
    
    
    return error;
} /* End of gvacqAberrationSensor */



int GetValidSubsWithABSLensletSize(int ABSLensletSize, int *validsubs, int *Nslopes){

    int error=CPL_ERROR_NONE;
    int i;
    int nsub = ABSNUMSPOTS; /* no of spots avilable with 9x9 SH */
    int dummy;
    if(ABSLensletSize<3) ABSLensletSize =3;

    if(ABSLensletSize>=9){/* 9 case */
    *Nslopes=136;    
    
    for (i = 0; i < nsub; i++) {
    if (i == 0 || i == 6 || i == 7 || i == 15 || i==38 || i == 61 || i == 69 || i == 70 || i==76) 
	{
	validsubs[i] = 0;
	}else validsubs[i] = 1;
    }
    
    }else if(ABSLensletSize == 8){/* 8 case */
    int Nzero=25;
    int Sub_Ap_number8[25]={1,2,3,5,6,7,8,9,15,16,17,25,39,53,61,62,63,69,70,71,72,73,75,76,77};
    
    for (i = 0; i < nsub; i++) validsubs[i]=1; 
    
    for (i = 0; i < Nzero; i++){
    dummy=Sub_Ap_number8[i]-1;
    validsubs[dummy] = 0;
    }
    *Nslopes=(nsub-Nzero)*2;

    }else if (ABSLensletSize == 7){/* 7 case */
    int Nzero=33;
    int Sub_Ap_number7[33]={1,2,3,4,5,6,7,8,9,15,16,17,25,26,34,35,39,43,44,52,53,61,62,63,69,70,71,72,73,74,
			    75,76,77}; 
    for (i = 0; i < nsub; i++) validsubs[i]=1;
    
    for (i = 0; i < Nzero; i++){
    dummy=Sub_Ap_number7[i]-1;
    validsubs[dummy] = 0;
    }
    
    *Nslopes=(nsub-Nzero)*2;
    }else if (ABSLensletSize == 6){/* 6 case */
    int Navilable=36;
    int Sub_Ap_number6[36]={11,12,13,19,20,21,22,23,27,28,29,30,31,32,33,36,37,38,40,41,42,45,46,47,48,49,
			     50,51,55,56,57,58,59,65,66,67};

    for (i = 0; i < nsub; i++) validsubs[i]=0;
    
    for (i = 0; i < Navilable; i++){
    dummy=Sub_Ap_number6[i]-1;
    validsubs[dummy] = 1;
    }
    *Nslopes=Navilable*2;

    }else if (ABSLensletSize == 5){/* 5 case */
    int Navilable=24;
    int Sub_Ap_number5[24]={19,20,21,22,23,28,29,30,31,32,37,38,40,41,46,47,48,49,50,55,56,57,58,59};
    
    for (i = 0; i < nsub; i++) validsubs[i]=0;
    
    for (i = 0; i < Navilable; i++){
    dummy=Sub_Ap_number5[i]-1;
    validsubs[dummy] = 1;
    } 
    *Nslopes=Navilable*2;
    }else if(ABSLensletSize == 4){/* 4 case */
    int Navilable=20;
    int Sub_Ap_number4[20]={20,21,22,28,29,30,31,32,37,38,40,41,46,47,48,49,50,56,57,58};
    
    for (i = 0; i < nsub; i++) validsubs[i]=0;
    
    for (i = 0; i < Navilable; i++){
    dummy=Sub_Ap_number4[i]-1;
    validsubs[dummy] = 1;
    }
    *Nslopes=Navilable*2;
    }else if(ABSLensletSize == 3){/* 3 case */
    int Navilable=16;
    int Sub_Ap_number3[16]={20,21,22,29,30,31,37,38,40,41,47,48,49,56,57,58};
    
    for (i = 0; i < nsub; i++) validsubs[i]=0;
    
    for (i = 0; i < Navilable; i++){
    dummy=Sub_Ap_number3[i]-1;
    validsubs[dummy] = 1;
    }
    *Nslopes=Navilable*2;

    }else{
    error=-1;
    }
    return error;
}

/* ----------------------------------------------------------------------*/
/**
 * @internal
 * @brief  Measures reference centroids for the four Shack-Hartmanns using instial guessed positions
 * @param  Ref_image   pointer to the  reference detector image
 * @param  ABSwindow aberation sensor window centers
 * @param  Aberration sensor window size
 * @param  ABSiStart starting positions to find the reference positions
 * @param  Extended_source source of operatiom, if point source = 0, extended = 1
 * @param  ABSRefPosition (output) reference centroids for the 4 telescopes
 * @param  ABSRefPosition (output) reference centroid errors for the 4 telescopes
 * @return error
 * @see    gvacqAberrationSensor, gvacqPupilTrackerCalib
 *
 * Basically, thye Shack-Hartmann should be illuminated with point source wavefront. 
 The centers of the spots are
 * measured by using standerd center of gravity centroid algorithm or multiple Gaussian Fit.  
 *
 * @par Error Handling:
 * CPL_ERROR_NULL_INPUT if(one of ) the input pointer (s) is NULL
 * CPL_ERROR_NULL
 */
cpl_error_code gvoacqAberrationSensorCalib(cpl_image * Ref_image, int * ABSwindow, 
					   int ABSwindowSize, double * ABSiStart,
					   int Extended_source, int ABSLensletSize, 
					   double * ABSRefPosition, double * ABSRefPositionError, 
					   double *FactoryRefPositions, double * RefPositionsMeasured) {
    /* test entries */
    
    if (Ref_image== NULL ){
        printf("Error occured at %s:%s:%d\n",  __FILE__, __FUNCTION__, __LINE__);
        return -1;
    }
    
    if(ABSwindow[0] <226-50 || ABSwindow[0] > 226+50 ||
       ABSwindow[1] <689-50 || ABSwindow[1] > 689+50 ||
       ABSwindow[2] <1160-50 || ABSwindow[2] >1160+50 ||
       ABSwindow[3] <1620-50 || ABSwindow[3] > 1620+50 ||
       ABSwindow[4] <980-50 || ABSwindow[4] > 980+50 ||
       ABSwindow[5] <980-50 || ABSwindow[5] > 980+50 ||
       ABSwindow[6] <980-50 || ABSwindow[6] > 980+50 ||
       ABSwindow[6] <980-50 || ABSwindow[6] > 980+50 ){
    
        printf("ERROR: The aberration sensor window positions read from database are wrong. \
Please check it!!! \n");
        return -1;
    }
    int i, ii;
    
    for(i=0; i<154*4; i++){
        if(ABSiStart[i] <=0 || ABSiStart[i] >= 200) {
            printf("ERROR: The aberration sensor spots factory positions read from database are wrong. \
Please check it!!!! \n");
            return -1;
        }
    } 
    
    if(ABSwindowSize != 176){
        printf("ERROR: The aberration window size expected 176, but given from database is %d \n", 
	       ABSwindowSize);
        return -1;
    }
    
    cpl_error_code error = CPL_ERROR_NONE;
    
    int nsub = ABSNUMSPOTS;
    /* int no_valid_subs = Nslopes_/2; */
   
    int window = 16;
    int yoffset = 0.0;
/*    double *SHimage;*/
    int ABSiLoc[ABSNUMSPOTS];
    int ABSjLoc[ABSNUMSPOTS];
    int no_of_slopes = Nslopes_;
    int nvalidsubs = 0;
    int vsc;
    int validsubs[nsub];  
    double centriod_array_ref[no_of_slopes];
    double centriod_array_ref_error[no_of_slopes];
    double ABSSpotsFlux[nsub];
    /* To choose the 68 valid sub spots  */
    
    for (i = 0; i < nsub; i++) {
    if (i == 0 || i == 6 || i == 7 || i == 15 || i==38 || i == 61 || i == 69 || i == 70 || i==76) 
	{
	validsubs[i] = 0;
	}else validsubs[i] = 1;
    }
  
    for (i = 0; i < nsub; i++) nvalidsubs += validsubs[i];
    
  
    
    for (i = 0; i < no_of_slopes; i++) {
        centriod_array_ref[i] = 0.0;
        centriod_array_ref_error[i] = 0.0;
    }
    int ABSErrorStatus[4];
    
  
    double RefPositionsMeasuredError[16];
    error=FocusWithAutoCorrelationCalib(Ref_image,
					ABSwindow,
					ABSwindowSize, 
					FactoryRefPositions,
					RefPositionsMeasured,
					RefPositionsMeasuredError,
					ABSErrorStatus);

    /*
      The no of pixels separation between two telescope pupil images is 472 pixels.
    */
    
    for (ii = 1; ii <= 4; ii++) {
        cpl_image *window_abs = cpl_image_extract(Ref_image, 
                                                  ABSwindow[ii-1], 
                                                  ABSwindow[ii+3], 
                                                  ABSwindow[ii-1]+ABSwindowSize-1, 
                                                  ABSwindow[ii+3]+ABSwindowSize-1);
        
	/*********************
	  Hot pixel estimation: START
	***********************/
	int nsigma=10;
	cpl_mask * Mask = NULL;
	double *image_data=cpl_image_get_data_double(window_abs);
	int npix_size_x = (int)cpl_image_get_size_y(window_abs);
	int npix=(int)cpl_image_get_size_x(window_abs)*(int)cpl_image_get_size_y(window_abs);
	
	cpl_image *bpm=cpl_image_new(npix_size_x, npix_size_x, CPL_TYPE_INT);
	int *bpm_data=cpl_image_get_data_int(bpm);

	int nbpm= BadPixDetection(image_data, bpm_data, nsigma, npix_size_x, npix);

	cpl_image *bpm_image= cpl_image_wrap_int (npix_size_x, npix_size_x, bpm_data);
	Mask=cpl_mask_threshold_image_create(bpm_image, 0, 2);

	if(nbpm < 3000)	error = cpl_image_reject_from_mask(window_abs, Mask);
	/* error=  cpl_detector_interpolate_rejected(window_abs); */
         

	cpl_mask_delete(Mask);
	cpl_image_unwrap(bpm_image);
	cpl_image_delete(bpm);
	

        /* for ABSiStart visit setABSrefCoordinates. These are the starting points to calibrate */
        for(i=0; i<ABSNUMSPOTS; i++) {
            ABSiLoc[i]= roundd(ABSiStart[(ii-1)*ABSNUMSPOTS*2 + i]); /* ABS i spots locations */
            ABSjLoc[i]= roundd(ABSiStart[(ii-1)*ABSNUMSPOTS*2 + ABSNUMSPOTS+ i]); /* ABS j spots locations */
        }
        
        if(window_abs !=NULL) {
            if(Extended_source)
                {
                    error = gvacqAberrationSensorSpots2Slopes(window_abs, ABSiLoc, ABSjLoc,
                                                              nsub, window, yoffset, validsubs, 
                                                              centriod_array_ref,
                                                              centriod_array_ref_error);
                }else
                { /*  point source Gauss fit*/
/*	SHimage=cpl_image_get_data_double(window_abs);*/
                    error = gvacqAberrationSensorSpotsMultiGauss(window_abs,  ABSiLoc, ABSjLoc, 
								 validsubs, centriod_array_ref, 
                                                                 centriod_array_ref_error, ABSSpotsFlux);
                    
                }
        }else printf("Error at %s:%s:%d\n", __FILE__, __FUNCTION__, __LINE__); /*end of window_abs !=NULL*/
        
        
        
        vsc=0;
        for (i = 0; i <ABSNUMSPOTS; i++) {
            if (validsubs[i]) {
                vsc++;
                ABSiStart[(ii-1)*ABSNUMSPOTS*2 + i] = centriod_array_ref[vsc-1]; 
                ABSiStart[(ii-1)*ABSNUMSPOTS*2 + ABSNUMSPOTS+i]=centriod_array_ref[nvalidsubs+vsc-1]; 
                
                ABSRefPositionError[(ii-1)*ABSNUMSPOTS*2 + i] = centriod_array_ref_error[vsc-1]; 
                ABSRefPositionError[(ii-1)*ABSNUMSPOTS*2 + ABSNUMSPOTS+i]=
		    centriod_array_ref_error[nvalidsubs+vsc-1]; 
            }
        }
        
/* free up the memory */
        cpl_image_delete(window_abs);
    }/*End of ii */
    
    for(i=0; i<ABSNUMSPOTS*2*4; i++)  ABSRefPosition[i]=ABSiStart[i];
    
    return error;
}
/* End of gvacqAberrationSensorCalib */


/* ----------------------------------------------------------------------*/
/**
 * @internal
 * @brief converts Shack-Hartmann array spots into positions using CPL_CENTROID (grid positions)
 * @param fimage  final image with spots
 * @param imistart2 vector of i starts of each image
 * @param imjstart2 vector of j starts of each image
 * @param nsubs  number of subapertures
 * @param binxy2  side size of the single spot image window extracted from the given image
 * @para  fimnx  the aberration sensor image X dimension
 * @param fimnx  the aberration sensor image Y dimension
 * @param yoff  y offset (to process only part of the image)
 * @param centroidw centroid weight vector for centroid computations, X & Y
 * @param validsubs valid subaps within the set of subap for which image is computed
 * @param mesvec final slope measurement vector
 * @param mesvec_error slope measurement error vector
 *
 *
 * pass one aberration sensor image buffer array
 * and a set of indices (start and end of each subapertures). Then it takes the each subaperture
 * spot and computes the weighted centriod and keeps the measurement in the slopes array vector.
 * flow:
 * 1. select the subaperture spot
 * 2. measure the X and Y centriod positions
 * 3. keep the centriod in the mesvec
 * 4. free memory
 *
 * @par Error Handling:
 * CPL_ERROR_NULL_INPUT if(one of ) the input pointer (s) is NULL
 * CPL_ERROR_NULL
 */

cpl_error_code gvacqAberrationSensorSpots2CPL_CENTROID(cpl_image *fimage,
		int *imistart2, int *imjstart2, int nsubs, int binxy2, int yoff,
		int *validsubs, double *mesvec, double *mesvec_error) {

    /* test entries */
/*
  cpl_ensure_code(fimage   || imistart2   || imjstart2  
  || validsubs   , CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (nsubs>0 ||binxy2>0 , CPL_ERROR_ILLEGAL_INPUT);
*/
    
    if (fimage == NULL || imistart2 == NULL || imjstart2 == NULL || validsubs == NULL 
        || nsubs <0 || binxy2 < 0){
        printf("ERROR::-> Inputs for the function are NULL at %s:%s:%d\n", 
               __FILE__, __FUNCTION__, __LINE__);
        
        return -1;
    }
    
    /* Declarations */
    cpl_error_code error = CPL_ERROR_NONE;
    int nvalidsubs, vsc, i, l;
    cpl_image *Live;
   
    vsc = -1;

    for(i=0; i<nsubs; i++){
        imistart2[i]=imistart2[i]-8;
        imjstart2[i]=imjstart2[i]-8;
    }
    
    cpl_image *ffimage2 = cpl_image_duplicate(fimage);
    
    
    nvalidsubs = 0;
    for (i = 0; i < nsubs; i++) nvalidsubs += validsubs[i];
    
    for (l = 0; l < nsubs; l++) {
        vsc++;
        Live = cpl_image_extract(ffimage2, imistart2[l] + 1, imjstart2[l] + 1,
                                 imistart2[l] + 16, imjstart2[l] + 16);
        
        /* cpl_ensure_code(Live  , CPL_ERROR_NULL_INPUT); */
        
        if (Live != NULL ) {
            mesvec[vsc] = imistart2[l] + cpl_image_get_centroid_x(Live) - 0.5 ; 
            mesvec[nsubs + vsc] = imjstart2[l] + cpl_image_get_centroid_y(Live) - 0.5;
        } else printf("Problem with Ref and Live in the aberration sensor at %s:%s:%d \n", 
                      __FILE__, __FUNCTION__, __LINE__);
        
        mesvec_error[vsc] = -1;
        mesvec_error[nsubs + vsc] = -1;
        
        cpl_image_delete(Live);
    }
    
    cpl_image_delete(ffimage2);
    return error;
}


/* ----------------------------------------------------------------------*/
/**
 * @internal
 * @brief Shack-Hartmann spots reference centers using Gaussian fit 
 * @param fimage  final image with spots
 * @param imistart2 vector of i starts of each image
 * @param imjstart2 vector of j starts of each image
 * @param Centroid final reference spot center measurement vector
 * @param Centroid_error reference spot center measurement error vector
 * 
 * The spot centers are measured using the Gaussian fit
 * 
 * @par Error Handling:
 * CPL_ERROR_NULL_INPUT if(one of ) the input pointer (s) is NULL
 * CPL_ERROR_NULL
 */
cpl_error_code gvacqAberrationSensorSpotsMultiGauss(cpl_image *fimage,
						    int *imistart2, int *imjstart2, int *validsubs,
						    double *Centroid, double *Centroid_error, double *Flux) {
    
    if (fimage == NULL || imistart2 == NULL || imjstart2 == NULL ){
        printf("ERROR::-> Inputs for the function are NULL at %s:%s:%d\n", 
               __FILE__, __FUNCTION__, __LINE__);
        
        return -1;
    }
    
    int i;
    int robustness = 9;
    int nvalidsubs = 0;
    cpl_error_code Cerror=CPL_ERROR_NONE;
    int vsc =-1;
    int GaussFitWindowSize=10; /* pixels */
    for (i = 0; i <ABSNUMSPOTS; i++) nvalidsubs += validsubs[i];
    
    cpl_matrix * ABSspotsLoc = cpl_matrix_new(2, nvalidsubs); /* 1 spot (x and y) */
    cpl_matrix * xy_centre = cpl_matrix_new(2, nvalidsubs); /* 1 spot (x and y) */
    cpl_matrix * xy_centre_err = cpl_matrix_new(2, nvalidsubs);
    cpl_matrix * xy_sigma = cpl_matrix_new(2, nvalidsubs);
    cpl_matrix * xy_sigma_err = cpl_matrix_new(2, nvalidsubs);
    cpl_matrix * xy_fwhm = cpl_matrix_new(2, nvalidsubs);
    cpl_matrix * xy_fwhm_err = cpl_matrix_new(2, nvalidsubs);
    cpl_array * all_error_codes = cpl_array_new(1, CPL_TYPE_INT);
    
  
    for(i=0; i<ABSNUMSPOTS; i++){
        if(validsubs[i]) {
            vsc++;
            cpl_matrix_set(ABSspotsLoc, 0, vsc, imistart2[i]); /* Guess spot position to fit Guassin */
            cpl_matrix_set(ABSspotsLoc, 1, vsc, imjstart2[i]); 
        }
    }
    Cerror=gvacqAberrationSensorSpotsFlux(fimage,imistart2, imjstart2, validsubs, Flux);
    Cerror = clipm_centroiding_multi_gauss(fimage, 
                                           ABSspotsLoc, 
                                           GaussFitWindowSize, 
                                           &xy_centre, 
                                           &xy_centre_err, 
                                           &xy_sigma, 
                                           &xy_sigma_err,
                                           &xy_fwhm, 
                                           &xy_fwhm_err, 
                                           NULL, 
                                           &all_error_codes,
                                           robustness);
    
    
    
   
    cpl_matrix *copy=cpl_matrix_duplicate (xy_centre);
    cpl_matrix *Orig=cpl_matrix_duplicate (xy_centre);
    Cerror = cpl_matrix_subtract(copy, ABSspotsLoc);	
    double *copy_double = cpl_matrix_get_data(copy);
    int size_x=(int) cpl_image_get_size_x(fimage);

    double avg_x=0,avg_y=0;
    for(i=0; i<nvalidsubs; i++){
    if(copy_double[i]<size_x) {
    avg_x+=copy_double[i];
    }else  avg_x+=0;
    
    if(copy_double[i+nvalidsubs]<size_x) {
    avg_y+=copy_double[i+nvalidsubs];
    }else  avg_y+=0;
    }
    avg_x=avg_x/((double)nvalidsubs);
    avg_y=avg_y/((double)nvalidsubs);
    

    double *Ref = cpl_matrix_get_data(ABSspotsLoc);
  
    for(i=0; i<=vsc; i++) {
    if(  fabs(copy_double[i]) > GaussFitWindowSize/2.0 ||
	 fabs(copy_double[i+vsc+1]) > GaussFitWindowSize/2.0 ){
	cpl_matrix_set(xy_centre, 0, i, Ref[i]+avg_x);
	cpl_matrix_set(xy_centre, 1, i, Ref[i+vsc+1]+avg_y);
    }
    }
   
    Cerror = cpl_matrix_subtract(Orig, xy_centre);	
    /*cpl_matrix_dump(Orig, NULL);*/

    double *centroid = cpl_matrix_get_data(xy_centre);
    double *centroid_error = cpl_matrix_get_data(xy_centre_err);
    
    for(i=0; i<nvalidsubs*2; i++) {
        Centroid[i] = centroid[i];
        Centroid_error[i]=centroid_error[i];
    }
    
    cpl_matrix_delete(Orig);
    cpl_matrix_delete(copy);
    cpl_matrix_delete(ABSspotsLoc);
    cpl_matrix_delete(xy_centre);
    cpl_matrix_delete(xy_centre_err);
    cpl_matrix_delete(xy_sigma);
    cpl_matrix_delete(xy_sigma_err);
    cpl_matrix_delete(xy_fwhm);
    cpl_matrix_delete(xy_fwhm_err);
    cpl_array_delete(all_error_codes);
    
    return Cerror;
}


int VerifySlopes(double *copy_double, int arraySize, double * AvgAndRMS){
    double RMSerrorX=0;
    double RMSerrorY=0;
    double avg_x=0,avg_y=0;
    int i;

    for(i=0; i<arraySize; i++){
    avg_x+=copy_double[i];
    avg_y+=copy_double[i+arraySize];
    }
    avg_x=avg_x/((double)arraySize);
    avg_y=avg_y/((double)arraySize);

    for(i=0; i<arraySize; i++){
    RMSerrorX+=(avg_x-copy_double[i])*(avg_x-copy_double[i]);
    RMSerrorY+=(avg_y-copy_double[i+arraySize])*(avg_y-copy_double[i+arraySize]);
    }

    RMSerrorX=sqrt(RMSerrorX)/(arraySize);
    RMSerrorY=sqrt(RMSerrorY)/(arraySize);

    AvgAndRMS[0]=avg_x;
    AvgAndRMS[1]=avg_y;
    AvgAndRMS[2]=RMSerrorX;
    AvgAndRMS[3]=RMSerrorY;
    
    return 0;
}

/* ----------------------------------------------------------------------*/
/**
 * @internal
 * @brief converts Shack-Hartmann array spots into slopes (tilts)
 * @param fimage  final image with spots
 * @param imistart2 vector of i starts of each image
 * @param imjstart2 vector of j starts of each image
 * @param nsubs  number of subapertures
 * @param binxy2  side size of the single spot image window extracted from the given image
 * @para  fimnx  the aberration sensor image X dimension
 * @param fimnx  the aberration sensor image Y dimension
 * @param yoff  y offset (to process only part of the image)
 * @param centroidw centroid weight vector for centroid computations, X & Y
 * @param validsubs valid subaps within the set of subap for which image is computed
 * @param mesvec final slope measurement vector
 * @param mesvec_error slope measurement error vector
 *
 *
 * pass one aberration sensor image buffer array
 * and a set of indices (start and end of each subapertures). Then it takes the each subaperture
 * spot and computes the weighted centriod and keeps the measurement in the slopes array vector.
 * flow:
 * 1. select the subaperture spot
 * 2. measure the X and Y centriod positions
 * 3. keep the centriod in the mesvec
 * 4. free memory
 *
 * @par Error Handling:
 * CPL_ERROR_NULL_INPUT if(one of ) the input pointer (s) is NULL
 * CPL_ERROR_NULL
 */

cpl_error_code gvacqAberrationSensorSpots2Slopes(cpl_image *fimage,
		int *imistart2, int *imjstart2, int nsubs, int binxy2, int yoff,
		int *validsubs, double *mesvec, double *mesvec_error) {

	/* test entries */
	/* cpl_ensure_code(fimage   || imistart2   || imjstart2 || validsubs   , CPL_ERROR_NULL_INPUT); */
    
    if (fimage == NULL || imistart2 == NULL || imjstart2 == NULL || validsubs == NULL 
        || nsubs <0 || binxy2 < 0){
        printf("ERROR::-> Inputs for the function are NULL at %s:%s:%d\n", 
               __FILE__, __FUNCTION__, __LINE__);
        
        return -2;
    }
    
    /* Declarations */
    cpl_error_code error = CPL_ERROR_NONE;
    int nvalidsubs, vsc, i, l;
    double xycen[4];
    cpl_image *Live, *ref;
    xycen[0] = 0.0;
    xycen[1] = 0.0;
    vsc = -1;
   
    cpl_image *template = cpl_image_new(binxy2, binxy2, CPL_TYPE_DOUBLE);
    cpl_image *ffimage1 = cpl_image_duplicate(fimage);
    cpl_image *ffimage2 = cpl_image_duplicate(fimage);
    
    for(i=0; i<nsubs; i++){
        /* printf("%d %d \n", imistart2[i], imjstart2[i]);*/
        imistart2[i]=imistart2[i]-8;
        imjstart2[i]=imjstart2[i]-8;
    }
    
    error = gvacqcorr_template_image_ext(ffimage1, imistart2, imjstart2, nsubs,
                                         validsubs, template);
    if(template == NULL) {
        printf("template calculation went wrong at %s:%s:%d \n", __FILE__, __FUNCTION__, __LINE__);
        return -2;
    }
   
    
    double off=0.5;
    double data_t[16*16];
    error = splineResample(template, 0, 1+off, binxy2+off, 
                           1+off, binxy2+off, binxy2, data_t);
    cpl_image * temp_shifted = cpl_image_wrap_double(binxy2, binxy2, data_t);
    
    if (error != CPL_ERROR_NONE)
        cpl_msg_error(cpl_func, "Error occured at %s:%d\n", __FILE__, __LINE__);
    
    /* Count the number of valid sub-apertures */
    nvalidsubs = 0;
    for (i = 0; i < nsubs; i++) nvalidsubs += validsubs[i];
    
    
    /* For each sub aperture, compute the center */
    for (l = 0; l < nsubs; l++) {
        if (validsubs[l]) vsc++;
    
        if (validsubs[l]) {
            ref = cpl_image_duplicate(template);
            Live = cpl_image_extract(ffimage2, imistart2[l] + 1, imjstart2[l] + 1,
                                     imistart2[l] + binxy2, imjstart2[l] + binxy2);
            /*
              cpl_ensure_code(ref, CPL_ERROR_NULL_INPUT);
              cpl_ensure_code(Live  , CPL_ERROR_NULL_INPUT);
            */
            
            if (ref != NULL || Live != NULL )
                {
                    /*  SSD  */
                    /* error = CorrCen(Live,  ref, 1, 1, xycen); */
                    
                    /*  FFT centroiding */
                    error = FFTCorrCenInterp(ref, temp_shifted, Live, off, xycen); 
                    
                    /* CCI */
                    /*error = CorrCen(ref, Live, 1, 1, xycen); */
                    
                } else {
                printf("Problem with Ref and Live in the aberration sensor \n");
                cpl_msg_error(cpl_func, "Detected at %s:%s:%d\n", 
                              __FILE__, __FUNCTION__, __LINE__);
            }
            
            if (error != CPL_ERROR_NONE) cpl_msg_error(cpl_func, "Detected at %s:%s:%d\n", 
                                                       __FILE__, __FUNCTION__, __LINE__);
            
            /* Copying to the return array */
            mesvec[vsc] = xycen[0]+imistart2[l];
            mesvec[nvalidsubs + vsc] = xycen[1]+imjstart2[l];
            
            /* Error bars  */
            mesvec_error[vsc] = xycen[2];
            mesvec_error[nvalidsubs + vsc] = xycen[3];
            
            /* Free up memory  */
            cpl_image_delete(Live);
            cpl_image_delete(ref);
        }
    }/* end of nsubs */
    
    cpl_image_delete(template);
    cpl_image_delete(ffimage1);
    cpl_image_delete(ffimage2);
    cpl_image_unwrap(temp_shifted);
    return error;
}

cpl_error_code gvacqAberrationSensorZ2SLoad(int ArraySize,
					    int* llx, int* lly, 
					    int* urx, int* ury){
    cpl_error_code error=CPL_ERROR_NONE;
    
    int S[7]={136, 104, 88, 72, 48, 40, 32};
    if(ArraySize<3) ArraySize=3;

    *lly=1;
    *ury=28;

    if(ArraySize>=9){
    *llx=1;
    *urx=S[0];
    }else if(ArraySize==8){
    *llx=S[0]+1;
    *urx=S[0]+S[1];
    }else if(ArraySize==7){
    *llx=S[0]+S[1]+1;
    *urx=S[0]+S[1]+S[2];
    }else if(ArraySize==6){
		*lly=1;
		*ury=24;
    *llx=S[0]+S[1]+S[2]+1;
    *urx=S[0]+S[1]+S[2]+S[3];
    }else if(ArraySize==5){
		*lly=1;	
		*ury=16;	
    *llx=S[0]+S[1]+S[2]+S[3]+1;
    *urx=S[0]+S[1]+S[2]+S[3]+S[4];
    }else if(ArraySize==4){
    *lly=1;
		*ury=13;
		*llx=S[0]+S[1]+S[2]+S[3]+S[4]+1;
    *urx=S[0]+S[1]+S[2]+S[3]+S[4]+S[5];
    }else if(ArraySize==3){
    *lly=1;
		*ury=9;
		*llx=S[0]+S[1]+S[2]+S[3]+S[4]+S[5]+1;
    *urx=S[0]+S[1]+S[2]+S[3]+S[4]+S[5]+S[6];
    }else{
    return -1;
    }
  
    return error;
}

/* ----------------------------------------------------------------------*/
/**
 * @internal
 * @brief  Outputs the Zernike coefficient array corresponding input slopes array.
 * @parame slopes slopes vector array corresponding the quasi aberrated wavefront
 * @parame slopes_error slopes error vector array
 * @parame Zernike_coefficient_array (output) Zernike coefficients corresponding the quasi aberrated wavefront
 * @parame Zernike_coefficient_array_error (output) Zernike coefficient error corresponding
 *          the quasi aberrated wavefront
 * @return error
 * @see    gvacqAberrationSensor,gvacqAberrationSensor_spots2slopes
 *
 * It measures exactly 68 Zernike coefficients here (if you want measure more Zernikes
 * need to change the the inverse of Zernike to slopes matrix )
 * @par Error Handling:
 * CPL_ERROR_NULL_INPUT if(one of ) the input pointer (s) is NULL
 * CPL_ERROR_NULL
 * CPL_ERROR_INCOMPATIBLE_INPUT  if AB_inv_Z2S_68_136Fits dims not match
 */

cpl_error_code gvacqAberrationSensorSlopes2Zernike(double *slopes,
                                                   double *slopes_error, 
                                                   char * AB_inv_Z2S_68_136Fits,
						   int Slopes_size,
						   int ABSlensletSize,
                                                   double *Zernike_coefficient_array, 
						   double *Zernike_coefficient_array_error) {
    
    /* test entries */
    cpl_error_code error = CPL_ERROR_NONE;
    int i;
	
    if (AB_inv_Z2S_68_136Fits == NULL){
        printf("ERROR::-> Please make sure database reading AB_inv_Z2S_68_136.fits correctly. \
Error occured at %s:%s:%d\n", 
               __FILE__, __FUNCTION__, __LINE__);
        
        return -1;
    }
    
    if (slopes == NULL || slopes_error == NULL){
        printf("ERROR::-> Input slopes are NULL. Error occured at %s:%s:%d\n", 
               __FILE__, __FUNCTION__, __LINE__);
        return -2;
    }

    int llx, lly, urx, ury;
    cpl_image * inv_of_Zernike_deriv_data = cpl_image_load(AB_inv_Z2S_68_136Fits,
                                                           CPL_TYPE_DOUBLE, 0, 0);


   
    if (inv_of_Zernike_deriv_data == NULL ){
    printf("ERROR:: -> AB_inv_Z2S_68_136.fits INCOMPATABLE input file from database.\n");
    printf("The file path submitted: %s \n", AB_inv_Z2S_68_136Fits);
    
    cpl_image_delete(inv_of_Zernike_deriv_data);
    return -2;
    }
    
    error=gvacqAberrationSensorZ2SLoad(ABSlensletSize,
				       &llx, &lly, &urx, &ury);
    
    cpl_image *buffer=  cpl_image_extract (inv_of_Zernike_deriv_data, llx, lly, urx, ury);
    if (buffer == NULL || error !=0 ){
    printf("ERROR: gvacqAberrationSensorZ2SLoad failed to measurement correct dimension\n");
    if(buffer != NULL) cpl_image_delete(buffer);
    cpl_image_delete(inv_of_Zernike_deriv_data);
    return -3;
    }
   
    int size_y = (int)cpl_image_get_size_y(buffer);
    int size_x = (int)cpl_image_get_size_x(buffer);
    
    double *inv_Z2S_double_Extracted = cpl_image_get_data_double(buffer);
    cpl_matrix *Z2Sinv = cpl_matrix_wrap(size_y, size_x, inv_Z2S_double_Extracted);
    
  
    if (Z2Sinv == NULL){
    printf("ERROR: gvacqAberrationSensorZ2SLoad failed to measurement correct dimension\n");
    
    if(Z2Sinv!=NULL) cpl_matrix_unwrap(Z2Sinv);
    cpl_image_delete(buffer);
    cpl_image_delete(inv_of_Zernike_deriv_data);
    return -3;
    }
  
    /* slopes are copying to a cpl_matrix formalism
       for multiplication requirement
    */
    cpl_matrix *mesvec_matrix = cpl_matrix_new(Slopes_size, 1);
    cpl_matrix *mesvec_error_matrix = cpl_matrix_new(Slopes_size, 1);
	
    /* fill the slopes matrix with zeros */
    for (i = 0; i < Slopes_size; i++) {
        cpl_matrix_set(mesvec_matrix, i, 0, slopes[i]);
        cpl_matrix_set(mesvec_error_matrix, i, 0, slopes_error[i]);
    } 
    
    /* Check whether the inverse Zernike derivative matrix compatible */
    if((int)cpl_matrix_get_ncol (Z2Sinv) == (int)cpl_matrix_get_nrow(mesvec_matrix)){
    cpl_matrix * Zernike_coeff_vec = cpl_matrix_product_create(Z2Sinv,
                                                               mesvec_matrix);
    
    cpl_matrix * Zernike_coeff_error_vec = cpl_matrix_product_create(Z2Sinv,
								     mesvec_error_matrix);
    
    /* Returing the Zernike coefficients to the calling function */
    double *Zresult_tmp = cpl_matrix_get_data(Zernike_coeff_vec);
    double *Zresult_error_tmp = cpl_matrix_get_data(Zernike_coeff_error_vec);
    
		for(i=0; i<ZernikeNum; i++)  Zernike_coefficient_array[i]=0.0;
    for (i = 0; i <ury; i++) {
        Zernike_coefficient_array[i] = Zresult_tmp[i];
        Zernike_coefficient_array_error[i] = Zresult_error_tmp[i];
    }
    
    cpl_matrix_delete(Zernike_coeff_vec);
    cpl_matrix_delete(Zernike_coeff_error_vec);

    }else {

    if(PRINTMSG) printf("--> S2Zer: %d X %d, slopes:%d X 1 %d\n", 
                        (int)cpl_matrix_get_nrow (Z2Sinv), (int)cpl_matrix_get_ncol (Z2Sinv), 
			(int)cpl_matrix_get_nrow(mesvec_matrix), Slopes_size); 
    error =-3;
    }
    

    /* free up the  memory */
    cpl_matrix_delete(mesvec_matrix);
    cpl_matrix_delete(mesvec_error_matrix);
    cpl_matrix_unwrap(Z2Sinv);
    cpl_image_delete(buffer);
    cpl_image_delete(inv_of_Zernike_deriv_data);

    return error;
}



/* ----------------------------------------------------------------------*/
/**
 * @internal
 * @brief  Outputs the Zernike coefficient array corresponding input slopes array.
 * @parame slopes slopes vector array corresponding the quasi aberrated wavefront
 * @parame slopes_error slopes error vector array
 * @parame Zernike_coefficient_array (output) Zernike coefficients corresponding the quasi aberrated wavefront
 * @parame Zernike_coefficient_array_error (output) Zernike coefficient error corresponding
 *          the quasi aberrated wavefront
 * @return error
 * @see    gvacqAberrationSensor,gvacqAberrationSensor_spots2slopes
 *
 * It measures exactly 68 Zernike coefficients here (if you want measure more Zernikes
 * need to change the the inverse of Zernike to slopes matrix )
 * @par Error Handling:
 * CPL_ERROR_NULL_INPUT if(one of ) the input pointer (s) is NULL
 * CPL_ERROR_NULL
 * CPL_ERROR_INCOMPATIBLE_INPUT  if AB_inv_Z2S_68_136Fits dims not match
 */

cpl_error_code gvacqAberrationSensorSlopes2ZernikeFocusOnly(double *slopes,
                                                   double *slopes_error, 
                                                   char * AB_inv_Z2S_68_136Fits,
						   int Slopes_size,
						   int ABSlensletSize,
                                                   double *Zernike_coefficient_array, 
						   double *Zernike_coefficient_array_error) {
    
    /* test entries */
    cpl_error_code error = CPL_ERROR_NONE;
    int i;
	
    if (AB_inv_Z2S_68_136Fits == NULL){
        printf("ERROR::-> Please make sure database reading AB_inv_Z2S_68_136.fits correctly. \
Error occured at %s:%s:%d\n", 
               __FILE__, __FUNCTION__, __LINE__);
        
        return -1;
    }
    
    if (slopes == NULL || slopes_error == NULL){
        printf("ERROR::-> Input slopes are NULL. Error occured at %s:%s:%d\n", 
               __FILE__, __FUNCTION__, __LINE__);
        return -2;
    }

    int llx, lly, urx, ury;
    cpl_image * inv_of_Zernike_deriv_data = cpl_image_load(AB_inv_Z2S_68_136Fits,
                                                           CPL_TYPE_DOUBLE, 1, 0); /* Loads smallr invZ2S matrix*/
    if (inv_of_Zernike_deriv_data == NULL ){
    printf("ERROR:: -> AB_inv_Z2S_68_136.fits INCOMPATABLE input file from database.\n");
    printf("The file path submitted: %s \n", AB_inv_Z2S_68_136Fits);
    
    cpl_image_delete(inv_of_Zernike_deriv_data);
    return -2;
    }
    
    error=gvacqAberrationSensorZ2SLoad(ABSlensletSize,
				       &llx, &lly, &urx, &ury);
 
    ury=4; /* Which allows to measure only 4 Zernike coefficients*/

    cpl_image *buffer=  cpl_image_extract (inv_of_Zernike_deriv_data, llx, lly, urx, ury);
    if (buffer == NULL || error !=0 ){
    printf("ERROR: gvacqAberrationSensorZ2SLoad failed to measurement correct dimension\n");
    if(buffer != NULL) cpl_image_delete(buffer);
    cpl_image_delete(inv_of_Zernike_deriv_data);
    return -3;
    }
   
    int size_y = (int)cpl_image_get_size_y(buffer);
    int size_x = (int)cpl_image_get_size_x(buffer);
    
    double *inv_Z2S_double_Extracted = cpl_image_get_data_double(buffer);
    cpl_matrix *Z2Sinv = cpl_matrix_wrap(size_y, size_x, inv_Z2S_double_Extracted);
    
  
    if (Z2Sinv == NULL){
    printf("ERROR: gvacqAberrationSensorZ2SLoad failed to measurement correct dimension\n");
    
    if(Z2Sinv!=NULL) cpl_matrix_unwrap(Z2Sinv);
    cpl_image_delete(buffer);
    cpl_image_delete(inv_of_Zernike_deriv_data);
    return -3;
    }
  
    /* slopes are copying to a cpl_matrix formalism
       for multiplication requirement
    */
    cpl_matrix *mesvec_matrix = cpl_matrix_new(Slopes_size, 1);
    cpl_matrix *mesvec_error_matrix = cpl_matrix_new(Slopes_size, 1);
	
    /* fill the slopes matrix with zeros */
    for (i = 0; i < Slopes_size; i++) {
        cpl_matrix_set(mesvec_matrix, i, 0, slopes[i]);
        cpl_matrix_set(mesvec_error_matrix, i, 0, slopes_error[i]);
    } 
    
    /* Check whether the inverse Zernike derivative matrix compatible */
    if((int)cpl_matrix_get_ncol (Z2Sinv) == (int)cpl_matrix_get_nrow(mesvec_matrix)){
    cpl_matrix * Zernike_coeff_vec = cpl_matrix_product_create(Z2Sinv,
                                                               mesvec_matrix);
    
    cpl_matrix * Zernike_coeff_error_vec = cpl_matrix_product_create(Z2Sinv,
								     mesvec_error_matrix);
    
    /* Returing the Zernike coefficients to the calling function */
    double *Zresult_tmp = cpl_matrix_get_data(Zernike_coeff_vec);
    double *Zresult_error_tmp = cpl_matrix_get_data(Zernike_coeff_error_vec);
    
		for(i=0; i<ZernikeNum; i++)  Zernike_coefficient_array[i]=0.0;
    for (i = 0; i <ury; i++) {
        Zernike_coefficient_array[i] = Zresult_tmp[i];
        Zernike_coefficient_array_error[i] = Zresult_error_tmp[i];
    }
    
    cpl_matrix_delete(Zernike_coeff_vec);
    cpl_matrix_delete(Zernike_coeff_error_vec);

    }else {

    if(PRINTMSG) printf("--> S2Zer: %d X %d, slopes:%d X 1 %d\n", 
                        (int)cpl_matrix_get_nrow (Z2Sinv), (int)cpl_matrix_get_ncol (Z2Sinv), 
			(int)cpl_matrix_get_nrow(mesvec_matrix), Slopes_size); 
    error =-3;
    }
    

    /* free up the  memory */
    cpl_matrix_delete(mesvec_matrix);
    cpl_matrix_delete(mesvec_error_matrix);
    cpl_matrix_unwrap(Z2Sinv);
    cpl_image_delete(buffer);
    cpl_image_delete(inv_of_Zernike_deriv_data);

    return error;
}



/* ----------------------------------------------------------------------*/
/**
 * @internal 
 * @brief  Estimates image shift of extended scene using correlation algorithm
 * @param  image target 16x16 pixel
 * @param  template reference 16x16 pixel
 * @param  InterpSampling level of sampling. ex. 1, 0.5, 0.1 .. 
 * @param  key type of algorithm to use: 1 CCI, 2 SSD, 3 ZSSD
 * @param  XYcen out put [x, y]
 * @return status error 
 */
cpl_error_code CorrCen(cpl_image *image, cpl_image * template,
                       double InterpSampling, int key, double *XYcen) {

    /* test entries */
    /* cpl_ensure_code(image  || template , CPL_ERROR_NULL_INPUT); */
    
    if (image== NULL || template == NULL){
        printf("Input args are NULL at %s:%s:%d\n",  __FILE__, __FUNCTION__, __LINE__);
        return -1;
    }
    
   
    
    cpl_error_code error = CPL_ERROR_NONE;
    int sizeX, sizeY, PadsizeX;
    int i, j, itr;
    double corr, step = 1;
    int pis_rejected = 0;
    double CenPeak[2];
    cpl_image* corrpeakResz;
    cpl_image * CorrPeak;
    double Centers[(int) InterpSampling * 2];
    int win;
    cpl_image * ZeroPadImage;
    cpl_image *copyImage;
    int offset = 1;
    double XYcenErr[2];
    
    for(i=0; i<4; i++) XYcen[i] = 0.0;
    
    sizeX = cpl_image_get_size_x(image);
    sizeY = cpl_image_get_size_x(template);
    win = sizeX / 2;
    
   
    if (sizeX != sizeY){
        printf("Works only with square images at %s:%s:%d\n",  __FILE__, __FUNCTION__, __LINE__);
        return -1;
    }
    /* cpl_ensure_code(sizeX==sizeY, CPL_ERROR_ILLEGAL_INPUT); */
    

    
    /* padding zeros around 16x16 px image */
    ZeroPadImage = cpl_image_new(sizeX + win, sizeX + win, CPL_TYPE_DOUBLE);
    PadsizeX = (int) cpl_image_get_size_x(ZeroPadImage);
    
    error = cpl_image_multiply_scalar(ZeroPadImage, 0.0);
    if (error != CPL_ERROR_NONE) cpl_msg_error(cpl_func, "Error occured at %s:%d\n", 
					       __FILE__, __LINE__);
    
    for (i = 1; i <= sizeX; i++) {
        for (j = 1; j <= sizeX; j++) {
            error = cpl_image_set(ZeroPadImage, i + win / 2, j + win / 2,
                                  cpl_image_get(image, i, j, &pis_rejected));
            
            if (error != CPL_ERROR_NONE) cpl_msg_error(cpl_func, "Error occured at %s:%d\n", 
						       __FILE__, __LINE__);
            
        }
    }
    
    for (itr = 0; itr < (int) InterpSampling; itr++) {
        CorrPeak = cpl_image_new(sizeX + win, sizeX + win, CPL_TYPE_DOUBLE);
        
        for (i = win + offset; i <= PadsizeX - (win + 1 + offset); i++) {
            for (j = win + offset; j <= PadsizeX - (win + 1 + offset); j++) {
                step = itr / InterpSampling;
                
                if (InterpSampling == 1 || itr == 0 ) {
                    copyImage = cpl_image_extract(ZeroPadImage, i - win + 1, j - win + 1,
                                                  i + win, j + win);
                    error = cpl_image_subtract(copyImage, template);
                    error = cpl_image_power(copyImage, 2.0);
                    corr = cpl_image_get_absflux(copyImage);
                    if (error != CPL_ERROR_NONE)
                        cpl_msg_error(cpl_func, "Error occured at %s:%d\n", __FILE__,
                                      __LINE__);
                    corr = 1. / (corr + 1e-100);
                    cpl_image_delete(copyImage);
                } else{
                    corr= corr_coeff_with_interp(ZeroPadImage, template, i, j, win,
                                                 step, key);
                }
                error = cpl_image_set(CorrPeak, i, j, corr); 
                if (error != CPL_ERROR_NONE)
                    cpl_msg_error(cpl_func, "Error occured at %s:%d\n", __FILE__,
                                  __LINE__);
            }/* i */
        }/* j */
        
        corrpeakResz = cpl_image_extract(CorrPeak, PadsizeX / 2 - win + 1, PadsizeX
                                         / 2 - win + 1, PadsizeX / 2 + win, PadsizeX / 2 + win);
        
        /* cpl_ensure_code(corrpeakResz , CPL_ERROR_NULL_INPUT); */
        
        error = ParabolaFit(corrpeakResz, CenPeak);
        if (error != CPL_ERROR_NONE)
            cpl_msg_error(cpl_func, "Error occured at %s:%d\n", __FILE__, __LINE__);
        
        Centers[itr] = CenPeak[0] + step;
        Centers[itr + (int) InterpSampling] = CenPeak[1] + step;
        cpl_image_delete(CorrPeak);
        cpl_image_delete(corrpeakResz);
    }/* itr */
    
    for (i = 0; i < (int) InterpSampling; i++) {
        XYcen[0] += Centers[i];
        XYcen[1] += Centers[i + (int) InterpSampling];
    }
    
    XYcen[0] = XYcen[0] / InterpSampling;
    XYcen[1] = XYcen[1] / InterpSampling;
    
    double Flux = cpl_image_get_absflux(ZeroPadImage);
    double ron = 12; /* detector read out noise*/
    error = gvacqCorr_center_error(ZeroPadImage, Flux, ron, XYcenErr);
    if (error != CPL_ERROR_NONE) {
        cpl_msg_error(cpl_func, "Error found at %s:%s:%d\n", __FILE__,
                      __FUNCTION__, __LINE__);
    }
    
    XYcen[2] = XYcenErr[0];
    XYcen[3] = XYcenErr[1];
    cpl_image_delete(ZeroPadImage);
    
    return error;
}

/* ----------------------------------------------------------------------*/
/**
 * @internal 
 * @brief  Estimates correlation coefficient between target and reference image
 * @param  image target 16x16 pixel
 * @param  template reference 16x16 pixel
 * @param  i ith pixel
 * @param  j jth pixel
 * @param  win window size
 * @param  step sampling size
 * @param  key type of algorithm to use: 1 CCI, 2 SSD, 3 ZSSD
 * @return correlation coefficient
 */
double corr_coeff_with_interp(cpl_image *Image, cpl_image *template, int i,
				      int j, int win, double step, int key) {
    
    
    if (Image == NULL || template == NULL ){
    printf("Error occured at %s:%s:%d\n",  __FILE__, __FUNCTION__, __LINE__);
    return -2;
    }

    /* initialize  */
    double corr_coeff = 0.0;
    cpl_image *InterpImage;
    int nx, ny, nxT, nyT;
    double fluxI, fluxT;
    cpl_error_code error=CPL_ERROR_NONE;
    
    nx = cpl_image_get_size_x(Image);
    ny = cpl_image_get_size_y(Image);
    nxT = cpl_image_get_size_x(template);
    nyT = cpl_image_get_size_y(template);
    double data[(int) nxT * (int) nyT];
    
    cpl_image *copytemplate = cpl_image_duplicate(template);
    cpl_image *copyImage = cpl_image_duplicate(Image);

    if (copytemplate == NULL || copyImage == NULL ){
        printf("Null input. Error occured at %s:%d\n", __FUNCTION__, __LINE__);
        
        return -2;
    }

    /*test entries */
    if ((int) nx != (int) ny && (int) nxT != (int) nyT) {
        cpl_msg_error(cpl_func, "It works only with only square images \n");
        
        return -3;
    }
    
    error = splineResample(copyImage, 0, i - win + 1 + step, i + win + step, 
                           j-win+1+step, j+win+ step, (int) nxT, data);
    
    if (error != CPL_ERROR_NONE)
        cpl_msg_error(cpl_func, "Error occured at %s:%d\n", __FILE__, __LINE__);
    
    InterpImage = cpl_image_wrap_double(nxT, nxT, data);
    
    if (key == 1) { /* corr image domain */
        error = cpl_image_multiply(InterpImage, template);
        corr_coeff = cpl_image_get_absflux(InterpImage);
    }
    
    if (key == 2) {/* SSD */
        fluxI = cpl_image_get_absflux(InterpImage);
        fluxT = cpl_image_get_absflux(copytemplate);
        error = cpl_image_multiply_scalar(copytemplate, fluxI / fluxT);
        error = cpl_image_subtract(InterpImage, copytemplate);
        error = cpl_image_power(InterpImage, 2.0);
        corr_coeff = cpl_image_get_absflux(InterpImage);
        corr_coeff = 1. / (corr_coeff + 1e-100);
    }
    
    if (key == 3) {/* ZSSD */
        fluxI = cpl_image_get_absflux(InterpImage);
        fluxT = cpl_image_get_absflux(copytemplate);
        error = cpl_image_multiply_scalar(copytemplate, fluxI / fluxT);
        error = cpl_image_subtract_scalar(InterpImage, cpl_image_get_mean(
                                              InterpImage));
        
        error = cpl_image_subtract_scalar(copytemplate, cpl_image_get_mean(
                                              copytemplate));
        error = cpl_image_subtract(InterpImage, copytemplate);
        error = cpl_image_power(InterpImage, 2.0);
        
        corr_coeff = cpl_image_get_absflux(InterpImage);
        corr_coeff = 1. / (corr_coeff + 1e-100);
    }
    
    /* free the buffers */
    cpl_image_unwrap(InterpImage);
    cpl_image_delete(copytemplate);
    cpl_image_delete(copyImage);
    
    return corr_coeff;
}/* corr_coeff_with_interp */


cpl_error_code image_roll_double(double * image, int im_x)
{
    cpl_error_code error=CPL_ERROR_NONE;
    int i;
    long ioff[2] = { im_x / 2 -1, im_x / 2 -1};
    long istd[2] = { 1, im_x };
    long ndims[2] = { im_x, im_x };
    long nn2[2] = { im_x, 1 };
    double ws[im_x];
    
    for (i = 0; i < im_x; i++) ws[i] = 0.0;
    
    for (i = 0; i < 2; i++){
        error = roll2d(image, ioff[i], istd[i], ndims[i], nn2[i], ws);
        if (error != CPL_ERROR_NONE) cpl_msg_error(cpl_func, "Error at roll2d \n");
    }
    
    return error;
}

/* @internal 
   @brief  Estimates image shift of extended scene using FFT correlation algorithm 
   with window shift
   @param temp template image
   @param template_shifted shifted template
   @param im live image
   @param off template image is offsetted
   @param XYCEN return image shifts w.r.t template
*/
cpl_error_code FFTCorrCenInterp(cpl_image *temp, cpl_image * template_shifted, 
				cpl_image * im, double off, double *XYCEN){
    cpl_error_code error;
    
    double XYcen[4];
    double XYcen1[4];
    
    error = FFTCorrCen(temp, im, XYcen);
    error = FFTCorrCen(template_shifted, im, XYcen1);
    
    XYCEN[0]=(XYcen[0]+XYcen1[0]-off)/2.0;
    XYCEN[1]=(XYcen[1]+XYcen1[1]-off)/2.0;

    XYCEN[2]=XYcen[2];
    XYCEN[3]=XYcen[3];

    return error;
}

/* ----------------------------------------------------------------------*/
/**
 * @internal 
 * @brief  Estimates image shift of extended scene using FFT correlation algorithm
 * @param  Ref reference 16x16 pixel
 * @param  image target 16x16 pixel
 * @param  XYcen out put [x, y]
 * @return status error 
 */
cpl_error_code FFTCorrCen(cpl_image *ref, cpl_image *im, double *XYcen) {
	int i;
	cpl_error_code error = CPL_ERROR_NONE;
	/* cpl_ensure_code(Ref  || image , CPL_ERROR_NULL_INPUT); */
	
	if (ref == NULL) {
        printf("Input argument is NULL \n");
	return -1;
	}
    
	int im_x = cpl_image_get_size_x(ref);
	int im_y = cpl_image_get_size_y(im);
    
	if (im_x != im_y) {
		cpl_msg_error(cpl_func,
                      "input image size should be squared: im_x = %d, im_y=%d\n", im_x, im_y);
		return -1;
	}

    /*Images are zero padded to cancel out aliasing errors */
	im_x = 2*im_x;
	im_y = 2*im_y;
	
	cpl_image *Ref=cpl_image_new(im_x, im_y, CPL_TYPE_DOUBLE);
	cpl_image *image= cpl_image_new(im_x, im_y, CPL_TYPE_DOUBLE);

	error = cpl_image_copy (Ref, ref, im_x/2-im_x/4, im_y/2-im_y/4);
	error = cpl_image_copy (image, im, im_x/2-im_x/4, im_y/2-im_y/4);
/* 
   The flow chart of finding the centroid
   template= real1 + j complex1
   image =   real2 + j coplex2
   
   here, real1= Ref, real2=image
   
   FFT(template) = real1 + j complex1 = cpl_fft_image(real1, complex1, CPL_FFT_DEFAULT);
   FFT(image) =    real2 + j complex2 = cpl_fft_image(real2, complex2, CPL_FFT_DEFAULT);
   The original data is replaced by the data after fourier transform in the above function

   conj(FFT(template)) = real1 - j complex1
   conj(FFT(template)) * FFT(image) = real1 X real2 + complex1 X complex2 +
   j (real1 X complex2 - real2 X complex1 );
   
   Re_part= real1Xreal2 + complex1Xcomplex2;
   Im_part  = real1Xcomplex2 -real2Xcomplex1;
   
   -1
   FFT  { conj(FFT(template))*FFT(image) } = cpl_fft_image(Re_part, Im_part, CPL_FFT_INVERSE);
   
   Correlation peak = roll2(Re_part);
   
   Finally find the parabolic fit on the Correlation peak
*/


	cpl_image * complex1 = cpl_image_new(im_x, im_y, CPL_TYPE_DOUBLE);
	cpl_image * complex2;
    
	for (i = 1; i <= im_x; i++) {
		error = cpl_image_set(complex1, i, i, 0.0);
	}
    
	complex2 = cpl_image_duplicate(complex1);
	/* cpl_ensure_code(complex2, CPL_ERROR_NULL_INPUT); */
    
	error = cpl_image_fft(Ref, complex1, CPL_FFT_DEFAULT); /* Forward tranform */
	error = cpl_image_fft(image, complex2, CPL_FFT_DEFAULT);/* Forward tranform */
    
	if (error != CPL_ERROR_NONE)
		cpl_msg_error(cpl_func, "Error occured at %s:%d\n", __FILE__, __LINE__);
    
	cpl_image * real1Xreal2 = cpl_image_multiply_create(Ref, image);
	cpl_image * complex1Xcomplex2 = cpl_image_multiply_create(complex1, complex2);
	cpl_image * real1Xcomplex2 = cpl_image_multiply_create(Ref, complex2);
	cpl_image * real2Xcomplex1 = cpl_image_multiply_create(image, complex1);
	cpl_image * Re_part = cpl_image_add_create(real1Xreal2, complex1Xcomplex2);
	cpl_image * Im_part = cpl_image_subtract_create(real1Xcomplex2, real2Xcomplex1);

	/* Inverse tranform */
	error = cpl_image_fft(Re_part, Im_part, CPL_FFT_INVERSE);

	if (error != CPL_ERROR_NONE)
		cpl_msg_error(cpl_func, "Error occured at %s:%d\n", __FILE__, __LINE__);
    
	/*
	 * Roll the image to the center
	 */
	double * data_corr = cpl_image_get_data_double(Re_part);
	error= image_roll_double(data_corr, im_x);
	cpl_image * corr_peak_im = cpl_image_wrap_double(im_x, im_y, data_corr);
	cpl_image *corr_peak = cpl_image_extract (corr_peak_im, 
                                              im_x/2-im_x/4+1, im_y/2-im_y/4+1, 
                                              im_x/2+im_x/4, im_y/2+im_y/4);

	/* Parabolic fit */
	error = ParabolaFit(corr_peak, XYcen);

	double Flux = cpl_image_get_absflux(image);
	double ron = RON; /* detector read out noise*/
	double XYcenErr[2];

	error = gvacqCorr_center_error(corr_peak, Flux, ron, XYcenErr);
	if (error != CPL_ERROR_NONE) {
		cpl_msg_error(cpl_func, "Error found at %s:%s:%d\n", __FILE__,__FUNCTION__,__LINE__);
	}

	XYcen[2] = XYcenErr[0];
	XYcen[3] = XYcenErr[1];

	/* Free buffer */
	cpl_image_unwrap(corr_peak_im);
	cpl_image_delete(corr_peak);
	cpl_image_delete(Ref);
	cpl_image_delete(image);
	cpl_image_delete(Re_part);
	cpl_image_delete(Im_part);
	cpl_image_delete(real1Xcomplex2);
	cpl_image_delete(real2Xcomplex1);
	cpl_image_delete(real1Xreal2);
	cpl_image_delete(complex1Xcomplex2);
	cpl_image_delete(complex1);
	cpl_image_delete(complex2);
    
	return error;
}/* FFTCorrCen */

/*
  It rolls selected dimensions of the image array.  
*/
cpl_error_code roll2d(double * a, long ioff, long istd, 
                      long n, long n2, double * ws) {
#define FUNCws_1(a1) ws[a1-1]
#define FUNCa_3(a1,a2,a3) a[a1-1+istd*(a2-1+n*(a3-1))]
    
	/*
	 *     a      double precision  FUNCa_3(istd,n,n2)
	 *            to be "rolled", in place, on its middle dimension
	 *
	 *     ioff   the offset by which to "roll" a_3 (iabs(ioff) .lt. n)
	 *            for ioff>0:
	 *               output  FUNCa_3(1+ioff)=  input  FUNCa_3(1)
	 *               output  FUNCa_3(2+ioff)=  input  FUNCa_3(2)
	 *                      ...         ...
	 *               output  FUNCa_3(n-1)=     input  FUNCa_3(n-ioff-1)
	 *               output  FUNCa_3(n)=       input  FUNCa_3(n-ioff)
	 *               output  FUNCa_3(1)=       input  FUNCa_3(n-ioff+1)
	 *               output  FUNCa_3(2)=       input  FUNCa_3(n-ioff+2)
	 *                      ...         ...
	 *               output  FUNCa_3(ioff-1)=  input FUNCa_3(n-1)
	 *               output  FUNCa_3(ioff)=    input FUNCa_3(n)
	 *            for ioff<0:
	 *               output  FUNCa_3(1)=        input FUNCa_3(1-ioff)
	 *               output  FUNCa_3(2)=        input FUNCa_3(2-ioff)
	 *                      ...          ...
	 *               output  FUNCa_3(n+ioff-1)= input FUNCa_3(n-1)
	 *               output  FUNCa_3(n+ioff)=   input FUNCa_3(n)
	 *               output  FUNCa_3(n+ioff+1)= input FUNCa_3(1)
	 *               output  FUNCa_3(n+ioff+2)= input FUNCa_3(2)
	 *                      ...          ...
	 *               output  FUNCa_3(n-1)=      input FUNCa_3(-ioff-1)
	 *               output  FUNCa_3(n)=        input FUNCa_3(-ioff)
	 *
	 *     istd   stride of values in a, that is, length of dimensions
	 *            of a before the one of length n
	 *
	 *     n2     length of dimensions of a after the one of length n
	 *
	 *     ws     double precision ws_1(n) of working space
	 **/

	long i, i2, j, joff;
	cpl_error_code error = CPL_ERROR_NONE;
    
	if (ioff == 0)
		return -1;
	if (ioff < 0) {
		joff = n + ioff;
		for (i2 = 1; i2 <= n2; i2 += 1) {
			for (i = 1; i <= istd; i += 1) {
				for (j = 1; j <= -ioff; j += 1) {
					FUNCws_1(j+joff) = FUNCa_3(i, j, i2);
				}
				for (j = 1; j <= joff; j += 1) {
					FUNCws_1(j) = FUNCa_3(i, j-ioff, i2);
				}
				for (j = 1; j <= n; j += 1) {
					FUNCa_3(i, j, i2) = FUNCws_1(j);
				}
			}
		}
	} else {
		joff = n - ioff;
		for (i2 = 1; i2 <= n2; i2 += 1) {
			for (i = 1; i <= istd; i += 1) {
				for (j = 1; j <= joff; j += 1) {
					FUNCws_1(j+ioff) = FUNCa_3(i, j, i2);
				}
				for (j = 1; j <= ioff; j += 1) {
					FUNCws_1(j) = FUNCa_3(i, j+joff, i2);
				}
				for (j = 1; j <= n; j += 1) {
					FUNCa_3(i, j, i2) = FUNCws_1(j);
				}
			}
		}
	}
	
	return error;
}/* roll2d */


/* ----------------------------------------------------------------------*/
/**
 * @internal
 * @brief  Parabolic fit using 5 pixels around the maximum 
 * @param  corr_peak   correlation peak image
 * @param  XYcen (output) subpixel centroid 
 * @return error
 *
 * Formulation adapted form Eq. 33,of S.Thomas et al., MNRS, 371,323-336 (2006).
 * 
 * @see Corr, gvacqAberrationSensor
 *
 * @par Error Handling:
 * CPL_ERROR_NULL_INPUT if(one of ) the input pointer (s) is NULL
 * CPL_ERROR_NULL
 */

cpl_error_code ParabolaFit(cpl_image *corr_peak, double *XYcen) {

    /* cpl_ensure_code(corr_peak  , CPL_ERROR_NULL_INPUT); */
    
    cpl_error_code error = CPL_ERROR_NONE;

#ifdef VLT2014
    cpl_size px, py;
#else 	 
	int px, py;
#endif

    
    int s = 0;
    double I01, I11, I21, I10, I12;

    if (corr_peak == NULL) {
        printf("Input argument is NULL \n");
        return -1;
    }

    int nx = cpl_image_get_size_x(corr_peak);
    int ny = cpl_image_get_size_y(corr_peak);
    
    if (nx != ny) {
        printf("Works only with  square sub images \n");
    return -1;
    }
    
    
    XYcen[0] = 0.0;
    XYcen[1] = 0.0;
    
   
    error = cpl_image_get_maxpos(corr_peak, &px, &py); /*px*/
    if (px == 1)
        I01 = cpl_image_get(corr_peak, nx, py, &s);
    else
        I01 = cpl_image_get(corr_peak, px - 1, py, &s);
    if (px == nx)
        I21 = cpl_image_get(corr_peak, 1, py, &s);
    else
        I21 = cpl_image_get(corr_peak, px + 1, py, &s);
    if (py == 1)
        I10 = cpl_image_get(corr_peak, px, ny, &s);
    else
        I10 = cpl_image_get(corr_peak, px, py - 1, &s);
    if (py == ny)
        I12 = cpl_image_get(corr_peak, px, 1, &s);
    else
        I12 = cpl_image_get(corr_peak, px, py + 1, &s);
    
    I11 = cpl_image_get(corr_peak, px, py, &s);
    XYcen[0] = px + 0.5 * (I01 - I21) / (I01 + I21 - 2 * I11);
    XYcen[1] = py + 0.5 * (I10 - I12) / (I10 + I12 - 2 * I11);
    
    return error;
}/* ParabolaFit */


/* ----------------------------------------------------------------------*/
/**
 * @internal
 * @brief  Measures centroid error of correlation algorithm 
 * @param  corr_peak   correlation peak image
 * @param  Flux flux of the spot
 * @param  ron detector Read Out Noise
 * @param  XYcenErr (output) centroid error
 * @return error
 *
 * Formulation adapted form Eq. 34, 35 of S.Thomas et al., MNRS, 371,323-336 (2006).
 * 
 * @see    Corr, gvacqAberrationSensor
 *
 * @par Error Handling:
 * CPL_ERROR_NULL_INPUT if(one of ) the input pointer (s) is NULL
 * CPL_ERROR_NULL
 */
cpl_error_code gvacqCorr_center_error(cpl_image *corr_peak, double Flux,
		double ron, double *XYcenErr) {
	/* test entries */
	/* cpl_ensure_code(corr_peak, CPL_ERROR_NULL_INPUT); */
    
	if (corr_peak == NULL) {
        printf("Input argument is NULL \n");
        return -1;
	}
    
#ifdef VLT2014
    cpl_size px, py;
#else 	 
	int px, py;
#endif
	cpl_error_code error = CPL_ERROR_NONE;
	double fwhm_x, fwhm_y;
	
    
	int im_size_x=cpl_image_get_size_x(corr_peak);
	
	if(im_size_x*16 /*RON + Background ?*/ > cpl_image_get_absflux(corr_peak))
	    {
            printf("No star data found in the image\n");
            
            return -3;
	    }
	
	error = cpl_image_get_maxpos(corr_peak, &px, &py); /*px*/
	error = cpl_image_get_fwhm(corr_peak, px, py, &fwhm_x, &fwhm_y); /*px*/
    
	if (error != CPL_ERROR_NONE)
		cpl_msg_error(cpl_func, "Error occured at %s:%d\n", __FILE__, __LINE__);
    
	double sigma_ron_x = 2 * fwhm_x * fwhm_x * ron / Flux; /*px*/
	double sigma_ron_y = 2 * fwhm_y * fwhm_y * ron / Flux;
	double sigma_ph = PI / sqrt(2 * log(2) * Flux); /*px*/
	XYcenErr[0] = sqrt(sigma_ron_x * sigma_ron_x + sigma_ph * sigma_ph); /*px*/
	XYcenErr[1] = sqrt(sigma_ron_y * sigma_ron_y + sigma_ph * sigma_ph); /*px*/
	
	return error;
}

/*
 * @internal
 * @brief  Generates a reference/template image from the aberration sensor image itself
 * @param  fimage shack-hartmann image
 * @param  imistart2 x image locations
 * @param  imjstart2 y image locations
 * @param  nsubs total number of subapertures
 * @param  validsubs valid subapertures
 * @param  template (output) genarated reference 
 * @return error
 */
cpl_error_code gvacqcorr_template_image(cpl_image *fimage, int *imistart2,
		int *imjstart2, int nsubs, int *validsubs, cpl_image *template) {

	/* test entries */
	/* cpl_ensure_code(fimage   || imistart2   || imjstart2  || validsubs   , CPL_ERROR_NULL_INPUT); */
    
	if (fimage==NULL   || imistart2==NULL   || imjstart2==NULL  || validsubs==NULL) {
        printf("Input argument is NULL \n");
        return -1;
	}
    
	/* Declarations */
	cpl_error_code error = CPL_ERROR_NONE;
	int i, l;
	int nvalidsubs;
	double Num = 0;
	double xy[2];
	cpl_image *tmp;
	int shiftx, shifty;
	nvalidsubs = 0;
	for (i = 0; i < nsubs; i++) nvalidsubs += validsubs[i];
    
	for (l = 0; l < nsubs; l++) {
        if (validsubs[l]) {
            Num += 1;
            tmp = cpl_image_extract(fimage, imistart2[l] + 1, imjstart2[l] + 1,
                                    imistart2[l] + 16, imjstart2[l] + 16);
            error =ParabolaFit(tmp, xy);
            shiftx=roundd(8-xy[0]);
            shifty=roundd(8-xy[1]);
            error = cpl_image_shift(tmp, shiftx, shifty);
           
            if(tmp != NULL) {
                error = cpl_image_add(template, tmp);
            }else printf("Data not found at %s:%s:%d\n", __FILE__, __FUNCTION__, __LINE__); 
            cpl_image_delete(tmp);
        }
        
	}
	error = cpl_image_multiply_scalar(template, 1.0 / Num);
	
	if (error != CPL_ERROR_NONE)
	    cpl_msg_error(cpl_func, "Error occured at %s:%d\n", __FILE__, __LINE__);
	
	return error;
}



cpl_error_code gvacqcorr_template_image_ext(cpl_image *fimage, int *imistart2,
		int *imjstart2, int nsubs, int *validsubs, cpl_image *template) {

	/* test entries */
	/* cpl_ensure_code(fimage   || imistart2   || imjstart2  || validsubs   , CPL_ERROR_NULL_INPUT); */

	if (fimage==NULL   || imistart2==NULL   || imjstart2==NULL  || validsubs==NULL) {
        printf("Input argument is NULL \n");
        return -1;
	}

	/* Declarations */
	cpl_error_code error = CPL_ERROR_NONE;
	int i, l;
	int nvalidsubs;
	double Num = 0;
	cpl_image *tmp;
	int shiftx, shifty;
	
#ifdef VLT2014
    cpl_size px, py;
#else 	 
	int px, py;
#endif
	nvalidsubs = 0;
	for (i = 0; i < nsubs; i++) nvalidsubs += validsubs[i];
    
	for (l = 1; l < nsubs; l++) {
        if (validsubs[l]) {
			Num += 1;
			tmp = cpl_image_extract(fimage, imistart2[l] + 1, imjstart2[l] + 1,
                                    imistart2[l] + 16, imjstart2[l] + 16);
			error=cpl_image_get_maxpos(tmp, &px, &py);
			shiftx=roundd(8-px);
			shifty=roundd(8-py);
		
			error = cpl_image_shift(tmp, shiftx, shifty);
			if(tmp != NULL) {
                error = cpl_image_add(template, tmp);
			}else printf("Data not found at %s:%s:%d\n", __FILE__, __FUNCTION__, __LINE__); 
			cpl_image_delete(tmp);
		}
	
	}
	error = cpl_image_multiply_scalar(template, 1.0 / Num);
	
	if (error != CPL_ERROR_NONE)
	    cpl_msg_error(cpl_func, "Error occured at %s:%d\n", __FILE__, __LINE__);
	
	return error;
}




/* ----------------------------------------------------------------------*/
/**
 * @internal
 * @brief converts Shack-Hartmann array spots into slopes (tilts)
 * @param fimage  final image with spots
 * @param imistart2 vector of i starts of each image
 * @param imjstart2 vector of j starts of each image
 * @param nsubs  number of subapertures
 * @param binxy2  side size of the single spot image window extracted from the given image
 * @para  fimnx  the aberration sensor image X dimension
 * @param fimnx  the aberration sensor image Y dimension
 * @param yoff  y offset (to process only part of the image)
 * @param centroidw centroid weight vector for centroid computations, X & Y
 * @param validsubs valid subaps within the set of subap for which image is computed
 * @param mesvec final slope measurement vector
 * @param mesvec_error slope measurement error vector
 *
 *
 * pass one aberration sensor image buffer array
 * and a set of indices (start and end of each subapertures). Then it takes the each subaperture
 * spot and computes the weighted centriod and keeps the measurement in the slopes array vector.
 * flow:
 * 1. select the subaperture spot
 * 2. measure the X and Y centriod positions
 * 3. keep the centriod in the mesvec
 * 4. free memory
 * This rotine adopted from yao library version 4.8.2
 * Copyrights with  F. Rigaut
 *
 * @par Error Handling:
 * CPL_ERROR_NULL_INPUT if(one of ) the input pointer (s) is NULL
 * CPL_ERROR_NULL
 */
int gvoacqAberrationSensorSpots2Slopes(double *fimage, int *imistart2,
                                       int *imjstart2, int nsubs, int binxy2, 
				       int fimnx, int fimny, int yoff,
                                       double *centroidw, double *threshold, 
				       int *validsubs, double *mesvec,
                                       double *mesvec_error)
    
{
	/* test entries */
	cpl_ensure_code(fimage && imistart2 && imjstart2, CPL_ERROR_NULL_INPUT);

	/* Declarations */
	double *bimage2;
	double *fnoise;
	double centx, centy, centi, tmpf, thback;
	long nb2;
	int i, j, k, l, koff, xyoff;
	int nvalidsubs, vsc; /* vsc = valid sub counter */
	double VarianceIntesnity[binxy2 * binxy2]; /* (sigma_I) */
	double XW_Xc, YW_Yc;
	double sigma_Xc = 0;
	double sigma2_Xc = 0;
	double sigma_Yc = 0;
	double sigma2_Yc = 0;
	double ron = 0.0;

	nb2 = binxy2 * binxy2;
	vsc = -1;
	xyoff = yoff * fimnx;

	cpl_ensure_code(fimage, CPL_ERROR_NULL_INPUT);

	/* Allocate memory for the input operands and check its availability. */
	bimage2 = (double*) malloc(nb2 * sizeof(double));
	fnoise = (double*) malloc(fimnx * fimny * sizeof(double));

	cpl_ensure_code( bimage2 || fnoise, CPL_ERROR_NULL_INPUT);
	/* if ( bimage2 == NULL || fnoise == NULL ) { return (1); } */

	/* how many valid subs? */
	nvalidsubs = 0;
	for (i = 0; i < nsubs; i++)
		nvalidsubs += validsubs[i];
    
	thback = (double) threshold[nsubs];
	for (i = xyoff; i < (xyoff + fimnx * fimny); i++)
		fimage[i] -= thback;
    
	for(i=0; i<nsubs; i++){
        imistart2[i]=imistart2[i]-8;
        imjstart2[i]=imjstart2[i]-8;
	}
    
	/* OK, now let's finally compute the slopes over all subaps */
	for (l = 0; l < nsubs; l++) {
		if (validsubs[l])
			vsc++;
        
		/* retrieve the final/integrated subaperture spot image from fimage */
		koff = imistart2[l] + imjstart2[l] * fimnx;
		for (j = 0; j < binxy2; j++) {
			for (i = 0; i < binxy2; i++) {
				k = koff + i + j * fimnx;
				bimage2[i + j * binxy2] = fimage[k];
			}
		}
        
		/* Apply threshold */
		tmpf = (double) threshold[l];
		for (i = 0; i < nb2; i++) {
			/* see above for thback meaning */
			bimage2[i] = bimage2[i] + thback - tmpf;
			if (bimage2[i] < 0.0f)
				bimage2[i] = 0.0f;
		}

		/* is that a valid subaperture? */
		if (validsubs[l]) {
			/* Compute centroids */
			centx = centy = centi = 0.0f;
			for (i = 0; i < binxy2; i++) {
				for (j = 0; j < binxy2; j++) {
					centx += centroidw[i] * bimage2[i + j * binxy2]; /* sum(xWI) */
					centy += centroidw[j] * bimage2[i + j * binxy2];
					centi += bimage2[i + j * binxy2]; /* sum(I) */
				}
			}
            
			for (i = 0; i < binxy2; i++) {
				for (j = 0; j < binxy2; j++) {
					VarianceIntesnity[i + j * binxy2] = bimage2[i + j * binxy2] + ron
                        * ron;
					/* sigma^2_Inoise = Iphoton + sigma^2_ron */
                    
					XW_Xc = centroidw[i] - centx / centi; /* xW -Xc */
					YW_Yc = centroidw[j] - centy / centi;

					sigma_Xc = VarianceIntesnity[i + j * binxy2] * (XW_Xc) / (centi);
					sigma2_Xc += sigma_Xc * sigma_Xc;
                    
					sigma_Yc = VarianceIntesnity[i + j * binxy2] * (YW_Yc) / (centi);
					sigma2_Yc += sigma_Yc * sigma_Yc;
				}
			}

			/* Compute measurements */
			if (centi > 0.0f) {
				mesvec[vsc] = centx / centi  + 8.0+ imistart2[l]; 
				mesvec[nvalidsubs + vsc] = (centy / centi)+ 8.0+ imjstart2[l];
				mesvec_error[vsc] = sqrt(sigma2_Xc);
				mesvec_error[nvalidsubs + vsc] = sqrt(sigma2_Yc);
				/* printf("%14.7e \n", sqrt(sigma2_Xc)); */
			} /* otherwise stay at zero (init value)
                 vsc++; increment valid subaperture counter */
		}
	}
	free(bimage2);
	free(fnoise);
	return (0);
}



/* double to round integer */
int roundd(double value)
{
    return (value >= 0) ? (int)(value + 0.5) : (int)(value - 0.5);
}



/**
 * @internal 
 * @brief  Measure of flux of Shack-Hartmann spots 
 * @param  fimage  pointer to the  Shack-Hartmann image
 * @param  imistart2 locations of spots (along x-axis)
 * @param  imjstart2 locations of spots (along y-axis)
 * @param  validsubs choose the spots to apply
 * @param  Flux Flux array for validsubs
*/
cpl_error_code gvacqAberrationSensorSpotsFlux(cpl_image *fimage,
					      int *imistart2, int *imjstart2, int *validsubs,
					      double *Flux) {
    
    if (fimage == NULL || imistart2 == NULL || imjstart2 == NULL ){
        printf("ERROR::-> Inputs for the function are NULL at %s:%s:%d\n", 
               __FILE__, __FUNCTION__, __LINE__);
        
        return -1;
    }
    
    int nsubs=ABSNUMSPOTS;
    int l, i;
    cpl_error_code error=CPL_ERROR_NONE;

    for(i=0; i<nsubs; i++){
    imistart2[i]=imistart2[i]-7;
    imjstart2[i]=imjstart2[i]-7;

    } 
   
    int nvalidsubs = 0;
    for (i = 0; i < nsubs; i++) nvalidsubs += validsubs[i];
    int vsc=-1;

    for (l = 0; l < nvalidsubs; l++) {
    vsc++;
    Flux[vsc] =  cpl_image_get_absflux_window(fimage, imistart2[l] + 1, imjstart2[l] + 1,
                                 imistart2[l] + 14, imjstart2[l] + 14);
   
    }
   
    return error;
}

/**
 * @internal 
 * @brief  Measuring effective lenslet size of Shack-Hartmann based on flux of spots
 * @param  detector_image   pointer to the  detector image
 * @param  ABSwindow aberation sensor window centers
 * @param  ABSwindowSize window size
 * @param  ABSRefSpotsCenters positions of spots, used to measure flux of spots
 * @param  ABSLensletSizeArray measured lenslet sizes for four tel.
 * @param  ABSErrorStatus error in the measurement for each telescope
*/
cpl_error_code gvacqAberrationSensorEstimateABSLensletSize(cpl_image *detector_image, 
							   int * ABSwindow,
							   int  ABSwindowSize,
							   double *ABSRefSpotsCenters,
							   int* ABSLensletSizeArray,
							   int * ABSErrorStatus){
    cpl_error_code error=CPL_ERROR_NONE;
    int nsub = ABSNUMSPOTS; /* no of spots avilable with 9x9 SH */
    int validsubs[ABSNUMSPOTS];
    int  ABSiLoc[ABSNUMSPOTS];
    int  ABSjLoc[ABSNUMSPOTS];
    double Flux[ABSNUMSPOTS];
    int GoodSpotsFluxCheck[ABSNUMSPOTS];
    int nvalidsubs;
  

    int i, ii;
    double Threshold_flux=100;
    for(i=0; i<4; i++){
       ABSErrorStatus[i]=0.0;
       ABSLensletSizeArray[i]=8;
    }

    
    for (i = 0; i < nsub; i++) {
       if (i == 0 || i == 6 || i == 7 || i == 15 || i==38 || i == 61 || i == 69 || i == 70 || i==76) validsubs[i] = 0;
       else validsubs[i] = 1;
    }

    nvalidsubs=0;
    for (i = 0; i < nsub; i++) nvalidsubs += validsubs[i];
    for(i=0; i<nvalidsubs; i++) GoodSpotsFluxCheck[i]=1;

    for (ii = 1; ii <= 4; ii++) {/* For each telescope */
    int GoodSpots=0;

    for(i=0; i<ABSNUMSPOTS; i++) {
       ABSiLoc[i]= roundd(ABSRefSpotsCenters[(ii-1)*ABSNUMSPOTS*2 + i]); /* ABS i spots locations */
       ABSjLoc[i]= roundd(ABSRefSpotsCenters[(ii-1)*ABSNUMSPOTS*2 + ABSNUMSPOTS+ i]); /* ABS j spots locations */
    }
      
    cpl_image *window_abs = cpl_image_extract(detector_image, 
					      ABSwindow[ii-1], 
					      ABSwindow[ii+3], 
					      ABSwindow[ii-1]+ABSwindowSize-1, 
					      ABSwindow[ii+3]+ABSwindowSize-1);
 

    for(i=0; i<nvalidsubs; i++) Flux[i]=1.0;
    error = gvacqAberrationSensorSpotsFlux(window_abs,  ABSiLoc, ABSjLoc, validsubs,  Flux);
    
    
    if(Flux[29]> Threshold_flux && Flux[37]> Threshold_flux && 
       Flux[39]> Threshold_flux && Flux[47] >  Threshold_flux){
    double Flux_ADU=(Flux[29]+Flux[37]+Flux[39]+Flux[47])/4.0;
 
    for(i=0; i<nvalidsubs; i++) {
    if(Flux_ADU/(Flux[i]+1e-2) > 2.5) GoodSpotsFluxCheck[i]=0;
    }
    
    for (i = 0; i < nvalidsubs; i++) {
    GoodSpots += GoodSpotsFluxCheck[i];
    }

    if(GoodSpots >= 60) {
    ABSLensletSizeArray[ii-1]=9;
    }else if (GoodSpots <60 && GoodSpots >=52 ){
    ABSLensletSizeArray[ii-1]=8;
    }else if (GoodSpots <52 && GoodSpots >=44 ){
    ABSLensletSizeArray[ii-1]=7;
    }else if (GoodSpots <44 && GoodSpots >=36 ){
    ABSLensletSizeArray[ii-1]=6;
    }else if (GoodSpots <36 && GoodSpots >=24 ){
    ABSLensletSizeArray[ii-1]=5;
    }else if (GoodSpots <24 && GoodSpots >=20 ){
    ABSLensletSizeArray[ii-1]=4;
    }else if (GoodSpots <20 && GoodSpots >=9 ){
    ABSLensletSizeArray[ii-1]=3;
    }else{
    ABSLensletSizeArray[ii-1]=-1;
    ABSErrorStatus[ii-1]=-4;
    }

    }else { /*End of if(Flux[37]> Threshold_flux) */
    ABSLensletSizeArray[ii-1]=-1;
    ABSErrorStatus[ii-1]=-4;
    }/*End of if(Flux[37] > Threshold_flux) */
    
    if(PRINTMSG)printf("Tel: %d, GoodSpots: %d \n", ii, GoodSpots);
    
    cpl_image_delete(window_abs);
    }/*end of if(ii = 1; ii <= 4; ii++) */
     
    return error;
}

/**
 * @internal 
 * @brief  Gaussian fits over multiple spots
 * @param  CorrPeak detector image
 * @param  locX rough x-positions of spots
 * @param  locY rough y-positions of spots
 * @param  Npeaks number of spots/peaks
 * @param  window Gaussian fit window size
 * @param  FitPosition measured positions of spots [x_1, x_2, -- x_n; y_1, y_2, -- y_n]
 * @param  FitPosition_error measured positions err of spots [x_1, x_2, -- x_n; y_1, y_2, -- y_n]
*/
cpl_error_code MultipleGaussFit(cpl_image * CorrPeak, double *locX, double *locY, int Npeaks, 
				int window, double *FitPosition, double *FitPosition_error){
    
    int i;
    cpl_error_code error=CPL_ERROR_NONE;
    cpl_matrix * PeaksLoc = cpl_matrix_new(2, Npeaks); /* 1 spot (x and y) */
    cpl_matrix * xy_centre = cpl_matrix_new(2, Npeaks); /* 1 spot (x and y) */
    cpl_matrix * xy_centre_err = cpl_matrix_new(2, Npeaks);
    cpl_matrix * xy_sigma = cpl_matrix_new(2, Npeaks);
    cpl_matrix * xy_sigma_err = cpl_matrix_new(2, Npeaks);
    cpl_matrix * xy_fwhm = cpl_matrix_new(2, Npeaks);
    cpl_matrix * xy_fwhm_err = cpl_matrix_new(2, Npeaks);
    cpl_array * all_error_codes = cpl_array_new(1, CPL_TYPE_INT);
    int robustness = 7;
  
    for(i=0; i<Npeaks; i++){
    cpl_matrix_set(PeaksLoc, 0, i, locX[i]); /* Guess spot position to fit Guassin */
    cpl_matrix_set(PeaksLoc, 1, i, locY[i]); 
    }
    
    error = clipm_centroiding_multi_gauss(CorrPeak, 
                                           PeaksLoc, 
                                           window, 
                                           &xy_centre, 
                                           &xy_centre_err, 
                                           &xy_sigma, 
                                           &xy_sigma_err,
                                           &xy_fwhm, 
                                           &xy_fwhm_err, 
                                           NULL, 
                                           &all_error_codes,
                                           robustness);
    
    cpl_matrix *copy=cpl_matrix_duplicate (xy_centre);
    error = cpl_matrix_subtract(copy, PeaksLoc);	
    double *copy_double = cpl_matrix_get_data(copy);

    double avg_x=0,avg_y=0;
    for(i=0; i<Npeaks; i++){
    avg_x+=copy_double[i];
    avg_y+=copy_double[i+Npeaks];
    }
    avg_x=avg_x/((double)Npeaks);
    avg_y=avg_y/((double)Npeaks);
    
   
    for(i=0; i<Npeaks; i++) {
    if( (copy_double[i] > window/3.0 || copy_double[i] < - window/3.0) || 
	copy_double[i+Npeaks] > window/3.0 || copy_double[i+Npeaks] < - window/3.0){
   
   
    cpl_matrix_set(xy_centre, 0, i, (double)locX[i]+avg_x);
    cpl_matrix_set(xy_centre, 1, i, (double)locY[i]+avg_y);
    }
    }
   


    double *Fit = cpl_matrix_get_data(xy_centre);
    double *Fit_error = cpl_matrix_get_data(xy_centre_err);
   
    
    for(i=0; i<2*Npeaks; i++) {
    FitPosition[i] = Fit[i];
    FitPosition_error[i]=Fit_error[i];
    }

    cpl_matrix_delete(copy);
    cpl_matrix_delete(PeaksLoc);
    cpl_matrix_delete(xy_centre);
    cpl_matrix_delete(xy_centre_err);
    cpl_matrix_delete(xy_sigma);
    cpl_matrix_delete(xy_sigma_err);
    cpl_matrix_delete(xy_fwhm);
    cpl_matrix_delete(xy_fwhm_err);
    cpl_array_delete(all_error_codes);
    return error;
}

/**
 * @internal 
 * @brief  Measurment of focus using autocorrelation method.
 * @param  detector_image   pointer to the  detector image
 * @param  ABSwindow aberation sensor window centers
 * @param  ABSwindowSize window size
 * @param  RefDefocusPositionsMeasured reference positions
 * @param  FocusArray measured focus values
 * @param  ABSErrorStatus errror in the measurement for each telescope
*/
cpl_error_code FocusWithAutoCorrelation(cpl_image *detector_image, 
					int * ABSwindow,
					int  ABSwindowSize,
					double *RefDefocusPositionsMeasured,
					double * FocusArray,
					double *TargetPositions,
					int * ABSErrorStatus){
    cpl_error_code error=CPL_ERROR_NONE;
    int offset=24;

    int ii, i;
    int N=64;
    double locX[3];
    double locY[3];
    double FittedX[3];
    double FittedY[3];
    double FitPositions[6];
    double FitPositions_error[6];
    double x_radius, y_radius;
    double x_shift, y_shift;

    for (ii = 1; ii <= 4; ii++) {/* For each telescope */
    cpl_image *window_abs = cpl_image_extract(detector_image, 
					      ABSwindow[ii-1]+offset, 
					      ABSwindow[ii+3]+offset, 
					      ABSwindow[ii-1]+ABSwindowSize-offset-1, 
					      ABSwindow[ii+3]+ABSwindowSize-offset-1);

 /*********************
	  Hot pixel estimation: START
    ***********************/
    cpl_mask * Mask = NULL;
    int nsigma=10;
    double *image_data=cpl_image_get_data_double(window_abs);
    int npix_size_x = (int)cpl_image_get_size_y(window_abs);
    int npix=(int)cpl_image_get_size_x(window_abs)*(int)cpl_image_get_size_y(window_abs);
    
    cpl_image *bpm=cpl_image_new(npix_size_x, npix_size_x, CPL_TYPE_INT);
    int *bpm_data=cpl_image_get_data_int(bpm);
    
    int nbpm= BadPixDetection(image_data, bpm_data, nsigma, npix_size_x, npix);
    
    cpl_image *bpm_image= cpl_image_wrap_int (npix_size_x, npix_size_x, bpm_data);
    Mask=cpl_mask_threshold_image_create(bpm_image, 0, 2);
    
    if(nbpm < 3000)	error = cpl_image_reject_from_mask(window_abs, Mask);
    error=  cpl_detector_interpolate_rejected(window_abs); 
   
    cpl_mask_delete(Mask);
    cpl_image_unwrap(bpm_image);
    cpl_image_delete(bpm);
    
    /*********************
	  Hot pixel estimation: END
    ***********************/


    locX[0]=N;
    locX[1]=N+RefDefocusPositionsMeasured[0+(ii-1)*4];
    locX[2]=N+RefDefocusPositionsMeasured[1+(ii-1)*4];

    locY[0]=N;
    locY[1]=N+RefDefocusPositionsMeasured[2+(ii-1)*4];
    locY[2]=N+RefDefocusPositionsMeasured[3+(ii-1)*4];

    cpl_image *corr_image= FFTCorr(window_abs, window_abs);
    error=MultipleGaussFit(corr_image, locX, locY, 3, 12, FitPositions, FitPositions_error);

    for(i=0; i<3; i++) {
    FittedX[i]=FitPositions[i];
    FittedY[i]=FitPositions[i+3];
    }
    
    x_radius= sqrt( (FittedX[1]-locX[0])*(FittedX[1]-locX[0]) + /* diff along x */
			   (FittedY[1]-locY[0])*(FittedY[1]-locY[0]) /* diff along y */); 
    
    y_radius= sqrt( (FittedX[2]-locX[0])*(FittedX[2]-locX[0]) + /* diff along x */
			   (FittedY[2]-locY[0])*(FittedY[2]-locY[0]) /* diff along y */); 

    x_shift= x_radius-RefDefocusPositionsMeasured [0+ (ii-1)*4];
    y_shift= y_radius-RefDefocusPositionsMeasured [3+ (ii-1)*4];
    
    FocusArray[ii-1]= (x_shift+y_shift)/2.0;

    TargetPositions[0+(ii-1)*4] = x_radius;
    TargetPositions[1+(ii-1)*4] = 0;
    TargetPositions[2+(ii-1)*4] = 0;
    TargetPositions[3+(ii-1)*4] = y_radius;

    cpl_image_delete(corr_image);    
    cpl_image_delete(window_abs);
    }/*end of if(ii = 1; ii <= 4; ii++) */

   
   
    return error;
}



/**
 * @internal 
 * @brief  Measurment of focus using autocorrelation method.
           This procedure computes reference positions for the autocorrelation peaks
 * @param  detector_image   pointer to the  detector image
 * @param  ABSwindow aberation sensor window centers
 * @param  ABSwindowSize window size
 * @param  FactoryRefPositions  rough positions from database
 * @param  RefPositionsMeasured measured positions 
 * @param  RefPositionsMeasuredError Gauss fit error for the measured positions 
 * @param  ABSErrorStatus errror in the measurement for each telescope
*/
cpl_error_code FocusWithAutoCorrelationCalib(cpl_image *detector_image,
					     int * ABSwindow,
					     int  ABSwindowSize, 
					     double* FactoryRefPositions,
					     double* RefPositionsMeasured,
					     double* RefPositionsMeasuredError,
					     int * ABSErrorStatus){
    cpl_error_code error=CPL_ERROR_NONE;
    int offset=24;
    int Nx =128;
    int ii,i;
    int N=Nx/2;
    double locX[3];
    double locY[3];
    double FittedX[3];
    double FittedY[3];
    double FitPositions[6];
    double FitPositions_error[6];
    double x_radius, y_radius;
    
    for (ii = 1; ii <= 4; ii++) {/* For each telescope */
    cpl_image *window_abs = cpl_image_extract(detector_image, 
					      ABSwindow[ii-1]+offset, 
					      ABSwindow[ii+3]+offset, 
					      ABSwindow[ii-1]+ABSwindowSize-offset-1, 
					      ABSwindow[ii+3]+ABSwindowSize-offset-1);

    
 /*********************
	  Hot pixel estimation: START
    ***********************/
    cpl_mask * Mask = NULL;
    int nsigma=10;
    double *image_data=cpl_image_get_data_double(window_abs);
    int npix_size_x = (int)cpl_image_get_size_y(window_abs);
    int npix=(int)cpl_image_get_size_x(window_abs)*(int)cpl_image_get_size_y(window_abs);
    
    cpl_image *bpm=cpl_image_new(npix_size_x, npix_size_x, CPL_TYPE_INT);
    int *bpm_data=cpl_image_get_data_int(bpm);
    
    int nbpm= BadPixDetection(image_data, bpm_data, nsigma, npix_size_x, npix);
    
    cpl_image *bpm_image= cpl_image_wrap_int (npix_size_x, npix_size_x, bpm_data);
    Mask=cpl_mask_threshold_image_create(bpm_image, 0, 2);
    
    if(nbpm < 3000)	error = cpl_image_reject_from_mask(window_abs, Mask);
    error=  cpl_detector_interpolate_rejected(window_abs); 
    
    
    cpl_mask_delete(Mask);
    cpl_image_unwrap(bpm_image);
    cpl_image_delete(bpm);
    
    /*********************
	  Hot pixel estimation: END
    ***********************/


    locX[0]=N;
    locX[1]=N+FactoryRefPositions[0+(ii-1)*4];
    locX[2]=N+FactoryRefPositions[1+(ii-1)*4];
    
    locY[0]=N;
    locY[1]=N+FactoryRefPositions[2+(ii-1)*4];
    locY[2]=N+FactoryRefPositions[3+(ii-1)*4];

    cpl_image *corr_image= FFTCorr(window_abs, window_abs);
    error=MultipleGaussFit(corr_image, locX, locY, 3, 14, FitPositions, FitPositions_error);
    
    for(i=0; i<3; i++) {
    FittedX[i]=FitPositions[i];
    FittedY[i]=FitPositions[i+3];
    }
    
    x_radius= sqrt( (FittedX[1]-locX[0])*(FittedX[1]-locX[0]) + /* diff along x */
			   (FittedY[1]-locY[0])*(FittedY[1]-locY[0]) /* diff along y */); 
    
    y_radius= sqrt( (FittedX[2]-locX[0])*(FittedX[2]-locX[0]) + /* diff along x */
			   (FittedY[2]-locY[0])*(FittedY[2]-locY[0]) /* diff along y */); 
    
    RefPositionsMeasured[0+(ii-1)*4]=x_radius;
    RefPositionsMeasured[1+(ii-1)*4]=0;
    RefPositionsMeasured[2+(ii-1)*4]=0;
    RefPositionsMeasured[3+(ii-1)*4]=y_radius;
    
    cpl_image_delete(corr_image);
    cpl_image_delete(window_abs);
    }/*end of if(ii = 1; ii <= 4; ii++) */
    return error;
}



/* ----------------------------------------------------------------------*/
/**
 * @internal 
 * @brief  Estimates image shift of extended scene using FFT correlation algorithm
 * @param  Ref reference 16x16 pixel
 * @param  image target 16x16 pixel
 * @param  XYcen out put [x, y]
 * @return status error 
 */
cpl_image *FFTCorr(cpl_image *ref, cpl_image *im) {

	cpl_error_code error = CPL_ERROR_NONE;
	/* cpl_ensure_code(Ref  || image , CPL_ERROR_NULL_INPUT); */
	
	if (ref == NULL) {
        printf("Input argument is NULL \n");
	return NULL;
	}
    
	int im_x = cpl_image_get_size_x(ref);
	int im_y = cpl_image_get_size_y(im);
    
	if (im_x != im_y ) {
		cpl_msg_error(cpl_func,
                      "input image size should be squared: im_x = %d, im_y=%d\n", im_x, im_y);
		return NULL;
	}
	


    /*Images are zero padded to cancel out aliasing errors */
	im_x = 2*im_x;
	im_y = 2*im_y;
	
	cpl_image *Ref=cpl_image_new(im_x, im_y, CPL_TYPE_DOUBLE);
	cpl_image *image= cpl_image_new(im_x, im_y, CPL_TYPE_DOUBLE);

	error = cpl_image_copy (Ref, ref, im_x/2-im_x/4, im_y/2-im_y/4);
	error = cpl_image_copy (image, im, im_x/2-im_x/4, im_y/2-im_y/4);
/* 
   The flow chart of finding cross-correlation
   template= real1 + j complex1
   image =   real2 + j coplex2
   
   here, real1= Ref, real2=image
   
   FFT(template) = real1 + j complex1 = cpl_fft_image(real1, complex1, CPL_FFT_DEFAULT);
   FFT(image) =    real2 + j complex2 = cpl_fft_image(real2, complex2, CPL_FFT_DEFAULT);
   The original data is replaced by the data after fourier transform in the above function

   conj(FFT(template)) = real1 - j complex1
   conj(FFT(template)) * FFT(image) = real1 X real2 + complex1 X complex2 +
   j (real1 X complex2 - real2 X complex1 );
   
   Re_part= real1Xreal2 + complex1Xcomplex2;
   Im_part  = real1Xcomplex2 -real2Xcomplex1;
   
   -1
   FFT  { conj(FFT(template))*FFT(image) } = cpl_fft_image(Re_part, Im_part, CPL_FFT_INVERSE);
   
   Correlation peak = roll2(Re_part);
*/


	cpl_image * complex1 = cpl_image_new(im_x, im_y, CPL_TYPE_DOUBLE);
	cpl_image * complex2;
	error=cpl_image_multiply_scalar(complex1, 0.0);
	complex2 = cpl_image_duplicate(complex1);
	/* cpl_ensure_code(complex2, CPL_ERROR_NULL_INPUT); */
    
	error = cpl_image_fft(Ref, complex1, CPL_FFT_DEFAULT); /* Forward tranform */
	error = cpl_image_fft(image, complex2, CPL_FFT_DEFAULT);/* Forward tranform */
    
	if (error != CPL_ERROR_NONE)
		cpl_msg_error(cpl_func, "Error occured at %s:%d\n", __FILE__, __LINE__);
    
	cpl_image * real1Xreal2 = cpl_image_multiply_create(Ref, image);
	cpl_image * complex1Xcomplex2 = cpl_image_multiply_create(complex1, complex2);
	cpl_image * real1Xcomplex2 = cpl_image_multiply_create(Ref, complex2);
	cpl_image * real2Xcomplex1 = cpl_image_multiply_create(image, complex1);
	cpl_image * Re_part = cpl_image_add_create(real1Xreal2, complex1Xcomplex2);
	cpl_image * Im_part = cpl_image_subtract_create(real1Xcomplex2, real2Xcomplex1);

	/* Inverse tranform */
	error = cpl_image_fft(Re_part, Im_part, CPL_FFT_INVERSE);

	if (error != CPL_ERROR_NONE)
		cpl_msg_error(cpl_func, "Error occured at %s:%d\n", __FILE__, __LINE__);
    
	/*
	 * Roll the image to the center
	 */
	double * data_corr = cpl_image_get_data_double(Re_part);
	error= image_roll_double(data_corr, im_x);
	cpl_image * corr_peak_im = cpl_image_wrap_double(im_x, im_y, data_corr);
	cpl_image *corr_peak = cpl_image_extract (corr_peak_im, 
                                              im_x/2-im_x/4+1, im_y/2-im_y/4+1, 
                                              im_x/2+im_x/4, im_y/2+im_y/4);

	cpl_image *return_image=cpl_image_duplicate(corr_peak);


	/* Free buffer */
	cpl_image_unwrap(corr_peak_im);
	cpl_image_delete(corr_peak);
	cpl_image_delete(Ref);
	cpl_image_delete(image);
	cpl_image_delete(Re_part);
	cpl_image_delete(Im_part);
	cpl_image_delete(real1Xcomplex2);
	cpl_image_delete(real2Xcomplex1);
	cpl_image_delete(real1Xreal2);
	cpl_image_delete(complex1Xcomplex2);
	cpl_image_delete(complex1);
	cpl_image_delete(complex2);
    
	return return_image;
}/* FFTCorrCen */



/**
 * @internal 
 * @brief  Measurment of focus using central spots
 * @param  detector_image   pointer to the  detector image
 * @param  ABSwindow aberation sensor window centers
 * @param  ABSwindowSize window size
 * @param  RefPos Reference centers
 * @param  FocusArray focus values measured
 * @param  ABSErrorStatus errror in the measurement for each telescope
*/
cpl_error_code SpecialFocus(cpl_image *detector_image, 
					int * ABSwindow,
					int  ABSwindowSize,
					double *RefPos,
					double* FocusArray,
					int * ABSErrorStatus){
    cpl_error_code error=CPL_ERROR_NONE;

    

    int ii, i,  n;
    double FitPositions[32];
    double FitPositions_error[32];
    double ShiftsX[16];
    double ShiftsY[16];
    double Defocus4Spots[4];
    double largeShift;
    int N=16;
    double avg_x, avg_y;

/* Horizental 8 spots going through center of SH and 
   Vertical 8 spots going through center of SH  */
    int Loc[16]={34, 35, 36, 37, 39, 40, 41, 42, 3, 11, 20, 29, 47, 56, 65, 73}; 
    double Total_spots_x[N];
    double Total_spots_y[N];
    int P;
    double Diff[N];
    
    
    cpl_mask * Mask = NULL;
    int nsigma=10;    
    for(i=0; i<4; i++) FocusArray[i]=0;

    for (ii = 1; ii <= 4; ii++) {/* For each telescope */
    cpl_image *window_abs = cpl_image_extract(detector_image, 
					      ABSwindow[ii-1], 
					      ABSwindow[ii+3], 
					      ABSwindow[ii-1]+ABSwindowSize-1, 
					      ABSwindow[ii+3]+ABSwindowSize-1);
    for(i=0; i<N; i++){
    P=Loc[i];
    Total_spots_x[i]=RefPos[P+(ii-1)*2*ABSNUMSPOTS];
    Total_spots_y[i]=RefPos[77+P+(ii-1)*2*ABSNUMSPOTS];
    }

        /*********************
	  Hot pixel estimation: START
	***********************/
	double *image_data=cpl_image_get_data_double(window_abs);
	int npix_size_x = (int)cpl_image_get_size_y(window_abs);
	int npix=(int)cpl_image_get_size_x(window_abs)*(int)cpl_image_get_size_y(window_abs);
	
	cpl_image *bpm=cpl_image_new(npix_size_x, npix_size_x, CPL_TYPE_INT);
	int *bpm_data=cpl_image_get_data_int(bpm);

	int nbpm= BadPixDetection(image_data, bpm_data, nsigma, npix_size_x, npix);

	cpl_image *bpm_image= cpl_image_wrap_int (npix_size_x, npix_size_x, bpm_data);
	Mask=cpl_mask_threshold_image_create(bpm_image, 0, 2);

	if(nbpm < 3000)	error = cpl_image_reject_from_mask(window_abs, Mask);
	/* error=  cpl_detector_interpolate_rejected(window_abs); */
         

	cpl_mask_delete(Mask);
	cpl_image_unwrap(bpm_image);
	cpl_image_delete(bpm);
	
      	/*********************
	  Hot pixel estimation: END
	***********************/

    error=MultipleGaussFit(window_abs, Total_spots_x, Total_spots_y, 
			   N, 10, FitPositions, FitPositions_error);

    avg_x=0;
    avg_y=0;

    for(i=0; i<N; i++) {
    ShiftsX[i]=FitPositions[i]-Total_spots_x[i];
    ShiftsY[i]=FitPositions[i+16]-Total_spots_y[i];

    avg_x += ShiftsX[i];
    avg_y += ShiftsY[i];
    }
   
    avg_x=avg_x/((double)N);
    avg_y=avg_y/((double)N);
    
    if(fabs(ShiftsX[0]) > 3.0) ShiftsX[0] = avg_x;
    if(fabs(ShiftsX[7]) > 3.0) ShiftsX[7] = avg_x;
    if(fabs(ShiftsY[8]) > 3.0) ShiftsY[8] = avg_y;
    if(fabs(ShiftsY[15]) >3.0) ShiftsY[15] = avg_y;

    largeShift=(fabs(ShiftsX[0])+fabs(ShiftsX[7])+ fabs(ShiftsY[8])+fabs(ShiftsY[15]))/4.0+1.0;

    for(i=0; i<4; i++){
    Diff[i]=ShiftsX[4+i]+ShiftsX[3-i];
    Diff[i+4]=ShiftsX[12+i]+ShiftsX[11-i];
    
    Diff[8+i]=ShiftsY[4+i]+ShiftsY[3-i];
    Diff[12+i]=ShiftsY[12+i]+ShiftsY[11-i];
    }
    

   
    n=0;
    for(i=0; i<4; i++){
    if(fabs(Diff[i]) > largeShift){
    n++;
    ShiftsX[3-i] = (ShiftsX[3-i]-ShiftsX[i+4]+ShiftsY[11-i]-ShiftsY[12+i])/4.0;
    ShiftsX[4+i] = (ShiftsX[i+4]-ShiftsX[3-i]-ShiftsY[11-i]+ShiftsY[12+i])/4.0;
    }

    if(fabs(Diff[i+4]) > largeShift){
    ShiftsX[12+i] = 0.0;
    ShiftsX[11-i] = 0.0;
    }

    if(fabs(Diff[i+8]) > largeShift){
    ShiftsY[4+i] = 0.0;
    ShiftsY[3-i] = 0.0;
    }

    if(fabs(Diff[i+12]) > largeShift){
    n++;
    ShiftsY[11-i] = (ShiftsX[3-i]-ShiftsX[i+4]+ShiftsY[11-i]-ShiftsY[12+i])/4.0;
    ShiftsY[12+i] = (ShiftsX[i+4]-ShiftsX[3-i]-ShiftsY[11-i]+ShiftsY[12+i])/4.0;
    }
    }

    for(i=0; i<4; i++) Defocus4Spots[i]=0;
    for(i=0; i<4; i++){ 
    Defocus4Spots[i] +=-ShiftsX[3-i];
    Defocus4Spots[i] +=ShiftsX[4+i];
    Defocus4Spots[i] +=-ShiftsY[11-i];
    Defocus4Spots[i] +=ShiftsY[12+i];
    }

    double LookUpTable[4]={0.372931, 0.743628, 1.069492, 1.445504};

    for(i=0; i<4; i++) {
    Defocus4Spots[i]=Defocus4Spots[i]/4.0;

    }

    for(i=1; i<3; i++) {
    FocusArray[ii-1] +=Defocus4Spots[i]/LookUpTable[i];
    }    
    FocusArray[ii-1]=FocusArray[ii-1]/2.0;

    cpl_image_delete(window_abs);
    }/*end of if(ii = 1; ii <= 4; ii++) */

    return error;
}

/**
 * @internal 
 * @brief  Filter out bad pixels using input sigma (pixels)
 * @param  ima input image
 * @param  bpm bad pixel preallocated array
 * @param  nsigma number of rms about the local mean out of which is
 * @param  image_size image size
 * @param  npix  total number of pixels in the image
 * @return status error 
 *
 *  Filter out the pixels that deviate from the local statistics.
 *  The mean and rms of the 8 (the minimum and maximum of these
 *  8 neighbors are excluded in the mean and rms computation) is
 *  computed. All pixels that deviates more than "nsigma" rms
 *  from the mean are flagged as bad pixels. The image and newly
 *  created bad pixel map are passed to the routine dpc()
 *  for correction. The processus can be iterated.
 *
 *  Taken from yorick-yao
*/

int BadPixDetection(double *ima, int *bpm, double nsigma, int stride, int npix)
{
    int i,k,n, nbpm;
    double avg1,avg2,val,std; 
    nbpm=0;

    for (i=0;i<npix;i++) {
    avg1 = avg2 = 0.0;
    n = 0;
    /* compute the avg and rms in a 3x3 box centered on pixel */
    /* k = i; avg1=ima[k]; avg2=ima[k]*ima[k]; n++; */
    
    k = i-1;        if ((k>=0)&&(k<npix)) { avg1+=ima[k]; avg2+=ima[k]*ima[k]; n++; }
    k = i+1;        if ((k>=0)&&(k<npix)) { avg1+=ima[k]; avg2+=ima[k]*ima[k]; n++; }
    k = i-stride;   if ((k>=0)&&(k<npix)) { avg1+=ima[k]; avg2+=ima[k]*ima[k]; n++; }
    k = i-stride-1; if ((k>=0)&&(k<npix)) { avg1+=ima[k]; avg2+=ima[k]*ima[k]; n++; }
    k = i-stride+1; if ((k>=0)&&(k<npix)) { avg1+=ima[k]; avg2+=ima[k]*ima[k]; n++; }
    k = i+stride;   if ((k>=0)&&(k<npix)) { avg1+=ima[k]; avg2+=ima[k]*ima[k]; n++; }
    k = i+stride-1; if ((k>=0)&&(k<npix)) { avg1+=ima[k]; avg2+=ima[k]*ima[k]; n++; }
    k = i+stride+1; if ((k>=0)&&(k<npix)) { avg1+=ima[k]; avg2+=ima[k]*ima[k]; n++; }
    
    avg1 = avg1 / (double)n;
    avg2 = avg2 / (double)n;
    std  = avg2 - avg1 * avg1;
    if (std<0.) std=0.0;
    else std = sqrt(std);
    val = fabs(ima[i]-avg1);
    if (val>(nsigma*std)){
    bpm[i]=1;
    nbpm++;
    } else bpm[i]=0;
    }
    
    return nbpm;
}
/*___oOo___*/

