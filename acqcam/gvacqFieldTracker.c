 /************************************************************************
 *   This file is part of the E.S.O. - VLTI GRAVITY project
 *
 *   
 *   who       when        what
 *   --------  ----------  ----------------------------------------------
 *   aa       2017-01-05  Introduced the authomatic detection of saturated stars. If flagued the FIT is modifyed
 *   aa       2016-07-15  On UTs the AO makes stars images less isotropic increased fabs(FWHM[0]-FWHM[1])<2.0
 *   aa       2016-07-10  Picking the second most luminous star if nearer to prev  
 *   na,aa,fe 2016-05-13  Debugged ADR code and made it work
 *   narsi    2016-05-13  Converting adr angle to pixel coordinates on the gvrtdACQ using new ZenitheAngle
 *   narsi    2016-05-12  Considering North and Zenith axes for four telescopes instead of one
 *   narsi    2016-04-13  Finding FWHM of defocussed star with large fit window size
 *   narsi    2016-03-28  Checking the ADR results for wrong measurements
 *   feisenha 2015-10-31  changed sign of paralactic angle and subtracted 90 degrees
 *                        for the moment ad hoc adjusted to get sign and direction of ADR right
 *   narsi    2015-09-20  spilt the gvacqProcessImage.c functions to gvacqFieldTracker.c
 *   narsi    2015-09-20  created
 */

/****************************************************************************
 *   NAME 
 *   gvacqFieldTracker - procedere to compute field guiding
 *  
 *   SYNOPSIS
 *   #include <gvoProcessImageAC.h>
 * 
 *   DESCRIPTION
 *   Provides the position of bright object in the field imager pixel coordinate system.
 *
 *   This function takes the pointer to the detector image and determines the
 *   the brightest object in the field in pixel coordinates.
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
#define PI 3.14159265359
#define USEPREVSTAR 1

/* Defaults are over written by the database */
static int SBA[8]={ 351, 821, 1250, 1729, 279, 261, 235, 295 }; /* FI window ceneters */
static double FIWindowAndSigma[2]={110, 2.0}; /* Field image window size/2, sigma */
static int FIwindowGauss=12; /*field image gauss fit window */
static int PRINTMSG=0;
#ifdef USEPREVSTAR
static double prevObjectPosition[8]={-1,-1,-1,-1,-1,-1,-1,-1};
#endif
static int FIStar2fromStar1[8]={0,0,0,0,0,0,0,0}; /* Star1 to Star2 delta seed coordinates x1 x2 x3 x4 y1 y2 y3 y4  */


/*-----------------------------------------------------------------------------
 Private function prototypes
 ----------------------------------------------------------------------------*/
cpl_error_code FIFindPosition(cpl_image *FI_image, double FIfwhm_minthreshold, int dx, int dy, double *StarPos, double *StarPos2, double *StarPosError, double *StarFWHM, double *StarFlux , double FISigma, int reject);
cpl_error_code FIFindFWHM(cpl_image *FI_image, double px, double py, double FIwindowGauss, double* FWHM);

/* Refractive index measurement algorithm */
double RefractiveIndex_Owens(double lambda, double P, double t, double Pw);
double RefractiveIndex_Schubert(double lambda, double p, double t, double p_w);

/* Estimate adr value in arcsec */
double ADR(double n, double nRef, double Zenith);
double partialPressure(double rh, double t);
double ADRshift(double lambda, double lambdaRef, double P, double t, double rh, double Zenith);
double polynomial_1d(double *p, double x);
double atanFit(double *a, double x);
cpl_error_code EffectiveWavelength(double H_K, double *ACQPolyCoeff, double *BCPolyCoeff, double *EffectiveLambda);
cpl_error_code ADRCorrectionOffset(double H_K, double *ACQPolyCoeff, double *BCPolyCoeff, 
				   double P, double t, double humidity, double Zenith, 
				   double * EffectiveLambda, double *ObjCorr);



/*
 *******************************************************************************
 *                                                                             *
 * Field imager                                                                *
 *                                                                             *
 *******************************************************************************
 */

/**
 * @internal
 * @brief  Provides the position of bright object in the pixel coordinate system.
 * @param  DetPointer   pointer to the detector image
 * @param  ObjectPosition  (output) The brightest object position (in pixels)
 * @param  ObjectPositionError (output)  The brightest object position error (in pixels)
 * @return status error  0: OK, else: Bad
 * @see   setCoordinates, setFIfitWindowSizeAndSigma
 *
 * This function takes the pointer to the detector image and determines the
 * the brightest object in the field in pixel coordinates.
 *
 * @note
 * 1. check the inputs
 * 2. find the brightest object position in the field by using
 *    cpl_apertures_extract by thresholding with sigma. Next sort the objects by flux.
 * 3. Fit the Gaussian (with normal weights now; not fixed the ADR elongated spots)
 *    around the positions obtained in step (2) with window size equal to 16px;
 * 4. Return the object pixel positions and corresponding errors
 *
 * Note: The input image buffer should be background cleaned
 * 
 * Function pre-required
 *  func setCoordinates sets field imager windows to global variable
 *  func setFIfitWindowSizeAndSigma sets field imager widnow size and 
 *       sigma in pixels to scan stars in the field
 *
 * @par Error Handling:
 * 
 * Return 0 on successful fitting
 *          -> The object in return after gaussfit
 * Return -1 for image is NULL and database reading problem
 *          -> The database object position is returned 
 * Return -2 if the no objects flux is found
 *          ->The database object position is returned
 * Return -3 if the object is near the edges of the window chosen
 *          -> The object which is scanned by the first algorithm is returned
 * 
 */
cpl_error_code gvacqFieldImager(cpl_image *DetPointer,
				double FIfwhm_minthreshold,
				 double * ObjectPosition, 
				 double * ObjectPosition2, 
				double * ObjectPositionError,
				double * FiFWHM,
				double * FiFlux,
				int * FIErrorStatus)
{


	/* initialising variables */
	cpl_error_code error = CPL_ERROR_NONE;
	int i, j;
	int SBAlow[8], SBAhigh[8];
	int Sub_window = (int) FIWindowAndSigma[0];
	double FISigma = FIWindowAndSigma[1];
	double StarPos[2];
	double StarPos2[2];
	double StarPosError[2];
	double StarFWHM;
	double StarFlux;
	int reject=1;
    
    /*test entries */
	if (DetPointer ==NULL || ObjectPosition == NULL  || ObjectPosition2 == NULL || ObjectPositionError == NULL ){ 
        printf("ERROR::The Det pointer or input parameters are NULL. Occured at %s:%s:%d\n", 
               __FILE__, __FUNCTION__,__LINE__);
        
        /* The database object position is returned  */
        for(i=0; i<4; i++){
		ObjectPosition[i] = -1;
		ObjectPosition[i + 4] = -1;
		ObjectPosition2[i] = -1;
		ObjectPosition2[i + 4] = -1;
		ObjectPositionError[i] = -1;
		ObjectPositionError[i + 4] = -1;
        }
        for(j=0; j<4; j++) FIErrorStatus[i]=-1;
        return -1;
	}
    
	for (i = 0; i < 8; i++) {
        SBAlow[i] = SBA[i] - Sub_window;
        SBAhigh[i] = SBA[i] + Sub_window-1;
	}
	
	for (i = 0; i < 4; i++) {/* tel */
        cpl_image *FI_image = cpl_image_extract(DetPointer, SBAlow[0 + i], SBAlow[4+ i], 
                                                SBAhigh[0 + i], SBAhigh[4 + i]);
#ifdef USEPREVSTAR
        /* feed prevObjPos into StartPos */
        if  (prevObjectPosition[i]==-1) {
           StarPos[0]=-1; 
           StarPos[1]=-1; 
        } else {
           StarPos[0]=prevObjectPosition[i]-(double)SBAlow[i];
           StarPos[1]=prevObjectPosition[i+4]-(double)SBAlow[i+4];
        }
        /* end  prevObjPos */
#endif

        error =  FIFindPosition(FI_image, FIfwhm_minthreshold, FIStar2fromStar1[i],  FIStar2fromStar1[i+4],  StarPos, StarPos2, StarPosError , 
				&StarFWHM, &StarFlux, FISigma, reject);
        cpl_image_delete(FI_image);
        
	FIErrorStatus[i]=error;
	

	if( fabs(StarPos[0]) >= 2*Sub_window  || fabs(StarPos[1]) >=  2*Sub_window){
	FIErrorStatus[i]=-2;
	}

        
        if(error !=0){
            if(error==-1 || error ==-2 || error ==-3){
                if(error==-1){
                    printf("Tel.%d: The Det pointer or input parameters are NULL. \n", i+1);
                }else if(error==-2){
                    printf("Tel.%d: No stars are found in the image. \n", i+1);
                }else printf( "Tel.%d: scanned for object data but not found [with condition: \
fwhm(object) > FIfwhm_minthreshold (= %f px)]  \n", 
                              i + 1, FIfwhm_minthreshold);
            }else {
                if(error==-4){
                    printf("Tel.%d: Bright star is found near the Gaussian window. Not able to fit. \n", i+1);
                }else printf("Tel.%d: The Gaussian fit is not successful. \
Results are not very accurate. \n", i+1);

	    }
	    
	    ObjectPosition[i] = -1;
	    ObjectPosition[i + 4] = -1;
	    ObjectPosition2[i] = -1;
	    ObjectPosition2[i + 4] = -1;
	    ObjectPositionError[i] = -1;
	    ObjectPositionError[i + 4] = -1;
            FiFWHM[i]=-1;
            FiFlux[i]=-1;
            
            
        }else {/* end of error handling */     

            ObjectPosition[i] = StarPos[0]+(double)SBAlow[i];
            ObjectPosition[i + 4] = StarPos[1]+(double)SBAlow[i+4];
            ObjectPosition2[i] = StarPos2[0]+(double)SBAlow[i];
            ObjectPosition2[i + 4] = StarPos2[1]+(double)SBAlow[i+4];
            ObjectPositionError[i] = StarPosError[0];
            ObjectPositionError[i + 4] = StarPosError[1];
            FiFWHM[i]=StarFWHM;
            FiFlux[i]=StarFlux;

	}

#ifdef USEPREVSTAR
        prevObjectPosition[i] =  ObjectPosition[i];
        prevObjectPosition[i+4] =  ObjectPosition[i+4];
#endif

	}/* End of i telescope */	
	return 0;
}

 /* End of Field Imager function */

/**
* @internal
 * @brief  Finds star position in  field imager window of a telescope
 * @param  FI_image field imager widnow of a telescope
 * @param  StarPos (output) star position(x,y)
 * @param  StarPos2 (output) second star position(x,y)
 * @param  StarPosError(output) star position error (x,y)
 * @param  FISigma sigma is used to scan the stars in the window
 * @return reject 1, it rejects saturated pixels in stars scanning, if 0 all pixels are considered 
 * @see    gvacqFieldImager
 */
cpl_error_code FIFindPosition(cpl_image *FI_image, double FIfwhm_minthreshold, int dx, int dy,
			      double *StarPos, double *StarPos2, double *StarPosError , double *StarFWHM, 
			      double *StarFlux, double FISigma, int reject)
{
    cpl_error_code error=CPL_ERROR_NONE;
#ifdef VLT2014
    cpl_size px, py, px2, py2,buf; 
#else  	 
    double px, py, px2, py2,buf;
#endif
    
    int RON=12;
    double max_center_error = 15.0; /*pixels*/
    double max_sigma_error = 30.0; /*pixels*/
    int robustness = 15; /*maximum # of fitting tries */
    cpl_mask * Mask = NULL;
    cpl_image * Mask_im = NULL;
    
/* ------ STEP 1: Check if the image is NULL input------------ */
    if (FI_image == NULL){
        return -1;
    }
    
    int im_size_x = cpl_image_get_size_x(FI_image);
    int im_size_y = cpl_image_get_size_y(FI_image);
    
/* ------ STEP 2: Check if the image has flux and stars ------------ */
    
    if(im_size_x*RON > cpl_image_get_absflux(FI_image) ) {
        return -2;
    }
    
    
/* ------ STEP 3: Check if the image has stars or only hot pixels------------ */
    error = Check_objects_found(FI_image, FIfwhm_minthreshold);
    if(error == -1 || error ==-2 ){
        return -3;
    }
    
/* ------ STEP 4: finally scan the stars in the image using sigma as threshold------------ */
    /*The threshold used for the detection is the median plus the avg distance to 
      the medeian times sigma*/
    cpl_apertures *objects = cpl_apertures_extract_sigma(FI_image, FISigma); 
    if(objects==NULL) {  
        return -2;
    }
    
    /* int Nobjects = cpl_apertures_get_size(objects);	 */
    error = cpl_apertures_sort_by_flux(objects); /* Sort the objects by flux */
    if (error != CPL_ERROR_NONE)
        cpl_msg_error(cpl_func, "Error occured at %s:%d\n", __FILE__, __LINE__);
    
    
    /* find the first moment of the brightest object (pixels)*/
    px = (double) cpl_apertures_get_centroid_x(objects, 1);
    py = (double) cpl_apertures_get_centroid_y(objects, 1);
    if (cpl_apertures_get_size(objects) >=2 ) { 

       /* default */
       px2 = (double) cpl_apertures_get_centroid_x(objects, 2);
       py2 = (double) cpl_apertures_get_centroid_y(objects, 2);

       /* if the dx, dy from star1 to star2 is defined, then look for aperture closer to (px+dx,py+dy) */
       if ( dx!=0 || dy!=0) { 
          int naper;
          double dist,px2temp,py2temp;
          dist==sqrt((px2-px-dx)*(px2-px-dx)+(py2-py-dy)*(py2-py-dy));
          for (naper=3; naper < cpl_apertures_get_size(objects); naper++) {
                px2temp = (double) cpl_apertures_get_centroid_x(objects, naper );
                py2temp = (double) cpl_apertures_get_centroid_y(objects, naper );
               if (sqrt((px2temp-px-dx)*(px2temp-px-dx)+(py2temp-py-dy)*(py2temp-py-dy))<dist) {
                   px2 = px2temp;
                   py2 = py2temp;
                   dist==sqrt((px2-px-dx)*(px2-px-dx)+(py2-py-dy)*(py2-py-dy));
               }
          }

       }
       /* printf("prev %f %f cent1 %f %f flux %f cent2 %f %f flux %f\n",StarPos[0]+1,StarPos[1]+1,
           px,py,cpl_apertures_get_flux(objects,1),px2,py2,cpl_apertures_get_flux(objects,2)); */
    } else {
       px2 = 0;
       py2 = 0;
    }

#ifdef USEPREVSTAR
    /* test is the second less luminous is closer to prev */
    if ( StarPos[0] != -1  &&  /* previous position already defined */
       cpl_apertures_get_size(objects) >=2  && /* a second aperture exists */ 
       cpl_apertures_get_flux(objects,2) >= 0.4*cpl_apertures_get_flux(objects,1)  && /* flux2>=0.4*flux1 */
       (StarPos[0]+1-px2)*(StarPos[0]+1-px2) +  (StarPos[1]+1-py2)*(StarPos[1]+1-py2) <  
       (StarPos[0]+1-px)*(StarPos[0]+1-px) + (StarPos[1]+1-py)*(StarPos[1]+1-py) /* second is closer to prev */ 
       ) { 
          buf= px;
          px = px2;
          px2= buf;
          buf= py;
          py = py2;
          py2= buf;
          printf("Choosing second brightest star\n");
    } 
    /* end dealing with prev */
#endif
    
    /* ------ STEP 5: reject/fix pixels if the found object is saturated star ------------ */    
    
     /* If the image is saturated
        reject==1->> The saturated pixels are rejected for the fitting
        else ->> The saturated pixels can be interpolated using 8 nearest neighbours
     */
    
    if(reject==1){
        Mask=cpl_mask_threshold_image_create(FI_image, -0.9999999, 0.000001);
        error = cpl_image_reject_from_mask (FI_image, Mask);
    }else{	
        Mask_im = cpl_image_new_from_mask (Mask);
        error=     cpl_detector_interpolate_rejected(FI_image);		
        cpl_image_delete(Mask_im);
    }
    
    /* ------ STEP 6: Check if the star found near the edge of the window given ------------ */
    /*  Checking if the object is near to the edges, before gaussian fit */
    if ((px - FIwindowGauss/ 2) <= 2 || (px + FIwindowGauss / 2) >= im_size_x - 2 || 
        (py- FIwindowGauss / 2) <= 2 || (py + FIwindowGauss / 2) >= im_size_y - 2) {
        
        StarPos[0] = px-1;
        StarPos[1] = py-1;
        StarPos2[0] = px2-1;
        StarPos2[1] = py2-1;
        StarPosError[0] = -1;
        StarPosError[1] = -1;
        
        cpl_apertures_delete(objects);
        cpl_mask_delete(Mask);
        return -4;
    } else {/* Gauss Fit*/

        /* ------ STEP 7: Create radial profile  to define saturated flag and mark pixels as bad -------------- */ 

#define PROFILE_MAXLEN 51
        double content,profile[PROFILE_MAXLEN];
        int i,j,r,pis_rejected,Nprofile[PROFILE_MAXLEN];
        int FIwindowGaussFactor=1;

        /* set profile to 0 */
        for (r=0; r<PROFILE_MAXLEN; r++) {
            profile[r]=0;
            Nprofile[r]=0;
        }

        /* loop over window and commulate profile around (px,py) */
        for (i=-FIwindowGauss ; i <= FIwindowGauss ; i++ ) 
          for (j=-FIwindowGauss ; j <= FIwindowGauss ; j++) {
              content=cpl_image_get(FI_image,px+i,py+j,&pis_rejected);
              r=(int)(sqrt(i*i+j*j)+0.5);

              if (!pis_rejected) {
                 profile[r]+=content;
                 ++Nprofile[r];
          }
        }

        /* normalize profile average on circular disks */
        for (r=0; r<FIwindowGauss+1; r++) {
           profile[r]=profile[r]/Nprofile[r];
        }

        /* find the profile maximum */
        for (r=0; r<FIwindowGauss && profile[r+1]>profile[r] ; r++) {}

        /* Flag saturated if large counts and either max more than 2 pixels away from the center or *
         * profile is not decreasing by more than 30% from the center to two pixels away            */

        if (profile[r]>10000 && ( r>1 || (profile[2]-profile[0])/profile[r]>-0.2 ) ) {

          if (r>1) {
          /* optimizing a bit px, py  by looking first neigbors for a lower value */
          double pxstore=px;
          double pystore=py;
          content=cpl_image_get(FI_image,px,py,&pis_rejected);
          for (i=-1; i<=1; i++) 
              for (j=-1; j<=1; j++) 
                   if ( cpl_image_get(FI_image,px+i,py+j,&pis_rejected) < content && !pis_rejected ) { 
                         content=cpl_image_get(FI_image,px+i,py+j,&pis_rejected);
                         pxstore=px+i;
                         pystore=py+j;
                   } 
          px=pxstore;
          py=pystore;
          }

          /* in clipm_centroiding_gauss() that is used in clipm_centroiding_multi_gauss()
           * @par Bad Pixel Handling:
           * Bad pixel maps are supported. This means that bad pixels are omitted
           * during computation of the marginal distributions, this means during
           * averaging.
           */

           for (i=-FIwindowGauss ; i <= FIwindowGauss ; i++ ) 
             for (j=-FIwindowGauss ; j <= FIwindowGauss ; j++) {
                 if (r > (int)(sqrt(i*i+j*j)+0.5) ) error = cpl_image_reject(FI_image, px+i, py+j );
             }

           /* also increase the window size if within boundaries */
          if ( (px - FIwindowGauss) > 2 && (px + FIwindowGauss ) < im_size_x - 2  && 
               (py - FIwindowGauss) > 2 && (py + FIwindowGauss ) < im_size_y - 2) 
                 FIwindowGaussFactor=2; 

           /* do not allow the fit window to decrease by too much */
           robustness=0;

        } /* end saturated changes */

        /* ------ STEP 7.1: Fit a Gaussian on the above found star ------------ */
        
        cpl_matrix * Objectlocations = cpl_matrix_new(2, 1); /* 1 spot (x and y) */
        cpl_matrix * xy_centre = cpl_matrix_new(2, 1); /* 1 spot (x and y) */
        cpl_matrix * xy_centre_err = cpl_matrix_new(2, 1);
        cpl_matrix * xy_sigma = cpl_matrix_new(2, 1);
        cpl_matrix * xy_sigma_err = cpl_matrix_new(2, 1);
        cpl_matrix * xy_fwhm = cpl_matrix_new(2, 1);
        cpl_matrix * xy_fwhm_err = cpl_matrix_new(2, 1);
        cpl_array * all_error_codes = cpl_array_new(1, CPL_TYPE_INT);
        
        cpl_matrix_set(Objectlocations, 0, 0, px); /* Guess spot position to fit Guassin */
        cpl_matrix_set(Objectlocations, 1, 0, py); 
        
        error = clipm_centroiding_multi_gauss(FI_image, 
                                              Objectlocations, 
                                              FIwindowGauss*FIwindowGaussFactor, 
                                              &xy_centre, 
                                              &xy_centre_err, 
                                              &xy_sigma, 
                                              &xy_sigma_err,
                                              &xy_fwhm, 
                                              &xy_fwhm_err, 
                                              NULL, 
                                              &all_error_codes,
                                              robustness);
        
        /* reset to unsaturated condition */
        FIwindowGaussFactor=1; 
        robustness=15;
        
        double xerr = cpl_matrix_get(xy_centre_err, 0, 0); /* x center error */
        double yerr = cpl_matrix_get(xy_centre_err, 1, 0); /* y center error */
        double sigmax = cpl_matrix_get(xy_sigma, 0, 0); /* spot sigma */
        double sigmay = cpl_matrix_get(xy_sigma, 1, 0);
        double SecondFWHM[2];
	double FWHM[2];
        /*
         * Error handling for failing fit
         * The brightets object centers and  fitting errors returning
         */
        if (error == CPL_ERROR_NONE && xerr < max_center_error
            && yerr < max_center_error && sigmax
            < max_sigma_error && sigmay < max_sigma_error) {
            /* The brightets object centers and  fitting errors returning*/
            
            StarPos[0] = cpl_matrix_get(xy_centre, 0, 0) - 1;
            StarPos[1] = cpl_matrix_get(xy_centre, 1, 0) - 1;
            
            StarPosError[0] = cpl_matrix_get(xy_centre_err, 0, 0);
            StarPosError[1] = cpl_matrix_get(xy_centre_err, 1, 0);

	    *StarFlux=cpl_apertures_get_flux(objects,1); 

	    FWHM[0]=cpl_matrix_get(xy_fwhm,0,0);
	    FWHM[1]=cpl_matrix_get(xy_fwhm,1,0);

	    if(FWHM[0] !=-1 && FWHM[1] !=-1 && FWHM[0] >= 2.0 && FWHM[1] >= 2.0 
	       && fabs( FWHM[0]- FWHM[1]) < 2.0){
	    *StarFWHM=(FWHM[0]+FWHM[1])/2.0;
	    }else {
	      FIFindFWHM(FI_image, StarPos[0], StarPos[1], FIwindowGauss*2, SecondFWHM);
	    if(SecondFWHM[0] !=-1 && SecondFWHM[1] !=-1 && SecondFWHM[0] >= 2.0 
	       && SecondFWHM[1] >= 2.0 && fabs( SecondFWHM[0]- SecondFWHM[1]) <=10){
	    *StarFWHM=(SecondFWHM[0]+SecondFWHM[1])/2.0;
	    }else *StarFWHM=-1;
	    }
     
            /* Second star fit */
            if (px2 == 0 || py2 ==0 || (px2 - FIwindowGauss/ 2) <= 2 || (px2 + FIwindowGauss / 2) >= im_size_x - 2 || 
                           (py2 - FIwindowGauss/ 2) <= 2 || (py2 + FIwindowGauss / 2) >= im_size_y - 2) {
                  StarPos2[0] = px2-1;
                  StarPos2[1] = py2-1;

            } else {

                  cpl_matrix_set(Objectlocations, 0, 0, px2); 
                  cpl_matrix_set(Objectlocations, 1, 0, py2); 

                  error = clipm_centroiding_multi_gauss(FI_image, 
                                              Objectlocations, 
                                              FIwindowGauss, 
                                              &xy_centre, 
                                              &xy_centre_err, 
                                              &xy_sigma, 
                                              &xy_sigma_err,
                                              &xy_fwhm, 
                                              &xy_fwhm_err, 
                                              NULL, 
                                              &all_error_codes,
                                              robustness);

                  if (error == CPL_ERROR_NONE ) {
                       StarPos2[0] = cpl_matrix_get(xy_centre, 0, 0) - 1;
                       StarPos2[1] = cpl_matrix_get(xy_centre, 1, 0) - 1;
                  } else {
                       StarPos2[0] = px2-1;
                       StarPos2[1] = py2-1;
                  }
            } /* end of second star fit */

        } else {

            StarPos[0] = px-1;
            StarPos[1] = py-1;
            StarPos2[0] = px2-1;
            StarPos2[1] = py2-1;
            StarPosError[0] = -1;
            StarPosError[1] = -1;
            *StarFWHM=-1;
	    *StarFlux=-1;

            /*free the buffers*/
            cpl_apertures_delete(objects);
            cpl_mask_delete(Mask);
            cpl_matrix_delete(Objectlocations);
            cpl_matrix_delete(xy_centre);
            cpl_matrix_delete(xy_centre_err);
            cpl_matrix_delete(xy_sigma);
            cpl_matrix_delete(xy_sigma_err);
            cpl_matrix_delete(xy_fwhm);
            cpl_matrix_delete(xy_fwhm_err);
            cpl_array_delete(all_error_codes);
            return -5;
        }
        
        /*free the buffers*/
        cpl_apertures_delete(objects);
        cpl_mask_delete(Mask);
        cpl_matrix_delete(Objectlocations);
        cpl_matrix_delete(xy_centre);
        cpl_matrix_delete(xy_centre_err);
        cpl_matrix_delete(xy_sigma);
        cpl_matrix_delete(xy_sigma_err);
        cpl_matrix_delete(xy_fwhm);
        cpl_matrix_delete(xy_fwhm_err);
        cpl_array_delete(all_error_codes);
    }/* End of else Gaussian Fit */
    
    return error;
}/* End of FIFindPosition */


cpl_error_code FIFindFWHM(cpl_image *FI_image, double px, double py, 
			  double FIwindowGauss, double* FWHM){
    
    int robustness=9;
    cpl_error_code error=CPL_ERROR_NONE;
    cpl_matrix * Objectlocations = cpl_matrix_new(2, 1); /* 1 spot (x and y) */
    cpl_matrix * xy_centre = cpl_matrix_new(2, 1); /* 1 spot (x and y) */
    cpl_matrix * xy_centre_err = cpl_matrix_new(2, 1);
    cpl_matrix * xy_fwhm = cpl_matrix_new(2, 1);
    cpl_matrix * xy_fwhm_err = cpl_matrix_new(2, 1);
    cpl_array * all_error_codes = cpl_array_new(1, CPL_TYPE_INT);
    
    cpl_matrix_set(Objectlocations, 0, 0, px); /* Guess spot position to fit Guassin */
    cpl_matrix_set(Objectlocations, 1, 0, py); 
   
    error = clipm_centroiding_multi_gauss(FI_image, 
					  Objectlocations, 
					  FIwindowGauss, 
					  &xy_centre, 
					  &xy_centre_err, 
					  NULL, 
					  NULL,
					  &xy_fwhm, 
					  &xy_fwhm_err, 
					  NULL, 
					  &all_error_codes,
					  robustness);
          
    FWHM[0]=cpl_matrix_get(xy_fwhm,0,0);
    FWHM[1]=cpl_matrix_get(xy_fwhm,1,0);
    
    cpl_matrix_delete(Objectlocations);
    cpl_matrix_delete(xy_centre);
    cpl_matrix_delete(xy_centre_err);
    cpl_matrix_delete(xy_fwhm);
    cpl_matrix_delete(xy_fwhm_err);
    cpl_array_delete(all_error_codes);
    
    return error;    
}


/* Atmospheric differential refraction (ADR) code */
/*
 * @internal
 * @brief computes Refractive index using lambda, pressure (P), temperature (T) 
 *        and partial pressure (Pw)
 * 
 * Formula used: Owens, App. optics 1967, Eq. 32. http://dx.doi.org/10.1364/AO.6.000051
 */
double RefractiveIndex_Owens(double lambda, double P, double t, double Pw)
{
    
    if(lambda > 2.5 || lambda < 1.0 ) {
        printf("Warning: lambda should be within 1.0-2.5 \n");
        return -999;
    }
    
    if(P<0 || t < 0 || Pw < 0 ) {
        printf("Pressure, Temperature and Partial pressure needs to be positive \n");
        return -999;
    } 
    
    double p,  f, a1, a2, a3, a4, a5, a6, a7;
    double sigma=1.0/lambda;
    
    p = 0.750061683*(P); /* convert to millibars to toricelli */
    f=0.750061683*Pw;/* convert to millibars to toricelli */
    t=t-273.15; /* convert K to Centigrade */
    
    a1=8342.13;
    a2=2406030/(130-sigma*sigma);
    a3=15997/(38.9-sigma*sigma);
    a4=p/720.775;
    a5=1+p*(0.817 - 0.0133*t)*1e-6;
    a6=1+0.0036610*t;
    a7=f*(5.722-0.0457*sigma*sigma);
    return 1 + 1e-8*( (a1+a2+a3)*a4*(a5/a6) - a7 );
}




/*
 * @internal
 * @brief computes Refractive index using lambda, pressure (P), temperature (T) 
 *        and partial pressure (Pw)
 * 
 * Formula used: Schubert, G. and Walterscheid, R. L. (2000). Earth. 
 * Allen’s astrophysical quantities.
 */
double RefractiveIndex_Schubert(double lambda, double p, double t, double p_w)
{
    
    if(lambda > 2.5 || lambda < 1.0 ) {
        printf("Warning: lambda should be within 1.0-2.5 \n");
        return -999;
    }
    
    if(p<0 || t < 0 || p_w < 0 ) {
        printf("Pressure, Temperature and Partial pressure needs to be positive \n");
        return -999;
    } 

    double p_s, Ts, A, B, C;   
    
    p_s = 1013.25; 
    Ts = 288.15; 
    A= 64.328 + 29498.1/(146.0 - 1/(lambda*lambda)) + 255.4/(41.0- 1/(lambda*lambda) );
    B= ((p*Ts)/(p_s*t))*1e-6;
    C= - 43.49*(1. - 7.956e-3/(lambda*lambda))*(p_w/p_s)*1e-6;

    return 1. + A*B  + C;
}

/*
 * @internal
 * @brief computes difference between true and apparent zenith distance 
 * using target refractive index (n)
 * and reference refractive index (nRef) and Zenith angle in degrees 
 * 
 * Formula used: Schubert, G. and Walterscheid, R. L. (2000). Earth. 
 * Allen’s astrophysical quantities.
 */
double ADR(double n, double nRef, double Zenith)
{
    
    if(n<0 || nRef < 0 ) {
        printf("refractive index needs to be positive \n");
        return -999;
    } 

    double n1, n2;
    n1= (n*n-1.0)/(2*n*n);
    n2= (nRef*nRef-1.0)/(2*nRef*nRef);
   
    return 206264.806247*(n2-n1)*tan(Zenith); /* 1rad= 206264.806247 as */
}

double partialPressure(double rh, double t)
{
    double ps;
    ps= -10474 + 116.43*t- 0.43284*t*t+ 0.00053840*t*t*t;
    return (rh/100.)*ps;
}


/*
 * @internal
 * @brief Top level ADR function which computes ADR induced shifts in arc sec
 * @param lambdaACQ Acquisition camera effective wavelegth (H-band)
 * @param lambdaBCI Beam combiner instrument effective wavelength (K-band)
 * @param P Pressure
 * @param T Temperature
 * @param rh humidity
 * @param Zenith Zenith angle in degrees
 * 
 * It returns ADR induced shifts in arc sec
 * 
 */
double ADRshift(double lambda, double lambdaRef, double P, double t, double rh, double Zenith)
{
    
    if(lambda > 2.5 || lambda < 1.0 ) {
        printf("Warning: lambda should be within 1.0-2.5 \n");
        return -999;
    }
    
    if(P<0 || t < 0 || rh < 0 || Zenith <0  || Zenith > 90) {
        printf("Pressure, Temperature (shall be in Kelvin), \
humidity percentage and Zenith angle (0-90 degrees) needs to be positive \n");
        return -999;
    } 

    double n, nRef;
   
    Zenith=Zenith*PI/180.0; /* convert defrees to radians */
    
    double Pw= partialPressure(rh,t);
    
    nRef=RefractiveIndex_Schubert(lambdaRef, P, t,  Pw);
    n= RefractiveIndex_Schubert(lambda, P, t,  Pw);    
    
    if(nRef==-999 || n ==-999 || nRef==n){
    printf("%s:%s:%d -----------------------------ADR errror \n",  
	   __FILE__, __FUNCTION__,  __LINE__);
    return 0;
    }
    return ADR(n, nRef, Zenith);
}

/*
 * @internal
 * @brief It evalutes effective wavelength from a 1d polynomial
 * @param polynomial coefficients of size 11
 * @param x Input value to estimate polynomial
 * 
 * return p[0] + x*( p[1] + x*( p[2] + x*( p[3] + x*( p[4] + x*( p[5] + x*( p[6] + x*( p[7] + x*( p[8] + x*( p[9] + x*p[10] ) ))))))));
 */
double polynomial_1d(double *p, double x)
{
    return p[0] + x*( p[1] + x*( p[2] + x*( p[3] + x*( p[4] + x*( p[5] + x*( p[6] + x*( p[7] + x*( p[8] + x*( p[9] + x*p[10] ) ))))))));   
}
/*
 * @internal
 * @brief It evalutes effective wavelength for H_K magnitude input using atan function
 * @param coefficients 
 * @param x Input H_K magnitude to estimate lambda
 * 
 * f(a, x)=a[0]+a[1]/PI*atan( a[2]*x-a[3] + PI/2 );
 * a=[1.37822,0.77782,0.730068,0.928234]
 */
double atanFit(double *a, double x){  
    return a[0]+a[1]/PI*atan( a[2]*x-a[3] + PI/2 );
}

/*
 * @internal
 * @brief It evalutes effective wavelength from a 1d polynomial
 * @param H_K star H_K magnitude
 * @param ACQPolyCoeff Acquisition camera H_K vs effective wavelength polynomial (H-band basically)
 * @param BCPolyCoeff  Beam combiner H_K vs effective wavelength polynomial (K-band basically)
 * @param EffectiveLambda (output) ACQ and BC effective wavelengths are returned 
 * return Error
 */
cpl_error_code EffectiveWavelength(double H_K, double *ACQPolyCoeff, 
				   double *BCPolyCoeff, double *EffectiveLambda)
{
    cpl_error_code error = CPL_ERROR_NONE;
    if( H_K >10 || H_K<-10) {
        printf("Expected H-K color is with -10 to 10, entered value is %f \n", H_K);
        return -1;
    }
    EffectiveLambda[0] = atanFit(ACQPolyCoeff, H_K); /* lambdaACQ  */
    EffectiveLambda[1] = atanFit(BCPolyCoeff, H_K); /* lambdaBC  */
    return error;
}
/*
 * @internal
 * @brief It evalutes effective wavelength from a 1d polynomial
 * @param H_K star H_K magnitude
 * @param ACQPolyCoeff Acquisition camera H_K vs effective wavelength polynomial (H-band basically)
 * @param BCPolyCoeff  Beam combiner H_K vs effective wavelength polynomial (K-band basically)
 * @param P Pressure
 * @param t Temperature
 * @param humidity humidity percentage
 * @param Zenith Zenith angle in degrees
 * @param ObjCorr (output) ADR biased error in arc sec
 * return Error
 */
cpl_error_code ADRCorrectionOffset(double H_K, double *ACQPolyCoeff, double *BCPolyCoeff, 
                                   double P, double t, double humidity, double Zenith, 
                                   double * EffectiveLambda, double *ObjCorr)
{
    cpl_error_code error = CPL_ERROR_NONE;
    double ObjCorrArcSec;
    if( H_K >10 || H_K<-10) {
        printf("ERROR:  Expected H-K color is with -10 to 10, entered value is %f \n", H_K);
        return -1;
    }
    error = EffectiveWavelength(H_K, ACQPolyCoeff, BCPolyCoeff, EffectiveLambda);
    if(error !=0) {
        printf("Error %s:%s:%d \n", __FILE__,__FUNCTION__, __LINE__);
        return error;
    }
    
    ObjCorrArcSec = ADRshift(EffectiveLambda[1], EffectiveLambda[0], P, t, humidity, Zenith);
    *ObjCorr=ObjCorrArcSec;

    return error;
}

/*
 * @internal
 * @brief It evalutes effective wavelength from a 1d polynomial
 * @param H_K star H_K magnitude
 * @param ACQPolyCoeff Acquisition camera H_K vs effective wavelength polynomial (H-band basically)
 * @param BCPolyCoeff  Beam combiner H_K vs effective wavelength polynomial (K-band basically)
 * @param P Pressure
 * @param t Temperature (degrees)
 * @param humidity humidity percentage
 * @param Zenith Zenith angle in degrees
 * @param PosAng Increasing Y-direction of detector array clockwise to the sky north vector
 * @param ParAng Angle from sky north vector counterclockwise to elevation up (Zenith) vector
 * @param XYoffset (output) Object XY offsets in pixels  
 * return Error
 */
cpl_error_code ADRCorrPixelsXY(double H_K, double *ACQPolyCoeff, double *BCPolyCoeff, 
			       double P, double t, 
				     double humidity, double Zenith, double*PosAng, 
			       double* ParAng, 
				     double FIPixelScale, double * EffectiveLambda, 
			       double *XYcorrArcSec, double *XYcorrPX)
{
    cpl_error_code error = CPL_ERROR_NONE;
    double ObjCorrArcSec;
    int i;
    double sumH=0;
    double sumK=0;
    if( H_K >10 || H_K<-10) {
        printf("ERROR:-> Expected H-K color is with -10 to 10, entered value is %f \n", H_K);
        return -1;
    }
    if ( ACQPolyCoeff == NULL || BCPolyCoeff == NULL ){
        printf("ERROR::-> Error in reading polynomial coefficients.\n");
        return -1;
    }
    
    if ( Zenith <0 || Zenith > 90 ){
        printf("ERROR::-> Error in reading Zenith angles, given value is  %f\n", Zenith);
    return -1;
    }

    if ( t <-20 || t > 40 ){
        printf("ERROR::-> Error in reading temperature, expected -20 to 40 degree, but given value is %f\n", t);
        return -1;
    }
    
    if ( P <600 || t > 900 ){
        printf("ERROR::-> Error in reading pressure, expected 600 to 900 millibars, but given value is %f\n", P);
    return -1;
    }

    for(i=0; i<11; i++){
        sumK +=BCPolyCoeff[i];
        sumH +=ACQPolyCoeff[i];
    }
    
    if(sumK == 0 || sumH == 0) {
     printf("ERROR::-> Error in reading ADR polynomials from database setup\n");
     return -1;
    }
    
    t=t+273.15;
    error = ADRCorrectionOffset(H_K, ACQPolyCoeff, BCPolyCoeff, P, t, 
				humidity, Zenith,  EffectiveLambda, &ObjCorrArcSec );
    if(error !=0) {
        printf("Error %s:%s:%d \n", __FILE__,__FUNCTION__, __LINE__);
        return error;
    }
    
    

/* FE changed sign of paralactic angle and subtract 90 degrees */
/* ErW: obviously Pos_ParAngle is not used anymore
    double Pos_ParAngle[4];
    for(i=0; i<4; i++){
      Pos_ParAngle[i] = (PosAng[i]+ParAng[i]-90)*PI/180.0;
    }
*/

    *XYcorrArcSec= ObjCorrArcSec;

    if (PRINTMSG) printf("ADR angle measured: %f \n", ObjCorrArcSec);

/* ErW: obviously ObjCorrPx is not used anymore
    double ObjCorrPx;
    if (FIPixelScale == 0) { 
        ObjCorrPx=0;
    } else {
        ObjCorrPx=ObjCorrArcSec/FIPixelScale;
    }
*/

/*
    for(i=0; i<4; i++){
    XYcorrPX [i]= ObjCorrPx*sin(Pos_ParAngle[i]);
    XYcorrPX [i+4]= ObjCorrPx*cos(Pos_ParAngle[i]);
    }
*/

 
    return error;
}


/* ----------------------------------------------------------------------*/
/**
 * @internal
 * @brief  Set the aberration sensor reference spot (5th spot) corrdinates
 * @param  FIcoordinates the spot window centre locations [x1,x2,x3,x4,y1,y2,y3,y4]
 * @see    setPTrefCoordinates, setABSrefCoordinates
 * 
 * Reads the coordinates from FIcoordinates and set to the SBA variable
 * 
 * @par Error Handling:
 */
void setFIrefCoordinates(int *FIcoordinates) {
	int i;
	for (i = 0; i < 8; i++) {
		SBA[i] = FIcoordinates[i];
	}
}


/* ----------------------------------------------------------------------*/
/**
 * @internal
 * @brief  Set the aberration sensor reference spot (5th spot corrdinates
 * @param  FIcoordinates the spot window centre locations [x1,x2,x3,x4,y1,y2,y3,y4]
 * @see    setPTrefCoordinates, setABSrefCoordinates
 * 
 * Reads the coordinates from FIcoordinates and set to the SBA variable
 * 
 * @par Error Handling:
 * CPL_ERROR_NULL_INPUT if(one of ) the input pointer (s) is NULL
 */

void setCoordinates(int* coordinates) {
    int i;
    for(i=0; i<8; i++) SBA[i]=coordinates[i];
    
    if(PRINTMSG){
        cpl_msg_info(cpl_func,"Set FI window Pos X (T1-T4):  %d %d %d %d\n", 
		     coordinates[0], coordinates[1], coordinates[2], coordinates[3]);	
        cpl_msg_info(cpl_func,"Set FI window Pos Y (T1-T4):  %d %d %d  %d  \n", 
		     coordinates[4], coordinates[5], coordinates[6], coordinates[7]);  
    }
}


/* ----------------------------------------------------------------------*/
/**
 * @internal
 * @brief  Sets the field imager field window size and the sigma value used to find the objects
 * @param  FIfieldWindowSizeAndSigma field imager window size (4" field around 220 pixels) and 
 *         sigma value used to find the objects cpl apertures
 * @see    setPTrefCoordinates, setABSrefCoordinates
 *  
 * @par Error Handling:
 */
void setFIfitWindowSizeAndSigma(double * FIfieldWindowSizeAndSigma) {
	int i;
	for (i = 0; i < 2; i++)
		FIWindowAndSigma[i] = FIfieldWindowSizeAndSigma[i];
	if(PRINTMSG) {cpl_msg_info(cpl_func,"FIWindowAndSigma %f %f \n", 
				   FIWindowAndSigma[0], FIWindowAndSigma[1]); }
}

/* ----------------------------------------------------------------------*/
/**
 * @internal
 * @brief  Sets the field imager relative position of star2 relative to star1
 * @param  tel, dx , dy
 * @see    setPTrefCoordinates, setABSrefCoordinates
 *  
 * @par Error Handling:
 */
void setFIStar2fromStar1(int tel, int dx, int dy) {

if ( (tel >=1) && (tel <=4) ) {
      FIStar2fromStar1[tel-1]=dx;
      FIStar2fromStar1[tel-1+4]=dy;
}

}


/*
  Once the bright object position is known roughly, this window size is used
  to fit Gaussian at the center of rough position
*/
void setFIwindowGauss(vltINT32 FIwindowGauss_in) {
    FIwindowGauss = FIwindowGauss_in;
    if(PRINTMSG) {cpl_msg_info(cpl_func,"FIwindowGauss = %d \n", FIwindowGauss);}
}
