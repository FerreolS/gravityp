/************************************************************************
 *   This file is part of the E.S.O. - VLTI GRAVITY project
 *
 *
 *   who       when        what
 *   --------  ----------  ----------------------------------------------
 *   aamorim  2016-08-16  Several Fixes in the PT counting of nspots...
 *   aamorim  2016-08-12  Fix problem if diffraction picture of two aligned spots  creates a connection between the two in the mask
 *   narsi    2016-05-07  Filtering pupil tracker spots with median filter
 *   narsi    2016-03-23  Add PT spots flux estimation proceders
 *   narsi    2016-02-20  modified to also work with only two subapertures illuminated
 *   narsi    2015-09-20  Splitted the gvacqProcessImageAC.c into gvacqPupilTracker.c
 *   narsi    2015-09-20  created
 */

/****************************************************************************
 *   NAME
 *   gvacqPupilTracker.c - Provides the lateral and longitudinal position of pupil,
 *         and measures reference spot position
 *
 *   SYNOPSIS
 *   #include <gvoProcessImageAC.h>
 *
 *   DESCRIPTION
 *   This function loads a acquistion camera fits file and determines the
 *   lateral and longitudinal  position of Unit Telescope pupil and write it to
 *   the database.
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
 **------------------------------------------------------------------------
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
 * Includes
 * ----------------------------------------------------------------------------*/

#include "gvoProcessImageAC.h" /* */
#include "cpl_apertures.h"


/*-----------------------------------------------------------------------------
 * Defines
 * ----------------------------------------------------------------------------*/
#define VLT2014
#define PI             3.14159265359
#define f_PT           14e-3    /* pupil tracker lenslet FL*/
#define subap_pupil    1.015e-3 /* STOP diameter of PT sub aperture */
#define f_lens         467e-3   /* folding optics lens FL */
#define Llambda        1.2e-6   /* laser diode wavelength */
#define Hlambda        1.65e-6  /* star beam wavelength */


/* These global variables are written by the database parameters */
static unsigned char   PTAlgo[64 * 4];
int           nTelAlgo;
static int    RefWindowPosPT[8];
static int    PTImageWindowSize          = 200; /* pixels, pupil tracker image window */
static int    PTwindowGauss              = 12;  /* pixels, Gauss fitting window size PT*/
static double PT_spots_scan_sigmaOrig[4] = { 10, 8, 6, 3.0 };
/* Positive, decreasing sigmas to apply for PT spots scan*/

static int    IMAGE_X                         = 2048; /* pixels */
static int    IMAGE_Y                         = 1536; /* pixels */
static int    RON                             = 12;
static double FIfwhm_minthreshold             = 1.99; /* pixels */
static int    PRINTMSG                        = 0;
static double TwoSubApertureLatPupilThreshold = 27;   /* px */


/*-----------------------------------------------------------------------------
 * Private function prototypes
 * ----------------------------------------------------------------------------*/
cpl_error_code gvacqPupilTrackerFCReferenceSpots(cpl_image *image,
                                                 double    *GuessCentre,                     /* Xc[0-3], Yc[4-7] */
                                                 double    *Fiber_coupler_reference_centers, /* Xcm[0-3], Ycm[4-7] */
                                                 double    *Fiber_coupler_reference_centers_error);

cpl_error_code gvacqPupilTrackerDetectedCentroidsOneTel(cpl_image *DetPointer,
                                                        int llx, int lly,
                                                        int urx, int ury,
                                                        double *Spot5RefTel,
                                                        double *PTspotsCen,
                                                        int *MissingSub_Apertures,
                                                        int *PredictedSpotsLocations);

cpl_error_code gvacqPupilTrackerComputeShifts(cpl_image *DetPointer,
                                              double *PTroughPositions,
                                              int llx, int lly,
                                              double *Spot5RefTel,
                                              int *MissingSubapertures,
                                              int *PredictedSpotsLocations,
                                              double *Barycentre,
                                              double *Barycentre_error,
                                              double *Pupil_positions,
                                              double *Pupil_positions_error,
                                              double *PTRefPosition_tmp,
                                              double *PTRefPositionError_tmp,
                                              double *Flux, int *Error_detected);

cpl_error_code gvacqPupilTrackerCentroids2BaryCenters(cpl_matrix *PTDetectedCenters,
                                                      cpl_matrix *PTDetectedCentersError,
                                                      int llx, int lly,
                                                      double *PTBarycenters,
                                                      double *PTBarycentersError);


cpl_error_code gvacqPupilTrackerBaryCenters2PupilPositions(double *Fiber_coupler_reference_centers,
                                                           double *PTBarycenters,
                                                           double *PTBarycentersError,
                                                           int    *MissingSubapertures,
                                                           double *PupilPosition,
                                                           double *PupilPositionError);
double gvacqPupilTrackerLongPostionPixels2Meters(double PositionPixels);
cpl_error_code gvacqPupilTrackerSpotsDetectionTwoSubs(cpl_image *Detector,
                                                      cpl_apertures *Spots,
                                                      int llx, int lly,
                                                      double *OneTelSpot5,
                                                      double *SpotsData,
                                                      int *MissingTelescopes,
                                                      int *PredictedSpotsLocations);


cpl_error_code gvacqPupilTrackerSpotsDetectionThreeORFourSubs(cpl_image *DetPointer,
                                                              cpl_apertures *PTspots,
                                                              int llx, int lly, int urx, int ury,
                                                              double *OneTelSpot5,
                                                              cpl_vector *sig,
                                                              double *PTspotsCen,
                                                              int *MissingSub_Apertures,
                                                              int *PredictedSpotsLocations);

cpl_error_code gvacqPupilShiftsMeasureFromCentroids(double *, cpl_matrix *,
                                                    cpl_matrix *,
                                                    double *, double *);


cpl_error_code cpl_apertures2XYdoubleArray(cpl_apertures *apert, int nx, double *DataDouble);

cpl_error_code Predict4thSpotWithThree(cpl_apertures *ptspots1, double *Data4Spots);

cpl_error_code Check_objects_found_PT(cpl_image *FI_image,
                                      double    fwhm_threshold_min,
                                      double    fwhm_threshold_max);

cpl_error_code PTSpotsFlux(cpl_image *fimage,
                           int       *imistart2,
                           int       *imjstart2,
                           int       nspots,
                           double    *Flux);

double minApertureDistance(cpl_apertures *PTspots);

int spotsInApertureMultiplicity(cpl_apertures *PTspots);

unsigned char *getrefPTAlgo() {
    return(PTAlgo);
}


/*
 *******************************************************************************
 *                                                                             *
 *      Pupil Tracker                                                          *
 *                                                                             *
 *******************************************************************************
 */
/* ----------------------------------------------------------------------*/

/**
 * @internal
 * @brief  Provides the lateral and longitudinal position of pupil,
 *         and measures reference spot position
 * @param  DetPointer   pointer to a  detector image
 * @param  Spot5RefData Fibe coupler laser reference positions
 * @param  FIfwhm_minthreshold fwhm thershold used to scan if there is spots or only hot pixels?
 * @param  PupilPosition (output) pupil positions for each telescope [4Tel.*3(X, Y, Z)]
 * @param  PupilPositionError (output) pupil position errors for each telescope [4Tel.*3(X, Y, Z)]
 * @param  barycentre (output) avg position of 4 spots (4 lasers imaged on lenslet) [4Tel.*4Spots*2(x,y)]
 * @param  barycentre (output) avg position of 4 spots error (4 lasers imaged on lenslet) [4Tel.*4Spots*2(x,y)]
 * @param  PTRefPosition (output) pupil tracker spot positions [4Tel.*4Spots*16(x,y)]
 * @param  PTRefPositionError (output) pupil tracker spot positions error [4Tel.*16Spots*2(x,y)]
 * @return error
 * @see    gvacqPupilTrackerCalib, gvacqPupilShiftsMeasureFromCentroids,  gvacqPupilTrackerFCReferenceSpots
 *
 * This function loads a acquistion camera fits file and determines the
 * lateral and longitudinal  position of Unit Telescope pupil and write it to
 * the database.
 *
 * FLOW:
 * 1. check inputs
 * 2. Fit the spots using Gaussian function. Calculate the spot shifts by comparing
 *    reference spot positions from func gvacqPupilTrackerCalib.
 * 3. Measure the slopes from the spot shifts and then fit the
 *    Zernike coefficients and then pupil shifts function func gvacqPupilShiftsMeasureFromCentroids.
 *
 * @par Description:
 * - The function clipm_centroiding_gauss() is invoked @a N times, for @a N
 *   different square windows of @a image. The windows are defined by
 *   their centers, given by @a locations, which must be of the form:
 *
 * Functions pre-required:
 * func setPTImageWindowSize sets pupil tracker window size (~ 220 pixels)
 * FUNC setRefWindowPosPT sets pupil tracker window centers for 4 telescopes
 * func setPTwindowGauss sets pupil tracker spots finding Gauss window (~ 8 pixels)
 * func setPTrefCoordinates sets pupil tracker reference positions(4 centers * 2(x,y) * 4telescopes)
 *
 *
 *  Here I am considering the spots in clockwise direction as mentioned
 * like  below
 * 2x2  lenslet sub-aperture  positions
 *
 * _______________________________
 *
 *
 *   b. (xb,yb)         c. (xc,yc)
 *
 *
 *
 *   a. (xa,yb)         d. (xd,yd)
 *
 *   _______________________________
 *
 * For the description of the above result  refer
 * Salas-Peimbert et al. 2005.
 *
 *
 * The lateral pupil positions [lx, ly] can be  measured as below
 * lx=[xa+xb+xc+xd]/4 in pixels
 * lx= lx*59.85 meters  at UT 8m beam
 *
 * ly=[ya+yb+yc+yd]/4 in pixels
 * ly= ly*59.85 meters  at UT 8m beam
 *
 * Where xa, xb, xc, xd are the spot shifts (distorted spot position - fiber coupler reference position)
 * caused by the aberrations.
 *
 * aberration defocus coefficient Ad = -[ (xa-xd) - (xc-xb) + (ya+yd) -(yb+yc) ]/8/f_PT radians
 * Where d is the diamter of lenslet
 *
 * @par Error Handling:
 * CPL_ERROR_NULL_INPUT if(one of ) the input pointer (s) is NULL
 * CPL_ERROR_NULL
 * return -1, NULL input image
 * return -2, pupil guiding lasers are not turned ON (no data in the image)
 * return -3, pupil tracker data consists huge background
 * All the above cases PupilPosition return as -100
 */
cpl_error_code gvacqPupilTracker(cpl_image *DetPointer,
                                 double    *Spot5RefData,
                                 double    FIfwhm_minthreshold,
                                 double    PTFitStateErrorThreshold,
                                 double    *PupilPosition,
                                 double    *PupilPositionError,
                                 double    *barycentre,
                                 double    *barycentre_error,
                                 double    *PTRefPosition,
                                 double    *PTRefPositionError,
                                 double    *SpotsFlux,
                                 int       *PTErrorStatus) {
    /* test entries */
    int i, j;

    for (i = 0; i < 32 * 4; i++) {
        PTRefPosition[i] = 0.0;
    }
    for (i = 0; i < 12; i++) {
        PupilPosition[i] = 0.0;
    }
    for (i = 0; i < 4 * 2 * 4; i++) {
        barycentre[i] = 0;
    }
    int win = 100;
    for (j = 0; j < 4; j++) {
        for (i = 0; i < 64; i++) {
            PTAlgo[4 * i + j] = 0;
        }
    }

    if ((DetPointer == NULL) ||
        (PupilPosition == NULL) || (PupilPositionError == NULL)) {
        printf("ERROR::The input parameters are NULL. Occured at %s:%s:%d\n",
               __FILE__, __FUNCTION__, __LINE__);
        for (i = 0; i < 4; i++) {
            PTErrorStatus[i] = -1;
        }
        return(-1);
    }

    if ((cpl_image_get_size_x(DetPointer) != IMAGE_X) ||
        (cpl_image_get_size_y(DetPointer) != IMAGE_Y)) {
        printf(":: Return: The input image dimensions are %dx%d. Expected image is 2048x1536 \n",
               (int)cpl_image_get_size_x(DetPointer),
               (int)cpl_image_get_size_y(DetPointer));

        for (i = 0; i < 4; i++) {
            PTErrorStatus[i] = -1;
        }
        return(-1);
    }

    if ((RefWindowPosPT [0] < 295 - win) || (RefWindowPosPT[0] > 295 + win) ||
        (RefWindowPosPT[1] < 770 - win) || (RefWindowPosPT[1] > 770 + win) ||
        (RefWindowPosPT[2] < 1230 - win) || (RefWindowPosPT[2] > 1230 + win) ||
        (RefWindowPosPT[3] < 1700 - win) || (RefWindowPosPT[3] > 1700 + win) ||
        (RefWindowPosPT[4] < 1380 - win) || (RefWindowPosPT[4] > 1380 + win) ||
        (RefWindowPosPT[5] < 1380 - win) || (RefWindowPosPT[5] > 1380 + win) ||
        (RefWindowPosPT[6] < 1380 - win) || (RefWindowPosPT[6] > 1380 + win) ||
        (RefWindowPosPT[7] < 1380 - win) || (RefWindowPosPT[7] > 1380 + win)) {
        for (i = 0; i < 8; i++) {
            printf("%d ", RefWindowPosPT [i]);
        }
        printf("ERROR: The pupil tracker window positions read from \
database are wrong. Please check it!!! \n");
        for (i = 0; i < 4; i++) {
            PTErrorStatus[i] = -1;
        }
        return(-1);
    }


    if ((PTImageWindowSize <= 110) || (PTImageWindowSize >= 350)) {
        printf("ERROR: PT window size expected 110-300 but given %d. \
Reading from database is wrong. Please check it!!! \n", PTImageWindowSize);
        for (i = 0; i < 4; i++) {
            PTErrorStatus[i] = -1;
        }
        return(-1);
    }

    if ((PTwindowGauss <= 3) || (PTwindowGauss >= 17)) {
        printf("ERROR: PT Gaussian window size expected 4-16 but given %d. \
Reading from database is wrong. Please check it!!! \n", PTwindowGauss);
        for (i = 0; i < 4; i++) {
            PTErrorStatus[i] = -1;
        }
        return(-1);
    }

    for (i = 0; i < 4; i++) {
        if ((PT_spots_scan_sigmaOrig[0] > 15) || (PT_spots_scan_sigmaOrig[0] < 1)) {
            printf("ERROR: PT spots scan sigma expected 2-15 but given %f. \
Reading from database is wrong. Please check it!!! \n", PT_spots_scan_sigmaOrig[0]);
            for (i = 0; i < 4; i++) {
                PTErrorStatus[i] = -1;
            }
            return(-1);
        }
    }
    if (PTwindowGauss >= 100) {
        printf("ERROR: PT Gaussian window size expected 0-100 but given %d. \
Reading from database is wrong. Please check it!!! \n", PTwindowGauss);
        for (i = 0; i < 4; i++) {
            PTErrorStatus[i] = -1;
        }
        return(-1);
    }

    if ((FIfwhm_minthreshold < 1.0) || (FIfwhm_minthreshold > 2.0)) {
        printf("ERROR: FIfwhm_minthreshold size expected 1-2 but given %f. \
Reading from database is wrong. Please check it!!! \n", FIfwhm_minthreshold);
        for (i = 0; i < 4; i++) {
            PTErrorStatus[i] = -1;
        }
        return(-1);
    }


    /* initialising variables */
    cpl_error_code error  = CPL_ERROR_NONE;
    int            Nspots = 16; /* # of spots avialble in PT mode */
    double         Pupil_positions[3];
    double         Pupil_positions_error[3];
    double         barycentre_temp[8];
    double         barycentre_error_temp[8];
    /* X-spots[0-15], Y-spots[16-32] for Tel.1 */
    cpl_matrix *UT1locations = cpl_matrix_new(2, Nspots);
    double     PTroughPositions[32];
    double     PTguessPositions[32];
    double     Spot5RefTel[8];
    int        Error_detected[3];
    double     PTRefPosition_tmp[32];
    double     PTRefPositionError_tmp[32];
    int        MissingSubapertures[4] = { 0, 0, 0, 0 };
    int        PredictedSpotsLocations[32];
    int        llx, urx, lly, ury;
    int        MissingSubs;
    double     FluxONorOFF[4];
    double     SpotsFluxPerTel[16];

    /* measured flux of pupil window; used for blink mode   */
    error = PupilWindowFlux(DetPointer,
                            RefWindowPosPT,
                            PTImageWindowSize,
                            FluxONorOFF);

    if (PRINTMSG) printf("Flux of PT window (subtracted): %f %f %f %f \n",
                         FluxONorOFF[0], FluxONorOFF[1],
                         FluxONorOFF[2], FluxONorOFF[3]);

    if (FluxONorOFF[1] < 0.0) {
        /*error = cpl_image_multiply_scalar(DetPointer,-1.0);*/
    }

    /* Estimate pupil positions for each telescope in a loop  */
    for (i = 0; i < 4; i++) { /* Do it for all four telescopes  */
        llx = RefWindowPosPT[i] - PTImageWindowSize / 2;
        lly = RefWindowPosPT[i + 4] - PTImageWindowSize / 2;
        urx = RefWindowPosPT[i] + PTImageWindowSize / 2;
        ury = RefWindowPosPT[i + 4] + PTImageWindowSize / 2;

        /*PTErrorStatus[i]=error;*/
        MissingSubs = 0;

        for (j = 0; j < 32; j++) {
            PredictedSpotsLocations[j] = 0.0;
        }
        for (j = 0; j < 4; j++) {
            Spot5RefTel[j]     = Spot5RefData[j + 4 * i];      /* x - centers  */
            Spot5RefTel[j + 4] = Spot5RefData[16 + j + 4 * i]; /* y - centers  */
        }

        cpl_image *DetPointer_dup = cpl_image_duplicate(DetPointer);

        /*
         * Estimate rough positions of pupil tracker spots using
         * cpl_apertures detection using sigma as threshold (PTroughPositions)
         */
        nTelAlgo = i;
        error    = gvacqPupilTrackerDetectedCentroidsOneTel(DetPointer_dup,
                                                            llx, /* llx */
                                                            lly, /* lly */
                                                            urx, /* urx */
                                                            ury, /* ury */
                                                            Spot5RefTel,
                                                            PTroughPositions,
                                                            MissingSubapertures,
                                                            PredictedSpotsLocations);
        for (j = 0; j < 4; j++) {
            if (MissingSubapertures[j] != 0) MissingSubs++;
        }
        cpl_image_delete(DetPointer_dup);

        PTErrorStatus[i] = error;


        if ((error == 0) && (MissingSubs <= 2)) {
            /* Extracted PT window  */
            cpl_image *PTs_im = cpl_image_extract(DetPointer,
                                                  llx, /* llx */
                                                  lly, /* lly */
                                                  urx, /* urx */
                                                  ury /* ury */);

            /* Saturated pixels rejection before Gauss fit in the next function*/
            cpl_mask *PTMask = cpl_mask_threshold_image_create(PTs_im, -0.9999999, 0.000001);

            cpl_mask *DetMask = cpl_mask_new(IMAGE_X, IMAGE_Y);

            error = cpl_mask_copy(DetMask, PTMask,
                                  RefWindowPosPT[i] - PTImageWindowSize / 2,
                                  RefWindowPosPT[i + 4] - PTImageWindowSize / 2);

            error = cpl_image_reject_from_mask(DetPointer, DetMask);

            /* End of saturated pixels rejection */
            for (j = 0; j < 16; j++) {
                PTguessPositions[j]      = PTroughPositions[j] - llx + 1;
                PTguessPositions[j + 16] = PTroughPositions[j + 16] - lly + 1;
            }


            /*
             * 1. Do Gauss fit on rough estimated positions
             * 2. Compared estimated positions with reference positions to estimate
             *   pupil positions
             */
            error = gvacqPupilTrackerComputeShifts(PTs_im, PTguessPositions,
                                                   llx, lly, Spot5RefTel,
                                                   MissingSubapertures,
                                                   PredictedSpotsLocations,
                                                   barycentre_temp,
                                                   barycentre_error_temp,
                                                   Pupil_positions,
                                                   Pupil_positions_error,
                                                   PTRefPosition_tmp,
                                                   PTRefPositionError_tmp,
                                                   SpotsFluxPerTel,
                                                   Error_detected);



            cpl_image_delete(PTs_im);
            cpl_mask_delete(PTMask);
            cpl_mask_delete(DetMask);

            for (j = 0; j < 3; j++) {
                if (Error_detected[j] != 0) PTErrorStatus[i] = Error_detected[j];
            }
            for (j = 0; j < 3; j++) {
                PupilPosition[j + i * 3]      = Pupil_positions[j];
                PupilPositionError[j + i * 3] = Pupil_positions_error[j];
            } /* j*/

            for (j = 0; j < 16; j++) {
                SpotsFlux[j + i * 16] = SpotsFluxPerTel[j];
            }

            for (j = 0; j < 4; j++) {
                barycentre[j + i * 4]      = barycentre_temp[j];
                barycentre[j + i * 4 + 16] = barycentre_temp[j + 4];

                barycentre_error[j + i * 4]      = barycentre_error_temp[j];
                barycentre_error[j + i * 4 + 16] = barycentre_error_temp[j + 4];
            } /* j*/

            for (j = 0; j < 16; j++) {
                PTRefPosition[j + i * 32]      = PTRefPosition_tmp[j];
                PTRefPosition[j + 16 + i * 32] = PTRefPosition_tmp[j + 16];

                PTRefPositionError[j + i * 32]      = PTRefPositionError_tmp[j];
                PTRefPositionError[j + 16 + i * 32] = PTRefPositionError_tmp[j];
            } /* j*/
        } else {
            PTErrorStatus[i] = -3;

            for (j = 0; j < 3; j++) {
                PupilPosition[j + i * 3]      = -999;
                PupilPositionError[j + i * 3] = -999;
            }
        } /* End of if error != -1 || error!= -2 || error !=-3; end of error handling */
    }     /* END (i=0; i<4; i++) */

    /* free up the  memory */
    cpl_matrix_delete(UT1locations);
    return(error);
} /* End of gvacqPupilTracker */


/**
 * @brief  Computes pupil positions/displacements
 * @param  PTimage pupil tracker window of a given telescope
 * @param  PTroughPositions pupil spots rough positions to apply Gaussian fits
 * @param  llx window low left position x
 * @param  lly window low left position y
 * @param  Spot5RefTel pupil tracker reference positions [4 centers x 2(x,y)]
 * @param  (output) Pupil_positions computed pupil positions
 * @param  (output) Pupil_positions_error computed pupil positions error
 * @param  (output) PTRefPosition_tmp pupil tracker spots positions [16 spots x 2 (x,y)]
 * @param  (output) PTRefPositionError_tmp pupil tracker spots positions error  [16 spots x 2 (x,y)]
 * @param  (output) Error_detected error detected during computations
 * @return   error
 */
cpl_error_code gvacqPupilTrackerComputeShifts(cpl_image *PTimage,
                                              double *PTroughPositions,
                                              int llx, int lly,
                                              double *Spot5RefTel,
                                              int *MissingSubapertures,
                                              int *PredictedSpotsLocations,
                                              double *Barycentre,
                                              double *Barycentre_error,
                                              double *Pupil_positions,
                                              double *Pupil_positions_error,
                                              double *PTRefPosition_tmp,
                                              double *PTRefPositionError_tmp,
                                              double *Flux, int *Error_detected) {
    cpl_error_code error  = CPL_ERROR_NONE;
    int            Nspots = 16;
    int            i, j;

    int        robustness          = 7;
    cpl_matrix *UT1locations       = cpl_matrix_new(2, Nspots);
    cpl_matrix *xy_centre          = cpl_matrix_new(2, Nspots);
    cpl_matrix *xy_centre_err      = cpl_matrix_new(2, Nspots);
    cpl_matrix *xy_sigma           = cpl_matrix_new(2, Nspots);
    cpl_matrix *xy_sigma_err       = cpl_matrix_new(2, Nspots);
    cpl_matrix *xy_fwhm            = cpl_matrix_new(2, Nspots);
    cpl_matrix *xy_fwhm_err        = cpl_matrix_new(2, Nspots);
    cpl_matrix *centre_intensities = cpl_matrix_new(1, Nspots);
    cpl_array  *all_error_codes    = cpl_array_new(Nspots, CPL_TYPE_INT);
    double     FCspotsData[8];
    int        imistart2[Nspots];
    int        imjstart2[Nspots];

    /* Copying positions of spots (x,y) to 2x16 matrix */
    for (i = 0; i < Nspots; i++) {
        cpl_matrix_set(UT1locations, 0, i, PTroughPositions[i]);          /*X-spots  Tel.1 */
        cpl_matrix_set(UT1locations, 1, i, PTroughPositions[i + Nspots]); /* Y-spots  Tel. 2*/
        imistart2[i] = (int)PTroughPositions[i];
        imjstart2[i] = (int)PTroughPositions[i + Nspots];
    }


    error = PTSpotsFlux(PTimage, imistart2, imjstart2, Nspots, Flux);

    /*-------------- STEP 2:: On the rough centers fit 2d Gaussian-----------*/

    error = clipm_centroiding_multi_gauss(PTimage,
                                          UT1locations,
                                          PTwindowGauss,
                                          &xy_centre,
                                          &xy_centre_err,
                                          &xy_sigma,
                                          &xy_sigma_err,
                                          &xy_fwhm,
                                          &xy_fwhm_err,
                                          &centre_intensities,
                                          &all_error_codes,
                                          robustness);
    Flux = cpl_matrix_get_data(centre_intensities);
    Error_detected[0] = error;

    /*When a spot is not there, Gauss fit give centroid error.
     * So reading directly the previoulsy predicted spot */

    for (j = 0; j < Nspots; j++) {
        if ((((2 > cpl_matrix_get(xy_fwhm, 0, j)) || (cpl_matrix_get(xy_fwhm, 0, j) > 8)) &&
             ((2 > cpl_matrix_get(xy_fwhm, 1, j)) || (cpl_matrix_get(xy_fwhm, 1, j) > 8))) &&
            (PredictedSpotsLocations[j] == 0)) {
            if (j < 4) MissingSubapertures[0] = 1;
            if ((j >= 4) && (j < 8)) MissingSubapertures[1] = 1;
            if ((j >= 8) && (j < 12)) MissingSubapertures[2] = 1;
            if ((j >= 12) && (j < 16)) MissingSubapertures[3] = 1;
        }

        if (PredictedSpotsLocations[j] == 1) {
            error = cpl_matrix_set(xy_centre, 0, j, PTroughPositions[j]);
            error = cpl_matrix_set(xy_centre, 1, j, PTroughPositions[j + Nspots]);
        }
    }

    if (error != CPL_ERROR_NONE)
        cpl_msg_warning(cpl_func, "Gauss fit fails at %s:%s:%d\n", __FILE__,
                        __FUNCTION__, __LINE__);

    /*--------------- STEP 3: Read the fiber coupler reference spots from instrument database-----*/
    for (j = 0; j < 4; j++) {
        FCspotsData[j]     = Spot5RefTel[j];     /* x - centers  */
        FCspotsData[j + 4] = Spot5RefTel[j + 4]; /* y - centers  */
    } /* j */


    /* ---------STEP 4: Measure the 4 barycenetrs of the pupil tracker spots ------*/
    error = gvacqPupilTrackerCentroids2BaryCenters(xy_centre,
                                                   xy_centre_err,
                                                   llx, lly,
                                                   Barycentre,
                                                   Barycentre_error);

    Error_detected[1] = error;
    /*---STEP 5: Using the fiber coupler refs and barycenters measure the pupil positions --------*/
    error = gvacqPupilTrackerBaryCenters2PupilPositions(FCspotsData,
                                                        Barycentre,
                                                        Barycentre_error,
                                                        MissingSubapertures,
                                                        Pupil_positions,
                                                        Pupil_positions_error);
    Error_detected[2] = error;
    /* Fill the arrays for database writing */
    for (j = 0; j < 16; j++) {
        PTRefPosition_tmp[j]      = cpl_matrix_get(xy_centre, 0, j) + llx - 1;
        PTRefPosition_tmp[j + 16] = cpl_matrix_get(xy_centre, 1, j) + lly - 1;

        PTRefPositionError_tmp[j]      = cpl_matrix_get(xy_centre_err, 0, j);
        PTRefPositionError_tmp[j + 16] = cpl_matrix_get(xy_centre_err, 1, j);
    } /* j*/


    /* free up the  memory */
    cpl_matrix_delete(UT1locations);
    cpl_matrix_delete(xy_centre);
    cpl_matrix_delete(xy_centre_err);
    cpl_matrix_delete(xy_sigma);
    cpl_matrix_delete(xy_sigma_err);
    cpl_matrix_delete(xy_fwhm);
    cpl_matrix_delete(xy_fwhm_err);
    cpl_matrix_delete(centre_intensities);
    cpl_array_delete(all_error_codes);

    return(error);
}


/*
 * This function is called when two lenslet sub-apertures found
 * in the pupil tracker image.
 */
cpl_error_code gvacqPupilTrackerSpotsDetectionTwoSubs(cpl_image *DetPointer,
                                                      cpl_apertures *PTspots,
                                                      int llx, int lly,
                                                      double *OneTelSpot5,
                                                      double *PTspotsCen,
                                                      int *MissingSub_Apertures,
                                                      int *PredictedSpots) {
    cpl_error_code error = CPL_ERROR_NONE;
    int            i, j, Nspots;
    double         psig[4];
    double         px, py;
    int            Left         = 0, Right = 0, Bottom = 0, Top = 0;
    double         Avg_2sub_x   = 0.0;
    double         Avg_2sub_y   = 0.0;
    double         Avg_sub_x[2] = { 0.0, 0.0 };
    double         Avg_sub_y[2] = { 0.0, 0.0 };
/*    int extract1=20; modify to 25 to cope with UT rotating in all directions */
    int extract1 = 25;
    int extract2 = 50;
    int N_a      = 0, N_b = 0, N_c = 0, N_d = 0;

#ifdef VLT2014
    cpl_size pos;
#else
    int pos;
#endif

    Nspots = cpl_apertures_get_size(PTspots);

    for (i = 0; i <= 3; i++) {
        psig[i] = PT_spots_scan_sigmaOrig[i] / 3.0;
    }
    cpl_vector *sig = cpl_vector_wrap(4, (double *)psig);

    for (j = 0; j < Nspots; j++) {
        px = (double)cpl_apertures_get_centroid_x(PTspots, j + 1);
        py = (double)cpl_apertures_get_centroid_y(PTspots, j + 1);

        Avg_2sub_x += px;
        Avg_2sub_y += py;
    }

    /* Find the ceneters of pupil window using avilable spots */
    Avg_2sub_x = Avg_2sub_x / (double)Nspots;
    Avg_2sub_y = Avg_2sub_y / (double)Nspots;

    Avg_2sub_x = Avg_2sub_x + llx - 1;
    Avg_2sub_y = Avg_2sub_y + lly - 1;

    /* storing the division of the different subapertures for debug */
    PTAlgo[30*4+nTelAlgo]=Avg_2sub_x-(RefWindowPosPT[nTelAlgo]-PTImageWindowSize/2);
    PTAlgo[31*4+nTelAlgo]=Avg_2sub_y-(RefWindowPosPT[nTelAlgo+4]-PTImageWindowSize/2);

    /*
     * Devide the pupil tracker window
     *
     * --------
     * . * *  .
     * .      . -->  ptspots8_b
     * . * *  .
     * -------
     * . * *  .
     * .      . -->  ptspots8_a
     * . * *  .
     * -------
     *
     *
     *          ------------------
     *           . * *  .   * *  .
     * ptspots8_c -->  .      .        .    <--- ptspots8_d
     *           . * *  .   * *  .
     *          ------------------
     */

    cpl_apertures *ptspots8_a = cpl_apertures_extract_window(DetPointer,
                                                             sig,
                                                             Avg_2sub_x - extract1,
                                                             Avg_2sub_y - extract2,
                                                             Avg_2sub_x + extract1,
                                                             Avg_2sub_y,
                                                             &pos);

    cpl_apertures *ptspots8_b = cpl_apertures_extract_window(DetPointer,
                                                             sig,
                                                             Avg_2sub_x - extract1,
                                                             Avg_2sub_y,
                                                             Avg_2sub_x + extract1,
                                                             Avg_2sub_y + extract2,
                                                             &pos);

    cpl_apertures *ptspots8_c = cpl_apertures_extract_window(DetPointer,
                                                             sig,
                                                             Avg_2sub_x - extract2,
                                                             Avg_2sub_y - extract1,
                                                             Avg_2sub_x,
                                                             Avg_2sub_y + extract1,
                                                             &pos);

    cpl_apertures *ptspots8_d = cpl_apertures_extract_window(DetPointer,
                                                             sig,
                                                             Avg_2sub_x,
                                                             Avg_2sub_y - extract1,
                                                             Avg_2sub_x + extract2,
                                                             Avg_2sub_y + extract1,
                                                             &pos);

    if (((ptspots8_a == NULL) && (ptspots8_b == NULL)) ||
        ((ptspots8_c == NULL) && (ptspots8_d == NULL))) {
        if (ptspots8_a != NULL) cpl_apertures_delete(ptspots8_a);
        if (ptspots8_b != NULL) cpl_apertures_delete(ptspots8_b);
        if (ptspots8_c != NULL) cpl_apertures_delete(ptspots8_c);
        if (ptspots8_d != NULL) cpl_apertures_delete(ptspots8_d);
        cpl_vector_unwrap(sig);
        return(-3);
    } else {
        if (ptspots8_a != NULL) N_a = (int)cpl_apertures_get_size(ptspots8_a);
        if (ptspots8_b != NULL) N_b = (int)cpl_apertures_get_size(ptspots8_b);
        if (ptspots8_c != NULL) N_c = (int)cpl_apertures_get_size(ptspots8_c);
        if (ptspots8_d != NULL) N_d = (int)cpl_apertures_get_size(ptspots8_d);
    }

    if (((N_a == -1) && (N_b == -1)) || ((N_c == -1) && (N_d == -1))) {
        if (ptspots8_a != NULL) cpl_apertures_delete(ptspots8_a);
        if (ptspots8_b != NULL) cpl_apertures_delete(ptspots8_b);
        if (ptspots8_c != NULL) cpl_apertures_delete(ptspots8_c);
        if (ptspots8_d != NULL) cpl_apertures_delete(ptspots8_d);
        cpl_vector_unwrap(sig);
        return(-3);
    }

    int    Pa         = 0, Pb = 0;
    int    Vertical   = 0;
    int    Horizental = 0;
    double Sub_Ap_spots1[8];
    double Sub_Ap_spots2[8];
    int    Error_a = 0, Error_b = 0, Error_c = 0, Error_d = 0;

    PTAlgo[4 * 17 + nTelAlgo] = N_a;
    PTAlgo[4 * 18 + nTelAlgo] = N_b;
    PTAlgo[4 * 19 + nTelAlgo] = N_c;
    PTAlgo[4 * 20 + nTelAlgo] = N_d;

    /* Checking for vertical alignment */
    if ((N_a + N_b) > (N_c + N_d)) {
        Vertical = 1;

        /*Read positions from cpl_apertures */
        if (N_a >= 4)
            error = cpl_apertures2XYdoubleArray(ptspots8_a, 4, Sub_Ap_spots1);
        else if (N_a == 3) {
            error = Predict4thSpotWithThree(ptspots8_a, Sub_Ap_spots1);
            Pa    = 1;
        } else Error_a = -1;


        if (N_b >= 4)
            error = cpl_apertures2XYdoubleArray(ptspots8_b, 4, Sub_Ap_spots2);
        else if (N_b == 3) {
            error = Predict4thSpotWithThree(ptspots8_b, Sub_Ap_spots2);
            Pb    = 1;
        } else Error_b = -1;

        /* Checking for horizental alignment */
    } else if ((N_a + N_b) < (N_c + N_d))  {
        Horizental = 1;

        if (N_c >= 4)
            /*Read positions from cpl_apertures */
            error = cpl_apertures2XYdoubleArray(ptspots8_c, 4, Sub_Ap_spots1);
        else if (N_c == 3) {
            error = Predict4thSpotWithThree(ptspots8_c, Sub_Ap_spots1);
            Pa    = 1;
        } else Error_c = -1;

        if (N_d >= 4)
            error = cpl_apertures2XYdoubleArray(ptspots8_d, 4, Sub_Ap_spots2);
        else if (N_d == 3) {
            error = Predict4thSpotWithThree(ptspots8_d, Sub_Ap_spots2);
            Pb    = 1;
        } else Error_d = -1;
    } else {
        Error_a = -1;
        Error_b = -1;
        Error_c = -1;
        Error_d = -1;
    }

    if ((Error_a == -1) || (Error_b == -1) || (Error_c == -1) || (Error_d == -1)) {
        printf("Error: Sub-apertures spots not detected properly(N_a, N_b, N_c, N_d)=\
(%d,%d,%d,%d); %s:%d\n", N_a, N_b, N_c, N_d, __FUNCTION__, __LINE__);
        if (ptspots8_a != NULL) cpl_apertures_delete(ptspots8_a);
        if (ptspots8_b != NULL) cpl_apertures_delete(ptspots8_b);
        if (ptspots8_c != NULL) cpl_apertures_delete(ptspots8_c);
        if (ptspots8_d != NULL) cpl_apertures_delete(ptspots8_d);
        cpl_vector_unwrap(sig);
        return(-3);
    }

    for (j = 0; j < 4; j++) {
        Avg_sub_x[0] += Sub_Ap_spots1[j];
        Avg_sub_y[0] += Sub_Ap_spots1[j + 4];
    }

    for (j = 0; j < 4; j++) {
        Avg_sub_x[1] += Sub_Ap_spots2[j];
        Avg_sub_y[1] += Sub_Ap_spots2[j + 4];
    }

    Avg_sub_x[0] = Avg_sub_x[0] / 4.0;
    Avg_sub_y[0] = Avg_sub_y[0] / 4.0;
    Avg_sub_x[1] = Avg_sub_x[1] / 4.0;
    Avg_sub_y[1] = Avg_sub_y[1] / 4.0; /*check if two spots are less than 3 */


    /*
     * Now check for whether the vertical alignment is avilable right side or left side ?
     * Now check for whether the horizental alignment is avilable bottom or top ?
     */
    if (Avg_sub_x[0] - OneTelSpot5[0] < TwoSubApertureLatPupilThreshold) Left = 1;
    if (Avg_sub_x[0] - OneTelSpot5[3] > -TwoSubApertureLatPupilThreshold) Right = 1;
    if (Avg_sub_y[0] - OneTelSpot5[4] < TwoSubApertureLatPupilThreshold) Bottom = 1;
    if (Avg_sub_y[0] - OneTelSpot5[5] > -TwoSubApertureLatPupilThreshold) Top = 1;


    /* If Vertical alignment  */
    if ((Vertical == 1) && (Horizental == 0)) {
/*
 * Now check for whether the vertical alignment is avilable right side or left side
 */

        if (Left) {
            if (PRINTMSG) printf("Executing-->Pattern: Vertical and Left \n");
            MissingSub_Apertures[2] = 3;
            MissingSub_Apertures[3] = 4;

            /* Assigning centers to array */
            for (j = 0; j < 4; j++) {
                PTspotsCen[j]      = Sub_Ap_spots1[j]; /* 0 sub-aperture */
                PTspotsCen[j + 16] = Sub_Ap_spots1[j + 4];

                PTspotsCen[j + 4]  = Sub_Ap_spots2[j]; /* 1 sub-aperture */
                PTspotsCen[j + 20] = Sub_Ap_spots2[j + 4];

                if (Pa == 1)
                    if (j == 3) {
                        PredictedSpots[j]      = 1;
                        PredictedSpots[j + 16] = 1;
                    }
                if (Pb == 1)
                    if (j == 3) {
                        PredictedSpots[j + 4]  = 1;
                        PredictedSpots[j + 20] = 1;
                    }

            } /*for*/
        } else if (Right)  {
            if (PRINTMSG) printf("Executing-->Pattern: Vertical and Right \n");
            MissingSub_Apertures[0] = 1;
            MissingSub_Apertures[1] = 2;

            for (j = 0; j < 4; j++) {
                PTspotsCen[j + 12] = Sub_Ap_spots1[j];
                PTspotsCen[j + 28] = Sub_Ap_spots1[j + 4]; /* 3 sub-aperture */

                PTspotsCen[j + 8]  = Sub_Ap_spots2[j]; /* 2 sub-aperture */
                PTspotsCen[j + 24] = Sub_Ap_spots2[j + 4];

                if (Pa == 1)
                    if (j == 3) {
                        PredictedSpots[j + 12] = 1;
                        PredictedSpots[j + 28] = 1;
                    }
                if (Pb == 1)
                    if (j == 3) {
                        PredictedSpots[j + 8]  = 1;
                        PredictedSpots[j + 24] = 1;
                    }

            } /* for*/
        }     /* Right */
    } else if ((Horizental == 1) && (Vertical == 0))         {
        if (Bottom) {
            if (PRINTMSG) printf("Executing-->Pattern: Horizental Bottom \n");
            MissingSub_Apertures[1] = 2;
            MissingSub_Apertures[2] = 3;

            for (j = 0; j < 4; j++) {
                PTspotsCen[j]      = Sub_Ap_spots1[j];
                PTspotsCen[j + 16] = Sub_Ap_spots1[j + 4]; /* 0 sub-aperture */

                PTspotsCen[j + 12] = Sub_Ap_spots2[j]; /* 3rd sub-aperture */
                PTspotsCen[j + 28] = Sub_Ap_spots2[j + 4];

                if (Pa == 1)
                    if (j == 3) {
                        PredictedSpots[j]      = 1;
                        PredictedSpots[j + 16] = 1;
                    }
                if (Pb == 1)
                    if (j == 3) {
                        PredictedSpots[j + 12] = 1;
                        PredictedSpots[j + 28] = 1;
                    }

            } /* for */
        } else if (Top)   {
            if (PRINTMSG) printf("Executing--> Pattern: Horizental Top\n");
            MissingSub_Apertures[0] = 1;
            MissingSub_Apertures[3] = 4;

            for (j = 0; j < 4; j++) {
                PTspotsCen[j + 4]  = Sub_Ap_spots1[j]; /* 1st sub-aperture */
                PTspotsCen[j + 20] = Sub_Ap_spots1[j + 4];

                PTspotsCen[j + 8]  = Sub_Ap_spots2[j]; /* 2rd sub-aperture */
                PTspotsCen[j + 24] = Sub_Ap_spots2[j + 4];

                if (Pa == 1)
                    if (j == 3) {
                        PredictedSpots[j + 4]  = 1;
                        PredictedSpots[j + 20] = 1;
                    }
                if (Pb == 1)
                    if (j == 3) {
                        PredictedSpots[j + 8]  = 1;
                        PredictedSpots[j + 24] = 1;
                    }

            }
        } /*End of Top */
    } else {
        cpl_apertures_delete(ptspots8_a);
        cpl_apertures_delete(ptspots8_b);
        cpl_apertures_delete(ptspots8_c);
        cpl_apertures_delete(ptspots8_d);
        cpl_vector_unwrap(sig);

        return(-3);
    }

    for (i = 0; i < 4; i++) {
      for (j = 0; j < 4; j++) {
        if (PTspotsCen[j+4*i]!= 0) 
            PTAlgo[(32+2*(j+i*4))*4+nTelAlgo]=PTspotsCen[j+4*i]-(RefWindowPosPT[nTelAlgo]-PTImageWindowSize/2);
        if (PTspotsCen[j+4*i+16]!=0) 
            PTAlgo[(33+2*(j+i*4))*4+nTelAlgo]=PTspotsCen[j+4*i+16]-(RefWindowPosPT[nTelAlgo+4]-PTImageWindowSize/2);
      }
    }

    cpl_apertures_delete(ptspots8_a);
    cpl_apertures_delete(ptspots8_b);
    cpl_apertures_delete(ptspots8_c);
    cpl_apertures_delete(ptspots8_d);
    cpl_vector_unwrap(sig);
    return(error);
}


/*
 * This function is called when three or four lenslet sub-apertures found
 * in the pupil tracker image.
 */
cpl_error_code gvacqPupilTrackerSpotsDetectionThreeORFourSubs(cpl_image *DetPointer,
                                                              cpl_apertures *PTspots,
                                                              int llx, int lly,
                                                              int urx, int ury,
                                                              double *OneTelSpot5,
                                                              cpl_vector *sig,
                                                              double *PTspotsCen,
                                                              int *MissingSub_Apertures,
                                                              int *PredictedSpots) {
    /* test entries */
    if ((DetPointer == NULL) || (PTspots == NULL) ||
        (llx <= 0) || (lly <= 0) || (urx <= 0) || (ury <= 0)) {
        printf("ERROR::The input parameters are NULL. Occured at %s:%s:%d\n",
               __FILE__, __FUNCTION__, __LINE__);
        return(-1);
    }

    /* initialise  */

#ifdef VLT2014
    cpl_size pos;
#else
    int pos;
#endif

    double         DataDouble[4 * 2];
    int            Nspots, nspots;
    int            i, j, k;
    int            offset = 28;              /* to devide the pupil tracker image into 4 parts  */
    int            Sub_aperture_window = 26; /* This size image is used to detect the spots */
    double         Cpx, Cpy;
    double         px, py;
    cpl_error_code error;
    int            Pattern_left = 0, Pattern_Right = 0, Pattern_bottom = 0, Pattern_Top = 0;

    int llx_s, lly_s, urx_s, ury_s;

    for (i = 0; i < 4; i++) {
        MissingSub_Apertures[i] = 0;
    }
    for (i = 0; i < 32; i++) {
        PTspotsCen[i] = 0.0;
    }


    Nspots = cpl_apertures_get_size(PTspots);
    error  = cpl_apertures_sort_by_flux(PTspots);
    if (Nspots > 16) Nspots = 16;
    /*  Nspots=Nspots-2; left from Narsi. He assumes that the two last spots with less flux might be wrong. */

    /* --------------------------------------------- */
    double PTspots_double[Nspots * 2];
    double PTspots_double_x[Nspots];
    double PTspots_double_y[Nspots];

    error = cpl_apertures2XYdoubleArray(PTspots, Nspots, PTspots_double);

    for (k = 0; k < Nspots; k++) {
        PTspots_double_x[k] = PTspots_double[k] + llx;
        PTspots_double_y[k] = PTspots_double[k + Nspots] + lly;
    }

    cpl_vector *Vec_x = cpl_vector_wrap(Nspots, (double *)PTspots_double_x);
    cpl_vector *Vec_y = cpl_vector_wrap(Nspots, (double *)PTspots_double_y);

#ifdef VLT2014
    error = cpl_vector_sort(Vec_x, CPL_SORT_ASCENDING);
    error = cpl_vector_sort(Vec_y, CPL_SORT_ASCENDING);
#else
    error = cpl_vector_sort(Vec_x, 1);
    error = cpl_vector_sort(Vec_y, 1);
#endif

    double PTspots_avg_x = cpl_vector_get_mean(Vec_x);
    double PTspots_avg_y = cpl_vector_get_mean(Vec_y);


    for (k = 0; k < Nspots; k++) {
        if (cpl_vector_get(Vec_x, k) < PTspots_avg_x)
            Pattern_left++;
        else if (cpl_vector_get(Vec_x, k) > PTspots_avg_x)
            Pattern_Right++;
        if (cpl_vector_get(Vec_y, k) < PTspots_avg_y)
            Pattern_bottom++;
        else if (cpl_vector_get(Vec_y, k) > PTspots_avg_y)
            Pattern_Top++;
    }

    int pattern_diff_x = Pattern_Right - Pattern_left;
    int pattern_diff_y = Pattern_Top - Pattern_bottom;

/*
 *  cpl_vector_dump(Vec_x, NULL);
 *  cpl_vector_dump(Vec_y, NULL);
 */
    cpl_vector_unwrap(Vec_x);
    cpl_vector_unwrap(Vec_y);

    Cpx = 0.0;
    Cpy = 0.0;


    if (Nspots < 16) {
        for (j = 0; j < Nspots; j++) {
            px = (double)cpl_apertures_get_centroid_x(PTspots, j + 1) + llx - 1;
            py = (double)cpl_apertures_get_centroid_y(PTspots, j + 1) + lly - 1;

            Cpx += px;
            Cpy += py;
        }
        Cpx = Cpx / Nspots;
        Cpy = Cpy / Nspots;


        if ((Nspots < 14) && (fabs(Pattern_Right - Pattern_left) > 2) && (fabs(Pattern_Top - Pattern_bottom) > 2)) {
            if (PRINTMSG) printf("--->Pattern: L -----\n");
            if ((Pattern_left < Nspots / 2) && (fabs(pattern_diff_x) >= 2))
                Cpx = Cpx - 10;
            else if ((Pattern_left > Nspots / 2) && (fabs(pattern_diff_x) >= 2))
                Cpx = Cpx + 10;
            if ((Pattern_bottom < Nspots / 2) && (fabs(pattern_diff_y) >= 2))
                Cpy = Cpy - 10;
            else if ((Pattern_bottom > Nspots / 2) && (fabs(pattern_diff_y) >= 2))
                Cpy = Cpy + 10;
        } else if ((Nspots < 14) && (fabs(fabs(Pattern_Right - Pattern_left) - fabs(Pattern_Top - Pattern_bottom)) >= 0))            {
            if (PRINTMSG) printf("--->Pattern: ------ Not  L ----\n");
            if ((pattern_diff_x > 0) && (fabs(pattern_diff_y) < 3))
                Cpx = Cpx - 10;
            else if ((pattern_diff_x < 0) && (fabs(pattern_diff_y) < 3))
                Cpx = Cpx + 10;
            if ((pattern_diff_y > 0) && (fabs(pattern_diff_x) < 3))
                Cpy = Cpy - 10;
            else if ((pattern_diff_y < 0) && (fabs(pattern_diff_x) < 3))
                Cpy = Cpy + 10;
        } /*Not L */
    } else { /* Nspots >= 16 */
        for (j = 0; j < 16; j++) {
            px = (double)cpl_apertures_get_centroid_x(PTspots, j + 1) + llx - 1;
            py = (double)cpl_apertures_get_centroid_y(PTspots, j + 1) + lly - 1;

            Cpx += px;
            Cpy += py;
        }
        Cpx = Cpx / 16.0;
        Cpy = Cpy / 16.0;
    }

    /* storing the division of the different subapertures for debug */
    PTAlgo[30*4+nTelAlgo]=Cpx-(RefWindowPosPT[nTelAlgo]-PTImageWindowSize/2);
    PTAlgo[31*4+nTelAlgo]=Cpy-(RefWindowPosPT[nTelAlgo+4]-PTImageWindowSize/2);

    for (i = 0; i < 4; i++) { /* 4 sub-apertures */
        /* STEP1: Prepare to extract sub-aperture images */
        /* here llx lly urx ury are copied for sub images */
        if (i == 0) { /* left bottom corner */
            llx_s = (int)Cpx - offset - Sub_aperture_window;
            urx_s = (int)Cpx - offset + Sub_aperture_window;
            lly_s = (int)Cpy - offset - Sub_aperture_window;
            ury_s = (int)Cpy - offset + Sub_aperture_window;
        }

        if (i == 1) { /* left up  */
            llx_s = (int)Cpx - offset - Sub_aperture_window;
            urx_s = (int)Cpx - offset + Sub_aperture_window;
            lly_s = (int)Cpy + offset - Sub_aperture_window;
            ury_s = (int)Cpy + offset + Sub_aperture_window;
        }

        if (i == 2) { /* right up  */
            llx_s = (int)Cpx + offset - Sub_aperture_window;
            urx_s = (int)Cpx + offset + Sub_aperture_window;
            lly_s = (int)Cpy + offset - Sub_aperture_window;
            ury_s = (int)Cpy + offset + Sub_aperture_window;
        }

        if (i == 3) { /* right bottom  */
            llx_s = (int)Cpx + offset - Sub_aperture_window;
            urx_s = (int)Cpx + offset + Sub_aperture_window;
            lly_s = (int)Cpy - offset - Sub_aperture_window;
            ury_s = (int)Cpy - offset + Sub_aperture_window;
        }

        /* STEP2: Extract sub-aperture images */
        cpl_image *ssimage = cpl_image_extract(DetPointer, llx_s, lly_s, urx_s, ury_s);

        int err = 0;
        if (ssimage == NULL) {
            printf("ERROR in extracting image %s:%s:%d\n", __FILE__, __FUNCTION__, __LINE__);
            err = -1;
        } else {
            if ((double)cpl_image_get_size_x(ssimage) * RON > cpl_image_get_absflux(ssimage)) {
                if (PRINTMSG) printf("No flux in sub-image %f > %f\n",
                                     (double)cpl_image_get_size_x(ssimage) * cpl_image_get_size_y(ssimage) * RON, cpl_image_get_absflux(ssimage));
                err = -1;
            } else err = Check_objects_found_PT(ssimage, FIfwhm_minthreshold, 10.0);
        }

        if (err != 0) {
            MissingSub_Apertures[i] = i + 1;
            if (PRINTMSG) printf("-----> Missing subaperture in PT = %d\n", MissingSub_Apertures[i]);
        } else {
            cpl_apertures *ptspots = cpl_apertures_extract_window(DetPointer,
                                                                  sig,
                                                                  llx_s, lly_s,
                                                                  urx_s, ury_s,
                                                                  &pos);

            if (ptspots == NULL) {
                MissingSub_Apertures[i] = i + 1;
                if (PRINTMSG) printf("-----> ptspots==NULL, Missing subaperture = %d\n",
                                     MissingSub_Apertures[i]);
            } else { /* vilability of spots */
                /* STEP 9:: Sort them via flux in order to select only bright 4 objects */
                error = cpl_apertures_sort_by_flux(ptspots);

                nspots = cpl_apertures_get_size(ptspots);

                PTAlgo[4 * (17 + i) + nTelAlgo] = nspots;

                if (PRINTMSG) { if (nspots != 4) printf("Sub-aperture: %d, spots %d \n", i + 1, nspots); else printf("\n"); }

                if ((nspots > 2) && (nspots <= 5)) {
                    if (nspots > 4) nspots = 4;
                    error = cpl_apertures_sort_by_flux(ptspots);

                    int wrongspot = 0;
                    for (k = 0; k < nspots; k++) {
                        if ((cpl_apertures_get_max(ptspots, k + 1) / cpl_apertures_get_max(ptspots, 1)) < 1 / 25.0) wrongspot++;
                    } /* end of for k=0; k<nspots; k++ */


                    if ((nspots >= 4) && (wrongspot <= 2)) /*  TOP level */
                        error = cpl_apertures2XYdoubleArray(ptspots, 4, DataDouble);
                    else { /*Predict 4th spot */
                        if ((nspots == 3) && (wrongspot == 0)) {
                            error = Predict4thSpotWithThree(ptspots, DataDouble);
                            for (k = 3; k < 4; k++) {
                                PredictedSpots[k + i * 4]      = 1;
                                PredictedSpots[16 + k + i * 4] = 1;
                            }
                        } else {
                            MissingSub_Apertures[i] = i + 1;
                            for (k = 0; k < 4; k++) {
                                DataDouble[k] = (llx_s + urx_s) / 2.0;
                                DataDouble[k] = (lly_s + ury_s) / 2.0;
                            }
                        }
                    } /* End of else if (nspots <4 ) */

                    cpl_apertures_delete(ptspots);
                } else {
                    MissingSub_Apertures[i] = i + 1;
                    if (PRINTMSG) printf("-----> Missing subaperture = %d\n", MissingSub_Apertures[i]);
                    cpl_apertures_delete(ptspots);
                }
            } /* if else spots-apertures exist */
        }     /* if else objects found */
        cpl_image_delete(ssimage);

        /* STEP 10:: fill the array with centers to return */
        for (k = 0; k < 4; k++) {
            PTspotsCen[k + i * 4]      = DataDouble[k];
            PTspotsCen[16 + k + i * 4] = DataDouble[k + 4];

            if (DataDouble[k]!=0) 
                PTAlgo[(32+2*(k+i*4))*4+nTelAlgo]=DataDouble[k]-(RefWindowPosPT[nTelAlgo]-PTImageWindowSize/2);
            if (DataDouble[k+4]!=0) 
                PTAlgo[(33+2*(k+i*4))*4+nTelAlgo]=DataDouble[k + 4]-(RefWindowPosPT[nTelAlgo+4]-PTImageWindowSize/2);
        }

    } /*END OF i or sub-aperture*/


    return(error);
}


/* ----------------------------------------------------------------------*/

/**
 * @internal
 * @brief  Detects the 1 telescope pupil tracker spots positions
 * @param  DetPointer   pointer to a  detector image
 * @param llx	Lower left x position (FITS convention, 1 for leftmost)
 * @param lly	Lower left y position (FITS convention, 1 for lowest)
 * @param urx	Upper right x position (FITS convention)
 * @param ury	Upper right y position (FITS convention)
 * @param PTspotsCen pupil tracker 16 spots positions [1Tel.*16Spots*2(x,y)]
 *
 * Requirement: Pupil tracker spots are under rotation (due to K-mirror), we don't know its positions at particular instance.
 * That requires the spots to be detected automatically for any pupil displacement and rotation without using
 * guess positions.
 *
 * Implementation:
 * 1. Find the center of pupil tracker spots C(x,y).
 * For this pupil tracker window is extracted from the instrument database
 * with input locations and scanned for spots with input sigma.
 * It will be searched for 16 spots. The detected spots does not tell us that correspond to particular lens-subaperture.
 * The avg of spot postions gives the C(x,y).
 *
 * 2. Using C(x,y), input pupil tracker is devided into 4 parts centered at C(x,y).
 * In each sub-image (corner image), existing four spots are searched and then return to top level function.
 * @return error
 */
cpl_error_code gvacqPupilTrackerDetectedCentroidsOneTel(cpl_image *DetPointer,
                                                        int llx, int lly, int urx, int ury,
                                                        double *OneTelSpot5,
                                                        double *PTspotsCen,
                                                        int *MissingSub_Apertures,
                                                        int *PredictedSpots) {
    /* test entries */
    if ((DetPointer == NULL) || (PTspotsCen == NULL) || (llx <= 0) || (lly <= 0) || (urx <= 0) || (ury <= 0)) {
        printf("ERROR::The input parameters are NULL. Occured at %s:%s:%d\n",
               __FILE__, __FUNCTION__, __LINE__);
        return(-1);
    }

    /* initialise  */

#ifdef VLT2014
    cpl_size pos;
#else
    int pos;
#endif


    int            Nspots;
    int            i, j;
    cpl_apertures  *PTspots;
    cpl_error_code error;
/*    double    psig[4]; */
    double PT_spots_scan_sigma[4];
    int    Maximum_check   = 8;
    double sigma_step_grow = 0.3;
    int    pis_rejected    = 0;

    for (i = 0; i < 4; i++) {
        MissingSub_Apertures[i] = 0;
    }
    for (i = 0; i < 32; i++) {
        PTspotsCen[i] = 0.0;
    }

    /* --------STEP 1:: Select the particular telescope image ---------*/
    cpl_image *PT_subimage = cpl_image_extract(DetPointer,
                                               llx, lly,
                                               urx, ury);

    /* low pass filtering  */
    cpl_matrix *Ker = cpl_matrix_new(1, 1);
    cpl_matrix_fill(Ker, 1.0);
    cpl_image *PT_subimage_filtered = cpl_image_filter_median(PT_subimage, Ker);
    cpl_matrix_delete(Ker);
    cpl_image_delete(PT_subimage);



    /* ------ STEP 2: Check if the image has  spots data ------------ */
    int size_x = cpl_image_get_size_x(PT_subimage_filtered);


    if (size_x * RON > cpl_image_get_absflux(PT_subimage_filtered)) {
        if (PRINTMSG) printf("ERROR:: Please make sure pupil guiding lasers are turned ON, %s:%s:%d \n",
                             __FILE__, __FUNCTION__, __LINE__);
        cpl_image_delete(PT_subimage_filtered);

        return(-2);
    }


    error = Check_objects_found(PT_subimage_filtered, FIfwhm_minthreshold);
    if ((error == -1) || (error == -2)) {
        if (error == -1)
            printf("The inputs are NULL %s:%d \n", __FUNCTION__, __LINE__);
        else if (PRINTMSG)
            printf("ERROR:: Pupil guiding lasers spots are NOT FOUND \n");

        cpl_image_delete(PT_subimage_filtered);
        return(-2);
    }

    /*copy pt spots scanning sigmas to local variables */
    for (i = 0; i <= 3; i++) {
        PT_spots_scan_sigma[i] = PT_spots_scan_sigmaOrig[i];
/*      psig[i]=PT_spots_scan_sigma[i]/3.0; */
    }
    cpl_vector *sigmas = cpl_vector_wrap(4, (double *)PT_spots_scan_sigma);


    /*-------- STEP 3:: Extract the spots in the pupil tracker image-------- */
    PTspots = cpl_apertures_extract(PT_subimage_filtered, sigmas, &pos);

    if (PTspots == NULL) {
        if (PRINTMSG)
            printf("ERROR:: Please make sure pupil guiding lasers are \
turned ON or make it fall under PT window; %s:%s:%d \n",
                   __FILE__, __FUNCTION__, __LINE__);
        cpl_image_delete(PT_subimage_filtered);
        cpl_vector_unwrap(sigmas);
        cpl_apertures_delete(PTspots);
        return(-2);
    }

    cpl_vector *sigmas_dup = cpl_vector_duplicate(sigmas);

    Nspots           = cpl_apertures_get_size(PTspots);
    PTAlgo[nTelAlgo] = Nspots;

    error = cpl_apertures_sort_by_flux(PTspots);

    /*------------ STEP 4: if Nspots>16, select only 16 spots by increasing sigma ------------- */
    if (Nspots > 16) {
        i = 0;
        do {
            error = cpl_vector_add_scalar(sigmas, sigma_step_grow);
            cpl_apertures_delete(PTspots);
            PTspots = cpl_apertures_extract(PT_subimage_filtered, sigmas, &pos);
            Nspots  = cpl_apertures_get_size(PTspots);
            i++;
        } while (Nspots > 16 && i < Maximum_check); /* end of do */
    }

    PTAlgo[4 + nTelAlgo] = Nspots;

    /*------------ STEP 5: if Nspots<16, select only 16 spots by increasing or decreasing sigma ------------- */
    /* by A. Amorim Aug 2016 */
    if ((Nspots < 16) && (Nspots >= 6)) {
        /*------------ see if it the diffraction of two aligned spots creates a connection between the two in the mask - */
        /*------------ try to increase the threshold to separate the two spots in the mask  ------------- */
        int origNspots = Nspots;
        i = 0;
        do {
            error = cpl_vector_add_scalar(sigmas, 2 * sigma_step_grow);
            cpl_apertures_delete(PTspots);
            PTspots = cpl_apertures_extract(PT_subimage_filtered, sigmas, &pos);
            Nspots  = cpl_apertures_get_size(PTspots);
            i++;
        } while (Nspots < 16 && i < Maximum_check); /* end of do */

        /* if failed revert changes and try decreasing */
        if ((Nspots <= origNspots) || (minApertureDistance(PTspots) < 6)) {
            for (i = 0; i < 4; i++) {
                cpl_vector_set(sigmas, i, PT_spots_scan_sigmaOrig[i]);
            }
            cpl_apertures_delete(PTspots);
            PTspots = cpl_apertures_extract(PT_subimage_filtered, sigmas, &pos);
            Nspots  = cpl_apertures_get_size(PTspots);

            /*--------- now try by decreasing siga and see if some spot was bellow the threshold  ------ */
            i = 0;
            do {
                error = cpl_vector_add_scalar(sigmas, -sigma_step_grow);
                cpl_apertures_delete(PTspots);
                PTspots = cpl_apertures_extract(PT_subimage_filtered, sigmas, &pos);
                Nspots  = cpl_apertures_get_size(PTspots);
                i++;
            } while (Nspots < 16 && i < Maximum_check); /* end of do */
            /* if failed revert changes*/
            if ((Nspots <= origNspots) || (minApertureDistance(PTspots) < 6)) {
                for (i = 0; i < 4; i++) {
                    cpl_vector_set(sigmas, i, PT_spots_scan_sigmaOrig[i]);
                }
                cpl_apertures_delete(PTspots);
                PTspots = cpl_apertures_extract(PT_subimage_filtered, sigmas, &pos);
                Nspots  = cpl_apertures_get_size(PTspots);
            }
        }
    }

    PTAlgo[2 * 4 + nTelAlgo] = Nspots;

    /*------ STEP 6:: still spots are < 13, use the original sigma to find the Nspots------  */
    if (Nspots < 8) {
        cpl_apertures_delete(PTspots);
        PTspots = cpl_apertures_extract(PT_subimage_filtered, sigmas_dup, &pos);
        Nspots  = cpl_apertures_get_size(PTspots);
    }
    cpl_vector_delete(sigmas_dup);


    /*-------STEP 6:: Sort them via flux in order to select only bright 16 objects----- */
    error = cpl_apertures_sort_by_flux(PTspots); /* Sort the objects by flux */


    /* debug the 2 spots with less flux before reducing using spotsInApertureMultiplicity */
    if (Nspots >2 ){
      PTAlgo[28*4+nTelAlgo]=cpl_apertures_get_centroid_x(PTspots, Nspots)+llx-(RefWindowPosPT[nTelAlgo]-PTImageWindowSize/2);
      PTAlgo[29*4+nTelAlgo]=cpl_apertures_get_centroid_y(PTspots, Nspots)+lly-(RefWindowPosPT[nTelAlgo+4]-PTImageWindowSize/2);
      PTAlgo[26*4+nTelAlgo]=cpl_apertures_get_centroid_x(PTspots, Nspots-1)+llx-(RefWindowPosPT[nTelAlgo]-PTImageWindowSize/2);
      PTAlgo[27*4+nTelAlgo]=cpl_apertures_get_centroid_y(PTspots, Nspots-1)+lly-(RefWindowPosPT[nTelAlgo+4]-PTImageWindowSize/2);
    }

    /* cuting one or two of the last apertures if they have multiplicity = 0 */
    Nspots=spotsInApertureMultiplicity(PTspots);

    PTAlgo[3 * 4 + nTelAlgo] = Nspots;

    for (i = 1; i <= urx - llx; i++) {
        for (j = 1; j <= ury - lly; j++) {
            error = cpl_image_set(DetPointer, llx + i, lly + j,
                                  cpl_image_get(PT_subimage_filtered, i, j, &pis_rejected));
        }
    }

    if (PRINTMSG) printf("NSPOTS after optimizing sigmas %d \n", Nspots);
    if ((Nspots > 6) && (Nspots <= 10)) {

        error = gvacqPupilTrackerSpotsDetectionTwoSubs(DetPointer,
                                                       PTspots,
                                                       llx, lly,
                                                       OneTelSpot5,
                                                       PTspotsCen,
                                                       MissingSub_Apertures,
                                                       PredictedSpots);
        PTAlgo[4 * 4 + nTelAlgo] = PTAlgo[4 * 4 + nTelAlgo] | 1;

    } else if ((Nspots > 10) && (Nspots <= 20)) { /* the spots with less flux will be dropped anyhow */

        cpl_vector *sig = cpl_vector_duplicate(sigmas);
        error = cpl_vector_multiply_scalar(sig, 1.0 / 3.0);
        error = gvacqPupilTrackerSpotsDetectionThreeORFourSubs(DetPointer,
                                                               PTspots,
                                                               llx, lly, urx, ury,
                                                               OneTelSpot5,
                                                               sig,
                                                               PTspotsCen,
                                                               MissingSub_Apertures,
                                                               PredictedSpots);
        cpl_vector_unwrap(sig);

    } else error = -3;

    /* Free the buffers */
    cpl_vector_unwrap(sigmas);
    cpl_apertures_delete(PTspots);
    cpl_image_delete(PT_subimage_filtered);
    return(error);
}


int Get_average_position_apertures(cpl_apertures *aperture, double *XY) {
    int    j, nx;
    double px, py, cpx, cpy;

    cpx = 0.0;
    cpy = 0.0;

    if (aperture == NULL) return(-1);

    nx = cpl_apertures_get_size(aperture);
    for (j = 0; j < nx; j++) {
        px = (double)cpl_apertures_get_centroid_x(aperture, j + 1);
        py = (double)cpl_apertures_get_centroid_y(aperture, j + 1);

        cpx += px;
        cpy += py;
    }

    cpx = cpx / nx;
    cpy = cpy / nx;

    XY[0] = cpx;
    XY[1] = cpy;

    return(0);
}


int Maximum(double a[], int n) {
    int i, max, index;

    max   = a[0];
    index = 0;

    for (i = 1; i < n; i++) {
        if (a[i] > max) {
            index = i;
            max   = a[i];
        }
    }

    return(index);
}


cpl_error_code Predict4thSpotWithThree(cpl_apertures *ptspots1,
                                       double        *Data4Spots) {
    int    i = 0, j = 0, nx = 0, a = 0, b = 0, c = 0, d = 0;
    double DataDouble[8];
    double px = 0, py = 0, error = 0;

    double cpx = 0.0;
    double cpy = 0.0;

    nx = cpl_apertures_get_size(ptspots1);

    if (nx < 2) error = -1;
    double SpotsThree[6];


/*  Find the spots centers*/
    for (j = 0; j < nx; j++) {
        px = (double)cpl_apertures_get_centroid_x(ptspots1, j + 1);
        py = (double)cpl_apertures_get_centroid_y(ptspots1, j + 1);

        SpotsThree[j]      = px;
        SpotsThree[j + nx] = py;

        cpx += px;
        cpy += py;
    }     /* j finished the positions of spots */

    cpx = cpx / nx;
    cpy = cpy / nx;

    double R[3]; /* radius */
    R[0] = sqrt(fabs(SpotsThree[0] - SpotsThree[1]) * fabs(SpotsThree[0] - SpotsThree[1]) +
                fabs(SpotsThree[3 + 0] - SpotsThree[3 + 1]) * fabs(SpotsThree[3 + 0] - SpotsThree[3 + 1]));

    R[1] = sqrt(fabs(SpotsThree[0] - SpotsThree[2]) * fabs(SpotsThree[0] - SpotsThree[2]) +
                fabs(SpotsThree[3 + 0] - SpotsThree[3 + 2]) * fabs(SpotsThree[3 + 0] - SpotsThree[3 + 2]));

    R[2] = sqrt(fabs(SpotsThree[1] - SpotsThree[2]) * fabs(SpotsThree[1] - SpotsThree[2]) +
                fabs(SpotsThree[3 + 1] - SpotsThree[3 + 2]) * fabs(SpotsThree[3 + 1] - SpotsThree[3 + 2]));

    int i_1, i_2;
    if (Maximum(R, 3) == 0) {
        i_1 = 0;
        i_2 = 1;
    } else if (Maximum(R, 3) == 1)     {
        i_1 = 0;
        i_2 = 2;
    } else if (Maximum(R, 3) == 2)     {
        i_1 = 1;
        i_2 = 2;
    }

    double diagonal[2];
    diagonal[0] = (SpotsThree[i_1] + SpotsThree[i_2]) / 2.0;
    diagonal[1] = (SpotsThree[i_1 + 3] + SpotsThree[i_2 + 3]) / 2.0;

    for (i = 0; i < 3; i++) {
        Data4Spots[i]     = SpotsThree[i]; /* copying good spots*/
        Data4Spots[i + 4] = SpotsThree[i + 3];
    }

    /*--------------------------------------- */
    /* -------  Step C: Spot manual guess --- */
    /* -------------------------------------- */
    /* which spot is missed ?*/
    for (j = 0; j < nx; j++) {
        px = (double)SpotsThree[j];
        py = (double)SpotsThree[j + nx];

        if ((cpx > px) && (cpy > py)) { DataDouble[0] = px; DataDouble[4] = py; a = 1; }
        if ((cpx > px) && (cpy < py)) { DataDouble[1] = px; DataDouble[5] = py; b = 1; }
        if ((cpx < px) && (cpy < py)) { DataDouble[2] = px; DataDouble[6] = py; c = 1; }
        if ((cpx < px) && (cpy > py)) { DataDouble[3] = px; DataDouble[7] = py; d = 1; }
    }    /* End of which spot is missed ? */


    /* PREDICT th location of a missing spot using existing spots info*/
    if (a + b + c + d == 3) {
/*
 * -------------------------------------------------
 |                          |                    |
 |  b=DataDouble(1,5)       |  c=DataDouble(2,6) |
 |                          |
 | -----------------------(cxx,cyy)----------------|
 |                          |                    |
 |  a=DataDouble(0,4)       |  d=DataDouble(3,7) |
 | -------------------------------------------------
 */
        double cxx, cyy;
        if ((a == 0) || (c == 0)) {
            cxx = (DataDouble[1] + DataDouble[3]) / 2.0; /* (cxx, cyy)= (b+d)/2 */
            cyy = (DataDouble[5] + DataDouble[7]) / 2.0;
        } else if ((b == 0) || (d == 0))           {
            cxx = (DataDouble[0] + DataDouble[2]) / 2.0;
            cyy = (DataDouble[4] + DataDouble[6]) / 2.0;
        }

        if (a == 0) {
            Data4Spots[3] = cxx + (cxx - DataDouble[2]);
            Data4Spots[7] = cyy + (cyy - DataDouble[6]);
        }

        if (b == 0) {
            Data4Spots[3] = cxx + (cxx - DataDouble[3]);
            Data4Spots[7] = cyy + (cyy - DataDouble[7]);
        }

        if (c == 0) {
            Data4Spots[3] = cxx + (cxx - DataDouble[0]);
            Data4Spots[7] = cyy + (cyy - DataDouble[4]);
        }

        if (d == 0) {
            Data4Spots[3] = cxx + (cxx - DataDouble[1]);
            Data4Spots[7] = cyy + (cyy - DataDouble[5]);
        }
    }


    if ((diagonal[0] - cpx) < 0.5)
        Data4Spots[3] = diagonal[0] + 3 * (diagonal[0] - cpx);
    else if (diagonal[0] - cpx > 0.5)
        Data4Spots[3] = diagonal[0] + 3 * (diagonal[0] - cpx);

    if ((diagonal[1] - cpy) < 0.5)
        Data4Spots[7] = diagonal[1] + 3 * (diagonal[1] - cpy);
    else if (diagonal[1] - cpy > 0.5)
        Data4Spots[7] = diagonal[1] + 3 * (diagonal[1] - cpy);
    if (PRINTMSG) printf("Predicted spot: %f, %f \n", Data4Spots[3], Data4Spots[7]);

    for (i = 0; i < 3; i++) {
        Data4Spots[i]     = SpotsThree[i]; /* copying good spots*/
        Data4Spots[i + 4] = SpotsThree[i + 3];
    }
    return(error);
}


/* ----------------------------------------------------------------------*/

/**
 * @internal
 * @brief  Determines the reference spots (5th spot refer to Fiber Coupler) from
 *         Guess grid positions for all 4 Telescope beams at a time
 * @param  image  pointer to the  detector image
 * @param  GuessCentre GuessCentre positions (first moment of 4 spots) [4*2]
 * @param  reference_spot (output) reference spots (5th) positions ([S0x,S1x,S2x,S3x,S0y,S1y,S2y,S3y,S4y])
 * @param  reference_spot_error (output)  reference spots positions error
 * @return status error
 * @see   gvacqPupilTracker
 *
 *
 * @par Error Handling:
 * CPL_ERROR_NULL_INPUT if(one of ) the input pointer (s) is NULL
 * CPL_ERROR_NULL
 */
cpl_error_code gvacqPupilTrackerFCspots(cpl_image *image,
                                        double    *GuessCentre, /* Xc[0-7], Yc[4-15] */
                                        double    *FCspots,     /* Xcm[0-7], Ycm[4-15] */
                                        double    *FCspots_error,
                                        int       *PTCalibErrorStatus) {
    /* test entries */

    if ((image == NULL) || (GuessCentre == NULL) || (FCspots == NULL) || (FCspots_error == NULL)) {
        printf("The input image buffer or parameters  NULL at func FCspots. Occured at %s:%s:%d\n",
               __FILE__, __FUNCTION__, __LINE__);
        return(-1);
    }

/*
 * cpl_ensure_code(image   || GuessCentre  , CPL_ERROR_NULL_INPUT);
 * cpl_ensure_code (cpl_image_get_size_x(image) ==IMAGE_X ||
 * cpl_image_get_size_y(image) ==IMAGE_Y,
 * CPL_ERROR_ILLEGAL_INPUT);
 */
    cpl_error_code error = CPL_ERROR_NONE;
    int            i, j, k;
    int            Npos = 2 * 4; /* 4 FC spots * 2(x,y) */
    int            NTel = 4;
    double         tmp_centers[8];
    double         tmp_measured_spots[Npos], tmp_measured_spots_error[Npos];
    double         DataX[16];
    double         DataY[16];


    for (i = 0; i < 16; i++) {
        DataX[i] = GuessCentre[i];
        DataY[i] = GuessCentre[16 + i];
    }

    /* 4 tel.  */
    for (i = 0; i < NTel; i++) {
        for (j = 0; j < 4; j++) {
            tmp_centers[j]     = DataX[j + i * 4];
            tmp_centers[j + 4] = DataY[j + i * 4];

            /* printf("%d %f %f %f %f \n", i, DataX[i*4+j], DataY[i*4+j],
             * tmp_centers[j], tmp_centers[j+4]); */
        }
        error = gvacqPupilTrackerFCReferenceSpots(image, tmp_centers,     /* Xc[0-3], Yc[4-7] */
                                                  tmp_measured_spots,     /* Xcm[0-3], Ycm[4-7] */
                                                  tmp_measured_spots_error);
        PTCalibErrorStatus[i] = error;
        if (error != 0) printf("T%d FC spots computing error, at %s:%s:%d\n", i + 1,
                               __FILE__, __FUNCTION__, __LINE__);
        for (k = 0; k < 4; k++) {
            FCspots[k + i * 4]      = tmp_measured_spots[k];
            FCspots[k + i * 4 + 16] = tmp_measured_spots[k + 4];

            FCspots_error[k + i * 4]      = tmp_measured_spots_error[k];
            FCspots_error[k + i * 4 + 16] = tmp_measured_spots_error[k + 4];
        }
        if (error != CPL_ERROR_NONE)
            cpl_msg_error(cpl_func, "Error found at %s:%s:%d\n",
                          __FILE__, __FUNCTION__, __LINE__);
    }    /* i  tel */

    /* for(i=0; i<8; i++) printf("%f %f %f %f \n", FCspots[i], FCspots[i+8], FCspots[i+16], FCspots[i+24]); */
    return(error);
} /* end of gvacqPupilTrackerFCReferenceSpots  */


/* ----------------------------------------------------------------------*/

/**
 * @internal
 * @brief  Determines the reference spots (5th spot refer to Fiber Coupler) from
 *         Guess grid positions for one Telescope beam at a time
 * @param  image  pointer to the  detector image
 * @param  GuessCentre GuessCentre positions (first moment of 4 spots) [4*2]
 * @param  reference_spot (output) reference spots (5th) positions ([S0x,S1x,S2x,S3x,S0y,S1y,S2y,S3y,S4y])
 * @param  reference_spot_error (output)  reference spots positions error
 * @return status error
 * @see   gvacqPupilTracker
 *
 * Functions pre-required:
 * setPTrefCoordinates
 *
 * @par Error Handling:
 * CPL_ERROR_NULL_INPUT if(one of ) the input pointer (s) is NULL
 * CPL_ERROR_NULL
 */
cpl_error_code gvacqPupilTrackerFCReferenceSpots(cpl_image *image,
                                                 double    *GuessCentre,    /* Xc[0-3], Yc[4-7] */
                                                 double    *reference_spot, /* Xcm[0-3], Ycm[4-7] */
                                                 double    *reference_spot_error) {
    /* test entries */

    if ((image == NULL) || (GuessCentre == NULL) || (reference_spot == NULL) ||
        (reference_spot_error == NULL)) {
        printf("The input parameters are NULL. Occured at %s:%s:%d\n",
               __FILE__, __FUNCTION__, __LINE__);
        return(-1);
    }

    cpl_error_code error = CPL_ERROR_NONE;
    int            i, j, k;
    int            robustness = 7;
    int            Nspots     = 4;

    cpl_matrix *SpotLocations      = cpl_matrix_new(2, Nspots);
    cpl_matrix *xy_centre          = cpl_matrix_new(2, Nspots);
    cpl_matrix *xy_centre_err      = cpl_matrix_new(2, Nspots);
    cpl_matrix *xy_sigma           = cpl_matrix_new(2, Nspots);
    cpl_matrix *xy_sigma_err       = cpl_matrix_new(2, Nspots);
    cpl_matrix *xy_fwhm            = cpl_matrix_new(2, Nspots);
    cpl_matrix *xy_fwhm_err        = cpl_matrix_new(2, Nspots);
    cpl_matrix *centre_intensities = cpl_matrix_new(1, Nspots);
    cpl_array  *all_error_codes    = cpl_array_new(Nspots, CPL_TYPE_INT);
    int        GwindowFit          = 10;

    for (i = 0; i < Nspots; i++) {
        cpl_matrix_set(SpotLocations, 0, i, GuessCentre[i]);          /* X-spots */
        cpl_matrix_set(SpotLocations, 1, i, GuessCentre[i + Nspots]); /* Y-spots  */
    }

    /* finding the centriod for the aberrated laser wavefront*/
    error = clipm_centroiding_multi_gauss(image, SpotLocations, GwindowFit,
                                          &xy_centre, &xy_centre_err, &xy_sigma, &xy_sigma_err, &xy_fwhm,
                                          &xy_fwhm_err, &centre_intensities, &all_error_codes, robustness);
    if (error != CPL_ERROR_NONE)
        cpl_msg_error(cpl_func, "Error found at %s:%s:%d\n", __FILE__,
                      __FUNCTION__, __LINE__);
         /* cpl_matrix_dump(xy_centre, NULL); */
    double *centroid       = cpl_matrix_get_data(xy_centre);
    double *centroid_error = cpl_matrix_get_data(xy_centre_err);

    for (k = 0; k < Nspots; k++) {
        for (j = 0; j < 2; j++) {
            reference_spot[k + j * Nspots]       = centroid[k + j * Nspots];       /*X=centroid[0,3]  */
            reference_spot[k + j * Nspots]       = centroid[k + j * Nspots];       /*X=centroid[4,7]   */
            reference_spot_error[k + j * Nspots] = centroid_error[k + j * Nspots];
            reference_spot_error[k + j * Nspots] = centroid_error[k + j * Nspots];
        }
    }

    /* free up the  memory */
    cpl_matrix_delete(SpotLocations);
    cpl_matrix_delete(xy_centre);
    cpl_matrix_delete(xy_centre_err);
    cpl_matrix_delete(xy_sigma);
    cpl_matrix_delete(xy_sigma_err);
    cpl_matrix_delete(xy_fwhm);
    cpl_matrix_delete(xy_fwhm_err);
    cpl_matrix_delete(centre_intensities);
    cpl_array_delete(all_error_codes);
    return(error);
} /* end of gvacqPupilTrackerFCReferenceSpots  */


/* ----------------------------------------------------------------------*/

/**
 * @internal
 * @brief  Determines the bary centers from the pupil tracker spots
 * @param  PTDetectedCenters pupil trcaker spots [16 spots *2(x,y)]
 * @param  PTDetectedCentersError pupil trcaker spots error [16 spots *2(x,y)]
 * @param  PTBarycenters (output) returns the barycenter positions [4  * 2(x,y)]
 * @param  PTBarycentersError (output) returns the barycenter positions error [4  * 2(x,y)]
 * @see   gvacqPupilTracker
 * @par Error Handling:
 * CPL_ERROR_NULL_INPUT if(one of ) the input pointer (s) is NULL
 * CPL_ERROR_NULL
 */
cpl_error_code gvacqPupilTrackerCentroids2BaryCenters(cpl_matrix *PTDetectedCenters,
                                                      cpl_matrix *PTDetectedCentersError,
                                                      int llx, int lly,
                                                      double *PTBarycenters,
                                                      double *PTBarycentersError) {
/* test entries */
    if ((PTDetectedCenters == NULL) || (PTDetectedCentersError == NULL) ||
        (PTBarycenters == NULL) || (PTBarycentersError == NULL)) {
        printf("The input parameters are NULL. Occured at %s:%s:%d\n",
               __FILE__, __FUNCTION__, __LINE__);
        return(-1);
    }

    int            NlopesX = 4;
    cpl_error_code error   = CPL_ERROR_NONE;
    int            i, j;



    /* initilise arrays  */
    for (i = 0; i < 8; i++) {
        PTBarycenters[i]      = 0.0;
        PTBarycentersError[i] = 0.0;
    }



    for (i = 0; i < NlopesX; i++) {
        for (j = 0; j < NlopesX; j++) {
            PTBarycenters[i]     += cpl_matrix_get(PTDetectedCenters, 0, j + i * NlopesX) + llx - 1;
            PTBarycenters[i + 4] += cpl_matrix_get(PTDetectedCenters, 1, j + i * NlopesX) + lly - 1;

            PTBarycentersError[i]           += cpl_matrix_get(PTDetectedCentersError, 0, j + i * NlopesX);
            PTBarycentersError[i + NlopesX] += cpl_matrix_get(PTDetectedCentersError, 1, j + i * NlopesX);
        } /* j */
    }     /* i */


    for (i = 0; i < 2 * NlopesX; i++) {
        PTBarycenters[i]      = PTBarycenters[i] / ((double)NlopesX);
        PTBarycentersError[i] = PTBarycentersError[i] / ((double)NlopesX);
    } /* i */

    return(error);
}


cpl_error_code gvacqPupilTrackerBaryCenters2PupilPositions(double *Fiber_coupler_reference_centers,
                                                           double *PTBarycenters,
                                                           double *PTBarycentersError,
                                                           int    *MissingSubapertures,
                                                           double *PupilPosition,
                                                           double *PupilPositionError) {
    /* test entries */
    if ((Fiber_coupler_reference_centers == NULL) || (PTBarycenters == NULL) ||
        (PTBarycentersError == NULL) || (PupilPosition == NULL) ||
        (PupilPositionError == NULL)) {
        printf("The input parameters are NULL. Occured at %s:%s:%d\n",
               __FILE__, __FUNCTION__, __LINE__);
        return(-1);
    }


    /* initialise */
    int            i;
    int            NSlopesX              = 4;
    cpl_error_code error                 = CPL_ERROR_NONE;
    double         lat_pupil_shift_x     = 0.0, lat_pupil_shift_y = 0.0;
    double         lat_pupil_shift_x_err = 0.0, lat_pupil_shift_y_err = 0.0;
    double         ShiftsArray[8];
    double         longi_pupil_shift_meter, longi_pupil_shift_meter_err;
    double         longi_pupil_shift_px, longi_pupil_shift_px_err;
    double         Wrong_threshold  = 5; /* Narsi suggest to lower this number */
    double         radius           = 57;
    double         radius_min_thres = radius - 2 * Wrong_threshold;
    double         radius_max_thres = radius + 2 * Wrong_threshold;

    int sub0 = 0, sub1 = 0, sub2 = 0, sub3 = 0;
    if (MissingSubapertures[0] == 0) sub0 = 1;
    if (MissingSubapertures[1] == 0) sub1 = 1;
    if (MissingSubapertures[2] == 0) sub2 = 1;
    if (MissingSubapertures[3] == 0) sub3 = 1;


    if (sub0 + sub1 + sub2 + sub3 < 2) error = -3;

    if (sub0 + sub1 + sub2 + sub3 == 3) {
        lat_pupil_shift_x     = 0.0;
        lat_pupil_shift_y     = 0.0;
        lat_pupil_shift_x     = 0.0;
        lat_pupil_shift_x_err = 0.0;
        lat_pupil_shift_y_err = 0.0;
        longi_pupil_shift_px  = 0.0;

        if (sub0 == 0) {
            PTBarycenters[0] = PTBarycenters[1];
            PTBarycenters[4] = PTBarycenters[7];
            if (PRINTMSG) printf("Missing 4th Barycenter predicted: %f %f  \n",
                                 PTBarycenters[0], PTBarycenters[4]);
        }

        if (sub1 == 0) {
            PTBarycenters[1] = PTBarycenters[0];
            PTBarycenters[5] = PTBarycenters[6];
            if (PRINTMSG) printf("Missing 4th Barycenter predicted: %f %f  \n",
                                 PTBarycenters[1], PTBarycenters[5]);
        }

        if (sub2 == 0) {
            PTBarycenters[2] = PTBarycenters[3];
            PTBarycenters[6] = PTBarycenters[5];
            if (PRINTMSG) printf("Missing 4th Barycenter predicted: %f %f  \n",
                                 PTBarycenters[2], PTBarycenters[6]);
        }

        if (sub3 == 0) {
            PTBarycenters[3] = PTBarycenters[2];
            PTBarycenters[7] = PTBarycenters[4];
            if (PRINTMSG) printf("Missing 4th Barycenter predicted: %f %f  \n",
                                 PTBarycenters[3], PTBarycenters[7]);
        }
    }


    /* Finding out the shifts between the refs and current measured barycenters   */
    for (i = 0; i < 8; i++) {
        ShiftsArray[i] = PTBarycenters[i] - Fiber_coupler_reference_centers[i];
    }

    /* Measure the lateral pupil shifts in pixels */
    for (i = 0; i < NSlopesX; i++) {
        lat_pupil_shift_x += ShiftsArray[i];            /* in pixels */
        lat_pupil_shift_y += ShiftsArray[i + NSlopesX]; /* in pixels */

        lat_pupil_shift_x_err += PTBarycentersError[i];
        lat_pupil_shift_y_err += PTBarycentersError[i + NSlopesX];
    } /* i */

    lat_pupil_shift_x     = lat_pupil_shift_x / 4.0;
    lat_pupil_shift_y     = lat_pupil_shift_y / 4.0;
    lat_pupil_shift_x_err = lat_pupil_shift_x_err / 4.0;
    lat_pupil_shift_y_err = lat_pupil_shift_y_err / 4.0;

    /* Measure the longitudinal pupil displacement */
    longi_pupil_shift_px = -((ShiftsArray[0] - ShiftsArray[3]) - (ShiftsArray[2] - ShiftsArray[1]) +
                             (ShiftsArray[4] + ShiftsArray[7]) - (ShiftsArray[6] + ShiftsArray[5])) / 8.0;
    /* (x4-x1 + x3-x2 + y2-y1 + y3-y4)/8 positive as enlarging */
    /* assuming that x4-x1=x3-x2=y2-y1=y3-y4 */
    /* 2 3 */
    /* 1 4 */

    /* Calculation of the errors */
    longi_pupil_shift_px_err = -((PTBarycentersError[0] - PTBarycentersError[3]) -
                                 (PTBarycentersError[2] - PTBarycentersError[1]) +
                                 (PTBarycentersError[4] + PTBarycentersError[7]) -
                                 (PTBarycentersError[6] + PTBarycentersError[5])) / 8.0;

    longi_pupil_shift_meter     = gvacqPupilTrackerLongPostionPixels2Meters(longi_pupil_shift_px);
    longi_pupil_shift_meter_err = gvacqPupilTrackerLongPostionPixels2Meters(longi_pupil_shift_px_err);



    if (sub0 + sub1 + sub2 + sub3 >= 3) {
        radius = sqrt((PTBarycenters[0] - PTBarycenters[3]) * (PTBarycenters[0] - PTBarycenters[3]) +
                      (PTBarycenters[4] - PTBarycenters[7]) * (PTBarycenters[4] - PTBarycenters[7]));

        if ((fabs(PTBarycenters[0] - PTBarycenters[1]) > Wrong_threshold) ||
            (fabs(PTBarycenters[2] - PTBarycenters[3]) > Wrong_threshold) ||
            (fabs(PTBarycenters[4] - PTBarycenters[7]) > Wrong_threshold) ||
            (fabs(PTBarycenters[5] - PTBarycenters[6]) > Wrong_threshold) ||
            (fabs(fabs(PTBarycenters[0] - PTBarycenters[3]) - fabs(PTBarycenters[4] - PTBarycenters[5])) > (3 + Wrong_threshold)) ||
            (fabs(fabs(PTBarycenters[1] - PTBarycenters[2]) - fabs(PTBarycenters[6] - PTBarycenters[7])) > (3 + Wrong_threshold)) ||
            radius<radius_min_thres || radius> radius_max_thres) {
            if (PRINTMSG) printf("%f %f %f %f %f %f, Error %s:%s:%d\n",
                                 fabs(PTBarycenters[0] - PTBarycenters[1]), fabs(PTBarycenters[2] - PTBarycenters[3]),
                                 fabs(PTBarycenters[4] - PTBarycenters[7]), fabs(PTBarycenters[5] - PTBarycenters[6]),
                                 fabs(fabs(PTBarycenters[0] - PTBarycenters[3]) - fabs(PTBarycenters[4] - PTBarycenters[5])),
                                 radius, __FILE__, __FUNCTION__, __LINE__);

            for (i = 0; i < 8; i++) {
                PTBarycenters[i] = 0;
            }
            error = -3;
        }
    }

    /* Two sub-apertures case */
    /* we assume that the longitudinal shift only moves the baricenters at 45deg relative to x,y */
    /* convention : suba1 suba2 */
    /*              suba0 suba3 */

    if (sub0 + sub1 + sub2 + sub3 == 2) {
        lat_pupil_shift_x     = 0.0;
        lat_pupil_shift_y     = 0.0;
        lat_pupil_shift_x     = 0.0;
        lat_pupil_shift_x_err = 0.0;
        lat_pupil_shift_y_err = 0.0;
        longi_pupil_shift_px  = 0.0;

        for (i = 0; i < 8; i++) {
            ShiftsArray[i] = PTBarycenters[i] - Fiber_coupler_reference_centers[i];
        }

        if ((sub0 == 1) && (sub1 == 1)) {
            if (PRINTMSG) printf("--------> Missing two sub-apertures 2, 3\n");


            longi_pupil_shift_px        = ((ShiftsArray[5] - ShiftsArray[4])) / (2.0);
            longi_pupil_shift_meter     = gvacqPupilTrackerLongPostionPixels2Meters(longi_pupil_shift_px);
            longi_pupil_shift_meter_err = gvacqPupilTrackerLongPostionPixels2Meters(longi_pupil_shift_px_err);

            lat_pupil_shift_x = (ShiftsArray[0] + ShiftsArray[1]) / 2.0 + longi_pupil_shift_px;

            /* adding longitudinal pupil shifts becuase the lateral shifts measured by the two sub-apertures strategy
             * has offsets with pupil defocus */
            lat_pupil_shift_y = (ShiftsArray[4] + ShiftsArray[5]) / 2.0;

            PTBarycenters[2] = 0;
            PTBarycenters[3] = 0;
            PTBarycenters[6] = 0;
            PTBarycenters[7] = 0;

            radius = sqrt((PTBarycenters[0] - PTBarycenters[1]) * (PTBarycenters[0] - PTBarycenters[1]) +
                          (PTBarycenters[4] - PTBarycenters[5]) * (PTBarycenters[4] - PTBarycenters[5]));

            if ((fabs(ShiftsArray[0] - ShiftsArray[1]) > Wrong_threshold) || radius<radius_min_thres ||
                                                                                    radius> radius_max_thres) {
                if (PRINTMSG) printf("TwoSubs case: fabs(ShiftsArray[0]-ShiftsArray[1]): %f, Error %s:%s:%d\n",
                                     fabs(ShiftsArray[0] - ShiftsArray[1]), __FILE__, __FUNCTION__, __LINE__);
                for (i = 0; i < 8; i++) {
                    PTBarycenters[i] = 0;
                }
                error = -3;
            }
        } else if ((sub0 == 1) && (sub3 == 1))           {
            if (PRINTMSG) printf("--------> Misssing  two sub-apertures 1, 2 \n");

            longi_pupil_shift_px        = ((ShiftsArray[3] - ShiftsArray[0])) / (2.0);
            longi_pupil_shift_meter     = gvacqPupilTrackerLongPostionPixels2Meters(longi_pupil_shift_px);
            longi_pupil_shift_meter_err = gvacqPupilTrackerLongPostionPixels2Meters(longi_pupil_shift_px_err);

            lat_pupil_shift_x = (ShiftsArray[0] + ShiftsArray[3]) / 2.0;
            lat_pupil_shift_y = (ShiftsArray[4] + ShiftsArray[7]) / 2.0 + longi_pupil_shift_px;

            PTBarycenters[1] = 0;
            PTBarycenters[2] = 0;
            PTBarycenters[5] = 0;
            PTBarycenters[6] = 0;

            radius = sqrt((PTBarycenters[0] - PTBarycenters[3]) * (PTBarycenters[0] - PTBarycenters[3]) +
                          (PTBarycenters[4] - PTBarycenters[7]) * (PTBarycenters[4] - PTBarycenters[7]));

            if ((fabs(ShiftsArray[4] - ShiftsArray[7]) > Wrong_threshold) || radius<radius_min_thres ||
                                                                                    radius> radius_max_thres) {
                if (PRINTMSG) printf("TwoSubs case: fabs(ShiftsArray[4]-ShiftsArray[7]): %f, Error %s:%s:%d\n",
                                     fabs(ShiftsArray[4] - ShiftsArray[7]), __FILE__, __FUNCTION__, __LINE__);
                for (i = 0; i < 8; i++) {
                    PTBarycenters[i] = 0;
                }
                error = -3;
            }
        } else if ((sub1 == 1) && (sub2 == 1))           {
            if (PRINTMSG) printf("--------> Misssing two sub-apertures  0, 3\n");

            longi_pupil_shift_px        = ((ShiftsArray[2] - ShiftsArray[1])) / (2.0);
            longi_pupil_shift_meter     = gvacqPupilTrackerLongPostionPixels2Meters(longi_pupil_shift_px);
            longi_pupil_shift_meter_err = gvacqPupilTrackerLongPostionPixels2Meters(longi_pupil_shift_px_err);

            lat_pupil_shift_x = (ShiftsArray[1] + ShiftsArray[2]) / 2.0;
            lat_pupil_shift_y = (ShiftsArray[5] + ShiftsArray[6]) / 2.0 - longi_pupil_shift_px;


            PTBarycenters[3] = 0;
            PTBarycenters[0] = 0;
            PTBarycenters[4] = 0;
            PTBarycenters[7] = 0;

            radius = sqrt((PTBarycenters[1] - PTBarycenters[2]) * (PTBarycenters[1] - PTBarycenters[2]) +
                          (PTBarycenters[5] - PTBarycenters[6]) * (PTBarycenters[5] - PTBarycenters[6]));

            if ((fabs(ShiftsArray[5] - ShiftsArray[6]) > Wrong_threshold) || radius<radius_min_thres ||
                                                                                    radius> radius_max_thres) {
                if (PRINTMSG) printf("TwoSubs case: fabs(ShiftsArray[5]-ShiftsArray[6]): %f, Error %s:%s:%d\n",
                                     fabs(ShiftsArray[5] - ShiftsArray[6]), __FILE__, __FUNCTION__, __LINE__);
                for (i = 0; i < 8; i++) {
                    PTBarycenters[i] = 0;
                }
                error = -3;
            }
        } else if ((sub2 == 1) && (sub3 == 1))           {
            if (PRINTMSG) printf("--------> Misssing two sub-apertures 0, 1 \n");

            longi_pupil_shift_px        = ((ShiftsArray[6] - ShiftsArray[7])) / (2.0);
            longi_pupil_shift_meter     = gvacqPupilTrackerLongPostionPixels2Meters(longi_pupil_shift_px);
            longi_pupil_shift_meter_err = gvacqPupilTrackerLongPostionPixels2Meters(longi_pupil_shift_px_err);

            lat_pupil_shift_x = (ShiftsArray[2] + ShiftsArray[3]) / 2.0 - longi_pupil_shift_px;
            lat_pupil_shift_y = (ShiftsArray[6] + ShiftsArray[7]) / 2.0;

            PTBarycenters[0] = 0;
            PTBarycenters[1] = 0;
            PTBarycenters[4] = 0;
            PTBarycenters[5] = 0;

            radius = sqrt((PTBarycenters[2] - PTBarycenters[3]) * (PTBarycenters[2] - PTBarycenters[3]) +
                          (PTBarycenters[6] - PTBarycenters[7]) * (PTBarycenters[6] - PTBarycenters[7]));

            if ((fabs(ShiftsArray[2] - ShiftsArray[3]) > Wrong_threshold) || radius<radius_min_thres ||
                                                                                    radius> radius_max_thres) {
                if (PRINTMSG) printf("TwoSubs case: %f %f, Error %s:%s:%d\n",
                                     fabs(ShiftsArray[2] - ShiftsArray[3]), radius,
                                     __FILE__, __FUNCTION__, __LINE__);
                for (i = 0; i < 8; i++) {
                    PTBarycenters[i] = 0;
                }
                error = -3;
            }
        } else error = -3;  /* End */

        if ((fabs(lat_pupil_shift_x) > TwoSubApertureLatPupilThreshold) ||
            (fabs(lat_pupil_shift_y) > TwoSubApertureLatPupilThreshold)) {
            if (PRINTMSG) printf("TwoSubs case: %f %f, Error %s:%s:%d\n",
                                 lat_pupil_shift_x, lat_pupil_shift_y, __FILE__, __FUNCTION__, __LINE__);
            error = -3;
        }
    }

    /*assiging values for return */
    PupilPosition[0]      = lat_pupil_shift_x;
    PupilPosition[1]      = lat_pupil_shift_y;
    PupilPosition[2]      = longi_pupil_shift_meter;
    PupilPositionError[0] = lat_pupil_shift_x_err;
    PupilPositionError[1] = lat_pupil_shift_y_err;
    PupilPositionError[2] = longi_pupil_shift_meter_err;
    return(error);
}


/**
 * @internal
 * @brief  Converts longitudinal pupil position from pixels to meters
 * @param  PositionPixels in pixel inits
 * @param  Flag_telescope which telescope is operating: for AT, UT 1 and 0
 * @return PositionMeters in meters
 * @see   gvacqPupilTracker
 */
double gvacqPupilTrackerLongPostionPixels2Meters(double PositionPixels) {
    double D_beam    = 18e-3; /* meter */
    double D_pixel   = 18e-6;
    double D_AT      = 1.8;
    double D_lenslet = 2 * subap_pupil;
    double Defocus_coeff_rms_rad, Defocus_coeff_rms_meter, Defocus_coeff_P2V_meter;
    double LongitudinalDefocusShift, PositionMeters;

    Defocus_coeff_rms_rad    = PositionPixels * D_pixel * D_beam / (f_PT * D_lenslet); /* Radians */
    Defocus_coeff_rms_meter  = Defocus_coeff_rms_rad * Llambda / (2 * PI);             /* meters */
    Defocus_coeff_P2V_meter  = 3.5 * Defocus_coeff_rms_meter;                          /* RMS to peak to valley */
    LongitudinalDefocusShift = 8 * (f_PT / D_lenslet) * (f_PT / D_lenslet) * Defocus_coeff_P2V_meter;
    PositionMeters           = f_lens * f_lens * LongitudinalDefocusShift / (f_PT + LongitudinalDefocusShift) /
                               f_PT * (D_AT / D_lenslet); /* Longitudinal pupil shift in meters */

    return(PositionMeters);
}


/* ----------------------------------------------------------------------*/

/**
 * @internal
 * @brief  Determines the lateral and longitudinal position of pupil using centiod postions
 * @param  fiber_coupler_reference (input)   Reference case of centoid data [2*4]
 * @param  xy_centre_distorted (input)  distorted case centoid data [2*16]
 * @param  xy_centre_distorted_error (input) distorted case centoid data [2*16]
 * @param  PupilPosition (output) measures the pupil positions for the centroid data [x,y,z]
 * @param  PupilPosition (output) measures the pupil position errors for the centroid data [sigma_x,sigma_y, sigma_z]
 * @return error
 * @see   gvacqPupilTracker
 *
 * @par Error Handling:
 * CPL_ERROR_NULL_INPUT if(one of ) the input pointer (s) is NULL
 * CPL_ERROR_ILLEGAL_INPUT if the dims of xy_centre, xy_centre_distorted not same
 * CPL_ERROR_NULL
 */
cpl_error_code gvacqPupilShiftsMeasureFromCentroids(double     *fiber_coupler_reference,
                                                    cpl_matrix *xy_centre_distorted,
                                                    cpl_matrix *xy_centre_distorted_error,
                                                    double     *PupilPosition,
                                                    double     *PupilPositionError) {
    /* test entries */
    if ((fiber_coupler_reference == NULL) || (xy_centre_distorted == NULL) ||
        (xy_centre_distorted_error == NULL)) {
        printf("The input parameters are NULL. Occured at %s:%s:%d\n",
               __FILE__, __FUNCTION__, __LINE__);
        return(-1);
    }


    int            NlopesX = 4;
    cpl_error_code error   = CPL_ERROR_NONE;
    int            i, j;
    double         lat_pupil_shift_x     = 0.0, lat_pupil_shift_y = 0.0;
    double         lat_pupil_shift_x_err = 0.0, lat_pupil_shift_y_err = 0.0, longi_pupil_shift_err = 0.0;
    double         BaryCenterArray[2 * NlopesX];
    double         BaryCenterArrayError[2 * NlopesX];
    double         ShiftsArray[8];

    /* cpl_matrix_dump(xy_centre_distorted, NULL);    */

    /* initilise arrays  */
    for (i = 0; i < 2 * NlopesX; i++) {
        BaryCenterArray[i]      = 0.0;
        BaryCenterArrayError[i] = 0.0;
    }

    for (i = 0; i < NlopesX; i++) {
        for (j = 0; j < NlopesX; j++) {
            BaryCenterArray[i]     += cpl_matrix_get(xy_centre_distorted, 0, j + i * NlopesX);
            BaryCenterArray[i + 4] += cpl_matrix_get(xy_centre_distorted, 1, j + i * NlopesX);

            BaryCenterArrayError[i]           += cpl_matrix_get(xy_centre_distorted_error, 0, j + i * NlopesX);
            BaryCenterArrayError[i + NlopesX] += cpl_matrix_get(xy_centre_distorted_error, 1, j + i * NlopesX);
        } /* j */
    }     /* i */

    for (i = 0; i < 2 * NlopesX; i++) {
        BaryCenterArray[i]      = BaryCenterArray[i] / ((double)NlopesX);
        BaryCenterArrayError[i] = BaryCenterArrayError[i] / ((double)NlopesX);
    }

    for (i = 0; i < 2 * NlopesX; i++) {
        ShiftsArray[i] = BaryCenterArray[i] - fiber_coupler_reference[i];
    }

    for (i = 0; i < NlopesX; i++) {
        /* printf("%f %f \n", BaryCenterArrayError[i], BaryCenterArrayError[i+NlopesX]); */
        lat_pupil_shift_x += ShiftsArray[i];           /* in pixels */
        lat_pupil_shift_y += ShiftsArray[i + NlopesX]; /* in pixels */

        lat_pupil_shift_x_err += BaryCenterArrayError[i];
        lat_pupil_shift_y_err += BaryCenterArrayError[i + NlopesX];
    }

    lat_pupil_shift_x     = lat_pupil_shift_x / 4.0;
    lat_pupil_shift_y     = lat_pupil_shift_y / 4.0;
    lat_pupil_shift_x_err = lat_pupil_shift_x_err / 4.0;
    lat_pupil_shift_y_err = lat_pupil_shift_y_err / 4.0;

    /* printf("%f %f \n", lat_pupil_shift_x_err, lat_pupil_shift_y_err); */


    /* The longitudinal pupil displacement */
    double longi_pupil_shift_px = -((ShiftsArray[0] - ShiftsArray[3]) - (ShiftsArray[2] - ShiftsArray[1]) +
                                    (ShiftsArray[4] + ShiftsArray[7]) - (ShiftsArray[6] + ShiftsArray[5]));

    /* Calculate the errors */
    longi_pupil_shift_err = -((BaryCenterArrayError[0] - BaryCenterArrayError[3]) -
                              (BaryCenterArrayError[2] - BaryCenterArrayError[1]) +
                              (BaryCenterArrayError[4] + BaryCenterArrayError[7]) -
                              (BaryCenterArrayError[6] + BaryCenterArrayError[5]));

    /*assiging values for return */
    PupilPosition[0] = lat_pupil_shift_x;
    PupilPosition[1] = lat_pupil_shift_y;
    PupilPosition[2] = longi_pupil_shift_px;

    PupilPositionError[0] = lat_pupil_shift_x_err;
    PupilPositionError[1] = lat_pupil_shift_y_err;
    PupilPositionError[2] = longi_pupil_shift_err;

    return(error);
} /* end of gvacqPupilShiftsMeasureFromCentroids */


/*
 * It appcepts cpl_aperture array and return the positions of apertures in  double array
 */
cpl_error_code cpl_apertures2XYdoubleArray(cpl_apertures *apert,
                                           int nx, double *DataDouble) {
    if ((apert == NULL) || (nx <= 0) || (DataDouble == NULL))
        printf("Inputs are not allocated. Occured at %s:%s:%d \n",
               __FILE__, __FUNCTION__, __LINE__);

    int i;
    for (i = 0; i < nx; i++) {
#ifdef VLT2014
        DataDouble[i]      = cpl_apertures_get_centroid_x(apert, i + 1);
        DataDouble[i + nx] = cpl_apertures_get_centroid_y(apert, i + 1);
#else
        DataDouble[i]      = cpl_apertures_get_centroid_x(apert, i + 1);
        DataDouble[i + nx] = cpl_apertures_get_centroid_y(apert, i + 1);
#endif
    }
    return(0);
}


/* Setting database initial points to the local global variabls */

/*
 * It is the pupil tracker box size center positions for 4 telescopes
 */
void setRefWindowPosPT(vltINT32 *RefWindowPosPT_in) {
    int i;

    if (PRINTMSG) {
        cpl_msg_info(cpl_func, "Set PT window Pos X (T1-T4):  %d %d %d %d\n",
                     RefWindowPosPT_in[0], RefWindowPosPT_in[1],
                     RefWindowPosPT_in[2], RefWindowPosPT_in[3]);
        cpl_msg_info(cpl_func, "Set PT window Pos Y (T1-T4):  %d %d %d  %d  \n",
                     RefWindowPosPT_in[4], RefWindowPosPT_in[5],
                     RefWindowPosPT_in[6], RefWindowPosPT_in[7]);
    }
    for (i = 0; i < 8; i++) {
        RefWindowPosPT[i] = RefWindowPosPT_in[i];
    }
}


/*
 * It is the pupil tracker box size where spots are searched
 */
void setPTImageWindowSize(vltINT32 PTImageWindowSize_in) {
    PTImageWindowSize = PTImageWindowSize_in;
    if (PRINTMSG) cpl_msg_info(cpl_func, "PTImageWindowSize= %d \n",
                               PTImageWindowSize);
}


/*
 * Once the pupil tracker spots known roughly, this window size is used
 * to fit Gaussian at the center of rough positions
 */
void setPTwindowGauss(vltINT32 PTwindowGauss_in) {
    PTwindowGauss = PTwindowGauss_in;
    if (PRINTMSG) cpl_msg_info(cpl_func, "PTwindowGauss= %d \n",
                               PTwindowGauss);
}


/*
 * To find pupil tracker spots, an array of sigmas used as threshold.
 * If first element fails second one used and so on ..
 */
void setPTspotsScanSigma(double *PTspotsScanSigma_in) {
    int i;

    for (i = 0; i < 4; i++) {
        PT_spots_scan_sigmaOrig[i] = PTspotsScanSigma_in[i];
    }
    if (PRINTMSG) cpl_msg_info(cpl_func, "PT_spots_scan_sigma= %f %f %f %f\n",
                               PTspotsScanSigma_in[0], PTspotsScanSigma_in[1],
                               PTspotsScanSigma_in[2], PTspotsScanSigma_in[3]);
}


/**
 * @brief  Quickly scans if there exits any objects in the input image
 * @param  fwhm_threshold Thershold fwhm is used to reject hot pixels and cosmic events
 * @return error
 *
 * Used for field imager and pupil tracker
 */
cpl_error_code Check_objects_found_PT(cpl_image *FI_image,
                                      double    fwhm_threshold_min,
                                      double    fwhm_threshold_max) {
    cpl_error_code error = CPL_ERROR_NONE;
    int            j, k;

#ifdef VLT2014
    cpl_size px1, py1;
#else
    int px1, py1;
#endif
    int Check_iteration = 50;

    int    im_size_x, im_size_y, status;
    double fwhm_xx, fwhm_yy;

    /*test entries */
    if (FI_image == NULL) return(-1);

    im_size_x = cpl_image_get_size_x(FI_image);
    im_size_y = cpl_image_get_size_y(FI_image);
    status    = 0;


    do {
        error = cpl_image_get_maxpos(FI_image, &px1, &py1); /* Get X and Y pixel positions*/

        if (error != CPL_ERROR_NONE)
            cpl_msg_error(cpl_func, "Detected at %s:%s:%d\n", __FILE__,
                          __FUNCTION__, __LINE__);


        /*To avoid the objects/pixels which are edge of the window */
        if ((im_size_x - px1 < 3) || (px1 < 3) || (im_size_y - py1 < 3) || (py1 < 3)) {
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

        if ((fwhm_xx <= FIfwhm_minthreshold) || (fwhm_yy <= FIfwhm_minthreshold))
            for (j = 1; j < 2; j++) {
                for (k = 1; k < 2; k++) {
                    error = cpl_image_set(FI_image, px1 - 1 + j, py1 - 1 + k, 0.0);
                    if (error != CPL_ERROR_NONE)
                        cpl_msg_error(cpl_func, "Detected at %s:%s:%d\n", __FILE__,
                                      __FUNCTION__, __LINE__);
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
        if (status >= Check_iteration)

            return(-2);

    } while (fwhm_xx <= fwhm_threshold_min ||
             fwhm_yy <= fwhm_threshold_min ||
             fwhm_xx >= fwhm_threshold_max ||
             fwhm_yy >= fwhm_threshold_max); /*end of do while  */


    return(error);
}


/*
 * Computes flux in the pupil tracker windows
 * using cpl_image_get_flux_window.
 */
cpl_error_code PupilWindowFlux(cpl_image *PupilDetector,
                               int       *PTWindow_coor,
                               int       WindowSize,
                               double    *Flux) {
    cpl_error_code error = CPL_ERROR_NONE;
    int            llx, lly, urx, ury, i;

    for (i = 0; i < 4; i++) { /* Do it for all four telescopes  */
        llx = PTWindow_coor[i] - WindowSize / 2;
        lly = PTWindow_coor[i + 4] - WindowSize / 2;
        urx = PTWindow_coor[i] + WindowSize / 2;
        ury = PTWindow_coor[i + 4] + WindowSize / 2;

        Flux[i] = cpl_image_get_flux_window(PupilDetector, llx, lly, urx, ury);
    }

    return(error);
}


/*
 * This script goes to the spot position and sums up the pixel intenities with in a window of 10 pixels.
 *
 * If the image is not preprocessed, the measurement also contains background counts
 */
cpl_error_code PTSpotsFlux(cpl_image *fimage,
                           int *imistart2, int *imjstart2, int nspots,
                           double *Flux) {
    if ((fimage == NULL) || (imistart2 == NULL) || (imjstart2 == NULL)) {
        printf("ERROR::-> Inputs for the function are NULL at %s:%s:%d\n",
               __FILE__, __FUNCTION__, __LINE__);

        return(-1);
    }


    int            l, i;
    cpl_error_code error = CPL_ERROR_NONE;
    i = -1;
    for (l = 0; l < nspots; l++) {
        i++;

        if ((imistart2[l] > 10) && (imistart2[l] > 10) && (imistart2[l] < 200) && (imjstart2[l] < 200))
            Flux[i] = cpl_image_get_absflux_window(fimage, imistart2[l] - 5, imjstart2[l] - 5,
                                                   imistart2[l] + 5, imjstart2[l] + 5);

    }

    return(error);
}


/* tool to provide the minimal distance from a set of apertures */
double minApertureDistance(cpl_apertures *PTspots) {
    int    i, j, Nspots;
    int    px[20], py[20];
    double d, dtemp;

    Nspots = cpl_apertures_get_size(PTspots);
    if (Nspots < 2) return(0);
    if (Nspots >20 ) Nspots=20;

    for (j = 0; j < Nspots; j++) {
        px[j] = cpl_apertures_get_centroid_x(PTspots, j + 1);
        py[j] = cpl_apertures_get_centroid_y(PTspots, j + 1);
    }
    d = sqrt((px[0] - px[1]) * (px[0] - px[1]) + (py[0] - py[1]) * (py[0] - py[1]));
    for (i = 0; i < Nspots; i++) {
        for (j = i + 1; j < Nspots; j++) {
            dtemp = sqrt((px[i] - px[j]) * (px[i] - px[j]) + (py[i] - py[j]) * (py[i] - py[j]));
            if (dtemp < d) d = dtemp;
        }
    }

    return(d);
}


/* tool to cut on spot multiplicity 0 */
int spotsInApertureMultiplicity(cpl_apertures *PTspots) {
    int    i, j, Nspots;
    int    px[20], py[20], npar[20];
    double d, dtemp;

    Nspots = cpl_apertures_get_size(PTspots);
    if (Nspots < 2) return(0);
    if (Nspots >20 ) Nspots=20;

    for (j = 0; j < Nspots; j++) {
        px[j] = cpl_apertures_get_centroid_x(PTspots, j + 1);
        py[j] = cpl_apertures_get_centroid_y(PTspots, j + 1);
        npar[j]=0;
    }
    for (i = 0; i < Nspots; i++) {
        for (j = i + 1; j < Nspots; j++) {
            if ((abs(px[i]-px[j])>52) && (abs(px[i]-px[j])<61) && (abs(py[i]-py[j])<4))  {++npar[i]; ++npar[j];}
            if ((abs(py[i]-py[j])>52) && (abs(py[i]-py[j])<61) && (abs(px[i]-px[j])<4))  {++npar[i]; ++npar[j];}
            if ((abs(px[i]-px[j])>52) && (abs(px[i]-px[j])<61) && (abs(py[i]-py[j])>52) && (abs(py[i]-py[j])<61)){++npar[i]; ++npar[j];}
        }
    }
    if ( npar[Nspots-1]==0) --Nspots;
    if ( npar[Nspots-1]==0) --Nspots;

    return(Nspots);
}


/*___oOo___*/
