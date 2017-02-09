/************************************************************************
 *   This file is part of the E.S.O. - VLTI GRAVITY project
 *
 *  "@(#) $Id: gvoProcessImageAC.c 249945 2013-12-18 07:08:07Z pgarcia $
 *   
 *   who       when        what
 *   --------  ----------  ----------------------------------------------
 *   narsi    2014-03-03  added bicubic spline libs
 *   narsi    2014-03-02  created
 */

/****************************************************************************
 *   NAME 
 *   gvacqBicubicSpline.c - procederes to compute bicubic spline interpolations
 *  
 *   SYNOPSIS
 *   #include <gvoProcessImageAC.h>
 * 
 *   DESCRIPTION
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
 #include <gvoProcessImageAC.h>
 * @endcode
 */



/*-----------------------------------------------------------------------------
                               Includes
 ----------------------------------------------------------------------------*/

#include <gvoProcessImageAC.h> /* */


/*-----------------------------------------------------------------------------
                                Defines 
 ----------------------------------------------------------------------------*/





/*
  library bicubic spline interpolation on image 
  
  
  The bicubic spline is calculated following way
  Firstly, measure N one-dimensional splines across the rows of 
  the 2d image.
  Secondly, generate N spline interpolations for an array values  
  y(x1i, x2), ranging from i=0, ...,M-1.
  Thirdly, calculate a one-dimensional spline through those array values. 
  Finally, spline interpolate to the desired value y(x1, x2)

  References: Numerical Recipes 3rd Edition: The Art of Scientific Computing, 
  William H. Press et al., (2007)
  
  The code is adopted from spline.c of https://github.com/frigaut/yorick-imutil 
  authored by Francois Rigaut.
*/

#include <gvoProcessImageAC.h>
/**
 *@internal
 *@brief Interpolates the data point "y" for  "x" using the input data [xa,ya]  
 *@param xa X array data
 *@param ya Y array data
 *@param n  number of data points
 *@param x  For the data point to construct interpolation
 *@param y (output) interpolated value for x
 *@return error
 */
cpl_error_code evalspline(double *xa, 
			  double *ya, 
			  double *y2a, 
			  double n, 
			  double x, 
			  double *y)
{
    long klo,khi,k;
    double h,b,a;
    cpl_error_code error=CPL_ERROR_NONE;
	
    klo=0;
    khi=n-1;
    while (khi-klo > 1){
    k=(khi+klo) >> 1;
    if (xa[k] > x) khi=k;
    else klo=k;
    }
    h=xa[khi]-xa[klo];
    if (h == 0.0) printf("Bad xa input to routine evalspline");
    a=(xa[khi]-x)/h;
    b=(x-xa[klo])/h;
    *y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0; 
    
    return error;
}

/*
  Measures the derivatives for the given array
*/
cpl_error_code splinf(double *x, 
		      double *y, 
		      long n, 
		      double *y2)
{
    long i,k;
    double p,qn,sig,un;
    double u[n-1];
    cpl_error_code error=CPL_ERROR_NONE;
    y2[0]=u[0]=qn=un=0.0;
    
    for (i=1;i<=n-2;i++) {
    sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
    p=sig*y2[i-1]+2.0;
    y2[i]=(sig-1.0)/p;
    u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
    u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
    }
    
    y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);   
    for (k=n-2;k>=0;k--) y2[k]=y2[k]*y2[k+1]+u[k];
    return error;
}

cpl_error_code
splin2(double *xin, 
       double *yin, 
       double *image, 
       double *deriv, 
       long nx, 
       long ny, 
       long *nvalidx, 
       double xout, 
       double yout, 
       double *res)
{
    long j,m;
    long n=0;
    double y2tmp[ny]; 
    double yytmp[ny];
    cpl_error_code error=CPL_ERROR_NONE;
    for (j=0;j<=ny-1;j++) {
    m = nvalidx[j];
    error=evalspline(&xin[n],&image[n],&deriv[n],m,xout,&yytmp[j]);
    n += m;
    }
    
    error=splinf(yin,yytmp,ny,y2tmp);
    error=evalspline(yin,yytmp,y2tmp,ny,yout,res);
    return error;
}

/*
  Measures the derivatives for the given image
*/
cpl_error_code
splie2( double *xin, 
        double *im, 
        long nx, 
        long ny,                                
        double *deriv, 
        long *nvalidx)
{
    long j,m;
    long n=0;
    cpl_error_code error=CPL_ERROR_NONE;
    for (j=0;j<=ny-1;j++) {
    m = nvalidx[j];
    error=splinf(&xin[n],&im[n],m,&deriv[n]);
    n += m;
    }
    return error;
}

/*
  Interpolate regularly sampled image using bicubic spline.
*/
cpl_error_code
spline2(double *xin, 
        double *yin, 
        double *im, 
        double *deriv, 
        long nx, 
        long ny,                                
        double *xout, 
        double *yout, 
        long npt, 
        long *nvalidx, 
        double *res)
{
  long i;
  cpl_error_code error=CPL_ERROR_NONE;
  for (i=0;i<=npt;i++) error=splin2(xin,yin,im,deriv,nx,ny,nvalidx,   
                              xout[i],yout[i],&res[i]);
  return error;
}

/*
  Interpolates an image to different resize
*/
cpl_error_code 
spline2grid(double*xin, 
            double *yin, 
            double *im, 
            double *deriv, 
            long nx, 
            long ny, 
            double *xout, 
            double *yout, 
            long nxout, 
            long nyout, 
            long *nvalidx,
            double *res)
{
    long j,ii,jj;
    double y2tmp[ny]; 
    double yytmp[ny];
    cpl_error_code error=CPL_ERROR_NONE;
    long n;
    long m;
    
    for (ii=0;ii<=nxout-1;ii++) {/* loop on out x */
    n=0;
    
    /* fill Y vector for xout(ii) */
    for (j=0;j<=ny-1;j++) {
    m = nvalidx[j];
    error=evalspline(&xin[n],&im[n],&deriv[n],m,xout[ii],&yytmp[j]);
    n += m;
    }
    /* find second derivative */
    splinf(yin,yytmp,ny,y2tmp);
    
    /* find and fill interpolated out Y vector for this xout(ii) */
    for (jj=0;jj<=nyout-1;jj++) {
    error=evalspline(yin,yytmp,y2tmp,ny,yout[jj],&res[jj*nxout+ii]);
    }
    }
    return error;
}

/*
  Outputs an array of N doubles equally spaced from START to STOP.
*/
cpl_error_code span(double start, double stop, int N, double *result)
{
    cpl_error_code error=CPL_ERROR_NONE;
    int i;
    if(start < 0 || stop < 0 || N < 0) printf("start:stop:N must be positive\n");  
    double diff = (stop-start)/(N-1);
    for(i=0; i< N; i++) result[i]=start+ i*diff;
    return error;
}

cpl_error_code span2Dx (double start, double stop, int N, double * result)
{
    cpl_error_code error =CPL_ERROR_NONE;
    int i,j;
    double xin_tmp[N];
    if(start < 0 || stop < 0 || N < 0) 
        printf("start:stop:N must be positive\n"); 
   
    error= span(start, stop, N, xin_tmp); 
    for(i=0; i<N; i++){
    for(j=0; j<=N; j++){
    result[i*N+j]=xin_tmp[j];
    }
    }
    return error;
}

cpl_error_code span2Dy (double start, double stop, int N, double *result)
{
    cpl_error_code error =CPL_ERROR_NONE;
    int i,j;
    double xin_tmp[N];
    if(start < 0 || stop < 0 || N< 0) 
	printf("start:stop:N must be positive\n"); 
    error= span(start, stop, N, xin_tmp);
    
    for(i=0; i<N; i++){
    for(j=0; j<=N; j++){
    result[i*N+j]=xin_tmp[i];
    }
    }
    return error;
}

/*
  Outputs of an array of ints running from "start" to "stop"
*/
cpl_error_code indgen(int start, int stop, double * result)
{
    cpl_error_code error=CPL_ERROR_NONE;
    int i;
    if(start < 0 || stop < 0 ) printf("start:stop must be positive\n"); 
    for(i=0; i<=stop-start; i++) result[i]=start+i;
    return error;
}

/*
  A sub image will be selected from the image with double sampling using
  bicubic splines
*/
cpl_error_code 
splineResample(cpl_image *image, long dimx, double is, double ie, 
	       double js, double je, int sampling, double *res)
{
    cpl_error_code error=CPL_ERROR_NONE;
    int nx = cpl_image_get_size_x(image);
    int ny = cpl_image_get_size_y(image);
    int i, j;
    int nxreb = nx*ny;
    double xin[nxreb];
    double yin[nx];
    double xreb[nxreb];
    double yreb[nxreb];     
    double imageCopy[nxreb];
    double deriv[nxreb];
    double *im= cpl_image_get_data_double(image); 
    
    if (nx != ny) {
    printf("Accepts only square image and square window interpolation %d %d %f %f\n", 
	   nx, ny, is, ie);
    return -1;
    }
    
    /*initialise the input arrays for interpolation */
    for(i=0; i<nxreb; i++) imageCopy[i]=im[i];
    for(i=0; i<ny; i++) yin[i]=i+1;
    error= span2Dx (1, nx, nx, xin);
    
    /* number of valid points in the interpolation */
    int imask[nxreb];
    for(i=0; i<nxreb; i++)  imask[i]= 1;
    
    long nvalidx[nx];                          
    int sum;
    for(i=0; i<nx; i++){
    sum=0;
    for(j=0; j<nx; j++){
    sum += imask[j];
    }
    nvalidx[i]=sum;  
    }
    splie2(xin,imageCopy,nx,ny,deriv,nvalidx);
    
    /* determine the 2D image grid x and y pixel positions to interpolate   */
    error= span2Dx (is, ie, sampling, xreb);
    error= span2Dy (js, je, sampling, yreb);
    
    int npt =sampling*sampling-1;
    
    /* do the interpolation and fill the res array */
    if(dimx==0){     
    for(i=0; i<sampling*sampling; i++)  res[i]= 0.0;
    spline2(xin, yin, imageCopy,deriv,nx,ny,xreb,yreb,npt,nvalidx,res);    
    }else{
    for(i=0; i<dimx*dimx; i++)  res[i]= 0.0;
    double ingen[dimx];
    double  xout[dimx];
    error=indgen(1, dimx, ingen );
    for(i=0; i<dimx; i++) xout[i]= (ingen[i]-1.0)/dimx*nx+1.0;
    error=spline2grid(xin,yin,im,deriv,nx,ny,xout,xout,dimx,dimx,nvalidx,res);
    }
    return error;
}




