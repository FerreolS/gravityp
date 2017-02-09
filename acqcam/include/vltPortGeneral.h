#ifndef VLT_GENERAL_H
#define VLT_GENERAL_H
/*******************************************************************************
* E.S.O. - VLT project
*
* "@(#) $Id: vltPortGeneral.h 268615 2015-04-30 15:31:15Z sfeyrin $"
*
* who       when      what
* --------  --------  ----------------------------------------------
* bgilli   2004/09/21 Created, extracting general interest defines from ccsLite.h
* sfeyrin  2009-01-11 Change vltINT8 to signed char
*/


/************************************************************************
 *                          VLT   Data  Types                           *
 ************************************************************************/

/* for VXWORKS signed to ensure portability for PPC */
typedef signed char        vltINT8;      /*  8 bits integers           */

typedef unsigned char      vltUINT8;     /*  8 bits unsigned integers  */
typedef short              vltINT16;     /* 16 bits integers           */
typedef unsigned short     vltUINT16;    /* 16 bits unsigned integers  */
typedef int                vltINT32;     /* 32 bits integers           */
typedef unsigned int       vltUINT32;    /* 32 bits unsigned integers  */

typedef unsigned char      vltLOGICAL;   /* logical (rtLogical)        */

typedef double             vltDOUBLE;    
typedef float              vltFLOAT;

/*  CCS time structure */

typedef struct {
    unsigned long  tv_sec;   /* seconds since midnight January 1,1970 UTC */
    long           tv_usec;  /* and microseconds		          */
} vltTIMEVAL;


typedef struct  {                     /* Polar data type       */
        vltDOUBLE   magnitude;
        vltDOUBLE   phase;
} vltPOLAR;  


typedef struct  {                     /* Rectangular data type */
        vltDOUBLE   real;
        vltDOUBLE   imaginary;
} vltRECTANGULAR;  

typedef unsigned char      vltBYTES4[4]; 
typedef unsigned char      vltBYTES8[8]; 
typedef unsigned char      vltBYTES12[12]; 
typedef unsigned char      vltBYTES16[16]; 
typedef unsigned char      vltBYTES20[20]; 
typedef unsigned char      vltBYTES32[32]; 
typedef unsigned char      vltBYTES48[48]; 
typedef unsigned char      vltBYTES64[64]; 
typedef unsigned char      vltBYTES80[80]; 
typedef unsigned char      vltBYTES128[128]; 
typedef unsigned char      vltBYTES256[256]; 

/* Data types defined to support string processing */

typedef char  vltSTRING4[4]; 
typedef char  vltSTRING8[8]; 
typedef char  vltSTRING12[12]; 
typedef char  vltSTRING16[16]; 
typedef char  vltSTRING20[20]; 
typedef char  vltSTRING32[32]; 
typedef char  vltSTRING48[48]; 
typedef char  vltSTRING64[64]; 
typedef char  vltSTRING80[80]; 
typedef char  vltSTRING128[128]; 
typedef char  vltSTRING256[256]; 
typedef char  vltSTRING512[512]; 


typedef vltTIMEVAL    vltDATE;
typedef vltTIMEVAL    vltTIME_OF_DAY;
typedef vltTIMEVAL    vltABS_TIME;

#if defined(i386) || defined(__i386)
/*
 * byte swapping macros needed for LINUX on PCs
 */

/* most significant byte of 2-byte integer */
#define CCS_MSB_INT16(x)        (((x) >>  8) & 0xff)
/* least significant byte of 2-byte integer */
#define CCS_LSB_INT16(x)         ((x)        & 0xff)
/* most significant word of 2-word integer */
#define CCS_MSW_INT32(x) (((x) >> 16) & 0xffff)
/* least signifcant byte of 2-word integer */
#define CCS_LSW_INT32(x)  ((x)        & 0xffff)

/* swap the MSB with the LSB of a 16 bit integer */
#define CCS_SWAP_INT16(x) (CCS_MSB_INT16(x) | (CCS_LSB_INT16(x) << 8))

/* swap all bytes of a 32 bit integer */
#define CCS_SWAP_INT32(x) (CCS_SWAP_INT16(CCS_LSW_INT32(x)) << 16) | \
                             CCS_SWAP_INT16(CCS_MSW_INT32(x))
/*void ccs_swap_2int32(int32_t* i1, int32_t* i2)
  { int32_t i=*i1; *i1=*i2; *i2=i; };*/
                          
/* swap all bytes of two consecutive 32 bit integers */
#define CCS_SWAP_2INT32(x_ptr,y_ptr) {int32_t i=*x_ptr;\
                                     *x_ptr=*y_ptr;\
                                     *y_ptr=i;}\
                           *x_ptr = CCS_SWAP_INT32(*x_ptr);         \
                           *y_ptr = CCS_SWAP_INT32(*y_ptr) 
/* courtesy alongino */			   
#define CCS_SWAP_FLOAT(v)  {vltINT32 *i1; i1=(vltINT32 *)&(v); \
                            *i1 = CCS_SWAP_INT32(*i1);} 
#define CCS_SWAP_DOUBLE(v) {vltINT32 *i1,*i2; i1=(vltINT32 *)&(v); \
                            i2=i1+1;CCS_SWAP_2INT32(i1,i2);} 
#endif


#endif /*!VLT_GENERAL_H*/
