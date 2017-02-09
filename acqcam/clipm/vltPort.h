/*************************************************************************
* E.S.O. - VLT project
*
* "@(#) $Id: vltPort.h.Linux,v 2.251 2009/03/18 09:55:04 vltsccm Exp $" 
*
* vltPort.h for Linux
*
* who        when      what
* --------  ---------- ----------------------------------------------
* eallaert  1998-06-08 created from vltPort.h.Solaris2
* ahuxley   1999-09-29 added definition of _SVID_SOURCE
* bgilli    1999-10-01 added timeval stuff just like on Solaris.
* bgilli    1999-10-06 #define _SVID_SOURCE transformed in #define _GNU_SOURCE
* bgilli    1999-10-06 added missing endif!
* psivera   2002-11-29 _GNU_SOURCE used only for C code, not C++
* bgilli    2004-09-22 added correct vltPortGeneral.h
* pkratzer  2009-06-12 ATTRIBUTE_UNUSED definition added
* pkratzer  2008-06-24 ATTRIBUTE_UNUSED moved to (now included) vltAttrib.h
* pkratzer  2009-02-27 if MAKE_VXWORKS is set include "vxWorks.h"
* sfeyrin   2009-03-18 move include "vxWorks.h" after "vltPortGeneral.h"
*
*/

/************************************************************************
*  vltPort.h - Include file to mask differences between platforms.
*              This file should be included in all source files.
*              It relies on macro definitions preceeding the
*              inclusion of this file.
*                
*  REMARK: This file belongs to the "vltMake" module.
*------------------------------------------------------------------------
*/

#ifndef VLTPORT_H
#define VLTPORT_H

#include <vltAttrib.h>

/*
* When it is used, vltPort.h MUST be the very first file included
* in ANSI ".c" files.
* Cause a syntax error if we detect that any other include file has been
* included before vltPort.h in an ANSI ".c" file.
*/
#if defined(__STDC__) && \
                         (defined(_H_STANDARDS) || \
                          defined(_SYS_STDSYMS_INCLUDED) || \
                          defined(_STANDARDS_H_))
#    error "vltPort.h MUST BE THE VERY FIRST FILE INCLUDED IN ANSI '.c' FILES"
#endif

/*
 * This file is used also by some VxWorks code.
 * To be compatible with existing code, SUN_COMP is defined for both
 * gcc and cc68k, but the following definitiond do not influence cc68k
 */
 
#define LINUX

/* 
 * at present, SELECT is defined in the code using it. It should be done here
 * for all. May be in the next release.
 */

/*
 * Adjust name-space information.
 */
#if defined(_ALL_SOURCE)
#    undef _POSIX_C_SOURCE
#endif

#if defined(_XOPEN_SOURCE)
#    undef _POSIX_C_SOURCE
#endif
  
/* for Thomas Ebert re. ipc.h AH 29/09/99. BGI:Modified in _GNU_SOURCE 06/10/99 */
#ifndef __cplusplus
#define _GNU_SOURCE
#endif

/* For CCS_Lite, just like on Solaris. */
#ifndef MAKE_VXWORKS
#include <stddef.h>
#include <sys/types.h>
#include <sys/time.h>
 
#if !defined(timercmp)
    struct timeval
        {
        long      tv_sec;         /* seconds */
        long      tv_usec;        /* and microseconds */
        };
#define crTIMEVAL_TIMEZONE_DEFINED
#endif
#endif /* MAKE_VXWORKS */
/* 
 * general interest defines
 */
#include "vltPortGeneral.h"

#if defined(__VXWORKS__) || defined(__vxworks)
/* -fno-builtin suppresses -Wformat checking of printf/scanf like functions which should be default with -Wall.
 * Enable these checks explicitely for standard functions
 */
#include <stdio.h>
#ifndef __RTP__
#include <private/stdioP.h> /* In Kernel vscanf, vsscanf, #define vfscanf __svfscanf not in visible API */
#endif
#include <time.h>
#define __printfport(fmtarg, firstvararg) \
        __attribute__((__format__ (__printf__, fmtarg, firstvararg)))
#define __scanfport(fmtarg, firstvararg) \
        __attribute__((__format__ (__scanf__, fmtarg, firstvararg)))
#define __strftimeport(fmtarg, firstvararg) \
        __attribute__((__format__ (__strftime__, fmtarg, firstvararg)))

#ifdef __cplusplus
extern "C" {
#endif

extern int  fprintf (FILE *, const char *, ...) __printfport(2,3);
extern int  printf (const char *, ...) __printfport(1,2);
extern int  sprintf (char *, const char *, ...) __printfport(2,3);
extern int  snprintf (char *, size_t, const char *, ...) __printfport(3,4);
extern int  fdprintf (int fd, const char *fmt, ...) __printfport(2,3);

extern int  vfprintf (FILE *, const char *, va_list) __printfport(2,0);
extern int  vprintf (const char *, va_list) __printfport(1,0);
extern int  vsprintf (char *, const char *, va_list) __printfport(2,0);
extern int  vfdprintf (int fd, const char *fmt, va_list ap) __printfport(2,0);

extern int  fscanf (FILE *, const char *, ...) __scanfport(2,3);
extern int  scanf (const char *, ...) __scanfport(1,2);
extern int  sscanf (const char *, const char *, ...) __scanfport(2,3);

extern int  vscanf (const char *, va_list) __scanfport(1,0);
extern int  vsscanf (const char *, const char *, va_list) __scanfport(2,0);
extern int  vfscanf (FILE *, const char *, va_list) __scanfport(2,0);

extern size_t strftime (char *_s, size_t _n, const char *_format, const struct tm *_tptr) __strftimeport(3,0);

#ifdef __cplusplus
}
#endif
#endif

#ifdef MAKE_VXWORKS
#include "vxWorks.h"
/* removed from sysVltLib.h and sysLib.h to not forcing dependency of apps to BSP internals because of incomplete std API */
#if !defined(_WRS_VXWORKS_MAJOR) || _WRS_VXWORKS_MAJOR <= 5
    extern int  snprintf ( char *str, size_t maxLen, const char *fmt, ...);
#endif
#endif

#endif /*!VLTPORT_H*/
