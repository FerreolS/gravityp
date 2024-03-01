# Process this file with autoconf to produce a configure script.

AC_INIT([GRAVITY Instrument Pipeline], [1.7.0b1], [https://support.eso.org],
[gravity])
AC_PREREQ([2.59])

AC_CONFIG_SRCDIR([Makefile.am])
AC_CONFIG_AUX_DIR([admin])
AC_CONFIG_MACRO_DIR([m4macros])

AC_CANONICAL_HOST
AC_CANONICAL_TARGET

AM_INIT_AUTOMAKE
AM_CONFIG_HEADER([config.h])
AM_SILENT_RULES([yes])

AM_MAINTAINER_MODE

# Immediately before every release do:
#-------------------------------------
#   if (the interface is totally unchanged from previous release)
#       REVISION++;
#   else {
#       /* interfaces have been added, removed or changed */
#       REVISION = 0;
#       CURRENT++;
#       if (any interfaces have been _added_ since last release)
#           AGE++;
#       if (any interfaces have been _removed_ or incompatibly changed)
#           AGE = 0;
#   }

# Order of arguments: VERSION, CURRENT, REVISION, AGE
GRAVI_SET_VERSION_INFO([$VERSION], [0], [0], [0])

# Checks for programs.
AC_PROG_CC
AC_PROG_CPP
AM_PROG_CC_C_O
AC_PROG_LN_S

ESO_ENABLE_DEBUG(no)
ESO_ENABLE_STRICT(no)
AC_PROG_CC_C99
AC_DEFINE([_POSIX_C_SOURCE], [200809L], [Minimun version of POSIX standard to use. Needed for mkstemp, snprintf, strndup, posix_fallocate, ftruncate])
AC_DEFINE([_XOPEN_SOURCE], [700], [Minimun version of SUS standard to use]. Needed for M_PI)
AC_SYS_LARGEFILE

ESO_CHECK_DOCTOOLS

AC_ENABLE_STATIC(no)
AC_ENABLE_SHARED(yes)

LT_INIT([disable-static])

# Checks for libraries.
AC_CHECK_LIB(m, pow, [LIBS="$LIBS -lm"])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE

# Checks for library functions.
AC_CHECK_FUNCS([floor pow sqrt isinf isnan])

# Check for CPL presence and usability
CPL_CHECK_LIBS

# Check for GSL 
CPL_CHECK_GSL

#Check for ERFA
ESO_CHECK_ERFA

# mark hdrl folder as containing a configurable package
AC_CONFIG_SUBDIRS([hdrl])
# check and define all required variables to build and
# link hdrl external located in the hdrl folder
HDRL_CHECK([hdrl])

#Check the presence of pkgconfig
#ESO_PROG_PKGCONFIG

# Check for <xx>clipm
# - can be called multiple times
# - replace the argument "xxclipm" by your library name,
#   e.g. with KMOS: XXCLIPM_CHECK_LIB(kmclipm)
# Comment in if required!
# XXCLIPM_CHECK_LIB(xxclipm)

# call after CPL_CHECK_LIBS and XXCLIPM_CHECK_LIB
GRAVI_SET_PATHS 
GRAVI_CREATE_SYMBOLS

MIRA_SET_PATHS

AM_WITH_DMALLOC


AC_CONFIG_FILES(Makefile
                doxygen/Doxyfile
                gravi/Makefile
                gravi/tests/Makefile
                irplib/Makefile
                irplib/tests/Makefile
                reflex/Makefile
                reflex/gravity_wkf.xml
                reflex/gravity_viscal.xml
                reflex/gravity_mastercalib.xml
                recipes/Makefile
                mira/Makefile)
AC_OUTPUT
