## Process this file with automake to produce Makefile.in

##   This file is part of the GRAVI Pipeline
##   Copyright (C) 2002,2003 European Southern Observatory

##   This library is free software; you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation; either version 2 of the License, or
##   (at your option) any later version.

##   This program is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.

##   You should have received a copy of the GNU General Public License
##   along with this program; if not, write to the Free Software
##   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

AUTOMAKE_OPTIONS = 1.8 foreign

SUBDIRS = . tests

DISTCLEANFILES = *~


if MAINTAINER_MODE

MAINTAINERCLEANFILES = Makefile.in

endif


AM_CPPFLAGS = -DMIRADIR=\"$(srcdir)/\" $(CPL_INCLUDES) $(GSL_CFLAGS) $(ERFA_CFLAGS)

INCLUDES = $(all_includes)

noinst_HEADERS =    gravi_calib.h \
                    gravi_wave.h \
                    gravi_utils.h \
                    gravi_ellipse.h \
                    gravi_pfits.h \
                    gravi_dfs.h \
                    gravi_data.h \
                    gravi_image.h \
                    gravi_metrology.h \
                    gravi_acqcam.h \
                    gravi-test.h \
                    gravi_cpl.h \
                    gravi_vis.h \
                    gravi_p2vmred.h \
                    gravi_tf.h \
                    gravi_p2vm.h \
                    gravi_preproc.h \
                    gravi_signal.h \
                    gravi_eop.h

pkginclude_HEADERS =

privatelib_LTLIBRARIES = libgravi.la

libgravi_la_SOURCES =    gravi_calib.c \
	                     gravi_wave.c \
                         gravi_utils.c \
                         gravi_ellipse.c \
                         gravi_pfits.c \
                         gravi_dfs.c \
                         gravi_data.c \
                         gravi_metrology.c \
                         gravi_acqcam.c \
                         gravi_image.c \
                         gravi-test.c \
                         gravi_cpl.c \
                         gravi_vis.c \
                         gravi_p2vmred.c \
                         gravi_tf.c \
                         gravi_p2vm.c \
                         gravi_preproc.c \
                         gravi_signal.c \
                         gravi_disp.c \
                         gravi_eop.c

libgravi_la_LDFLAGS = $(CPL_LDFLAGS) $(GSL_LDFLAGS) -version-info $(LT_CURRENT):$(LT_REVISION):$(LT_AGE)
libgravi_la_LIBADD = $(LIBIRPLIB) $(LIBCPLDFS) $(LIBCPLUI) $(LIBCPLDRS) $(LIBCPLCORE) $(GSL_LIBS) $(ERFA_LIBS)
libgravi_la_DEPENDENCIES = $(LIBIRPLIB)
