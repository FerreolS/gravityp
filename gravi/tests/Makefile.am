## Process this file with automake to produce Makefile.in

##   This file is part of the GRAVI Pipeline Library
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

AUTOMAKE_OPTIONS = 1.6 foreign

DISTCLEANFILES = *~


if MAINTAINER_MODE

MAINTAINERCLEANFILES = Makefile.in

endif


AM_CPPFLAGS = -DDATADIR=\"$(srcdir)/\" $(all_includes)

check_PROGRAMS = gravi_dfs-test \
                 gravi_data-test \
                 gravi_vis-test \
                 gravi_utils-test \
                 gravi_calib-test \
                 gravi_strehl-test \
                 gravi_image-test \
                 gravi_eop-test

gravi_dfs_test_SOURCES = gravi_dfs-test.c
gravi_dfs_test_LDFLAGS = $(CPL_LDFLAGS) $(GSL_LDFLAGS) $(HDRL_LDFLAGS)
gravi_dfs_test_LDADD = $(LIBCPLUI) $(LIBCPLCORE) $(LIBGRAVI) $(GSL_LIBS) $(HDRL_LIBS)
gravi_dfs_test_DEPENDENCY = $(LIBHDRL)

gravi_data_test_SOURCES = gravi_data-test.c
gravi_data_test_LDFLAGS = $(CPL_LDFLAGS) $(HDRL_LDFLAGS)
gravi_data_test_LDADD = $(LIBCPLUI) $(LIBCPLCORE) $(LIBGRAVI) $(GSL_LIBS) $(HDRL_LIBS)
gravi_data_test_DEPENDENCY = $(LIBHDRL)

gravi_calib_test_SOURCES = gravi_calib-test.c
gravi_calib_test_LDFLAGS = $(CPL_LDFLAGS) $(GSL_LDFLAGS) $(HDRL_LDFLAGS)
gravi_calib_test_LDADD = $(LIBCPLUI) $(LIBCPLCORE) $(LIBGRAVI) $(GSL_LIBS) $(HDRL_LIBS)
gravi_calib_test_DEPENDENCY = $(LIBHDRL)

gravi_vis_test_SOURCES = gravi_vis-test.c
gravi_vis_test_LDFLAGS = $(CPL_LDFLAGS) $(GSL_LDFLAGS) $(HDRL_LDFLAGS)
gravi_vis_test_LDADD = $(LIBCPLUI) $(LIBCPLCORE) $(LIBGRAVI) $(GSL_LIBS) $(HDRL_LIBS)
gravi_vis_test_DEPENDENCY = $(LIBHDRL)

gravi_utils_test_SOURCES = gravi_utils-test.c
gravi_utils_test_LDFLAGS = $(CPL_LDFLAGS) $(GSL_LDFLAGS) $(HDRL_LDFLAGS)
gravi_utils_test_LDADD = $(LIBCPLUI) $(LIBCPLCORE) $(LIBGRAVI) $(GSL_LIBS) $(HDRL_LIBS)
gravi_utils_test_DEPENDENCY = $(LIBHDRL)

gravi_strehl_test_SOURCES = gravi_strehl-test.c
gravi_strehl_test_LDFLAGS = $(CPL_LDFLAGS) $(GSL_LDFLAGS) $(HDRL_LDFLAGS)
gravi_strehl_test_LDADD = $(LIBCPLUI) $(LIBCPLCORE) $(LIBGRAVI) $(GSL_LIBS) $(HDRL_LIBS)
gravi_strehl_test_DEPENDENCY = $(LIBHDRL)

gravi_image_test_SOURCES = gravi_image-test.c
gravi_image_test_LDFLAGS = $(CPL_LDFLAGS) $(GSL_LDFLAGS) $(HDRL_LDFLAGS)
gravi_image_test_LDADD = $(LIBCPLUI) $(LIBCPLCORE) $(LIBGRAVI) $(GSL_LIBS) $(HDRL_LIBS)
gravi_image_test_DEPENDENCY = $(LIBHDRL)

gravi_eop_test_SOURCES = gravi_eop-test.c
gravi_eop_test_LDFLAGS = $(CPL_LDFLAGS) $(GSL_LDFLAGS) $(HDRL_LDFLAGS)
gravi_eop_test_LDADD = $(LIBCPLUI) $(LIBCPLCORE) $(LIBGRAVI) $(GSL_LIBS) $(HDRL_LIBS)
gravi_eop_test_DEPENDENCY = $(LIBHDRL)

# Be sure to reexport important environment variables.
TESTS_ENVIRONMENT = MAKE="$(MAKE)" CC="$(CC)" CFLAGS="$(CFLAGS)" \
        CPPFLAGS="$(CPPFLAGS)" LD="$(LD)" LDFLAGS="$(LDFLAGS)" \
        LIBS="$(LIBS)" LN_S="$(LN_S)" NM="$(NM)" RANLIB="$(RANLIB)" \
        OBJEXT="$(OBJEXT)" EXEEXT="$(EXEEXT)"

TESTS = $(check_PROGRAMS)

# Expected failures
XFAIL_TESTS =

# We need to remove any files that the above tests created.
clean-local:
	-rm test_data_*.fits
	

EXTRA_DIST = data1.oifits \
Dark.fits \
Flat_1.fits \
Flat_2.fits \
Flat_3.fits \
Flat_4.fits \
Wave.fits \
gravi_badpix_map.fits \
gravi_profile_map.fits \
test_compute_dark.fits \
gravi_dark_map.fits \
gravi_wave_map.fits \
gravi_p2vm_map.fits \
p2vm_map.fits \
vis_data.fits

