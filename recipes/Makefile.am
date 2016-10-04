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

DISTCLEANFILES = *~


if MAINTAINER_MODE

MAINTAINERCLEANFILES = Makefile.in

endif


INCLUDES = $(all_includes)

noinst_HEADERS = 

plugin_LTLIBRARIES = gravity_dark.la \
						gravity_p2vm.la \
						gravity_vis.la \
						gravity_postprocess.la \
						gravity_viscal.la \
						gravity_eop.la \
						gravity_image.la \
						gravity_nab.la \
						gravity_disp.la \
						gravity_piezo.la \
						gravity_vis_from_p2vmred.la \
						gravity_badpix.la \
						gravity_biasmask.la \
						gravity_wavelamp.la

gravity_dark_la_SOURCES = gravity_dark.c 
gravity_dark_la_LIBADD = $(LIBGRAVI)
gravity_dark_la_LDFLAGS = -module -avoid-version
gravity_dark_la_DEPENDENCIES = $(LIBGRAVI)

gravity_disp_la_SOURCES = gravity_disp.c 
gravity_disp_la_LIBADD = $(LIBGRAVI)
gravity_disp_la_LDFLAGS = -module -avoid-version
gravity_disp_la_DEPENDENCIES = $(LIBGRAVI)

gravity_piezo_la_SOURCES = gravity_piezo.c 
gravity_piezo_la_LIBADD = $(LIBGRAVI)
gravity_piezo_la_LDFLAGS = -module -avoid-version
gravity_piezo_la_DEPENDENCIES = $(LIBGRAVI)

cgravity_wavelamp_la_SOURCES = gravity_wavelamp.c 
gravity_wavelamp_la_LIBADD = $(LIBGRAVI)
gravity_wavelamp_la_LDFLAGS = -module -avoid-version
gravity_wavelamp_la_DEPENDENCIES = $(LIBGRAVI)

gravity_p2vm_la_SOURCES = gravity_p2vm.c 
gravity_p2vm_la_LIBADD = $(LIBGRAVI)
gravity_p2vm_la_LDFLAGS = -module -avoid-version
gravity_p2vm_la_DEPENDENCIES = $(LIBGRAVI)

gravity_postprocess_la_SOURCES = gravity_postprocess.c 
gravity_postprocess_la_LIBADD = $(LIBGRAVI)
gravity_postprocess_la_LDFLAGS = -module -avoid-version
gravity_postprocess_la_DEPENDENCIES = $(LIBGRAVI)

gravity_viscal_la_SOURCES = gravity_viscal.c 
gravity_viscal_la_LIBADD = $(LIBGRAVI)
gravity_viscal_la_LDFLAGS = -module -avoid-version
gravity_viscal_la_DEPENDENCIES = $(LIBGRAVI)

gravity_vis_la_SOURCES = gravity_vis.c 
gravity_vis_la_LIBADD = $(LIBGRAVI)
gravity_vis_la_LDFLAGS = -module -avoid-version
gravity_vis_la_DEPENDENCIES = $(LIBGRAVI)

gravity_vis_from_p2vmred_la_SOURCES = gravity_vis_from_p2vmred.c
gravity_vis_from_p2vmred_la_LIBADD = $(LIBGRAVI)
gravity_vis_from_p2vmred_la_LDFLAGS = -module -avoid-version
gravity_vis_from_p2vmred_la_DEPENDENCIES = $(LIBGRAVI)

gravity_badpix_la_SOURCES = gravity_badpix.c
gravity_badpix_la_LIBADD = $(LIBGRAVI)
gravity_badpix_la_LDFLAGS = -module -avoid-version
gravity_badpix_la_DEPENDENCIES = $(LIBGRAVI)

gravity_biasmask_la_SOURCES = gravity_biasmask.c
gravity_biasmask_la_LIBADD = $(LIBGRAVI)
gravity_biasmask_la_LDFLAGS = -module -avoid-version
gravity_biasmask_la_DEPENDENCIES = $(LIBGRAVI)

gravity_nab_la_SOURCES = gravity_nab.c 
gravity_nab_la_LIBADD = $(LIBGRAVI)
gravity_nab_la_LDFLAGS = -module -avoid-version
gravity_nab_la_DEPENDENCIES = $(LIBGRAVI)

gravity_image_la_SOURCES = gravity_image.c 
gravity_image_la_LIBADD = $(LIBGRAVI)
gravity_image_la_LDFLAGS = -module -avoid-version
gravity_image_la_DEPENDENCIES = $(LIBGRAVI)

gravity_eop_la_SOURCES = gravity_eop.c 
gravity_eop_la_LIBADD = $(LIBGRAVI)
gravity_eop_la_LDFLAGS = -module -avoid-version
gravity_eop_la_DEPENDENCIES = $(LIBGRAVI)
