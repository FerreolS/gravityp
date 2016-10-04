/* $Id: gravi_dfs.h,v 1.9 2011/04/31 06:10:40 nazouaoui Exp $
 *
 * This file is part of the GRAVI Pipeline
 * Copyright (C) 2002,2003 European Southern Observatory
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef GRAVI_IMAGE_H
#define GRAVI_IMAGE_H

/*-----------------------------------------------------------------------------
                                   Define
 -----------------------------------------------------------------------------*/

//#define MIRA_SCRIPT         "/home/vincent/Developpement/EclipseWorkspace/dev_mira/mira/mira-script.i"
//#define MIRA_SCRIPT2        GRAVI_PLUGIN_PATH "/" PACKAGE "-" PACKAGE_VERSION "/mira/mira-script.i"
//#define YORICK_BIN			"/usr/bin/yorick"
//#define YORICK_BIN2			GRAVI_PLUGIN_PATH "/" PACKAGE "-" PACKAGE_VERSION "/yorick/bin/yorick"

/*-----------------------------------------------------------------------------
                                Public prototypes
 -----------------------------------------------------------------------------*/

cpl_image * gravi_image( const cpl_frame * input_frame, const cpl_parameterlist * input_param);

#endif
