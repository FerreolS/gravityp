/* $Id: gravi_ellipse.c,v 1.10 2011/05/31 06:10:40 nazouaoui Exp $
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

#ifndef GRAVI_ELLIPSE_H
#define GRAVI_ELLIPSE_H

/*-----------------------------------------------------------------------------
                                   Includes
 -----------------------------------------------------------------------------*/

#include <cpl.h>

/*-----------------------------------------------------------------------------
                                Macros
 -----------------------------------------------------------------------------*/

/* Use linear ellipse fitting or lmvq fit */
#define USE_LINEAR_ELLIPSE 1

/*-----------------------------------------------------------------------------
                              Public prototypes
 -----------------------------------------------------------------------------*/

cpl_vector * gravi_ellipse_phase_create (cpl_vector * vectCA,
                                         cpl_vector * vectDB,
                                         cpl_vector * envelope_vector);

cpl_vector * gravi_ellipse_phase_create_fast (cpl_vector * vectCA,
                                              cpl_vector * vectDB);

cpl_vector * gravi_ellipse_meanopd_create (cpl_table * spectrum_table,
                                           cpl_table * detector_table,
                                           cpl_table ** oiwave_tables,
                                           cpl_vector * guess_vector,
                                           int base);

#endif
