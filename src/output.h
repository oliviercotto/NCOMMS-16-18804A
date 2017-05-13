/** $Id: output.h,v 1.6.2.2 2016-04-28 13:07:30 fred Exp $
*
*  @file output.h
*  Nemo2
*
*   Copyright (C) 2006-2015 Frederic Guillaume
*   frederic.guillaume@ieu.uzh.ch
*
*   This file is part of Nemo
*
*   Nemo is free software; you can redistribute it and/or modify
*   it under the terms of the GNU General Public License as published by
*   the Free Software Foundation; either version 2 of the License, or
*   (at your option) any later version.
*
*   Nemo is distributed in the hope that it will be useful,
*   but WITHOUT ANY WARRANTY; without even the implied warranty of
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*   GNU General Public License for more details.
*
*   You should have received a copy of the GNU General Public License
*   along with this program; if not, write to the Free Software
*   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*
*  created on @date 08.09.05
* 
*  @author fred
*/

#ifndef OUTPUTS_H
#define OUTPUTS_H

#include <cstdlib>

#ifdef _R_OUTPUT_
#include <R.h>
#endif

extern bool SILENT_RUN;

extern void message (const char* message, ...);
extern void warning (const char* message, ...);
extern int  error   (const char* message, ...);
extern void fatal   (const char* message, ...);

#endif
