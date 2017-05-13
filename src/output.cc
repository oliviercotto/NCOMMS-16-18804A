/** $Id: output.cc,v 1.8.2.2 2016-03-09 15:19:24 fred Exp $
*
*  @file output.cc
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

#include <iostream>
#include <string>
#include <stdarg.h>
#include <stdio.h>
#include "output.h"
#include "MPImanager.h"

bool SILENT_RUN = false;


void message (const char* message, ...)
{
  
  if(SILENT_RUN) return;
  
  va_list(ap);
  
  va_start(ap, message);
#ifdef _R_OUTPUT_
  Rvprintf(message, ap);
#else
  vprintf(message,ap);
#endif
  va_end(ap);
}

void warning (const char* str, ...)
{
  if(SILENT_RUN) return;
  
  va_list(ap);
  
  std::cout<<"\n>>>> WARNING::";
   
  va_start(ap, str);
#ifdef _R_OUTPUT_
  Rvprintf(str, ap);
#else
  vprintf(str,ap);
#endif
  va_end(ap);
}

int error (const char* str, ...)
{
  va_list(ap);
  
  std::cerr<<"***ERROR*** ";
 
  va_start(ap, str);
#ifdef _R_OUTPUT_
  REvprintf(str, ap);
#else
  vfprintf(stderr,str, ap);
#endif
  va_end(ap);
  
  return 0;
}

void fatal (const char* str, ...)
{
  va_list(ap);
  
  std::cerr<<"***ERROR*** ";

  va_start(ap, str);
#ifdef _R_OUTPUT_
  REvprintf(str, ap);
#else
  vfprintf(stderr,str,ap);
#endif
  va_end(ap);

  MPIenv::abort(1);
}
