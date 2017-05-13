/** $Id: stathandler.cc,v 1.8.2.1 2014-04-29 18:11:45 fred Exp $
*
*  @file stathandler.cc
*  Nemo2
*
*   Copyright (C) 2006-2011 Frederic Guillaume
*   frederic.guillaume@env.ethz.ch
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
*  created on @date 05.12.2005
*
*  @author fred
*/

#include <iostream>
#include "stathandler.h"
#include "statrecorder.h"
#include "metapop.h"
#include "simenv.h"

// ----------------------------------------------------------------------------------------
// init
// ----------------------------------------------------------------------------------------
void StatHandlerBase::init()
{
  _pop = get_pop_ptr();

  if(SIMenv::getGenerations() == 0)
    fatal("cannot initialize stat recorders, generation nbr is null\n");
  
  if(SIMenv::getReplicates() == 0) 
    fatal("cannot initialize stat recorders, replicate nbr is null\n");

}
// ----------------------------------------------------------------------------------------
// reset
// ----------------------------------------------------------------------------------------
void StatHandlerBase::reset ( )
{
  STAT_IT IT = _stats.begin();

  while(IT != _stats.end()) {
    delete (*IT);
    IT++;
  }

  _stats.clear();

  clear();
}

