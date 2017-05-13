/**  $Id: main.cc,v 1.3 2014-01-24 10:57:21 fred Exp $
*
*  @file main.cc
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
*  created on @date 19.07.2004
* 
*  @author fred
*/


#include "simulation.h"
#include "simenv.h"
#include "metapop.h"

using namespace std;

//--------------------------------------------------------------------
int main (int argc, char **argv)
{
  //--------------------------------------------------------------------
  
  Metapop *thePop = new Metapop;
  SimRunner *theSim = SIMenv::setMainSim( thePop );
  
  if( !theSim->run(argc, argv) ) fatal(" could not run the simulation!\n");
  
  delete theSim;
  delete thePop;
  
  return 0;
}
