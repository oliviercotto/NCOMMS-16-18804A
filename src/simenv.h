/**  $Id: simenv.h,v 1.6.2.2 2014-04-29 18:16:55 fred Exp $
 *
 *  @file simenv.h
 *  Nemo2
 *
 *   Copyright (C) 2008-2011 Frederic Guillaume
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
 *  created on @date 31.12.2008
 * 
 *  @author fred
 */

#ifndef SIMENV_H
#define SIMENV_H

#include "simulation.h"

/**Global class exposing the current SimRunner throughout the code. **/
class SIMenv {
  
  SIMenv(){}
  
public:
  
  static SimRunner* MainSim;
  
  static void loadDefaultComponents(SimRunner* sim);
  
  static SimRunner* setMainSim ( Metapop *pop ) {
    MainSim = new SimRunner( pop );
    loadDefaultComponents( MainSim );
    return MainSim;
  }
  
  static SimRunner* getNewSimulation ( ) {
    SimRunner *sim = new SimRunner( new Metapop() );
    loadDefaultComponents( sim );
    return sim;
  }
  
  static void setMainSim( SimRunner* sim) {MainSim = sim;}
  
  static unsigned int getCurrentGeneration() {return MainSim->getCurrentGeneration();}
  static unsigned int getGenerations() {return MainSim->getGenerations();}
  static unsigned int getCurrentReplicate() {return MainSim->getCurrentReplicate();}
  static unsigned int getReplicates() {return MainSim->getReplicates();}
  static int getCurrentRankInLifeCycle() {return MainSim->getCurrentRankInLifeCycle();}
};

#endif
