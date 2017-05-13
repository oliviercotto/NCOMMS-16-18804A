/**  $Id: updaterservices.h,v 1.6 2014-01-24 10:53:00 fred Exp $
 *
 *  @file updaterservices.h
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
 *  created on @date 30.12.2008
 *
 *  @author fred
 */

#ifndef UPDATERSERVICES_H
#define UPDATERSERVICES_H

#include <map>
#include <list>
#include "service.h"
#include "param.h"

using namespace std;

class SimComponent;
/**Class to update the simulation components' state during a simulation. It stores the
   ParamUpdater objects. Each ParamUpdater stores a list of parameters from a SimComponent
   that use a single updating function from that component. That updating function is called
   to update the state of the SimComponent that declares it at a given generation during a
   simulation.
   When attaching an updater to the UpdaterServices, it is added to the _componentUpdater map
   and the list of its parameters is stored to the _paramUpdater map. Each map has the generation
   number at which the parameter/component must be updated as a key.
   The updaters are attached during the simulation setup (SimRunner::register_component). 
   Each parameter set/simulation component is scanned for the presence of temporal parameters 
   and is asked to add its updaters to the stack. The updaters will appear only once in the 
   _updaters list as they store a single updating function although several parameters may upload
   the same updater. This allows several parameters to be attached to a single updater and use a 
   single updating function from their simulation component.
   
   The updating procedure is as follows: the parameters are first updated with their temporal
   argument value at the right generation. The components are then updated by using the function
   pointers stored in the updaters that will fetch the parameters values and update its state. 
   Those updating functions must be implemented seperately for each component and set of parameters. 
   The default SimComponent::setParameters function may suffice most of the time. That function is
   called only once and thus will update the component's state at the generations given in input by
   by the user.
*/
class UpdaterServices : public Service {

private:
  
  map< unsigned int, list< Param* > > _paramUpdater;
  map< unsigned int, list< ParamUpdaterBase* > > _componentUpdater;
  list< ParamUpdaterBase* > _updaters;
  
public:
  
  UpdaterServices () { }
  virtual ~UpdaterServices () {reset();}

  virtual bool init ();
  virtual void notify () {}
  /**The generation time is passed to the updaters.*/
  virtual void notify (unsigned int generation);
  /**Loads the updaters in case a component's parameter has temporal arguments.*/
  virtual void load ( SimComponent* sc );
  /**Attach the updater to the stack.*/
  virtual void attach ( Handler* H );
  /**Clears the containers.*/
  virtual void reset();
  /**Calls the parameters updating procedure. The map is searched for the occurence of the generation.*/
  void update_params (unsigned int generation);
  /**Update the components using the functions stored in the updaters.*/
  void update_components (unsigned int generation);
  /**Checks if any updater has been uploaded.*/
  bool hasTemporals () {return !_paramUpdater.empty();}
};

#endif
