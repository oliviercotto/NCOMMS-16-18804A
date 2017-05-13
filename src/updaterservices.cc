/**  $Id: updaterservices.cc,v 1.6 2014-01-24 10:53:01 fred Exp $
 *
 *  @file updaterservices.cc
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

#include <deque>
#include "updaterservices.h"
#include "simcomponent.h"
#include "lifecycleevent.h"
#include "param.h"

//-----------------------------------------------------------------------------------------
// init
// ----------------------------------------------------------------------------------------
bool UpdaterServices::init ( )
{
  //consolidate the lists of updaters by removing the duplicates (if they exist)
  for(map<unsigned int, list<ParamUpdaterBase*> >::iterator it = _componentUpdater.begin();
      it != _componentUpdater.end(); it++)
    it->second.unique();
  
  for(map<unsigned int, list<Param*> >::iterator it = _paramUpdater.begin();
      it != _paramUpdater.end(); it++)
    it->second.unique();
    
  _updaters.unique();
  
  for(list<ParamUpdaterBase*>::iterator it = _updaters.begin();
      it != _updaters.end(); it++)
    (*it)->init();
  
  return true;
}
//-----------------------------------------------------------------------------------------
// reset
// ----------------------------------------------------------------------------------------
void UpdaterServices::reset ( )
{
  map<unsigned int, list<Param*> >::iterator pit = _paramUpdater.begin();
  
  for(;pit != _paramUpdater.end(); pit++) {
    pit->second.clear();
  }
  _paramUpdater.clear();
  
  map<unsigned int, list<ParamUpdaterBase*> >::iterator scit = _componentUpdater.begin();
  
  for(;scit != _componentUpdater.end(); scit++) {
    scit->second.clear();
  }
  _componentUpdater.clear();
  
  _updaters.clear();
}
//-----------------------------------------------------------------------------------------
// load
// ----------------------------------------------------------------------------------------
void UpdaterServices::load ( SimComponent* sc )
{
  sc->loadUpdaters(this);
}
//-----------------------------------------------------------------------------------------
// attach
// ----------------------------------------------------------------------------------------
void UpdaterServices::attach ( Handler* H )
{
  ParamUpdaterBase* PU = dynamic_cast<ParamUpdaterBase*> (H);
  
  if(!PU) return;
    
  list< Param* > params = PU->getParams();
  
  deque< unsigned int > dates;
    
  for(list<Param*>::const_iterator pit = params.begin();
      pit != params.end(); pit++) {
    
    dates = (*pit)->getUpdatingDates();
    
    for(unsigned int i = 0; i < dates.size(); i++) {
      _paramUpdater[ dates[i] ].push_back((*pit));
      _componentUpdater[ dates[i] ].push_back(PU);
    }
  }
  
  _updaters.push_back(PU);
}
//-----------------------------------------------------------------------------------------
// notify
// ----------------------------------------------------------------------------------------
void UpdaterServices::notify (unsigned int generation)
{
  update_params(generation);
  update_components(generation);
}
//-----------------------------------------------------------------------------------------
// update_params
// ----------------------------------------------------------------------------------------
void UpdaterServices::update_params (unsigned int generation)
{
  map<unsigned int, list<Param*> >::iterator mit = _paramUpdater.find(generation);
  
  if(mit == _paramUpdater.end()) return;
  
  for(list<Param*>::iterator pit = mit->second.begin(); 
      pit != mit->second.end(); pit++)
    (*pit)->update(generation);
  
}
//-----------------------------------------------------------------------------------------
// update_components
// ----------------------------------------------------------------------------------------
void UpdaterServices::update_components (unsigned int generation)
{
  map<unsigned int, list<ParamUpdaterBase*> >::iterator mit = _componentUpdater.find(generation);
  
  if(mit == _componentUpdater.end()) return;
  
  for(list<ParamUpdaterBase*>::iterator pit = mit->second.begin(); 
      pit != mit->second.end(); pit++) 
  {
    if( !(*pit)->update(generation) )
      warning("could not update component %s at generation %i.\n", 
              (*pit)->getComponent()->get_name().c_str(), generation);
    }
}

