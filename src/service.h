/** $Id: service.h,v 1.6 2014-01-24 10:53:01 fred Exp $
 *
 *  @file service.h
 *  Nemo2
 *
 *  Copyright (C) 2006-2011 Frederic Guillaume
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
 *  created on @date 16.06.2005
 * 
 *  @author fred
 */

#ifndef SERVICE_H
#define SERVICE_H

#include <string>
#include <list>
#include "handler.h"

class SimComponent;

/**Interface for the simulation services (files and stats).
 * Implements the observer pattern. Notify the observers (Handler) to update their state.
 * Contains the observer list. Provides interface to attach the observers to the service.
 **/
class Service {
  
public:
  
  Service( ) { }
  
  virtual ~Service( ) { }
  /**Inits internals. */
  virtual bool init ( ) = 0;
  /**Notifies all observers to update their state. */
  virtual void notify ( ) = 0;
  /**Interface used by a simulation component to load its obervers onto a service provider.*/
  virtual void load ( SimComponent* sc ) = 0;
  /**Adds an observer to the list. */
  virtual void attach ( Handler* h ) = 0;
  /**Clears the observers list. */
  virtual void reset ( ) = 0;
  
};

#endif //SERVICE_H
