/** $Id: stathandler.h,v 1.8.2.2 2014-04-29 18:11:11 fred Exp $
 *
 *  @file stathandler.h
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
 *  created on @date 07.04.2005
 * 
 *  @author fred
 */

#ifndef STATHANDLER_H
#define STATHANDLER_H

#include <list>
#include <string>
#include <iostream>
#include "handler.h"
#include "statservices.h"

class Metapop;

class StatRecBase;
/**Base class of the StatHandler class, implements the Handler interface.
 * This class stores the Handler state and the list of links to the StatRecBase.
 * This list is duplicated in the StatHandler as a list of the StatRecorder templates. 
 * It allows the StatService to access the stat recorders without knowledge of the
 * actual StatHandler and StatRecorder template type. 
 */
class StatHandlerBase : public Handler {
  
private:
  /**Link to the StatRecorder list elements in the StatHandler derived class.*/
  std::list<StatRecBase*> _stats;

  /**Link to the StatService.*/
  StatServices*      _service;
  
  typedef std::list< StatRecBase* >::iterator STAT_IT;
  
protected:
  /**Link to the current population, set through the link to the StatService.*/
  Metapop* _pop;
  
public:
  
  StatHandlerBase ( ) :  _service(0) {}
  
  virtual           ~StatHandlerBase     ( ) { }

  ///@name Accessors
  ///@{
  Metapop*          get_pop_ptr      ( )                 {return _service->get_pop_ptr();}
  
  void              set_service      (StatServices* srv) {_service = srv;}
  
  StatServices*     get_service      ( )                 {return _service;}
  
  unsigned int      getOccurrence    ( )                 {return _service->getOccurrence();}
  
  unsigned int      getNumOccurrences( )                 {return _service->getNumOccurrences();}
  
  unsigned int      getCurrentOccurrence ( )             {return _service->getCurrentOccurrence();}
  
  unsigned int      getNbRecorders   ( )                 {return _stats.size();}
    
  std::list<StatRecBase*>& getStats  ( )                 {return _stats;}
  
  virtual void      add              (StatRecBase* rec)  {_stats.push_back(rec);}
  ///@}
  /**Empties the _stats list and calls clear() (defined in the derived class).*/
  virtual void      reset            ( );
  
  ///@name Handler implementation
  ///@{
  virtual void      init             ( );
  /**This function is left empty as the StatServices calls StatRecorder::setVal directly.*/
  virtual void      update           ( ) { }
  ///@}
  ///@name StatHandler interface declaration
  ///@{
  
  virtual bool      setStatRecorders (std::string& token) = 0;
  
  virtual void      clear            ( ) = 0;
  ///@}
};

/**A class to compute and store the summary statistics associated with a SimComponent.
 * The template type must be the type of the class that declares the methods linked into the
 * StatRecorder elements. */
template <class SH> class StatHandler: public StatHandlerBase {
  
protected:
  /**The list of stat recorders.*/
  std::list<StatRecorder<SH>*> _recorders;
  
  typedef typename std::list< StatRecorder<SH>* >::iterator REC_IT;
  
public:
  
  StatHandler( ) { }
  
  virtual ~StatHandler ( ) { }
  
  /**Empties the _recorders list, they are destroyed in StatHandlerBase::reset().*/
  virtual void clear   ( )  {_recorders.clear();}

  /**Adds a StatRecorder to the list, it is also added to the StatHandlerBase::_stats list.
   Two types of function variables are passed to this function. The "getter" and "setter".
   A "getter" returns a double value that will be stored in the StatRecorder structure. It 
   may or may not take an argument. Only one getter should be passed to the stat recorder.
   The setter is used to set variables in the SH class which are then read by the getter.
   A setter and a getter may be given together but a setter alone will issue an error at
   runtime as no getter is present in the stat recorder.
   * @param Title the stat title (as shown in the R frontend fro eg.)
   * @param Name the stat name (headers in the stat output file)
   * @param AGE the age class for which the stat will be recorded
   * @param ARG1 the first argument to pass to the SH function
   * @param ARG2 the second argument to pass to the SH function
   * @param getStatNoArg function ptr to a SH getter that takes no argument
   * @param getStatOneArg function ptr to a SH getter with a single integer argument
   * @param getStatTwoArg function ptr to a SH getter with two integer arguments
   * @param setStat function ptr to a SH setter
   */
  virtual StatRecorder<SH>* add (std::string Title, std::string Name, age_t AGE, unsigned int ARG1, unsigned int ARG2,
                                 double(SH::* getStatNoArg)(void), double(SH::* getStatOneArg)(unsigned int),
                                 double(SH::* getStatTwoArg)(unsigned int, unsigned int), void(SH::* setStat)(void))
  {
    StatRecorder<SH>* new_rec = new StatRecorder<SH>();

    new_rec->set(Title, Name, AGE, ARG1, ARG2, 
                 getStatNoArg, getStatOneArg, getStatTwoArg, setStat);
    
    new_rec->setHandler(dynamic_cast<SH*> (this));
    
    _recorders.push_back(new_rec);
    
    StatHandlerBase::add(new_rec);
    
    return new_rec;
  }
};

// TraitStatHandler
//
/**Template class for the trait's StatHandler. Constructor links to a given trait prototype.
 The pointer and index of the trait prototype are stored in the class. They
 can be accessd through the _SHLinkedTrait and _SHLinkedTraitIndex members to get the stats.*/
template <class TP, class SH> class TraitStatHandler : public StatHandler<SH> {
protected:
  /**Pointer to a TraitProtoype object.*/
  TP* _SHLinkedTrait;
  /**Index of the trait in the Individual::Traits table.*/
  int _SHLinkedTraitIndex;
public:
  TraitStatHandler (TP* trait_proto);
  virtual ~TraitStatHandler () { }
};

template<class TP, class SH> TraitStatHandler<TP,SH>::TraitStatHandler(TP* trait_proto)
{
  _SHLinkedTrait = trait_proto;
  _SHLinkedTraitIndex = trait_proto->get_index();
}

// EventStatHandler
//
/**Template class for the LCEs StatHandler classes. Constructor links to a give LCE.
 The LCE can be accessed through the _SHLinkedEvent member to get the stats.*/
template <class LCE, class SH> class EventStatHandler : public StatHandler<SH> {
protected:
  /**Pointer to the linked LCE.*/
  LCE* _SHLinkedEvent;
public:
  EventStatHandler (LCE* lce);
  virtual ~EventStatHandler () { }
};

template<class LCE, class SH> EventStatHandler<LCE,SH>::EventStatHandler(LCE* lce)
{
  _SHLinkedEvent = lce;
}


#endif
