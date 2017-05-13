/**  $Id: lifecycleevent.h,v 1.6 2013-02-01 10:11:32 fred Exp $
*
*  @file lifecycleevent.h
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
*  created on @date 07.07.2004
*
*  @author fred
*/

#ifndef LIFECYCLEEVENT_H
#define LIFECYCLEEVENT_H

#include <string>
#include <list>
#include <map>
#include "param.h"
#include "simcomponent.h"
#include "metapop.h"

/**Base class of the Life Cycle Events, declares the LCE interface.
 * A Life Cycle Event (LCE) is a population operator that modifies the state of the population at each iteration of the life cycle.
 * It declares the execute() method that is called during the life cycle loop (Metapop::Cycle()).
 * No interface is given to manage the periodicity of the event, the execute() method is called at each iteration.
 * The periodicity check must thus be implemented in the derived class, within the exectute() method.
 *
 * Each LCE has a link to the current Metapop instance and a link to a particular trait which are set by LifeCycleEvent::init()
 * at simulation initialization time (called by Metapop::setLifeCycle() through Metapop::init() and SimRunner::setup()).
 * The linked trait can be specified in the constructor by its type identifier (a character string literal) or later,
 * at runtime. In the first case, the trait's index will be set by a call to IndFactory::getTraitIndex() at initialization
 * (see LifeCycleEvent::init()). In the later case, the trait to link with must be specified at initialization time and LifeCycleEvent::init()
 * must be overloaded. The linked trait type (or name) can then be specified by the user in the init file and set in the derived class init()
 * method before explicitly calling LifeCycleEvent::init() which will set the linked trait's index (see LCE_Breed_Selection::init() for an example).
 * That index is then used to easily access the given trait through the Individual interface (see Individual::getTrait()).
 *
 * Each LCE also has a name by which it is called from the user's defined parameter file to specify when in the life cycle that particular event
 * has to be executed. The position in the life cycle is the LCE's rank (i.e. the value given in input to the parameter "name"). As only one 
 * rank value is allowed per LCE, it is executed only once in the life cycle, this is a limitation that might be removed in future 
 * versions if need be!
 *
 * The main action of an LCE is often to move individuals between populations, between age-classes or both. The changes done to the population
 * age structure by a particular LCE are tracked through the removeAgeClass() and addAgeClass() interface methods called after each execution
 * of an LCE. The metapopulation age flag is thus updated by these two functions.
 *
 * \b Note: Overloading LifeCycleEvent::init().
 * If more than one trait is used by your LCE, think about overloading the init() function to allow the 
 * registration of several trait links within your derived LCE. The init() method will also be overloaded if you 
 * add additional parameters to the LCE's ParamSet. In any cases, LifeCycleEvent::init() MUST be called in your init() function
 * to properly set the metapop pointer (unless you copy the code into your function...).
 
 * \b Note: the LCEs \a must be inited \a after a call to IndFactory::makePrototype() has been issued so that the trait links can be set!
 */
class LifeCycleEvent :  public virtual SimComponent
{
protected:  
  /**The param name to be read in the init file.*/
  std::string _event_name;
  /**The ptr to the current Metapop. */
  Metapop* _popPtr;
  /**The name of the linked trait.*/
  std::string _LCELinkedTraitType;
  /**The index in the individual's trait table of the linked trait. A value of -1 means the link is broken.*/
  int _LCELinkedTraitIndex;
  
public: 
  /**Cstor. 
    @param name the name of the LCE as it must appear in the parameter input file
    @param trait_link the name of the linked trait used by this LCE
    @callgraph
  */
  LifeCycleEvent (const char* name, const char* trait_link) 
  : _popPtr(0), _LCELinkedTraitType(trait_link), _LCELinkedTraitIndex(-1)
  {
//    cout << "calling LifeCycleEvent("<<name<<","<<trait_link<<")\n";
    set_event_name(name);
  }

  virtual ~LifeCycleEvent ( ) { }
  /**Sets the pointer to the current Metapop and the trait link if applicable.
   * DEV: Don't forget to explicitly call this function (\c LifeCycleEvent::init()) when you overload it!
   * @param popPtr the pointer to the current instance of \c Metapop
   * @callgraph
   */
  virtual void init(Metapop* popPtr)
  { 
    _popPtr = popPtr;
    if(!attach_trait(_LCELinkedTraitType)) fatal("bailing out\n");
    if(!setParameters()) fatal("bailing out\n");
  }
  
  virtual bool attach_trait (string trait)
  {
    _LCELinkedTraitType = trait;
    
    if(_LCELinkedTraitType.size() != 0) {
      _LCELinkedTraitIndex = _popPtr->getTraitIndex(_LCELinkedTraitType.c_str());
      if(_LCELinkedTraitIndex == -1) {
        error("cannot attach trait \"%s\" to life cycle event \"%s\", trait has not been initiated.\n",
              _LCELinkedTraitType.c_str(), _event_name.c_str());
        return false;
      }
    }
    return true;
  }
  
  virtual void set_paramset (std::string name, bool required, SimComponent* owner )
  {
    SimComponent::set_paramset(name, required, owner);
    add_parameter(name.c_str(),INT,1,0,0,0,0);
  }
  
  /**Set the name of the event (name of the ParamSet) and add the corresponding parameter to the set.   */
  virtual void set_event_name (std::string& name)
  {
    _event_name = name;
    set_paramset(name, 0, this);
  }
  
  virtual void set_event_name (const char* name)
  {
    _event_name = name;
    set_paramset(name, 0, this); 
  }
  
  /**Accessor to the LCE's name.*/
  virtual string& get_event_name ( ) 
  { return _event_name; }
  
  /**Accessor to the LCE rank in the life cycle.*/

  virtual int get_rank ( )
  { return (int)get_parameter_value(_event_name.c_str()); }
  
  /**Accessors for the population pointer
   * @param popPtr The pointer to the current Metapop */
  virtual void set_pop_ptr (Metapop* popPtr)
  {_popPtr = popPtr;}
  
  virtual Metapop*   get_pop_ptr ( ) 
  { return _popPtr; }
  
  ///@name LCE interface
  ///@{
  /**Execute the event on the pop. */
  virtual void  execute () = 0;
  /**Cloning interface. */
  virtual LifeCycleEvent*  clone () = 0;
  /**Removes the returned age-class flag(s) from the current Metapop age-class flags.*/
  virtual age_t removeAgeClass () = 0;
  /**Adds the returned age-class flag(s) to the current Metapop age-class flags.*/
  virtual age_t addAgeClass () = 0;
  /**Specifies what age-classes are required by the LCE to execute.*/
  virtual age_t requiredAgeClass () = 0;
  ///@}
};

#endif //LIFECYCLEEVENT_H

