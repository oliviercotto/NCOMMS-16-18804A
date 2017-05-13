/** $Id: basicsimulation.h,v 1.7.2.2 2014-04-29 16:13:14 fred Exp $
*
*  @file basicsimulation.h
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
*  Created on @date 30.10.2005
*  @author fred
**/

#ifndef BASICSIMULATION_H
#define BASICSIMULATION_H

#include <list>
#include <map>
#include <string>
#include "ttrait.h"
#include "lifecycleevent.h"

using namespace std;
/* /_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/ */
/** Class to manage the simulation components and prototypes. 
 *  This class stores and provides accessors to the simulation components. It also stores
 *  the trait prototype and life cycle event template instances. 
 **/
/* \_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\ */

class ComponentManager {

public:
  ComponentManager ( )  {}
  ~ComponentManager ( ) {/*message("ComponentManager::~ComponentManager\n");*/}
  /**Clears and builds the list of all components from the lists of trait prototype templates and life cycle event templates. **/
  void                             build_component_list   ( );
  /**Push a component at the back of the component list. **/
  void                             add_component          (SimComponent* cmpt)     {_components.push_back(cmpt);}
  /**Add a trait prototype to the template and component lists. **/
  void                             add_trait     (TraitPrototype* trait) {_TTrait_Templates.push_back(trait); _components.push_back(trait);}
  /**Add a life cycle event to the template and component lists. **/
  void                             add_LCE       (LifeCycleEvent* event) {_LCE_Templates.push_back(event);_components.push_back(event);}
  /**Search for component with "name" in the trait prototype list. 
   * @return NULL if component with "name" not found, a pointer to that component otherwise
   **/
  TraitPrototype*                  get_trait     (string name); 
  /**Search for component with "name" in the life cycle events list. 
    * @return NULL if component with "name" not found, a pointer to that component otherwise
    **/
  LifeCycleEvent*                  get_LCE       (string name); 
  
protected:
  
  /**List of all the simulation components. **/
  list< SimComponent* >            _components;
  /**List of all trait prototypes of the simulation, a subset of _components list. */
  list< TraitPrototype* >          _TTrait_Templates;
  /**List of all the life-cycle events of the simulation, a subset of _components list. */
  list< LifeCycleEvent* >          _LCE_Templates;
};

/* /_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/ */
/**Class to manage the sets of parameters of the simulation components. 
 *  This class performs parameters setting and checking for the whole set of the simulation
 *  components. Provides access to derived classes to the complete list of parameter sets. 
 *  Also sets the list of simulations parameters in case of sequential parameters found in 
 *  input.
 *  It stores and builds the simulation parameters set.
 *
 *  @see ParamsParser and derived classes. These are the input parameters providers.
 **/
/* \_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\ */
class ParamManager: public ComponentManager {
public:
  /**Cstor. Builds the simulation PramaSet. **/
  ParamManager ( );
  ~ParamManager ( ) {/*message("ParamManager::~ParamManager\n");*/}
  /**Adds a ParamSet to the list of the parameter sets of the simulation. **/
  void                         add_paramset           (ParamSet* paramset) {_allParams.push_back(paramset);}
  /**Looks for paramset with "name" in the list of parameter sets. 
   * @return NULL if parameter set "name" not found
  **/
  ParamSet*                    get_paramset           (string& name); 
  /**Clears and fills the _allParams list with the ParamSet's of the simulation components. **/
  void                             build_allParams        ( );
  /**Accessor of the whole list of the parameter sets. 
    *@return the reference to the list of ParamSet
   **/
  list<ParamSet*>&             get_allParams          ( )                      {return _allParams;}
  
  /**Sets the parameters of the simulation with the argument strings given in input. 
    *Scans the _allParams list to set the parameters present in the input map simparams. Each ParamSet checks 
    *internally for the presence of a Param with the given name string and sets its value with the given argument, if present.
    
    \b Note: all ParamSet owning a Param with the same name will use the same argument string. The input map is not a
    * multimap, each param name is present only once.
    *
    *@param simparams a map containing the parameter names and their argument string
    *@param silent will be silent about parameters that could not be set.
    *@return the status of the \c param_consistency_check() function
    *@callgraph
    **/
  bool                             set_parameters         (map< string,string >& simparams, bool silent);
  
  /**Checks if all the mandatory parameters are set so that the simulation can be launched. 
    *@return TRUE if all \c ParamSet::check_consistency() returned true
    *@callgraph
    **/
  bool                             param_consistency_check ( );
  /**Builds the list of simulation parameters from the parsed input file(s). @callgraph**/
  void                             build_records          (map< string, vector<string> >& initParams);
  /**Accessor to the simulations parameter list. **/
  list< map< string,string > >&    get_simRecords         ( )                      {return _simRecords;}
  /**Accessor to the first element in the simulations parameter list. **/
  map< string,string >&            get_firstRecord        ( )                      {return (*_simRecords.begin());}
  /**Accessor to the size of the simulations parameter list, i.e. the number of simulations to perform. **/
  int                              get_nbSims             ( )                      {return _simRecords.size();}
  
protected:
  /**A list of all the parameter sets of all the simulation components loaded in the _component list of the ComponentManager. **/
  list< ParamSet* >              _allParams;
  /**A map of the parameters and their arguments of the current (running) simulation. **/
  map< string, string >          _inputParams;
  /**Lists of parameters to be updated during a simulation indexed by generation update time.*/
  map< unsigned int, list < pair< string, string> > >   _temporalParams;
  /**Sets of parameters of all the simulations to perform. **/
  list< map< string, string > >  _simRecords;
  /**The ParamSet param set of the simulation. **/
  ParamSet                       _paramSet;
  
private:
  
  string                           setFilename            (string& fstring, unsigned int sim, vector<string>& args, vector<unsigned int>& arg_no, bool check_arg_no);
  string                           stripFormatString      (string& str, unsigned int& index);
  string                           setArgString           (string& fmt, string& arg, unsigned int arg_pos);
  string                           lowercase              (string& input);

};
/* /_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/ */
/**Provides methods to build the user's selected set of life cycle events and traits from the parameters.
 *  This class implements methods to build the lists of selected traits and life cycle events from the user's defined parameters.
 *  Each simulation component that has its ParamSet in the "set" state is elligible to be part of the current simulation.
 *  Accessors to these components are provided. This class does however not provide a runnable simulation object.
 **/
/* \_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\ */
class SimBuilder: public ParamManager {
public:
  SimBuilder ( ) { }
  /**copy cstor.*/
  SimBuilder (const SimBuilder& SB);
  ~SimBuilder ( );
  /**Builds the list of parameters from user's defined input parameters.
   * @param simparams Hashtable of the parsed input parameters
   * @return TRUE if parameters are consitently set
   **/
  bool                             build_currentParams    (map< string,string >& simparams);
  /**Selects the trait prototypes that have their parameters set.
   * @return Hashtable of the trait prototypes
   **/
  map< trait_t,TraitPrototype* >&  build_currentTraits    ( );
  
  /**Selects the life cycle events that have their parameters set.*/
  void                             build_LifeCycle ( );
  
  /**Accessor to the list of current trait prototypes. 
   * @param type the trait type
   * @return ptr to the prototype
   **/
  TraitPrototype*                  get_current_trait      (trait_t type);
  /**Accessor to the list of current LCEs. 
   * @param name the name of the LCE
   * @return ptr to the LCE
   **/
  LifeCycleEvent*                  get_current_event      (string& name);
  /**Accessor to the list of the selected parameter sets. 
   * @return list of ParamSet
   **/
  list<ParamSet*>&             get_currentParams      ( )     {return _currentParams;}
  
  age_t                        getFirstRequiredAgeInLifeCycle ( );
  
protected:
  /**List of the selected simulation components from the user defined input parameters. **/
  list< ParamSet* >            _currentParams;
  /**List of the selected trait prototypes from the user defined input parameters. **/
  map< trait_t, TraitPrototype* >  _currentTraits;
  /**List of the selected life cycle events from the user defined input parameters. **/
  map< int, LifeCycleEvent* >      _LifeCycle;
  
  typedef map< int, LifeCycleEvent* >::const_iterator LCE_ITER;
  typedef map< trait_t, TraitPrototype* >::const_iterator TRAIT_ITER;
  
};
#endif

