/** $Id: simcomponent.h,v 1.5.2.1 2014-04-29 18:26:04 fred Exp $ 
*
*  @file simcomponent.h
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
*  created on @date 16.06.2005
* 
*  @author fred
*/

#ifndef SIMCOMPONENT_H
#define SIMCOMPONENT_H

#include <string>
#include "fileservices.h"
#include "statservices.h"
#include "updaterservices.h"
#include "param.h"
#include "binarystoragebuffer.h"

/**Interface to all basic components of a simulation (traits, life cycle events, pop, etc.\ ). 
 * Implements the interface to handle the files and stats handler, the parameter updaters,
 * and the parameters container.
 * Contains the component's parameters and declares the interface to interact with its parameters.
 */
class SimComponent {
protected:
  /**The parameters container. */
  ParamSet* _paramSet;
  
public:
  
  SimComponent ( ) : _paramSet(0) {}//{ cout << "calling SimComponent()" << endl;}
  //**Dstor. Deletes the parameter container. */
  virtual ~SimComponent ( ) {if(_paramSet != NULL) delete _paramSet;}
  
  ///@name Services interface
  ///@{
  /**Loads the component's FileHandler onto the FileServices. 
    @param loader the file service */
  virtual void loadFileServices ( FileServices* loader ) = 0;
  
  /**Loads the component's StatHandler onto the StatServices. 
     @param loader the stat service */
  virtual void loadStatServices ( StatServices* loader ) = 0;
  
  /**Loads the parameters and component updater onto the updater manager.
     @param loader the updater manager */
  virtual void loadUpdaters  ( UpdaterServices* loader ) {
    list<ParamUpdaterBase*> updaters = _paramSet->getUpdaters();
    //remove duplicates:
    updaters.unique();
    for (list<ParamUpdaterBase*>::iterator u = updaters.begin(); u != updaters.end(); u++) {
      loader->attach( (*u) );
    }
  }
  ///@}

  ///@name Parameters handling
  ///@{  
  /** Default interface needed to initialize the component's variables from its input 
      parameters value. Formerly called 'init'.*/
  virtual bool setParameters () = 0;
  
  /**Sets the ParamSet member.
     @param paramset a pointer to a ParamSet object */
  virtual void set_paramset(ParamSet * paramset)
  { _paramSet = paramset; }
  
  /**Sets a new ParamSet and name it. 
    @param name the name of the parameters container
    @param required tag whether the component is required to run a simulation
    @param owner a reference to the owner of the ParamSet (should be this).**/
  virtual void set_paramset (std::string name, bool required, SimComponent* owner ) {
    if(_paramSet != NULL) delete _paramSet;
    _paramSet = new ParamSet(); 
    _paramSet->setName(name);
    _paramSet->setIsRequired(required);
    _paramSet->setOwner(owner);
  }
  
  /**Reset the set of parameters from a another set.
    @param PSet the parameter set to copy parameters from */
  virtual void set_paramsetFromCopy (const ParamSet &PSet) {
    if(_paramSet != NULL) delete _paramSet;
    _paramSet = new ParamSet(PSet); 
  }
  /**ParamSet accessor. **/
  virtual ParamSet*  get_paramset () {return _paramSet;}
  
  /**Interface to add a parameter to the set.
   * @param param The parameter to add to the set */
  virtual void add_parameter (Param* param) {_paramSet->add_param(param);}
  
  /**Interface to add a parameter to the set.
   * @param Name The param string name as read in the init file
   * @param Type The type of the argument (DBL=double, INT=integer, BOOL=boolean, STR=string, etc., see type.h)
   * @param isRequired True if the parameter is mandatory
   * @param isBounded True if the parameter takes bounded values as argument
   * @param low_bnd The lower value the argument can take
   * @param up_bnd The upper value the argument can take */
  virtual void add_parameter (std::string Name, param_t Type, bool isRequired, bool isBounded,
                              double low_bnd, double up_bnd)
  {_paramSet->add_param(Name, Type, isRequired, isBounded, low_bnd, up_bnd, 0);}
    
  /**Interface to add a parameter and its updater to the set.
   * @param Name The param string name as read in the init file
   * @param Type The type of the argument (DBL=double, INT=integer, BOOL=boolean, STR=string, etc., see type.h)
   * @param isRequired True if the parameter is mandatory
   * @param isBounded True if the parameter takes bounded values as argument
   * @param low_bnd The lower value the argument can take
   * @param up_bnd The upper value the argument can take
   * @param updater a pointer to a ParamUpdaterBase object used to reset its value during a simulation */
  virtual void add_parameter (std::string Name, param_t Type, bool isRequired, bool isBounded,
                              double low_bnd, double up_bnd, ParamUpdaterBase* updater)
  {_paramSet->add_param(Name, Type, isRequired, isBounded, low_bnd, up_bnd, updater);}
  
  /**Param getter. Asks the ParamSet to return a pointer to Param \a name.
    @see ParamSet::get_param */
  virtual Param*    get_parameter       (std::string name)  {return _paramSet->get_param(name);}
  
  /**Param value getter. Asks the ParamSet to return the value of Param \a name.
    @see ParamSet::getValue */
  virtual double    get_parameter_value (std::string name) {return _paramSet->getValue(name);}
  
  /**Returnd the name of the ParamSet, i.e. the name of the component.*/
  virtual string    get_name () {return _paramSet->getName();}
  ///}
  
};

/**Provides an interface to binary data saving and uploading.
   This interface is used to save and load data from binary data files (serialization).
   @see BinaryDataSaver and BinaryDataLoader
*/
class StorableComponent {
 
public:
  /**Interface to store the component data (e.g.\ gene values) into a binary buffer. **/
  virtual void store_data    ( BinaryStorageBuffer* saver  ) = 0;
  /**Interface to retrieve the same data from the binary buffer. **/
  virtual bool retrieve_data ( BinaryStorageBuffer* reader ) = 0;
  
  virtual ~StorableComponent ( ) { }
  
};

#endif  //SIMCOMPONENT_H

