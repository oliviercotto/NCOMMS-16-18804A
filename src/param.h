/** $Id: param.h,v 1.8.2.3 2016-03-10 14:56:47 fred Exp $
 *
 *  @file param.h
 *  Nemo2
 *
 *   Copyright (C) 2006-2015 Frederic Guillaume
 *   frederic.guillaume@ieu.uzh.ch
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
 *  created on @date 27.08.2004
 * 
 *  @author fred
 */

#ifndef PARAM_H
#define PARAM_H

#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <assert.h>
#include <map>
#include <deque>
#include <list>
#include <vector>
#include "handler.h"
#include "tmatrix.h"
#include "types.h"
#include "output.h"

using namespace std;
class SimComponent;
class ParamUpdaterBase;
/**This structure stores one parameter, its definition and its string argument.
 Parameters are aggregated into a ParamSet.
 */
class Param
  {
    /**The name of the parameter as read in the init file.*/
    string          _name;
    /**The argument string, set by the ParamsParser upon initialization.*/
    string          _arg;
    /**The input argument string, kept untouched as a backup for logging.*/
    string          _input_arg;
    /**The type of the argument, must one of DBL, INT, BOOL, STR or MAT (see types.h).*/
    param_t         _type;
    /**Flag set once the parameter has recieved the right argument.*/
    bool            _isSet;
    /**Flag telling if the parameter has a bounded argument value.*/
    bool            _isBounded;
    /**Flag telling if this parameter is mandatory. A SimComponent with a mandatory parameter in the 'unset' state
     will not be selected to be part of a running simulation component.*/
    bool            _isRequired;
    /**The argument value boundaries.*/
    double          _bounds[2];
    /**Generation at which the parameter has been set/updated*/ 
    unsigned int    _setAtGeneration;
    /**Pointer to the component that declared this parameter.*/
    SimComponent*   _myOwner;
    /**Pointer to an ParamUpdater object.*/
    ParamUpdaterBase* _myUpdater;
    
    /**External argument file.*/
    vector< string > _externalFile;
    bool   _hasExternalFile;
    
    /**The temporal arguments.*/
    map< unsigned int, string > _temporalArgs;
    /**Flag telling if this parameter has temporal arguments.*/
    bool _isTemporal;
    
    /**Multiple arguments*/
    vector< string > _multiArgs;
    bool _hasMultipleArgs;
        
  public:
    /**Constructor.
     * @param Name the name of the parameter as read in the init file 
     * @param Type the type of the parameter argument (see types.h), used to convert the argument string into a value
     * @param mandatory specifies whether this parameter is mandatory for the ParamSet owning it
     * @param bounded specifies whether the values this parameter can take are bounded
     * @param low_bnd the lower bound of the parameter value
     * @param up_bnd the upper bound
     * @param owner pointer to the SimComponents that owns this parameter
     * @param updater the param updater object used to update the parameter's state during a simulation
     **/
    Param  (string& Name, param_t Type, bool mandatory, bool bounded, double low_bnd, double up_bnd, 
            SimComponent* owner, ParamUpdaterBase* updater);
    
    /**Copy cstor, provides shallow copy of the parameter definition, do not copy arguments.*/
    Param (const Param& P);
    
    virtual ~Param ();
    
    /**Clears the _isSet flag and the argument string.**/
    void    reset   ();
    
    //accessors:
    /**Sets the parameter's name**/
    void            setName             (string value)          {_name = value;}
    /**Sets the parameter's argument**/
    void            setArg              (string value)          {_arg = value;}
    /**Sets the parameter's type (see types.h)**/
    void            setType             (param_t value)         {_type = value;}
    /**Sets the _isSet flag**/
    void            setIsSet            (bool value)            {_isSet = value;}
    /**Sets the _isBounded flag**/
    void            setIsBounded        (bool value)            {_isBounded = value;}
    /**Sets the bounds**/
    void            setBounds   (double low_bnd,double up_bnd)  {_bounds[0]=low_bnd; _bounds[1]=up_bnd;}
    /**Sets the pointer to owner**/
    void            setOwner            (SimComponent *owner)   {_myOwner = owner;}
    /**Sets the pointer to the updater object*/
    void            setUpdater          (ParamUpdaterBase* updater) {_myUpdater = updater;}
    
    void            setAtGeneration     (unsigned int generation) {_setAtGeneration = generation;}
    
    
    string          getName             ()                      {return _name;}
    string          getArg              ()                      {return _input_arg;}
    param_t         getType             ()                      {return _type;}
    bool            isSet               ()                      {return _isSet;}
    bool            isBounded           ()                      {return _isBounded;}
    bool            isRequired          ()                      {return _isRequired;}
    bool            isTemporal          ()                      {return _isTemporal;}
    bool            hasMultipleArgs     ()                      {return _hasMultipleArgs;}
    bool            hasExternalFile     ()                      {return _hasExternalFile;}
    double          getBound            (unsigned int i)        {return _bounds[i];}
    SimComponent*   getOwner            ()                      {return _myOwner;}
    ParamUpdaterBase* getUpdater        ()                      {return _myUpdater;}
    
    deque< unsigned int >  getUpdatingDates ();  
    
    deque< string > getTemporalArgs();
    
    vector< string > getMultiArgs();
    
    vector< string > getExternalFiles() {return _externalFile;}

    /**Sets the _isSet flag to true and _arg to arg if the arg is of the right type and whithin the bounds.
     Called at simulation setup (i.e. generation 0). For parameter update during the simulation, check the
     Param::update member function.**/
    bool            set                  (string arg, string& errmsg);
    
    /**Updates the parameter value at a given generation during the simulation.*/
    bool            update               (unsigned int generation);
    
    /**Returns the argument value according to its type.
     *@return -1 if the parameter is not set or not of a the right type
     **/
    double          getValue            ();
    
    /**Checks if the argument is of matrix type.**/
    bool            isMatrix            () 
    {  return (_type == MAT || (_arg.size() != 0 ? _arg[0] == '{' : false) ); }
    
    /**Sets the matrix from the argument string if the parameter is set and of matrix type. Passes the param 'mat' to
     *function Param::parse_matrix.
     *@param mat a TMatrix ptr, mat dimensions and values will be reset to the values read in the init file.
     **/
    void            getMatrix           (TMatrix* mat);
    
    void            getVariableMatrix   (vector< vector <double> >* mat);
   
    /**Parses the matrix from the argument string.
     *@param mat a TMatrix ptr, mat dimensions and values will be reset to the values read in the init file.
     **/
    void            parse_matrix         (TMatrix* mat);
    void            parse_variable_matrix         (vector< vector <double> >* mat);
    bool            parseArgument        (string& arg);
    bool           parseTemporalArgument (const string& arg);
    bool           parseAgeSpecArgument  (const string& arg);
    bool           parseSubParamArgument (const string& arg);
    string         checkArgumentForExpansion (string arg);
    string         getArgumentFromFile   (string file);
    
    /**Print state to stdout.**/
    void            show_up              ();
    
  };


/**Parameters container, implemented in each SimComponent. 
 A SimComponent is added to the set of active components of a simulation only if all
 its required parameters are set (isSet = true).
 */
class ParamSet
  {
    
    string                  _name;
    bool                    _isSet;
    bool                    _isRequired;
    map<string, Param*>     _params;
    /**Pointer to the component that declared this parameter.*/
    SimComponent*           _myOwner;
    
  public:
    
    ParamSet                             ( ) : _isSet(0), _isRequired(0) { }
    ParamSet                             (const ParamSet& PS);
    ~ParamSet                            ( );
    /**Put the container in the unset state, reset each Param it contains.**/
    void            reset                ( );
    
    /**Empties the parameter containers (no delete).*/
    void            clear                ( ) {_params.clear();}
    /**Checks for the status of the required parameters.
     *@return TRUE if all required parameters are or if nothing is set and the container is not required
     **/
    bool            check_consistency    ( );
    /**print info to stdout.**/
    void            show_up              ( );
    /**print all set parameters to the outpout file stream**/
    void            print                (ofstream& FILE);
    /**Returns the number of parameters contained**/
    int             size                 ( )                     {return _params.size();}
    /**Returns the complete list of parameters**/
    map<string, Param*>& getAllParams    ( )                     {return _params;}
    ///@name Accessors to Param members.
    ///@{
    /**Adds the param argument to the list**/
    void            add_param            (Param* param)         {_params[param->getName()] = param;}
    /**Adds a new param specified by arguments to the list.
     *@param Name the name of the parameter
     *@param Type the type of the parameter
     *@param mandatory specifies if this parameter is required and must be set for the container to gain the "set" status
     *@param isBounded specified whether this parameter is bounded
     *@param low_bnd the lower value the parameter can take, used if isBounded is true
     *@param up_bnd the upper value the parameter can take, used if isBounded is true
     **/
    void            add_param            (string Name, param_t Type, bool mandatory, bool isBounded,
                                          double low_bnd, double up_bnd)
    {add_param(Name, Type, mandatory, isBounded, low_bnd, up_bnd, 0);}
    
    void            add_param            (string Name, param_t Type, bool mandatory, bool isBounded,
                                          double low_bnd, double up_bnd, ParamUpdaterBase* updater);
    /**Look for a param named "Name" and try to set it with the "Arg" argument string.
     *@return TRUE if param Name has been found and set with Arg
     *@return FALSE otherwise
     *@param Name the name of the parameter to find in the list
     *@param Arg the argument string as found in the init params
     **/
    bool            set_param            (string Name, string Arg);
    /**Look for a param "name" in its parameters list.
     *@return NULL if no Param with _name = name exists
     **/
    Param*          get_param            (string name);
    /**Calls the updating procedure of each of its Param.*/
    bool            update_param         (string Name, unsigned int generation);
    /**Sets the container's name**/
    void            setName              (string value)          {_name = value;}
    /**Sets the _isRequired flag meaning this container is mandatory and must be set in order to run a simulation**/
    void            setIsRequired        (bool value)            {_isRequired = value;}
    /**Sets the pointer to the SimComponents that owns this set.*/
    void            setOwner             (SimComponent* owner)   {_myOwner = owner;}
    /**Accessor to the status flag. **/
    bool            isSet                ()                      {return _isSet;}
    /**Name accessor.**/
    string          getName              ()                      {return _name;}
    /**Accessor to the mandatory flag. **/
    bool            isRequired           ()                      {return _isRequired;}
    /**Accessor to the parameters status flag. **/
    bool            isSet                (string name)           {return (get_param(name))->isSet();}
    /**Check if the parameter "name" is of matrix type. **/
    bool            isMatrix             (string name)           {return (get_param(name))->isMatrix();}
    /**Check if the parameter "name" has temporal arguments. **/
    bool            isTemporal           (string name)           {return (get_param(name))->isTemporal();}
    /**Accessor to the parameters argument string.**/
    string          getArg               (string name)           {return (get_param(name))->getArg();}
    /**Accessor the parameters value.**/
    double          getValue             (string name)           {return (get_param(name))->getValue();}
    /**Accessor to the parameters matrix.**/
    void            getMatrix   (string name, TMatrix* mat)      {return (get_param(name))->getMatrix(mat);}
    /**Collects the parameter updaters from the set of parameters*/
    list<ParamUpdaterBase*> getUpdaters();
    ///@}
    ParamSet&       operator= (const ParamSet& PS);
  };


//class ParamUpdaterBase
//
/**Base class of the ParamUpdater class used to handle the temporal parameter argument values*/
class ParamUpdaterBase : public Handler {
  
protected:
  /**List of the parameters affected by this updater.*/
  list< Param* > _params;
  
public:
  
  typedef list< Param* >::iterator PIT;
  
public:
  
  ParamUpdaterBase() { }
  
  ParamUpdaterBase (const ParamUpdaterBase& PU)
  {
    _params.assign(PU._params.begin(), PU._params.end());
  }
  
  virtual ~ParamUpdaterBase() {}
  
  virtual void init   () = 0;
  virtual void update () {}
  /**Updating procedure.*/
  virtual bool update (unsigned int generation) = 0;
  
  virtual SimComponent* getComponent () = 0;
  /**Adds a parameter to the stack.*/
  virtual void addParam        (Param* param)       {_params.push_back( param );}
  /**Clears the parameters stack.*/
  virtual void reset           ( )                  {_params.clear();}
  /**Returns the list of parameters.*/
  list< Param* > getParams     ( )                  {return _params;}
  
};


//class ParamUpdater
//
/**Implementation of the ParamUpdaterBase interface. Stores the pointers to the SimComponent and its
   updating function. The template argument allows to pass SimComponent class type necessary to
   call the updating function.
   The updating procedure is as follow: the parameters are first updated with their generational 
   values (if the generation matches) and the component's is then updated using its updating function.
   The updating process is implemented in the UpdaterServices class.
   Several parameters may want to use the same updating function and thus this class allows to run
   that function only once, after updating all the parameters that use it. */
   
template <class SC> class ParamUpdater : public virtual ParamUpdaterBase {
  
  bool (SC::*_myUpdaterFunction) (void);
  SC* _myComponent;
  
public:
  
  ParamUpdater(bool (SC::*updateFuncPtr) (void) )
  {
    _myUpdaterFunction = updateFuncPtr;
  }
  
  ParamUpdater (const ParamUpdater< SC >& PU)
  {
    _params.assign(PU._params.begin(), PU._params.end());
    _myUpdaterFunction = PU._myUpdaterFunction;
    _myComponent = PU._myComponent;
  }
  
  virtual ~ParamUpdater ()
  {
    for(PIT pit = _params.begin(); pit != _params.end(); pit++)
      if((*pit)->getUpdater() == this) (*pit)->setUpdater(0);
  }
  
  /**Sets the pointer to the SimComponent. It is deduced from the parameters' owner.*/
  virtual void init ()
  {
    _myComponent = dynamic_cast<SC*> ( (*_params.begin())->getOwner() );
    
    for(PIT pit = ++_params.begin(); pit != _params.end(); pit++)
      assert(_myComponent == (*pit)->getOwner());
  }
  /**Calls the SimComponent's updating function using its pointer.*/
  virtual bool update (unsigned int generation)
  {        
    return ((( dynamic_cast<SC*> (_myComponent))->*_myUpdaterFunction) ());
  }
  /**Sets the pointer to the updating function.*/
  void setFuncPtr   ( bool (SC::*updateFuncPtr) (void))
  {
    _myUpdaterFunction = updateFuncPtr;
  }
  /**Accessor to the SimCimponent.*/
  virtual SC* getComponent () {return _myComponent;}
};

#endif
