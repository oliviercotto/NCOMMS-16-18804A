/** $Id: param.cc,v 1.8.2.5 2016-03-10 14:56:47 fred Exp $
 *
 *  @file param.cc
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

//#include <stdlib.h>
//#include <iostream>
#include <vector>
#include <sstream>
#include "paramsparser.h"
#include "param.h"
#include "output.h"
#include "tstring.h"

// ----------------------------------------------------------------------------------------
// Param
// ----------------------------------------------------------------------------------------
Param::Param (string& Name,param_t Type,bool mandatory,bool bounded,double low_bnd,double up_bnd,
              SimComponent* owner, ParamUpdaterBase* updater)
: _name(Name),_arg(""),_type(Type),_isSet(0),_isBounded(bounded),_isRequired(mandatory),
  _setAtGeneration(0),_myOwner(owner), _myUpdater(updater),_hasExternalFile(0),_isTemporal(0),
  _hasMultipleArgs(0)

{
  _bounds[0] = low_bnd;
  _bounds[1] = up_bnd;
}
// ----------------------------------------------------------------------------------------
// Param
// ----------------------------------------------------------------------------------------
Param::Param (const Param& P)
: _name(P._name),_arg(P._arg),_type(P._type),_isSet(P._isSet),_isBounded(P._isBounded),
_isRequired(P._isRequired),_setAtGeneration(P._setAtGeneration),_myOwner(P._myOwner),
_myUpdater(0),_hasExternalFile(P._hasExternalFile),_hasMultipleArgs(P._hasMultipleArgs)
{
  _bounds[0] = P._bounds[0];
  _bounds[1] = P._bounds[1];
}
// ----------------------------------------------------------------------------------------
// ~Param
// ----------------------------------------------------------------------------------------
Param::~Param ( ) 
{
  if(_myUpdater) delete _myUpdater;
}
// ----------------------------------------------------------------------------------------
// reset
// ----------------------------------------------------------------------------------------
void Param::reset ()
{
  _arg = "";
  _isSet = 0;
  _setAtGeneration = 0;
  _multiArgs.clear();
  _externalFile.clear();
}
// ----------------------------------------------------------------------------------------
// getUpdatingDate
// ----------------------------------------------------------------------------------------
deque< unsigned int >  Param::getUpdatingDates ()
{
  deque< unsigned int > dates;
  
  for(map< unsigned int, string >::iterator tmp_it = _temporalArgs.begin();
      tmp_it != _temporalArgs.end(); tmp_it++)
    dates.push_back( tmp_it->first );
  
  return dates;
}
// ----------------------------------------------------------------------------------------
// getTemporalArgs
// ----------------------------------------------------------------------------------------
deque< string > Param::getTemporalArgs()
{
  deque< string > args;
  
  for(map< unsigned int, string >::iterator tmp_it = _temporalArgs.begin();
      tmp_it != _temporalArgs.end(); tmp_it++)
    args.push_back( tmp_it->second );
  
  return args;
}
// ----------------------------------------------------------------------------------------
// getTemporalArgs
// ----------------------------------------------------------------------------------------
vector< string > Param::getMultiArgs()
{
  if (!_hasMultipleArgs) {
    _multiArgs.clear();
    _multiArgs.push_back(_arg);
  }
  return _multiArgs;
}
// ----------------------------------------------------------------------------------------
// set
// ----------------------------------------------------------------------------------------
bool Param::set (string arg, string& errmsg)
{  
//  cout << "Param::set ("<<arg<<")\n";

  _input_arg = arg; //backup the input parameter argument string

  _externalFile.clear(); //ext. file names will be pushed here
  //check if argument is in an external file, and read it from that file
  if (arg[0] == '&') {
      //read the argument from an external file, file name is stored within the call.
      string expanded = checkArgumentForExpansion(arg);
      arg = expanded;
    }

  if( arg[0] == '(' ) {
    if( !parseArgument( arg ) ) {
      errmsg = "could not parse the argument string";
      return false;
    }
  }

  //integer or decimal: 
  if(_type == INT || _type == DBL) {
    //!!a matrix may also be specified!! in that case, we don't check values here
    if(arg[0] != '{') {
      
      if(!tstring::isanumber(arg)) {
        errmsg = "argument is not a number";
        return false;
      }
      
      double val = tstring::str2dble(arg);
      
      if(_isBounded && (val < _bounds[0] || val > _bounds[1]) ){
        errmsg = "argument is out of bounds";
        return false;
      }
    }
    
    //boolean:
  } else if(_type == BOOL) {
    
    if( !tstring::isanumber(arg) ){
      errmsg = "argument is not a boolean";
      return false;
    }
    
  } else if(_type == MAT) {
    if(arg[0] != '{'){
      errmsg = "argument is not a matrix";
      return false;
    }
  }
  // else if(_type == P || _type == S): no conditions to check
  
  _isSet = (_type == BOOL && tstring::str2int(arg) == 0 ? false : true);

  //now we set the value of the parameter to the value of the argument read in the file
  _arg = arg;

  return true;
}
// ----------------------------------------------------------------------------------------
// parseArgument
// ----------------------------------------------------------------------------------------
bool Param::parseArgument (string& arg)
{/**Parse an argument string enclosed in '( )', might hold temporal or multiple values.*/

  arg = tstring::removeEnclosingChar(arg, '(', ')');
  
  vector< string > args = tstring::splitExcludeEnclosedDelimiters(arg);
  
  _temporalArgs.clear();
  _multiArgs.clear();
  
//  cout<< "Param::parseArgument: " << endl;

  
  for(unsigned int i = 0; i < args.size(); ++i) {
//    cout << args[i] << endl;
    if(args[i][0] == '@') {
      if(args[i][1] == 'g') {if(!parseTemporalArgument(args[i])) return false;}
      else if(args[i][1] == 'a') {if(!parseAgeSpecArgument(args[i])) return false;}
      else {
        error("parameter \"%s\" has unrecognized specifier in its argument string \"%s\".\n", _name.c_str(), args[i].c_str());
        return false;
      }
    } else { // (arg1, arg2, arg3, ...)
      //values might still be in external files
      _multiArgs.push_back( checkArgumentForExpansion(args[i]) );
    }
  }
  
  if(_multiArgs.size() != 0) _hasMultipleArgs = true;
  else _hasMultipleArgs = false;
  
  if(_temporalArgs.size() != 0) {
    
    //cout << "temporal arg values:\n";
    //    for(map< unsigned int , string >::iterator iter = _temporalArgs.begin();
    //        iter != _temporalArgs.end();
    //        iter++)
    //        cout << iter->first <<" "<< iter->second << endl;
    
    if(_myUpdater != 0) {
      
      if(_temporalArgs.find(0) != _temporalArgs.end()) {
        _isTemporal = true;
        //set the value for the first generation
        arg = _temporalArgs.find(0)->second;
      }
      else {
        error("first generation argument value for temporal parameter \"%s\" is missing; no \"@g0\".\n", _name.c_str());
        return false;
      }
      
    } else {
      warning("trying to pass a temporal argument to a non temporal parameter (\"%s\").\n", _name.c_str());
      _isTemporal = false;
      _temporalArgs.clear();
    }
    
  } else _isTemporal = false;
  
  
  return true;
}
// ----------------------------------------------------------------------------------------
// parseTemporalArgument
// ----------------------------------------------------------------------------------------
bool Param::parseTemporalArgument (const string& arg)
{
  vector< string > args = tstring::split(arg, ' ', true);
  
  if(args.size() != 2) {
    error("Param::parseTemporalArgument:: missing argument value in \"%s %s\".\n", _name.c_str(), arg.c_str());
    return false;
  }
  
  unsigned int gen = tstring::str2uint(tstring::removeFirstCharOf(tstring::removeFirstCharOf(args[0], '@'), 'g'));
  
  args[1] = checkArgumentForExpansion(args[1]);
  
  _temporalArgs[gen] = args[1]; //still need to process the argument string for replacements!!
  
  return true;
}
// ----------------------------------------------------------------------------------------
// parseArgument
// ----------------------------------------------------------------------------------------
bool Param::parseAgeSpecArgument (const string& arg)
{
  return true;
}
// ----------------------------------------------------------------------------------------
// parseSubParamArgument
// ----------------------------------------------------------------------------------------
bool Param::parseSubParamArgument (const string& arg)
{
  return true;
}
// ----------------------------------------------------------------------------------------
// checkArgumentForExpansion
// ----------------------------------------------------------------------------------------
string Param::checkArgumentForExpansion (string arg)
{
  string expanded;
  
  if (arg[0] == '&') { //this is the char indicating an external argument file
    expanded = getArgumentFromFile(arg);
  } else
    expanded = arg;
  
  return expanded;
}
// ----------------------------------------------------------------------------------------
// getArgumentFromFile
// ----------------------------------------------------------------------------------------
string Param::getArgumentFromFile (string file)
{
  StreamParser Parser(tstring::removeFirstCharOf(file, '&').c_str());
  int dummyCnt;
  string arg;
  
  ifstream EXT(tstring::removeFirstCharOf(file, '&').c_str());
  if(!EXT) fatal("External parameter file '%s' could not be found!\n", file.c_str());
  while(Parser.readArguments(EXT, dummyCnt, arg));
  EXT.close();
  
  //add the external file name to the list
  _externalFile.push_back(file);
  _hasExternalFile = true;

  return arg;
}
// ----------------------------------------------------------------------------------------
// getValue
// ----------------------------------------------------------------------------------------
double Param::getValue ()
{
  if( !(isMatrix() || _type == STR) && _isSet)
    return atof(_arg.c_str());
  else
    return -1.0;
}
// ----------------------------------------------------------------------------------------
// getMatrix
// ----------------------------------------------------------------------------------------
void Param::getMatrix (TMatrix* mat)
{
  if( isMatrix() && _isSet ){
    
    parse_matrix(mat);
    
  } else
    warning("param \"%s\" is not a matrix!\n",_name.c_str());
}

// ----------------------------------------------------------------------------------------
// parse_matrix
// ----------------------------------------------------------------------------------------
void Param::parse_matrix (TMatrix* mat)
{
  std::vector< std::vector<double> > tmpMat;
  std::istringstream IN;
  
  IN.str(_arg);
  
  unsigned int cols = 0;
  double elmnt;
  char c;
  
  
  int rows = -1, pos = -1;
  //we count the number of rows
  do {
    pos = _arg.find("{", pos + 1);
    //message("pos %i, arg %s\n",pos,_arg.c_str());
    rows++;
  }while(pos != (int)string::npos);
  //decrement by 1, the first doesn't count for a row
  rows--;
  
  for(int i = 0; i < rows; i++)
    tmpMat.push_back( vector<double>()  );
  
  //remove the first enclosing bracket
  IN>>c;
  //then read the rows
  for(unsigned int i = 0; i < tmpMat.size(); i++) {
    
    cols = 0;
    
    //read a row enclosed by {...}:
    while(IN) {
      
      //first character:
      IN>>c;
      
      if(c == '{' || c == ',' || c == ';') {
        //read a row element:
        IN>>elmnt;
        cols++;
        tmpMat[i].push_back(elmnt);
        
      } else if(c == '}')
        //go to next row
        break;
      
    }
  }
  //check for matrix coherence:
  for(unsigned int i = 0; i < tmpMat.size(); i++) {
    if(tmpMat[i].size() != cols)
      fatal("%s: not same number of elements in all rows of matrix! (%i, cols %i)\n",_name.c_str(),tmpMat[i].size(),cols);
    cols = tmpMat[i].size();
  }
  
  //copy to input TMatrix:
  mat->reset(rows, cols);
  for(int i = 0; i < rows; ++i)
    for(unsigned int j = 0; j < cols; ++j)
      mat->set(i,j,tmpMat[i][j]);
  
//  for(int i = 0; i < rows; i++) {
//    delete &(tmpMat[i]);
//  }
  tmpMat.clear();
}
// ----------------------------------------------------------------------------------------
// getVariableMatrix
// ----------------------------------------------------------------------------------------
void Param::getVariableMatrix (vector< vector< double > >* mat)
{
  if( isMatrix() && _isSet ){
    
    parse_variable_matrix(mat);
    
  } else
    warning("param \"%s\" is not a (variable) matrix!\n",_name.c_str());
}
// ----------------------------------------------------------------------------------------
// parse_variable_matrix
// ----------------------------------------------------------------------------------------
void Param::parse_variable_matrix (vector< vector< double > >* mat)
{
  //  std::vector< std::vector<double> > tmpMat;
  std::istringstream IN;
  
  IN.str(_arg);
  
  double elmnt;
  char c;
  
  //purge the input matrix:
  mat->clear();
  
  int rows = -1, pos = -1;
  //we count the number of rows
  do {
    pos = _arg.find("{", pos + 1);
    //message("pos %i, arg %s\n",pos,_arg.c_str());
    rows++;
  }while(pos != (int)string::npos);
  //decrement by 1, the first doesn't count for a row
  rows--;
  
  for(int i = 0; i < rows; i++)
    mat->push_back( vector<double>()  );
  
  //remove the first enclosing bracket
  IN>>c;
  //then read the rows
  for(unsigned int i = 0; i < mat->size(); i++) {
    
    //read a row enclosed by {...}:
    while(IN) {
      
      //first character:
      IN>>c;
      
      if(c == '{' || c == ',' || c == ';') {
        //read a row element:
        IN>>elmnt;
        (*mat)[i].push_back(elmnt);
        
      } else if(c == '}')
        //go to next row
        break;
    }
  }
}
// ----------------------------------------------------------------------------------------
// update
// ----------------------------------------------------------------------------------------
bool Param::update (unsigned int generation)
{
  string error_msg;
  if(_temporalArgs.find( generation ) != _temporalArgs.end()) {
    if( !set( _temporalArgs[generation], error_msg ) ){
      error("could not set \"%s\": %s", _name.c_str(), error_msg.c_str());
      return false;
    } else
      return true;
  } else
    error("parameter \"%s\" does not have an argument value for generation %i.\n", _name.c_str(), generation);
  return false; 
}
// ----------------------------------------------------------------------------------------
// show_up
// ----------------------------------------------------------------------------------------
void Param::show_up ()
{
  message("%s = %s",_name.c_str(),_arg.c_str());
  if(_isSet)
    message(" (set at gen %i)\n",_setAtGeneration);
  else
    message(" (not set)\n");
}
// ------------------------------------------------------------------------------

//                             ParamSet

// ----------------------------------------------------------------------------------------
// ParamSet
// ----------------------------------------------------------------------------------------
ParamSet::ParamSet (const ParamSet& PS) 
: _name(PS._name), _isSet(PS._isSet), _isRequired(PS._isRequired)
{
  map<string, Param*>::iterator param;
  
  if( _params.size() ) {
    param =  _params.begin();
    while(param != _params.end()) {
      delete param->second;
      param++;
    }
    _params.clear();
  }
  map<string, Param*> PS_params = PS._params;
  param =  PS_params.begin();
  while(param != PS_params.end()) {
    _params[param->first] = new Param( *param->second );
    param++;
  }
  
}

// ----------------------------------------------------------------------------------------
// ~ParamSet
// ----------------------------------------------------------------------------------------
ParamSet::~ParamSet ()
{
  map<string, Param*>::iterator param =  _params.begin();
  while(param != _params.end()) {
    delete param->second;
    param++;
  }
}
// ----------------------------------------------------------------------------------------
// reset
// ----------------------------------------------------------------------------------------
void ParamSet::reset ()
{
  _isSet = 0;
  map<string, Param*>::iterator param =  _params.begin();
  while(param != _params.end()) {
    param->second->reset();
    param++;
  }
}
// ----------------------------------------------------------------------------------------
// add_param
// ----------------------------------------------------------------------------------------
void ParamSet::add_param (string Name, param_t Type, bool mandatory, bool isBounded,
                          double low_bnd, double up_bnd, ParamUpdaterBase* updater)
{
  //  if(!_name.empty()) Name = _name + "_" + Name; //not yet...
  Param* param = new Param(Name, Type, mandatory, isBounded, low_bnd, up_bnd, _myOwner, updater);
  if(updater) updater->addParam(param);
  _params[Name] = param;
}
// ----------------------------------------------------------------------------------------
// getUpdaters
// ----------------------------------------------------------------------------------------
list<ParamUpdaterBase*> ParamSet::getUpdaters()
{
  list<ParamUpdaterBase*> updaters;
  map<string, Param*>::iterator param = _params.begin();
  for(;param != _params.end(); param++) {
    if(param->second->getUpdater() != 0)
      updaters.push_back(param->second->getUpdater());
  }
  return updaters;
}
// ----------------------------------------------------------------------------------------
// set_param
// ----------------------------------------------------------------------------------------
bool ParamSet::set_param (string Name, string Arg)
{
  map<string, Param*>::iterator param = _params.find(Name);
  string error_msg;
  
  if(param == _params.end()) {
    //    error("could not set \"%s\": parameter not found.\n",Name.c_str());
    //exist silently
    return false;
  }
  
  if( !param->second->set(Arg, error_msg) ) {
    error("could not set \"%s\": %s.\n", Name.c_str(), error_msg.c_str());
    return false;
  } else
    return true;
}
// ----------------------------------------------------------------------------------------
// get_param
// ----------------------------------------------------------------------------------------
Param* ParamSet::get_param (string Name)
{
  map<string, Param*>::iterator param = _params.find(Name);
  
  if(param != _params.end())
    return param->second;
  else
    fatal("parameter \"%s\" is not a member of \"%s\".\n", Name.c_str(), _name.c_str());
  //return NULL;
  return NULL;
}
// ----------------------------------------------------------------------------------------
// update_param
// ----------------------------------------------------------------------------------------
bool ParamSet::update_param (string Name, unsigned int generation)
{
  Param* param = get_param(Name);
  if(param) return param->update( generation );
  return false;
}
// ----------------------------------------------------------------------------------------
// check_consistency
// ----------------------------------------------------------------------------------------
bool ParamSet::check_consistency ()
{
  map<string, Param*>::iterator param = _params.begin();
  bool isOK = true;
  bool touched = false;
  
  while(param != _params.end()) {
    //check if all required field have been set properly
    if(param->second->isRequired())
      isOK &= param->second->isSet();
    //else we don't care..
    //check if at least one param has been set
    touched |= param->second->isSet();
    param++;
  }
  _isSet = isOK;
  //return isOk or check if _isRequired in case no params are set (untouched params)
  return ( isOK | (!_isRequired));// & !touched) );
}
// ----------------------------------------------------------------------------------------
// operator=
// ----------------------------------------------------------------------------------------
ParamSet& ParamSet::operator= (const ParamSet& PS) 
{
  if(this != &PS) {
    _name = PS._name;
    _isSet = PS._isSet;
    _isRequired = PS._isRequired;
    
    map<string, Param*>::iterator param;
    
    if( _params.size() ) {
      param =  _params.begin();
      while(param != _params.end()) {
        delete param->second;
        param++;
      }
      _params.clear();
    }
    map<string, Param*> PS_params = PS._params;
    param =  PS_params.begin();
    while(param != PS_params.end()) {
      _params[param->first] = new Param( *param->second );
      param++;
    }
  }
  return *this;
}
// ----------------------------------------------------------------------------------------
// show_up
// ----------------------------------------------------------------------------------------
void ParamSet::show_up ()
{
  message("%s\n",_name.c_str());
  map<string, Param*>::iterator param = _params.begin();
  while(param != _params.end()) {
    param->second->show_up();
    param++;
  }
}
// ----------------------------------------------------------------------------------------
// print
// ----------------------------------------------------------------------------------------
void ParamSet::print (ofstream& FILE)
{
  map<string, Param*>::iterator paramRec = _params.begin();
  Param* param;
  
  while(paramRec != _params.end()) {
    
    param = paramRec->second;
    
    if(param->isSet()){
      //    FILE.width(20);
      //    FILE.setf(ios::left,ios::adjustfield);
      
      FILE<<param->getName()<<" ";
      
//      if (param->isTemporal()) {
//
//        deque< unsigned int > gen = param->getUpdatingDates();
//        deque< string > args = param->getTemporalArgs();
//
//        FILE<<"(";
//
//        for (unsigned int i = 0, size = gen.size(); i < size; i++) {
//
//          FILE<<"@g"<<gen[i]<<" "<<args[i];
//
//          if (i < size-1) {
//            FILE << ", ";
//          }
//          else FILE << ")";
//        }
//
//      } else if (param->hasExternalFile()) {
//
//    	vector<string> ext_files = param->getExternalFiles();
//
//    	if(ext_files.size() == 1) {
//
//          FILE<<ext_files[0];
//
//    	} else {
//
//          FILE<<"(";
//          for (unsigned int i = 0; i < ext_files.size() -1; ++i) {
//            FILE<<ext_files[i]<<", ";
//          }
//          FILE<<ext_files[ ext_files.size() -1 ]<<")";
//    	}
//
//      } else
          FILE<<param->getArg();
      
      FILE<<endl;
    }
    
    paramRec++;
  }
}
