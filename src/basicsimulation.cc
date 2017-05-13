/** $Id: basicsimulation.cc,v 1.13.2.3 2016-04-28 13:01:32 fred Exp $
*
*  @file basicsimulation.cc
*  Nemo2
*
*  Copyright (C) 2006-2015 Frederic Guillaume
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
*  created on @date 21.11.2005
* 
*  @author fred
*/
#include <iostream>
#include <sstream>
#include <cmath>
#include "basicsimulation.h"
#include "output.h"

extern MPIenv *_myenv;

//----------------------------------------------------------------------------------------
// build_component_list
// ----------------------------------------------------------------------------------------
void ComponentManager::build_component_list ()
{
  _components.clear();
  
  list< TraitPrototype* >::iterator TT_it = _TTrait_Templates.begin();
  
  while(TT_it != _TTrait_Templates.end()) {
    _components.push_back( (*TT_it) );
    TT_it++;
  }
  
  list< LifeCycleEvent* >::iterator LCE_it = _LCE_Templates.begin();
  
  while(LCE_it != _LCE_Templates.end()) {
    _components.push_back( (*LCE_it) );
    LCE_it++;
  }
}
//----------------------------------------------------------------------------------------
// get_trait
// ----------------------------------------------------------------------------------------
TraitPrototype* ComponentManager::get_trait (string name)
{
  list< TraitPrototype* >::iterator TT = _TTrait_Templates.begin();
  
  while(TT != _TTrait_Templates.end()) {
	
	if( (*TT)->get_paramset()->getName().compare(name) == 0 ) return (*TT);
	
	TT++;
  }
  
  return NULL;
}
//----------------------------------------------------------------------------------------
// get_LCE
// ----------------------------------------------------------------------------------------
LifeCycleEvent* ComponentManager::get_LCE (string name)
{
  list< LifeCycleEvent* >::iterator LCE = _LCE_Templates.begin();
  
  while(LCE != _LCE_Templates.end()) {
    
    if( (*LCE)->get_event_name().compare(name) == 0 ) return (*LCE);
    
    LCE++;
  }
  
  return NULL;
}
/* /_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/ */

//                             ******ParamManager*******

//----------------------------------------------------------------------------------------
// cstor
// ---------------------------------------------------------------------------------------
ParamManager::ParamManager()
{
  _paramSet.setName("simulation");
  _paramSet.setIsRequired(true);
  _paramSet.setOwner(NULL);
  _paramSet.add_param("filename",STR,true,false,0,0);
  _paramSet.add_param("root_dir",STR,false,false,0,0);
  _paramSet.add_param("logfile",STR,false,false,0,0);
  _paramSet.add_param("postexec_script",STR,false,false,0,0);
  _paramSet.add_param("postexec_args",STR,false,false,0,0);
  _paramSet.add_param("random_seed",INT,false,false,0,0);
  _paramSet.add_param("replicates",INT,true,false,0,0);
  _paramSet.add_param("generations",INT,true,false,0,0);
  _paramSet.add_param("run_mode",STR,false,false,0,0);
}
//----------------------------------------------------------------------------------------
// build_allParams
// ----------------------------------------------------------------------------------------
void ParamManager::build_allParams ()
{
  list< SimComponent* >::iterator cmpt = _components.begin();

  _allParams.clear();
  
  _allParams.push_back(&_paramSet);
  
  while(cmpt != _components.end()) {
    _allParams.push_back( (*cmpt)->get_paramset() );
    cmpt++;
  }
}
//----------------------------------------------------------------------------------------
// param_consistency_check
// ----------------------------------------------------------------------------------------
bool ParamManager::param_consistency_check ()
{
  list<ParamSet*>::iterator current_paramset = _allParams.begin();
  
  bool check = true;
  
  while(current_paramset != _allParams.end()){
    if(!(*current_paramset)->check_consistency()){
      error("ParamManager::param_consistency_check::consistency not satisfied for \"%s\"\n",
            (*current_paramset)->getName().c_str());
      check = false;
    }
    current_paramset++;
  }
  return check;
}
//----------------------------------------------------------------------------------------
// set_parameters
// ----------------------------------------------------------------------------------------
bool ParamManager::set_parameters (map< string,string >& simparams, bool silent)
{
  _inputParams = simparams;
  list<ParamSet*>::iterator current_paramset = _allParams.begin();
  map< string,string >::iterator IT = _inputParams.begin();
  bool set = false;
  
  while(IT != _inputParams.end()) {
    
    current_paramset = _allParams.begin();
    
    while(current_paramset != _allParams.end()){
      set |= (*current_paramset)->set_param((string&)IT->first,IT->second);
      //is true if at least one paramset could set the param
      current_paramset++;
    }
    
    if(!set && !silent){
      //error("ParamManager::could not set param \"%s\"\n",IT->first.c_str());
      //    return false;
    }
    set = false;
    IT++;
  }
  
  return param_consistency_check();
}
//----------------------------------------------------------------------------------------
// get_paramset
// ----------------------------------------------------------------------------------------
ParamSet* ParamManager::get_paramset (string& name)
{
  list<ParamSet*>::iterator pset = _allParams.begin();
  
  while(pset != _allParams.end()) {
    
    if( (*pset)->getName().compare(name) == 0)
      return (*pset);
    
    pset++;
  }
  
  return NULL;
}

//----------------------------------------------------------------------------------------
// build_records
// ----------------------------------------------------------------------------------------
void ParamManager::build_records (map< string, vector<string> >& initParams)
{
  map< string, string > params;
  map< string, string > paramsToExpand;
  map< string, string >::iterator param_iter;
  map< string, vector<string> >::iterator Pit;
  unsigned int RecNb = 1, ArgNb, SeqParam = 0, BlockSize, seq_pos, num_seq_arg = 0;
  vector<unsigned int> sequence; //stores the number of args of each sequence parameters
  vector<string> currSeqArg;
  vector<unsigned int> currSeqPos;
  vector<string> currCombinArg;
  vector<unsigned int> currCombinPos;
  string NAME;
  string ARG;
  bool SEQUENCE = 0;
    
  //find sequential and combinatorial parameters:
  for(Pit = initParams.begin(); Pit != initParams.end(); ++Pit) {
    if(Pit->second.size() > 1) {
      //fetch the number of the current sequence parameter:
      sequence.push_back(Pit->second.size());
      num_seq_arg++;
      //increase the total number of simulations records:
      RecNb *= Pit->second.size();
    } else if(Pit->second[0].find_first_of('[') != string::npos) {
      //process sequential arg
    }
    
  }
  
  if(RecNb > 1)
    //we have a sequence of simulations:
    SEQUENCE = true;
  
  for(unsigned int i = 0; i < RecNb; ++i) {
    //now build the simulation records with the right params!
    //the map 'param' will get all the params used for one simulation
    //it is then added to the list of simulations' parameters map
    
    SeqParam = 0;//used as index of the 'sequence' vector
    
    for(Pit = initParams.begin(); Pit != initParams.end(); ++Pit) {
      
      if(!(Pit->first.compare("filename")==0) &&
         !(Pit->first.compare("stat")==0) &&
         !(Pit->first.compare("postexec_args")==0)) {
        
        //get the number of arguments for the current parameter:
        ArgNb = Pit->second.size();
        
        if(ArgNb > 1) {
          //the current param is a sequence param
          //increase the index of the sequence parameter
          SeqParam++;
          //then compute the right argument to give to the current simulation record:
          BlockSize = RecNb;
          
          for(unsigned int j = 0; j < SeqParam; ++j)
            BlockSize /= sequence[j];
          
          seq_pos = (i/BlockSize) % ArgNb;
          //we store the argument position in its sequence and the actual argument string to build the filename:
          currSeqPos.push_back( seq_pos );
          currSeqArg.push_back(Pit->second[ seq_pos ]);
          
          //assign the right argument to the parameter:
          //params[Pit->first] = Pit->second[ seq_pos ];
          
          ARG =  Pit->second[ seq_pos ];
          //we can expand a string argument of a sequential parameter as well
          if(ARG.find_first_of('%') != string::npos)
            paramsToExpand[Pit->first] = ARG;
          else 
            params[Pit->first] = ARG;
          
        } else if (ArgNb == 1) {
          //the current param isn't a sequence param but gets an argument
          //we might have to do some name expansion:
          ARG =  Pit->second[0];
          
          if(ARG.find_first_of('%') != string::npos)
            paramsToExpand[Pit->first] = ARG;
          else 
            params[Pit->first] = ARG;
          
        } else
          //current param has no argument (bool type param) we give it value 1 (true)
          params[Pit->first] = "1";
        
      } else if (Pit->first.compare("stat") == 0){
        
        params["stat"].assign("");
        //build a string with all the stats arguments
        for(unsigned int k = 0; k < Pit->second.size(); ++k)
          params["stat"] += Pit->second[k] + " ";
        
      } else if (Pit->first.compare("filename")==0) {
        
        if(SEQUENCE)
          NAME = Pit->second[0];
        else
          params["filename"] = Pit->second[0];
        
      } else if (Pit->first.compare("postexec_args") == 0){
        
        ARG = "";

        //build a string with all the postexec arguments
        for(unsigned int k = 0; k < Pit->second.size(); ++k)
          ARG += Pit->second[k] + " ";
          
        params["postexec_args"] = ARG;
        
#ifdef NEMOSUB
        if(ARG.find_first_of('%') != string::npos)
          paramsToExpand["postexec_args"] = ARG;
#endif
        
      }
    }
    
    if(SEQUENCE) { //perform expansion of param arguments containing the '%' marker
      
      params["filename"] = setFilename(NAME, i+1, currSeqArg, currSeqPos, true);
      
      for(param_iter = paramsToExpand.begin(); param_iter != paramsToExpand.end(); param_iter++)
        params[param_iter->first] = setFilename(param_iter->second, i+1, currSeqArg, currSeqPos, false); 
    }
    
    //add all the params previously computed to the main simulation recorder list (of maps)
    _simRecords.push_back(params);
    
    currSeqArg.clear();
    currSeqPos.clear();
  }
  /*
   RecNb = 1;
   for(RecIt = SimRecorder.begin();RecIt != SimRecorder.end();++RecIt){
     message("\nSimulation "<<RecNb++;
             for(P2 = RecIt->begin();P2 != RecIt->end();++P2)
             message("\n  "<<P2->first<<"\t"<<P2->second;
   }
*/
}
// ----------------------------------------------------------------------------------------
// setFilename
// ----------------------------------------------------------------------------------------
string ParamManager::setFilename(string& fstring, unsigned int sim,
                                 vector<string>& args,
                                 vector<unsigned int>& arg_no,
                                 bool check_arg_no)
{
  string out, tail, fmt;
  string::size_type pos = fstring.find_first_of('%');
  string::size_type next;
  unsigned int index, nstr_args = (pos != string::npos ? 1 : 0);
  bool add_sim_no = false;
  static bool has_warned = false;
  
  while(pos != string::npos)
    nstr_args += ( (pos = fstring.find('%', pos+1)) != string::npos);
  
  if(check_arg_no && nstr_args < args.size()) {
    if(!has_warned) warning("missing sequential arguments in filename parameter, adding simulation number to filename.\n");
    has_warned = true;
    add_sim_no = true;
  }
  
  pos = fstring.find_first_of('%');
  
  if(pos != string::npos) {
    
    if(pos > 0) out = fstring.substr(0, pos);
    
    tail = fstring.substr(pos+1, string::npos);
    
    fmt = stripFormatString(tail, index);
    
    if(index > args.size()) fatal("too many sequential arguments in \"%s\"\n",fstring.c_str());
    
    out += setArgString(fmt, args[index-1], arg_no[index-1]);
    
    next = tail.find_first_of('%');
    
    while(next != string::npos){
    
      out += tail.substr(0, next);
      
      tail = tail.substr(next+1, string::npos);
      
      fmt = stripFormatString(tail, index);
      
      if(index > args.size()) fatal("too many sequential arguments in \"%s\"\n",fstring.c_str());

      out += setArgString(fmt, args[index-1], arg_no[index-1]);
      
      next = tail.find_first_of('%');
    }
    
    out += tail.substr(0, next);

  } else
    out = fstring;
  
  if(add_sim_no) {
    
    ostringstream ostr;
    ostr << sim;    
    
    out += "-" + ostr.str();
  }
  
  return out;
}
// ----------------------------------------------------------------------------------------
// stripFormatString
// ----------------------------------------------------------------------------------------
string ParamManager::stripFormatString(string& str, unsigned int& index)
{
  string fmt;
  unsigned int digits;
  size_t fmt_end;
  
  //check for the presence of a format string, enclose by two \'
  if(str[0] == '\'') {
    
    fmt_end = str.find('\'', 1);
    
    if(fmt_end == string::npos) fatal("format string not closed in filename parameter\n");
    
    fmt = str.substr(1, fmt_end-1);
    
    str = str.substr(fmt_end +1, string::npos);
    
  } else {
    
    fmt = "";
    
  }
  
  index = (unsigned int) strtol(str.c_str(), NULL, 10);
  
  digits = (unsigned int) log10((double)index) +1;
  
  str = str.substr(digits, string::npos);
  
  return fmt;
  
}
// ----------------------------------------------------------------------------------------
// setArgString
// ----------------------------------------------------------------------------------------
string ParamManager::setArgString(string& fmt, string& arg, unsigned int arg_pos)
{
  unsigned int width, digits;
  double value;
  bool is_dotless = 0;
  string out, arg_str;
  size_t pos = 0;
  
  if(fmt.size() != 0) {
    
    if(fmt[0] == '.') {
      is_dotless = true;
      fmt = fmt.substr(1, string::npos);
    }
    
    if( !isdigit(fmt[0]) ) {
      error("syntax error in format string for filename parameter:\n");
      fatal("first character within \' \' must be a number\n");
    }
    
    width =  (unsigned int) strtol(fmt.c_str(), NULL, 10);
    
    digits = (unsigned int) log10((double)width) +1;
    
    fmt = fmt.substr(digits, string::npos);
        
    if(fmt[0] == '[') {
      
      if( (pos = fmt.find(']')) == string::npos){
        error("syntax error in format string for filename parameter:\n");
        fatal("enclosing ']' was not found\n");
      }
      
      arg_str = fmt.substr(1, pos-1);
      
      if(arg_str[0] == '+') {
        
        ostringstream ostr;
        ostr.width(width);
        ostr.fill('0');
        ostr << arg_pos+1;
        
        out = ostr.str();
        
      } else {
        
        if (arg_str.size() < (arg_pos+1)*width) {
          error("syntax error in format string for filename parameter:\n");
          fatal("the number of characters within [ ] is not sufficient.\n");
        }
        
        for(unsigned int i = 0, start = arg_pos*width; i < width; i++)
          out += arg_str[start + i];
      }
      
    } else {
      
      value = strtod(arg.c_str(), 0);
      
      if(is_dotless) {
        //take the decimal part:
        double dummy;
        value = modf(value, &dummy);
        value *= pow(10.0, (double)width);
      }
      
      ostringstream ostr;
      ostr.width(width);
      ostr.fill('0');
      ostr << value;
      
      out = ostr.str();
    }
    
  } else {
    
    out = arg;
    
  }
  
  return out;
}
// ----------------------------------------------------------------------------------------
// lowercase
// ----------------------------------------------------------------------------------------
string ParamManager::lowercase(string& input)
{
  for(unsigned int i=0;i<input.size();++i)
    input[i] = tolower(input[i]);
  return input;
}
/* /_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/ */

//                             ******SimBuilder*******

//----------------------------------------------------------------------------------------
// copy cstor
// ---------------------------------------------------------------------------------------
SimBuilder::SimBuilder (const SimBuilder& SB)
{  
  list< TraitPrototype* >::const_iterator TT = SB._TTrait_Templates.begin();
  
  for(;TT != SB._TTrait_Templates.end(); TT++)  add_trait( (*TT)->clone() );
  
  list< LifeCycleEvent* >::const_iterator LCE = SB._LCE_Templates.begin();
  
  for(;LCE != SB._LCE_Templates.end(); LCE++)   add_LCE( (*LCE)->clone() );
    
  
  _paramSet = SB._paramSet;
  
  build_allParams();
}
//----------------------------------------------------------------------------------------
// destructor
// ---------------------------------------------------------------------------------------
SimBuilder::~SimBuilder ()
{
#ifdef _DEBUG_
  message("SimBuilder::~SimBuilder\n");
#endif
  
  _components.clear();
  
  list< TraitPrototype* >::iterator TT_it = _TTrait_Templates.begin();
  
  while(TT_it != _TTrait_Templates.end()) {
    delete (*TT_it);
    TT_it++;
  }
  
  list< LifeCycleEvent* >::iterator LCE_it = _LCE_Templates.begin();
  
  while(LCE_it != _LCE_Templates.end()) {
    delete (*LCE_it);
    LCE_it++;
  }
    
}
//----------------------------------------------------------------------------------------
// get_current_trait
// ----------------------------------------------------------------------------------------
TraitPrototype* SimBuilder::get_current_trait (trait_t type)
{
  TRAIT_ITER trait = _currentTraits.find(type);
  
  if(trait != _currentTraits.end())
    return trait->second;
  else
    return NULL;
}
//----------------------------------------------------------------------------------------
// get_current_event
// ----------------------------------------------------------------------------------------
LifeCycleEvent* SimBuilder::get_current_event (string& name)
{
  LCE_ITER LCE = _LifeCycle.begin();
  
  while(LCE != _LifeCycle.end()) {
    if( LCE->second->get_event_name().compare(name) == 0)
      return LCE->second;
    LCE++;
  }
  
  return NULL;
}
//----------------------------------------------------------------------------------------
// build_currentParams
// ----------------------------------------------------------------------------------------
bool SimBuilder::build_currentParams (map< string,string >& simparams)
{
  if(! this->set_parameters(simparams, false) ) return false;
  
  //empty the current parameters container:
  _currentParams.clear();
  
  list<ParamSet*>::iterator current_paramset = this->_allParams.begin();
  
  //add all set parameters to the list of the current params:
  while(current_paramset != this->_allParams.end()){
    if((*current_paramset)->isSet()) {
#ifdef _DEBUG_
      message("SimBuilder::build_currentParams:adding paramset %s\n",
              (*current_paramset)->getName().c_str());
#endif
      _currentParams.push_back(*current_paramset);
    }
    current_paramset++;
  }
  
  return true;
}
//----------------------------------------------------------------------------------------
// build_currentTraits
// ----------------------------------------------------------------------------------------
map< trait_t,TraitPrototype* >& SimBuilder::build_currentTraits()
{
  //build the list of Traits for this simulation:
  list< TraitPrototype* >::iterator trait = this->_TTrait_Templates.begin();
  
  _currentTraits.clear();
  
  while(trait != this->_TTrait_Templates.end()) {
    
    if( (*trait)->get_paramset()->isSet() ){
      _currentTraits[(*trait)->get_type()] = (*trait);
    }
    
    trait++;
  }
  
  return _currentTraits;
}
//----------------------------------------------------------------------------------------
// build_currentLifeCycle
// ----------------------------------------------------------------------------------------
void SimBuilder::build_LifeCycle()
{
  //build the list of the life cycle events for this simulation:
  list< LifeCycleEvent* >::iterator LCE = this->_LCE_Templates.begin();
  
  _LifeCycle.clear();
  
  
  while(LCE != this->_LCE_Templates.end()) {
    if( (*LCE)->get_paramset()->isSet() ){
#ifdef _DEBUG_
    cout <<"adding LCE "<< (*LCE)->get_event_name() <<" : "<< (*LCE)->get_rank() << endl;
#endif
      _LifeCycle[(int)(*LCE)->get_rank()] = (*LCE);
    }
    LCE++;
  }

}
//----------------------------------------------------------------------------------------
// getFirstRequiredAgeInLifeCycle
// ----------------------------------------------------------------------------------------
age_t SimBuilder::getFirstRequiredAgeInLifeCycle ( )
{
  age_t requiredAgeClass = NONE;

  for(LCE_ITER IT = _LifeCycle.begin();
      IT != _LifeCycle.end() && requiredAgeClass == NONE;
      IT++)
    requiredAgeClass = IT->second->requiredAgeClass();
  
  return requiredAgeClass;
}
