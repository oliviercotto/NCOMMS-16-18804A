/**  $Id: simulation.cc,v 1.11.2.4 2016-02-09 14:04:19 fred Exp $
 *
 *  @file simulation.cc
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
 *  created on @date 21.07.2004
 *
 *  @author fred
 */

#include <iostream>
#include <iomanip>
#include <sstream>
#include <time.h>
#include <cerrno>
#include "simulation.h"
#include "metapop.h"
#include "MPImanager.h"
#include "fileservices.h"
#include "statservices.h"
#include "servicenotifiers.h"
#include "Uniform.h"
#include "output.h"
#include "version.h"
#include "lifecycleevent.h"
#include "tstring.h"

MPIenv *_myenv = 0;

#ifndef HAS_SPRNG
#ifdef HAS_GSL
gsl_rng * RAND::r = 0;
#else
long RAND::Seed1 = 0;
long RAND::Seed2 = 98280582;
#endif
#endif

//----------------------------------------------------------------------------------------
// init
// ----------------------------------------------------------------------------------------
bool SimRunner::init()
{
  if( !(_paramSet.isSet("filename")) ) {
	error("parameter \"filename\" is not properly set!\n");
	return false;
  }
  
  _FileServices.setBasename(_paramSet.getArg("filename"));
  
  if(_paramSet.isSet("root_dir")) {
    _FileServices.setRootDir(_paramSet.getArg("root_dir"));
  }
  
  if(_paramSet.isSet("run_mode")) {
    
    _modeArg = _paramSet.getArg("run_mode");
    
    if(_modeArg == "run"){
      _FileServices.setMode(0);
      _doRun = true;
      _mode = 0;
      SILENT_RUN = false;  //declared in output.h
    } else if( _modeArg == "overwrite" ) {
      _FileServices.setMode(1);
      _doRun = true;
      _mode = 1;
      SILENT_RUN = false;
    } else if(_modeArg == "skip"){
      _FileServices.setMode(2);
      _doRun = true;
      _mode = 2;
      SILENT_RUN = false;
    } else if(_modeArg == "dryrun"){
      _FileServices.setMode(3);
      _doRun = false;
      _mode = 3;
      SILENT_RUN = false;
    } else if(_modeArg == "create_init"){
      _FileServices.setMode(4);
      _doRun = false;
      _mode = 4;
      SILENT_RUN = false;
    } else if(_modeArg == "silent_run"){
      _FileServices.setMode(0);
      _doRun = true;
      _mode = 0;
      SILENT_RUN = true;
    } else {
      error("simulation run mode \"%s\" unknown\n",_modeArg.c_str());
      return false;
    }
  } else {
    _modeArg = "run";
    _doRun = true;
    _FileServices.setMode(0);
    SILENT_RUN = false;
  }
  
  _replicates = (unsigned int)_paramSet.getValue("replicates");
  
  _generations = (unsigned int)_paramSet.getValue("generations");
  
  
  if( _paramSet.isSet("logfile") ) {
    _logfile = _paramSet.getArg("logfile");
  } else 
    _logfile = "nemo.log";
  
  printLogHeader();
  
  if( _paramSet.isSet("postexec_script") ) {
    _postexec_script = _paramSet.getArg("postexec_script");
    _do_postexec = true;
    
    if( _paramSet.isSet("postexec_args") )
      _postexec_args = _paramSet.getArg("postexec_args");
    else
      _postexec_args = "";

  } else
    _do_postexec = false;
  
  return true;
}
//----------------------------------------------------------------------------------------
// dstor
// ----------------------------------------------------------------------------------------
SimRunner::~SimRunner ()
{
#ifdef _DEBUG_
  message("SimRunner::~SimRunner\n");
#endif
  
  reset();
  
}
//----------------------------------------------------------------------------------------
// reset
// ----------------------------------------------------------------------------------------
void SimRunner::reset()
{
#ifdef _DEBUG_
  message("SimRunner::reset\n");
#endif
  //reset all the parameter to the "unset" state
  list<ParamSet*>::iterator current_paramset = this->_allParams.begin();
  while(current_paramset != this->_allParams.end()) {
    (*current_paramset)->reset();
    current_paramset++;
  }
  
  reset_services();
}
//----------------------------------------------------------------------------------------
// reset_services
// ----------------------------------------------------------------------------------------
void SimRunner::reset_services()
{
  _FileServices.reset();
  _StatServices.reset();
  _ParamUpdaterManager.reset();
#ifdef _DEBUG_
  message("SimRunner::reset_services::done\n");
#endif
}
//----------------------------------------------------------------------------------------
// register_services
// ----------------------------------------------------------------------------------------
void SimRunner::register_component (SimComponent* cmpt)
{
  _FileServices.load(cmpt);
  _StatServices.load(cmpt);
  _ParamUpdaterManager.load(cmpt);
}
//----------------------------------------------------------------------------------------
// register_components
// ----------------------------------------------------------------------------------------
void SimRunner::register_component_handlers ()
{
  register_component(_thePop);
  
  for(TRAIT_ITER trait = _currentTraits.begin(); trait != _currentTraits.end(); trait++)
    register_component(trait->second);
  
  for(LCE_ITER LCE = _LifeCycle.begin(); LCE != _LifeCycle.end(); LCE++)
    register_component(LCE->second);  
}
// ----------------------------------------------------------------------------------------
// setLifeCycle
// ----------------------------------------------------------------------------------------
void SimRunner::setLifeCycle( )
{
  build_LifeCycle();

  for( LCE_ITER LCE = _LifeCycle.begin(); LCE != _LifeCycle.end(); LCE++)
    LCE->second->init(_thePop);
}
//----------------------------------------------------------------------------------------
// init_components
// ----------------------------------------------------------------------------------------
bool SimRunner::init_components (map< string,string >& simparams)
{
#ifdef _DEBUG_
  message("SimRunner::init_components\n");
#endif  
  //first reset all paramSets and services (clear the handlers' lists)
  reset();
  
  _FileServices.set_pop_ptr(_thePop);
  _StatServices.set_pop_ptr(_thePop);
  
#ifdef _DEBUG_
  message("SimRunner::init_components: building current params\n");
#endif
  
  //build the list of active component parameters from the simulation record:
  if(!this->build_currentParams(simparams)){
    error("SimRunner::init_components:couldn't build current params\n");
    return false;
  }
  
  //initialize the random generator's seed:
  if(_myenv->isMaster()) init_random_seed();
  
  //init the sim and pop (set parameters from input values)
  if( !init() || !_thePop->init() ) return false;
  
  //propagate replicate and generation numbers to metapop, often used by LCEs
  _thePop->setReplicates(_replicates);
  _thePop->setGenerations(_generations);

  //the traits and LCEs components are set only after the simulation and 
  //population parameters are set
  
  //build the Individual prototype, and init the TTraits
  //---> SimComponents::setParameters() is called here, this will also build the genetic map
  _thePop->makePrototype( build_currentTraits() );
  
  //build the list of life cycle events, and init LCEs 
  //---> SimComponent::setParameters() is called here
  setLifeCycle( );
  
  //load the stat-, file-, and updater-handlers of the simulation components:
  register_component_handlers();
  
  //StatServices: build the lists of stat recorders:
  if( !_StatServices.init() ) return false;
  
  if( !_ParamUpdaterManager.init() ) return false;
  
  if( _ParamUpdaterManager.hasTemporals() ) {
    LCE_ParamUpdaterNotifier* updater = new LCE_ParamUpdaterNotifier();
    updater->setManager( &_ParamUpdaterManager );
    updater->init(_thePop);
    _LifeCycle[ -1 ] = updater;
  }
  
  return true;
}
//----------------------------------------------------------------------------------------
// init_random_seed
// ----------------------------------------------------------------------------------------
void SimRunner::init_random_seed ()
{
    
  if(this->_paramSet.isSet("random_seed")) {
    _random_seed = (unsigned long) this->_paramSet.getValue("random_seed");
    message( "setting random seed from input value: %i\n", _random_seed);
  } else
    _random_seed =  2*time(0)+1;
  
  RAND::init(_random_seed);
}
//----------------------------------------------------------------------------------------
// run
// ----------------------------------------------------------------------------------------
bool SimRunner::run ( int ARGC, char **ARGV )
{  
  //initialize the MPI environment:  
  _myenv = new MPIenv( ARGC, ARGV, _my_mpi_manager );    
  
  //--------------------------------------------------------------------
  if(_myenv->isMaster()) {
    message("\n  N E M O %i.%i.%i%s %s\n", MAIN_VERSION, MINOR_VERSION, REVISION, RELEASE, VERSION_DATE);
    message("\n  Copyright (C) 2006-2014 Frederic Guillaume");
    message("\n  This is free software; see the source for copying");
    message("\n  conditions. There is NO warranty; not even for ");
    message("\n  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.");
    message("\n  http://nemo2.sourceforge.net");
    message("\n------------------------------------------------\n");
  }
  
  //--------------------------------------------------------------------  
  //check for the presence of a metapop:
  if(_thePop == NULL) fatal("SimRunner::run: no population is attached at startup!\n");
  
  //--------------------------------------------------------------------  
  //build the list of params from the components list
  this->build_allParams();
  
  //--------------------------------------------------------------------  
  //get the input parameters
  FileParser Reader("");
  if (ARGC == 1)
    build_records(Reader.getParsedParameters("Nemo2.ini"));
  else
    for (int i = 1; i < ARGC; ++i)
      build_records(Reader.getParsedParameters(ARGV[i]));
  
  //--------------------------------------------------------------------  
  //run the simulation
  bool status = run();
  
  _myenv->finish(_my_mpi_manager);
  
  delete _myenv;
  
  return status;
}
//----------------------------------------------------------------------------------------
// run
// ----------------------------------------------------------------------------------------
bool SimRunner::run ( )
{
  time_t t;
  unsigned int sim = 0, simnbre = _simRecords.size();
  
  list< map< string,string > >::iterator currentSim = _simRecords.begin();
  
  //first loop: perform all simulations contained in _simRecords:
  //---------------------------------------------------------------------------------------
  while(currentSim != _simRecords.end()) {
    
    sim++;
    
    //clear and build the lists of traits, LCEs, stats and files handlers, random seed
    //create trait prototypes and the genetic map within IndFactory::makePrototype
    if(!init_components(*currentSim))  return false;
    
    //output a few info:
    if ( _myenv->isMaster() ) {
      t = time(NULL);
      strftime(_startTime, 20, "%d-%m-%Y %H:%M:%S", localtime(&t));
      message("\n--- SIMULATION %i/%i ---- [ %s ]\n\n",sim,simnbre,_FileServices.getBaseFileName().c_str());
      message("    start: %s\n",_startTime);
      message("    mode: %s\n",_modeArg.c_str());
      message("    traits: ");
      
      map<trait_t, TraitPrototype* >::iterator trait = _currentTraits.begin();
      while(trait != _currentTraits.end()) {
        message("%s",trait->first.c_str());
        trait++;
        if(trait != _currentTraits.end()) message(", ");
      }
      message("\n    LCEs: ");
      
      for(LCE_ITER LCE = _LifeCycle.begin(); LCE != _LifeCycle.end(); LCE++) {
        message("%s(%i)",LCE->second->get_event_name().c_str(),LCE->first);
        if(LCE != _LifeCycle.end()) message(", ");
      }
      message("\n");
    }
    
    //init the file services:
    //->files check, false means the user wants to skip this simulation.
    if( !(_FileServices.init(this->_currentParams)) ) {
      currentSim++;
      _thePop->clearPrototype();
      continue;
    }
    
    if(_thePop->isSourceLoad()) {
      message("    loading population from: %s_*%s\n",_thePop->getSourceName().c_str(),
              _thePop->getSourceFileType().c_str());
      //should check if trait with map before issuing following warning...
      if(_thePop->getIndividualProtoype()->getTraitNumber() != 0) {
        warning("while resetting genetic map, the locus positions saved in source population are ignored");
        warning("loci positions are set from parameters in the init file only");
        warning("please check manually for mismatch\n");
      }
    }
    message("\n");
    
    //->save simparameters in log files, add the current random seed
    _FileServices.log_simparams();
    _FileServices.log("#random_seed "+tstring::ulong2str( _random_seed ));
    
    clock_t start = clock();
    
    //run the simulation
    //-------------------------------------------------------------------------------------
    if(_doRun) Replicate_LOOP( );
    
    clock_t stop = clock();
    
    t = time(NULL);
    
    strftime(_endTime, 20, "%d-%m-%Y %H:%M:%S", localtime(&t));
    
    _simElapsedTime = setElapsedTime(stop - start);
    
    if ( _myenv->isMaster() && _mode != 3 && _mode != 4) {
      printLog();
      message("\n\n    end: %s\n",_endTime);
      message("--- done (CPU time: %ss)\n",_simElapsedTime.c_str());
      //don't log this when in 'dryrun' or 'create_init' mode:
      _FileServices.log("\n# simulation finished " + string(_endTime));
      _FileServices.log("# CPU time used: " + _simElapsedTime + "s");
    }
     
    _thePop->clearPrototype();
    
    currentSim++;
  }
  //---------------------------------------------------------------------------------------
  if( _do_postexec && _myenv->isMaster() ) runPostExec();
  
  return true;
}
// ----------------------------------------------------------------------------------------
// Replicate_LOOP
// ----------------------------------------------------------------------------------------
void SimRunner::Replicate_LOOP()
{
  
#ifndef USE_MPI
  
  clock_t start;
  clock_t stop;
  time_t t;
  char t_buff[10];
  
  _meanReplElapsedTime = 0;
  _meanGenLength = 0;
  _currentGeneration = 0;
  
  //---------------------------------- REPLICATE LOOP -------------------------------------
  for(_currentReplicate = 1; !(_currentReplicate > _replicates); ++_currentReplicate) {
    
    
    // -------------------------- PRINT STUFF ------------------------------
    t=time(0);
    //print some output:
    strftime(t_buff,10,"%H:%M:%S",localtime(&t));
    
#ifdef LOW_VERBOSE
    message("    replicate %i/%i [%s] \n",_currentReplicate, _replicates, t_buff);
#else
    message("\r    replicate %i/%i [%s] 1/%i                ",_currentReplicate, _replicates, t_buff, _generations);
    fflush(stdout);
#endif
    
    //---------------- SET POPULATION FOR FIRST GENERATION ------------------
    
    setForFirstGeneration();
     
    //--------------------------- GENERATION LOOP ---------------------------   
    start = clock();
    
    Cycle(t_buff);
    
    stop = clock();
    //-----------------------------------------------------------------------
    
    _meanReplElapsedTime += (stop - start);
        
    _meanGenLength += _currentGeneration - 1;
    
    //call the file services to print the replicate stats in case of pop extinction 
    if( !_thePop->isAlive() && _currentGeneration-1 < _generations ) {
      _currentGeneration = _generations;
      _thePop->setCurrentGeneration(_generations);
      _FileServices.notify();
    }
    
  }
  //--------------------------------- /REPLICATE LOOP -------------------------------------
  
  if( !_thePop->isAlive() && _currentGeneration < _generations ) {
    _currentGeneration = _generations;
    _currentReplicate  = _replicates;
    _thePop->setCurrentGeneration(_generations);
    _thePop->setCurrentReplicate(_replicates);
    _FileServices.notify();
  }
  
  _meanReplElapsedTime /= _replicates;
  _meanGenLength /= _replicates;

#else
  //else ifndef USE_MPI
  //------------------------------- MPI REPLICATE LOOP ------------------------------------
  _currentReplicate = _my_mpi_manager->init( &_StatServices );
  
  while ( _currentReplicate <= _replicates )
    _my_mpi_manager->iterate( this, &_StatServices, &_currentGeneration, &_currentReplicate );
 
  _my_mpi_manager->finish( &_StatServices, &_currentGeneration, _currentReplicate );
  
  if( _myenv->isMaster() ) {
    _currentGeneration = _generations;
    _currentReplicate  = _replicates;
    _thePop->setCurrentGeneration(_generations);
    _thePop->setCurrentReplicate(_replicates);
    
    //write the stat files, only if "save_stats" LCE is part of the life cycle
    if(get_LCE("save_stats")->get_paramset()->isSet()) {
      FileHandler& statWriter = dynamic_cast<LCE_StatServiceNotifier*> ( get_LCE("save_stats") )->getFH();
      for (unsigned int i = 1; i <= _replicates; ++i) {
        _currentReplicate = i;
        statWriter.FHwrite();
      }
    }
  }
  //------------------------------ /MPI REPLICATE LOOP ------------------------------------
#endif
  //endif USE_MPI
  
  //delete all individuals present in the population and delete the patches:
  _thePop->clear();
}
// ----------------------------------------------------------------------------------------
// setForFirstGeneration
// ----------------------------------------------------------------------------------------
void SimRunner::setForFirstGeneration()
{
  _thePop->setCurrentGeneration(0);
  
  //do some memory clean-up:
  _thePop->purgeRecyclingPOOL();

  //reset temporal parameters/components to their first generation value/state
  _ParamUpdaterManager.notify(0);
  
  //reset stat recording iterator to first occurrence
  _StatServices.resetCurrentOccurrence();
  
  //build metapopulation for the current replicate (build first generation)
  _thePop->setPopulation(_currentReplicate, _replicates);  
}
// ----------------------------------------------------------------------------------------
// Cycle
// ----------------------------------------------------------------------------------------
void SimRunner::Cycle(char* startTime)
{
  
  // ------------------------------ CYCLE --------------------------------
  
  for(_currentGeneration = 1; !(_currentGeneration > _generations); _currentGeneration++) {
	
    // -------------------------- PRINT STUFF ------------------------------
#if !defined(LOW_VERBOSE) && !defined(USE_MPI)
    if( !(_currentGeneration % 100) || _currentGeneration < 100){
      message("\r    replicate %i/%i [%s] %i/%i", _currentReplicate, _replicates
			  ,startTime, _currentGeneration, _generations);
	  fflush(stdout);
	}
#endif
    
#ifdef _DEBUG_
    message("____Generation %i/%i____\n", _currentGeneration, _generations);
#endif    
    
    // -------------------------- STEP ONE GEN ------------------------------
    _thePop->setCurrentGeneration(_currentGeneration);

    // do some memory clean-up along the way, to avoid keeping too many individuals in the pool
    if( !(_currentGeneration % 500)) _thePop->purgeRecyclingPOOL();
    
    // do one iteration of the life cycle:
    if(_thePop->isAlive())
      
      step(1);
    
    else {

#if !defined(LOW_VERBOSE) && !defined(USE_MPI)
      message("\r    replicate %i/%i [%s] %i/%i -> Pop extinction !\n", _currentReplicate, _replicates
              ,startTime, _currentGeneration, _generations); 
#endif

      _currentGeneration++;

      break;
    }

  }
  // --------------------------- END OF CYCLE --------------------------
}
//----------------------------------------------------------------------------------------
// run_event
// ----------------------------------------------------------------------------------------
bool SimRunner::run_event (string& name)
{
  //before calling that fction, first set the params and call init_components()!!
  
  LifeCycleEvent* event = this->get_current_event(name);
  
  if(event == NULL) {
	error("SimRunner::run_event:event \"%s\" not found in set events!\n",name.c_str());
	return false;
  }
  
  if( !(event->get_paramset()->isSet()) ) return false;//highly unlikely!!
  
  event->execute();
  
  return true;
}
//----------------------------------------------------------------------------------------
// step
// ----------------------------------------------------------------------------------------
void SimRunner::step (unsigned int nb_gen)
{
  LCE_ITER LCE = _LifeCycle.begin();
  
  while(LCE != _LifeCycle.end()) {
    _currentRankInLifeCycle = LCE->second->get_rank();
     LCE->second->execute();
    _thePop->setCurrentAge(LCE->second);
    LCE++;
  }
}
//----------------------------------------------------------------------------------------
// displayElapsedTime
// ----------------------------------------------------------------------------------------
std::string SimRunner::setElapsedTime(clock_t time)
{
  int e_time = time / CLOCKS_PER_SEC;
  int hour =  e_time / 3600;
  int minute =  ((e_time % 3600) / 60);
  int sec = (e_time % 3600) % 60;
  
  std::ostringstream out(ios::out);
  
  out.fill('0');
  out.width(2);
  out<<hour<<":";
  out.width(2);
  out<<minute<<":";
  out.precision(2);
  out.fill('0');
  out.width(2);
  out<<sec;
  
  return(out.str());
}

// ----------------------------------------------------------------------------------------
// SimRunner::printLogHeader()
// ----------------------------------------------------------------------------------------
void SimRunner::printLogHeader()
{
  ofstream FH;
  ifstream IF;
  //check is the logfile already exists:
  IF.open(_logfile.c_str(),ios::in);
  if(IF){
    IF.close();
    return;
  }
  IF.close();
  
  FH.open(_logfile.c_str(),ios::out);
  if(!FH) {
    error("could not create simulation logfile \"%s\"\n",_logfile.c_str());
    return;
  }
  
  FH<<"--- N E M O ---\n"
  <<"    LOGFILE\n\n\n";
  FH<<"| basename                                |      start time     |      stop time      | e-time CPU |"
  <<" repl done | rpl e-time | mean gen | version               | hostname             | output files \n";
  
  FH.close();
}
// ----------------------------------------------------------------------------------------
// SimRunner::printLog()
// ----------------------------------------------------------------------------------------
void SimRunner::printLog()
{
  ofstream FH(_logfile.c_str(),ios::app);
  
  if(!FH.is_open()){
    error("could not open simulation logfile \"%s\"\n",_logfile.c_str());
    return;
  }
  
  FH<<"| ";
  FH.width(40);
  FH.setf(ios::left,ios::adjustfield);
  FH << _FileServices.getBaseFileName( );
  
  FH <<"| "<< _startTime <<" | "<< _endTime <<" | ";
  
  FH.width(10);
  FH.setf(ios::right,ios::adjustfield);
  FH<< _simElapsedTime <<" | ";
  
  FH.width(9);
  FH << _replicates <<" | ";
  
  FH.width(10);
  FH << setElapsedTime( _meanReplElapsedTime ) << " | ";
  
  FH.width(8);
  FH << _meanGenLength <<" | ";
  
  FH<<" "<<MAIN_VERSION<<"."<<MINOR_VERSION<<"."<<REVISION<<RELEASE
  <<" "<<VERSION_DATE;
  
  FH<<" | ";
  FH.width(20);
  FH.setf(ios::left,ios::adjustfield);
  char* host;
  if( (host = getenv("HOST")) != NULL )
    FH << host << " |";
  else if ( (host = getenv("HOSTNAME")) != NULL )
    FH << host << " |";
  else
    FH << "-" << " |";
  
  FileServices::file_it file = _FileServices.getFirstWriter(), last = _FileServices.getLastWriter() ;
  
  for(;file != last; file++)
    FH << " \"" << (*file)->get_extension() << "\":" << (*file)->get_path();
  
  FH << std::endl;
  
  FH.close();
}
// ----------------------------------------------------------------------------------------
// SimRunner::runPostExec()
// ----------------------------------------------------------------------------------------
void SimRunner::runPostExec ( )
{
  ifstream script(_postexec_script.c_str(),ios::in);
  string cmd;
  
  if(!script.is_open()) {
	error("could not open post simulation shell script!\n");
	return;
  }
  
  message("Executing shell script \"%s\" ",_postexec_script.c_str());
  fflush(stdout);
  
  cmd = "sh " + _postexec_script + " " + _postexec_args;
  
  if(system(cmd.c_str()) < 0){
	error("execution of `sh %s %s' failed: %s\n",_postexec_script.c_str(), _postexec_args.c_str(),strerror(errno));
	return;
  }
  
  message("...done\n");
}
