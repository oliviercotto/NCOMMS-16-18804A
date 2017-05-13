/**  $Id: simulation.h,v 1.9.2.1 2014-04-29 18:14:57 fred Exp $
*
*  @file simulation.h
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


#ifndef SIMULATION_H
#define SIMULATION_H

#include <list>
#include <map>
#include <vector>
#include <string>
#include <time.h>
#include "basicsimulation.h"
#include "fileservices.h"
#include "statservices.h"
#include "updaterservices.h"
#include "metapop.h"


/**Performs the setup of the Metapop and SimComponents and runs the simulation. 
 * This class implements the two main loops of a simulation, the replicate and the generation loops. The replicate
 * loop iterates the generation loop which itself iterates the life cycle loop composed of the LCEs selected by the user. 
 * A SimRunner brings together the basic simulation components and a metapopulation on which they act. 
 * It perfoms the setups necessary to have a Metapop ready for the simulation and runs the different simulations
 * stored in its ParamManager base class. Also hosts the file and stat services. 
 **/
class SimRunner: public SimBuilder {
private:
    
  MPImanager*                      _my_mpi_manager;
  
  Metapop*                         _thePop;
  
  std::string                      _logfile;
  
  char                             _startTime[20], _endTime[20];
  
  std::string                      _simElapsedTime;
  /**Clock counter, for logging.**/
  clock_t _meanReplElapsedTime;
  /**Generation counter, for logging.**/
  unsigned int _meanGenLength;
  
  //parameters:
  /**Number of generations to iterate.*/
  unsigned int _generations;
  /**Number of replicates to iterate.*/
  unsigned int _replicates;
  /**The current replicate in the replicate loop, starts at 1.*/
  unsigned int _currentReplicate;
  /**The current generation in the generation loop, starts at 1.*/
  unsigned int _currentGeneration;
  /**The current rank in the life cycle, corresponds to the rank of the current LCE, before it executes.*/
  int _currentRankInLifeCycle;
  /**The startup random seed of the random generator*/
  unsigned long _random_seed;
  /**The run mode ('overwrite', 'run', 'skip', 'dryrun', 'create_init').*/
  std::string                      _modeArg;
  /**The run mode code (0 = run, 1 = overwrite, 2 = skip, 3 = dryrun, 4 = create_init)*/
  unsigned int                     _mode;
  /**Boolean set when not in dryrun mode.*/
  bool                             _doRun;
  /**The path to the script to be executed after last simulation.*/
  std::string                      _postexec_script;
  /**Arguments to pass to the post-exec script.*/
  std::string                      _postexec_args;
  /**Boolean set when a post-exec script must be executed.*/
  bool                             _do_postexec;

public:
    
  FileServices					   _FileServices;
  
  StatServices					   _StatServices;

  UpdaterServices          _ParamUpdaterManager;
  
public:
  
  SimRunner   (Metapop* pop) : _my_mpi_manager(0)
  {
    attach_pop(pop);
  }
  /**Dstror.*/
                                   ~SimRunner            ( );
  /**Checks simulation parameters and init the FileServices with the base filename. @callgraph**/
  bool                             init                   ( );
  /**Performs the initialization of the different components of the simulation.
    *Builds the list of the simulation parameters and load the components that have their ParamSet in the set state.
    *Builds the population, the TraitPrototype, and life cycle, register the various services and init 
    *the StatServices and ParameterUpdaterManager.
    *@param simparams the hashtable containing the parameters and their arguments parsed from the init file
    */
  bool                             init_components        (map< string,string >& simparams);
  
  /**Sets the list of LifeCyckeEvent's currently active.*/
  void                             setLifeCycle           ( );
  
  /**Sets the population and the services ready for the first generation of a new replicate.*/
  void                             setForFirstGeneration  ( );
  
  /**Initialize the seed of the random generator. */
  void                             init_random_seed       ( );

  /**Resets all the parameters to the unset state, resets the services.*/
  void                             reset                  ( );
  
  /**Compute and print the simulation's elapsed time to stdout. 
    @param time elapsed time in ticks count*/
  std::string                      setElapsedTime         (clock_t time);
  
  /**Register the different Handler's attached to a SimComponent.
   @param cmpt a SimComponent**/
  void                             register_component      (SimComponent* cmpt);
  /**Register all the Handlers of the currently active simulation components */
  void                             register_component_handlers ( );
  /**Resets the FileServices and StatServices.*/
  void                             reset_services         ( );
  /**Returns the FileServices. **/
  FileServices*			           get_FileServices       ( )                      {return &_FileServices;}
  /**Returns the StatServices. **/
  StatServices*			           get_StatServices       ( )                      {return &_StatServices;}
  /**Returns the complete list of the stat recorders loaded after parameters initialization. **/
  list<StatRecBase*>               get_allRegisteredStats ( )                      {return _StatServices.getAllStats();}
  
  /**First loop of the simulation, performs the simulations stored in the ParamManager base class. */
  bool                             run                    ( int ARGC, char **ARGV );  
  /**First loop of the simulation, performs the simulations stored in the ParamManager base class. */
  bool                             run                    ( );
  /**Execute one specific life cycle event, if present in the list of current events.
    @param name the name of the (paramSet of the) LifeCycleEvent to execute.*/
  bool                             run_event              (string& name);
  /**Iterates the life cycle.
    @param nb_gen number of iterations to perform*/
  void                             step                   (unsigned int nb_gen);
  
  /**Attach a pop to the simulation. Adds it to the components list.
   * @param pop ptr to the pop object
   **/
  void                             attach_pop             (Metapop* pop) {_thePop = pop;this->_components.push_back(_thePop);}
  /**Accessor to the pop ptr. 
   * @return the pop ptr
   **/
  Metapop*                         get_pop                ( )                      {return _thePop;}
  /**Calls the Metapop init procedure with current traits and LCEs.  \callgraph*/
    /*!\callgraph 
    */
  bool                             build_pop              ( );
  
  ///@name Main loops
  ///@{
  /**Replicate loop, iterates the life cycle \a _replicates times.   */
  void Replicate_LOOP();
  
  /**Life cycle loop, executes the list of LCEs \a _generations times.
   @param startTime the starting time of the current replicate.
   */
  void Cycle(char* startTime);
  ///@}
  
  void         setCurrentGeneration (unsigned int gen)  {_currentGeneration = gen; _thePop->setCurrentGeneration(gen);}
  unsigned int getCurrentGeneration ()                  {return _currentGeneration;}
  void         setGenerations       (unsigned int gen)  {_generations = gen; _thePop->setGenerations(gen);}
  unsigned int getGenerations       ()                  {return _generations;}
  void         setCurrentReplicate  (unsigned int repl) {_currentReplicate = repl; _thePop->setCurrentReplicate(repl);}
  unsigned int getCurrentReplicate  ()                  {return _currentReplicate;}
  void         setReplicates        (unsigned int repl) {_replicates = repl; _thePop->setReplicates(repl);}
  unsigned int getReplicates        ()                  {return _replicates;}
  int     getCurrentRankInLifeCycle ()                  {return _currentRankInLifeCycle;}
  
  void printLogHeader ();
  void printLog ();
  void runPostExec();
};


#endif

