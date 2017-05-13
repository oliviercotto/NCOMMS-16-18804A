/**  $Id: statservices.cc,v 1.8.2.2 2014-04-29 17:49:23 fred Exp $
 *
 *  @file statservices.cc
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

#include <iostream>
#include <iomanip>
#include <sstream>
#include <list>
#include "statservices.h"
#include "stathandler.h"
#include "metapop.h"
#include "simenv.h"
#include "output.h"

extern MPIenv *_myenv;

//-----------------------------------------------------------------------------------------
// init
// ----------------------------------------------------------------------------------------
bool StatServices::init() 
{
  stat_it HIT;
  
  HIT = _statHandlers.begin();
  
  //first init the stat handlers -> will define the dims of the stat tables
  while(HIT != _statHandlers.end()) {
    (*HIT)->init();
    HIT++;
  }
  
  istringstream stats(_statArg);
  std::string arg;
  std::list<std::string> args;
  
  //extract the stat options from the string read in the init file
  while(stats>>arg) args.push_back(arg);
	
  std::list<std::string>::iterator token = args.begin();
  
  bool is_set = false;
  
  while(token != args.end()) {
    
    HIT = _statHandlers.begin();
    
    is_set = false;
    
    while(HIT != _statHandlers.end() && !is_set) {
      is_set |= (*HIT)->setStatRecorders( (*token) );
      HIT++;
    }
    
    if(!is_set) {
      error("the string \"%s\" is not a valid stat option\n",(*token).c_str());
      return false;
    }
    
    token++;
  }

  //---- set the list of stat recorders:
  _statRecorders.clear();
  
  for (HIT = _statHandlers.begin(); HIT != _statHandlers.end(); ++HIT) {

    _statRecorders.splice(_statRecorders.end(), (*HIT)->getStats());
    //this removes the recorders from the handler's list (in principle)
    //this means we need to delete the recorders within the StatServices!!
  } 
  
  _numRecorders = _statRecorders.size();

  
  //---- set iterator to first generation occurrence to record
  _current_occurrence = _occurrences.begin(); 
  //this is safe as the occurrence list is set before init() is called.
  
//  //reset output option for per-gen mean stats
//  _printAverages = 1;  //this resets user-defined option!!
  
  //---- set the stat table:
  reset_stat_table();
  
  for (unsigned int i = 0; i < SIMenv::getReplicates(); ++i) {
    _statValues.push_back(vector<double*>());
  }
  
  return true;
}
//-----------------------------------------------------------------------------------------
// setOccurrences
// ----------------------------------------------------------------------------------------
void StatServices::setOccurrences (std::map<unsigned int, unsigned int> timeTable)
{
  std::map<unsigned int, unsigned int>::const_iterator time = timeTable.begin(), next;
  
  _occurrences.clear();
  _occurrences.push_front(0);
  
  if (time->first == 0 && time->second != 0) { //this is a recursive logtime
    
    unsigned int next_term;
    
    for (;time != timeTable.end(); ++time) {
      
      next = time;
      next++;
      next_term = ( next != timeTable.end() ? 
                   next->first : SIMenv::getGenerations() );
      
      while( _occurrences.back() + time->second < next_term)
      {
        _occurrences.push_back( _occurrences.back() + time->second );
      }
      _occurrences.push_back(next_term);
    }
    
    
  } else if (time->second == 0) { //non-recursive logtime, only records specific times
    
    for (;time != timeTable.end(); ++time) _occurrences.push_back(time->first);
    
  } else { 
    
    fatal("something went wrong while setting stat recording times (first=%i, second=%i)\n",
          time->first, time->second);
  }
  
  //clean up:
  while (_occurrences.back() > SIMenv::getGenerations()) _occurrences.pop_back();
  
  if(_occurrences.back() != SIMenv::getGenerations()) //automatically record last gen
    _occurrences.push_back(SIMenv::getGenerations());
  
  if (_occurrences.front() == 0) _occurrences.pop_front();
  if (_occurrences.front() != 1) _occurrences.push_front(1); //automatically record first gen
}

//-----------------------------------------------------------------------------------------
// load
// ----------------------------------------------------------------------------------------
void StatServices::load ( SimComponent* sc )
{
  sc->loadStatServices(this);
}
//-----------------------------------------------------------------------------------------
// attach
// ----------------------------------------------------------------------------------------
void StatServices::attach ( Handler* H)
{
  StatHandlerBase* SH = dynamic_cast<StatHandlerBase*> (H);
  
  _statHandlers.push_back(SH);
  
  SH->set_service(this);
}
//-----------------------------------------------------------------------------------------
// notify
// ----------------------------------------------------------------------------------------
void StatServices::notify ()
{
#ifdef _DEBUG_
  message("StatServices::notify (%i/%i): ", SIMenv::getCurrentGeneration(),
          *_current_occurrence);
  fflush(stdout);
#endif
  
  if( *_current_occurrence == SIMenv::getCurrentGeneration() ) {
    recordStats( *_current_occurrence );
    _current_occurrence++;
  }
  
#ifdef _DEBUG_
  message("\n");
#endif
}
//-----------------------------------------------------------------------------------------
// reset
// ----------------------------------------------------------------------------------------
void StatServices::reset ()
{
  stat_it HIT = _statHandlers.begin();
  
  while(HIT != _statHandlers.end()) {
    (*HIT)->reset();
    HIT++;
  }
  
  _statHandlers.clear();
  
  //the stat recorders must be freed here, they have been removed from the StatHandlers
  rec_it REC = _statRecorders.begin();
  
  for(; REC != _statRecorders.end(); ++REC) delete (*REC);
  
  _statRecorders.clear();
  
  reset_stat_table();
}
//-----------------------------------------------------------------------------------------
// reset_stat_table
// ----------------------------------------------------------------------------------------
void StatServices::reset_stat_table ()
{
  for (unsigned int i = 0; i < _statValues.size(); ++i) {
    for (unsigned int j = 0; j < _statValues[i].size(); ++j) {
      delete [] _statValues[i][j];
    }
    _statValues[i].clear();
  }
  _statValues.clear();
  
  _currentStatValues = 0;
}
//-----------------------------------------------------------------------------------------
// recordStats
// ----------------------------------------------------------------------------------------
void StatServices::recordStats (unsigned int gen)
{
  _currentStatValues = new double [ _numRecorders + 2]; //+ replicate + generation
  
  _currentStatValues[0] = SIMenv::getCurrentReplicate();
  _currentStatValues[1] = gen;
  
  unsigned int i = 2, tab_idx = SIMenv::getCurrentReplicate() - 1;
  
  for (rec_it REC = _statRecorders.begin(); REC != _statRecorders.end(); ++REC) 
  {

    assert(i < _numRecorders+2 );
    
    _currentStatValues[i++] = (*REC)->setVal( _popPtr->getCurrentAge() ); 
    
  }
  
  _statValues[ tab_idx ].push_back(_currentStatValues);
  
  assert( _statValues[ tab_idx ].size() <= getNumOccurrences() );
}
//-----------------------------------------------------------------------------------------
// print_headers
// ----------------------------------------------------------------------------------------
void StatServices::printStatHeaders(std::ofstream& FH)
{
  for (rec_it IT = _statRecorders.begin(); IT != _statRecorders.end(); ++IT) {
    FH<<_separator;
    FH.width(_width);
    FH<<(*IT)->getName();
  }
}
//-----------------------------------------------------------------------------------------
// printStatValue
// ----------------------------------------------------------------------------------------
void StatServices::printStatValue(std::ofstream& FH, unsigned int repl_idx)
{
  
  FH.width(_width);
  FH.setf(ios::left,ios::adjustfield);
  FH<<setprecision(_precision);

  if(SIMenv::getCurrentReplicate() == 1) {
        
    //print the first row with stat headers:
    FH<<"replicate"<<_separator;
    FH.width(_width);
    FH<<"generation";
    
    //print the stat names:
    printStatHeaders(FH);
    
    FH<<endl;
  }
    
  for (unsigned int gen = 0, i = 2; gen < _statValues[repl_idx].size(); ++gen) {
    //for all generations recorded, print stuff:

    assert( repl_idx + 1 == _statValues[repl_idx][gen][0]); //check if replicate is ok
    
    FH.width(_width);
    FH<< _statValues[repl_idx][gen][0]; //replicate number
    FH<<_separator;
    FH.width(_width);
    FH<< _statValues[repl_idx][gen][1]; //generation recorded
    
    i = 2;
    //print all stats on one line
    for (rec_it IT = _statRecorders.begin(); IT != _statRecorders.end(); ++IT) {
      
      FH<<_separator;
      FH.width(_width);
      FH<< _statValues[repl_idx][gen][i];
      
      i++;
      
    }
    
    FH<<endl;
    
    assert( i == _numRecorders + 2 );
  }
  
}
//-----------------------------------------------------------------------------------------
// printStatAverage
// ----------------------------------------------------------------------------------------
void StatServices::printStatAverage(std::ofstream& FH)
{
  vector< double > means;
  unsigned int num_repl = SIMenv::getReplicates(), rep_cntr;
  unsigned int num_gen  = getNumOccurrences();
  unsigned int current_gen = 0;
  
  //print the first row with stats headers:
  FH.width(_width);
  FH.setf(ios::left,ios::adjustfield);
  FH<<"generation";
  
  //print the stat names:
  printStatHeaders(FH);
  
  FH<<_separator;
  FH.width(_width);
  FH<<"alive.rpl";
  FH<<endl;
  
  FH<<setprecision(_precision);
  
  for (unsigned int gen = 0; gen < num_gen; ++gen) {
    
    rep_cntr = 0;
    
    means.assign(_numRecorders,0);
    
    for (unsigned int rep = 0; rep < num_repl; ++rep) {
      
      //verify if this replicate was not extinct at this generation
      if ( _statValues[ rep ].size() > gen ) rep_cntr++;
      else continue; //skip this replicate

      for (unsigned int j = 0; j < _numRecorders; ++j) {
        means[j] += _statValues[rep][gen][j+2]; //two first are replicate and generation
      }
      current_gen = _statValues[rep][gen][1];
    }
    
    if(rep_cntr != 0) {
      
      FH.width(_width);
      FH.setf(ios::left,ios::adjustfield);
      //FH<<_statValues[0][gen][1];
      FH<<current_gen;
      
      for (unsigned int j = 0; j < _numRecorders; ++j) {
        FH<<_separator;
        FH.width(_width);
        FH<< means[j] / rep_cntr;
      }
      //print number of alive replicates 'til this generation
      FH<<_separator;
      FH.width(_width);
      FH<<rep_cntr;

      FH<<endl;
    }
    
  }
  
  
}
//-----------------------------------------------------------------------------------------
// getAllStats
// ----------------------------------------------------------------------------------------
list<StatRecBase*> StatServices::getAllStats ( )
{
  list< StatHandlerBase* >::iterator HIT = _statHandlers.begin();
  
  list<StatRecBase*> allStats, stats;
  
  while(HIT != _statHandlers.end()) {
    stats = (*HIT)->getStats();
    allStats.insert(allStats.end(), stats.begin(), stats.end());
    HIT++;
  }
  
  return allStats;
}
