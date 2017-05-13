/** $Id: statservices.h,v 1.9.2.2 2014-05-01 15:45:58 fred Exp $
*
*  @file statservices.h
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

#ifndef STATSERVICES_H
#define STATSERVICES_H

#include <list>
#include <string>
#include <fstream>
#include <map>
#include <vector>
#include <assert.h>
#include <cstring> //for memcpy
#include "service.h"
#include "statrecorder.h"

using namespace std;

class Metapop;

class StatHandlerBase;

/**The Service class used to manage the StatHandler objects.*/
class StatServices : public Service {
  
private:
	
  Metapop* _popPtr;
  
  /** List of stat handlers declared by currently active simulation components.*/
  list< StatHandlerBase* > _statHandlers;
  
  /** List of stat recorders. */
  list< StatRecBase *> _statRecorders;
  
  /** Number of stats to record. */
  unsigned int _numRecorders;
  
  /** Table containing all recorded stats, replicate x generation x (num recorders + 2) . */
  vector< vector< double*> > _statValues;
  
  /** Pointer to the last recorded stats.*/
  double* _currentStatValues;
  
  /** The string argument of the 'stat' input option. */
  string _statArg;
  
  /** Deprecated.*/
  unsigned int _occurrence;
  
  /** List of all generations to record.*/
  list< unsigned int > _occurrences;
  
  /** Iterator pointing to the current generation to record. 
      Is incremented once all stats have been recorded for the current generation.
   */
  list< unsigned int >::const_iterator _current_occurrence;
  
  //output format:
  unsigned int _width, _precision;
  
  unsigned char _separator;
  
  bool _printAverages;
  
public:
    
  typedef list< StatHandlerBase* >::const_iterator stat_it;
	
  typedef list< StatRecBase* >::const_iterator rec_it;
  
  StatServices ( ) : _popPtr(0), _occurrence(0), _currentStatValues(0), 
  _numRecorders(0), _width(12), _precision(4), _separator('\t'), _printAverages(1)
  { }
  
  virtual ~StatServices ( ) { reset_stat_table(); }
  
  virtual bool init ( );
  
  Metapop* get_pop_ptr          ( )      {return _popPtr;}
  
  void set_pop_ptr              (Metapop* pop) {_popPtr=pop;}
  
  void set                      (string& str, unsigned int occ) 
  {_statArg = str; _occurrence = occ;}
  
  void setStatOptions           (string& str) {_statArg = str;}
  
  string& getStatArg       ( )  {return _statArg;}
  
  //OCURRENCE
  unsigned int getOccurrence    ( ) {return _occurrence;}
  
  void setOcccurrence           (unsigned int value) {_occurrence = value;}
  /**Sets the list of generation for which statistics must be recorded during a run
     @param timeTable a map containing the generations at which stats must be recorded
   */
  void setOccurrences           (map<unsigned int, unsigned int> timeTable);
  /**Returns the maximum number of generation records per replicate.*/
  unsigned int getNumOccurrences( ) {return _occurrences.size();}
  
  /**Returns the number of generation records present in the stat table for a replicate.
     @param replicate the replicate number (not the index in the stat table).
   */
  unsigned int getNumOccurrences(unsigned int replicate) 
  {
    return _statValues[replicate-1].size();
  }
  /**Returns the last generation recorded for current replicate*/
  unsigned int getCurrentOccurrence () {return *_current_occurrence;}
  /**Resets the occurrence iterator to the beginning of the list of generation occurrences.*/
  void resetCurrentOccurrence () {_current_occurrence = _occurrences.begin();}
  
  //PRINT
  void printStatHeaders         (ofstream& FH);
  void doPrintAverages     ( )  {_printAverages = 1;}
  void cancelPrintAverages ( )  {_printAverages = 0;}
  bool getPrintAveragesOpt ( )  {return _printAverages;}
  
  /**Prints the stat values to the '.txt' output file.
    @param FH the file output stream
    @param repl_idx the replicate index in the stat tables */
  void printStatValue           (ofstream& FH, unsigned int repl_idx);
  
  void printStatAverage         (ofstream& FH);
  
  void setCompactOutputFormat   () {_width = 0; _separator = ' ';}
  
  void setFieldWidth            (unsigned int val) {_width = val;}
  
  void setFieldPrecision        (unsigned int val) {_precision = val;}
  
  void setFieldSeparator        (unsigned char c)  {_separator = c;}
  
  void setDefaultOutputFormat   () {_width = 12; _precision = 4; _separator = '\t';}
  
  //STAT RECORDS
  list<StatRecBase*> getAllStats ( );
  
  unsigned int getNumStats ( ) {return _statRecorders.size();}
  
  stat_it getFirst () {return _statHandlers.begin();}
  
  stat_it getLast () {return _statHandlers.end();}
  
  vector<double*>* getReplicateStatRecords (unsigned int replicate) 
  {
    return &_statValues[ replicate-1 ];
  }
    
  double* getGenerationStatValues (unsigned int replicate, unsigned int occurence) const
  {
    assert(occurence < _occurrences.size());
    return _statValues[ replicate-1 ][ occurence ];
  }
  
  void copyGenerationStatValues (unsigned int replicate, unsigned int occurence, double* values, unsigned int size)
  {
    assert( _statValues[ replicate-1 ].size() < _occurrences.size() );
    
    assert( size == _numRecorders+2);
    
    double* new_record = new double[size];
    
    memcpy(new_record, values, size*sizeof(double));
    
    _statValues[ replicate-1 ].push_back(new_record);
  }
  /**
   * record stat values in stat table by calling all stat recorders.
   */
  void recordStats              (unsigned int gen);
  
  virtual void notify ();
  /** 
	*  tell the SimComponent to load its stat handlers
	*  @param sc the SimComponent
	*/
  virtual void load ( SimComponent* sc );
  /** 
	*  attach the StatHandler to the current list (_statHandlers) of the StatServices
	*  @param H the StatHandler
	*/
  virtual void attach ( Handler* H);
  /** 
	*  clear the list of StatHandler
	*/
  virtual void reset ( );
  
  /** Deletes the stat tables.*/
  void reset_stat_table();
};
#endif //STATSERVICES_H

