/**  $Id: fileservices.cc,v 1.9.2.1 2014-04-29 18:26:23 fred Exp $
*
*  @file fileservices.cc
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
#include <time.h>
#include <cmath>
#include "simcomponent.h"
#include "fileservices.h"
#include "filehandler.h"
#include "metapop.h"
#include "output.h"

using namespace std;

extern MPIenv *_myenv;

 FileServices::FileServices() : _popPtr(0), _rep_filename(""), _basename(""), _root_dir("") 
{
  _logWriter = new FHLogWriter();
  attach(_logWriter);
}

FileServices::~FileServices (){
  if(_logWriter) delete _logWriter;
}
// ----------------------------------------------------------------------------------------
// attach
// ----------------------------------------------------------------------------------------
void FileServices::attach ( Handler* H )
{
//  Service::attach(FH);
  FileHandler* FH = dynamic_cast<FileHandler*> (H);
  _writers.push_back(FH);
  
  if(FH->get_isInputHandler()) attach_reader(FH);
  
  FH->set_service(this);
}
// ----------------------------------------------------------------------------------------
// attach_reader
// ----------------------------------------------------------------------------------------
void FileServices::attach_reader ( FileHandler* FH )
{
  _readers.push_back(FH);
  FH->set_service(this);
}
// ----------------------------------------------------------------------------------------
// load
// ----------------------------------------------------------------------------------------
void FileServices::load ( SimComponent* sc ) {sc->loadFileServices(this);}
// ----------------------------------------------------------------------------------------
// init
// ----------------------------------------------------------------------------------------
bool FileServices::init (list< ParamSet* >&  params)
{
  char yn;
  bool ok;
  file_it HIT;

  _params = params;
  
  _logWriter->set(0, 0, 0, 0, 0, "");

  //first, build path and filename for each file handler
  HIT = _writers.begin();
  message("    outputs: %s{",_root_dir.c_str());
  while(HIT != _writers.end()) {
    message("%s%s*%s",
            (*HIT)->get_path().c_str(),
            ((*HIT)->get_path().size() !=0? "/":""),
            (*HIT)->get_extension().c_str());
    (*HIT)->init();
    HIT++;
    if(HIT != _writers.end()) message(", ");
  }
  message("}\n");
  //then check if any file already exists
check:

  ok = true;
  HIT = _writers.begin();
  //we have to skip the first writer, that is the log-writer
  while(++HIT != _writers.end()) {
    ok &= (*HIT)->ifExist();

  }

  if(!ok) {

#if defined(_R_OUTPUT_) || defined(LOW_VERBOSE) || defined(USE_MPI)
    yn = 's';
    message("\nPlease choose an other base filename or move the existing file(s)\n");
#else
    //file mode is set during simulation settup, in SimRunner::init()
    if(_mode == 0) {
      message("\nDo you want to overwrite all the files that use it ? (y/n/s(kip)): ");
      cin>>yn;
    } else if( _mode == 1 ) {
      warning("Overwriting existing files.\n");
      yn = 'y';
    } else if( _mode == 2 ) {
      message("\nPlease choose another base filename or move the existing file(s)\n");
      yn = 's';
    } else if( _mode == 3 || _mode == 4) {
      //dryrun, pretend overwriting, silently
      yn = 'y';      
    } else {
      fatal("Run mode not properly set, check simulation parameter \"run_mode\".\n");
    }
#endif

    switch(yn) {
      case 'y':
        break;
      case 's':
        message("Skipping this simulation\n");
        return false;
      default: {
        message("Please give a new output filename : ");
        cin>>_basename;
        goto check;
      }
    }
  }
  
  return true;
}
// ----------------------------------------------------------------------------------------
// notify
// ----------------------------------------------------------------------------------------
void FileServices::notify()
{
  for (file_it file = _writers.begin(); file != _writers.end(); ++file) {
    (*file)->update();
  }
}
// ----------------------------------------------------------------------------------------
// reset
// ----------------------------------------------------------------------------------------
void FileServices::reset()
{
  //Service::reset(); 
  _writers.clear(); _writers.push_back(_logWriter);
  _readers.clear();
}
// ----------------------------------------------------------------------------------------
// log_simparams
// ----------------------------------------------------------------------------------------
void FileServices::log_simparams() 
{
  if(_mode == 4)
    _logWriter->createInitFile(_params);
  else
    _logWriter->save_simparams(_params);
}
// ----------------------------------------------------------------------------------------
// log
// ----------------------------------------------------------------------------------------
void FileServices::log (string message) 
{
  _logWriter->log_message(message);
}
// ----------------------------------------------------------------------------------------
// getFirstReplicateFileName
// ----------------------------------------------------------------------------------------
string& FileServices::getFirstReplicateFileName ()
{
  ostringstream rpl;

  rpl.fill('0');
  rpl.width( (int)log10((double)_popPtr->getReplicates()) + 1);
  rpl<< max( (unsigned)1, _myenv->slaveRank() );

  _rep_filename = _basename + "_" + rpl.str();

  return _rep_filename;
}
// ----------------------------------------------------------------------------------------
// getReplicateFileName
// ----------------------------------------------------------------------------------------
string& FileServices::getReplicateFileName ()
{
  _rep_filename = _basename + "_" + getReplicateCounter();

  return _rep_filename;
}
// ----------------------------------------------------------------------------------------
// getGenerationReplicateFileName
// ----------------------------------------------------------------------------------------
string FileServices::getGenerationReplicateFileName () 
{
  string name = _basename + "_" + getGenerationCounter() + "_" + getReplicateCounter();
  return name;
}
// ----------------------------------------------------------------------------------------
// getBaseFileName
// ----------------------------------------------------------------------------------------
string& FileServices::getBaseFileName () 
{
  return _basename;
}
// ----------------------------------------------------------------------------------------
// getRootDir
// ----------------------------------------------------------------------------------------
std::string& FileServices::getRootDir ()
{
  return _root_dir;
}
// ----------------------------------------------------------------------------------------
// getReplicateCounter
// ----------------------------------------------------------------------------------------
string FileServices::getReplicateCounter ()
{
  ostringstream rpl;
  
  rpl.fill('0');
  rpl.width( (int)log10((double)_popPtr->getReplicates()) + 1);
  rpl<<_popPtr->getCurrentReplicate();
  
  return rpl.str();
}
// ----------------------------------------------------------------------------------------
// getGenerationCounter
// ----------------------------------------------------------------------------------------
string FileServices::getGenerationCounter ()
{
  ostringstream gen;
  
  gen.fill('0');
  gen.width( (int)log10((double)_popPtr->getGenerations()) + 1);
  gen<<_popPtr->getCurrentGeneration();
  
  return gen.str();
}
// ----------------------------------------------------------------------------------------
// setBasename
// ----------------------------------------------------------------------------------------
void FileServices::setBasename (string name) 
{
    _basename = name;
}
// ----------------------------------------------------------------------------------------
// setRootDir
// ----------------------------------------------------------------------------------------
void FileServices::setRootDir (string name) 
{
  _root_dir = name;
 
  if(_root_dir.size() != 0 && _root_dir[_root_dir.length()-1] != '/')
	_root_dir += "/";
  
}
// ----------------------------------------------------------------------------------------
// getReader
// ----------------------------------------------------------------------------------------
FileHandler* FileServices::getReader (string& type)
{
  file_it file = getFirstReader(), last = getLastReader() ;
    
  for(;file != last; file++) {
    if( (*file)->get_extension().compare( type ) == 0) {
      return (*file);
    }
  }
  return 0;
}
