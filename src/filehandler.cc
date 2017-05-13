/**  $Id: filehandler.cc,v 1.9.2.2 2014-04-29 18:27:27 fred Exp $
*
*  @file filehandler.cc
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
#include <cstdlib>
#include <cmath>
#include <dirent.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <dirent.h>
#include <sys/stat.h>
#include "output.h"
#include "filehandler.h"
#include "LCEmisc.h"
#include "metapop.h"
#include "version.h"

using namespace std;

extern MPIenv *_myenv;

// ----------------------------------------------------------------------------------------
// init
// ----------------------------------------------------------------------------------------
void FileHandler::init ()
{
  _pop = _service->get_pop_ptr();
  //check if the occurrences exceed the pop parameters:
  if(_GenerationOccurrence > _pop->getGenerations())
    _GenerationOccurrence = _pop->getGenerations();

  if(_ReplicateOccurrence > _pop->getReplicates())
    _ReplicateOccurrence = _pop->getReplicates();

  set_path();
  
//  if ( _myenv->isMaster() != _isMasterExec ) return;

//  message("     %s*%s\n",get_path().c_str(),get_extension().c_str());
}
// ----------------------------------------------------------------------------------------
// ifExist
// ----------------------------------------------------------------------------------------
bool FileHandler::ifExist()
{
  //check if the basefilename is already used on disk:
  string filename = _path + _service->getBaseFileName();

  ostringstream rpl, gen;
    
  gen.fill('0');
  gen.width( (int)log10((double)_pop->getGenerations()) + 1);
  gen<<_GenerationOccurrence;
  
  rpl.fill('0');
  rpl.width( (int)log10((double)_pop->getReplicates()) + 1);
  rpl<<_ReplicateOccurrence;    
  
  if(!_isReplicatePeriodic)
    filename += _extension;
  else if(!_isGenerationPeriodic)
    filename += "_" + rpl.str() + _extension;
  else
    filename += "_" + gen.str() + "_" + rpl.str() + _extension;
  
  ifstream ifXist;
  ifXist.setstate(ios::failbit);
  
  ifXist.open(filename.c_str(),ios::in);
  
  if(ifXist.is_open()) {
    warning("filename \"%s\" used by \"%s\"\n",_service->getBaseFileName().c_str(),filename.c_str());
    ifXist.close();
    return false;
  }
  ifXist.close();

  return true;
}
// ----------------------------------------------------------------------------------------
// set_path
// ----------------------------------------------------------------------------------------
void FileHandler::set_path ( )
{
  string root_dir;
  
  root_dir = _service->getRootDir();
  
  if(_path.size() != 0 && _path[_path.length()-1] != '/')
	_path += "/";
  
  if(root_dir.size() != 0)
    _path = root_dir + _path;
  
  if(_path.size() != 0) {
    
    DIR *dir = opendir(_path.c_str());
    
    //check if we have to create the directory:
    if(_path.size() != 0 && dir == 0) {
      
#ifdef _WINDOWS_
      mode_t mode = S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH;
      if((mkdir(_path.c_str(), mode)) == -1)
#else
        string cmd = "mkdir -p " + _path;
      if ( system( cmd.c_str() ) < 0 )
#endif
      {
        error("could not create directory \"%s\", saving in current directory.\n",_path.c_str());
        _path = "";
      }
    } else {
      if(closedir(dir) == -1)
        warning("FileHandler::set_path::closedir: %s\n",strerror(errno));
    }
  }
}
// ----------------------------------------------------------------------------------------
// get_filename
// ----------------------------------------------------------------------------------------
string& FileHandler::get_filename ()
{
  if(_isReplicatePeriodic)
    if(_isGenerationPeriodic)
      _current_filename = _path + _service->getGenerationReplicateFileName() + _extension;
    else
      _current_filename = _path + _service->getReplicateFileName() + _extension;
  else
    _current_filename = _path + _service->getBaseFileName() + _extension;
  
  return _current_filename;
}
// ----------------------------------------------------------------------------------------
// update
// ----------------------------------------------------------------------------------------
void FileHandler::update()
{
  _current_replicate = _pop->getCurrentReplicate();
  
  _current_generation = _pop->getCurrentGeneration();
  
  // reset the generation iterator to the first in the list for new replicates
  if(_genITER == _generations.end() && _current_generation == 1)
    _genITER = _generations.begin();
  
  if ( ( (_isReplicatePeriodic && !(_current_replicate % _ReplicateOccurrence))
        || (_current_replicate == _ReplicateOccurrence) ) ) {
    if(_isGenerationPeriodic) {
      if(_GenerationOccurrence != 0) {
        if(!(_current_generation % _GenerationOccurrence) &&  _myenv->isMaster() == _isMasterExec) 
          FHwrite();
      } else if(_current_generation == (unsigned)(*_genITER)  &&  _myenv->isMaster() == _isMasterExec) {
        FHwrite();_genITER++;
      }
    } else if(_current_generation == _GenerationOccurrence &&  _myenv->isMaster() == _isMasterExec)
      FHwrite();
  }
}
// ----------------------------------------------------------------------------------------
// FHLogWriter::save_simparams
// ----------------------------------------------------------------------------------------
void FHLogWriter::save_simparams(list< ParamSet* >&  params)
{
  ofstream FH;
  
  open_logfile4writing(FH);
  
  list< ParamSet* >::iterator Pit = params.begin();
  char t_buff[20];
  time_t t = time(NULL);
  strftime(t_buff, 20, "%d-%m-%Y %H:%M:%S", localtime(&t));
  
  if(FH) {
    FH<<"# Nemo    v"<<MAIN_VERSION<<"."<<MINOR_VERSION<<"."<<REVISION<<RELEASE
    <<"  "<<VERSION_DATE<<endl<<"#"<<endl;
    FH<<"# simulation started on "<<t_buff<<endl<<endl;
    
    while(Pit != params.end()) {
      FH<<endl;
      (*Pit)->print(FH);
      Pit++;
    }
    FH<<endl;
    
    FH.close();
  }
}
// ----------------------------------------------------------------------------------------
// FHLogWriter::createInitFile
// ----------------------------------------------------------------------------------------
void FHLogWriter::createInitFile(list< ParamSet* >&  params)
{
  string file = get_service()->getBaseFileName() + ".ini"; // get_filename();
  ofstream FH;
  
  FH.open(file.c_str(), ios_base::out);
  
  list< ParamSet* >::iterator Pit = params.begin();
  char t_buff[20];
  time_t t = time(NULL);
  strftime(t_buff, 20, "%d-%m-%Y %H:%M:%S", localtime(&t));
  
  if(FH) {
    FH<<"# Nemo    v"<<MAIN_VERSION<<"."<<MINOR_VERSION<<"."<<REVISION<<RELEASE
    <<"  "<<VERSION_DATE<<endl<<"#"<<endl;
    FH<<"# init file created on "<<t_buff<<endl<<endl;
    
    while(Pit != params.end()) {
      FH<<endl;
      // !!! we have to change the run mode in the init file !!!
      if ((*Pit)->getName() == "simulation") {
        (*Pit)->get_param("run_mode")->setArg("run");
      }
      
      (*Pit)->print(FH);
      
      Pit++;
    }
    
    FH<<endl;
    FH.close();
  }
}
// ----------------------------------------------------------------------------------------
// FHLogWriter::log_message
// ----------------------------------------------------------------------------------------
void FHLogWriter::log_message (string& logstr)
{
  ofstream FH;
  
  open_logfile4writing(FH, ios_base::app);

  if(FH) {
    FH<<logstr<<endl;
    FH.close();
  }
}
// ----------------------------------------------------------------------------------------
// FHLogWriter::open_logfile4writing
// ----------------------------------------------------------------------------------------
void FHLogWriter::open_logfile4writing (ofstream& FH, ios_base::openmode flag)
{
  string file = get_filename();
  
  FH.open(file.c_str(), flag);
  
  if(!FH) {
    error("FileServices::init: could not output sim parameters to\"%s\"\n",file.c_str());
    FH.close();
  }
}
