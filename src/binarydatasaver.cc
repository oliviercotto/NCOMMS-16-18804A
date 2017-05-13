/** $Id: binarydatasaver.cc,v 1.9.2.4 2016-02-09 13:55:42 fred Exp $
 *
 *  @file binarydatasaver.cc
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
 *  Created on @date 18.10.2005
 *  @author fred
 */

#include <dirent.h>
#include <sys/stat.h>
#include <sys/wait.h>
#include <unistd.h>
#include <cstdlib>
#include <cctype>
#include <cerrno>
#include <string.h>
#include <fcntl.h>
#include <sstream>
#include "metapop.h"
#include "individual.h"
#include "ttrait.h"
#include "binarydatasaver.h"
#include "binarystoragebuffer.cc"
#include "output.h"
#include "version.h"
#include "simenv.h"

extern MPIenv *_myenv;
pid_t BinaryDataSaver::PID = 0;

// ----------------------------------------------------------------------------------------
// BinaryDataSaver::cstor
// ----------------------------------------------------------------------------------------
BinaryDataSaver::BinaryDataSaver()
: LifeCycleEvent("store",""), FileHandler(".bin"), _isPeriodic(0)
{
  ParamUpdater< BinaryDataSaver > * updater =
  new ParamUpdater< BinaryDataSaver > (&BinaryDataSaver::setParameters);

  add_parameter("store_generation",INT,true,false,0,0,updater);
  add_parameter("store_recursive",BOOL,false,false,0,0,updater);
  add_parameter("store_dir",STR,false,false,0,0,updater);
  add_parameter("store_nocompress",BOOL,false,false,0,0,updater);
  add_parameter("store_noarchive",BOOL,false,false,0,0,updater);
  add_parameter("store_compress_cmde",STR,false,false,0,0,updater);
  add_parameter("store_compress_extension",STR,false,false,0,0,updater);
  add_parameter("store_uncompress_cmde",STR,false,false,0,0,updater);
  add_parameter("store_archive_cmde",STR,false,false,0,0,updater);
  add_parameter("store_archive_extension",STR,false,false,0,0,updater);
}
// ----------------------------------------------------------------------------------------
// BinaryDataSaver::setParameters
// ----------------------------------------------------------------------------------------
bool BinaryDataSaver::setParameters()
{
  Param* param;

  _generation = (unsigned int)get_parameter_value("store_generation");

  //make sure we at least save the last generation:
  if(_generation > _popPtr->getGenerations()) _generation = _popPtr->getGenerations();

  param = get_parameter("store_dir");
  if(param->isSet())
    _dir = param->getArg();
  else
    _dir = "";

  param = get_parameter("store_recursive");
  if(param->isSet())
    _isPeriodic = true;
  else
    _isPeriodic = false;

  param = get_parameter("store_compress_cmde");
  if(param->isSet())
    _comp_cmd = param->getArg();
  else
    _comp_cmd = "bzip2"; //default compressor used

  param = get_parameter("store_compress_extension");
  if(param->isSet())
    _comp_ext = param->getArg();
  else
    _comp_ext = ".bz2"; //default

  param = get_parameter("store_uncompress_cmde");
  if(param->isSet())
    _uncomp_cmd = param->getArg();
  else
    _uncomp_cmd = "unbzip2"; //default

  param = get_parameter("store_archive_cmde");
  if(param->isSet())
    _tar_cmd = param->getArg();
  else
    _tar_cmd = "tar"; //default

  param = get_parameter("store_archive_extension");
  if(param->isSet())
    _tar_ext = param->getArg();
  else
    _tar_ext = ".tar"; //default

  param = get_parameter("store_nocompress");
  if(param->isSet()) {_comp_cmd = ""; _comp_ext ="";} //compression process disabled

  param = get_parameter("store_noarchive");
  if(param->isSet()) {_tar_cmd = ""; _tar_ext = "";} //archiving processed disabled

  //set the file handler, will be called by FileServices, will output data at each gen
  FileHandler::set(true, _isPeriodic, 1, _generation, get_rank(), _dir);

  _offsetDataStart = 0L; //this stores the offset of the stored data in the file, positioned after
                     // the simulation parameters stored in plain text.

  _occurrences.clear();

  if(_isPeriodic) {
      unsigned int num_occurrences = _popPtr->getGenerations() / _generation;

      for(unsigned int i = 1; i <= num_occurrences; ++i)
	_occurrences.push_back( i* _generation);
  } else
    _occurrences.push_back(_generation);

  _current_occurrence = _occurrences.begin();


  return true;
}
// ----------------------------------------------------------------------------------------
// BinaryDataSaver::execute
// ----------------------------------------------------------------------------------------
void BinaryDataSaver::execute()
{
  //set the buffer ready to be filled
  if(_popPtr->getCurrentGeneration() == 1) {
      _current_occurrence = _occurrences.begin();
  }

  if( _popPtr->getCurrentGeneration() == *_current_occurrence ) {

      _buff.set_buff(this);
      printHeader();
      storeData();
      printOffsetTable();
      finish();

      _current_occurrence++;
  }

//  //if the generation to store is > than the total nbr of generations,
//  //correct this and store the last one
//  if(!_isPeriodic && _generation == _popPtr->getCurrentGeneration()) {
//    storeData();
//  } else
//	  //if periodic, and generation match, store the data in the buffer
//	  if( _isPeriodic && !(_popPtr->getCurrentGeneration() % _generation)){
//		  storeData();
//	  }
//
//  //if reached last generation, write data to file,
//  //this empties the storage buffer with a call to printData() within FHwrite()
//  //thus prevents further calls to FHwrite to rewrite data to disc
//  if (_popPtr->getCurrentGeneration() == _popPtr->getGenerations()) {
//    //make sure the last generation is always stored before writing data to disc:
//    //catch if tot num gene is not a multiple of _generation
//    if( _popPtr->getCurrentGeneration() % _generation ||
//       //catch if tot num gen is a multiple but _generation < num gen and non periodic
//       (!_isPeriodic && _generation != _popPtr->getGenerations() &&
//        !(_popPtr->getCurrentGeneration() % _generation)))
//    {
//      storeData();
//    }
//    printData(); //this calls _buff.clear()
//    finish();  //this should call printOffsetTable()
//  }
}
// ----------------------------------------------------------------------------------------
// BinaryDataSaver::FHwrite
// ----------------------------------------------------------------------------------------
void BinaryDataSaver::FHwrite()
{
  //need to make sure data is saved even though last generation has not been reached
  //in case of population extinction, the file manager issues a notify() call after
  //the simulation exited from the generations loop, we try to catch this call here

  //keep in mind that this function is always called at the last generation

  //if data is already written to disc, the data buffer is empty
  //if the 'store' LCE is after the 'save_files' LCE in the life cycle, we don't write
  //--> issues two calls to printData, this one would be before storing last generation
  //if 'store' is before 'save_files', the data buffer should be empty by now

  if(_buff.getBuffLength() > 0 && SIMenv::getCurrentRankInLifeCycle() >= LifeCycleEvent::get_rank()) {
    printHeader();
    printOffsetTable(); //this calls _buff.clear()
    finish();
  }

}
// ----------------------------------------------------------------------------------------
// BinaryDataSaver::printHeader
// ----------------------------------------------------------------------------------------
void BinaryDataSaver::printHeader()
{
  //first, record the parameters in text mode:
  ofstream FILE (get_filename().c_str(), ios::out);

  if(!FILE) fatal("could not open Binary output file!!\n");

  list< ParamSet* > current_params = get_service()->get_params();
  list< ParamSet* >::iterator Pit = current_params.begin();

  FILE<<"#NEMO "<<MAIN_VERSION<<" "<<MINOR_VERSION<<" "<<REVISION<<" "<<RELEASE<<" "<<VERSION_DATE<<endl;

#ifdef _DEBUG_
  message("BinaryDataSaver::printHeader:storing parameters\n");
#endif
  while(Pit != current_params.end()) {
    (*Pit)->print(FILE);
    Pit++;
  }
  FILE.close();
}
// ----------------------------------------------------------------------------------------
// BinaryDataSaver::storeData
// ----------------------------------------------------------------------------------------
void BinaryDataSaver::storeData()
{
#ifdef _DEBUG_
  message("BinaryDataSaver::storeData\n");
  unsigned int byte_count = _buff.getBuffLength();
#endif

  unsigned int  generation = _popPtr->getCurrentGeneration();
  unsigned char separator[2] = {'@','G'}; //generation separator = '@G'

  //store the data, begin with generation separator and number:
  _buff.store(&separator, 2 * sizeof(unsigned char));

  _buff.store(&generation, sizeof(unsigned int));

  //open the output file;
  setFileDescriptor();

  //store the position of the generation separator in the buffer into the offset table:
  //_offsetDataStart = _buff.getTotByteRecorded();

  _offsetDataStart = lseek(_fdesc,0,SEEK_END);

//  cout<<"BinaryDataSaver::storeData::offset of start of record: "<<_offsetDataStart
//      <<"; data recorder so far: "<<_buff.getTotByteRecorded()<<"B\n";
  //start storing information, the data buckets will automatically be written to file if >500MB

  //Metapop is a StorableComponent, call store_data:
  //stores all individual info, including trait sequences
  _popPtr->store_data(&_buff);

#ifdef _DEBUG_
  message("BinaryDataSaver::storeData::stored %iB\n",(_buff.getBuffLength()-byte_count));
  byte_count = _buff.getBuffLength();
#endif
}
// ----------------------------------------------------------------------------------------
// BinaryDataSaver::printData
// ----------------------------------------------------------------------------------------
void BinaryDataSaver::setFileDescriptor()
{
  mode_t mode = S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH;
  int flag = O_WRONLY | O_APPEND;

  //open in write and append mode with permission flag set to rw-r--r--
  _fdesc = open(get_filename().c_str(), flag, mode);

  if(_fdesc == -1)
    fatal("BinaryDataSaver::printData::open %s:%s\n",get_filename().c_str(),strerror(errno));

}
// ----------------------------------------------------------------------------------------
// BinaryDataSaver::printData
// ----------------------------------------------------------------------------------------
void BinaryDataSaver::printData()
{
  //we print what is in the buffer and empty it for the next bucket
  if((write(_fdesc,_buff.getBuffer(),_buff.getBuffLength())) == -1)
    fatal("BinaryDataSaver::printData::write %s\n",strerror(errno));

  _buff.emptyBuffer();
}
// ----------------------------------------------------------------------------------------
// BinaryDataSaver::printData
// ----------------------------------------------------------------------------------------
void BinaryDataSaver::printOffsetTable()
{

  //get the offset position of the EOF, will be the nbr of bytes to add to the generation offsets
  // (data has already been written, the header plus any previous data)
  off_t pos = lseek(_fdesc,0,SEEK_END);

  if(pos == -1) fatal("BinaryDataSaver::printOffsetTable::lseek failed %s\n",strerror(errno));

//  cout<<"+++BinaryDataSaver::printOffsetTable::EOF at "<<pos;

  //offset of the offset table:
//  off_t off_table = _buff.getBuffLength() + pos;

//  cout<<"; adding "<<_buff.getBuffLength()<<"B of in-buff data +";

  //add the offset info at the end of the _buff:
  unsigned char separator[3] = {'@','O','T'};
  _buff.store(&separator, 3);

//  std::map<unsigned int, unsigned int>::iterator IT = _offset_table.begin();
//  while(IT != _offset_table.end()) {
//    IT->second += pos;
//    _buff.store((int*)&IT->first, sizeof(int));
//    _buff.store((int*)&IT->second, sizeof(off_t));
//    IT++;
//  }

  //_offset_table = pos; //EOF of current file, position of the start of the generation record

  unsigned int dummy = _popPtr->getCurrentGeneration();

  _buff.store(&dummy, sizeof(int));

  _buff.store(&_offsetDataStart, sizeof(off_t)); //previously stored offset of the generation separator

  //finally, record the number of generations recorded and the offset of start of the offset table
//  int nb_recgen = _offset_table.size();
//  int nb_recgen = 1;  //always store only one generation
//  _buff.store(&nb_recgen, sizeof(int));
//  _buff.store(&off_table, sizeof(off_t));

//  cout<<" "<<4*sizeof(int)+3<<"B for the offset table = "<<_buff.getBuffLength()<<"B\n";

  //now write the data to the file, that is, anything remaining in the buffer:
#ifdef _DEBUG_
  message("BinaryDataSaver::printOffsetTable:writing %iB of data ",_buff.getBuffLength());
#endif

  if((write(_fdesc,_buff.getBuffer(),_buff.getBuffLength())) == -1)
    fatal("BinaryDataSaver::printOffsetTable::write %s\n",strerror(errno));

  if((close(_fdesc)) == -1)
    error("BinaryDataSaver::printOffsetTable::close %s\n",strerror(errno));

  //empty the buffer, get ready for next record
  _buff.clear();

#ifdef _DEBUG_
  message(" [ok]\n");
#endif
}
// ----------------------------------------------------------------------------------------
// BinaryDataSaver::finish
// ----------------------------------------------------------------------------------------
void BinaryDataSaver::finish ()
{
  bool doComp = (_comp_cmd.size() > 0), doTar = (_tar_cmd.size() > 0),
  first = (_popPtr->getCurrentReplicate()
           == max((unsigned)1,_myenv->slaveRank()));

  if (!(doComp || doTar)) return;

  stringstream sysCmd;

  if (doComp)
    sysCmd << _comp_cmd << " " << get_filename() << ";";

  if (doTar) {
    sysCmd << _tar_cmd;
    if (first) sysCmd << " c"; //first replicate, create the tar archive:
    else 	   sysCmd << " r"; //next replicates, append files to it:
    sysCmd << "f " << get_path() << get_service()->getBaseFileName();
    if (!_myenv->isMaster()) sysCmd << _myenv->slaveRank();
    sysCmd << _tar_ext << " --remove-files " << get_filename() << _comp_ext;
  }


  //wait for last child:
//  if (!first) waitpid(PID,NULL,0);
//
//  pid_t tPID = fork();
//
//  if (tPID == 0) {

    if (system(sysCmd.str().c_str()) < 0)
      error("BinaryDataSaver::finish system cmd \"%s\" failed: %s\n",
            sysCmd.str().c_str(),strerror(errno));

//    _exit(EXIT_FAILURE);

//  } else if (tPID < 0) //fork failed :
//    error("BinaryDataSaver::finish::could not fork new process: %s\n",
//          strerror(errno));
//
//  PID = tPID;
}




