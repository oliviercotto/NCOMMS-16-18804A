/** $Id: binarydatasaver.h,v 1.6.2.2 2014-04-29 16:20:13 fred Exp $
* 
*  @file binarydatasaver.h
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
*  Created on @date 18.10.05
*  
*  @author fred
*/
#ifndef BINARYDATASAVER_H
#define BINARYDATASAVER_H

#include <string>
#include <map>
#include "lifecycleevent.h"
#include "filehandler.h"
#include "binarystoragebuffer.h"

/**A class to handle simulation data saving in binary format.
It is both an LCE and a FileHandler but the writing is executed through the 
LCE interface (i.e. LifeCycleEvent::execute()). Its inheritance from the 
FileHandler class is used to access the basic file services 
(i.e. basename and current replicate filename, etc.).
*/
class BinaryDataSaver: public LifeCycleEvent, public FileHandler {

private:

  BinaryStorageBuffer _buff;
  
  std::string _comp_cmd, _comp_ext, _uncomp_cmd, _tar_cmd, _tar_ext, _dir;
  bool _isPeriodic;
  unsigned int _generation;
  
//  std::map<unsigned int, off_t> _offset_table;

  off_t _offsetDataStart;

  vector< unsigned int > _occurrences;

  vector< unsigned int >::const_iterator _current_occurrence;

  int _fdesc;

  void setFileDescriptor ();
  int  getFileDescriptor () { return _fdesc;}

  void printHeader();
  void storeData();
  void printOffsetTable();
  void finish();
  
public:
  
  BinaryDataSaver ();
  ~BinaryDataSaver () {}

  static pid_t PID;


  void printData();

  //implements LifeCycleEvent:
  virtual bool setParameters ();
  virtual void execute ();
  virtual BinaryDataSaver* clone ( ) {return new BinaryDataSaver();}
  //implements FileHandler:
  virtual void FHwrite (); 
  virtual void FHread (string& filename) {}
  //implements SimComponent 
  virtual void  loadFileServices ( FileServices* loader ) {loader->attach(this);}
  virtual void  loadStatServices ( StatServices* loader ) {}
  virtual age_t removeAgeClass ( ) {return 0;}
  virtual age_t addAgeClass ( ) {return 0;}
  virtual age_t requiredAgeClass () {return 0;}
};

#endif
