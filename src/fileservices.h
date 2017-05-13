/** $Id: fileservices.h,v 1.8 2014-02-07 16:59:05 fred Exp $
*
*  @file fileservices.h
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



#ifndef FILESERVICES_H
#define FILESERVICES_H

#include <list>
#include <string>
#include <map>
#include "service.h"
#include "param.h"

class Metapop;

class FileHandler;
class FHLogWriter;
/**A class to manage the files associated with each components of the simulation.
   
   Implements the Observer design pattern (is the concrete subject), stores the base filename of the simulation and 
   updates the replicate filenames. It also performs files checking and saves the simulation parameters on init.
*/

class FileServices : public Service {

private:
  /**a pointer to the current Metapop*/
  Metapop*    _popPtr;
  /**a FileHandler used to save the simulation parameters on disk.*/
  FHLogWriter* _logWriter;
  
  /**the list of the FileHandler's registered by the SimComponent in output mode*/
  list< FileHandler* > _writers;
  
  /**the list of the FileHandler's registered by the SimComponent in input mode*/
  list< FileHandler* > _readers;  
  
  /**the file name associated with the current simulation replicate*/
  string _rep_filename;
  
  /**the base file name of the simulation, read from the init file (param "filename") */
  string _basename;
  
  /**the root directory for the simulation's results, read from the init file (param "root_dir") */
  string _root_dir;
  
  /**the list of the parameters of the current simulation.
   The list is created when init() is called from SimRunner::run().*/
  list< ParamSet* > _params;
  
  /**File mode, sets behavior when file must be overwritten or not.*/
  unsigned int _mode;
  
public:
	
  typedef std::list< FileHandler* >::const_iterator file_it;
    
  FileServices ( );
  
  virtual ~FileServices ( );
  
  virtual bool init ( ) {return false;}
 
  virtual void notify ();
  /**Checks if files with _basename already exist and save the simulation parameters in log files.
   * @param params a ref to the list of the current parameters of the simulation
   * @return true if the files check is ok
   * @return false if the user wants to skip this simulation
   */
  bool init (list< ParamSet* >&  params);
  
  /**Accessor to the pointer to the main population.
   * @return the pointer to the current Metapop as set during simulation setup.*/
  virtual Metapop*   get_pop_ptr ( )      {return _popPtr;}
  
  /**Sets the Metapop reference.*/
  virtual void set_pop_ptr (Metapop* pop) {_popPtr = pop;}
  
  /**Mode setter, determines if file will get overwritten or not.*/
  void setMode(unsigned int m) {_mode = m;}
  
  /**Writting mode getter.*/
  unsigned int getMode() {return _mode;}
  
  /**Sets the base file name of the simulation*/
  void setBasename (string name);
  
  /**Sets the root directory of the simulation*/
  void setRootDir (string name);
  
  /**Saves the current simulation parameters to the default parameter logfile. */
  void log_simparams();
  
  /**Write to the parameter logfile.*/
  void log (string message);
  
  /**Accessor to the list of the current parameters of the simulation.*/
  list< ParamSet* >& get_params() {return _params;};
  
  /**Accessor to first element of the list of output FileHandlers.*/
  file_it getFirstWriter() {return _writers.begin();}
  
  /**Accessor to last element of the list of output FileHandlers.*/
  file_it getLastWriter () {return _writers.end();}  
  
  /**Accessor to first element of the list of input FileHandlers.*/
  file_it getFirstReader() {return _readers.begin();}
  
  /**Accessor to last element of the list of input FileHandlers.*/
  file_it getLastReader () {return _readers.end();}
  
  /**Accessor to a specific file handler specified by its extension string.*/
  FileHandler* getReader (string& type);
  
  /**Accessor to the base file name of the simulation.*/
  string&  getBaseFileName ();
  
  /**Accessor to the name of the simulation's root output directory.   */
  string&  getRootDir ();
  
  /**Accessor to the current replicate counter string.*/
  string  getReplicateCounter ();
  
  /**Accessor to the first replicate file name.*/
  string&  getFirstReplicateFileName ();
  
  /**Accessor to the current replicate file name.*/
  string&  getReplicateFileName ();
  
  /**Accessor to the current generation counter string.*/
  string  getGenerationCounter ();
  
  /**Accessor to the current file name with generation and replicate counters added.*/
  string  getGenerationReplicateFileName ();
  
  /**Tells the SimComponent to load its file handlers.
   *  @param sc the SimComponent*/
  virtual void load ( SimComponent* sc );
  
  /**Attaches the FileHandler to the current list (_writers) of the FileServices.
   *  @param FH the FileHandler*/
  virtual void attach ( Handler* FH );
  
  /**Attaches the FileHandler to the current list (_readers) of the FileServices.
    *  @param FH the FileHandler*/
  virtual void attach_reader ( FileHandler* FH );
  
  /**Clears the list of FileHandlers. */
  virtual void reset ( );

};
#endif //FILESERVICES_H

