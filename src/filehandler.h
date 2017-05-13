/** $Id: filehandler.h,v 1.8 2014-02-07 16:59:05 fred Exp $
*
*  @file filehandler.h
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

#ifndef FILEHANDLER_H
#define FILEHANDLER_H

#include <string>
#include "handler.h"
#include "fileservices.h"
#include "tmatrix.h"

/** Interface to handle file input/output for any SimComponent.
 *  Stores the periodicity parameters and the file path and extension. The replicate file name
 *  is given by the FileServices. A file handler might be set to write output at a specific generation
 *  of a specific replicate or at some periodic time during the simulation.
 *  The default writing mode of a FileHandler is in output and is added as a 'writer' to the File Services
 *  as such using the FileServices::attach mode. The _isInputHandler flag is thus set to FALSE by default. 
 *  To set a FileHandler as a reader instead of a writer, the method FileHandler::set_isInputHandler(TRUE)
 *  must be called prior to attaching it to the FileServices using the FileServices::attach_reader method.
 *  Furthermore, a single instance of a FileHandler can be attached as both a writer and a reader by setting
 *  the _isInputHandler flag to TRUE and then using the FileServices::attach method to add it to the list of
 *  both the reader and writers.
 *  Note on memory management: the FileServices class only deals with its lists of FileHandlers pointers
 *  and will never try to delete any one of those. All memory management must thus be taken care of by the
 *  user class (typically the TraitPrototype...).
 */
class FileHandler : public Handler {
  
private:
  /**Link to the files manager.*/
  FileServices* _service;
  /**Flag telling if the file must be written by the master or the slave node.
    *An example of master-only file is the stats files (.txt and _bygen.txt) that
    *aggregate the stats from each node. Binary files or FSTAT files are written
    *by the slave nodes, each file beeing written for one replicate only.
    *The default behaviour is to be TRUE in the non-MPI version and FALSE in the
    *MPI version.
    */
  bool _isMasterExec;
  /**Writing mode flag. Must be true */
  bool _isInputHandler;
  /**Flag telling if the file should be written each replicate.*/
  bool _isReplicatePeriodic;
  /**Flag setting the per generation periodicity behaviour of the file.*/
  bool _isGenerationPeriodic;
  /**Tells every what replicate the file should be written. Set by the set() method below.*/
  unsigned int _ReplicateOccurrence;
  /**Tells every what generation the file should be written. Set by the set() method below.*/
  unsigned int _GenerationOccurrence;
  /**The current replicate number of the simulation. Set by FileHandler::update().*/
  unsigned int _current_replicate;
  /**The current generation number of the simulation. Set by FileHandler::update().*/
  unsigned int _current_generation;
  
  /**Multiple occurences.*/
  TMatrix _multipleOccurences;

  list< int >::const_iterator _genITER;
  list< int > _generations;
  
  /**unused... yet*/
  unsigned int _ExecRank;
  /**File path as set during initialization. Usually user-defined.*/
  std::string  _path;
  /**File extension, should be specific to the implementation.
    *Is set during construction.*/
  std::string  _extension;
  /**The current filename as set by FileHandler::get_filename().
    *Is composed by putting together the following strings: FileHandler::_path, 
    *FileServices::_basename (+ replicate and generation counters if needed), and 
    *FileHandler::_extension.
    */
  std::string  _current_filename;
  
protected:
    /**Pointer to the current metapop, set during initialization within the init function. */
    Metapop* _pop;
  
public:
	
    FileHandler (const char* ext) 
    : _service(0), 
#ifdef USE_MPI
    _isMasterExec(false), 
#else
    _isMasterExec(true), 
#endif
    _isInputHandler(0),
    _isReplicatePeriodic(0), _isGenerationPeriodic(0),
    _ReplicateOccurrence(0), _GenerationOccurrence(0),
    _current_replicate(0), _current_generation(0), _ExecRank(0),
    _path(), _extension(ext), _current_filename(), _pop(0)
    {}
  
  virtual ~FileHandler ( ) { }
  /**Called by notifier during simulation setup, performs file checking.*/
  virtual void  init ( );
  /**Checks if any file associated with the current file name already exists on disk. 
   * Only checks wether the first replicate file is present. 
   *@return false if filename already exists on disk, true otherwise.*/
  virtual bool ifExist( );
  ///@name Accessors
  ///@{
  /**Returns the pointer to the current metapop through the FileServices interface.*/
  Metapop*      get_pop_ptr ( )      {return _pop;}
  
  void          set_pop_ptr (Metapop* pop_ptr) { _pop = pop_ptr;}
  /**Returns pointer to the FileServices.*/
  FileServices* get_service () {return _service;}
  
  void          set_service (FileServices* srv) {_service = srv;}

  std::string&  get_path () {return _path;}
  
  void          set_path ( );
  
  std::string&  get_extension ( )  {return _extension;}
  
  void          set_extension (const char* ext) {_extension = ext;}
  /**Builds and returns the current file name depending on the periodicity of the file.*/
  std::string&  get_filename ();
  
  
  bool get_isInputHandler () {return _isInputHandler;}
  void set_isInputHandler (bool val) {_isInputHandler = val;}
  
  bool get_isReplicatePeriodic ()         {return _isReplicatePeriodic;}
  void set_isReplicatePeriodic (bool val) {_isReplicatePeriodic = val;}
  
  unsigned int  get_ReplicateOccurrence  ()         {return _ReplicateOccurrence;}
  void set_ReplicateOccurrence  (unsigned int val)  {_ReplicateOccurrence = val;}
  
  bool get_isGenerationPeriodic ()         {return _isGenerationPeriodic;}
  void set_isGenerationPeriodic (bool val) {_isGenerationPeriodic = val;}
  
  unsigned int  get_GenerationOccurrence  ()         {return _GenerationOccurrence;}
  void set_GenerationOccurrence  (unsigned int val)  {_GenerationOccurrence = val;}
  /**unused yet...*/
  unsigned int  get_ExecRank ()        {return _ExecRank;}
  void set_ExecRank (int val) {_ExecRank = val;}
  
  TMatrix* get_OccMatrix () {return &_multipleOccurences;}
  
  void set_OccMatrix (TMatrix* occ) {
    _multipleOccurences.reset(occ->nrows(), occ->ncols(), occ->getValArray());
    _generations.clear();
    for(unsigned int i = 0; i < _multipleOccurences.ncols(); i++)
      _generations.push_back((int)_multipleOccurences.get(0, i));
    _genITER = _generations.begin();
  }
  
  bool get_isMasterExec() {return _isMasterExec;}
#ifdef USE_MPI
  void set_isMasterExec(bool is) {_isMasterExec = is;}
#else
  void set_isMasterExec(bool is) {_isMasterExec = true;}
#endif
  ///@}
  /**Sets the hanlder parameters.
    @param rpl_per replicate periodicity
    @param gen_per generation periodicity
    @param rpl_occ replicate occurence
    @param gen_occ generation occurence
    @param rank the rank in the life cycle, actualy unused...
    @param path the file path
    */
  virtual void set (bool rpl_per, bool gen_per, int rpl_occ, int gen_occ, int rank, string path) {
	set_isReplicatePeriodic(rpl_per); set_isGenerationPeriodic(gen_per); set_ReplicateOccurrence(rpl_occ);
	set_GenerationOccurrence(gen_occ); set_ExecRank(rank); _path = path;}
  
  virtual void set_multi (bool rpl_per, bool gen_per, int rpl_occ, TMatrix* Occ, string path) {
    set_isReplicatePeriodic(rpl_per); set_isGenerationPeriodic(gen_per); set_ReplicateOccurrence(rpl_occ);
    set_GenerationOccurrence(0); set_OccMatrix(Occ); _path = path;}
  
  /**Default behavior of the class, called by Handler::update().**/
  virtual void  FHwrite () = 0;
  
  /**Default input function. Loads a pop from the genotypes read from the input file.*/
  virtual void  FHread (string& filename) = 0;
  /**Updates the inner replicate and generation counters and calls FHwrite if needed by the
    *the periodicity of the file.*/
  virtual void update();
  
};

//CLASS TraitFileHandler
//
/**Template class for the trait's FileHandler. Constructor links to a given trait prototype.
The pointer and index of the trait prototype are stored in the class. They
can be accessd through the _FHLinkedTrait and _FHLinkedTraitIndex members to get the stats. */
template <class TP> class TraitFileHandler : public FileHandler {
protected:
  TP* _FHLinkedTrait;
  int _FHLinkedTraitIndex;
public:
  TraitFileHandler(TP* trait_proto, const char* ext);
  virtual ~TraitFileHandler () {}
  
  virtual void  FHwrite () = 0;
    
  virtual void  FHread (string& filename) = 0;

  virtual void set (bool rpl_per, bool gen_per, int rpl_occ, int gen_occ, int rank, string path, TP* trait_proto);
};

template <class TP> TraitFileHandler<TP>::TraitFileHandler(TP* trait_proto, const char* ext) :
FileHandler(ext), _FHLinkedTrait(trait_proto), _FHLinkedTraitIndex(trait_proto->get_index())
{ }

template <class TP> void TraitFileHandler<TP>::set(bool rpl_per, bool gen_per, int rpl_occ, int gen_occ, int rank, string path, TP* trait_proto)
{
  FileHandler::set(rpl_per, gen_per, rpl_occ, gen_occ,rank, path);
  _FHLinkedTrait = trait_proto;
  _FHLinkedTraitIndex = trait_proto->get_index();
}

// EventFileHandler
//
/**Template class for the LCEs StatHandler classes. Constructor links to a give LCE.
The LCE can be accessed through the _SHLinkedEvent member to get the stats.*/
template <class LCE> class EventFileHandler : public FileHandler {
protected:
  LCE* _FHLinkedEvent;
public:
    EventFileHandler(LCE* event, const char* ext);
  virtual ~EventFileHandler () {}
  
  virtual void  FHwrite () = 0;
  
  virtual void FHread (string& filename) = 0;

  virtual void set (bool rpl_per, bool gen_per, int rpl_occ, int gen_occ, int rank, string path, LCE* event);
};

template <class LCE> EventFileHandler<LCE>::EventFileHandler(LCE* event, const char* ext) :
FileHandler(ext), _FHLinkedEvent(event)
{ }

template <class LCE> void EventFileHandler<LCE>::set(bool rpl_per, bool gen_per, int rpl_occ, int gen_occ, int rank, string path, LCE* event)
{
  FileHandler::set(rpl_per,gen_per,rpl_occ,gen_occ,rank,path);
  _FHLinkedEvent = event;
}

/**File Handler used to save the simulation parameters to a log file.*/
class FHLogWriter : public FileHandler {
public:
  FHLogWriter() : FileHandler(".log"){};
  virtual ~FHLogWriter(){}
  
  virtual void FHwrite (){}//@TODO implement log writer procedures
  virtual void FHread  (string& filename) {}
  
  
  void createInitFile(list< ParamSet* >&  params);
  
  void save_simparams(list< ParamSet* >&  params);

  void log_message (string& logstr);
  
  void open_logfile4writing (ofstream& FH, ios_base::openmode flag = ios_base::out);
};
#endif //FILEHANDLER_H

