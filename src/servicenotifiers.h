/** $Id: servicenotifiers.h,v 1.9.2.2 2014-04-29 18:19:11 fred Exp $
*
*  @file servicenotifiers.h
*  Nemo2
*
*  Copyright (C) 2008-2011 Frederic Guillaume
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
*  Created on @date 31.12.2008
*  @author fred
*/

#include "types.h"
#include "lifecycleevent.h"
#include "fileservices.h"
#include "statservices.h"
#include "updaterservices.h"
#include "filehandler.h"
#include "stathandler.h"


//LCE_ParamUpdaterNotifier
//
/**Calls the UpdaterServices to notify its components of a generation change.*/
class LCE_ParamUpdaterNotifier : public LifeCycleEvent {
  
private:
  UpdaterServices* _manager;
  
public:
  
  LCE_ParamUpdaterNotifier () : LifeCycleEvent("param_updater",""), _manager(0) { }
  virtual ~LCE_ParamUpdaterNotifier() {}
  
  void setManager (UpdaterServices* mng) {_manager = mng;}
  virtual void execute ();
  virtual bool setParameters () {return true;}
  
  virtual LifeCycleEvent* clone ( ) {return new LCE_ParamUpdaterNotifier();}
  
  //SimComponent overrides:
  virtual void  loadFileServices ( FileServices* loader ) {}
  virtual void  loadStatServices ( StatServices* loader ) {}
  virtual void  loadUpdaters ( UpdaterServices* loader) {_manager = loader;}
  virtual age_t removeAgeClass ( ) {return 0;}
  virtual age_t addAgeClass ( ) {return 0;}
  virtual age_t requiredAgeClass () {return 0;}
};

// LCE_FileServicesNotifier
//
/**Event used to notify all file handlers to update their state through the FileServices::notify() interface.*/
class LCE_FileServicesNotifier: public LifeCycleEvent
  {  
    
    FileServices* _service;
    
  public:
    LCE_FileServicesNotifier() 
    : LifeCycleEvent("save_files",""), _service(0) { }
    virtual ~LCE_FileServicesNotifier( ) { }
    
    virtual void execute ();
    virtual bool setParameters () {return true;}
    
    virtual LifeCycleEvent* clone ( ) {return new LCE_FileServicesNotifier();}
    
    //SimComponent overrides:
    virtual void  loadFileServices ( FileServices* loader ) {_service = loader;}
    virtual void  loadStatServices ( StatServices* loader ) {}
    virtual age_t removeAgeClass ( ) {return 0;}
    virtual age_t addAgeClass ( ) {return 0;}
    virtual age_t requiredAgeClass () {return 0;}
  };

// LCE_StatFH
//
/**FileHandler of the LCE_StatServiceNotifier class, writes the recorded stats to txt files.*/
class LCE_StatFH : public FileHandler {
  
  StatServices* _statService;
  
public:
  
  LCE_StatFH () 
  : FileHandler(".txt"), _statService(0) 
  {
    FileHandler::set(false,false,0,0,0,"");
    set_isMasterExec(true);
  }
  
  ~LCE_StatFH() { };
  
  virtual bool ifExist();
  
  void set_statService(StatServices* srv) {_statService = srv;}
  
  virtual void FHwrite();
  virtual void FHread (string& filename) {}
  
  void PrintStat_byGen ( );
};


// LCE_StatServiceNotifier
//
/**Initiates the StatServices' parameters (log time) when registering, calls StatServices::notify() when executing.
 * Registers the file handler used to save the stats to the '.txt' and '_bygen.txt' files.*/
class LCE_StatServiceNotifier: public LifeCycleEvent
  { 
    
    StatServices* _service;
    
    unsigned int _occurrence;
    
    string _arg, _dir;
    
    LCE_StatFH _fileHandler;
    
  public:
    
    LCE_StatServiceNotifier ( );
    
    virtual ~LCE_StatServiceNotifier ( ) { }
    
    FileHandler& getFH() {return _fileHandler;}
    bool setOccurence ();
    bool dummyUpdate  (){return true;}
    
    virtual bool setParameters ();
    virtual void  execute ();
    
    virtual LCE_StatServiceNotifier* clone ( ) {return new LCE_StatServiceNotifier();}
    
    //SimComponent overrides:
    virtual void loadFileServices ( FileServices* loader ) {loader->attach(&_fileHandler);}
    virtual void loadStatServices ( StatServices* loader );
    virtual age_t removeAgeClass ( ) {return 0;}
    virtual age_t addAgeClass ( ) {return 0;}
    virtual age_t requiredAgeClass () {return 0;}
    
  };
