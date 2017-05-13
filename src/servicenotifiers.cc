/** $Id: servicenotifiers.cc,v 1.9 2014-02-07 11:01:59 fred Exp $
*
*  @file servicenotifiers.cc
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

#include <iostream>
#include <iomanip>
#include "servicenotifiers.h"
#include "metapop.h"
#include "Uniform.h"
#include "output.h"
#include "simenv.h"
#include "tstring.h"

/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

//                         ******** LCE_ParamUpdaterNotifier ********/

// ----------------------------------------------------------------------------------------
// LCE_ParamUpdaterNotifier::execute
// ----------------------------------------------------------------------------------------
void LCE_ParamUpdaterNotifier::execute ()
{
#ifdef _DEBUG_
  message("LCE_ParamUpdaterNotifier::execute (gen: %i rpl: %i)\n",
		  this->get_pop_ptr()->getCurrentGeneration(),
		  this->get_pop_ptr()->getCurrentReplicate());
#endif
  _manager->notify( _popPtr->getCurrentGeneration() );
}

/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

//                         ******** LCE_FileServicesNotifier ********/

// ----------------------------------------------------------------------------------------
// LCE_FileServicesNotifier::execute
// ----------------------------------------------------------------------------------------
void LCE_FileServicesNotifier::execute ()
{
#ifdef _DEBUG_
  message("LCE_FileServicesNotifier::execute (gen: %i rpl: %i)\n",
		  this->get_pop_ptr()->getCurrentGeneration(),
		  this->get_pop_ptr()->getCurrentReplicate());
#endif
  _service->notify();
}

/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

//                         ******** LCE_StatServiceNotifier ********/

// ----------------------------------------------------------------------------------------
// LCE_StatServiceNotifier
// ----------------------------------------------------------------------------------------
LCE_StatServiceNotifier::LCE_StatServiceNotifier ()
: LifeCycleEvent("save_stats",""), _service(0)
{
  ParamUpdater<LCE_StatServiceNotifier> * updater;

  updater = new ParamUpdater<LCE_StatServiceNotifier> (&LCE_StatServiceNotifier::dummyUpdate);
  
  add_parameter("stat",STR,true,false,0,0);
  
  add_parameter("stat_log_time",INT,true,false,0,0,updater);
  
  add_parameter("stat_dir",STR,false,false,0,0);
  
  add_parameter("stat_output_compact",BOOL,false,false,0,0);
  
  add_parameter("stat_output_CSV",BOOL,false,false,0,0);
  
  add_parameter("stat_output_width",INT,false,false,0,0);
  
  add_parameter("stat_output_precision",INT,false,false,0,0);
  
  add_parameter("stat_output_no_means",BOOL,false,false,0,0);
  
  _service = SIMenv::MainSim->get_StatServices();
}
//-----------------------------------------------------------------------------
// setParameters
//-----------------------------------------------------------------------------
bool LCE_StatServiceNotifier::setParameters ( )
{  
  _arg = _paramSet->getArg("stat");

  setOccurence();
  
  _dir = _paramSet->getArg("stat_dir");
  
  _fileHandler.set(true, false, 1,
                   SIMenv::getGenerations(), //this is when data is written
                   this->get_rank(),
                   _dir);
  
  //output format:  
  bool is_compact = false;
  if (_paramSet->isSet("stat_output_compact")) {

    _service->setCompactOutputFormat();
    is_compact = true;
    
  } else if (_paramSet->isSet("stat_output_CSV")) {
    
    _service->setCompactOutputFormat();
    _service->setFieldSeparator(',');
    is_compact = true;

  } else {
    
    _service->setDefaultOutputFormat();
    
  }

  if (_paramSet->isSet("stat_output_width") && !is_compact) {
    _service->setFieldWidth((unsigned int)_paramSet->getValue("stat_output_width"));
  }
  
  if (_paramSet->isSet("stat_output_precision")) {
    _service->setFieldPrecision((unsigned int)_paramSet->getValue("stat_output_precision"));
  }
  return true;
}
//-----------------------------------------------------------------------------
// setOccurence
//-----------------------------------------------------------------------------
bool LCE_StatServiceNotifier::setOccurence ( )
{
  map<unsigned int, unsigned int> atTimes;
  
  Param* logtime = _paramSet->get_param("stat_log_time");
  
  
  if (!(logtime->isMatrix() || logtime->isTemporal())) {
  
    _occurrence = (unsigned int)logtime->getValue();
    
    atTimes[0] = _occurrence;

    
  } else if (logtime->isMatrix()) {
    
    
    //a matrix argument may be used to specify a set of generations
    //where results will be saved, it is not recursive
    
    TMatrix tmp;
    logtime->getMatrix(&tmp);
    
    if (tmp.nrows() > 1) {
      fatal("parameter \"stat_log_time\" only accepts one-dimensional arrays.\n");
      return false;
    }
    
    for (unsigned int i = 0; i < tmp.ncols(); ++i)
      atTimes[ tmp.get(0, i) ] = 0; //zero is used to indicate no recursivity
    
    
  } else if (logtime->isTemporal()) {
    
    
    //the occurence time changes during the course of a replicate
    //get the temporal args:
    
    deque< unsigned int > tempoDates = logtime->getUpdatingDates();
    deque< string > tempoArgs = logtime->getTemporalArgs();
    
    for (unsigned int i = 0; i < tempoDates.size(); i++)
      atTimes[ tempoDates[i] ] = tstring::str2int( tempoArgs[i] );
    
  } 
    
  _service->setOccurrences(atTimes);
  
  return true;
}
//-----------------------------------------------------------------------------
// loadStatServices
//-----------------------------------------------------------------------------
void LCE_StatServiceNotifier::loadStatServices ( StatServices* loader )
{
//  _service = loader; //already set in cstor
  
  _fileHandler.set_statService(loader);
  
  if(!get_paramset()->isSet())
    error("LCE_StatServiceNotifier::loadStatServices:stat params are not set!\n");

  loader->setStatOptions(_arg); //occurence has been set previously
  
  if ( _paramSet->isSet("stat_output_no_means") )
    loader->cancelPrintAverages();
  else
    loader->doPrintAverages();

}
//-----------------------------------------------------------------------------
// LCE_StatServiceNotifier::execute
//-----------------------------------------------------------------------------
void LCE_StatServiceNotifier::execute ()
{
#ifdef _DEBUG_
  message("LCE_StatServiceNotifier::execute (occurrence %i)\n",_service->getOccurrence());
#endif
  _service->notify();
}

/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/
//
//                               ******** LCE_StatFH ********/
//
// ----------------------------------------------------------------------------------------
// LCE_StatFH::ifExist
// ----------------------------------------------------------------------------------------
bool LCE_StatFH::ifExist()
{  
  bool status = true;
  ifstream ifExist1, ifExist2;
  ifExist1.setstate(ios::failbit);
  ifExist2.setstate(ios::failbit);
  
  //check if the basefilename is already used on disk:
  string filename = get_path() + get_service()->getBaseFileName() + ".txt";
  
  ifExist1.open(filename.c_str(),ios::in);
  if(ifExist1.is_open()) {
    warning("filename \"%s\" used by \"%s\"\n",get_service()->getBaseFileName().c_str(),filename.c_str());
    status = false;
  }
  ifExist1.close();
  
  filename = get_path() + get_service()->getBaseFileName() + "_bygen.txt";
  
  ifExist2.open(filename.c_str(),ios::in);
  if(ifExist2.is_open()) {
    warning("filename \"%s\" used by \"%s\"\n",get_service()->getBaseFileName().c_str(),filename.c_str());
    status = false;
  }
  ifExist2.close();
  
  return status;
}
// ----------------------------------------------------------------------------------------
// LCE_StatFH::FHwrite
// ----------------------------------------------------------------------------------------
void LCE_StatFH::FHwrite()
{
  ofstream FH;
  
  string filename = get_path() + get_service()->getBaseFileName() + get_extension();
  
#ifdef _DEBUG_
  message("LCE_StatFH::FHwrite (%s)\n",filename.c_str());
#endif
  
  if(SIMenv::getCurrentReplicate() == 1) {
    
    FH.open(filename.c_str(),ios::trunc);    
    
    if(!FH) fatal("LCE_StatFH::FHwrite:: could not open stat output file \"%s\"\n",filename.c_str());

  } else {
    
    FH.open(filename.c_str(),ios::app);
    
    if(!FH) fatal("LCE_StatFH::FHwrite:: could not open stat output file \"%s\"\n",filename.c_str());

  }
  
  _statService->printStatValue(FH, SIMenv::getCurrentReplicate() - 1);
  
  FH.close();
  
  //write the stat averages to the 'bygen.txt' file if at last replicate.
  if(SIMenv::getCurrentReplicate() == SIMenv::getReplicates() 
     && SIMenv::getReplicates() > 1
     && _statService->getPrintAveragesOpt())
    PrintStat_byGen();
}
// ----------------------------------------------------------------------------------------
// LCE_StatFH::PrintStat_byGen
// ----------------------------------------------------------------------------------------
void LCE_StatFH::PrintStat_byGen()
{
  ofstream FH;
  
  string filename = get_path() + get_service()->getBaseFileName() + "_bygen.txt";
  
  FH.open(filename.c_str(),ios::trunc);
  
  if(!FH) fatal("PrintStat:: could not open stat output file\n");
  
  _statService->printStatAverage(FH);
  
  FH.close();
}
