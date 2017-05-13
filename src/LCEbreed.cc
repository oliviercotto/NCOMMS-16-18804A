/**  $Id: LCEbreed.cc,v 1.9.2.6 2015-03-16 20:17:01 fred Exp $
 *
 *  @file LCEbreed.cc
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
 *  created on @date 07.07.2004
 *
 *  @author fred
 */
#include <deque>
#include <algorithm>
#include "output.h"
#include "LCEbreed.h"
#include "individual.h"
#include "metapop.h"
#include "Uniform.h"
#include "simenv.h"

/* _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/ */

//                             ******** LCE_Breed_base ********/

// ----------------------------------------------------------------------------------------
// LCE_Breed_base
// ----------------------------------------------------------------------------------------
LCE_Breed_base::LCE_Breed_base ()
: LifeCycleEvent("", ""),
_mating_system(0), 
_mating_proportion(1), 
_mean_fecundity(0),
_LeslieMatrix(0),
//_growthRates(0),
MatingFuncPtr(0),
DoBreedFuncPtr(0), 
FecundityFuncPtr(0), 
CheckMatingConditionFuncPtr(0), 
GetOffsprgSex(0)
//,GetPatchFecundityFuncPtr(0)
{
  ParamUpdater<LCE_Breed_base> * updater = new ParamUpdater<LCE_Breed_base> (&LCE_Breed_base::setMatingSystem);
  
  add_parameter("mating_system",INT,true,true,1,6, updater);
  add_parameter("mating_proportion",DBL,false,true,0,1,updater);
  add_parameter("mating_males",INT,false,false,0,0, updater);
  //  add_parameter("growth_model", INT, false, true, 1, 7, updater);
  //  add_parameter("growth_rate", DBL, false, false, 0, 0, updater);
  
  updater = new ParamUpdater< LCE_Breed_base > (&LCE_Breed_base::setFecundity);
  add_parameter("mean_fecundity",DBL,false,false,0,0, updater);
  add_parameter("fecundity_dist_stdev",DBL,false,false,0,0, updater);
  add_parameter("fecundity_distribution",STR,false,false,0,0, updater);
  
  updater = new ParamUpdater< LCE_Breed_base > (&LCE_Breed_base::setSexRatio);
  add_parameter("sex_ratio_mode",STR,false,false,0,0, updater);
}
// ----------------------------------------------------------------------------------------
// LCE_Breed_base::setParameters
// ----------------------------------------------------------------------------------------
bool LCE_Breed_base::setParameters ()
{
  return ( setMatingSystem() && setFecundity() && setSexRatio() );
}
// ----------------------------------------------------------------------------------------
// LCE_Breed_base::setMatingSystem
// ----------------------------------------------------------------------------------------
bool LCE_Breed_base::setMatingSystem ()
{
  _mating_system = (int)this->get_parameter_value("mating_system");
  
  if(get_parameter("mating_proportion")->isSet())
    _mating_proportion = this->get_parameter_value("mating_proportion");
  else
    _mating_proportion = 1;
  
  if(_paramSet->isSet("mating_males"))
    _mating_males = (int)_paramSet->getValue("mating_males");
  else
    _mating_males = 1;
  
  //set the mating functions ptr:
  CheckMatingConditionFuncPtr = &LCE_Breed_base::checkNoSelfing;
  DoBreedFuncPtr = &LCE_Breed_base::breed;
  _do_inherit = true; //is true unless the mating system is cloning
  
  switch(_mating_system) {
      //random mating:
    case 1:
    {
      MatingFuncPtr = &LCE_Breed_base::RandomMating;
      break;
    }
      //polygyny:
    case 2:
    {
      if(_mating_proportion == 1)
        if(_mating_males == 1)
          MatingFuncPtr = &LCE_Breed_base::fullPolyginy;
        else
          MatingFuncPtr = &LCE_Breed_base::fullPolyginy_manyMales;
        else
          if(_mating_males == 1)
            MatingFuncPtr = &LCE_Breed_base::partialPolyginy;
          else
            MatingFuncPtr = &LCE_Breed_base::partialPolyginy_manyMales;
      
      CheckMatingConditionFuncPtr = &LCE_Breed_base::checkPolygyny;
      break;
    }
      //monogamy:
    case 3:
    {
      if(_mating_proportion == 1)
        MatingFuncPtr = &LCE_Breed_base::fullMonoginy;
      else 
        MatingFuncPtr = &LCE_Breed_base::partialMonoginy;
      
      break;
    }
      //selfing:
    case 4:
    {
      if(_mating_proportion == 1)
        MatingFuncPtr = &LCE_Breed_base::fullSelfing;
      else
        MatingFuncPtr = &LCE_Breed_base::partialSelfing;
      
      CheckMatingConditionFuncPtr = &LCE_Breed_base::checkSelfing;
      break;
    }
      //cloning
    case 5:
    {
      if(_mating_proportion == 1)
        MatingFuncPtr = &LCE_Breed_base::fullSelfing;
      else
        MatingFuncPtr = &LCE_Breed_base::partialSelfing;
      
      CheckMatingConditionFuncPtr = &LCE_Breed_base::checkCloning;
      DoBreedFuncPtr = &LCE_Breed_base::breed_cloning;
      //      _do_inherit = false;
      break;
    }
      //random mating in Wright-Fisher model with hermaphrodites:
    case 6:
    {
      MatingFuncPtr = &LCE_Breed_base::random_hermaphrodite;
      CheckMatingConditionFuncPtr = &LCE_Breed_base::checkSelfing;
      break;
    }
      
  }
  
  //Growth model
  //  unsigned int model;
  //  if(_paramSet->isSet("growth_model"))
  //    model = get_parameter_value("growth_model");
  //  else 
  //    model = 1;
  //  
  //  switch (model) {
  //    case 1:
  //      GetPatchFecundityFuncPtr = &LCE_Breed_base::instantGrowth;
  //      break;
  //    case 2:
  //      GetPatchFecundityFuncPtr = &LCE_Breed_base::logisticGrowth;
  //      break;
  //    case 3:
  //      GetPatchFecundityFuncPtr = &LCE_Breed_base::stochasticLogisticGrowth;
  //      break;
  //    case 4:
  //      GetPatchFecundityFuncPtr = &LCE_Breed_base::conditionalLogisticGrowth;
  //      break;
  //    case 5:
  //      GetPatchFecundityFuncPtr = &LCE_Breed_base::conditionalStochasticLogisticGrowth;
  //      break;
  //    case 6:
  //      GetPatchFecundityFuncPtr = &LCE_Breed_base::fixedFecundityGrowth;
  //      break;
  //    case 7:
  //      GetPatchFecundityFuncPtr = &LCE_Breed_base::stochasticFecundityGrowth;
  //      break;
  //    default:
  //      GetPatchFecundityFuncPtr = &LCE_Breed_base::instantGrowth;
  //      break;
  //  }
  //  
  //  //growth rate
  //  if (model > 1 && model < 6) {
  //    
  //    if(!_paramSet->isSet("growth_rate")) {
  //      error("parameter \"growth_rate\" needs to be set\n");
  //      return false;
  //    }
  //    
  //    if(_growthRates) delete [] _growthRates;
  //    _growthRates = new double [ _popPtr->getPatchNbr() ];
  //    
  //    if(_paramSet->isMatrix("growth_rate")) {
  //      
  //      TMatrix tmp;
  //      
  //      _paramSet->getMatrix("growth_rate", &tmp);
  //      
  //      if(tmp.getNbCols() != _popPtr->getPatchNbr()){
  //        error("matrix argument to \"growth_rate\" has wrong number of elements,\
  //              must equal the number of patches.\n");
  //        return false;
  //      }
  //      
  //      for (unsigned int i = 0; i < _popPtr->getPatchNbr(); i++) {
  //        _growthRates[i] = tmp.get(0, i);
  //      }
  //      
  //      
  //    } else { //not a matrix
  //      
  //      _growthRates[0] = get_parameter_value("growth_rate");
  //      
  //      for (unsigned int i = 1; i < _popPtr->getPatchNbr(); i++) {
  //        _growthRates[i] = _growthRates[0];
  //      }
  //    }
  //    
  //    
  //  }
  
  return true;
}
// ----------------------------------------------------------------------------------------
// LCE_Breed_base::setSexRatio
// ----------------------------------------------------------------------------------------
bool LCE_Breed_base::setSexRatio ()
{
  // SEX-RATIO
  
  if(get_parameter("sex_ratio_mode")->isSet()) {
    
    if(get_parameter("sex_ratio_mode")->getArg().compare("fixed") == 0) 
      
      GetOffsprgSex = &LCE_Breed_base::getOffsprgSexFixed;
    
    else if(get_parameter("sex_ratio_mode")->getArg().compare("random") != 0) {
      
      error("\"sex_ratio_mode\" parameter argument must be either \"fixed\" or \"random\".");
      return false;
    } else
      GetOffsprgSex = &LCE_Breed_base::getOffsprgSexRandom;
  } else {
    
    switch(_mating_system) {
      case 4: GetOffsprgSex = &LCE_Breed_base::getOffsprgSexSelfing;
        break;
      case 5: GetOffsprgSex = &LCE_Breed_base::getOffsprgSexCloning;
        break;
      case 6: GetOffsprgSex = &LCE_Breed_base::getOffsprgSexSelfing;
        break;
      default: GetOffsprgSex = &LCE_Breed_base::getOffsprgSexRandom;
    }
  }
  return true;
}
// ----------------------------------------------------------------------------------------
// LCE_Breed_base::setFecundity
// ----------------------------------------------------------------------------------------
bool LCE_Breed_base::setFecundity ()
{
  // FECUNDITY
  _mean_fecundity = get_parameter_value("mean_fecundity");
  
  if(get_parameter("fecundity_distribution")->isSet()) {
    
    string dist = get_parameter("fecundity_distribution")->getArg();
    
    if( dist.compare("fixed") == 0 )
      
      FecundityFuncPtr = &LCE_Breed_base::getFixedFecundity;
    
    else if( dist.compare("poisson") == 0 )
      
      FecundityFuncPtr = &LCE_Breed_base::getPoissonFecundity;
    
    else if( dist.compare("normal") == 0 ) {
      
      FecundityFuncPtr = &LCE_Breed_base::getGaussianFecundity;
      
      if(get_parameter("fecundity_dist_stdev")->isSet())
        _sd_fecundity = get_parameter_value("fecundity_dist_stdev");
      else {
        error("parameter \"fecundity_dist_stdev\" is missing!\n");
        return false;
      }
      
    } else {
      error("unknown fecundity distribution parameter's argument!\n");
      return false;
    }
    
  } else { //default distribution is Poisson:
    FecundityFuncPtr = &LCE_Breed_base::getPoissonFecundity;
  }
  return true;
}
// ----------------------------------------------------------------------------------------
// LCE_Breed_base::getOffsprgSexFixed
// ----------------------------------------------------------------------------------------
sex_t LCE_Breed_base::getOffsprgSexFixed  ()
{
  static bool sex = RAND::RandBool();
  sex ^= 1;
  return (sex_t)sex;
}

// ----------------------------------------------------------------------------------------
// LCE_Breed_base::makeOffspring
// ----------------------------------------------------------------------------------------
Individual* LCE_Breed_base::makeOffspring(Individual* ind)
{
//  unsigned int cat = ind->getPedigreeClass();  //MOD SLIM NEMO 7.11
//
//  ind->getMother()->DidHaveABaby(cat);
//  if(cat!=4) ind->getFather()->DidHaveABaby(cat);
  
  return ind->create(doInheritance(), true);
}
// ----------------------------------------------------------------------------------------
// LCE_Breed_base::breed
// ----------------------------------------------------------------------------------------
Individual* LCE_Breed_base::breed(Individual* mother, Individual* father, unsigned int LocalPatch)
{
  return _popPtr->makeNewIndividual(mother, father, getOffsprgSex(), LocalPatch);
}
// ----------------------------------------------------------------------------------------
// LCE_Breed_base::breed_cloning
// ----------------------------------------------------------------------------------------
Individual* LCE_Breed_base::breed_cloning(Individual* mother, Individual* father, unsigned int LocalPatch)
{
  Individual *newind;
  
  if(mother == father) {
    newind = _popPtr->getNewIndividual();
    //cloning:
    (*newind) = (*mother); 
    
//    newind->reset_counters();  //MOD SLIM NEMO 7.11
    newind->setFather(NULL);
//    newind->setFatherID(0);//MOD SLIM NEMO 7.11
    newind->setMother(mother);
//    newind->setMotherID(mother->getID());//MOD SLIM NEMO 7.11
//    newind->setIsSelfed(true);  //MOD SLIM NEMO 7.11
    newind->setHome(LocalPatch);
    newind->setAge(0);
    _do_inherit = false;  
  } else {
    newind =  _popPtr->makeNewIndividual(mother, father, getOffsprgSex(), LocalPatch);
    _do_inherit = true;
  }
  
  return newind;
}
// ----------------------------------------------------------------------------------------

//                                     LCE_Breed

LCE_Breed::LCE_Breed ( ) : LifeCycleEvent("breed",""),//constructor inherit from lifecycleevent
_fecundities(0)
{
	add_parameter("breed_competition_coefficient",DBL,false,false,0,10,0);//add parameter
  //special syntax with ifset() see above if the parameter is not mandatory
  
}
// ----------------------------------------------------------------------------------------
// LCE_Breed::init
// ----------------------------------------------------------------------------------------
bool LCE_Breed::setParameters ( )
{ 
	//@TODO it is not a mandatory paramter, and it doesn't have a default value!
	_coeff_competition = get_parameter_value("breed_competition_coefficient");

	if(! _LeslieMatrix) _LeslieMatrix = new TMatrix(*_popPtr->getLeslieMatrix());

	unsigned int nb_class = _popPtr->getLeslieMatrix()->getNbCols();
  
	if(_fecundities) delete [] _fecundities;
  
	_fecundities = new double [nb_class];
  
	_LeslieMatrix->getRowView(0, nb_class, _fecundities);
  
	for(unsigned int i = 0; i < nb_class; ++i) _fecundities[i] = round(_fecundities[i]);
  
	//_carrying_capacity = get_parameter_value("Carrying_capacity_perpatch");
  
  return LCE_Breed_base::setParameters();
}
// ----------------------------------------------------------------------------------------
// LCE_Breed::execute
// ----------------------------------------------------------------------------------------
void LCE_Breed::execute()
{
  
	//declarations
  
  Patch* patch;
  Individual* mother;
  Individual* father;
  Individual* NewOffsprg;
  unsigned int nbBaby;
//  float  nb_alive_baby;
//  unsigned int entire_nb_alive_baby;
  unsigned int pop_size;
  double compet;
  unsigned int sum;
  age_idx current_age;
  
  unsigned int nb_class = _popPtr->getLeslieMatrix()->getNbCols();
  
#ifdef _DEBUG_
  message("LCE_Breed::execute (Patch nb: %i offsprg nb: %i adlt nb: %i)\n"
          ,_popPtr->getPatchNbr(),_popPtr->size( OFFSPRG ),_popPtr->size( ADULTS ));
#endif
  
  if(_popPtr->size(OFFSPRG) != 0 and _LeslieMatrix->get(0, 0)==0) {
    
    warning("offspring containers not empty at time of breeding, flushing.\n");

    _popPtr->flush(OFFSx);

  }
  
  for(unsigned int i = 0; i < _popPtr->getPatchNbr(); i++) {
    
	  patch = _popPtr->getPatch(i);
    
	  pop_size = patch->size(ADULTS) + patch->size(OFFSPRG);
    
    compet = 1/(1 + pop_size * _coeff_competition);
    
    //check
    //cout<<"patch-size"<<patch->size(FEM, ADLTx)<<endl;
    //
    
	  if( !checkMatingCondition(patch) ) continue;
    
	  for(int j = (int)nb_class-1; j >= 0; j--){
      
      current_age = static_cast<age_idx> (j);

      // MOD F.G. 01.05.14
      // fecunditites are now stored during initialisation in _fecundities.
      
//      fec = _LeslieMatrix->get(0, current_age);
//      
//      fec = round(fec);
      
      for(unsigned int size = patch->size(FEM, current_age), indexOfMother = 0;
          indexOfMother < size;
          indexOfMother++)
      {
        mother = patch->get(FEM, current_age, indexOfMother);
                
        nbBaby = (unsigned int)mother->setFecundity( _fecundities[current_age] ) ;//nb babies produced
        
        
        ///check
        // cout<<nbBaby<<endl;
        //  test=test+1;
        //	cout<<test<<endl;
        //////
        
        ///need to add a density dependent term for nb babies alive
        // it is a multiplicative term witch depend on the number of individual in the patch in the idelal case
        //including offspring
        
        // if( get_parameter("competition_coefficient")->isSet()){
        
        //compet = 1/(1 + pop_size * _coeff_competition);
        
        //}
        
        //if(get_parameter("Carrying_capacity_perpatch")->isSet()){
        
        //compet = (_carrying_capacity - pop_size)/_carrying_capacity;}
        
        //else{compet = 1;}
        
        sum=0;
        
        for( unsigned int k=0; k < nbBaby; k++ )
          sum += RAND::Bernoulli2(compet);
        
        //      nb_alive_baby = sum;
        
        nbBaby = sum;
        
        ////check
        //cout<<"nb Baby = "<<nbBaby<<endl;
        ////
        
        
        //-----------------------------------------------------------------------
        while(nbBaby != 0) {
          
          father = this->getFatherPtr(patch, mother, indexOfMother);
          
          NewOffsprg = makeOffspring( do_breed(mother, father, i) );
          
          patch->add(NewOffsprg->getSex(), OFFSx, NewOffsprg);
          
          nbBaby--;
        }//_END__WHILE
        
        
      }
      
      
      
    }


    
#ifdef _DEBUG_
    cout<<"make new offspring"<<endl;
#endif
  }
  //check//////////////////////////////////////////////////////////////
  
  //cout<<"pop_size: "<<pop_size<<endl;
  
  //message("LCE_Breed::execute (Patch nb: %i offsprg nb: %i adlt nb: %i fem off: %i)\n"
  //        ,_popPtr->getPatchNbr(),_popPtr->size( OFFSPRG ),_popPtr->size( ADULTS), _popPtr->size(FEM, OFFSx) );
  
  //cout<<"get fec: "<<getFecundity()<<endl;
  
  
  
  
  
  /////////////////////////////////////////////
  
}

// _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/ 

//                        ******** LCE_Breed_Wolbachia ********/

//----------------------------------------------------------------------------------------
// LCE_Breed_Wolbachia
// ----------------------------------------------------------------------------------------
LCE_Breed_Wolbachia::LCE_Breed_Wolbachia () : LifeCycleEvent("breed_wolbachia",WOLB),
_incomp_cost(0), 
_fec_cost(0),
_infected_fec(0),
_inoculum_size(0),
_inoculum_time(0),
_writer(0)
{
  ParamUpdater< LCE_Breed_Wolbachia > * updater =
  new ParamUpdater< LCE_Breed_Wolbachia > (&LCE_Breed_Wolbachia::setParameters);
  add_parameter("wolbachia_fecundity_cost",DBL,true,true,-1,1, updater);
  add_parameter("wolbachia_incompatibility_cost",DBL,true,true,-1,1, updater);
  add_parameter("wolbachia_inoculum_size",INT,true,false,0,0, updater);
  add_parameter("wolbachia_inoculum_time",INT,true,false,0,0, updater);
  add_parameter("wolbachia_model",INT,false,true,1,2, updater);
  add_parameter("wolbachia_output_dir", STR, false, false, 0, 0);
}

LCE_Breed_Wolbachia::~LCE_Breed_Wolbachia ( ) 
{
  if(_writer) delete _writer; 
}
// ----------------------------------------------------------------------------------------
// LCE_Breed_Wolbachia::setParameters
// ----------------------------------------------------------------------------------------
bool LCE_Breed_Wolbachia::setParameters ()
{
  if(!LCE_Breed_base::setParameters( )) return false;
  
  _fec_cost = this->get_parameter_value("wolbachia_fecundity_cost");
  _infected_fec = this->getMeanFecundity() * (1 - _fec_cost);
  _incomp_cost = this->get_parameter_value("wolbachia_incompatibility_cost");
  _inoculum_size = (unsigned int)this->get_parameter_value("wolbachia_inoculum_size");
  _inoculum_time = (unsigned int)this->get_parameter_value("wolbachia_inoculum_time");
  
  if (get_parameter("wolbachia_model")->isSet()) {
    _model = (unsigned int)get_parameter_value("wolbachia_model");
  } else
    _model = 1;
  
  
  switch(_model) {
    case 1: 
      _breed_func_ptr = &LCE_Breed_Wolbachia::wolbachia_model_1;
      break;
    case 2:
      _breed_func_ptr = &LCE_Breed_Wolbachia::wolbachia_model_2;
      break;
    default:  _breed_func_ptr = &LCE_Breed_Wolbachia::wolbachia_model_1;
  }
  
  
  return true;
}
// ----------------------------------------------------------------------------------------
// LCE_Breed_Wolbachia::loadFileServices
// ----------------------------------------------------------------------------------------
void LCE_Breed_Wolbachia::loadFileServices ( FileServices* loader )
{
  if(_writer != NULL) delete _writer;
  
  _writer = new TTWolbachiaFH(this);
  
  _writer->set(false, false, SIMenv::getReplicates(), SIMenv::getGenerations(), 0,
               get_parameter("wolbachia_output_dir")->getArg(), this);
  
  loader->attach(_writer);
}
//----------------------------------------------------------------------------------------
// LCE_Breed_Wolbachia::execute
// ----------------------------------------------------------------------------------------
void LCE_Breed_Wolbachia::execute ()
{
  
#ifdef _DEBUG_
  message("LCE_Breed_Wolbachia::execute\n");
#endif
  double infection_status;
  
  if( _popPtr->getCurrentGeneration() == _inoculum_time) inoculate_wolbachia();
  
  if( _popPtr->getCurrentGeneration() > _inoculum_time) {
    
    infection_status = hasInfectedFemale();
    
    if( infection_status == 0 || infection_status == 1) {
      _writer->record(_popPtr->getCurrentReplicate(), _popPtr->getCurrentGeneration(), infection_status);
      _popPtr->reset();
      return;
    }}
  
  
  if(_popPtr->size(OFFSPRG) != 0) {
    warning("offspring containers not empty at time of breeding, flushing.\n");
    _popPtr->flush(OFFSx);
  }
  
  (this->*_breed_func_ptr)();
  
  if (_popPtr->getCurrentGeneration() == SIMenv::getGenerations() ) {
    _writer->record(_popPtr->getCurrentReplicate(), _popPtr->getCurrentGeneration(), hasInfectedFemale());
  }
  
}
//----------------------------------------------------------------------------------------
// LCE_Breed_Wolbachia::wolbachia_model_1
// ----------------------------------------------------------------------------------------
void LCE_Breed_Wolbachia::wolbachia_model_1 ()
{
  Patch* current_patch;
  unsigned int indexOfMother, nbBaby;
  Individual* FatherPtr;
  Individual* MotherPtr;
  Individual* NewOffsprg;
  
  for(unsigned int i = 0; i < _popPtr->getPatchNbr(); i++) {
    
    current_patch = _popPtr->getPatch(i);
    
    if(current_patch->size(FEM, ADLTx) == 0 || current_patch->size(MAL, ADLTx) == 0) continue;
    
    for(indexOfMother = 0; indexOfMother < current_patch->size(FEM, ADLTx); indexOfMother++) {
      
      MotherPtr = current_patch->get(FEM, ADLTx, indexOfMother);
      
      if(*(bool*)MotherPtr->getTraitValue(_LCELinkedTraitIndex))
        nbBaby = (unsigned int)MotherPtr->setFecundity( getFecundity( _infected_fec ) );
      else
        nbBaby = (unsigned int)MotherPtr->setFecundity( getFecundity() );
      //-----------------------------------------------------------------------
      while(nbBaby != 0) {
        
        FatherPtr = getFatherPtr(current_patch, MotherPtr, indexOfMother);
        
        NewOffsprg = _popPtr->makeOffsprg(MotherPtr,FatherPtr,(sex_t)RAND::RandBool(),i);
        
        if(!(*(bool*)NewOffsprg->getTraitValue(_LCELinkedTraitIndex)) &&
           *(bool*)FatherPtr->getTraitValue(_LCELinkedTraitIndex)) {
          if(RAND::Uniform() < _incomp_cost) {
            _popPtr->recycle(NewOffsprg);
          } else
            current_patch->add( NewOffsprg->getSex(), OFFSx, NewOffsprg );
        } else
          current_patch->add( NewOffsprg->getSex(), OFFSx, NewOffsprg );
        
        nbBaby--;
      }//_END_WHILE nbBaby
    }//end_for indexOfMother
  }//end_for patch
}
//----------------------------------------------------------------------------------------
// LCE_Breed_Wolbachia::wolbachia_model_2
// ----------------------------------------------------------------------------------------
void LCE_Breed_Wolbachia::wolbachia_model_2 ()
{
  Patch* current_patch;
  unsigned int indexOfMother, nbBaby;
  double fec;
  Individual* FatherPtr;
  Individual* MotherPtr;
  Individual* NewOffsprg;
  
  for(unsigned int i = 0; i < _popPtr->getPatchNbr(); i++) {
    
    current_patch = _popPtr->getPatch(i);
    
    if(current_patch->size(FEM, ADLTx) == 0 || current_patch->size(MAL, ADLTx) == 0) continue;
    
    for(indexOfMother = 0; indexOfMother < current_patch->size(FEM, ADLTx); indexOfMother++) {
      
      MotherPtr = current_patch->get(FEM, ADLTx, indexOfMother);
      
      FatherPtr = getFatherPtr(current_patch, MotherPtr, indexOfMother);
      
      if(*(bool*)MotherPtr->getTraitValue(_LCELinkedTraitIndex)) {
        if(*(bool*)FatherPtr->getTraitValue(_LCELinkedTraitIndex))
          fec = _infected_fec * (1 - _incomp_cost);
        else
          fec = _infected_fec;
      } else if (*(bool*)FatherPtr->getTraitValue(_LCELinkedTraitIndex)) {
        fec = getMeanFecundity() * (1 - _incomp_cost);
      } else {
        fec = getMeanFecundity();
      }
      
      nbBaby = (unsigned int)MotherPtr->setFecundity( getFecundity( fec ) );
      
      //-----------------------------------------------------------------------
      for(;nbBaby != 0;nbBaby--) {
        
        NewOffsprg = _popPtr->makeOffsprg(MotherPtr,FatherPtr,(sex_t)RAND::RandBool(),i);
        
        current_patch->add( NewOffsprg->getSex(), OFFSx, NewOffsprg );
        
      }//end_for nbBaby
    }//end_for indexOfMother
  }//end_for patch
}
//----------------------------------------------------------------------------------------
// LCE_Breed_Wolbachia::inoculate_wolbachia
// ----------------------------------------------------------------------------------------
void LCE_Breed_Wolbachia::inoculate_wolbachia ()
{
  Patch* current_patch;
  bool T = 1, ok = false;
  
  for (unsigned int i = 0; i < _popPtr->getPatchNbr(); i++) {
    
    current_patch = _popPtr->getPatch(i);
    
    if(current_patch->get_isExtinct()) continue;
    
    if(_inoculum_size <= current_patch->size(FEM, ADLTx) && _inoculum_size <= current_patch->size(MAL, ADLTx)){
      for (unsigned int j = 0; j < _inoculum_size; j++) {
        current_patch->get(FEM, ADLTx, j)->setTrait(_LCELinkedTraitIndex,&T);
        //        current_patch->get(MAL, ADLTx, j)->setTrait(_LCELinkedTraitIndex,&T);
      }
      ok = true;
    }
    break;
  }
  if(!ok) fatal("could not inoculate wolbachia, check the inoculum size!\n");
}
//----------------------------------------------------------------------------------------
// LCE_Breed_Wolbachia::getFemaleInfection
// ----------------------------------------------------------------------------------------
double LCE_Breed_Wolbachia::hasInfectedFemale ()
{
  double indNbr = 0, mean = 0, size;
  Patch* crnt_patch;
  
  for(unsigned int i = 0; i < _popPtr->getPatchNbr(); ++i) {
    crnt_patch = _popPtr->getPatch(i);
    size = crnt_patch->size(FEM, ADLTx);
    indNbr += size;
    for(unsigned int j = 0; j < size; ++j)
      mean += (double)*(bool*)crnt_patch->get(FEM, ADLTx, j)->getTraitValue(_LCELinkedTraitIndex);
  }
  
  return (indNbr != 0 ? mean/indNbr : nanf("NULL"));
}

/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

//                             ******** TTWolbachiaFH ********

// ----------------------------------------------------------------------------------------
// record
// ----------------------------------------------------------------------------------------
void TTWolbachiaFH::record (unsigned int repl, unsigned int gen, double infection)
{
  _times[repl] = gen;
  _rate.push_back(infection);
}
// ----------------------------------------------------------------------------------------
// FHwrite
// ----------------------------------------------------------------------------------------
void TTWolbachiaFH::FHwrite ()
{
  string filename = get_filename();
  
  ofstream FILE (filename.c_str(), ios::out); 
  
  if(!FILE) fatal("could not open wolbachia output file \"%s\"\n", filename.c_str());
  
  FILE << "replicate\t" << "generation\t" << "infection_rate\n";
  
  unsigned int i = 0;
  
  map<unsigned int, unsigned int>::iterator IT = _times.begin();
  
  while (IT != _times.end() && i < _rate.size() ) {
    FILE << IT->first <<"\t"<< IT->second << "\t" << _rate[i] << endl;
    i++;
    IT++;
  }
  
  FILE.close();
}

//####################################################################################

//										CLONAGE										//

//####################################################################################

LCE_clonage::LCE_clonage():LifeCycleEvent ("cloning",""), _clonal_rate(0),
		_coeff_competition_cloning(0), _age_target(ADLTx) {
  
	add_parameter("cloning_rate",DBL,true,false,0,0);
  
	add_parameter("cloning_competition_coefficient",DBL,false,false,0,0);
  
    add_parameter("cloning_add_to_age_class",INT,false,false,0,0);
}

bool LCE_clonage::setParameters (){
  
  if(get_parameter("cloning_rate")->isSet())
	  _clonal_rate = get_parameter_value("cloning_rate");

  if(get_parameter("cloning_competition_coefficient")->isSet())
	  _coeff_competition_cloning = get_parameter_value("cloning_competition_coefficient");//clonage specifique
  
  if(get_parameter("cloning_add_to_age_class")->isSet())
	  _age_target = static_cast<age_idx>((unsigned int)get_parameter_value("cloning_add_to_age_class"));
  else
	  _age_target = ADLTx;
  
  return true;
}

void LCE_clonage::execute ()
{
  
	//cout<<"clonage::execute"<<endl;
  
  //declaration
  Individual* newind;
  Patch* patch;
  Individual* mother;
  unsigned int pop_size;
  double number_clones;
  double compet_clone;
  
  
  //loop for each patch, for each adult female
  
  for(unsigned int i = 0; i < _popPtr->getPatchNbr(); i++) {
    
    patch = _popPtr->getPatch(i);
    
    pop_size = patch->size(ADULTS); // + patch->size(OFFSPRG);
    
    for(unsigned int size = patch->size(FEM, ADULTS),
        indexOfMother = 0;indexOfMother < size; indexOfMother++)
    {
      
      //number of individuals to be generated from a poisson distribution
      mother = patch->get(FEM, ADULTS, indexOfMother);
      
      compet_clone = 1/(1 + pop_size * _coeff_competition_cloning);
      
      //cout<<"_clonal_rate "<<_clonal_rate<<endl;
      
      //cout<<"compet "<<compet<<endl;
      
      
      number_clones = compet_clone * _clonal_rate;
      //cout<<"number_clones1 "<<number_clones<<endl;
      
      number_clones = round(RAND::Poisson(number_clones));
      
      //cout<<"number_clones2 "<<number_clones<<endl;
      
      
      for (int k = 0 ; k < number_clones; k++){
        
        newind = _popPtr->getNewIndividual();
        
        //cloning:
        
        (*newind) = (*mother);
//        newind->reset_counters();  //MOD SLIM NEMO 7.11
        newind->setFather(NULL);
//        newind->setFatherID(0);//MOD SLIM NEMO 7.11
        newind->setMother(mother);
//        newind->setMotherID(mother->getID());//MOD SLIM NEMO 7.11
//        newind->setIsSelfed(true);  //MOD SLIM NEMO 7.11
        newind->setHome(i);
        newind->setAge(_age_target);
        
        patch->add(FEM, _age_target, newind->create(false,true));
        
        //cout<<k<<endl;
      }
    }
  }
  
  
  // we add new individuals in the offspring age class
  
  
  
  //cout<<"end_execute"<<endl;
  
}

