/**  $Id: LCEcomposite.cc,v 1.10.2.2 2015-03-16 20:17:01 fred Exp $
 *
 *  @file LCEcomposite.cc
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
 *  created on @date 09.10.2007
 *
 *  @author fred
 */

#include "LCEcomposite.h"

/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

//                             ******** LCE_Breed_Disperse ********/

// ----------------------------------------------------------------------------------------
// LCE_Breed_Disperse::LCE_Breed_Disperse
// ----------------------------------------------------------------------------------------
LCE_Breed_Disperse::LCE_Breed_Disperse() : LifeCycleEvent("breed_disperse",""),_breed_disperse(0),
_num_colonizers(0),_growthRates(0),_make_offspring(0),_get_numFemOffspring(0),
_get_numMalOffspring(0),_get_patchFecundity(0)
{
  ParamUpdater<LCE_Breed_Disperse>* updater = new ParamUpdater<LCE_Breed_Disperse>(&LCE_Breed_Disperse::setParameters);
  
  LCE_Disperse_base::addParameters("breed_disperse", updater);
  
  add_parameter("breed_disperse_dispersing_sex", STR, false, false, 0, 0, updater);
  add_parameter("breed_disperse_colonizers", INT,false,false,0,0, updater);
  add_parameter("breed_disperse_growth_model", INT, false, true, 1, 7, updater);
  add_parameter("breed_disperse_growth_rate", DBL, false, false, 0, 0, updater);
}
// ----------------------------------------------------------------------------------------
// LCE_Breed_Disperse::setParameters
// ----------------------------------------------------------------------------------------
bool LCE_Breed_Disperse::setParameters ()
{
  if(!LCE_Breed_base::setParameters()) return false;
  
  LCE_Disperse_base::set_isForward(false);
  
  if(!LCE_Disperse_base::setBaseParameters("breed_disperse")) return false;
  
  _num_colonizers = (int)get_parameter_value("breed_disperse_colonizers");
  
  _get_numFemOffspring = (_num_colonizers == -1 ?
                          &LCE_Breed_Disperse::numFemOffspring : 
                          &LCE_Breed_Disperse::numFemOffspring_colonizers);
  _get_numMalOffspring = &LCE_Breed_Disperse::numMalOffspring_notrandom;
  
  if(get_parameter("breed_disperse_dispersing_sex")->isSet()) {
    
    if (get_parameter("breed_disperse_dispersing_sex")->getArg() == "female") {
      _dispersing_sex = FEM;
    } else if (get_parameter("breed_disperse_dispersing_sex")->getArg() == "male") {
      _dispersing_sex = MAL;
    } else {
      error("the \"breed_disperse_dispersing_sex\" parameter only takes \"female\" or \"male\" as argument\n");
    }
    
  } else {
    _dispersing_sex = FEM;
  }
  
  //check mating system:  
  unsigned int model = this->getMatingSystem();
  
  if(model == 2 || model == 3) {
    error("Polygyny and Monogamy are not implemented within the breed_disperse LCE.\n");
    return false;
  } else if(model == 1) {
    _make_offspring = &LCE_Breed_Disperse::mate_random;
    _get_numMalOffspring = (_num_colonizers == -1 ?
                            &LCE_Breed_Disperse::numMalOffspring_random : 
                            &LCE_Breed_Disperse::numMalOffspring_random_colonizers);
  } else if(model == 4)
    if(this->getMatingProportion() != 1)
      _make_offspring = &LCE_Breed_Disperse::mate_selfing;
    else
      _make_offspring = &LCE_Breed_Disperse::mate_full_selfing;
    else if(model == 5)
      _make_offspring = &LCE_Breed_Disperse::mate_cloning;
    else if(model == 6)
      _make_offspring = &LCE_Breed_Disperse::mate_random_hermaphrodite;
  
  //check dispersal model:
  model = this->getDispersalModel();
  if(model == 2)
    _breed_disperse = &LCE_Breed_Disperse::do_breed_disperse_propagule;
  else
    _breed_disperse = &LCE_Breed_Disperse::do_breed_disperse;
  
  //growth model
  if(_paramSet->isSet("breed_disperse_growth_model"))
    model = get_parameter_value("breed_disperse_growth_model");
  else 
    model = 1;
  
  switch (model) {
    case 1:
      _get_patchFecundity = &LCE_Breed_Disperse::instantGrowth;
      break;
    case 2:
      _get_patchFecundity = &LCE_Breed_Disperse::logisticGrowth;
      break;
    case 3:
      _get_patchFecundity = &LCE_Breed_Disperse::stochasticLogisticGrowth;
      break;
    case 4:
      _get_patchFecundity = &LCE_Breed_Disperse::conditionalLogisticGrowth;
      break;
    case 5:
      _get_patchFecundity = &LCE_Breed_Disperse::conditionalStochasticLogisticGrowth;
      break;
    case 6:
      _get_patchFecundity = &LCE_Breed_Disperse::fixedFecundityGrowth;
      break;
    case 7:
      _get_patchFecundity = &LCE_Breed_Disperse::stochasticFecundityGrowth;
      break;
    default:
      _get_patchFecundity = &LCE_Breed_Disperse::instantGrowth;
      break;
  }
  
  //growth rate
  if (model > 1 && model < 6) {
    
    if(!_paramSet->isSet("breed_disperse_growth_rate")) {
      error("parameter \"breed_disperse_growth_rate\" needs to be set\n");
      return false;
    }
    
    if(_growthRates) delete [] _growthRates;
    _growthRates = new double [ _popPtr->getPatchNbr() ];
    
    if(_paramSet->isMatrix("breed_disperse_growth_rate")) {
      
      TMatrix tmp;
      
      _paramSet->getMatrix("breed_disperse_growth_rate", &tmp);
      
      if(tmp.getNbCols() != _popPtr->getPatchNbr()){
        error("matrix argument to \"breed_disperse_growth_rate\" has wrong number of elements,\
              must equal the number of patches.\n");
        return false;
      }
      
      for (unsigned int i = 0; i < _popPtr->getPatchNbr(); i++) {
        _growthRates[i] = tmp.get(0, i);
      }
      
      
    } else { //not a matrix
      
      _growthRates[0] = get_parameter_value("breed_disperse_growth_rate");
      
      for (unsigned int i = 1; i < _popPtr->getPatchNbr(); i++) {
        _growthRates[i] = _growthRates[0];
      }
    }
    
    
  }
  
  return true;
}
// ----------------------------------------------------------------------------------------
// LCE_Breed_Disperse::execute
// ----------------------------------------------------------------------------------------
void LCE_Breed_Disperse::execute()
{
#ifdef _DEBUG_
  message("LCE_Breed_Disperse::execute (Patch nb: %i offsprg nb: %i adlt nb: %i "
          ,_popPtr->getPatchNbr(), _popPtr->size( OFFSPRG ), _popPtr->size( ADULTS ));
#endif
  
  if(_npatch != _popPtr->getPatchNbr()) {
    _npatch = _popPtr->getPatchNbr();
    if(!updateDispMatrix()) fatal("bailing out\n");
  }
  
  LCE_Disperse_base::reset_counters();
  
  (this->*_breed_disperse)();
  
#ifdef _DEBUG_
  unsigned int a = 0, b = 0, c = 0;
  for(unsigned int i = 0; i < _popPtr->getPatchNbr(); i++){
    a += _popPtr->getPatch(i)->nbEmigrant;
    b += _popPtr->getPatch(i)->nbImigrant;
    c += _popPtr->getPatch(i)->nbPhilopat;
  }
  message("immigrate: %f, emigrants: %i, imigrants: %i)\n",(double)b/(b+c), a, b);
#endif
}
// ----------------------------------------------------------------------------------------
// LCE_Breed_Disperse::do_breed_disperse
// ----------------------------------------------------------------------------------------
void LCE_Breed_Disperse::do_breed_disperse()
{
  Patch* patch;
  unsigned int nfem, nmal;
  
  for(unsigned int i = 0; i < _popPtr->getPatchNbr(); i++) {
    
    patch = _popPtr->getPatch(i);
    
    if(patch->size(OFFSx) != 0) patch->flush(OFFSx, _popPtr);
    
    nfem = (this->* _get_numFemOffspring)(patch);
    nmal = (this->* _get_numMalOffspring)(patch);
    
    for(unsigned int j = 0; j < nfem; j++) 
      patch->add(FEM, OFFSx, 
                 LCE_Breed_base::makeOffspring( (this->*_make_offspring)(FEM, patch, i) )
                 );
    
    for(unsigned int j = 0; j < nmal; j++) 
      patch->add(MAL, OFFSx, 
                 LCE_Breed_base::makeOffspring( (this->*_make_offspring)(MAL, patch, i) )
                 );
  }
}
// ----------------------------------------------------------------------------------------
// LCE_Breed_Disperse::do_breed_disperse_propagule
// ----------------------------------------------------------------------------------------
void LCE_Breed_Disperse::do_breed_disperse_propagule()
{
  this->setIsland_PropagulePool_Matrix();
  do_breed_disperse();
}
// ----------------------------------------------------------------------------------------
// LCE_Breed_Disperse::mate_random
// ----------------------------------------------------------------------------------------
Individual* LCE_Breed_Disperse::mate_random(sex_t SEX, Patch *patch, unsigned int LocalPatch)
{
  return _popPtr->makeNewIndividual(get_parent(FEM, FEM, patch, LocalPatch), 
                                    get_parent(MAL, MAL, patch, LocalPatch),
                                    SEX, LocalPatch);
}
// ----------------------------------------------------------------------------------------
// LCE_Breed_Disperse::mate_random_hermaphrodite
// ----------------------------------------------------------------------------------------
Individual* LCE_Breed_Disperse::mate_random_hermaphrodite(sex_t SEX, Patch *patch,
                                                          unsigned int LocalPatch)
{
  return _popPtr->makeNewIndividual(get_parent(FEM, FEM, patch, LocalPatch), 
                                    get_parent(FEM, _dispersing_sex, patch, LocalPatch),
                                    FEM, LocalPatch);
}
// ----------------------------------------------------------------------------------------
// LCE_Breed_Disperse::mate_selfing
// ----------------------------------------------------------------------------------------
Individual* LCE_Breed_Disperse::mate_selfing(sex_t SEX, Patch *patch, unsigned int LocalPatch)
{
  Individual *mom, *dad;
  //  cout << "--- choosing mum in patch "<<LocalPatch<<" ("<<patch->size(FEM, ADLTx)<<")\n";
  mom = get_parent(FEM, FEM, patch, LocalPatch);
  
  //  cout << "--- choosing dad \n";
  
  if( RAND::Uniform() > this->getMatingProportion() ) {
    //random mating
    do{
      dad = get_parent(FEM, _dispersing_sex, patch, LocalPatch);
    }while(dad == mom && _popPtr->getPatchNbr() != 1 && patch->size(FEM,ADLTx) != 1);
  } else
    dad = mom;
  
  //  cout << "--- make offspring by crossing "<<mom->getID()<<" with "<<dad->getID()<<endl;
  
  
  return _popPtr->makeNewIndividual(mom, dad, FEM, LocalPatch);
}
// ----------------------------------------------------------------------------------------
// LCE_Breed_Disperse::mate_full_selfing
// ----------------------------------------------------------------------------------------
Individual* LCE_Breed_Disperse::mate_full_selfing(sex_t SEX, Patch *patch, unsigned int LocalPatch)
{
  Individual *mom = get_parent(FEM, FEM, patch, LocalPatch);
  
  return _popPtr->makeNewIndividual(mom, mom, FEM, LocalPatch);
}
// ----------------------------------------------------------------------------------------
// LCE_Breed_Disperse::mate_cloning
// ----------------------------------------------------------------------------------------
Individual* LCE_Breed_Disperse::mate_cloning(sex_t SEX, Patch *patch, unsigned int LocalPatch)
{    
  return LCE_Breed_base::breed_cloning(get_parent(FEM, FEM, patch, LocalPatch), NULL, LocalPatch);
}
// ----------------------------------------------------------------------------------------
// LCE_Breed_Disperse::get_parent
// ----------------------------------------------------------------------------------------
Individual* LCE_Breed_Disperse::get_parent(sex_t SEX, sex_t DispSex, Patch *local_patch,
                                           unsigned int LocalPatch)
{
  unsigned int SourcePatch;
  Patch* src_patch;
  //  cout << "  getting "<<(DispSex==1?"female":"male")<<" gamete to "<<LocalPatch<<flush;
  do {
    
    SourcePatch = LCE_Disperse_base::getMigrationPatchBackward(DispSex, LocalPatch);
    src_patch = _popPtr->getPatchPtr(SourcePatch);
    
  } while (src_patch->size( SEX, ADLTx ) == 0); //redraw if source patch is empty
  //  cout << " from "<<SourcePatch<<endl;
  
  if(LocalPatch != SourcePatch) {
    src_patch->nbEmigrant++;
    
    if( local_patch->size(ADLTx) == 0 )
      local_patch->nbKolonisers++;
    else
      local_patch->nbImigrant++;
    
  } else
    src_patch->nbPhilopat++;
  
  return src_patch->get(SEX, ADLTx, RAND::Uniform( src_patch->size( SEX, ADLTx ) ) );
}
/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

//                     ******** LCE_Breed_Selection_Disperse ********/

// ----------------------------------------------------------------------------------------
// LCE_Breed_Selection_Disperse::LCE_Breed_Selection_Disperse
// ----------------------------------------------------------------------------------------
LCE_Breed_Selection_Disperse::LCE_Breed_Selection_Disperse()
:  LifeCycleEvent("breed_selection_disperse","")
{
  add_parameter("breed_selection_disperse_fitness_threshold", DBL, false, true, 0, 1, 
                new ParamUpdater<LCE_Breed_Selection_Disperse> (&LCE_Breed_Selection_Disperse::setParameters) );
}
// ----------------------------------------------------------------------------------------
// LCE_Breed_Selection_Disperse::setParameters
// ----------------------------------------------------------------------------------------
bool LCE_Breed_Selection_Disperse::setParameters()
{
  if(!LCE_Breed_Disperse::setParameters()) return false;
  if(!LCE_Selection_base::setParameters()) return false;
  
  if(get_parameter("breed_selection_disperse_fitness_threshold")->isSet())
    _base_fitness = get_parameter_value("breed_selection_disperse_fitness_threshold");
  else
    _base_fitness = 0.05;
  
  unsigned int model = this->getDispersalModel();
  if(model == 2)
    _breed_selection_disperse = &LCE_Breed_Selection_Disperse::breed_selection_disperse_propagule;
  else
    _breed_selection_disperse = &LCE_Breed_Selection_Disperse::breed_selection_disperse;
  
  return true;
}
// ----------------------------------------------------------------------------------------
// LCE_Breed_Selection::loadStatServices
// ----------------------------------------------------------------------------------------
void LCE_Breed_Selection_Disperse::loadStatServices ( StatServices* loader )
{
  LCE_Selection_base::loadStatServices(loader);
}
// ----------------------------------------------------------------------------------------
// LCE_Breed_Selection_Disperse::LCE_Breed_Selection_Disperse
// ----------------------------------------------------------------------------------------
void LCE_Breed_Selection_Disperse::execute()
{
  int ind_count = 0;
  
#ifdef _DEBUG_
  message("LCE_Breed_Selection_Disperse::execute (mean fitness(adlt): %f, fit scale: %f, offsprg nb: %i, adlt nb: %i, "
          ,getMeanFitness(ADLTx), _scaling_factor, _popPtr->size( OFFSPRG ), _popPtr->size( ADULTS ));
  fflush(stdout);
#endif
  
  //  _scaling_factor = (_mean_fitness >= _base_fitness ? 1 : _base_fitness/_mean_fitness);
  
  if(_npatch != _popPtr->getPatchNbr()) {
    _npatch = _popPtr->getPatchNbr();
    if(!updateDispMatrix()) fatal("bailing out\n");
  }
  
  LCE_Selection_base::resetCounters();  
  LCE_Disperse_base::reset_counters();                                                  
  
  (this->*_breed_selection_disperse)(&ind_count);
  
  LCE_Selection_base::setMeans(ind_count);
  
#ifdef _DEBUG_
  unsigned int a = 0, b = 0, c = 0;
  for(unsigned int i = 0; i < _popPtr->getPatchNbr(); i++){
    a += _popPtr->getPatch(i)->nbEmigrant;
    b += _popPtr->getPatch(i)->nbImigrant;
    c += _popPtr->getPatch(i)->nbPhilopat;
  }
  message("immigrate: %f, emigrants: %i, imigrants: %i, mean fitness(offsprg): %f)\n",(double)b/(b+c), a, b, _mean_fitness);
#endif
}  
// ----------------------------------------------------------------------------------------
// LCE_Breed_Selection_Disperse::breed_selection_disperse
// ----------------------------------------------------------------------------------------
void LCE_Breed_Selection_Disperse::breed_selection_disperse (int* counter)
{
  Patch* patch;
  unsigned int nfem, nmal;
  unsigned int patchNbr = _popPtr->getPatchNbr();
  
  for(unsigned int i = 0; i < patchNbr; i++) {
    
    patch = _popPtr->getPatch(i);
    
    if(patch->size(OFFSx) != 0) patch->flush(OFFSx, _popPtr);
    
    _mean_fitness = getMeanPatchFitness(ADLTx, i);
    _scaling_factor = (_mean_fitness >= _base_fitness ? 1 : _base_fitness/_mean_fitness);
    
    nfem = (this->* _get_numFemOffspring)(patch);
    nmal = (this->* _get_numMalOffspring)(patch);
    
    do_breed(FEM, nfem, counter, patch, i);
    do_breed(MAL, nmal, counter, patch, i);
    
  } 
}
// ----------------------------------------------------------------------------------------
// LCE_Breed_Selection_Disperse::do_breed
// ----------------------------------------------------------------------------------------
void LCE_Breed_Selection_Disperse::breed_selection_disperse_propagule (int* counter)
{
  this->setIsland_PropagulePool_Matrix();
  breed_selection_disperse(counter);
}
// ----------------------------------------------------------------------------------------
// LCE_Breed_Selection_Disperse::do_breed
// ----------------------------------------------------------------------------------------
void LCE_Breed_Selection_Disperse::do_breed(sex_t SEX, unsigned int size, int* counter, Patch* patch,
                                            unsigned int cpatch)
{
  //  int viab_count;
  Individual *new_ind;
  for(unsigned int j = 0; j < size; j++) {
    //    viab_count = (*counter);
    do{
      new_ind = makeOffspringWithSelection( (this->*_make_offspring)(SEX, patch, cpatch), cpatch);
      (*counter)++;
      //      if((*counter) - viab_count > _max_try) {
      //        warning("it took more than 500 trials to get a fit offspring.");
      //      }
    } while (new_ind == NULL);
    
    patch->add(SEX, OFFSx, new_ind);
  }
}
/* _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/ */

//                        ******** LCE_Breed_Selection ********/

// ----------------------------------------------------------------------------------------
// LCE_Breed_Selection
// ----------------------------------------------------------------------------------------
LCE_Breed_Selection::LCE_Breed_Selection ( ) : LifeCycleEvent("breed_selection",""), _breed_selection(0)
{
  add_parameter("breed_selection_fecundity_fitness", BOOL, false, false, 0, 0, 0);
}
// ----------------------------------------------------------------------------------------
// LCE_Breed_Selection::setParameters
// ----------------------------------------------------------------------------------------
bool LCE_Breed_Selection::setParameters ( )
{
  if( !LCE_Breed_base::setParameters( ) ||  !LCE_Selection_base::setParameters( ) )
    return false;
  
  if(getMeanFecundity() <= 0) {
    error("\"mean_fecundity\" is not set!\n");
    return false;
  }
  
  if(get_parameter("breed_selection_fecundity_fitness")->isSet()){
    
    _breed_selection = &LCE_Breed_Selection::do_breed_selection_FecFitness ;
    
  } else {
    _breed_selection = &LCE_Breed_Selection::do_breed_selection_OffSurvival ;
  }

  
  return true;
}
// ----------------------------------------------------------------------------------------
// LCE_Breed_Selection::loadStatServices
// ----------------------------------------------------------------------------------------
void LCE_Breed_Selection::loadStatServices ( StatServices* loader )
{
  LCE_Selection_base::loadStatServices(loader);
}
// ----------------------------------------------------------------------------------------
// LCE_Breed_Selection::breed
// ----------------------------------------------------------------------------------------
void LCE_Breed_Selection::execute ()
{
  Patch* current_patch;
  unsigned int ind_count;
  
  if(_popPtr->size(OFFSPRG) != 0) {
    warning("offspring containers not empty at time of breeding, flushing.\n");
    _popPtr->flush(OFFSx);
  }
  
#ifdef _DEBUG_
  message("LCE_Breed_Selection::execute\n");
#endif
  
  _scaling_factor = 1;
  _mean_fitness = 0;
  
  LCE_Selection_base::resetCounters();
  
  for(unsigned int home = 0; home < _popPtr->getPatchNbr(); home++) {
    
    current_patch = _popPtr->getPatch(home);
    
    if( !checkMatingCondition(current_patch) ) continue;
    
    (this->* _breed_selection) (current_patch, home, &ind_count);
    
  }//end_for home
  
  LCE_Selection_base::setMeans(ind_count);
  
#ifdef _DEBUG_
  message("   (offsprg nb: %i adults nb: %i, mean fitness: %f)\n",_popPtr->size(OFFSPRG),
          _popPtr->size(ADULTS), _mean_fitness);
#endif
  
}
// ----------------------------------------------------------------------------------------
// LCE_Breed_Selection::do_breed_selection_FecFitness
// ----------------------------------------------------------------------------------------
void LCE_Breed_Selection::do_breed_selection_OffSurvival (Patch* patch, unsigned int patchID, unsigned int *cntr)
{
  unsigned int nbBaby;
  //double mean_fec;
  Individual* FatherPtr;
  Individual* MotherPtr;
  Individual* NewOffsprg;
  
  for(unsigned int size = patch->size(FEM, ADLTx), indexOfMother = 0;
      indexOfMother < size;
      indexOfMother++) {
    
    MotherPtr = patch->get(FEM, ADLTx, indexOfMother);
    
    nbBaby = (unsigned int)MotherPtr->setFecundity( getFecundity() );
    
    (*cntr) += nbBaby;
    //-----------------------------------------------------------------------
    while(nbBaby != 0) {
      
      FatherPtr = this->getFatherPtr(patch, MotherPtr, indexOfMother);
      
      NewOffsprg = makeOffspringWithSelection( LCE_Breed_base::do_breed(MotherPtr, FatherPtr, patchID), patchID );
      
      if(NewOffsprg != NULL)
        patch->add(NewOffsprg->getSex(), OFFSx, NewOffsprg);
      
      nbBaby--;
      
    }//_END__WHILE;
  }//end_for indexOfMother
  
}
// ----------------------------------------------------------------------------------------
// LCE_Breed_Selection::do_breed_selection_OffSurvival
// ----------------------------------------------------------------------------------------
void LCE_Breed_Selection::do_breed_selection_FecFitness (Patch* patch, unsigned int patchID, unsigned int *cntr)
{
  unsigned int nbBaby;
  //double mean_fec;
  Individual* FatherPtr;
  Individual* MotherPtr;
  Individual* NewOffsprg;
  unsigned int cat;

  for(unsigned int size = patch->size(FEM, ADLTx), indexOfMother = 0;
      indexOfMother < size;
      indexOfMother++) {
    
    (this->*_setScalingFactor)(ADLTx, patchID); //from within LCE_selection_base
    
    MotherPtr = patch->get(FEM, ADLTx, indexOfMother);
        
    nbBaby = (unsigned int)MotherPtr->setFecundity( getFecundity(getMeanFecundity()*getFitness(MotherPtr, patchID,0)) );
        
    (*cntr) += nbBaby;
    //-----------------------------------------------------------------------
    while(nbBaby != 0) {
      
      
      FatherPtr = this->getFatherPtr(patch, MotherPtr, indexOfMother);
      
      NewOffsprg = LCE_Breed_base::makeOffspring( LCE_Breed_base::do_breed(MotherPtr, FatherPtr, patchID) );
      
      patch->add(NewOffsprg->getSex(), OFFSx, NewOffsprg);
      
      nbBaby--;
      
      cat = NewOffsprg->getPedigreeClass(NewOffsprg->getMother(), NewOffsprg->getFather());
      
      _fitness[cat] += getFitness(NewOffsprg, patchID,0);
      _ind_cntr[cat]++;
      _survival[cat]++;
      
    }//_END__WHILE;
  }//end_for indexOfMother
  
}
// ----------------------------------------------------------------------------------------
// LCE_Breed_Selection::makeOffspringWithSelection
// ----------------------------------------------------------------------------------------
Individual* LCE_Breed_Selection::makeOffspringWithSelection (Individual* ind,
                                                             unsigned int natalpatch)
{
  if(!ind->getMother() || !ind->getFather()) fatal("found NULL parent ptr in offspring (LCE_Breed_Selection::makeOffspringWithSelection)\n");
  if(_LCELinkedTraitIndex != -1)
    ind->createTrait(_LCELinkedTraitIndex, LCE_Breed_base::doInheritance(), true);
  
  register double fitness = getFitness(ind, natalpatch,0);
  register unsigned int cat = ind->getPedigreeClass(ind->getMother(), ind->getFather());
  
  _fitness[cat] += fitness;
  _ind_cntr[cat]++;
  
  if(RAND::Uniform() > fitness * _scaling_factor ) {
    //this one dies
    _popPtr->recycle(ind);
    
//    ind->getMother()->addMating(cat);
//    if(cat!=4) ind->getFather()->addMating(cat);
    
    return NULL;
    
  } else {
    
    //update counters
    _survival[cat]++;
    
//    ind->getMother()->DidHaveABaby(cat);
//    if(cat!=4) ind->getFather()->DidHaveABaby(cat);
  }
  
  //compute inheritance and mutation of the other traits:
  for(int i = 0; i < (int)_nb_trait; i++)
    if(i != _LCELinkedTraitIndex)
      ind->createTrait(i, LCE_Breed_base::doInheritance(), true);
  
  return ind;
}


