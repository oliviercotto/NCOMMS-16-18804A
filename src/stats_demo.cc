/** $Id: stats_demo.cc,v 1.7.2.1 2014-04-29 18:06:51 fred Exp $
*
*  @file stats_demo.cc
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
*   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA.
*
*  Created on @date 23.01.2004
*
*  @author fred
*/

#include <cmath>
#include "MPStatHandler.h"
#include "metapop.h"


using namespace std;

/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

//                         ******** Migrants analysis ********
// ----------------------------------------------------------------------------------------
// getMeanEmigrantPerPatch
// ----------------------------------------------------------------------------------------
double MPStatHandler::getMeanEmigrantPerPatch()
{
  unsigned int meanM = 0;

  for(unsigned int i = 0; i < _pop->getPatchNbr(); i++)
    meanM += _pop->getPatch(i)->nbEmigrant;
  meanEmigrant = (double)meanM/_pop->getPatchNbr();
  return meanEmigrant;
}
// ----------------------------------------------------------------------------------------
// getMeanImigrantPerPatch
// ----------------------------------------------------------------------------------------
double MPStatHandler::getMeanImigrantPerPatch()
{
  double meanM = 0,nbpatch = 0;
  Patch* current_patch;
  for(unsigned int i = 0; i < _pop->getPatchNbr(); i++) {
    current_patch = _pop->getPatch(i);
    if(!current_patch->get_isExtinct()){
      meanM += current_patch->nbImigrant;
      nbpatch++;
    }
  }
  meanImigrant = (nbpatch != 0 ? meanM/nbpatch : nanf("NULL"));
  return meanImigrant;
}
// ----------------------------------------------------------------------------------------
// getMeanMigrantRatio
// ----------------------------------------------------------------------------------------
double MPStatHandler::getMeanMigrantRatio()
{
  return ((meanImigrant+meanResidant) != 0 ? meanImigrant/(meanImigrant+meanResidant) : 0);
}
// ----------------------------------------------------------------------------------------
// getMeanResidantPerPatch
// ----------------------------------------------------------------------------------------
double MPStatHandler::getMeanResidantPerPatch()
{
  double meanR = 0,nbpatch = 0;
  Patch* current_patch;
  for(unsigned int i = 0; i < _pop->getPatchNbr(); i++){
    current_patch = _pop->getPatch(i);
    if(!current_patch->get_isExtinct()){
      meanR += current_patch->nbPhilopat;
      nbpatch++;
    }
  }
  meanResidant = (nbpatch != 0 ? meanR/nbpatch : nanf("NULL"));
  return meanResidant;
}
// ----------------------------------------------------------------------------------------
// getMeanMigrantProportion
// ----------------------------------------------------------------------------------------
double MPStatHandler::getMeanKolonisersProportion()
{
  double nbpatch = 0, mean = 0;
  Patch* current_patch;
  for(unsigned int i = 0; i < _pop->getPatchNbr(); i++){
    current_patch = _pop->getPatch(i);
    if(current_patch->nbKolonisers >= 0) {
      mean += fmin((double)current_patch->nbKolonisers / current_patch->get_K(), 1.0);
      nbpatch++;
    }
  }

  return  (nbpatch != 0 ? mean/nbpatch : nanf("NULL"));
}
// ----------------------------------------------------------------------------------------
// getMeanKolonisersPerPatch
// ----------------------------------------------------------------------------------------
double MPStatHandler::getMeanKolonisersPerPatch()
{
  double mean = 0,nbpatch = 0;
  Patch* current_patch;
  for(unsigned int i = 0; i < _pop->getPatchNbr(); i++){
    current_patch = _pop->getPatch(i);
    if(current_patch->nbKolonisers >= 0) {
      mean += current_patch->nbKolonisers;
      nbpatch++;
    }
  }
  meanKolonisers = (nbpatch != 0 ? mean/nbpatch : nanf("NULL"));
  return meanKolonisers;
}

double MPStatHandler::getEmigrantInPatch  (unsigned int i)
{
  Patch* patch = _pop->getPatch(i);
  return (patch != 0 ? patch->nbEmigrant : 0);
}

double MPStatHandler::getResidantInPatch  (unsigned int i)
{
  Patch* patch = _pop->getPatch(i);
  return (patch != 0 ? patch->nbPhilopat : 0);
}

double MPStatHandler::getImigrateInPatch  (unsigned int i)
{
  Patch* patch = _pop->getPatch(i);
  return (patch != 0 ? (double)patch->nbImigrant/(double)(patch->nbImigrant + patch->nbPhilopat)
          : nanf("NULL"));
}


double MPStatHandler::getKolonisersInPatch(unsigned int i) 
{
  Patch* patch = _pop->getPatch(i);
  int colon = (patch != 0 ? patch->nbKolonisers : 0);
  return (colon != -1 ? colon : nanf("NULL"));
}

/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

//                     ****** Patch extinction analysis ********

// ----------------------------------------------------------------------------------------
// setObsrvdExtinctionRate
// ----------------------------------------------------------------------------------------
void MPStatHandler::setObsrvdExtinctionRate ()
{
  ObservedExtinctionRate = 0;
  
  for(unsigned int i = 0; i < _pop->getPatchNbr(); ++i)
	ObservedExtinctionRate += _pop->getPatch(i)->isEmpty();
  
  ObservedExtinctionRate /= _pop->getPatchNbr();
}

double MPStatHandler::getPatchAge (unsigned int i) 
{
  Patch* patch = _pop->getPatch(i);
  return (patch != 0 ? patch->get_age() : 0);
}

double MPStatHandler::getMeanPatchAge () 
{
  int mean = 0;
  for(unsigned int i = 0; i < _pop->getPatchNbr(); ++i)
    mean += _pop->getPatch(i)->get_age();
  return mean/_pop->getPatchNbr();
}

/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

//                   ****** Population saturation analysis ********

// ----------------------------------------------------------------------------------------
// getMeanOffsprgFemNumber
// ----------------------------------------------------------------------------------------
double MPStatHandler::getMeanOffsprgFemNumber ()
{
  double mean = 0;
  unsigned int nb_patch = 0;

  for(unsigned int i = 0; i < _pop->getPatchNbr(); i++) {
    nb_patch += (_pop->size(FEM,OFFSPRG,i) != 0);
    mean += _pop->size(FEM,OFFSPRG,i);
  }

  return (nb_patch != 0 ? mean/nb_patch : nanf("NULL"));
}
// ----------------------------------------------------------------------------------------
// getMeanOffsprgMalNumber
// ----------------------------------------------------------------------------------------
double MPStatHandler::getMeanOffsprgMalNumber ()
{
  double mean = 0;
  unsigned int nb_patch = 0;

  for(unsigned int i = 0; i < _pop->getPatchNbr(); i++) {
    nb_patch += (_pop->size(MAL,OFFSPRG,i) != 0);
    mean += _pop->size(MAL,OFFSPRG,i);
  }

  return (nb_patch != 0 ? mean/nb_patch : nanf("NULL"));
}
// ----------------------------------------------------------------------------------------
// getMeanOffsprgPopSatLevel
// ----------------------------------------------------------------------------------------
double MPStatHandler::getMeanOffsprgPopSatLevel()
{
  double mean = 0;
  unsigned int nb_patch = 0;
  for(unsigned int i = 0; i < _pop->getPatchNbr(); i++) {
    nb_patch += (_pop->size(OFFSPRG,i) != 0);
    mean += _pop->getPatch(i)->getPatchOffsprgSatLevel();
  }

  return (nb_patch != 0 ? mean/nb_patch : nanf("NULL"));
}
// ----------------------------------------------------------------------------------------
// getMeanOffsprgPopSatVar
// ----------------------------------------------------------------------------------------
double MPStatHandler::getMeanOffsprgPopSatVar()
{
  unsigned int i,nb_patch;
  vector<double> stat;
  double var = 0, mean = 0, val;

  for(i = 0; i < _pop->getPatchNbr(); i++) {
    if(_pop->size(OFFSPRG,i) != 0){
      val = _pop->getPatch(i)->getPatchOffsprgSatLevel();
      stat.push_back(val);
      mean += val;
    }
  }
  nb_patch = stat.size();
  mean = (nb_patch != 0 ? mean/nb_patch : 0);

  for(i = 0; i < nb_patch; i++)
    var += pow((stat[i]-mean),2);

  return (nb_patch != 0 ? var/nb_patch : nanf("NULL"));
}
// ----------------------------------------------------------------------------------------
// getMeanAdultsFemNumber
// ----------------------------------------------------------------------------------------
double MPStatHandler::getMeanAdultsFemNumber()
{
  double mean = 0;
  unsigned int nb_patch = 0;
  for(unsigned int i = 0; i < _pop->getPatchNbr(); i++) {
    nb_patch += (_pop->size(FEM,ADULTS,i) != 0);
    mean += _pop->size(FEM,ADULTS,i);
  }
  return (nb_patch != 0 ? mean/nb_patch : nanf("NULL"));
}
// ----------------------------------------------------------------------------------------
// getMeanAdultsMalNumber
// ----------------------------------------------------------------------------------------
double MPStatHandler::getMeanAdultsMalNumber()
{
  double mean = 0;
  unsigned int nb_patch = 0;
  for(unsigned int i = 0; i < _pop->getPatchNbr(); i++) {
    nb_patch += (_pop->size(MAL,ADULTS,i) != 0);
    mean += _pop->size(MAL,ADULTS,i);
  }
  return (nb_patch != 0 ? mean/nb_patch : nanf("NULL"));
}
// ----------------------------------------------------------------------------------------
// getMeanAdultsPopSatLevel
// ----------------------------------------------------------------------------------------
double MPStatHandler::getMeanAdultsPopSatLevel()
{
  double mean = 0;
  unsigned int nb_patch = 0;
  for(unsigned int i = 0; i < _pop->getPatchNbr(); i++) {
    nb_patch += (_pop->size(ADULTS,i) != 0);
    mean += _pop->getPatch(i)->getPatchAdultsSatLevel();
  }
  return (nb_patch != 0 ? mean/nb_patch : nanf("NULL"));
}
// ----------------------------------------------------------------------------------------
// getMeanAdultsPopSatVar
// ----------------------------------------------------------------------------------------
double MPStatHandler::getMeanAdultsPopSatVar()
{
  unsigned int i,nb_patch;
  vector<double> stat;
  double var = 0, mean = 0, val;

  for(i = 0; i < _pop->getPatchNbr(); i++) {
    if(_pop->size(ADULTS,i) != 0){
      val = _pop->getPatch(i)->getPatchAdultsSatLevel();
      stat.push_back(val);
      mean += val;
    }
  }

  nb_patch = stat.size();
  mean = (nb_patch != 0 ? mean/nb_patch : 0);

  for(i = 0; i < nb_patch; i++)
    var += pow((stat[i]-mean),2);

  return (nb_patch != 0 ? var/nb_patch : nanf("NULL"));
}
/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

//                   ****** Fecundity and matings analysis ********

// ----------------------------------------------------------------------------------------
// getMeanAssignedFecundity
// ----------------------------------------------------------------------------------------
double MPStatHandler::getMeanAssignedFecundity(unsigned int sex)
{
  double mean = 0, sum = 0;
  unsigned int nbpatch = 0;
  Patch* crnt_patch;

  if((bool)sex) {
    for(unsigned int i = 0; i < _pop->getPatchNbr(); ++i) {
	  crnt_patch = _pop->getPatch(i);
      if(crnt_patch->size(FEM, ADLTx) != 0) {
        nbpatch++;
        sum = 0;
        for(unsigned int j = 0; j < crnt_patch->size(FEM, ADLTx);++j)
          sum += crnt_patch->get(FEM, ADLTx, j)->getFecundity();
        mean += sum/crnt_patch->size(FEM, ADLTx);
      }
    }
  } else {
    for(unsigned int i = 0; i < _pop->getPatchNbr(); ++i) {
	  crnt_patch = _pop->getPatch(i);
      if(crnt_patch->size(MAL, ADLTx) != 0) {
        nbpatch++;
        sum = 0;
        for(unsigned int j = 0; j < crnt_patch->size(MAL, ADLTx);++j)
          sum += crnt_patch->get(MAL, ADLTx, j)->getFecundity();
        mean += sum/crnt_patch->size(MAL, ADLTx);
      }
    }
  }
  return (nbpatch != 0 ? mean/nbpatch : nanf("NULL"));
}
// ----------------------------------------------------------------------------------------
// getMeanMatings
// ----------------------------------------------------------------------------------------
double MPStatHandler::getMeanMatings(unsigned int sex)
{
  double mean = 0, sum = 0;
  unsigned int nbpatch = 0;
  Patch* crnt_patch;

  if((bool)sex) {
    for(unsigned int i = 0; i < _pop->getPatchNbr(); ++i) {
	  crnt_patch = _pop->getPatch(i);
      if(crnt_patch->size(FEM, ADLTx) != 0) {
        nbpatch++;
        sum = 0;
        for(unsigned int j = 0; j < crnt_patch->size(FEM, ADLTx);++j)
//          sum += crnt_patch->get(FEM, ADLTx, j)->getTotMatings(); //MOD SLIM NEMO 7.11
        mean += sum/crnt_patch->size(FEM, ADLTx);
      }
    }
  } else {
    for(unsigned int i = 0; i < _pop->getPatchNbr(); ++i) {
	  crnt_patch = _pop->getPatch(i);
      if(crnt_patch->size(MAL, ADLTx) != 0) {
        nbpatch++;
        sum = 0;
        for(unsigned int j = 0; j < crnt_patch->size(MAL, ADLTx);++j)
//          sum += crnt_patch->get(MAL, ADLTx, j)->getTotMatings(); //MOD SLIM NEMO 7.11
        mean += sum/crnt_patch->size(MAL, ADLTx);
      }
    }
  }
  return (nbpatch != 0 ? mean/nbpatch : nanf("NULL"));
}
// ----------------------------------------------------------------------------------------
// getMeanRealizedFecundity
// ----------------------------------------------------------------------------------------
double MPStatHandler::setReproductiveStats(unsigned int sex)
{
  double var = 0, mean = 0;
  unsigned int tot_size = _pop->size( (sex_t)sex, ADULTS ), v = 0;
  double *stat;
  Patch *crnt_patch;

  if(tot_size == 0) {
    _var_reprod_success = nanf("NULL");
    return nanf("NULL");
  }
  
  stat = new double [tot_size];
  
  //females
  if((bool)sex) {
    
    for(unsigned int i = 0; i < _pop->getPatchNbr(); ++i) {
	  crnt_patch = _pop->getPatch(i);
      for(unsigned int j = 0, size = crnt_patch->size(FEM, ADLTx);
          j < size;
          ++j)
      {
//        stat[v] = crnt_patch->get(FEM, ADLTx, j)->getTotRealizedFecundity(); //MOD SLIM NEMO 7.11
        mean += stat[v];
        v++;
      }
    }
    
    mean /= tot_size;
    
    for(unsigned int i = 0; i < tot_size; ++i)
      var += pow( (stat[i] - mean), 2);
   
    var /= tot_size;
  //males
  } else {
    
    for(unsigned int i = 0; i < _pop->getPatchNbr(); ++i) {
	  crnt_patch = _pop->getPatch(i);
      for(unsigned int j = 0, size = crnt_patch->size(MAL, ADLTx);
          j < size;
          ++j)
      {
//        stat[v] = crnt_patch->get(MAL, ADLTx, j)->getTotRealizedFecundity(); //MOD SLIM NEMO 7.11
        mean += stat[v];
        v++;
      }
    }
    
    mean /= tot_size;
    
    for(unsigned int i = 0; i < tot_size; ++i)
      var += pow( (stat[i] - mean), 2);
    
    var /= tot_size;
    
  }
  
  _var_reprod_success = var;
  
  delete [] stat;
  
  return mean;
}
/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

//                             ****** kinship analysis ********

// ----------------------------------------------------------------------------------------
// setKinship
// ----------------------------------------------------------------------------------------
void MPStatHandler::setKinship()
{
  unsigned int i,j,k;
  unsigned int patchNbr = this->_pop->getPatchNbr(), Msize=0, Fsize=0;
  Patch* current_patch;
  Individual *I1,*I2;
  
  //counters initialization
  for(i = 0; i < 5; ++i) _sib_prop[i] = 0.0;
  
  for( i = 0; i < patchNbr; ++i) {
    
    current_patch = _pop->getPatch(i);
    
    //male-male
	if ( (Msize = current_patch->size(MAL, OFFSx)) != 0) {
      
	  for(j = 0; j < Msize -1; ++j) {
        
        I1 = current_patch->get(MAL, OFFSx, j);
        
		for(k = j+1; k < Msize; ++k) {
          
          I2 = current_patch->get(MAL, OFFSx, k);
          
          setKinClassCounter(I1, I2);
          
		} //end for k < size
        
        //selfed offspring counter:
//        if(I1->getIsSelfed()) _sib_prop[4]++; //MOD SLIM NEMO 7.11
                
	  } //end for j < size-1    
      
      //don't forget the last one!
//      if(current_patch->get(MAL, OFFSx, Msize -1)->getIsSelfed()) _sib_prop[4]++; //MOD SLIM NEMO 7.11
      
	}//endif
    
    //female-female
	if ( (Fsize = current_patch->size(FEM, OFFSx)) != 0) {
      
	  for(j = 0; j < Fsize -1; ++j) {
        
        I1 = current_patch->get(FEM, OFFSx, j);
        
		for(k = j+1; k < Fsize; ++k) {
          
          I2 = current_patch->get(FEM, OFFSx, k);
          
          setKinClassCounter(I1, I2);
          
		} //end for k < size  
        
          //selfed offspring counter:
//        if(I1->getIsSelfed()) _sib_prop[4]++;
        
	  } //end for j < size-1
      
      //don't forget the last one!
//      if(current_patch->get(FEM, OFFSx, Fsize -1)->getIsSelfed()) _sib_prop[4]++; //MOD SLIM NEMO 7.11

	}//endif
    
    //male-female
    for(j = 0; j < Msize; ++j) {
      
      I1 = current_patch->get(MAL, OFFSx, j);
      
      for(k = 0; k < Fsize; ++k) {
        
        I2 = current_patch->get(FEM, OFFSx, k);
        
        setKinClassCounter(I1, I2);        
      } //end for k
    } //end for j

  } //end for i < patchNbr
  
  //total number of pairwise comparisons:
  double tot = _sib_prop[0] + _sib_prop[1] + _sib_prop[2] + _sib_prop[3];
  
  for(i = 0 ; i < 4; ++i) _sib_prop[i] /= tot;

  _sib_prop[4] /= _pop->size(OFFSPRG);
}
// ----------------------------------------------------------------------------------------
// setKinClassCounter
// ----------------------------------------------------------------------------------------
void MPStatHandler::setKinClassCounter(Individual *I1, Individual *I2)
{

  //MOD SLIM NEMO 7.11//MOD SLIM NEMO 7.11//MOD SLIM NEMO 7.11//MOD SLIM NEMO 7.11//MOD SLIM NEMO 7.11
  //non sibs
//  if((I1->getMotherID() != I2->getMotherID()) && (I1->getFatherID() != I2->getFatherID()))
//    _sib_prop[0]++;
//
//  //maternal half sibs
//  else if((I1->getMotherID() == I2->getMotherID()) && (I1->getFatherID() != I2->getFatherID()))
//    _sib_prop[1]++;
//
//  //paternal half sibs
//  else if((I1->getMotherID() != I2->getMotherID()) && (I1->getFatherID() == I2->getFatherID()))
//    _sib_prop[2]++;
//
//  //full sibs
//  else if((I1->getMotherID() == I2->getMotherID()) && (I1->getFatherID() == I2->getFatherID()))
//    _sib_prop[3]++;
}  

// ----------------------------------------------------------------------------------------
// setPedegreeCount
// ----------------------------------------------------------------------------------------
void MPStatHandler::setPedegreeCount()
{
  
  unsigned int i,j;
  unsigned int patchNbr = _pop->getPatchNbr();
  Patch* patch;
  
  //counters initialization
  for(i = 0; i < 5; ++i) _ped_prop[i] = 0.0;
  
  for(i = 0; i < patchNbr; ++i) {
    
    patch = _pop->getPatch(i);
    
    //males
    for(j = 0; j < patch->size(MAL, OFFSx); ++j) {
      
//      _ped_prop[ patch->get(MAL, OFFSx, j)->getPedigreeClass() ]++; //MOD SLIM NEMO 7.11
      
    }  
    
    //females
    for(j = 0; j < patch->size(FEM, OFFSx); ++j) {
      
//      _ped_prop[ patch->get(FEM, OFFSx, j)->getPedigreeClass() ]++; //MOD SLIM NEMO 7.11
      
    }
    
    
  } //end for i < patchNbr
  
  //total:
  double tot = _ped_prop[0] + _ped_prop[1] + _ped_prop[2] + _ped_prop[3] + _ped_prop[4];
  
  for(i = 0 ; i < 5; ++i) _ped_prop[i] /= tot;
  
}

double MPStatHandler::getOffsprgSexRatio () {return (_pop->size(MAL, OFFSx)!= 0 ? (double)_pop->size(FEM, OFFSx)/_pop->size(MAL, OFFSx) : 0);}
double MPStatHandler::getAdultSexRatio   () {return (_pop->size(MAL, ADULTS) != 0 ? (double)_pop->size(FEM, ADULTS)/_pop->size(MAL, ADULTS) : 0);}
double MPStatHandler::getOffsprgNumber   () {return _pop->size(OFFSPRG);}
double MPStatHandler::getOffspringNumber (unsigned int i){return _pop->size(OFFSx, i);}
double MPStatHandler::getOffFemNumber    (unsigned int i){return _pop->size(FEM, OFFSx,i);}
double MPStatHandler::getOffMalNumber    (unsigned int i){return _pop->size(MAL, OFFSx,i);}
double MPStatHandler::getFemNumber       () {return _pop->size(FEM, ADULTS);}
double MPStatHandler::getMalNumber       () {return _pop->size(MAL, ADULTS);}
double MPStatHandler::getAdultsNumber    () {return _pop->size(ADULTS);}
double MPStatHandler::getAdultsNumber    (unsigned int i){return _pop->size(ADULTS, i);}
double MPStatHandler::getFemNumber       (unsigned int i){return _pop->size(FEM, ADULTS, i);}
double MPStatHandler::getMalNumber       (unsigned int i){return _pop->size(MAL, ADULTS, i);}
double MPStatHandler::getFemNumber       (unsigned int age, unsigned int i){return _pop->size(FEM, static_cast<age_idx> (age), i);}
double MPStatHandler::getMalNumber       (unsigned int age, unsigned int i){return _pop->size(MAL, static_cast<age_idx> (age), i);}

