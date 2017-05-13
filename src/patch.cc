/**  $Id: patch.cc,v 1.5.2.2 2014-04-29 18:20:32 fred Exp $
*
*  @file patch.cc
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
*  @author fred
*/

#include <iostream>
#include "metapop.h"
#include "Uniform.h"
#include "output.h"
#include "simenv.h"

using namespace std;

// ----------------------------------------------------------------------------------------
// init
// ----------------------------------------------------------------------------------------
Patch* Patch::init(unsigned int nb_age_classes, unsigned int nbfem, unsigned int nbmal, unsigned int id)
{
  _ID = id;
  _KFem = nbfem;
  _KMal = nbmal;
  _K = _KFem + _KMal;
  _isExtinct = false;
  _age = 0;

  if(_nb_age_class != nb_age_classes) reset_containers(nb_age_classes);
  
  _nb_age_class = nb_age_classes;
  
  //in case age structure didn't change and patch is not empty:
  if(size(ALL) != 0) flush(ALL, SIMenv::MainSim->get_pop()); 
  
  reset_counters();
  
  return this;
}
// ----------------------------------------------------------------------------------------
// reset_counters
// ----------------------------------------------------------------------------------------
void Patch::reset_counters()
{
  nbEmigrant = 0;
  nbImigrant = 0;
  nbPhilopat = 0;
  nbKolonisers = 0;
}
// ----------------------------------------------------------------------------------------
// reset_containers
// ----------------------------------------------------------------------------------------
void Patch::reset_containers(unsigned int new_nb_classes)
{  
  
  for(unsigned int i = 0; i < 2; i++) {
    
    if(_containers[i]) {
      
      flush(ALL, SIMenv::MainSim->get_pop());
      
      for(unsigned int j = 0; j < _nb_age_class; j++) {
        
        delete [] _containers[i][j];
      }
      
      delete [] _containers[i];
    }
    
    _containers[i] = new Individual** [new_nb_classes];
    
    if (_sizes[i]) {
      delete [] _sizes[i];
    }
     
    _sizes[i] = new unsigned int  [new_nb_classes];
    
    if (_capacities[i]) {
      delete [] _capacities[i];
    }
    
    _capacities[i] = new unsigned int  [new_nb_classes];
  }
    
  for(unsigned int i = 0; i < new_nb_classes; i++) {
    
    _containers[MAL][i] = new Individual* [ _KMal ];
    
    for(unsigned int j = 0; j < _KMal; ++j)
      _containers[MAL][i][j] = 0;
    
    _containers[FEM][i] = new Individual* [ _KFem ];
    
    for(unsigned int j = 0; j < _KFem; ++j)
      _containers[FEM][i][j] = 0;
    
    _sizes[MAL][i] = 0;
    _sizes[FEM][i] = 0;
    _capacities[MAL][i] = _KMal;
    _capacities[FEM][i] = _KFem;
  }
}
// ----------------------------------------------------------------------------------------
// setNewGeneration
// ----------------------------------------------------------------------------------------
void Patch::setNewGeneration(age_t AGE, Metapop* pop)
{
  unsigned int mask = 1;

  for(unsigned int i = 0; i < _nb_age_class; i++) {
    if( (mask & AGE) != 0) setNewGeneration(static_cast<age_idx>(i), pop);
    mask<<=1;
  }
}
// ----------------------------------------------------------------------------------------
// setNewGeneration
// ----------------------------------------------------------------------------------------
void Patch::setNewGeneration(age_idx AGE, Metapop* pop)
{  
  Individual *new_ind;
  //--------------------------------------------------------------------
  //if too much females in the Patch, flush them into the RecyclingPOOL
  if(size(FEM, AGE) > 0) flush(FEM, AGE, pop);
  
  for(unsigned int i = 0; i < _KFem; i++) {
    new_ind = pop->makeNewIndividual(0,0,FEM,_ID);
    new_ind->create_first_gen();
    new_ind->setAge(pop->getAgeStructure()->get(0, AGE));
    add(FEM, AGE, new_ind);
  }
  
  //--------------------------------------------------------------------
  //males: same as for the females....
  if(size(MAL, AGE) > 0) flush(MAL, AGE, pop);
  
  for(unsigned int i = 0; i < _KMal; i++) {
    new_ind = pop->makeNewIndividual(0,0,MAL,_ID);
    new_ind->create_first_gen();
    new_ind->setAge(pop->getAgeStructure()->get(0, AGE));
    add(MAL, AGE, new_ind);
  }
  
}
// ----------------------------------------------------------------------------------------
// ~Patch
// ----------------------------------------------------------------------------------------
Patch::~Patch()
{
#ifdef _DEBUG_
//  message("Patch::~Patch\n");
#endif
  
  for (unsigned int i = 0; i < 2; ++i) {
    for(unsigned int j = 0; j < _nb_age_class; ++j) {
      for(unsigned int k = 0; k < _sizes[i][j] ; ++k)
        if(_containers[i][j][k] != 0) delete _containers[i][j][k];
      delete [] _containers[i][j];
    }
    delete [] _containers[i];
    delete [] _sizes[i];
    delete [] _capacities[i];
  }
  

}
// ----------------------------------------------------------------------------------------
// show_up
// ----------------------------------------------------------------------------------------
void Patch::show_up()
{
  message("Patch %i:\n  age: %i; K: %i, K_fem: %i; K_mal: %i\n",_ID, _age, _K, _KFem, _KMal);
  for(unsigned int j = 0; j < _nb_age_class; ++j)
    message("   age class %i: females: %i; males: %i\n", j, _sizes[FEM][j], _sizes[MAL][j]);
}
