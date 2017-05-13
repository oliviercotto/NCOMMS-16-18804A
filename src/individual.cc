/**  $Id: individual.cc,v 1.5.2.3 2014-05-01 07:51:16 fred Exp $
*
*  @file individual.cc
*  NEMO
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
*  Created on @date 05.08.2004
*
*  @author fred
*/

#include <iostream>
#include <string.h>
#include "individual.h"
#include "ttrait.h"
#include "output.h"
#include "binarystoragebuffer.cc"

unsigned long Individual::currentID = 0;

// ----------------------------------------------------------------------------------------
// individual
// ----------------------------------------------------------------------------------------
Individual::Individual ( ) 
: _age(0), _sex(MAL),_mother(NULL),_father(NULL),
_home(0),_fecundity(0),_trait_nb(0)  //_motherID(0),_fatherID(0),   //MOD SLIM NEMO 7.11
{
  _id = currentID++;
//  _pedigreeClass == 0;   //MOD SLIM NEMO 7.11
//  for(unsigned int i = 0; i < 5; i++) {
//    _matings[i] = 0; _realizedFecundity[i] = 0;
//  }
}
// ----------------------------------------------------------------------------------------
// init
// ----------------------------------------------------------------------------------------
Individual * Individual::init ()
{
  _age = 0;
//  _pedigreeClass = 0;  //MOD SLIM NEMO 7.11
//  _motherID = 0;
//  _fatherID = 0;
  _mother = NULL;
  _father = NULL;
  _home = 0;
//  for(unsigned int i = 0; i < 5; i++) {
//    _matings[i] = 0; _realizedFecundity[i] = 0;  //MOD SLIM NEMO 7.11
//  }
  
  if(_trait_nb != Traits.size()){
    error("Individual::init: trait counter and table size differ, resetting\n");
    _trait_nb = Traits.size();
  }
  
  for(unsigned int i = 0; i < _trait_nb; i++)
    Traits[i]->init();

  return this;
}
// ----------------------------------------------------------------------------------------
// reset
// ----------------------------------------------------------------------------------------
void Individual::reset ()
{
  _id = 0;
  _age = 0;
  _sex = MAL;
//  _pedigreeClass = 0;  //MOD SLIM NEMO 7.11
//  _motherID = 0;  //MOD SLIM NEMO 7.11
//  _fatherID = 0;  //MOD SLIM NEMO 7.11
  _mother = NULL;
  _father = NULL;
  _home = 0;
//  for(unsigned int i = 0; i < 5; i++) {
//    _matings[i] = 0; _realizedFecundity[i] = 0;  //MOD SLIM NEMO 7.11
//  }
   
  if(_trait_nb != Traits.size()){
    warning("Individual::reset: trait counter and table size differ, resetting\n");
    _trait_nb = Traits.size();
  }
}
// ----------------------------------------------------------------------------------------
// store_data
// ----------------------------------------------------------------------------------------
void Individual::store_data ( BinaryStorageBuffer* saver )
{
  saver->store(&_id, sizeof(unsigned long));
  saver->store(&_age, sizeof(unsigned short)); //not in previous versions
  unsigned long dum = 0;   //MOD SLIM NEMO 7.11
  saver->store(&dum, sizeof(unsigned long));
  saver->store(&dum, sizeof(unsigned long));
//  saver->store(&_motherID, sizeof(unsigned long));
//  saver->store(&_fatherID, sizeof(unsigned long));
  saver->store(&_sex, sizeof(sex_t));
  saver->store(&_home, sizeof(unsigned short));
  unsigned short dummy[5] = {0,0,0,0,0};  // MOD NEMO SLIM 7.11.16
  saver->store(&dummy[0], 5*sizeof(unsigned short));
  saver->store(&dummy[0], 5*sizeof(unsigned short));
  saver->store(&dummy[0], 1);
}
// ----------------------------------------------------------------------------------------
// retrieve_data
// ----------------------------------------------------------------------------------------
void Individual::retrieve_data ( BinaryStorageBuffer* reader )
{
  reader->read(&_id, sizeof(unsigned long));
  reader->read(&_age, sizeof(unsigned short));
  unsigned long dum;   //MOD SLIM NEMO 7.11
  reader->read(&dum, sizeof(unsigned long));
  reader->read(&dum, sizeof(unsigned long));
//  reader->read(&_motherID, sizeof(unsigned long));
//  reader->read(&_fatherID, sizeof(unsigned long));
  reader->read(&_sex, sizeof(sex_t));
  reader->read(&_home, sizeof(unsigned short));
  unsigned short dummy[5] = {0,0,0,0,0};  // MOD NEMO SLIM 7.11.16
  reader->read(&dummy, 5*sizeof(unsigned short));
  reader->read(&dummy, 5*sizeof(unsigned short));
  reader->read(&dummy, 1);
}
// ----------------------------------------------------------------------------------------
// show_up
// ----------------------------------------------------------------------------------------
void Individual::show_up ()
{
  message("\n Individual ID: %i\n\
           age: %i\n\
           sex: %i\n\
        mother: 0\n\
        father: 0\n\
pedigree class: NA\n\
          home: %i\n\
 traits values: \n",_id,_age,_sex, _home); //_motherID,_fatherID,_pedigreeClass, //MOD SLIM NEMO 7.11

for(unsigned int i = 0; i < _trait_nb; i++)
  Traits[i]->show_up();
}
// ----------------------------------------------------------------------------------------
// clone
// ----------------------------------------------------------------------------------------
Individual* Individual::clone ()
{
  Individual* myClone = new Individual();
  
  for(unsigned int i = 0; i < _trait_nb; i++)
    myClone->addTrait(Traits[i]->clone(), i);
  
  return myClone;
}
// ----------------------------------------------------------------------------------------
// operator=
// ----------------------------------------------------------------------------------------
Individual& Individual::operator=(const Individual& i)
{
  if(this != &i) {  
    
    if(Traits.size() != i.Traits.size())
      fatal("Individual::operator=:not same number of traits in left and right sides of assignment\n");

    if(_trait_nb != i._trait_nb) {
      error("Individual::operator=:trait counters differ, resetting\n");
      _trait_nb = i._trait_nb;
    }
    
    _sex = i._sex;
    _age = i._age;
//    _motherID = i._motherID;
//    _fatherID = i._fatherID;
    _mother = i._mother;
    _father = i._father;
    _home = i._home;
//    _pedigreeClass = i._pedigreeClass;  //MOD SLIM NEMO 7.11
    _fecundity = i._fecundity;
//    for(unsigned int j = 0; j < 5; j++) {  //MOD SLIM NEMO 7.11
//      _matings[j] = i._matings[j];
//      _realizedFecundity[j] = i._realizedFecundity[j];
//    }

    for(unsigned int t = 0; t < _trait_nb; t++){ 
      if(Traits[t]->get_type().compare(i.Traits[t]->get_type()) != 0) 
        fatal("Individual::operator=: not same kinds of traits on left and right sides of assignment\n");
      (*Traits[t]) = (*i.Traits[t]);
    }
  }
  return *this;
}
// ----------------------------------------------------------------------------------------
// operator==
// ----------------------------------------------------------------------------------------
bool Individual::operator==(const Individual& i)
{
  if(this != &i) {
    if(Traits.size() != i.Traits.size()) return false;

    for(unsigned int t = 0; t < Traits.size(); t++)
      if((*Traits[t]) != (*i.Traits[t])) return false;

  //if(_sex != i._sex) return false;
  }
  return true;
}
// ----------------------------------------------------------------------------------------
// operator!=
// ----------------------------------------------------------------------------------------
bool Individual::operator!=(const Individual& i)
{
  if(!((*this) == i))
    return true;
  else
    return false;
}
// ----------------------------------------------------------------------------------------
// getPedigreeClass
// ----------------------------------------------------------------------------------------
unsigned int Individual::getPedigreeClass (Individual* mother, Individual* father)
{
  if(mother == father) return 4; //selfed

  if(father == 0 && mother != 0) return 4; //selfed through cloning, with hermaphrodites

  if(father != 0 && mother == 0) return 4; //selfed through cloning, not hermaphrodites

  if(mother->getHome() != father->getHome()) return 0; //outbred between patch
  
  Individual* mm, *mf,*fm,*ff; //MOD SLIM NEMO 7.11
  //mother's parents:
  mm = mother;//->getMotherID();//MOD SLIM NEMO 7.11
  mf = mother;//->getFatherID();//MOD SLIM NEMO 7.11
  //father's parents:
  fm = father;//->getMotherID();//MOD SLIM NEMO 7.11
  ff = father;//->getFatherID();//MOD SLIM NEMO 7.11

  if(mm != fm && mf != ff) return 1; //outbred within patch
    
  else if((mm == fm && mf != ff) || (mm != fm && mf == ff)) return 2; //half sibs
    
  else return 3; //full sibs
}

