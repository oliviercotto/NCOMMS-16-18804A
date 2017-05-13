/** $Id: indfactory.cc,v 1.10.2.2 2014-05-01 13:34:42 fred Exp $
*
*  @file indfactory.cc
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
*  created on @date 23.11.2005
* 
*  @author fred
*/

#include "indfactory.h"
#include "Uniform.h"
#include "output.h"

using namespace std;

 IndFactory::~IndFactory ( )
{
  purgeRecyclingPOOL(); 
  _protoIndividual.clearTraits();
  _protoTraits.clear(); //the proto traits are deallocated within ~SimBuilder
  _TraitsIndex.clear();
}
// ----------------------------------------------------------------------------------------
// makePrototype
// ----------------------------------------------------------------------------------------
void IndFactory::makePrototype(map< trait_t, TraitPrototype* > TTlist)
{
#ifdef _DEBUG_
  message("IndFactory::makePrototype\n");
#endif

  //first, reset the ID counter:
  _protoIndividual.setCurrentID(0);
  
  //store the traits list for future use:
  _protoTraits = TTlist;

  //then add the traits from the trait prototypes:
  map< trait_t, TraitPrototype* >::iterator trait = TTlist.begin();

  _TraitsIndex.clear(); //kept here but also done in clearPrototype()
  
  int i = 0;
  
  while(trait != TTlist.end()) {
    
#ifdef _DEBUG_
    message("IndFactory::makePrototype::addTrait: %s\n",trait->first.c_str());
    cout << "trait proto ptr:"<<trait->second <<endl;
    cout << "trait proto type:"<<trait->second->get_type() <<endl;
//    TTrait* tt=trait->second->hatch();
//    tt->init();
//    tt->init_sequence();
//    tt->show_up();
//    delete tt;
#endif
    
    //init the trait prototype, set attribute from input file parameters value
    if(!trait->second->setParameters()) {
      error("initialization of prototype for trait \"%s\" failed\n",trait->first.c_str());
      fatal("bailing out from individual prototype initialization\n");
    }
    //set the prototype's index value
    trait->second->set_index(i);

    //create the trait and add it to the proto-individual
    _protoIndividual.addTrait(trait->second->hatch(),i);
    
    //store the trait's index into the index map
    _TraitsIndex[trait->first] = i++;
    
    trait++;
  }
  
}
// ----------------------------------------------------------------------------------------
// makePrototype
// ----------------------------------------------------------------------------------------
void IndFactory::clearPrototype ( )
{
#ifdef _DEBUG_
  message("IndFactory::clearPrototype\n");
#endif
  
  //clear the ttraits:
  _protoIndividual.clearTraits();
  
  map< trait_t, TraitPrototype* >::iterator trait = _protoTraits.begin();

  while(trait != _protoTraits.end()) {
     
    trait->second->reset();
    
    trait++;
  }
  
  _protoTraits.clear();
  _TraitsIndex.clear();
  TTProtoWithMap::_map.clear();
}
// ----------------------------------------------------------------------------------------
// getTraitIndex
// ----------------------------------------------------------------------------------------
int IndFactory::getTraitIndex (trait_t type)
{
  map< trait_t, int >::iterator trait = _TraitsIndex.find(type);
  
  if(trait == _TraitsIndex.end())
    return -1;
  
  return trait->second;
}
// ----------------------------------------------------------------------------------------
// getTraitPrototype
// ----------------------------------------------------------------------------------------
TraitPrototype* IndFactory::getTraitPrototype (trait_t type)
{
  map< trait_t, TraitPrototype* >::iterator trait = _protoTraits.find(type);
  
  if(trait == _protoTraits.end()) return NULL;
  
  return trait->second;
  
}
// ----------------------------------------------------------------------------------------
// makeNewIndividual
// ----------------------------------------------------------------------------------------
Individual* IndFactory::makeNewIndividual(Individual* mother, Individual* father, sex_t sex, unsigned short homepatch)
{
  Individual* newind;
  
  if(RecyclingPOOL.empty()) {
    //create new Individual by cloning the proto-individual
    newind = _protoIndividual.clone();
    //allocate memory for the traits' sequences:
    newind->init();
  } else {
    //recycle an Individual from the POOL
    newind = RecyclingPOOL[0];
    RecyclingPOOL.pop_front();
    newind->reset();
    newind->setID(Individual::currentID++);
  }
  
  newind->setSex(sex);
  newind->setHome(homepatch);
  if(mother != NULL && father != NULL) {
//    newind->setFatherID(father->getID());//MOD SLIM NEMO 7.11
//    newind->setMotherID(mother->getID());//MOD SLIM NEMO 7.11
    newind->setFather(father);
    newind->setMother(mother);
//    newind->setPedigreeClass(mother, father);  //MOD SLIM NEMO 7.11
  }
  
  return newind;
}
// ----------------------------------------------------------------------------------------
// makeOffsprg
// ----------------------------------------------------------------------------------------
Individual* IndFactory::makeOffsprg(Individual* mother, Individual* father, sex_t sex, unsigned short homepatch)
{  
  Individual* offspring = makeNewIndividual(mother,father,sex,homepatch);
//  unsigned int cat = offspring->getPedigreeClass();  //MOD SLIM NEMO 7.11
//  mother->DidHaveABaby(cat);
//  if(cat!=4) father->DidHaveABaby(cat);
  return offspring->create(mother, father);
}
