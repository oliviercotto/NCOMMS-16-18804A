/** $Id: indfactory.h,v 1.7.2.1 2014-04-29 18:27:44 fred Exp $
*
*  @file indfactory.h
*  Nemo2
*
*  Copyright (C) 2006-2011 Frederic Guillaume
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
*  Created on @date 20.08.2004
*  @author fred
*/

#ifndef FACTORY_H
#define FACTORY_H

#include <map>
#include <deque>
#include "individual.h"
#include "types.h"
#include "ttrait.h"

/**Factory of Individual, stores the individual prototype and the trait prototypes, manages the individual garbage collector.
 * Provides methods to generate new individuals within the Metapop. Each new individual is created by cloning a prototype itself
 * created at the simulation setup. New individuals are decorated with the appropriate traits as set by the trait prototypes and 
 * receives a unique ID (unique within a simulation).
 **/
class IndFactory {
protected:
  /**Map of the trait prototypes.*/
  std::map< trait_t,TraitPrototype* > _protoTraits;
  /**Table containing the index of each trait.*/
  std::map< trait_t, int > _TraitsIndex;

  /**The individuals prototype used to create any new individual in a simulation.*/
  Individual _protoIndividual;
  
  /**Garbage collector for unused Individual's.*/
  std::deque<Individual*> RecyclingPOOL;
  
public:
  
  IndFactory ( ) { }; 
  virtual ~IndFactory ( );
  
  /**Put an individual in the recycling pool.*/
  void recycle(Individual* ind) 
  { if(ind == NULL) error("IndFactory::recycle:ind is NULL!!\n"); else RecyclingPOOL.push_back(ind); }
  
  /**Empty the recycling pool.*/
  void purgeRecyclingPOOL ( )
  { for(unsigned int i=0; i < RecyclingPOOL.size(); ++i) delete RecyclingPOOL[i]; RecyclingPOOL.clear(); }
  
  /**Creates the individuals prototype from the selected trait prototypes.
    Resets the individual's ID counter to 0 and sets the traits index table.
    @param TTlist the list of the current trait prototype selected from the current simulation parameters.
    **/
  void                    makePrototype               (map< trait_t,TraitPrototype* > TTlist);
  
  /**Reset the trait prototypes, mostly done to unregister the genetic maps.*/
  void clearPrototype ( );
 
  /**Creates a blank individual which has to be "decorated" later.
   * ID is set and new traits are allocated but no genetic data is created (i.e. TTrait::init_sequence() is not called).
   * Sex has to be set later too.
   **/
  Individual*             getNewIndividual() {return makeNewIndividual(NULL,NULL,MAL,0);}
  
  /**Creates an individual with pointers to parents, sex and home ID set but no genetic data.
   * No inheritance or mutations on the trait sequences are done.
   * Sets the pedigree class of the individual.
   * Calls Individual::init() to allocate traits' sequences memory if individual is cloned from the prototype. Otherwise,
   * calls Individual::reset() if the new individual is coming from the recycling pool.
   * @param mother ptr to the mother
   * @param father ptr to the father
   * @param sex gender of the individual
   * @param homepatch ID of the Patch where this individual is born, usually the current position in the Metapop::vPatch array.
   **/
  Individual*             makeNewIndividual           (Individual* mother, Individual* father, sex_t sex, unsigned short homepatch);
  
  /**Completely creates an individual with inheritance and mutations on all traits.
   * Calls makeNewIndividual() to get the new offspring. 
   * @param mother ptr to the mother
   * @param father ptr to the father
   * @param sex gender of the individual
   * @param homepatch ID of the Patch where this individual is born, usually the current position in the Patch array
   **/
  Individual*             makeOffsprg                 (Individual* mother, Individual* father, sex_t sex, unsigned short homepatch);
  
  /**Individual prototype accessor.*/
  Individual*       getIndividualProtoype       ( )   {return &_protoIndividual;}
  
  /**Accessor to a TraitPrototype.
    *@param type the trait name
    **/
  TraitPrototype* getTraitPrototype (trait_t type);
  
  /**Accessor to the list of TraitPrototype's.*/
  std::map< trait_t,TraitPrototype* >& getTraitPrototypes ( )  {return _protoTraits;}
  
  /**Gives the index of trait with \a type.
    @param type the type of the trait (i.e. its "name")*/
  int getTraitIndex (trait_t type);
  
};

#endif

