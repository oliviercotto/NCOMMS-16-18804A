/**  $Id: individual.h,v 1.8.2.1 2014-04-29 18:29:24 fred Exp $
*
*  @file individual.h
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
*  created on @date Thu Aug 5 2004
*  @author fred
*/


#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H
#include <map>
#include <deque>
#include <limits>
#include "types.h"
#include "ttrait.h"
#include "ttrait_with_map.h"
#include "binarystoragebuffer.h"

/**This class contains traits along with other individual information (sex, pedigree, etc.\ ). 
 * The Individual class can be view as a trait container. It encapsulates the basic interface to manipulate
 * the traits an individual is carrying like initialization, inheritance, mutation and phenotype getter. It also
 * stores the basic individual info like its sex, pedigree (home Patch, mather and father id, etc.) and its
 * mating and fecundity values. It is not a StorableComponent but declares a similar interface to store and retrieve
 * all the previous info to/from a binary file.
 * All individuals in a simulation are instantiated by a call to the IndFactory class (through its derived class Metapop)
 * which contains the individual and traits prototypes.
 */
class Individual {
private:
  /**ID tag, unique for one simulation.*/
  unsigned long _id;
  /**Age.*/
  unsigned short _age;
  /**Sex tag.*/
  sex_t _sex;
  /**Parents ID tags.*/
//  unsigned long _motherID, _fatherID;   //MOD SLIM NEMO 7.11
  /**Parents pointers.*/
  Individual *_mother, *_father;
  /**Natal Patch tag.*/
  unsigned short _home;
  /**Pedigree class of the individual. 
    - 0: parents are from different demes.\
    - 1: parents are from the same deme but unrelated.\
    - 2: parents are half-sib\
    - 3: parents are full-sib\
    - 4: parents are the same individual (i.e. selfing).\
  */
//  unsigned char _pedigreeClass;  //MOD SLIM NEMO 7.11
  /**Assigned fecundity*/
  double _fecundity;
//  /**Mating counter. The relationship*/
//  unsigned short _matings[5];  //MOD SLIM NEMO 7.11
//  /**Number of surviving offspring from the different mating categories (see matings).*/
//  unsigned short _realizedFecundity[5];  //MOD SLIM NEMO 7.11
  /**Number of traits in the table*/
  unsigned int _trait_nb;
  
public:
  
  /**The ID counter, reset at the beginning of each simulation.*/
  static unsigned long currentID;
  
  typedef int IDX;
  
  /**The traits table.*/
  std::deque<TTrait*> Traits;
  
  Individual ();
  ~Individual () {clearTraits();}
  /**Inits parameters and traits.
   Called by IndFactory::makeNewIndividual() to allocate the traits' sequences memory.
   **/
  Individual*     init                    ();
  /**Resets parameters and traits values.
   * Does not de-allocate the traits' sequences memory (i.e. is not calling TTrait::reset() on each trait).
   * @callergraph
   **/
  void            reset                   ();
  
  ///@name Setters
  ///@{
  void            setID                  (unsigned long value)  {_id = value;}
  void            setAge                 (unsigned short value) {_age = value;}
  void            Aging                  ()                     {_age++;}
//  void            setFatherID            (unsigned long value)  {_fatherID = value;}//MOD SLIM NEMO 7.11
//  void            setMotherID            (unsigned long value)  {_motherID = value;}//MOD SLIM NEMO 7.11
  void            setFather              (Individual* f)        {_father = f;}
  void            setMother              (Individual* m)        {_mother = m;}
  void            setHome                (unsigned short value) {_home = value;}
  void            setSex                 (sex_t sex)            {_sex = sex;}
  void            setCurrentID           (unsigned long value)  {currentID = value;}
//  void            setIsSelfed            (bool s)               {_pedigreeClass = (s ? 4 : 0);}
//  void            setPedigreeClass       (Individual* mother, Individual* father) {_pedigreeClass = getPedigreeClass(mother,father);}
//  void            setPedigreeClass       (unsigned char ped)    {_pedigreeClass = ped;}

  ///@}
  
  ///@name Getters
  ///@{
  unsigned long   getID                  ()                      {return _id;}
  unsigned short  getAge                 ()                      {return _age;}
//  unsigned long   getFatherID            ()                      {return _fatherID;}//MOD SLIM NEMO 7.11
//  unsigned long   getMotherID            ()                      {return _motherID;}//MOD SLIM NEMO 7.11
  Individual*     getFather              ()                      {return _father;}
  Individual*     getMother              ()                      {return _mother;}
  unsigned short  getHome                ()                      {return _home;}
  sex_t           getSex                 ()                      {return _sex;}
  bool            isFemale               ()                      {return (_sex == FEM);}
//  bool            getIsSelfed            ()                      {return (_pedigreeClass == 4);}
  unsigned long   getcurrentID           ()                      {return currentID;}
  double          getFecundity           ()                      {return _fecundity;}
  /**Gives the number of matings that individual had with mates from a given pedigree class.
    @param cat the mating category (= the pedigree class of the offspring):
    - 0 = mates are from different patches
    - 1 = mates are from the same patch
    - 2 = mates are half sib
    - 3 = mates are full sib
    - 4 = selfed mating
  */
//  unsigned short  getMatings (unsigned int cat) //MOD SLIM NEMO 7.11
//  {
//    return _matings[cat];
//  }
//  /**Gives the number of times an individual mated with an individual from the same patch.*/
//  unsigned short  getLocalMatings        ( ) //MOD SLIM NEMO 7.11
//  {
//    return _matings[1]+_matings[2]+_matings[3]+_matings[4];
//  }
//  /**Gives the total number of matings of an individual.*/
//  unsigned int    getTotMatings          ( )  //MOD SLIM NEMO 7.11
//  {
//    return _matings[0] + getLocalMatings();
//  }
  /**Gives the number of surviving offspring for a given pedigree class of mating.
    @param cat the mating category:
    - 0 = mates are from different patches
    - 1 = mates are from the same patch
    - 2 = mates are half sib
    - 3 = mates are full sib
    - 4 = selfed mating
  */
//  unsigned short  getRealizedFecundity   (unsigned int cat)  //MOD SLIM NEMO 7.11
//  {
//    return _realizedFecundity[cat];
//  }
//
//  /**Gives the total number of surviving offspring when mating occures with mates of the same patch*/
//  unsigned int    getLocalRealizedFecundity( )  //MOD SLIM NEMO 7.11
//  {
//    return _realizedFecundity[1]+_realizedFecundity[2]+_realizedFecundity[3]+_realizedFecundity[4];
//  }
//
//  /**Gives the total number of surviving offspring for all categories of mating.*/  //MOD SLIM NEMO 7.11
//  unsigned int    getTotRealizedFecundity( )        {return _realizedFecundity[0] + getLocalRealizedFecundity( );}
//
//  /**Returns the pedigree class of the individual, as set during offspring creation.*/
//  unsigned int    getPedigreeClass       ( )  //MOD SLIM NEMO 7.11
//  {
//    return _pedigreeClass;
//  }
  /**Returns the pedigree class of two individuals.
    @param mother the first individual (e.g.the mother of this individual)
    @param father the second individual (e.g. the father of this individual)
    @return
             - 0 = mother and father are from different patches
             - 1 = mother and father are from the same patch
             - 2 = mother and father are half sib
             - 3 = mother and father are full sib
             - 4 = mother and father are the same individual (selfed offspring)
    */
  unsigned int    getPedigreeClass       (Individual* mother, Individual* father);
  ///@}
  
  ///@name implementation
  ///@{
  void store_data    ( BinaryStorageBuffer* saver  );
  
  void retrieve_data ( BinaryStorageBuffer* reader );
  ///@}
  
  ///@name Matings and Fecundity
  ///@{
  /**Sets the fecundity to the value given and returns it.
    @param value the fecundity*/
  double          setFecundity            (double value)         {_fecundity = value; return value;}

  /**Resets the mating and fecundity counters.*/
//  void            reset_counters         ( ) //MOD SLIM NEMO 7.11
//  {
//    for(unsigned int i = 0; i < 5; i++) {
//      _matings[i] = 0; _realizedFecundity[i] = 0;
//    } }

  /**Increments the mating counter according to the pedigree class of the offspring.
   * @param category the pedigree class of the offspring
   */ 
//  void   addMating  (unsigned int category) //MOD SLIM NEMO 7.11
//  { _matings[category]++; }

  /**Increments the mating and realized fecundity counters according to the pedigree class of the offspring.
   * @param category the pedigree class of the offspring
   */
//  void   DidHaveABaby  (unsigned int category) //MOD SLIM NEMO 7.11
//  { _matings[category]++; _realizedFecundity[category]++; }

  /**Returns the proportion of succesfull local matings (i.e. those that survived selection).*/
//  double getFecWithHomePatchMate  () //MOD SLIM NEMO 7.11
//  {
//    unsigned short mate = getLocalMatings();
//    return (mate != 0 ? (double) getLocalRealizedFecundity()/mate : 0.0);
//  }
  /**Returns the proportion of successfull remote matings.
    *@return _realizedFecundity[0]/_matings[0]
   */
//  double getFecWithOtherPatchMate () //MOD SLIM NEMO 7.11
//  {
//    return (_matings[0] != 0 ? (double) _realizedFecundity[0]/_matings[0] : 0.0);
//  }
  ///@}

  ///@name Trait interface
  ///@{
  
  /**Accessor to the size of the traits table.*/
  unsigned int getTraitNumber() { return _trait_nb;}
  
  /**Accessot to the traits table itself.
    * @return the traits table (a deque)*/
  std::deque<TTrait *>& getTraits() {return Traits;}
  
  /**Sets the phenotype/value of a trait to a particular value.
    * @param T the trait's index in the traits table
    * @param value the value passed to the trait (using the TTrait interface)*/
  void*  setTrait (IDX T, void* value) 
  { return getTrait(T)->set_trait(value); }
  
  /**Calls the value setting procedure of a particular trait.
    * @param T the trait's index in the traits table*/
  void   setTraitValue (IDX T)
  {  getTrait(T)->set_value(); }
  
  /**Calls the value setting procedure of all traits present in an individual.*/
  void   setTraitValue ()    
  { for(unsigned int i = 0; i < _trait_nb; i++) Traits[i]->set_value(); }
  
  /**Accessor to the value (phenotype) of a particular trait.
    * @param T the trait's index in the traits table
    * @return the trait's value, using the TTrait interface*/
  void*  getTraitValue (IDX T)
  { return getTrait(T)->getValue(); }
  
  /**Trait accessor.
    * @param T the trait's index in the traits table
    * @return the pointer to the trait*/
  TTrait*  getTrait (IDX T) 
  { if( T == -1 || !(T < (int)_trait_nb) )
    fatal("Individual::Trying to access a trait not present in the traits table (at %i, size %i)\n",T,_trait_nb);
    return Traits[T]; 
  }
  
  /**Adds a trait to the table.
   * @param theTrait pointer to the trait to add
     @param pos the position where the trait should be added in the traits table. 
                Used to check that the index of the trait in the table is correctly set.
     @see IndFactory::makePrototype*/
  void  addTrait (TTrait* theTrait, IDX pos)
  { if((int)_trait_nb != pos)
    fatal("Individual::adding a trait to the wrong position (at %i, size %i)!\n",pos,_trait_nb);
    Traits.push_back(theTrait); _trait_nb++;
  }
  
  /**Removes a trait from the table.
    * @param T the trait's index in the traits table*/
  void  removeTrait (IDX T)
  { delete Traits[T]; Traits.erase(Traits.begin() + T); _trait_nb--;}
  
  /**Clears the traits container.**/
  void  clearTraits ( )
  { if(_trait_nb != 0) {for(unsigned int i = 0; i < Traits.size(); ++i) delete Traits[i];
    Traits.clear(); _trait_nb = 0;}
  }
  
  /**Calls the inheritance procedure of a particular trait.
    * @param T the trait's index in the traits table
    * @param mother the mother
    * @param father the father */
  void  inheritTrait (IDX T, Individual* mother, Individual* father) 
  {  
    recombine(_id);
    getTrait(T)->inherit(mother->getTrait(T), father->getTrait(T)); 
  }
  
  /**Calls the mutation procedure of a particular trait.
    * @param T trait's index in the traits table*/
  void  mutateTrait (IDX T)
  {  getTrait(T)->mutate(); }
  
  /**Sets a particular trait's genotype and phenotype values from the two parents.
    Includes its inheritance, mutation and value setting.
    @param i the index of the trait to create in the traits table
    @param mother the mother
    @param father the father**/
  void createTrait (IDX i, Individual* mother, Individual* father)
  { if(!(mother && father)) 
      fatal("Individual::create::received null pointer!!!\n");
    TTrait* T = Traits[i];
    recombine(_id);
    T->inherit(mother->getTrait(i),father->getTrait(i));
    T->mutate();  
    T->set_value();
  }
  
  /**Creates an individual's genotypes and phenotypes with optional recombination or mutations on one trait only.
    @param i the index of the trait in the traits' table
    @param do_inherit if true, does inheritance and recombination of the trait
    @param do_mutate  if true, performs mutation of the trait**/
  Individual*  createTrait (IDX i, bool do_inherit, bool do_mutate)
  {
    if(do_inherit) inheritTrait(i, _mother, _father);
    if(do_mutate)  mutateTrait(i);
    setTraitValue(i);
    return this;
  } 
  
  /**Creates an individual's genotypes and phenotypes for first generation.**/
  Individual*  create_first_gen ()
  {
    for(unsigned int i = 0; i < _trait_nb; i++) {
      Traits[i]->init_sequence();
      Traits[i]->set_value();
    }
    //we have to set the parents ids otherwise the first offspring generation will be made of full sibs only.
    static unsigned long ID = std::numeric_limits< unsigned long >::max();
//    _motherID = ID--;//MOD SLIM NEMO 7.11
//    _fatherID = ID--;//MOD SLIM NEMO 7.11
    return this;
  }  
  
  /**Creates an individual's genotypes and phenotypes with recombination and mutations.**/
  Individual*  create ()
  { 
    inherit(_mother, _father);
    mutate();
    setTraitValue();


    return this;
  } 
  /**Creates an individual's genotypes and phenotypes with optional recombination or mutations.
    *@param do_inherit if true, do traits inheritance and recombination
    *@param do_mutate  if true, performs mutation of the traits**/
  Individual*  create (bool do_inherit, bool do_mutate)
  { 

    if(do_inherit) inherit(_mother, _father);

//#ifdef _DEBUG_
 // cout<<"inherit"<<endl;
//#endif

    if(do_mutate)  mutate();

//#ifdef _DEBUG_
 // cout<<"mutate"<<endl;
//#endif
    setTraitValue();

//#ifdef _DEBUG_
  //cout<<"set trait value"<<endl;
//#endif

    return this;
  } 
  
  /**Creates an individual, inherit, mutate and set all its trait's values.
    * @param mother the mother
    * @param father the father**/
  Individual*  create (Individual* mother, Individual* father)
  { 
    if(!(mother && father))
      fatal("Individual::create::received null parents pointer!!!\n");
    TTrait *TT;
    recombine(_id);
    for(unsigned int i = 0; i < _trait_nb; i++) {
      TT = Traits[i];
      TT->inherit(mother->getTrait(i), father->getTrait(i));
      TT->mutate();
      TT->set_value();
    }
    return this;
  }
  
  /**Calls the inheritance procedure of all the traits present in the individual.    
    * @param mother the mother
    * @param father the father**/ 
  void  inherit (Individual* mother, Individual* father)
  { 
    recombine(_id);
    for(unsigned int i = 0; i < _trait_nb; i++) 
      Traits[i]->inherit(mother->getTrait(i), father->getTrait(i)); 
  }
  
  /**Calls the mutation procedure of all the traits present in the individual. */
  void mutate ()
  { for(unsigned int i = 0; i < _trait_nb; i++) Traits[i]->mutate(); }
  
  void recombine (unsigned long ID)
  {  TTProtoWithMap::recombine(ID); }
  ///@}
  
  /**Write some info to stdout.*/
  void  show_up();
  
  /**Cloning procedure, clones all the traits present in the individual.*/
  Individual*  clone ();
  
  ///@name Operators
  ///@{
  /**Assignment, make a deep copy of the parameter values and traits.*/
  Individual& operator=(const Individual& i);
  /**Only checks for traits equivalence. */
  bool operator==(const Individual& i);
  bool operator!=(const Individual& i);
  ///@}
};

#endif //INDIVIDUAL_H

