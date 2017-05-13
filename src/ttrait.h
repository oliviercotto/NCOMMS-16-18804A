/**  $Id: ttrait.h,v 1.7.2.4 2014-05-01 14:48:08 fred Exp $
*
*  @file ttrait.h
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
*  created on @date 05.08.2004
*  @author fred
*/
 

#ifndef TTRAIT_H
#define TTRAIT_H

#include "types.h"
#include "simcomponent.h"

/**Interface for all trait types, declares all basic trait operations.
 * This pure abstract class declares the traits intreface. The precise genetic architecture 
 * of the trait is not defined and is up to the designer of the trait. It uses void pointers to
 * allow users to define their own structure. Existing traits use various types for their genes
 * sequence, like \c char, \c double or \c bitset. It is up to the trait's user to know what kind
 * of structure they expect.
 * The trait objects are contained in the Individual class. Their parameters are set by their 
 * TraitPrototype.
 **/

class TTrait: public StorableComponent {

public:
  /**Called to allocate the trait's genotypic sequences. Called each time a new Individual is created (\c Individual::init())**/
  virtual   void            init () = 0;
  
  /**Called at the start of each replicate, sets the initial genotypes. Called by \c Individual::create(). **/
  virtual   void            init_sequence () = 0;
  
  /**Called at the end of each simulation/replicate, deallocates sequence memory. **/
  virtual   void            reset () = 0;
  /** Inheritance procedure, creates a new trait from mother's and father's traits
   * @param mother the mother's trait
   * @param father the father's trait
   **/
  virtual   void            inherit (TTrait* mother, TTrait* father) = 0;
  /** Mutation procedure, perform mutations on the genes sequence. **/
  virtual   void            mutate () = 0;
  /** Called to set the phenotypic to a particular value or to give context-dependant value(s) to the trait.
    * @param value the value passed to the trait
    * @return the argument passed
    **/
  virtual   void*           set_trait (void* value) = 0;
  /** Called to set the sequence pointer to an existing trait
    * @param seq the existing sequence pointer
    **/
  virtual   void            set_sequence (void** seq) = 0;
  /**Tells the trait to set its phenotype from genotype, should be used instead of getValue().**/
  virtual   void            set_value () = 0;
  /** Genotype to phenotype mapper.
	* @return the phenotype computed from the genotype
    **/
  virtual   void*           getValue () const = 0;
  /** type accessor.
	* @return the trait's type
    **/
  virtual   trait_t         get_type () const = 0;
  /** sequence accessor.
	* @return the sequence pointer
    **/
  virtual   void**          get_sequence () const = 0;
  /** Called to read one allele value at a particular locus.
   * @return the allelic value at position 'all' at locus 'loc'
   * @param loc locus position in the sequence
   * @param all which allele we want to read the value from
   **/
  virtual   void*           get_allele   (int loc, int all) const = 0;
  /** Writes some info to stdout. **/
  virtual   void            show_up  () = 0;
  /** Returns a copy of itself. 
    \b Note: call the copy constructor of the trait which should only copy the parameters values
     not the complete state of the trait (i.e. shallow copy). The copy of the sequence data is
     made through the assignement operator!**/
  virtual   TTrait*         clone () = 0;
  ///@name Operators
  ///@{
  /**Copies the complete state of the trait from right to left side of the operator, sequence
     data included.*/
  virtual   TTrait& operator= (const TTrait&) = 0;
  /**Checks for parameters equivalence, not genetic equivalence.*/
  virtual   bool operator== (const TTrait&) = 0;
  virtual   bool operator!= (const TTrait&) = 0;
  ///@}
  virtual   ~TTrait ( ) { }
};

/**TTrait setter.
 * Encapsulates the methods to set the traits parameters and to generate traits. Also stores
 * the posistion of the trait in the individuals trait table.
 * This class manages the file and stat handlers through its inheritance of the SimComponent
 * interface.
*/
class TraitPrototype : public StorableComponent, public SimComponent {

protected:
  /**The trait index in the Individual traits table.*/
  int _index;
  
public:
  /**Called at the end of a simulation to reset the Traits' prototypes (e.g. unregister genetic maps).
     Called by \c IndFactory::clearPrototype. */
  virtual   void            reset () = 0;
  /**Creates the trait of which it is the prototype, called by \c IndFactory::makePrototype().*/
  virtual   TTrait*         hatch () = 0;
  /**Returns a copy of itself. \b Note: calls the copy constructor and only copy the parameters state.*/
  virtual   TraitPrototype* clone () = 0;
  /**Type accessor.
	* @return the trait's type
    **/
  virtual   trait_t         get_type () const = 0;
  /**Sets the traits index. Called by \c IndFactory::makePrototype().*/
  virtual void set_index(int idx) {_index = idx;}
  /**Index getter.*/
  virtual int  get_index() {return _index;}

};
#endif //TTRAIT_H

