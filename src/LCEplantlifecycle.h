/*
 * LCEplantlifecycle.h
 *
 *  Created on: Oct 12, 2016
 *      Author: fred
 */

#ifndef LCEPLANTLIFECYCLE_H_
#define LCEPLANTLIFECYCLE_H_

#include <cmath>
#include "LCEbreed.h"
#include "LCEdisperse.h"
#include "LCEselection.h"
#include "LCEmisc.h"         // ?
#include "lifecycleevent.h"

//CLASS LCE_Breed_Disperse_Clone_Regulate
//
/** Performs breeding, migration, cloning and regulation in one.
    The aim is to reduce the memory demand by producing only that amount of offspring that is needed after dispersal and regulation.
    Inherits parameters from LCE_Breed_base, LCE_Disperse_base, LCE_clonage.
*/

class LCE_Breed_Disperse_Clone_Regulate : public virtual LCE_Breed_base, public virtual LCE_Disperse_base, public virtual LCE_clonage//, public virtual LCE_Regulation
{ // it inherits from LCE_Breed, LCE_Disperse_base, LCE_clonage

private: // cannot be accessed from inherited classes or from outside

  double _competition_regulation;                // private in LCE_Regulation -> call header file ??
  vector<vector<double> > _reducedDispMatBack[2];

  // function output vectors and matrices
  vector<unsigned int> _postBreedOffspring; // offspring per patch after breeding       - vector[patch]                - get_postBreedOffspring()
  vector<unsigned int> _postDispOffspring ; // offspring per patch after dispersal      - vector[patch]                - get_get_postDispOffspring()
  vector< vector<unsigned int> > _newSeedlingsFromClones; //ordered: [patch][reverse adults age class without seedlings]
  vector<unsigned int> _newSeedlingsFromNewSeeds;
  vector< vector< Individual* > > _newSeedlingsFromSeedBank;
  vector< vector< Individual* > > _agedSeedlings; //temporary container before moving to pre-adults

  vector< vector<unsigned int> > _reducedpostDispOffspring[2];
  // how offspring is dispersing              - matrix[SEX][donor][receiver] - get_reducedpostDispOffspring()
  vector< vector<unsigned int> > _reducedpostDispOffspringBack[2];
  // how offspring is dispersing              - matrix[SEX][receiver][donor] - get_reducedpostDispOffspringBack()
  vector< vector<unsigned int> > _reducedpostRegulateOffspringBack[2];
  // offspring after dispersal and regulation - matrix[SEX][receiver][donor] - get_reducedpostRegulateOffspring()
  vector< vector<unsigned int> > _postSurvivalAdults;
  //number of surviving adults in each patch and age class, set in do_agingAdults_withoutMoving

  void (LCE_Breed_Disperse_Clone_Regulate::* _life_cycle)(void); //function pointer

public:
  LCE_Breed_Disperse_Clone_Regulate();

  virtual ~LCE_Breed_Disperse_Clone_Regulate (){}

  //void do_breed_disperse_clone_regulate             ();
  //void do_breed_disperse_clone_regulate_propagule   ();

  //life cycles:
  void aging_germination_cloning_regulation();
  void cloning_regulation();

  bool get_reducedDispMatBack();
  void get_postBreedOffspring();
  void get_postDispOffspring();
  void get_postDispOffspringBack();

  //life cycle 1:
  void do_agingAdults_withoutMoving();
  void do_SeedSurvivalAndGermination();
  void do_cloning();
  void do_SeedlingRegulation();
  void breed_newSeedlingsFromSeeds();
  void breed_newSeedlingsFromClones();
  void moveAgedAdults();

  //life cycle 2:
  void get_postRegulateOffspring();
  void Breed();
  void pro_Clone();


  virtual bool setParameters (); // A virtual member is a member function that can be redefined in a derived class, while preserving its calling properties through references.
  virtual void execute ();
  virtual LifeCycleEvent* clone () {return new LCE_Breed_Disperse_Clone_Regulate();}
  virtual void loadFileServices ( FileServices* loader ) {}
  virtual void loadStatServices ( StatServices* loader ) {}
  virtual age_t removeAgeClass ( ) {return NONE;}
  virtual age_t addAgeClass ( ) {return OFFSPRG;}
  virtual age_t requiredAgeClass () {return ADULTS;}


};


#endif /* LCEPLANTLIFECYCLE_H_ */
