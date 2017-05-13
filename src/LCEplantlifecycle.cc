/*
 * LCEplantlifecycle.cc
 *
 *  Created on: Oct 12, 2016
 *      Author: fred
 */

#include "LCEplantlifecycle.h"

/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

//                             ******** LCE_Breed_Disperse_Clone_Regulate ********/

// ----------------------------------------------------------------------------------------
// LCE_Breed_Disperse_Clone_Regulate::LCE_Breed_Disperse_Clone_Regulate
// ----------------------------------------------------------------------------------------
LCE_Breed_Disperse_Clone_Regulate::LCE_Breed_Disperse_Clone_Regulate() :
LifeCycleEvent("breed_disperse_clone_regulate",""), _competition_regulation(0), _life_cycle(0)
{  /**Interface to add a parameter and its updater to the set.
 * @param Name The param string name as read in the init file
 * @param Type The type of the argument (DBL=double, INT=integer, BOOL=boolean, STR=string, etc., see type.h)
 * @param isRequired True if the parameter is mandatory
 * @param isBounded True if the parameter takes bounded values as argument
 * @param low_bnd The lower value the argument can take
 * @param up_bnd The upper value the argument can take
 * @param updater a pointer to a ParamUpdaterBase object used to reset its value during a simulation */

  ParamUpdater<LCE_Breed_Disperse_Clone_Regulate> * updater = new ParamUpdater<LCE_Breed_Disperse_Clone_Regulate> (&LCE_Breed_Disperse_Clone_Regulate::setParameters);

  LCE_Disperse_base::addParameters("breed_disperse_clone_regulate", updater);

  add_parameter("competition_regulation",DBL,false,false,0,0);

  add_parameter("breed_disperse_clone_regulate_life_cycle",INT, true, true, 1, 2);

}
// ----------------------------------------------------------------------------------------

// ----------------------------------------------------------------------------------------
// LCE_Breed_Disperse_Clone_Regulate::setParameters
// ----------------------------------------------------------------------------------------
bool LCE_Breed_Disperse_Clone_Regulate::setParameters ()
{
  // if(!LCE_Regulation::setParameters ()) return false; // causes problems -> _popPtr is ambigous
//cout<<"LCE_Breed_Disperse_Clone_Regulate::setParameters : \n";

  if(!LCE_Breed_base::setParameters()) return false;
//cout<<"LCE_Breed_base::setParameters::OK\n";

  if(!LCE_Breed::setParameters()) return false;
  //cout<<"LCE_Breed::setParameters::OK\n";

  if(!LCE_Disperse_base::setBaseParameters("breed_disperse_clone_regulate")) return false; // determine prefix for disperse
  //cout<<"LCE_Disperse_base::setBaseParameters::OK\n";

  if(!LCE_clonage::setParameters ()) return false;
  //cout<<"LCE_clonage::setParameters::OK\n";

  if(!get_reducedDispMatBack()) return false;
  //cout<<"get_reducedDispMatBack::OK\n";

  if(get_parameter("competition_regulation")->isSet())
      _competition_regulation = get_parameter_value("competition_regulation");


  if(get_parameter_value("breed_disperse_clone_regulate_life_cycle") == 1)

      _life_cycle = &LCE_Breed_Disperse_Clone_Regulate::aging_germination_cloning_regulation;

  else if(get_parameter_value("breed_disperse_clone_regulate_life_cycle") == 2)

      _life_cycle = &LCE_Breed_Disperse_Clone_Regulate::cloning_regulation;
  else
     return error("\"breed_disperse_clone_regulate_life_cycle\" not properly set, should be 1 or 2\n");

//cout<<"OK\n";
  return true;
}

// ----------------------------------------------------------------------------------------
// create backwards migration matrix:                                                          (output: _reducedDispMatBack[matrix])
// ----------------------------------------------------------------------------------------    (demand: _reducedDispMat[matrix])
bool LCE_Breed_Disperse_Clone_Regulate::get_reducedDispMatBack()
{
  _reducedDispMatBack[0].clear();

  //cout<<"+++ get_reducedDispMatBack: size of reduced disp mat[0]: "<<_reducedDispMat[0].size()
    //  <<"; reduced disp mat[1]: "<<_reducedDispMat[1].size()<<endl;

  for (unsigned int i = 0; i < _reducedDispMat[0].size(); ++i)
    _reducedDispMatBack[0].push_back(vector<double>());

  for (unsigned int i = 0; i < _reducedDispMat[0].size(); ++i) {         // loop through every donor patch number (i)

    for(unsigned int j = 0; j < _reducedDispMat[0][i].size(); ++j){       // loop through each receiver patch position (j)

      _reducedDispMatBack[0][ _reducedDispMat[0][i][j] -1 ].push_back(i + 1);

    }
  }

  //cout<<"+++ get_reducedDispMatBack: OK\n";
  return true;
}

// ----------------------------------------------------------------------------------------
// LCE_Breed_Disperse_Clone_Regulate::execute
// ----------------------------------------------------------------------------------------
void LCE_Breed_Disperse_Clone_Regulate::execute()
{
#ifdef _DEBUG_
  message("LCE_Breed_Disperse_Clone_Regulate::execute (Patch nb: %i offsprg nb: %i adlt nb: %i "
          ,_popPtr->getPatchNbr(), _popPtr->size( OFFSPRG ), _popPtr->size( ADULTS ));
#endif

  if(_npatch != _popPtr->getPatchNbr()) {
    //cout << "_npatch updated" << endl;
    _npatch = _popPtr->getPatchNbr();
    if(!updateDispMatrix()) fatal("bailing out\n");
  }

//cout<<endl;

  // LIFE CYCLE:
  (this->*_life_cycle)();

}
// ----------------------------------------------------------------------------------------
// cloning_regulation --> life cycle option 2 (previous, wrong, life cycle, optimized)
// ----------------------------------------------------------------------------------------
void LCE_Breed_Disperse_Clone_Regulate::cloning_regulation()
{
    // 1. REPRODUCTION - count number of new seeds, don't instantiate any new individuals
    get_postBreedOffspring();

    // 2. SEED DISPERSAL - count number of emigrants from each patch to each connected patch
    get_postDispOffspring();

    // 2.1 reverse the postDispOffspring count table to get the number of immigrants per connected patch
    get_postDispOffspringBack();

    // 3. Cloning
    pro_Clone();

    // 4. offspring regulation
    get_postRegulateOffspring();

    // 5. actual offspring production
    Breed();
}
// ----------------------------------------------------------------------------------------
// aging_germination_cloning_regulation --> life cycle option 1
// ----------------------------------------------------------------------------------------
void LCE_Breed_Disperse_Clone_Regulate::aging_germination_cloning_regulation()
{
  // 1. REPRODUCTION - count number of new seeds, don't instantiate any new individuals
    get_postBreedOffspring();

    // 2. SEED DISPERSAL - count number of emigrants from each patch to each connected patch
    get_postDispOffspring();

    // 2.1 reverse the postDispOffspring count table to get the number of immigrants per connected patch
    get_postDispOffspringBack();


    // --- END OF CURRENT SEASON/YEAR --


    // 3. ADULT SURVIVAL AND AGING -- COUNT NUM ADULTS SURVIVING INTO THE NEXT SEASON

    // count surviving adults in each stage, but keep all individuals in place, they are
    // needed to create individuals from the new seeds counted in step 1.
    // the only ones to move are the seedlings, put in a temp container, seedling container is emptied

    do_agingAdults_withoutMoving();

    // --- NEW SEASON --- START WITH GERMINATION FROM SURVING SEEDS AND CLONING FROM SURVING ADULTS


    // 4. OFFSPRING SURVIVAL AND GERMINATION -- add/keep seeds to the seed bank and do germination from them

    // new seeds surviving to the seed bank are created from parents and added to offspring container
    // the new seeds are created from parents of the previous season (i.e. pre-survival)

    do_SeedSurvivalAndGermination();


    // 5. CLONING FROM SURVIVING ADULTS -- cloning only from pre-adults and adults, not from seedlings

    // don't create them now, will be created at the end, once the new adult cohort is known
    // we count the number of clones produced from each each class

    do_cloning();


    // 6. SEEDLING COMPETITION REGULATION --

    // count surviving seedlings from: new seeds, seeds in seed bank, and clones

    do_SeedlingRegulation();


    // 7. BREED/CREATE INDIVIDUALS THAT SURVIVED FROM SEEDS CREATED IN PREVIOUS SEASON (NEW SEEDS)

    // new individuals and germinated seeds from the seed bank are added to the seedling container
    // for this event, we need the pre-survival adults (from the previous season) still in their container

    breed_newSeedlingsFromSeeds();


    // 8. PROCEED WITH ADULTS AGING

    // create new adult 'cohort' for this season
    // we now randomly select the adults that survived from previous season and move them to their container

    moveAgedAdults();


    // 9. CREATE CLONES

    //effectively create seedlings from clones of the post-survival adults that are now in place

    breed_newSeedlingsFromClones();
}

// ----------------------------------------------------------------------------------------
// breed*:     calculate how much offspring (age class 0 = seeds) would be added by breeding   (output: _postBreedOffspring[patch])
// ----------------------------------------------------------------------------------------    (demand: _LeslieMatrix, _popPtr, _coeff_competition)
void LCE_Breed_Disperse_Clone_Regulate::get_postBreedOffspring()                         //    (comments: density dependent on adult + offspring number)
{
  //declaration
  Patch* patch;
  age_idx current_age;
  unsigned int nbBaby;

//  _postBreedOffspring.clear();
//
//  for (unsigned int i = 0; i < _reducedDispMat[0].size(); ++i)
//    _postBreedOffspring.push_back(0);

  // initialize counter
  _postBreedOffspring.assign(_popPtr->getPatchNbr(), 0);

  unsigned int nb_class = _popPtr->getLeslieMatrix()->getNbCols();        // get number of age classes

  #ifdef _DEBUG_                                                          // debug message
  message("LCE_Breed_Disperse_Clone_Regulate::get_postBreedOffspring (Patch nb: %i offsprg nb: %i adlt nb: %i)\n"
      ,_popPtr->getPatchNbr(),_popPtr->size( OFFSPRG ),_popPtr->size( ADULTS ));
  #endif

  if(_popPtr->size(OFFSPRG) != 0 and _LeslieMatrix->get(0, 0)==0) {      // check if offspring container is empty

    warning("offspring containers not empty at time of breeding, flushing.\n");

    _popPtr->flush(OFFSx);

  }

  //cout<<"+++ reproduction\n";

  for(unsigned int i = 0; i < _popPtr->getPatchNbr(); ++i) {   // loop through each patch

    patch = _popPtr->getPatch(i);                            // get patch pointer

    if( !checkMatingCondition(patch) ) continue;             // check carrying capacity

    for(int j = (int)nb_class-1; j > 0; j--){               // loop through age classes except offspring

      current_age = static_cast<age_idx> (j);

      if(_LeslieMatrix->get(0, j) > 0){

        for(unsigned int size = patch->size(FEM, current_age), indexOfMother = 0;
            indexOfMother < size;                                          // loop through each female adult
            indexOfMother++)
        {

          nbBaby = (unsigned int)round(RAND::Poisson(_LeslieMatrix->get(0,j))) ;

          _postBreedOffspring[i] += nbBaby ;

        }
      }
    }

    //cout<<"  patch "<<i<<" > actual state: |0: "<<patch->size(FEM, OFFSx)<<" ("<<_postBreedOffspring[i]<<")"
      //  <<" |1: "<<patch->size(FEM, ADLTx)
      //  <<" |2: "<<patch->size(FEM, static_cast<age_idx> (2))
      //  <<" |3: "<<patch->size(FEM, static_cast<age_idx> (3))<<endl;
  }
}

// ----------------------------------------------------------------------------------------
// dispersal*: create vector and matrix of how offspring (seeds) would be distributed          (output: _postDispOffspring[patch], _reducedpostDispOffspring[matrix])
// ----------------------------------------------------------------------------------------    (demand: _reducedDispMatProba)
void LCE_Breed_Disperse_Clone_Regulate::get_postDispOffspring()
{
  _reducedpostDispOffspring[0].clear();                                   // _reducedDispMatProba[SEX][LocalPatch][AimedPatch]

  for (unsigned int i = 0; i < _reducedDispMat[0].size(); ++i){

    _reducedpostDispOffspring[0].push_back(vector<unsigned int>());

    for (unsigned int j = 0; j < _reducedDispMat[0][i].size(); ++j){

      _reducedpostDispOffspring[0][i].push_back(0);

    }
  }

//  _postDispOffspring.clear();
//
//  for (unsigned int i = 0; i < _reducedDispMat[0].size(); ++i)
//    _postDispOffspring.push_back(0);      

  // initialize counter
  _postDispOffspring.assign(_popPtr->getPatchNbr(), 0);

  Patch* patch;                                                           // clear every generation + initiate with zeros?
  double sum = 0, random;
  unsigned int AimedPatch = 0;

  for (unsigned int i = 0; i < _reducedDispMatProba[0].size(); i++) {     // loop through every donor patch (i)

    patch = _popPtr->getPatch(i);

    for(unsigned int k = 0; k < _postBreedOffspring[i]; k++){           // loop through each offspring (k) of this donor patch i

      sum = 0;
      random = RAND::Uniform();                       // pick a random number between 0 and 1
      AimedPatch = 0;

      if(random > 0.999999) random = 0.999999;                        // this to avoid overflows when random == 1

      sum =  _reducedDispMatProba[0][i][AimedPatch];                  // get dispersal probability

      while (random > sum) {                                          // while random number is larger than the summed probability continue

        AimedPatch++;

        sum +=  _reducedDispMatProba[0][i][AimedPatch];

      }

      _reducedpostDispOffspring[0][i][AimedPatch] += 1;               // add offspring to the aimed patch (= stochastic process)

      _postDispOffspring[ _reducedDispMat[0][i][AimedPatch]-1 ] += 1;
      // the -1 is because patch number are [1 - num patch], given in input, to be changed!
    }
  }
}

// ----------------------------------------------------------------------------------------
// create backwards migration matrix:                                                          (output: _reducedpostDispOffspringBack[matrix])
// ----------------------------------------------------------------------------------------    (demand: _reducedpostDispOffspring[matrix], _reducedDispMatBack[matrix])
void LCE_Breed_Disperse_Clone_Regulate::get_postDispOffspringBack()
{
  
  // reverse the count table of post-dispersal offspring to get
  // the number of immigrants from each connected patch to the focal patch
  
  unsigned int rpn;  // receiver patch number
  unsigned int dpos; // donor patch number

  //[FG] recreating the whole table would be necessary only if the number of patches changes
  if(_reducedpostDispOffspringBack[0].size() != _popPtr->getPatchNbr()) {

    _reducedpostDispOffspringBack[0].clear();
  
    for (unsigned int i = 0; i < _reducedDispMatBack[0].size(); ++i){
    //add a container for each patch
      _reducedpostDispOffspringBack[0].push_back(vector<unsigned int>());
    }
  }

  
  for (unsigned int i = 0; i < _reducedDispMatBack[0].size(); ++i)
    _reducedpostDispOffspringBack[0][i].clear();
  

  for (unsigned int i = 0; i < _reducedpostDispOffspring[0].size(); i++) {
    // loop through every donor patch (i)

    for(unsigned int j = 0; j < _reducedpostDispOffspring[0][i].size(); j++){
      // loop through each receiver patch position (j)

      rpn = _reducedDispMat[0][i][j] -1 ; // subtract one because patch number in input are given [1 - num patch]!!!
      // get receiver patch number (rpn)

      dpos = _reducedpostDispOffspringBack[0][rpn].size() ;
      // get position to add offspring

      if(_reducedDispMatBack[0][rpn][dpos] -1 == i ){
        // check if it is added to the correct position

        _reducedpostDispOffspringBack[0][rpn].push_back( _reducedpostDispOffspring[0][i][j] ) ;
        // clean _reducedpostDispOffspringBack before each "refilling"!

      } //@TODO [FG] what if not?? -> error message or a fallback is necessary
    }
  }
}

// ----------------------------------------------------------------------------------------
// clone:     create clones (age class 1 = seedlings) based on adult + offspring* number       (output:  adds seedlings)
// ----------------------------------------------------------------------------------------    (demand: _popPtr, _coeff_competition_cloning, _clonal_rate)
void LCE_Breed_Disperse_Clone_Regulate::do_cloning()                                      
//    (comments: density dependent on adult + offspring
{
  Patch* patch;
  
  unsigned int last = _popPtr->getNumAgeClasses() - 1; //will be 3

  age_idx current_age;
  
  //we exclude the seedling and offspring stages -> will be 2 stages contributing clones for our case
  unsigned int num_contributing_stages = _popPtr->getNumAgeClasses() - 2;
  unsigned int num_last_adlt;
  
  
  //cout<<"+++ do cloning\n";

// inits
  if(_newSeedlingsFromClones.size() != _popPtr->getPatchNbr()) {

    _newSeedlingsFromClones.clear();

    for(unsigned int i = 0; i < _popPtr->getPatchNbr(); ++i) {

      _newSeedlingsFromClones.push_back(vector<unsigned int>());
    }
  }

  for(unsigned int i = 0; i < _popPtr->getPatchNbr(); ++i)
    _newSeedlingsFromClones[i].assign( num_contributing_stages , 0);


  for(unsigned int i = 0; i < _popPtr->getPatchNbr(); i++) {               // loop through each patch

    patch = _popPtr->getPatch(i);

    
    ///////// LAST AGE CLASS -- ADULTS //////////////////////////////////////////////
    //cloning is happening from adults that survived through the new season
    
    //the number of adults in the last age class is the sum of the ones that survived there
    //and the ones that aged from the pre-adults:
    num_last_adlt = _postSurvivalAdults[i][last] + _postSurvivalAdults[i][last-1];

    for(unsigned int k = 0; k < num_last_adlt; ++k)           // loop adult each adult female
    {
      _newSeedlingsFromClones[i][last - num_contributing_stages] += (unsigned int)round(RAND::Poisson(_clonal_rate));                   
      // get number clones
    }

    ///////// WHICH AGE CLASS DOES CLONING??? -> EXCLUDE SEEDLINGS, BUT NOT PRE-ADULTS
    for( int j = (int)last -1; j > 1; j--) {

      //we have to consider the survivors of the previous age class in the new season:
      for(unsigned int k = 0; k < _postSurvivalAdults[i][j-1]; ++k)           // loop adult each adult female
      {
         _newSeedlingsFromClones[i][j - num_contributing_stages] += (unsigned int)round(RAND::Poisson(_clonal_rate));                       
        // get number clones
      }
    }

    //cout<<"  patch "<<i<<" > new clones: "
      //  <<_newSeedlingsFromClones[i][0]<<" (from "<<_postSurvivalAdults[i][last-2]
	  //  <<": "<<(double)_newSeedlingsFromClones[i][0]/_postSurvivalAdults[i][last-2]<<") + "
          //  <<_newSeedlingsFromClones[i][1]
    //  <<" (from "<<_postSurvivalAdults[i][last]+_postSurvivalAdults[i][last-1]
    //  <<": "<<(double)_newSeedlingsFromClones[i][1]/(_postSurvivalAdults[i][last]+_postSurvivalAdults[i][last-1])<<")"
    //  <<endl;
  }
}

// ----------------------------------------------------------------------------------------
// regulate*:  calculate effect of regulation (= only on offspring)                            (output: _reducedpostRegulateOffspring[matrix])
// ----------------------------------------------------------------------------------------    (demand: _popPtr, _reducedDispMatProba, _competition_regulation, _reducedDispMat)
void LCE_Breed_Disperse_Clone_Regulate::get_postRegulateOffspring()                      //    (comments: density dependent on adults)
{

  Patch* patch;                    // clear every generation ??
  double capacityRatio = 0;

  //@TODO [FG] probably not necessary to recreate the whole table every generation!
  _reducedpostRegulateOffspringBack[0].clear();
  for (unsigned int k = 0; k < _reducedDispMatBack[0].size(); ++k){

    _reducedpostRegulateOffspringBack[0].push_back(vector<unsigned int>());

  }

  unsigned int cntr=0;

  //cout << "+++ regulation\n";

  for (unsigned int i = 0; i < _reducedpostDispOffspringBack[0].size(); i++) {
    // loop through each receiving patch (i)

    patch = _popPtr->getPatch(i);

    unsigned int pop_size = patch->size(ADULTS);

    double num_off  = _postDispOffspring[i];

    //@TODO [FG] is the following really necessary?
    double sum = 0;
    for(unsigned int a = 0; a < _reducedpostDispOffspringBack[0][i].size(); ++a)
      sum += _reducedpostDispOffspringBack[0][i][a];

    if(sum==num_off){
      // check if back migration matrix has been filled up correctly

      //@TODO  [FG] what if not?? no fallback or error message!

      double compet = 1/(1 + pop_size * _competition_regulation);
      // get competition coefficient

      unsigned int num_off_surviving = RAND::Binomial2(compet,num_off);
      // get surviving offspring

//      capacityRatio = 0;

      capacityRatio = num_off_surviving / num_off;
      // calculate ratio pre/post-regulation

      //cout<<"  patch "<<i<<" > surv offspring: "<<num_off_surviving<<" (from: "<<sum<<"; compet: "<<compet<<")";
      cntr = 0;
      if(capacityRatio < 1){
        unsigned int num;
        for(unsigned int j=0; j < _reducedpostDispOffspringBack[0][i].size(); j++) {
          // loop through reduced backdispersal matrix
          num = RAND::Binomial2( capacityRatio, _reducedpostDispOffspringBack[0][i][j]);
          cntr += num;
          _reducedpostRegulateOffspringBack[0][i].push_back(num);
        }
      }

      //cout<<" > effective num surv off: "<<cntr<<endl;

    } else
       fatal("something went wrong in get_postRegulateOffspring()\n");

    //cout<<"  patch "<<i<<" > actual state: |0: "<<patch->size(FEM, OFFSx)<<" ("<<_postBreedOffspring[i]<<" -> "<<cntr<<")"
      //  <<" |1: "<<patch->size(FEM, ADLTx)
      //  <<" |2: "<<patch->size(FEM, static_cast<age_idx> (2))
      //  <<" |3: "<<patch->size(FEM, static_cast<age_idx> (3))<<endl;
  }
}
// ----------------------------------------------------------------------------------------
// clone:     create clones (age class 1 = seedlings) based on adult + offspring* number       (output:  adds seedlings)
// ----------------------------------------------------------------------------------------    (demand: _popPtr, _coeff_competition_cloning, _clonal_rate)
void LCE_Breed_Disperse_Clone_Regulate::pro_Clone()                                      //    (comments: density dependent on adult + offspring
{

  //declaration
  Individual* newind;
  Individual* mother;
  Patch* patch;
  unsigned int pop_size;
  double number_clones;
  unsigned int clones;
  double compet_clone;

  unsigned int cntr = 0;

  //cout<<"+++ do cloning (rate = "<<_clonal_rate<<", coeff compet = "<<_coeff_competition_cloning<<")\n";

  for(unsigned int i = 0; i < _popPtr->getPatchNbr(); i++) {               // loop through each patch

    patch = _popPtr->getPatch(i);

    pop_size = patch->size(ADULTS); // + _postDispOffspring[i];                // get population size within patch i

    compet_clone = 1/(1 + pop_size * _coeff_competition_cloning);          // get competition coefficient (Beverton-Holt function)

    number_clones= compet_clone * _clonal_rate;                            // get number of clones

    //cout<<"  patch "<<i<<" > mean clones: "<<number_clones<<" (compet : "<<compet_clone<<") > new clones: "<<flush;

    cntr = 0;

    for(unsigned int size = patch->size(FEM, ADULTS),
        indexOfMother = 0;indexOfMother < size; indexOfMother++)           // loop adult each adult female
    {
        mother = patch->get(FEM, ADULTS, indexOfMother);                   // get mother

        clones= round(RAND::Poisson(number_clones));                       // get offspring number

        cntr += clones;

        for (unsigned int k = 0 ; k < clones; k++){                                 // loop through each created clone and exucete cloning

          newind = _popPtr->getNewIndividual();

          //cloning:
          (*newind) = (*mother);
//          newind->reset_counters();
          newind->setFather(NULL);
//          newind->setFatherID(0);//MOD SLIM NEMO 7.11
          newind->setMother(mother);
//          newind->setMotherID(mother->getID());//MOD SLIM NEMO 7.11
//          newind->setIsSelfed(true);
          newind->setHome(i);
          newind->setAge(1);                                                 // cloning adds seedlings !!

          patch->add(FEM,ADLTx,newind->create(false,true));
        }
    }

    //cout<<cntr<<endl;
  }


}

// ----------------------------------------------------------------------------------------
// breed:      actually breed individuals and disperse them
// ----------------------------------------------------------------------------------------
//(demand: _popPtr, _reducedDispMat)
void LCE_Breed_Disperse_Clone_Regulate::Breed()
{
  //declarations
  Patch* patch;
  Individual* mother;
  Individual* father;
  Individual* NewOffsprg;
  age_idx current_age;

  unsigned int nb_class = _popPtr->getLeslieMatrix()->getNbCols();

  #ifdef _DEBUG_
  message("LCE_Breed::breed_effectiveOffspring (Patch nb: %i offsprg nb: %i adlt nb: %i)\n"
      ,_popPtr->getPatchNbr(),_popPtr->size( OFFSPRG ),_popPtr->size( ADULTS ));
  #endif

  if(_popPtr->size(OFFSPRG) != 0 and _LeslieMatrix->get(0, 0)==0) {                         // check if offspring container is empty

    warning("offspring containers not empty at time of breeding, flushing.\n");

    _popPtr->flush(OFFSx);

  }

//cout<<"+++ offspring production\n";

  for(unsigned int i = 0; i < _reducedpostRegulateOffspringBack[0].size(); i++) {                                    // loop trough each receiver patch (number i)

    patch = ( _popPtr->getPatch(i));

    if( !checkMatingCondition(patch) ) continue;

    for(int j = (int)nb_class-1; j > 0; j--){                                                // loop through each age classe (as it potentially contributes offspring)

      current_age = static_cast<age_idx> (j);

      if(_LeslieMatrix->get(0, current_age) > 0){

        for(unsigned int k = 0; k < _reducedpostRegulateOffspringBack[0][i].size() ; ++k){    // loop through each donor patch (position k)

          patch = ( _popPtr->getPatch(_reducedDispMatBack[0][i][k] - 1));

          for(unsigned int l = 0; l < _reducedpostRegulateOffspringBack[0][i][k] ; ++l){    // loop through each offspring to be created

            int indexOfMother = RAND::Uniform(patch->size( current_age )-1);              // randomly pick a female

            mother = patch->get(FEM, current_age, indexOfMother);

            father = this->getFatherPtr(patch, mother, indexOfMother);

            NewOffsprg = makeOffspring( do_breed(mother, father, i) );

            _popPtr->getPatch(i) -> add(NewOffsprg->getSex(), OFFSx, NewOffsprg);

          }
        }
      }
    }
  }
}

// ----------------------------------------------------------------------------------------
// do_agingAdults_withoutMoving
// ----------------------------------------------------------------------------------------
void LCE_Breed_Disperse_Clone_Regulate::do_agingAdults_withoutMoving ()
{

  Individual* current_ind;
  Patch* current_patch;

  unsigned int age, nb_class, ind;
  double surv_rate;

  TMatrix *survival = _popPtr->getLeslieMatrix();

  nb_class = survival->getNbCols(); //that should be the same as _popPtr->getNumAgeClasses()

  unsigned int last = nb_class - 1; //index of the last age class

  age_idx last_idx = static_cast<age_idx> (nb_class -1);//!!age_idx is a type. type conversion.

  age_idx current_age;

  // init containers -----------------------------------------------------------
  if(_postSurvivalAdults.size() != _popPtr->getPatchNbr()) {

    _postSurvivalAdults.clear();
    _agedSeedlings.clear();

    for(unsigned int i = 0; i < _popPtr->getPatchNbr(); ++i) {
      _postSurvivalAdults.push_back(vector<unsigned int>());
      _agedSeedlings.push_back(vector<Individual*>());
    }
  }

  for(unsigned int i = 0; i < _popPtr->getPatchNbr(); ++i) {
    
    _postSurvivalAdults[i].assign(_popPtr->getNumAgeClasses(), 0);
    //this also assigns 0 for offspring even if no offspring stage considered
    //=> avoids subtracting 1 from stage indices below

    _agedSeedlings[i].clear(); //empty the backup container before re-filling it
  }

  //cout<<"+++ adult survival - stay in place\n";

  //for each patch
  for(unsigned int i = 0; i < _popPtr->getPatchNbr(); i++)
  {
    current_patch = _popPtr->getPatch(i);

//////////////////////// LAST AGE CLASS ////////////////////////////////////////////////////

    if( (surv_rate = survival->get(last, last)) != 0) //get the last row and column of the leslie matrix

    {

      //individual female in that class may stay if they survive

      _postSurvivalAdults[i][last] += RAND::Binomial2(survival->get(last, last), current_patch->size(FEM, last_idx));

    }



////////////////////////AGING IN OTHER AGE-CLASSES - EXCEPT OFFSPRING///////////////////////

    for( int j = (int)last -1; j > 0; j--) {

      current_age = static_cast<age_idx> (j); //this is the index in the containers

      //we capitalize on the fact that pre-adults and seedlings stay so only one year
      
      // WE DON'T CHECK FOR THE AGE OF THE INDIVIDUALS!!! this works only our life cycle

      _postSurvivalAdults[i][j] += RAND::Binomial2(survival->get(j+1, j),
                                                     current_patch->size(FEM, current_age));

      // SEEDLINGS - special treatment
      if(j == 1) { 
        
        //move survivors to backup container to free the container that will receive new seedlings later
        //we need to keep the pre-adults in their container for when creating the new offspring
        //all surviving adults from the previous year are later moved to their destination stage in 'moveAgedAdults()'

        for(unsigned int s = 0; s < _postSurvivalAdults[i][j] &&
                               current_patch->size(FEM, current_age) != 0;
                               ++s)
        {

          ind = RAND::Uniform(current_patch->size(FEM, current_age));

          current_ind = current_patch->remove(FEM, current_age, ind);

          _agedSeedlings[i].push_back(current_ind);

          current_ind->Aging();//still need to update the age tag

        }
        //trash the losers
        current_patch->flush(FEM, current_age, _popPtr);

      }
    } //end_for j age classes

    //cout<<"  patch "<<i<<" > age structure: 0: "<<_postSurvivalAdults[i][0]
      //  <<"; 1: "<<_postSurvivalAdults[i][0]<<" ("<<_agedSeedlings[i].size()<<" to 2)"
      //  <<"; 2: "<<_postSurvivalAdults[i][1]<<" ("<<_postSurvivalAdults[i][2]<<" to 3)"
      //  <<"; 3: "<<_postSurvivalAdults[i][2]+_postSurvivalAdults[i][3]<<" ("<<_postSurvivalAdults[i][3]<<" stay in 3)"<<endl;

    //cout<<"  patch "<<i<<" > actual state: 0: "<<current_patch->size(FEM, OFFSx)
      //  <<"; 1: "<<current_patch->size(FEM, ADLTx)<<" (is 0?)"
      //  <<"; 2: "<<current_patch->size(FEM, static_cast<age_idx> (2))
      //  <<"; 3: "<<current_patch->size(FEM, static_cast<age_idx> (3))<<endl;

  } //end_for i patchNbr

}
// ----------------------------------------------------------------------------------------
// do_SeedSurvivalAndGermination
// ----------------------------------------------------------------------------------------
void LCE_Breed_Disperse_Clone_Regulate::do_SeedSurvivalAndGermination ()
{
  Individual* current_ind;
  Patch* current_patch;

  unsigned int age, nb_class;
  double surv_rate;

  TMatrix *survival = _popPtr->getLeslieMatrix();

  nb_class = survival->ncols();

  age_idx last_idx = static_cast<age_idx> (nb_class -1);

  // inits -----------------------------------------------------------------------------------
  if(_reducedpostRegulateOffspringBack[0].size() != _popPtr->getPatchNbr()) {
  
    //that table first counts number of seedlings produced from new seeds
    //from local and immigrant seeds (need to know the origin for breeding later)
    //will be updated during competition regulation with num surviving
    _reducedpostRegulateOffspringBack[0].clear();
    
    assert(_reducedDispMatBack[0].size() == _popPtr->getPatchNbr());
    
    for (unsigned int k = 0; k < _reducedDispMatBack[0].size(); ++k){
      
      _reducedpostRegulateOffspringBack[0].push_back(vector<unsigned int>());
    }
  }
  
  for (unsigned int k = 0; k < _reducedDispMatBack[0].size(); ++k)
    _reducedpostRegulateOffspringBack[0][k].assign(_reducedDispMatBack[0][k].size(), 0);
  
  _newSeedlingsFromNewSeeds.assign(_popPtr->getPatchNbr(), 0);
  
  if(_newSeedlingsFromSeedBank.size() != _popPtr->getPatchNbr()) {
    
    _newSeedlingsFromSeedBank.clear();
    
    for(unsigned int i = 0; i < _popPtr->getPatchNbr(); ++i)
      _newSeedlingsFromSeedBank.push_back(vector<Individual*>());
    
  } else
    for(unsigned int i = 0; i < _popPtr->getPatchNbr(); ++i)
      _newSeedlingsFromSeedBank[i].clear();
  
  
  //cout<<"+++ seed survival and germination\n";


  //for each patch -------------------------------------------------------------------------
  for(unsigned int i = 0; i < _popPtr->getPatchNbr(); i++)
  {
    current_patch = _popPtr->getPatch(i);

    //SURVIVAL IN OFFSPRING WITH A SEED BANK
    //the seed bank is simply the offspring container, with a constant survival rate in it

    //1. deal with the offspring already present in the seedbank ---------------------------
    // note: in absence of a seedbank, the patch should currently not contain any offspring

    //cout<<"  patch "<<i<<" > seeds in bank: "<<current_patch->size(FEM, OFFSx)<<flush;

    for(int k = (int)current_patch->size(FEM, OFFSx) - 1; k >= 0; --k) {

      current_ind = current_patch->get(FEM, OFFSx, k);

      //if absence of a survival rate, check for survival and move them
      //this is silly, if survival == 0 there is no seed bank
      if(survival->get(0,0) == 0){

        if(RAND::Uniform() > survival->get(1, 0)) {
          _popPtr->recycle(current_ind);
          current_patch->remove(FEM, OFFSx, k);
        }

        else {
          current_ind->Aging();
//          current_patch->move(FEM, OFFSx, ADLTx, k);  //move to seedlings
          
          _newSeedlingsFromSeedBank[i].push_back(current_ind);
          
          current_patch->remove(FEM, OFFSx, k);
        }

      } else {

        //note that in this case survival off-adlt is the maturation rate
        //here, it is the germination rate of individuals in the seed bank
        //need to check for the maximum age in the seed bank

        age = (unsigned int)current_ind->getAge(); //this is age in number of cycles (years)

        
        //first check if individual must move to seedlings because of age limit in the seed bank
        
        if( age + 1 == _popPtr->getAgeStructure()->get(0, 1) ){

          //check for germination
          if(RAND::Uniform() > survival->get(1,0)) {

            _popPtr->recycle(current_ind);
            current_patch->remove(FEM, OFFSx, k);

          } else {

            current_ind->Aging();
//            current_patch->move(FEM, OFFSx, ADLTx, k); //don't move now, will be added later, only if survived
            
            //add to temporary container, will be emptied during competition regulation seedlings
            _newSeedlingsFromSeedBank[i].push_back(current_ind);
            
            current_patch->remove(FEM, OFFSx, k);

          }

        } else {
            //this individual may stay in the seed bank or germinate

            //check for germination
            if(RAND::Uniform() > 1 - survival->get(1,0)) {

              current_ind->Aging();
//                current_patch->move(FEM, OFFSx, ADLTx, k);
              
              _newSeedlingsFromSeedBank[i].push_back(current_ind);
              
              current_patch->remove(FEM, OFFSx, k);

            //if it doesn't germinate, check if it survives in the seed bank
            // >>> this means the effective survival rate is surv * (1-germination rate)
            } else if(RAND::Uniform() > survival->get(0,0)) {

              _popPtr->recycle(current_ind);
              current_patch->remove(FEM, OFFSx, k);

            } else {
              current_ind->Aging();
            }
        }
      }
    }//end for female

    //cout<<"; germinated: "<<_newSeedlingsFromSeedBank[i].size()
      //  <<"; survived: "<<current_patch->size(FEM, OFFSx)<<endl;

  }//end for patch

  //   2. deal with new seeds: check for germination and survival in the seed bank
  
  //   seeds surviving in the seed bank need be instantiated from parents >> we create their genetics here
  
  unsigned int num_off, num_off_germinating, num_off_surviving, new_seedlings, new_seedbank;
  Patch* patch_of_origin;
  unsigned int indexOfMother;
  Individual *mother, *father, *NewOffsprg;

  for (unsigned int i = 0; i < _reducedpostDispOffspringBack[0].size(); i++) {
    
     // loop through each receiving patch (i) and apply survival and germination rates
     // to each number of individuals received from the different 'donor' patches (immigrants)

     current_patch = _popPtr->getPatch(i);

     num_off  = _postDispOffspring[i]; //this is sum(_reducedpostDispOffspringBack[0][i][...])

     num_off_germinating = 0;

     num_off_surviving = 0;

     //cout<<"  patch "<<i<<" > new seeds: "<< num_off<<flush;

     //for the current patch, cycle through all connected patches to create new offspring
     for(unsigned int k = 0, n_connect = _reducedpostDispOffspringBack[0][i].size(); k < n_connect; ++k) {

       if( _reducedpostDispOffspringBack[0][i][k] == 0) continue; //skip empty cells

       // use the 'reversed' dispersal matrix to find which is the patch of origin of the migrants
       // we subtract one because the patch numbers in the connectivity matrix are given in range [1 - num patch] in input!!!
       patch_of_origin = _popPtr->getPatch(_reducedDispMatBack[0][i][k] - 1);

       new_seedlings = RAND::Binomial2(survival->get(1,0), _reducedpostDispOffspringBack[0][i][k]);

       //effective survival is surv * (1-germination rate)
       new_seedbank = RAND::Binomial2(survival->get(0,0),
                                      _reducedpostDispOffspringBack[0][i][k] - new_seedlings);

       num_off_germinating += new_seedlings;
       
       num_off_surviving += new_seedbank;

       
       // check whether we keep numbers within total number of seeds received+produced in the patch
       assert(num_off_germinating + num_off_surviving <= num_off);


       
       // new seedlings; assumes only the last age class reproduces ----------------
       ////// FOR NOW JUST COUNT THEM - THEY WILL BE CREATED AFTER REGULATION in 'breed_newSeedlingsFromSeeds()'
       _newSeedlingsFromNewSeeds[i] += new_seedlings;

       _reducedpostRegulateOffspringBack[0][i][k] = new_seedlings;

       
       
       // MOVE SEEDS TO SEED BANK --------------------------------------------------
       // create the seeds that will stay in the seedbank
       for(unsigned int l = 0; l < new_seedbank; ++l) {

         indexOfMother = RAND::Uniform(patch_of_origin->size(FEM, last_idx )-1);

         mother = patch_of_origin->get(FEM, last_idx, indexOfMother);

         // breed assuming hermaphroditic plants (i.e. father is same as mother)
         NewOffsprg = makeOffspring( do_breed(mother, mother, i) );

         current_patch-> add(FEM, OFFSx, NewOffsprg); //stay in stage 0 = seeds

       }

     }// end for connected patches to focal patch

     //cout<<"; germinated: "<<num_off_germinating
       //  <<"("<<(double)num_off_germinating/num_off
       //  <<"); survived > seed bank: "<<num_off_surviving<<"("<<(double)num_off_surviving/num_off<<")"
       //  <<"; removed: "<<num_off - num_off_surviving - num_off_germinating<<"\n";

   }//end for patch

}
// ----------------------------------------------------------------------------------------
// do_SeedlingRegulation
// ----------------------------------------------------------------------------------------
void LCE_Breed_Disperse_Clone_Regulate::do_SeedlingRegulation()
{

  // compute the number of surviving seedlings after density-depdt regulation
  // seedlings are from three origins: seed bank, new seeds, and clones
  // for strength of competition, we need to calculate the current, post-ageing
  // population size; calculate tot pop size from post survival adult counters
  // and seedling counters
  
  unsigned int pop_size, at;

  double compet;

  unsigned int num_surviving_clones, num_surviving_newSeed, num_surviving_seedbank;

  unsigned int cnt;

  Patch* patch;

  //cout<<"+++ seedling regulation\n";

// patch loop ------------------------------------------------------------------------------
  for(unsigned int i = 0; i < _popPtr->getPatchNbr(); ++i) {

    pop_size = 0;

    cnt = 0; //only used for the //cout's below

    patch = _popPtr->getPatch(i);

    //count new seedlings from clones, there are two stages not contributing clones, seedlings and seeds
    for(unsigned int a = 0; a < _popPtr->getNumAgeClasses() - 2; ++a) {
      pop_size += _newSeedlingsFromClones[i][a];
    }

    cnt+=pop_size;
    
    //count surviving adults from each non-seed age class from previous 'year'
    for(unsigned int a =  _popPtr->getNumAgeClasses() - 1; a > 0; --a)
      pop_size += _postSurvivalAdults[i][a];

    //also need to count new seedlings in total pop size for competition
    pop_size += _newSeedlingsFromSeedBank[i].size() + _newSeedlingsFromNewSeeds[i];

    cnt +=_newSeedlingsFromSeedBank[i].size() + _newSeedlingsFromNewSeeds[i];
    
    // strength of competition in Beverton-Holt:
    compet = 1/(1 + pop_size * _competition_regulation);


    //cout<<"  patch "<<i<<" > pop size: "<<pop_size<<"; coef compet: "<<compet
      //  <<"; seedlings: "<<cnt<<endl;

    unsigned int counter = 0, counter_clones = 0; //counters for the //cout's below

    
    
    // survival within clones --------------------------------------------------------------
    for(unsigned int a = 0; a < _popPtr->getNumAgeClasses() - 2; ++a) {

       num_surviving_clones = RAND::Binomial2(compet, _newSeedlingsFromClones[i][a]);

       _newSeedlingsFromClones[i][a] = num_surviving_clones;

       counter_clones += num_surviving_clones;
    }

    counter += counter_clones;

    
    
    // survival within seedlings from the seedbank -----------------------------------------
    num_surviving_seedbank = RAND::Binomial2(compet, _newSeedlingsFromSeedBank[i].size());

    vector<Individual*>::iterator IT;

    // remove non-surviving seedlings from the temp container and delete them
    while(_newSeedlingsFromSeedBank[i].size() - num_surviving_seedbank > 0) {

      IT = _newSeedlingsFromSeedBank[i].begin();

      at = RAND::Uniform(_newSeedlingsFromSeedBank[i].size());

      delete _newSeedlingsFromSeedBank[i][at];
      //a call to metapop::recycle would be more adequate...

      _newSeedlingsFromSeedBank[i].erase(IT + at);
    }

    //what remains are the surviving ones
    counter += _newSeedlingsFromSeedBank[i].size();

    
    
    // survival within seedlings from fresh seeds ------------------------------------------
    
    _newSeedlingsFromNewSeeds[i] = 0; //reset to hold the post-regulation numbers
    
    for(unsigned int j = 0; j < _reducedpostRegulateOffspringBack[0][i].size(); ++j) {

      //the numbers within table are updated directly - we don't create a new table
      //that table was initialized in 'do_SeedSurvivalAndGermination'
      _reducedpostRegulateOffspringBack[0][i][j] = RAND::Binomial2(compet,
                                            _reducedpostRegulateOffspringBack[0][i][j]);

      _newSeedlingsFromNewSeeds[i] += _reducedpostRegulateOffspringBack[0][i][j];

    }

    counter += _newSeedlingsFromNewSeeds[i];

    //cout<<"          > survivors: "<< counter<<" ("<<counter/pop_size<<") (clones: "<<counter_clones
      //  <<" + seed bank: "<<_newSeedlingsFromSeedBank[i].size()
      //  <<" + new seeds: "<<_newSeedlingsFromNewSeeds[i]<<")"<<endl;
  }

}
// ----------------------------------------------------------------------------------------
// breed_newSeedlingsFromSeeds
// ----------------------------------------------------------------------------------------
void LCE_Breed_Disperse_Clone_Regulate::breed_newSeedlingsFromSeeds()
{

  // at this point the seedlings container should be empty and the number of surviving
  // seedlings known from the three sources: seed bank, new seeds, and clones
  // to create new seedlings from new seeds (from last year), we need the pre-aging adults
  // in the adult containers; the adult containers are updated after this breeding episode

  Patch* patch, *patch_of_origin;

  unsigned int indexOfMother;

  Individual *mother, *NewOffsprg;

  age_idx last_idx = static_cast<age_idx> (_popPtr->getNumAgeClasses() -1);


  //cout<<"+++ create new seedlings from seeds: \n";

  for(unsigned int i = 0; i < _popPtr->getPatchNbr(); ++i) {

    patch = _popPtr->getPatch(i);

    //the seedlings container should be empty by now
    assert(patch->size(FEM, ADLTx) == 0);

    // 1. move individuals from the seedbank

    for(unsigned int j = 0; j < _newSeedlingsFromSeedBank[i].size(); ++j) {

      patch->add(FEM, ADLTx, _newSeedlingsFromSeedBank[i][j]);

    }

    _newSeedlingsFromSeedBank[i].clear();

    // 2. create new seedlings from fresh seeds -- only the last age class reproduces

    for(unsigned int j = 0; j < _reducedpostRegulateOffspringBack[0][i].size(); ++j) {

      patch_of_origin = _popPtr->getPatch(_reducedDispMatBack[0][i][j] - 1);

      for(unsigned int l = 0; l < _reducedpostRegulateOffspringBack[0][i][j]; ++l) {

         indexOfMother = RAND::Uniform( patch_of_origin->size(FEM, last_idx )-1);

         mother = patch_of_origin->get(FEM, last_idx, indexOfMother);

         // breed assuming hermaphroditic plants
         NewOffsprg = makeOffspring( do_breed(mother, mother, i) );

         NewOffsprg->setAge(1);//seedling stage = age 1

         patch-> add(FEM, ADLTx, NewOffsprg); //adding to the seedling stage = first adult stage

       }

    }
    //cout<<"  patch "<<i<<" : "<<patch->size(FEM, ADLTx)<<endl;
  }
}
// ----------------------------------------------------------------------------------------
// moveAgedAdults
// ----------------------------------------------------------------------------------------
void LCE_Breed_Disperse_Clone_Regulate::moveAgedAdults()
{

  // we proceed with replacing the individuals in adult containers (
  Individual* current_ind;
  Patch* current_patch;

  unsigned int age, nb_class, ind;
  double surv_rate;

  TMatrix *survival = _popPtr->getLeslieMatrix();

  nb_class = survival->getNbCols();

  unsigned int last = nb_class - 1;

  age_idx last_idx = static_cast<age_idx> (nb_class -1);//!!age_idx is a type. type conversion.

  age_idx current_age, next_stage;


  //cout<<"+++ move aged adults to next stage\n";

  //for each patch
  for(unsigned int i = 0; i < _popPtr->getPatchNbr(); i++)
  {
    current_patch = _popPtr->getPatch(i);

////////////// LAST AGE CLASS ////////////////////////////////////////////////////

    //cout<<"  patch "<<i<<" > remove ("<<current_patch->size(FEM, last_idx)
    //	<<" - "<<_postSurvivalAdults[i][last]<<")"<<flush;
    
    //need to randomly remove non-surviving adults
    while(current_patch->size(FEM, last_idx) > _postSurvivalAdults[i][last])
    {
      ind = RAND::Uniform(current_patch->size(FEM, last_idx));

      _popPtr->recycle( current_patch->remove(FEM, last_idx, ind) );

    }

//    assert(current_patch->size(FEM, last_idx) == _postSurvivalAdults[i][last]);

    // survivors stay in last stage container
    for(unsigned int j = 0; j < _postSurvivalAdults[i][last]; ++j) {

      current_patch->get(FEM, last_idx, j)->Aging(); //this one survived, increase its age

    }
    
    
////////////////////////AGING IN OTHER AGE-CLASSES - EXCEPT OFFSPRING///////////////////////
    for( int j = (int)last - 1; j > 0; j--) {

      current_age = static_cast<age_idx> (j); //this is the index in the containers
      next_stage = static_cast<age_idx> (j+1);

      if(j == 1){ //seedlings - last age class we consider in this loop
        
        // random selection of survivors has already happened in seedling regulation stage
        
        for(unsigned int k = 0; k < _agedSeedlings[i].size(); ++k)
          current_patch->add(FEM, next_stage, _agedSeedlings[i][k]); //already aged
        
        _agedSeedlings[i].clear();
        
      }
      else {

        unsigned int non_surviving = current_patch->size(FEM, current_age) - _postSurvivalAdults[i][j];

        // randomly select survivors and move them to next stage container:
        
        while(current_patch->size(FEM, current_age) > non_surviving){

          ind = RAND::Uniform(current_patch->size(FEM, current_age));

          current_patch->get(FEM, current_age, ind)->Aging(); //this one survived, increase its age

          current_patch->move(FEM, current_age, next_stage, ind);

        }
        
        current_patch->flush(FEM, current_age, _popPtr); //ok here, remove remaining non-surviving individuals
      
      }
    }

    //cout<<"  patch "<<i<<" > age structure: 0: "<<current_patch->size(FEM, OFFSx)
      //  <<"; 1: "<<current_patch->size(FEM, ADLTx)
      //  <<"; 2: "<<current_patch->size(FEM, static_cast<age_idx> (2))
      //  <<"; 3: "<<current_patch->size(FEM, static_cast<age_idx> (3))<<endl;

  }

}
// ----------------------------------------------------------------------------------------
// breed_newSeedlingsFromClones
// ----------------------------------------------------------------------------------------
void LCE_Breed_Disperse_Clone_Regulate::breed_newSeedlingsFromClones()
{
  
  // create the surviving seedlings arising from clonal reproduction
  // for this we need the post-survival adults, thus called after 'moveAgedAdults()'
  // as survival is random respective of the phenotype of the adults, we don't care
  // which adult produced which clone during cloning stage, cloning done randomly here
  
  Individual* newind;
  Individual* mother;
  Patch* patch;
  unsigned int last = _popPtr->getNumAgeClasses() - 1;
  age_idx current_age;

  //calculate how many stages contribute clones
  //we exclude the seedling stage -> will be 2 stages only for our case
  unsigned int num_contributing_stages = _popPtr->getNumAgeClasses() - 2;

  //cout<<"+++ create new seedlings from clones\n";

  unsigned int cnt;

  for(unsigned int i = 0; i < _popPtr->getPatchNbr(); i++) {               // loop through each patch

    patch = _popPtr->getPatch(i);

    cnt = patch->size(FEM, ADLTx);

  /////WHICH AGE CLASS DOES CLONING??? -> EXCLUDE SEEDLINGS, BUT NOT PRE-ADULTS?
    for( int j = (int)last; j > 1; j--) {

      current_age = static_cast<age_idx> (j);

      for(unsigned int k = 0, stage = j - num_contributing_stages;
          k < _newSeedlingsFromClones[i][stage]; ++k)
      {
        //randomly select mother plant in current patch
        mother = patch->get(FEM, current_age, RAND::Uniform(patch->size(FEM,current_age)));

        newind = _popPtr->getNewIndividual();

        //cloning:
        (*newind) = (*mother);
//        newind->reset_counters();//MOD SLIM NEMO 7.11
        newind->setFather(NULL);
//        newind->setFatherID(0);//MOD SLIM NEMO 7.11
        newind->setMother(mother);
//        newind->setMotherID(mother->getID());//MOD SLIM NEMO 7.11
//        newind->setIsSelfed(true);//MOD SLIM NEMO 7.11
        newind->setHome(i);
        newind->setAge(1);

        patch->add(FEM,ADLTx,newind->create(false,true));
      }
    }

    //cout<<"  patch "<<i<<" : "<<patch->size(FEM, ADLTx)<<" ("<<patch->size(FEM, ADLTx)-cnt<<")"<<endl;
    //cout<<"  patch "<<i<<" > age structure: 0: "<<patch->size(FEM, OFFSx)
      //  <<"; 1: "<<patch->size(FEM, ADLTx)
      //  <<"; 2: "<<patch->size(FEM, static_cast<age_idx> (2))
      //  <<"; 3: "<<patch->size(FEM, static_cast<age_idx> (3))<<endl;

  }
}
