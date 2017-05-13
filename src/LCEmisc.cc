/** $Id: LCEmisc.cc,v 1.15.2.12 2016-04-28 13:04:42 fred Exp $
 *
 *  @file LCEmisc.cc
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
 *  created on @date 16.06.2005
 * 
 *  @author fred
 */

#include <cmath>
#include <deque>
#include "LCEmisc.h"
#include "metapop.h"
#include "Uniform.h"
#include "output.h"
#include "tstring.h"

/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

//                           ******** LCE_Aging ********

// ----------------------------------------------------------------------------------------
// LCE_Aging::execute
// ----------------------------------------------------------------------------------------
void LCE_Aging::execute ()
{  
#ifdef _DEBUG_
  message("LCE_Aging::execute (Patch nb: %i offsprg nb: %d adlt nb: %d; ",_popPtr->getPatchNbr()
          ,_popPtr->size( OFFSPRG ), _popPtr->size( ADULTS ));
#endif
  unsigned int nbInd = 0;
  Patch *patch;
  for(unsigned int i = 0; i < _popPtr->getPatchNbr(); i++) {
    
    patch = _popPtr->getPatch(i);
    
    patch->flush(ADLTx, _popPtr);
    
    nbInd = 0;
    
    while(nbInd++ < patch->get_KFem() && patch->size(FEM, OFFSx) != 0) 
      patch->move(FEM, OFFSx, ADLTx, RAND::Uniform( patch->size(FEM, OFFSx) ) );
    
    nbInd = 0;
    
    while(nbInd++ < patch->get_KMal() && patch->size(MAL, OFFSx) != 0)
      patch->move(MAL, OFFSx, ADLTx, RAND::Uniform( patch->size(MAL, OFFSx) ) );
    
    patch->flush(OFFSx, _popPtr);
    
    //set the Patch extinction and age tags:
    if(patch->size(ADULTS) == 0) {
      patch->set_isExtinct(true);
      patch->set_age(0);
    } else {
      patch->set_isExtinct(false);
      patch->set_age( patch->get_age() + 1 );
    }
  }
#ifdef _DEBUG_
  message("after: %i)\n",_popPtr->size( ));
#endif
}
/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

//                           ******** LCE_Aging_Multi ********/

// ----------------------------------------------------------------------------------------
bool LCE_Aging_Multi::setParameters ( )
{
  
  //if(_survival != NULL) delete _survival;
  
  // !! do not do this here, it is a source of errors if the leslie matrix is changed
  //using the pop parameter, e.g. with a temporal parameter, this copy will not be updated!
  
  //_survival = new TMatrix(*_popPtr->getLeslieMatrix());
  
  return true;
}
// ----------------------------------------------------------------------------------------
// LCE_Aging_Multi::execute
// ----------------------------------------------------------------------------------------
void LCE_Aging_Multi::execute ()
{  
  
  //define variables of the function
  Individual* current_ind;
  Patch* current_patch;
  
  unsigned int age, nb_class;
  double surv_rate;
  
  //this corrects the problem mentioned in the param setter above:
  _survival = _popPtr->getLeslieMatrix();
  
  nb_class = _survival->getNbCols();
  
  unsigned int last = nb_class - 1;
  
  age_idx last_idx = static_cast<age_idx> (nb_class -1);//!!age_idx is a type. type conversion.
  
  //check
  //cout<<"last id = "<<last_idx<<endl;
  //
  
  //  cout<<"+++ aging\n";
  
  age_idx current_age;
  
#ifdef _DEBUG_
  message("LCE_Aging_Multi::execute (Patch nb: %i offsprg nb: %d adlt nb: %d)\n",_popPtr->getPatchNbr()
          ,_popPtr->size( OFFSPRG ), _popPtr->size( ADULTS ));
  _popPtr->show_up();
  cout << "survival matrix (Leslie matrix):\n";
  _survival->show_up();
  cout << "age structure matrix:\n";
  _popPtr->getAgeStructure()->show_up();
#endif
  
  
  //for each patch
  for(unsigned int i = 0; i < _popPtr->getPatchNbr(); i++)
  {
    current_patch = _popPtr->getPatch(i);
    
#ifdef _DEBUG_
    message("LCE_Aging_Multi:: in patch %i (%i offsp, %i adlts)\n",
            current_patch->getID(), 
            current_patch->size(OFFSPRG),
            current_patch->size(ADULTS));
#endif
    
    ///start with the last age class, check if aggregated class:
    
    //cout<<current_patch->size(FEM, OFFSx)<<endl;
    
    //cout<<current_patch->size(FEM, last_idx+1)<<endl;
    
    // cout<<last_idx<<endl;
    
    //cout<<"age_structure0 1 2 3 4 = "<<current_patch->size(FEM, static_cast<age_idx>(last_idx-3)) <<
    //	" "<<current_patch->size(FEM, static_cast<age_idx>(last_idx-2))<<" "<<current_patch->size(FEM, static_cast<age_idx>(last_idx-1 ))<<" "
    //<<current_patch->size(FEM, static_cast<age_idx>(last_idx))<<endl;
    
    if( (surv_rate = _survival->get(last, last)) != 0) //get the last row and column of the leslie matrix
      
    {
      
      //individuals female in that class may stay if they survive
      
      for(unsigned int k = 0; k < current_patch->size(FEM, last_idx); k++)
      {
        current_ind = current_patch->get(FEM, last_idx, k);
        //check for survival
        if(RAND::Uniform() > surv_rate) 
        {
          _popPtr->recycle( current_patch->remove(FEM, last_idx, k) );
          k--;
        } else {
          current_ind->Aging(); //this one survived, increase its age
          // current_ind->mutate(); //individuals may accumulate somatic mutations
        }
        
      } //for females
      
      
      //individual male in that class may stay if they survive
      
      for(unsigned int j=0; j < current_patch->size(MAL, last_idx);j++)
      {
        current_ind=current_patch->get(MAL,last_idx,j);
        //check for survival
        if(RAND::Uniform()>surv_rate)
        {
          _popPtr->recycle(current_patch->remove(MAL,last_idx,j));
          j--;
        } else {
          current_ind->Aging();//increase age if survived}
        }
      } //for males 
      
      
    } else { //all individuals in the last age class die
      
      //cout<<"age_structureI 1 2 3 4 = "<<current_patch->size(FEM, last_idx - 2) <<
      //" "<<current_patch->size(FEM, last_idx - 1)<<" "<<current_patch->size(FEM, last_idx )<<" "
      //<<current_patch->size(FEM, last_idx+1)<<endl;
      
      
#ifdef _DEBUG_
      message("LCE_Aging_Multi::flushing last adult age class (%i)\n",last_idx);
#endif
      current_patch->flush(last_idx, _popPtr);\
      
      //cout<<"age_structureII 1 2 3 4 = "<<current_patch->size(FEM, static_cast<age_idx>(last_idx-3)) <<
      //" "<<current_patch->size(FEM, static_cast<age_idx>(last_idx-2))<<" "<<current_patch->size(FEM, static_cast<age_idx>(last_idx-1 ))<<" "
      //<<current_patch->size(FEM, static_cast<age_idx>(last_idx))<<endl;
      
    }
    
    
    
    ///Regulation in other age-classes - EXCEPT OFFSPRING
    
    for( int j = (int)nb_class - 2; j > 0; j--) {
      
      current_age = static_cast<age_idx> (j); //this is the index in the containers
      
#ifdef _DEBUG_
      message("  aging in age class (%i) ",current_age);
      unsigned int move_cntr = 0, aging_cntr = 0, candidates = 0, mean_age = 0, ind_cntr = 0;
#endif
      
      
      //for each female///
      
      for(int k = current_patch->size(FEM, current_age) -1; k >= 0; --k) {
        
        current_ind = current_patch->get(FEM, current_age, k);
        
        age = (unsigned int)current_ind->getAge(); //this is age in generation/cycle number
        
        
#ifdef _DEBUG_
        if(j>0) {mean_age += age;ind_cntr++;}
#endif 
        
        
        //first check if individual must move to upper age class
        if( age + 1 >= _popPtr->getAgeStructure()->get(0, j+1) )
          
        {
        	//check for survival
          
#ifdef _DEBUG_
          candidates++;
#endif
          
          
          if(RAND::Uniform() > _survival->get(j+1, j)) {
            
            _popPtr->recycle(current_ind);
            current_patch->remove(FEM, current_age, k);
            
          } else {
            current_ind->Aging();
            current_patch->move(FEM, current_age, static_cast<age_idx> (current_age + 1), k);
#ifdef _DEBUG_
            move_cntr++; aging_cntr++;
#endif
          }
        }
        
        else   //individuals may stay in the same age class if they survive
          
        {
          
          if(_survival->get(j, j) == 0){
            
            //that means the individual hasn't reached the max age, but cannot survive in the same stage
            //we move it to the next stage if it survives (instead of killing it instantly):
            
            if(RAND::Uniform() > _survival->get(j+1, j)) {
              
              _popPtr->recycle(current_ind);
              current_patch->remove(FEM, current_age, k);
              
            } else {
              current_ind->Aging();
              current_patch->move(FEM, current_age, static_cast<age_idx> (current_age + 1), k);
            }
            
            //otherwise, the individual can stay in the same stage and survive there
          } else if(RAND::Uniform() > _survival->get(j, j)) {
            
            _popPtr->recycle(current_ind);
            current_patch->remove(FEM, current_age, k);
            
          } else {
            current_ind->Aging();
#ifdef _DEBUG_
            aging_cntr++;
#endif
          }
        }
        
      }
      
      
      ///for each male
      
      for(int k = current_patch->size(MAL, current_age) -1; k >= 0; --k) {
        
        current_ind = current_patch->get(MAL, current_age, k);
        
        age = (unsigned int)current_ind->getAge(); //this is age in generation/cycle number
        
        //first check if individual must move to upper age class
        if( age + 1 >= _popPtr->getAgeStructure()->get(0, j+1) ) {
          
          //check for survival
          if(RAND::Uniform() > _survival->get(j+1, j)) {
            
            _popPtr->recycle(current_ind);
            current_patch->remove(MAL, current_age, k);
            
          } else {
            current_ind->Aging();
            current_patch->move(MAL, current_age, static_cast<age_idx> (current_age + 1), k);
          }
        }
        else   //individuals may stay in the same age class if they survive
        {
          if(_survival->get(j, j) == 0){
            
            //that means the individual hasn't reached the max age, but cannot survive in the same stage
            //we move it to the next stage if it survives (instead of killing it instantly):
            
            if(RAND::Uniform() > _survival->get(j+1, j)) {
              
              _popPtr->recycle(current_ind);
              current_patch->remove(MAL, current_age, k);
              
            } else {
              current_ind->Aging();
              current_patch->move(MAL, current_age, static_cast<age_idx> (current_age + 1), k);
            }
            
            //otherwise, the individual can stay in the same stage and survive there
          } else if(RAND::Uniform() > _survival->get(j, j)) {
            
            _popPtr->recycle(current_ind);
            current_patch->remove(MAL, current_age, k);
            
          } else {
            current_ind->Aging();
          }
        }
      }
      
#ifdef _DEBUG_
      message("(%f): moved %i(%i) ind from %i to %i (%i aged)\n", 
              (double)mean_age/ind_cntr, move_cntr, candidates, 
              current_age, current_age+1, aging_cntr);
#endif
    } //end_for j age classes
    
    
    //REGULATION IN OFFSPRING
    //(allows maturation rate - some offspring may stay offspring)
    
    //Option to regulate to carrying capacity - offspring survive according to density
    //occurs before aging
    
    // DO WE KEEP THIS ??
    double survival_K;
    
    if (get_parameter("regulation_before_aging")->isSet()){
      
    	unsigned int adult_size=current_patch->size(ADULTS);
      
    	double K = current_patch->get_K(FEM) * (nb_class-1); //other classes than offspring
      
      
    	if(adult_size > K) survival_K = 0;
      
    	else survival_K = (K - adult_size) / K;
      
    	/*cout<<"K "<<K<<endl;
       cout<<"Asize "<<adult_size<<endl;*/
    }
    
    else survival_K = 1;
    
    /*cout<<"S_K "<<survival_K<<endl;
     cout<<current_patch->size(FEM, OFFSx)<<endl;
     cout<<"gen"<<_popPtr->getCurrentGeneration();
     cout<<"\n"<<endl;*/
    
    ///end of Option
    
    //@TODO aging of MAL offspring
    
    //  cout<<"  patch "<<i<<" > offspring: "<<current_patch->size(FEM, OFFSx)<<flush;
    
    unsigned int cntr = 0, aged = 0, surv_in_bank = 0;
    
    for(int k = (int)current_patch->size(FEM, OFFSx) - 1; k >= 0; --k) {
      
      current_ind = current_patch->get(FEM, OFFSx, k);
      
      ///// NO SEED BANK /////
      
      if(_survival->get(0,0) == 0){
        
        if(RAND::Uniform() > _survival->get(1, 0)) {
          _popPtr->recycle(current_ind);
          current_patch->remove(FEM, OFFSx, k);
        }
        
        // that shouldn't be necessary anymore...
        else if(RAND::Uniform() > survival_K) {
          //cout<<"dead"<<endl;
          _popPtr->recycle(current_ind);
          current_patch->remove(FEM, OFFSx, k);
        }
        
        else {
          current_ind->Aging();
          current_patch->move(FEM, OFFSx, ADLTx, k);  //move to first adult class
          cntr++;
        }
        
      }
      
      ///// SEED BANK /////
      
      if( _survival->get(0,0) != 0) {
        
        age = (unsigned int)current_ind->getAge(); //this is age in number of cycles (years)
        
        //first check if individual must move to seedlings because of age limit in the seed bank
        if( age + 1 == _popPtr->getAgeStructure()->get(0, 1) ){
          
          //check for survival to next age class
          if(RAND::Uniform() > _survival->get(1,0)) {
            
            _popPtr->recycle(current_ind);
            current_patch->remove(FEM, OFFSx, k);
            
          } else {
            
            current_ind->Aging();
            current_patch->move(FEM, OFFSx, ADLTx, k);
            current_patch->remove(FEM, OFFSx, k);
            aged++;
          }
          
        } else {
          
          //note that in this case Survival off-adlt is the maturation rate//
          
          if(RAND::Uniform() > 1 - _survival->get(1, 0)  ){
            current_ind->Aging();
            current_patch->move(FEM, OFFSx, ADLTx, k); //move to next stage
            cntr++;
          }
          // remove ind if not survive
          else if(RAND::Uniform() > _survival->get(0, 0)){
            _popPtr->recycle(current_ind);
            current_patch->remove(FEM,OFFSx, k);
          } else {
            // otherwise it just stays there, but we need to age it
            current_ind->Aging();
            surv_in_bank++;
          }
        }
      }
    }
    
    //  cout<<"; germinated: "<<cntr<<"; aged in bank: "<<surv_in_bank<<"; aged out of bank: "<<aged<<"; in bank: "<<current_patch->size(FEM, OFFSx)<<endl;
    
    
    //for(unsigned int k = 0; k < current_patch->size(FEM, OFFSx); ++k) {
    
    //current_ind = current_patch->get(FEM, OFFSx, k);
    
    //if(RAND::Uniform() > survival_K) {
    // _popPtr->recycle(current_ind);
    // current_patch->remove(FEM, OFFSx, k);
    // k--;}}
    
    
    
    
    
#ifdef _DEBUG_
    message("LCE_Aging_Multi::finished aging in patch %i (%i offsp, %i adlts)\n",
            current_patch->getID(),
            current_patch->size(OFFSPRG),
            current_patch->size(ADULTS));
#endif
    
    //ceiling regulation !!! this need fixing for errors !!! not the most efficient way
    //if(current_patch->size(FEM, ADULTS) > current_patch->get_K(FEM)) {
    
#ifdef _DEBUG_
    message("LCE_Aging_Multi::ceiling regulation on adults in patch %i (%i fem)\n",
            current_patch->getID(), current_patch->size(FEM, ADULTS));
#endif
    
    //unsigned int size = current_patch->size(FEM, ADULTS);
    //unsigned int ceiling = current_patch->get_K(FEM);
    
    //double age_prop[nb_class - 1]; //the offspring are not regulated here
    //double sum, rand;
    
    //for(unsigned int a = 1; a < nb_class; a++) //only consider adult classes:
    
    //proportion of adult with age a in the pop. index starts at 0
    
    //age_prop[a-1] = (double)current_patch->size(FEM, static_cast<age_idx> (a)) / size;
    
    //while the number of individual in the focal patch is above the carrying capacity
    
    // while(current_patch->size(FEM, ADULTS) > ceiling) {
    
    //  sum = age_prop[0]; //initialize the counter
    // age = 0;
    // rand = RAND::Uniform();
    
    //problem here when the age class targeted has been emptied
    //but the proportion has not been updated to 0
    //==>will produce an error at runtime
    
    // while(rand > sum && age < nb_class -1) {
    //  age++;
    //  sum += age_prop[ age ]; //equivalent to sum=sum + age_prop
    //}
    
    //chose randomly an individual within the focal age-class and patch and remove it
    
    //current_age = static_cast<age_idx> (age+1);
    //ind = RAND::Uniform( current_patch->size( FEM, current_age ) );
    //_popPtr->recycle( current_patch->remove(FEM, current_age, ind) );
    
    
    //}
    // }
    
    //set the Patch extinction and age tags:
    if(current_patch->size(ALL) == 0) {
      current_patch->set_isExtinct(true);
      current_patch->set_age(0);
      
    } else {
      current_patch->set_isExtinct(false);
      current_patch->set_age( current_patch->get_age() + 1 );
    }
    
    // cout<<"age_structureIII 1 2 3 4 = "<<current_patch->size(FEM, static_cast<age_idx>(last_idx-3)) <<
    //" "<<current_patch->size(FEM, static_cast<age_idx>(last_idx-2))<<" "<<current_patch->size(FEM, static_cast<age_idx>(last_idx-1 ))<<" "
    //<<current_patch->size(FEM, static_cast<age_idx>(last_idx))<<endl;
    
    
  } //end_for i patchNbr
  
#ifdef _DEBUG_
  cout<<"   age classes after regulation: ";
  for(unsigned int i =0; i < nb_class; i++)
    cout<<i<<":"<<_popPtr->size(FEM,static_cast<age_idx> (i))<<" | ";
  cout<<endl;
  _popPtr->show_up();
#endif
  
}

/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

//                             ******** LCE_Regulation ********

// ----------------------------------------------------------------------------------------
// LCE_Regulation::setParameters
// ----------------------------------------------------------------------------------------
bool LCE_Regulation::setParameters ()
{
 
  if(! get_parameter("regulation_by_competition")->isSet() && 
     ! get_parameter("regulation_carrying_capacity")->isSet())
    return error("regulation mode not set for LCE regulation\n");
  
  _regulate_ptrs.clear();
  
  // defaults:
  _competition_regulation = 0;
  _age_flag = ADULTS; // = 31
  
  // set the parameters for regulation by competition
	if(get_parameter("regulation_by_competition")->isSet()) {

		_competition_regulation = get_parameter_value("regulation_by_competition");
	
    // regulation can act on different age classes, it acts on the offspring by default
    if(get_parameter("regulation_by_competition_affected_age")->isSet()) {
      
      _age_target = static_cast<age_idx>((unsigned int)get_parameter_value("regulation_by_competition_affected_age"));
      
      if(_age_target > _popPtr->getNumAgeClasses() -1)
        return(error("parameter \"regulation_by_competition_affected_age\" targets a non-existing age class\n"));
      
    } else
      _age_target = OFFSx;
    
    cout<<"DOING REGULATION ON "<<_age_target<<" (offsx: "<<OFFSx<<"; adltx: "<<ADLTx<<")\n";
    
    // the competition is exerted on the targeted age class by individuals in all other age classes,
    // we need to set the age flag corresponding to those age classes if different from all adult
    // age classes; age flag 31 = 0xFFFFFFFE, is all age classes except the offspring
    // (the max is 32 age classes, see the definition of the ADULTS age flag in types.h)
    if(get_parameter("regulation_by_competition_count_age_flag")->isSet()) {
      
      _age_flag = (unsigned int)get_parameter_value("regulation_by_competition_count_age_flag");
      
      if(_age_flag > 31)
        return(error("in \"regulation_by_competition\", the age flag is > 31\n"));
      
    }
    
    _regulate_ptrs.push_back( &LCE_Regulation::regulatePatch );
    
    
  } // end regulation_by_competition
   
  // for ceiling regulation, we add the pointer to the set of regulation function
  if (get_parameter("regulation_carrying_capacity")->isSet()) {
    
    _regulate_ptrs.push_back( &LCE_Regulation::regulatePatchCeiling );
  
  }
  
	return true;
}
// ----------------------------------------------------------------------------------------
// LCE_Regulation::execute
// ----------------------------------------------------------------------------------------
void LCE_Regulation::execute ()
{
#ifdef _DEBUG_
  message("LCE_Regulation::execute (Patch nb: %i, offsprg: %i, adults: %i ",_popPtr->getPatchNbr()
          ,_popPtr->size( OFFSPRG ), _popPtr->size( ADULTS ));
#endif
  Patch* patch;
  
  _nb_class = _popPtr->getNumAgeClasses();
  
  unsigned int pop_size;
  
  for(unsigned int i = 0; i < _popPtr->getPatchNbr(); i ++) {
    
    patch = _popPtr->getPatch(i);
    
    pop_size = patch->size(_age_flag); //total adult size matters for the strength of competition
    
    //if(!pop_size) continue; //go to next patch if this one is empty
    
    for (unsigned int i = 0; i < _regulate_ptrs.size(); ++i) {
      (this->* _regulate_ptrs[i])(patch, FEM, pop_size);
      (this->* _regulate_ptrs[i])(patch, MAL, pop_size);
    }
    
  }
  
#ifdef _DEBUG_
  message("after: offsprg: %i, adults: %i)\n",_popPtr->size( OFFSPRG ), _popPtr->size( ADULTS ));
#endif
  
}
// ----------------------------------------------------------------------------------------
// LCE_Regulation::regulatePatch
// ----------------------------------------------------------------------------------------
void LCE_Regulation::regulatePatch (Patch* patch, sex_t sex, unsigned int pop_size)
{  
  unsigned int num_affected = patch->size(sex, _age_target);
  
  if(!num_affected) return; //exit if sex-age class container is empty
  
  //option by default: beverton holt regulation
  
  if(_competition_regulation != 0) {
    
    double compet = 1/(1 + pop_size * _competition_regulation);
    
    unsigned int num_surviving = RAND::Binomial2(compet,num_affected);
    
    while(patch->size(sex, _age_target) > num_surviving){
      
      _popPtr->recycle( patch->remove(sex, _age_target,
                                      RAND::Uniform( patch->size(sex, _age_target) ) ) );
      
    }
  }

}//end function
// ----------------------------------------------------------------------------------------
// LCE_Regulation::regulatePatchCeiling
// ----------------------------------------------------------------------------------------
void LCE_Regulation::regulatePatchCeiling (Patch* patch, sex_t sex, unsigned int pop_size)
{
  // adjust to carrying capacity randomly removing individuals
  _K = patch->get_K(sex);
  
  age_idx age;
  
  if(_K == 0) {
    
    if(patch->size(sex, ALL) != 0) {
      
        patch->flush(sex, _popPtr);
    
    }
  }
  else {
    while( patch->size(ALL) > _K ) {
      
      age = static_cast<age_idx>(RAND::Uniform(_nb_class));
      
      if(patch->size(sex, age) !=0 ){
        
        _popPtr->recycle( patch->remove(sex, age,
                                        RAND::Uniform( patch->size(sex, age) ) ) );
        
      }
    }
  }
  
}//end function

/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

//                         ******** LCE_Patch_Extinction ********

//-----------------------------------------------------------------------------
// LCE_Patch_Extinction
//-----------------------------------------------------------------------------
LCE_Patch_Extinction::LCE_Patch_Extinction( ) : LifeCycleEvent("extinction",""), _Xtion_rate(0),
_harvest_size(0), _harvest_proportion(0), _harvest_size_varies(0), _by_size(0), 
_by_proportion(0), _harvest_dist_stdev(0), _extinction_threshold(0), _rand_size_fct(0)
{
  ParamUpdater< LCE_Patch_Extinction > * updater =
  new ParamUpdater< LCE_Patch_Extinction > (&LCE_Patch_Extinction::setParameters);
  add_parameter("extinction_rate", DBL, 0, 1, 0, 1, updater);
  add_parameter("extinction_threshold", DBL, 0, 1, 0, 1, updater);
  add_parameter("extinction_size", INT, 0, 0, 0, 0, updater);
  add_parameter("extinction_proportion", DBL, 0, 1, 0, 1, updater);
  add_parameter("extinction_size_distribution", STR, 0, 0, 0, 0, updater);
  add_parameter("extinction_size_dist_stdev", DBL, 0, 0, 0, 0, updater);
  add_parameter("extinction_size_dist_shape", DBL, 0, 0, 0, 0, updater);
}

//-----------------------------------------------------------------------------
// LCE_Patch_Extinction::init
//-----------------------------------------------------------------------------
bool LCE_Patch_Extinction::setParameters ()
{
  
  if(get_parameter("extinction_rate")->isSet()) {
    
    if(!_Xtion_rate) _Xtion_rate = new TMatrix();
    
    if(!set_matrix_param(_Xtion_rate, "extinction_rate")) return false;
    
  } else {
    if(_Xtion_rate) delete _Xtion_rate; _Xtion_rate = 0;
  }
  
  if(get_parameter("extinction_size")->isSet()) {
    
    if(!_harvest_size) _harvest_size = new TMatrix();
    
    if(!set_matrix_param(_harvest_size, "extinction_size")) return false;
    
    _by_size = true;
    
  } else {
    if(_harvest_size) delete _harvest_size; _harvest_size = 0;
    _by_size = false;
  }
  
  if(get_parameter("extinction_proportion")->isSet()) {
    if(!_harvest_proportion) _harvest_proportion = new TMatrix();
    
    if(!set_matrix_param(_harvest_proportion, "extinction_proportion")) return false;
    
    _by_proportion = true;
    
  } else {
    if(_harvest_proportion) delete _harvest_proportion; _harvest_proportion = 0;
    _by_proportion = false;
  }
  
  _extinction_threshold = get_parameter_value("extinction_threshold");
  
  if( !_Xtion_rate && !_by_size && !_by_proportion) {
    error("Please give one of the following parameter: \"extinction_rate\", \"extinction_size\", or \"extinction_proportion\".\n");
    return false;
  }
  else if(_by_size && _by_proportion) {
    warning("Both \"extinction_size\" and \"extinction_proportion\" are set, using sizes only.\n");
    _by_proportion = false;
  }
  
  if(get_parameter("extinction_size_distribution")->isSet()) {
    
    if(!_by_size) {
      error("\"extinction_size_distribution\" is set but the \"extinction_size\" parameter is not!\n");
      return false;
    }
    
    _harvest_distribution = _paramSet->getArg("extinction_size_distribution");
    _harvest_size_varies = true;
    _harvest_dist_stdev = get_parameter_value("extinction_size_dist_stdev");
    
    if(_harvest_distribution.compare("poisson") == 0)
      
      _rand_size_fct = &LCE_Patch_Extinction::rand_poisson;
    
    else if(_harvest_distribution.compare("uniform") == 0)
      
      _rand_size_fct = &LCE_Patch_Extinction::rand_uniform;
    
    else if(_harvest_distribution.compare("normal") == 0) {
      
      _rand_size_fct = &LCE_Patch_Extinction::rand_gaussian;
      
      if(_harvest_dist_stdev == -1) {
        error("Standard deviation of the normal distribution for the harvesting size distribution is missing!\n");
        return false;
      }
      
    } else if(_harvest_distribution.compare("lognormal") == 0) {
      
      _rand_size_fct = &LCE_Patch_Extinction::rand_lognormal;
      
      if(_harvest_dist_stdev == -1) {
        error("Standard deviation of the lognormal distribution for the harvesting size distribution is missing!\n");
        return false;
      }
      
    } else if(_harvest_distribution.compare("exponential") == 0)
      
      _rand_size_fct = &LCE_Patch_Extinction::rand_exp;
    
    else {
      error("Distribution \"%s\" is not a valid option for \"harvest_size_distribution\"\n",
            _harvest_distribution.c_str());
      return false;
    }
    
  }
  return true;
}
//-----------------------------------------------------------------------------
// LCE_Patch_Extinction::set_matrix_param
//-----------------------------------------------------------------------------
bool LCE_Patch_Extinction::set_matrix_param (TMatrix *mat, string name)
{
  double value;
  Param* param = get_parameter(name);
  
  if(param->isMatrix()) {
    
    param->getMatrix(mat);
    
    if(mat->getNbRows() > 1) {
      error("The \"%s\" matrix must be a one-dimensional array.\n", name.c_str());
      return false;
    }
    
    if(mat->getNbCols() != _popPtr->getPatchNbr()) {
      error("The length of the \"%s\" array must be equal to the number of patches.\n", name.c_str());
      return false;
    }
    
  } else {
    value = param->getValue();
    mat->reset(1, _popPtr->getPatchNbr());
    mat->assign(value);
  }
  return true; 
}
//-----------------------------------------------------------------------------
// LCE_Patch_Extinction::execute
//-----------------------------------------------------------------------------
void LCE_Patch_Extinction::execute ()
{
#ifdef _DEBUG_
  message("LCE_Patch_Extinction::execute ");
  unsigned int cnt = 0;
#endif
  Patch *patch;
  for(unsigned int i = 0; i < _popPtr->getPatchNbr(); i++) {
    
    patch = _popPtr->getPatch(i);
    
    if(_by_size || _by_proportion) {
      do_remove(OFFSx, patch);
      do_remove(ADLTx, patch);
    } else if(_Xtion_rate)
      if( RAND::Uniform() < _Xtion_rate->get(0, i) )
        do_flush(patch);
    
    if(_extinction_threshold != -1) {
      if( _extinction_threshold < 1 && (double)patch->size(ALL)/patch->get_K() < _extinction_threshold ) 
        do_flush(patch);
      else if( patch->size(ALL) < _extinction_threshold )
        do_flush(patch);
    }
#ifdef _DEBUG_
    cnt += (patch->get_isExtinct());
#endif
  }
#ifdef _DEBUG_
  message("(%i extinct patches)\n",cnt);
#endif  
}
//-----------------------------------------------------------------------------
// LCE_Patch_Extinction::do_flush
//-----------------------------------------------------------------------------
void LCE_Patch_Extinction::do_flush (Patch *patch)
{
  patch->flush(_popPtr);
  patch->set_isExtinct(true);
  patch->set_age(0);
}
//-----------------------------------------------------------------------------
// LCE_Patch_Extinction::do_remove
//-----------------------------------------------------------------------------
void LCE_Patch_Extinction::do_remove(age_idx AGE, Patch* patch)
{
  unsigned int remove_size;
  sex_t sex;
  
  //check if probability of event is set, and if removal will happen
  if(_Xtion_rate) if( RAND::Uniform() > _Xtion_rate->get(0, patch->getID()) ) return;
  
  if(patch->size(AGE) != 0) {
    
    remove_size = get_harvest_size(AGE, patch);
    
    for(unsigned int i = 0; i < remove_size; ++i) {
      
      sex = (sex_t)RAND::RandBool(); 
      
      if(patch->size(sex, AGE) != 0)
        _popPtr->recycle( patch->remove( sex, AGE, (unsigned int)RAND::Uniform(patch->size(sex, AGE)) ) );
      else //we already know here that the patch is not empty
        _popPtr->recycle( patch->remove( (sex_t)!sex, AGE, (unsigned int)RAND::Uniform(patch->size( (sex_t)!sex, AGE) ) ) );
      
      if(patch->size(AGE) == 0) break;
    }
    //    cout<<"--removed "<<remove_size<<" individuals in age class "<<AGE<<", patch "<< patch->getID()<<" size = "<<patch->size(AGE)<<endl;
  }
}
//-----------------------------------------------------------------------------
// LCE_Patch_Extinction::get_harvest_size
//-----------------------------------------------------------------------------
unsigned int LCE_Patch_Extinction::get_harvest_size (age_idx AGE, Patch *patch)
{
  
  if( _by_size ) {
    
    if(_harvest_size_varies)  
      return (this->*_rand_size_fct) (_harvest_size->get(0, patch->getID()));
    
    else return (unsigned int)_harvest_size->get(0, patch->getID());
    
  } else if( _by_proportion ) {
    
    return (unsigned int)(_harvest_proportion->get(0, patch->getID()) * patch->size(AGE));
    
  }
  
  return 0;
}

/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

//                             ******** LCE_Cross ********
// ----------------------------------------------------------------------------------------
// LCE_Cross::execute
// ----------------------------------------------------------------------------------------
LCE_Cross::LCE_Cross() : LifeCycleEvent("cross",""), _nSire(0), _nDam(0), _nOffspring(0), 
_atGeneration(0)
{
  ParamUpdater< LCE_Cross > * updater = new ParamUpdater< LCE_Cross > (&LCE_Cross::setParameters);
  add_parameter("cross_num_sire", INT, 1, 0, 0, 0, updater);
  add_parameter("cross_num_dam", INT, 1, 0, 0, 0, updater);
  add_parameter("cross_num_offspring", INT, 1, 0, 0, 0, updater);
  add_parameter("cross_at_generation", INT, 1, 0, 0, 0, updater);
  add_parameter("cross_do_among_pop", BOOL, 0, 0, 0, 0, updater);
  add_parameter("cross_do_within_pop", BOOL, 0, 0, 0, 0, updater);
  add_parameter("cross_with_replacement", BOOL, 0, 0, 0, 0, updater);
}
// ----------------------------------------------------------------------------------------
// LCE_Cross::init
// ----------------------------------------------------------------------------------------
bool LCE_Cross::setParameters ( )
{
  _nSire = (unsigned int)get_parameter_value("cross_num_sire");
  _nDam = (unsigned int)get_parameter_value("cross_num_dam");
  _nOffspring = (unsigned int)get_parameter_value("cross_num_offspring");
  _atGeneration = (unsigned int)get_parameter_value("cross_at_generation");
  
  if(get_parameter("cross_do_among_pop")->isSet())
    _doAmongPop = (unsigned int)get_parameter_value("cross_do_among_pop");
  else
    _doAmongPop = 0;
  
  string arg = get_parameter("cross_do_within_pop")->getArg();
  if (tstring::str2int(arg) == 1) {
    _doWithinPop = 1;
  } else if (tstring::str2int(arg) == 0) {
    _doWithinPop = 0;
  } else {
    _doWithinPop = 1;
  }
  
  //  if(get_parameter("cross_do_within_pop")->isSet())
  //    _doWithinPop = (unsigned int)get_parameter_value("cross_do_within_pop");
  //  else
  //    _doWithinPop = 1;
  
  if(get_parameter("cross_with_replacement")->isSet())
    _doReplace = (unsigned int)get_parameter_value("cross_with_replacement");
  else 
    _doReplace = 0;
  
  return true;
}
// ----------------------------------------------------------------------------------------
// LCE_Cross::execute
// ----------------------------------------------------------------------------------------
void LCE_Cross::execute ()
{
  Patch* patch;
  unsigned int nsire;
  deque<Individual*> males, females;
  
#ifdef _DEBUG_
  message("LCE_Cross::execute\n");
#endif
  
  if(_popPtr->getCurrentGeneration() != _atGeneration) return;
  
  if(_popPtr->size(OFFSPRG) != 0) {
    warning("Offspring are still present at time of crossing, flushing\n");
    _popPtr->flush(OFFSPRG);
  }
  
  for(unsigned int i = 0; i < _popPtr->getPatchNbr(); i ++) {
    
    patch = _popPtr->getPatch(i);
    
    if(patch->size(MAL, ADLTx) == 0) continue;
    
    if(_nSire > patch->size(MAL, ADLTx)) {
      warning("LCE_Cross:: num_sire greater than actual number of males in patch %i, reseting.\n", patch->size(MAL, ADLTx));
      nsire = patch->size(MAL, ADLTx);
    } else
      nsire = _nSire;
    
    males.clear();
    patch->getCopy(MAL, ADLTx, males);
    
    if(_doAmongPop) sampleAmongPop(patch, males, nsire);
    
    males.clear();
    females.clear();
    
    if(patch->size(FEM, ADLTx) == 0) continue;
    
    patch->getCopy(MAL, ADLTx, males);
    patch->getCopy(FEM, ADLTx, females);
    
    if(_doWithinPop) sampleWithinPop(patch, males, females, nsire);
  }
  
}
// ----------------------------------------------------------------------------------------
// LCE_Cross::sampleAmongPop
// ----------------------------------------------------------------------------------------
void LCE_Cross::sampleAmongPop(Patch* patch, deque<Individual*>& males, unsigned int nsire)
{
  unsigned int at, deme, npatch = _popPtr->getPatchNbr();
  Individual *sire, *dam;
  Patch* aimedPatch;
  
  for(unsigned int s = 0; s < nsire; ++s) {
    
    //sampling parents
    at = (unsigned int)RAND::Uniform( males.size() );
    
    sire = males[at];
    
    if( !_doReplace ) males.erase(males.begin() + at);
    
    for(unsigned int d = 0; d < _nDam; ++d) {
      //sample a deme
      do {
        deme = (unsigned int)RAND::Uniform(npatch);
      } while(deme == patch->getID() || patch->size(FEM, ADLTx) == 0);
      
      aimedPatch = _popPtr->getPatchPtr(deme);
      //sample a female
      dam = aimedPatch->get( FEM, ADLTx, (unsigned int)RAND::Uniform( aimedPatch->size(FEM, ADULTS) ) );
      
      if(dam == NULL) {d--; continue;}
      
      //mate and breed
      for(unsigned int k = 0; k < _nOffspring; ++k) {
        
        patch->add(FEM, OFFSx, _popPtr->makeOffsprg(dam, sire, FEM, patch->getID()) );
        
      }
    }
  }
}
// ----------------------------------------------------------------------------------------
// LCE_Cross::sampleWithinPop
// ----------------------------------------------------------------------------------------
void LCE_Cross::sampleWithinPop(Patch* patch, deque<Individual*>& males, deque<Individual*>& females, unsigned int nsire)
{  
  unsigned int ndam, at;
  Individual *sire, *dam;
  
  if(!_doReplace && nsire * _nDam > patch->size(FEM,ADLTx)) {
    warning("LCE_Cross:: total num of dam greater than available females in patch %i, reseting.\n",patch->getID());
    ndam = patch->size(FEM,ADLTx)/nsire;
  } else
    ndam = _nDam;
  
  if(nsire < 2 || ndam == 0) return;
  
  for(unsigned int s = 0; s < nsire; ++s) {
    
    //sampling parents
    at = (unsigned int)RAND::Uniform( males.size() );
    
    sire = males[at];
    
    if( !_doReplace ) males.erase(males.begin() + at);
    
    for(unsigned int d = 0; d < ndam; ++d) {
      
      at = (unsigned int)RAND::Uniform( females.size() );
      
      dam = females[at];
      
      if( !_doReplace ) females.erase(females.begin() + at);
      
      //now mate and breed
      for(unsigned int k = 0; k < _nOffspring; ++k) {
        
        patch->add(FEM, OFFSx, _popPtr->makeOffsprg(dam, sire, FEM, patch->getID()) );
        
      }
    }
  }
}
/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

//                             ******** LCE_Resize ********/
// ----------------------------------------------------------------------------------------
// LCE_Resize::execute
// ----------------------------------------------------------------------------------------
LCE_Resize::LCE_Resize() : LifeCycleEvent("resize",""), _atGeneration(0),
_do_flush(0), _do_fill(0), _do_regulate(0)
{
  ParamUpdater< LCE_Resize > * updater = new ParamUpdater< LCE_Resize > (&LCE_Resize::setParameters);
  add_parameter("resize_at_generation", INT, 1, 0, 0, 0, updater);
  updater = new ParamUpdater< LCE_Resize > (&LCE_Resize::updateParameters);
  add_parameter("resize_patch_number", INT, 0, 0, 0, 0, updater);
  add_parameter("resize_patch_capacity", INT, 0, 0, 0, 0, updater);
  add_parameter("resize_female_capacity", INT, 0, 0, 0, 0, updater);
  add_parameter("resize_male_capacity", INT, 0, 0, 0, 0, updater);
  add_parameter("resize_age_class", STR, 0, 0, 0, 0, updater);
  add_parameter("resize_do_flush", BOOL, 0, 0, 0, 0, updater);
  add_parameter("resize_do_regulate", BOOL, 0, 0, 0, 0, updater);
  add_parameter("resize_do_fill", BOOL, 0, 0, 0, 0, updater);
  add_parameter("resize_keep_patch", MAT, 0, 0, 0, 0, updater);
}
// ----------------------------------------------------------------------------------------
// LCE_Resize::setParameters
// ----------------------------------------------------------------------------------------
bool LCE_Resize::setParameters ( )
{
  
  _generations.clear();
  
  if(_paramSet->isMatrix("resize_at_generation") ) {
    
    TMatrix tmp;
    
    _paramSet->getMatrix("resize_at_generation", &tmp);
    
    for(unsigned int i = 0; i < tmp.getNbCols(); i++)
      _generations.push_back((int)tmp.get(0, i));
    
  } else
    _generations.push_back((int)get_parameter_value("resize_at_generation"));
  
  _genITER = _generations.begin();
  _atGeneration =  (*_genITER);
  
  return updateParameters();
}
// ----------------------------------------------------------------------------------------
// LCE_Resize::updateParameters
// ----------------------------------------------------------------------------------------
bool LCE_Resize::updateParameters ( )
{
  string option = _paramSet->getArg("resize_age_class");
  
  if(option.length() == 0) _setAge = ALL;
  else if(option.compare("OFFSPRG") == 0 || option.compare("offspring") == 0 ||
          option.compare("0") == 0){
    _setAge = OFFSPRG;
  } else if(option.compare("ADULTS") == 0 || option.compare("adults") == 0 ||
            option.compare("1") == 0) {
    _setAge = ADULTS;
  } else if(option.compare("ALL") == 0 || option.compare("all") == 0) {
    _setAge = ALL;
  } else {
    warning("\"%s\" is not a valid option for parameter \"resize_age_class\".\n",option.c_str());
    _setAge = ALL;
  }
  
  if( !get_parameter("resize_patch_number")->isSet() &&
     !get_parameter("resize_patch_capacity")->isSet() &&
     !get_parameter("resize_female_capacity")->isSet() &&
     !get_parameter("resize_male_capacity")->isSet() &&
     !get_parameter("resize_keep_patch")->isSet()) {
    error("LCE_Resize:: at least one of the population size parameters must be specified.\n");
    return false;
  }
  
  _do_regulate = get_parameter("resize_do_regulate")->isSet();
  _do_flush = get_parameter("resize_do_flush")->isSet();
  _do_fill = get_parameter("resize_do_fill")->isSet();
  
  if(get_parameter("resize_keep_patch")->isSet()) {
    _paramSet->getMatrix("resize_keep_patch", &_patch2keep);
    if(_patch2keep.getNbCols() > _popPtr->getPatchNbr() && 
       _popPtr->getCurrentGeneration() == (unsigned)_atGeneration) {
      error("LCE_Resize:: more patches to keep than existing patches in the population!\n");
      return false;
    }
  } else
    _patch2keep.reset();
  return true;
}
// ----------------------------------------------------------------------------------------
// LCE_Resize::execute
// ----------------------------------------------------------------------------------------
void LCE_Resize::execute ()
{
  
  if(_atGeneration > 0)
    if(_popPtr->getCurrentGeneration() != (unsigned)_atGeneration) return;
  
  if(++_genITER != _generations.end())
    _atGeneration = (*_genITER);
  else {
    _genITER = _generations.begin();
    _atGeneration = (*_genITER);
  }
  
  // if(_everyGeneration > 0)
  //    if(_atGeneration > 0) {
  //      //in this case, the _atGeneration value is taken as the starting point
  //      if( _popPtr->getCurrentGeneration() < (unsigned)_atGeneration ||
  //         (_popPtr->getCurrentGeneration() - _atGeneration) % _everyGeneration ) return;
  //    } else if(_popPtr->getCurrentGeneration() % _everyGeneration) return;
  
#ifdef _DEBUG_
  message("\nLCE_Resize::execute at %i (Patch nb: %i offsprg nb: %i adlt nb: %i)"
          ,_popPtr->getCurrentGeneration(),_popPtr->getPatchArraySize(),_popPtr->size( OFFSPRG ),_popPtr->size( ADULTS ));
#endif
  
  if( !(_popPtr->getCurrentAge() & _setAge) || _popPtr->getCurrentAge() == NONE ) {
    warning("LCE_Resize::execute:: required aged class not present in population, exiting.\n");
    return;
  }
  //first reset the metapop size parameters:
  ParamSet *originalPSet = _popPtr->get_paramset();
  //make a backup copy of the initial parameters state to restore them later:
  ParamSet *pSet = new ParamSet(*originalPSet);
  
  if( get_parameter("resize_patch_number")->isSet() )
    pSet->set_param("patch_number", _paramSet->getArg("resize_patch_number"));
  
  if( get_parameter("resize_keep_patch")->isSet() ) {
    if(_patch2keep.getNbCols() <= _popPtr->getPatchNbr()) {
      pSet->set_param("patch_number", tstring::int2str(_patch2keep.getNbCols()));
    } else if(_patch2keep.getNbCols() > _popPtr->getPatchNbr() && 
              _popPtr->getCurrentGeneration() == (unsigned)_atGeneration)
      //this is a problem if one want to resize AND reorder... unless temporal arg used in pop params
      fatal("LCE_Resize::execute:: more patches to keep than present in the population!\n");
  }
  
  if( get_parameter("resize_patch_capacity")->isSet() )
    pSet->set_param("patch_capacity", _paramSet->getArg("resize_patch_capacity"));
  
  if( get_parameter("resize_female_capacity")->isSet() )
    pSet->set_param("patch_nbfem", _paramSet->getArg("resize_female_capacity"));
  
  if( get_parameter("resize_male_capacity")->isSet() )
    pSet->set_param("patch_nbmal", _paramSet->getArg("resize_male_capacity"));
  
  _popPtr->set_paramset(pSet);
  //just reset the patch number and patch capacities, doesn't touch the patch structure yet
  _popPtr->setPopulationParameters();
  
  //now restore the original param set from the backup, will be needed to start the next replicate:
  _popPtr->set_paramset(originalPSet);
  
  delete pSet;
  
  if(_do_flush) {
    
    //new populations will be empty and supernumerary patches are flushed and deleted
    buildNewPatchArrayNoBackup();
    
  } else {
    
    if(_patchBackup) {
      _patchBackup->clear();
      delete _patchBackup;
    }
    _patchBackup = new Patch();
    
    _patchBackup->init(_popPtr->getNumAgeClasses(), _popPtr->getPatchKFem(), _popPtr->getPatchKMal(), 0);
    
    //should be used when existing patches are merged together (fusion) or split (fission)
    buildNewPatchArrayWithBackup();
    
  }
  
  // updatePatchCapacities();
  
  if(_do_regulate)
    if(_do_flush)  regulate( &LCE_Resize::regulateAgeClassNoBackup );
    else           regulate( &LCE_Resize::regulateAgeClassWithBackup );
  
  if(_do_fill)
    if(_do_flush) fillPop( &LCE_Resize::fillPatchNoBackup );
    else          fillPop( &LCE_Resize::fillPatchWithBackup );
  
#ifdef _DEBUG_
  message(" ---> (Patch nb: %i offsprg nb: %i adlt nb: %i, do_fill %i, do_flush %i)\n"
          ,_popPtr->getPatchNbr(),_popPtr->size( OFFSPRG ),_popPtr->size( ADULTS ), _do_fill, _do_flush);
#endif
}
// ----------------------------------------------------------------------------------------
// LCE_Resize::updatePatchCapacities
// ----------------------------------------------------------------------------------------
void LCE_Resize::updatePatchCapacities()
{
  TMatrix *cap = _popPtr->getPatchCapacities();
  
  if(cap->getNbCols() != _popPtr->getPatchNbr())
    fatal("LCE_Resize:: more patch capacities than patches in the population!\n");
  
  for (unsigned int i = 0; i < _popPtr->getPatchNbr(); i++) {
    _popPtr->getPatchPtr(i)->set_KFem((unsigned int)cap->get(FEM, i));
    _popPtr->getPatchPtr(i)->set_KMal((unsigned int)cap->get(MAL, i));
  }
}
// ----------------------------------------------------------------------------------------
// LCE_Resize::buildNewPatchArrayNoBackup
// ----------------------------------------------------------------------------------------
void LCE_Resize::buildNewPatchArrayNoBackup()
{
  if(_patch2keep.getNbCols() != 0)  removeDesignatedPatch(false);
  
  //
  //add or remove patches until the array size matches Metapop::_patchNbr
  //up to this point, the number of patches to keep could not be != from that nbre
  _popPtr->updatePatchArray();
  
}
// ----------------------------------------------------------------------------------------
// LCE_Resize::buildNewPatchArrayWithBackup
// ----------------------------------------------------------------------------------------
void LCE_Resize::buildNewPatchArrayWithBackup()
{
  unsigned int patchNbr =  _popPtr->getPatchNbr();
  
  if(_patch2keep.getNbCols() != 0)  removeDesignatedPatch(true);
  
  //reset the right number of patches to match new parameters
  if(_popPtr->getPatchArraySize() > patchNbr) {
    
    while(_popPtr->getPatchArraySize() > patchNbr) {
      _popPtr->getPatchPtr(0)->copy2patch(_patchBackup);
      _popPtr->getPatchPtr(0)->clear();
      _popPtr->removePatch(0);
    }
  } else while(_popPtr->getPatchArraySize() < patchNbr) _popPtr->addPatch(new Patch());
  //set the patch ID and Ks:
  _popPtr->updatePatchState();
  //regulate the patches, backup the supernumerary individuals, will be used in empty patches
  if(_do_fill) regulate( &LCE_Resize::regulateAgeClassWithBackup );
}
// ----------------------------------------------------------------------------------------
// LCE_Resize::removeDesignatedPatch
// ----------------------------------------------------------------------------------------
void LCE_Resize::removeDesignatedPatch(bool do_backup)
{
  deque< Patch* > to_keep;
  Patch* patch;
  unsigned int patchNbr = _popPtr->getPatchArraySize() ;
  unsigned int i, nKeep = _patch2keep.getNbCols();
  unsigned int id;
  
  if(patchNbr < nKeep) fatal("LCE_Resize::more patches to keep than available in the population\n");
  
  for(i = 0; i < nKeep; i++)
    to_keep.push_back( _popPtr->getPatchPtr((unsigned int)_patch2keep.get(0, i)-1 ) );
  
  //remove the patches we want to keep from the patch array, without touching the contents
  for(i = 0; i < nKeep; i++) {
    id = to_keep[i]->getID();
    for(unsigned int j = 0; j < _popPtr->getPatchArraySize(); ++j)
      if(_popPtr->getPatchPtr(j)->getID() == id){ _popPtr->removePatch(j); break;}
  }
  //delete the patches we don't want to keep, remove the individuals and backup them if needed
  for(i = 0; i < _popPtr->getPatchArraySize(); i++) {
    
    patch = _popPtr->getPatchPtr(i);
    
    if(do_backup) { patch->copy2patch(_patchBackup); patch->clear(); }
    
    _popPtr->deletePatch(i); //this deletes the patch and the individuals it contains
    
  }
  //rebuild the patch array with the patches we want, in the right order
  for(i = 0; i < nKeep; i++)
    _popPtr->addPatch(to_keep[i]);
  
  
}
// ----------------------------------------------------------------------------------------
// LCE_Resize::fillPop
// ----------------------------------------------------------------------------------------
void LCE_Resize::fillPop ( void (LCE_Resize:: *fillFuncPtr) (unsigned int p, age_idx age))
{
  age_t active_age = (_popPtr->getCurrentAge() & _setAge);
  
  if(active_age & OFFSPRG) { 
    for(unsigned int i=0; i < _popPtr->getPatchNbr(); ++i) 
      (this->*fillFuncPtr)(i, OFFSx); 
  }
  
  if(active_age & ADULTS) {
    for(unsigned int i=0; i < _popPtr->getPatchNbr(); ++i) 
      (this->*fillFuncPtr)(i, ADLTx); 
  }
}
// ----------------------------------------------------------------------------------------
// LCE_Resize::fillPatchNoBackup
// ----------------------------------------------------------------------------------------
void LCE_Resize::fillPatchNoBackup(unsigned int p, age_idx age)
{
  Patch *patch = _popPtr->getPatchPtr(p);
  while (patch->size(FEM, age) != patch->get_KFem()) {
    patch->add(FEM, age, _popPtr->makeNewIndividual(NULL, NULL, FEM, patch->getID()));
  }
  while (patch->size(MAL, age) != patch->get_KMal()) {
    patch->add(MAL, age, _popPtr->makeNewIndividual(NULL, NULL, MAL, patch->getID()));
  }
}
// ----------------------------------------------------------------------------------------
// LCE_Resize::fillPatchWithBackup
// ----------------------------------------------------------------------------------------
void LCE_Resize::fillPatchWithBackup(unsigned int p, age_idx age)
{
  Patch *patch = _popPtr->getPatchPtr(p);
  while (patch->size(FEM, age) != patch->get_KFem() && _patchBackup->size(FEM, age) != 0) {
    patch->add(FEM, age, _patchBackup->get(FEM, age, 0));
    _patchBackup->remove(FEM, age, 0);
  }
  while (patch->size(MAL, age) != patch->get_KMal() && _patchBackup->size(MAL, age) != 0) {
    patch->add(MAL, age, _patchBackup->get(MAL, age, 0));
    _patchBackup->remove(MAL, age, 0);
  }
}
// ----------------------------------------------------------------------------------------
// LCE_Resize::regulateWithBackup
// ----------------------------------------------------------------------------------------
void LCE_Resize::regulate( void (LCE_Resize::* regFuncPtr) (Patch *patch, age_idx age))
{
  Patch *patch;
  for(unsigned int i = 0; i < _popPtr->getPatchNbr(); i++) {
    patch = _popPtr->getPatchPtr(i);
    (this->*regFuncPtr)(patch, OFFSx);
    (this->*regFuncPtr)(patch, ADLTx);
  }
}
// ----------------------------------------------------------------------------------------
// LCE_Resize::regulateAgeClassWithBackup
// ----------------------------------------------------------------------------------------
void LCE_Resize::regulateAgeClassWithBackup(Patch *patch, age_idx age)
{
  unsigned int ind;
  
  while(patch->size(FEM, age) > patch->get_KFem()) {
    ind = RAND::Uniform( patch->size(FEM, age) ) ;
    _patchBackup->add(FEM, age, patch->get(FEM, age, ind));
    patch->remove(FEM, age, ind);
  }
  
  while(patch->size(MAL, age) > patch->get_KMal()) {
    ind = RAND::Uniform( patch->size(MAL, age) ) ;
    _patchBackup->add(MAL, age, patch->get(MAL, age, ind));
    patch->remove(MAL, age, ind);
  }
  
}
// ----------------------------------------------------------------------------------------
// LCE_Resize::regulateAgeClassNoBackup
// ----------------------------------------------------------------------------------------
void LCE_Resize::regulateAgeClassNoBackup(Patch *patch, age_idx age)
{
  unsigned int ind;
  
  while(patch->size(FEM, age) > patch->get_KFem()) {
    ind = RAND::Uniform( patch->size(FEM, age) ) ;
    _popPtr->recycle( patch->get(FEM, age, ind) );
    patch->remove(FEM, age, ind);
  }
  
  while(patch->size(MAL, age) > patch->get_KMal()) {
    ind = RAND::Uniform( patch->size(MAL, age) ) ;
    _popPtr->recycle( patch->get(MAL, age, ind) );
    patch->remove(MAL, age, ind);
  }
  
}

