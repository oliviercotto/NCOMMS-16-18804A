/**  $Id: LCEmisc.h,v 1.7.2.6 2016-04-28 13:05:05 fred Exp $
*
*  @file LCEmisc.h
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
*  created on @date 06.16.2005
* 
*  @author fred
*/

#ifndef LCEMISC_H
#define LCEMISC_H

#include <iostream>
#include <string>
#include <map>
#include "types.h"
#include "lifecycleevent.h"
#include "fileservices.h"
#include "statservices.h"
#include "filehandler.h"
#include "stathandler.h"
#include "Uniform.h"

// LCE_Regulation
//
/** Regulates the patches to their carrying capacity, acts on each age class separately.
 *  Randomly removes individuals in each patch until each age class has reached the 
 *  (sex-specific) carrying capacities of the patch.
**/
class LCE_Regulation: public LifeCycleEvent
{  
private:

	double _competition_regulation;
	unsigned int _K;
//	TMatrix* _survival_regulation;
	unsigned int _nb_class;
	age_idx _age_target;
	age_t _age_flag;
  
  // pointer to function used for regulation
  vector<void (LCE_Regulation::*) (Patch*, sex_t, unsigned int)> _regulate_ptrs;

public:
  
  LCE_Regulation( ) : LifeCycleEvent("regulation",""), _competition_regulation(0), _K(0),
  	  	  	  	  	  _nb_class(0), _age_target(OFFSx), _age_flag(ADULTS)
  {
	  add_parameter("regulation_by_competition",DBL,false,false,0,0);
	  add_parameter("regulation_by_competition_affected_age",INT,false,false,0,0);
	  add_parameter("regulation_by_competition_count_age_flag",INT,false,true,0,32);
	  add_parameter("regulation_carrying_capacity",BOOL,false,false,0,0);
  }
  
  virtual ~LCE_Regulation( ) {}
  
  // regulate by competition:
  void regulatePatch (Patch* patch, sex_t sex, unsigned int pop_size);
  
  // ceiling regulation to carrying capacity:
  void regulatePatchCeiling (Patch* patch, sex_t sex, unsigned int pop_size);
  
  //implementations:
  virtual bool setParameters ();// {return true;}
  virtual void execute ();
  virtual LCE_Regulation* clone ( ) {return new LCE_Regulation();}
  virtual void loadFileServices ( FileServices* loader ) {}
  virtual void loadStatServices ( StatServices* loader ) {}
  virtual age_t removeAgeClass ( ) {return NONE;}
  virtual age_t addAgeClass ( ) {return NONE;}
  virtual age_t requiredAgeClass () {return NONE;}
};

// LCE_Aging
//
/**Removes all adults from the patches and randomly moves the offspring to the adults age class.
 * Patch regulation is performed at the same time, leaving the patch at carrying capacity (if enough
 * offspring individuals were present).
 * This is the only LCE that actually removes the adults from the patches.
 * Also checks whether the patch is filled and sets the extinction flag accordingly.*/
class LCE_Aging: public LifeCycleEvent
{
public:
  
  LCE_Aging( ) : LifeCycleEvent("aging","") {}
  
  virtual ~LCE_Aging( ) { }
  
  //implementations:
  virtual bool setParameters () {return true;}
  virtual void execute ();
  virtual LCE_Aging* clone ( ) {return new LCE_Aging();}
  virtual void loadFileServices ( FileServices* loader ) {}
  virtual void loadStatServices ( StatServices* loader ) {}
  virtual age_t removeAgeClass ( ) {return OFFSPRG;}
  virtual age_t addAgeClass ( ) {return ADULTS;}
  virtual age_t requiredAgeClass () {return OFFSPRG;}
};

// LCE_Aging_multi
//
/**Performs aging with overlapping generations.
 * Uses the Leslie matrix to decide on individual survival during aging.
 * Some parameters are specific to an older, now unused model.
 * Will need clean-up!!*/
class LCE_Aging_Multi: public LifeCycleEvent
{
  
  TMatrix* _survival;
  
public:
  
  LCE_Aging_Multi( ) : LifeCycleEvent("aging_multi",""),
  _survival(0)
  {add_parameter("regulation_before_aging",BOOL,false,false,0,0,0);}
  
  virtual ~LCE_Aging_Multi( ) {}
  
  
  //implementations:
  virtual bool setParameters ();
  virtual void execute ();
  
  virtual LCE_Aging_Multi* clone ( ) {return new LCE_Aging_Multi();}
  
  virtual void loadFileServices ( FileServices* loader ) {}
  virtual void loadStatServices ( StatServices* loader ) {}
  virtual age_t removeAgeClass ( ) {return OFFSPRG;}
  virtual age_t addAgeClass ( ) {return ADULTS;}
  virtual age_t requiredAgeClass () {return ADULTS | OFFSPRG;}
};

// LCE_Patch_Extinction
//
/**Randomly removes individuals from the patches according to the extinction rate parameter.
 * Sets the patches extinction flag accordingly.**/
class LCE_Patch_Extinction: public LifeCycleEvent
{
  /**Patch extinction probability.*/
  TMatrix* _Xtion_rate;
  /**Number of individual to remove per patch.*/
  TMatrix* _harvest_size;
  /**Proportion of the patch size to remove.*/
  TMatrix* _harvest_proportion;
  /**Name of the distribution to use*/
  string _harvest_distribution;
  /**Flags*/
  bool _harvest_size_varies, _by_size, _by_proportion;
  /**Standard deviate to use with the Gaussian dist.*/
  double _harvest_dist_stdev;
  /**shape variable to use with the gamma dist.*/
  double _harvest_dist_shape;
  /**Patch extinction threshold in % of total size of the patch*/
  double _extinction_threshold;
  
  unsigned int (LCE_Patch_Extinction:: *_rand_size_fct) (double);
  
  unsigned int rand_uniform   (double max) {return RAND::Uniform((unsigned int)max);}
  unsigned int rand_poisson   (double mean){return (unsigned int)RAND::Poisson(mean);}
  unsigned int rand_gaussian  (double mean){return (unsigned int)abs(mean + RAND::Gaussian(_harvest_dist_stdev));}
  unsigned int rand_exp       (double mean){return (unsigned int)(-1.0 * mean * log(1.0-RAND::Uniform()));}
  unsigned int rand_lognormal (double mean){return (unsigned int)RAND::LogNormal(mean, _harvest_dist_stdev);}
  
public:
	
    LCE_Patch_Extinction( ) ;  
  virtual ~LCE_Patch_Extinction( ) 
  {
    if(_Xtion_rate) delete _Xtion_rate;
    if(_harvest_size) delete _harvest_size; 
    if(_harvest_proportion) delete _harvest_proportion;
  }
  
  bool set_matrix_param (TMatrix* mat, string name);
  void do_flush (Patch *patch);
  void do_remove (age_idx AGE, Patch *patch);
  unsigned int get_harvest_size (age_idx AGE, Patch *patch);
  
  //LifeCycleEvent interface:
//  virtual void init (Metapop* popPtr);
  virtual bool setParameters ();
  virtual void execute ();
  virtual LifeCycleEvent* clone ( ) {return new LCE_Patch_Extinction();}

  //SimComponent interface:
  virtual void loadFileServices ( FileServices* loader ) {}
  virtual void loadStatServices ( StatServices* loader ) {}
  virtual age_t removeAgeClass ( ) {return 0;}
  virtual age_t addAgeClass ( ) {return 0;}
  virtual age_t requiredAgeClass () {return 0;}
};

//CROSSING    
//
/**A class to perform crosses within patches, implements the NC1 mating design.*/
class LCE_Cross : public LifeCycleEvent
{
  unsigned int _nSire, _nDam;
  unsigned int _nOffspring;
  unsigned int _atGeneration;
  bool _doAmongPop, _doWithinPop, _doReplace;
  
  
public:
    
    LCE_Cross( );
  virtual ~LCE_Cross( ) { }
    
  void sampleAmongPop(Patch* patch, deque<Individual*>& males, unsigned int nsire);
  void sampleWithinPop(Patch* patch, deque<Individual*>& males, deque<Individual*>& females, unsigned int nsire);
  
  //LifeCycleEvent interface:
//  virtual void init (Metapop* poapPtr);
  virtual bool setParameters ();
  virtual void execute ();
  virtual LifeCycleEvent* clone ( ) {return new LCE_Cross();}
  
  //SimComponent interface:
  virtual void loadFileServices ( FileServices* loader ) {}
  virtual void loadStatServices ( StatServices* loader ) {}
  virtual age_t removeAgeClass ( ) {return 0;}
  virtual age_t addAgeClass ( ) {return OFFSPRG;}
  virtual age_t requiredAgeClass () {return ADULTS;}
  
};

//RESIZE POP    
//
/**A class to change the size of the population/patches during a simulation.*/
class LCE_Resize : public LifeCycleEvent
  {
    list< int >::const_iterator _genITER;
    int _atGeneration;
    list< int > _generations;
    bool _do_flush, _do_fill, _do_regulate;
    TMatrix _patch2keep;
    age_t _setAge;
    Patch* _patchBackup;
  public:
    
    LCE_Resize( );
    virtual ~LCE_Resize( ) { }
    
    void buildNewPatchArrayNoBackup();
    void buildNewPatchArrayWithBackup();
    void removeDesignatedPatch(bool do_backup);
    void updatePatchCapacities();
    void fillPop ( void (LCE_Resize:: *fillFuncPtr) (unsigned int p, age_idx age));
    void fillPatchNoBackup(unsigned int p, age_idx age);
    void fillPatchWithBackup(unsigned int p, age_idx age);
    void regulate( void (LCE_Resize::* regFuncPtr) (Patch *patch, age_idx age));
    void regulateAgeClassWithBackup(Patch *patch, age_idx age);
    void regulateAgeClassNoBackup(Patch *patch, age_idx age);
    
    //LifeCycleEvent interface:
    virtual bool setParameters ();
    bool         updateParameters ();
    virtual void execute ();
    virtual LifeCycleEvent* clone ( ) {return new LCE_Resize();}
    
    //SimComponent interface:
    virtual void loadFileServices ( FileServices* loader ) {}
    virtual void loadStatServices ( StatServices* loader ) {}
    virtual age_t removeAgeClass ( ) {return 0;}
    virtual age_t addAgeClass ( ) {return 0;}
    virtual age_t requiredAgeClass () {return 0;}
    
  };

#endif //LCEMISC_H
