/** $Id: LCEdisperse.h,v 1.6.2.3 2015-04-13 15:11:55 fred Exp $
 *
 *  @file LCEdisperse.h
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
 *  Created on @date 07.08.2004
 *  @author fred
 */

#ifndef LCEDISPERSE_H
#define LCEDISPERSE_H
#include <vector>
#include "lifecycleevent.h"
#include "param.h"

class Patch;
/**The base class of the dispersal LCEs, all events move offspring to the post-dispersal patch containers.
 * Stores the dispersal matrices and dispersal model parameters and interface.
 */
class LCE_Disperse_base: public virtual LifeCycleEvent
{
  
  int    _disp_model;
  double _disp_propagule_prob;
  vector<unsigned int> _PropaguleTargets;
  double _fem_rate, _mal_rate;
  bool _isForward;
  
  /**The sex-specific dispersal matrices, [0] for males, [1] for females, might be used as connectivity matrix as well*/
  TMatrix* _DispMatrix[2]; 
  TMatrix* _patchConnected;
  TMatrix* _pDisp;

  string _prefix;
  
  friend class LCE_Disperse_ConstDisp;
  friend class LCE_Disperse_EvolDisp;
  
protected:
  
  unsigned int _npatch;

  vector<Patch*> _postDispPop;
  vector< vector<double> > _reducedDispMat[2];
  vector<vector<double> > _reducedDispMatProba[2];
  
public:
  
  LCE_Disperse_base();
  
  /**Deallocates the disp matrix.*/
  virtual ~LCE_Disperse_base(); 
  
  bool setBaseParameters(string prefix);
  void setParamPrefix (string pref) {_prefix = pref;}
  void addParameters  (string prefix, ParamUpdaterBase* updater);
  void setPostDispPop ();
  void resetPostDispPop ();
  
  ///@name Dispersal Matrix
  ///@{
  void set_isForward(bool val) {_isForward = val;}
  bool checkForwardDispersalMatrix    (TMatrix* mat);
  bool checkBackwardDispersalMatrix   (TMatrix* mat);
  void allocateDispMatrix  (sex_t sex, unsigned int dim);
  bool setDispMatrix();
  bool updateDispMatrix();
  bool setReducedDispMatrix();
  bool setIsland_MigrantPool_Matrix();
  bool setIsland_PropagulePool_Matrix();
  bool setSteppingStone1DMatrix();
  bool setLatticeMatrix();
  bool setBasicLatticeMatrix(unsigned int side, double phi_mal, double phi_fem, double disp_mal, double disp_fem);
  bool setLatticeTorrusMatrix(unsigned int side, double disp_mal, double disp_fem);
  bool setLatticeAbsorbingMatrix( );
  bool setLatticeReflectingMatrix(unsigned int side);
  ///@}  
  unsigned int getMigrationPatchForward  (sex_t SEX, unsigned int LocalPatch);
  unsigned int getMigrationPatchBackward (sex_t SEX, unsigned int LocalPatch);
  void setPropaguleTargets ( );
  void swapPostDisp ( );
  void reset_counters ( );
  ///@name Accessors
  ///@{
  unsigned int getDispersalModel  ( ) {return _disp_model;}
  double       getPropaguleProb   ( ) {return _disp_propagule_prob;}
  unsigned int getPropaguleTarget (unsigned int home) {return _PropaguleTargets[home];}
  ///@}
  ///@name Implementations
  ///@{
  virtual void loadFileServices ( FileServices* loader ) {}
  virtual void loadStatServices ( StatServices* loader ) {}
  virtual age_t removeAgeClass ( ) {return NONE;}
  virtual age_t addAgeClass ( ) {return NONE;}
  virtual age_t requiredAgeClass () {return OFFSPRG;}
  ///@}
};

/**Dispersal event with constant dispersal rates.
 * Sets and uses the dispersal matrices according to the dispersal model chosen. 
 * Dispersal models implemented so far are:
 * <ul> 
 * <li>1: Island Model with migrant pool
 * <li>2: Island Model with propagule pool
 * <li>3: Stepping stone model in 1 dimension (ring population)
 * <li>4: Lattice model (stepping stone in 2D)
 * </ul>
 **/
class LCE_Disperse_ConstDisp: public virtual LCE_Disperse_base 
{
  void (LCE_Disperse_ConstDisp::* doMigration) (void);
  void (LCE_Disperse_ConstDisp::* doPatchMigration) (sex_t SEX, unsigned int local_patch);
  
public:
	
  LCE_Disperse_ConstDisp (); // : LifeCycleEvent ("disperse","") {}  
  virtual ~LCE_Disperse_ConstDisp(){}

  bool setParameters (string prefix);
  
  void Migrate ( );
  void Migrate_propagule ( );
  void MigratePatch (sex_t SEX, unsigned int LocalPatch);
  void MigratePatch_AbsorbingBorder (sex_t SEX, unsigned int LocalPatch);
  ///@name Implementations
  ///@{
  virtual bool setParameters () {return setParameters("dispersal");}
  virtual void execute ();
  virtual LifeCycleEvent* clone () {return new LCE_Disperse_ConstDisp();}
  ///@}
};



class LCE_SeedDisp : public LCE_Disperse_ConstDisp {
  

public:
  LCE_SeedDisp ();
  virtual ~LCE_SeedDisp(){}
  virtual bool setParameters () {return LCE_Disperse_ConstDisp::setParameters("seed_disp");}
  virtual LifeCycleEvent* clone () {return new LCE_SeedDisp();}
};



/**Dispersal event with an evolving dispersal rate given by the "disp" trait.
 * The dispersal models implemented so far are:
 * <ul> 
 * <li>1: Island Model with migrant pool
 * <li>2: Island Model with propagule pool
 * <li>3: Stepping stone model in 1 dimension (ring population)
 * </ul>
 * 
**/
class LCE_Disperse_EvolDisp: public virtual LCE_Disperse_base
{
  
  double _fem_cost, _mal_cost, _fixed_disp_rate;
  unsigned int (LCE_Disperse_EvolDisp::* getAimedPatch) (unsigned int);
  void (LCE_Disperse_EvolDisp::* exec) ();
  
  unsigned int Migrate_Island (unsigned int home);
  unsigned int Migrate_Island_Propagule (unsigned int home);
  unsigned int Migrate_SteppingStone1D (unsigned int home);
  void exec_evolmale ();
  void exec_evolfemale ();
  void exec_evol2sex ();
  void evoldisp (sex_t SEX, int trait_link, double cost);
  void fixdisp (sex_t SEX, double rate, double cost);
  
  int _fdisp_trait_link, _mdisp_trait_link;
  
//  typedef LifeCycleEvent< LCE_Disperse_EvolDisp > LCEDispEvol;
public:
	
    LCE_Disperse_EvolDisp();
  
  virtual ~LCE_Disperse_EvolDisp(){}
    
  ///@name Implementations
  ///@{
  virtual bool setParameters ();
  virtual void execute ();
  virtual LifeCycleEvent* clone () {return new LCE_Disperse_EvolDisp();}
  ///@}
};

#endif //LCEDISPERSE_H
