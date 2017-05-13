/** $Id: LCEdisperse.cc,v 1.14.2.6 2015-08-06 06:46:06 fred Exp $
 *
 *  @file LCEdisperse.cc
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
 *  Created on @date 08.07.2004
 *  @author fred
 */

#include <stdlib.h>
#include <string>
#include <cmath>
#include "LCEdisperse.h"
#include "metapop.h"
#include "individual.h"
#include "Uniform.h"

using namespace std;
/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

//                             ******** LCE_Disperse_base ********/

// ----------------------------------------------------------------------------------------
LCE_Disperse_base::LCE_Disperse_base () 
: LifeCycleEvent("", ""), _disp_model(-1), _disp_propagule_prob(-1.0), _PropaguleTargets()
, _fem_rate (-1), _mal_rate(-1), _isForward(1), _patchConnected(0), _pDisp(0)
{
  _DispMatrix[0] = NULL;
  _DispMatrix[1] = NULL;
//  cout << "calling LCE_Disperse_base()\n";
}
// ----------------------------------------------------------------------------------------
LCE_Disperse_base::~LCE_Disperse_base ()
{
  if(NULL != _DispMatrix[0])
    delete _DispMatrix[0];
  
  if(NULL != _DispMatrix[1])
    delete _DispMatrix[1];
  
  for (unsigned int i = 0; i < _postDispPop.size(); ++i) {
    if(_postDispPop[i]) delete _postDispPop[i]; //!! careful when deleting individuals here...
  }
  _postDispPop.clear();
}
// ----------------------------------------------------------------------------------------
void LCE_Disperse_base::addParameters (string prefix, ParamUpdaterBase* updater)
{
  add_parameter(prefix + "_model",INT,false,true,1,4,updater);
  add_parameter(prefix + "_border_model",INT,false,true,1,3,updater);
  add_parameter(prefix + "_lattice_range",INT,false,true,1,2,updater);
  add_parameter(prefix + "_propagule_prob",DBL,false,true,0,1,updater);
  add_parameter(prefix + "_matrix",MAT,false,false,0,0,updater);
  add_parameter(prefix + "_matrix_fem",STR,false,false,0,0,updater);
  add_parameter(prefix + "_matrix_mal",STR,false,false,0,0,updater);
  add_parameter(prefix + "_rate",DBL,false,true,0,1,updater);
  add_parameter(prefix + "_rate_fem",DBL,false,true,0,1,updater);
  add_parameter(prefix + "_rate_mal",DBL,false,true,0,1,updater);
  add_parameter(prefix + "_connectivity_matrix",MAT,false,false,0,0,updater);
  add_parameter(prefix + "_reduced_matrix",MAT,false,false,0,0,updater);
  
  //for backward compatibility with previous Olivier's version
  add_parameter(prefix + "_patchConnected",MAT,false,false,0,0,updater);
  add_parameter(prefix + "_pDisp",MAT,false,false,0,0,updater);
}
// ----------------------------------------------------------------------------------------
// LCE_Disperse_base::setParameters
// ----------------------------------------------------------------------------------------
bool LCE_Disperse_base::setBaseParameters(string prefix)
{
  
  //cout<<"set base parametre"<<endl;
  
  _prefix = prefix;
  
  _npatch = _popPtr->getPatchNbr();
  
  setPostDispPop();
  
  _disp_model = (int)_paramSet->getValue(prefix + "_model");
  
  _disp_propagule_prob = _paramSet->getValue(prefix + "_propagule_prob");
  

  
  if(_paramSet->isSet(prefix + "_connectivity_matrix") && 
     _paramSet->isSet(prefix + "_reduced_matrix")){
    
    _disp_model=0;
        
    get_parameter(prefix + "_connectivity_matrix")->getVariableMatrix(&_reducedDispMat[0]);
    
    get_parameter(prefix + "_reduced_matrix")->getVariableMatrix(&_reducedDispMatProba[0]);
    
    //need to copy matrices for the other sex (a single matrix is given in input)
    
    _reducedDispMat[1].clear();
    _reducedDispMatProba[1].clear();
    
    if (_reducedDispMat[0].size() != _reducedDispMatProba[0].size()) {
      error("The connectivity and reduced dispersal matrices don't have same number of rows\n");
      return false;
    }
    
    for (unsigned int i = 0; i < _reducedDispMat[0].size(); ++i) {
      
      if (_reducedDispMat[0][i].size() != _reducedDispMatProba[0][i].size()) {
        error("Row %i of the connectivity and reduced dispersal matrices are not of same size\n", i+1);
        return false;
      }
      
      _reducedDispMat[1].push_back(vector<double>());
      
      _reducedDispMatProba[1].push_back(vector<double>());
      
      double row_sum = 0;
      //@TODO we should here deduce 1 from the patch ID in the connectivity matrix
      // this is because we assume IDs be given in range [1 - num patch] in input
      // this forces us to deduce 1 when we read IDs during migration, silly!
      for (unsigned int j = 0; j < _reducedDispMat[0][i].size(); ++j) {
        _reducedDispMat[1][i].push_back( _reducedDispMat[0][i][j] );
        _reducedDispMatProba[1][i].push_back( _reducedDispMatProba[0][i][j] );
        row_sum += _reducedDispMatProba[0][i][j];
      }
      
      if(row_sum < 0.999999 || row_sum > 1.000001) {
        error("the elements of row %i of the reduced dispersal matrix do not sum to 1!\n",i+1);
        return false;
      }
    }
    
    
  }
  
  else if ( (_paramSet->isSet(prefix + "_connectivity_matrix") && 
             !_paramSet->isSet(prefix + "_reduced_matrix")) || 
            (!_paramSet->isSet(prefix + "_connectivity_matrix") && 
             _paramSet->isSet(prefix + "_reduced_matrix")))
    
  {
    error("both \"%s_connectivity_matrix\" and \"%s_reduced_matrix\" must be set together\n", prefix.c_str(), prefix.c_str());
    return false;
  }
  
  else
  
  {
    
    if(_paramSet->isSet(prefix + "_matrix")) {
      
      if(_DispMatrix[0] == NULL) _DispMatrix[0] = new TMatrix();
      
      _paramSet->getMatrix(prefix + "_matrix",_DispMatrix[0]);
      
      if(_isForward) {
        if(!checkForwardDispersalMatrix(_DispMatrix[0])) return false;
      } else { 
        if(!checkBackwardDispersalMatrix(_DispMatrix[0])) return false;
      }
      
      if(_DispMatrix[1] != NULL) delete _DispMatrix[1];
      
      _DispMatrix[1] = new TMatrix(*_DispMatrix[0]);
      
    } else {
      if(_paramSet->isSet(prefix + "_matrix_fem")) {
        
        if(_DispMatrix[FEM] == NULL) _DispMatrix[FEM] = new TMatrix();
        
        _paramSet->getMatrix(prefix + "_matrix_fem",_DispMatrix[FEM]);
        
        if(_isForward) {
          if(!checkForwardDispersalMatrix(_DispMatrix[FEM])) return false;
        } else { 
          if(!checkBackwardDispersalMatrix(_DispMatrix[FEM])) return false;
        }
      }
      
      if(_paramSet->isSet(prefix + "_matrix_mal")) {
        
        if(_DispMatrix[MAL] == NULL) _DispMatrix[MAL] = new TMatrix();
        
        _paramSet->getMatrix(prefix + "_matrix_mal",_DispMatrix[MAL]);
        
        if(_isForward) {
          if(!checkForwardDispersalMatrix(_DispMatrix[MAL])) return false;
        } else {
          if(!checkBackwardDispersalMatrix(_DispMatrix[MAL])) return false;
        }
      }
    }
    
    if( _paramSet->isSet(prefix + "_matrix") || 
       ( _paramSet->isSet(prefix + "_matrix_fem") && _paramSet->isSet(prefix + "_matrix_mal") )  )
    {
      if(  ( _paramSet->isSet(prefix + "_rate") ||
            (_paramSet->isSet(prefix + "_rate_fem") &&  _paramSet->isSet(prefix + "_rate_mal")) )
         || _paramSet->isSet(prefix + "_model") )
        warning("parameter \"dispersal_matrix\" takes precedence over parameters \"dispersal_rate\" and \"dispersal_model\"\n");
      
      _disp_model = 0;
      
    } else {
      
      if(!_paramSet->isSet(prefix + "_model")) {
        error("Dispersal model not set!\n");
        return false;
      }
      
      if(_paramSet->isSet(prefix + "_rate")) 
        
      {
        _fem_rate = _mal_rate = _paramSet->getValue(prefix + "_rate");
        
        if(!setDispMatrix()) return false;
      }
      
      else if(  _paramSet->isSet(prefix + "_rate_fem") &&  _paramSet->isSet(prefix + "_rate_mal")  ) 
        
      {
        _fem_rate = _paramSet->getValue(prefix + "_rate_fem");
        
        _mal_rate = _paramSet->getValue(prefix + "_rate_mal");
        
        if(!setDispMatrix()) return false;
      }
      
      else {
        error("Dispersal rate parameters not set!\n");
        return false;
      }
      
    }  
    
  }
  return true;
}
// ----------------------------------------------------------------------------------------
// LCE_Disperse_base::updateDispMatrix
// ----------------------------------------------------------------------------------------
bool LCE_Disperse_base::updateDispMatrix()
{

//cout<<"check1"<<endl;

  if ( getDispersalModel() == 0) {
    error("cannot update the dispersal matrix provided in input when number of populations changes.\n");
    return false;
  }
  
  return setDispMatrix();
}
// ----------------------------------------------------------------------------------------
// LCE_Disperse_base::setPostDispPop
// ----------------------------------------------------------------------------------------
void LCE_Disperse_base::setPostDispPop ( )
{

  Patch* new_patch;
  
  if (_postDispPop.size() != 0) {
    resetPostDispPop();
  }
  
  for (unsigned int i = 0; i < _npatch; ++i) {
    new_patch = new Patch();
    new_patch->init(1, _popPtr->getPatchKFem(), _popPtr->getPatchKMal(), i);
    _postDispPop.push_back(new_patch);
  }
  
}
// ----------------------------------------------------------------------------------------
// LCE_Disperse_base::setPostDispPop
// ----------------------------------------------------------------------------------------
void LCE_Disperse_base::resetPostDispPop ( )
{




  for (unsigned int i = 0; i < _postDispPop.size(); ++i) {
    if(_postDispPop[i]->size() != 0) _postDispPop[i]->clear(); //active pointers are lost.
    delete _postDispPop[i];//clear must be called before delete otherwise individuals are lost
  }
  _postDispPop.clear();
}
// ----------------------------------------------------------------------------------------
// LCE_Disperse_base::swapPostDisp
// ----------------------------------------------------------------------------------------
void LCE_Disperse_base::swapPostDisp ( )
{




  Patch *patch;
  
  for(unsigned int i = 0; i < _npatch; i++) {
    patch = _popPtr->getPatch(i);
    patch->slice(FEM, OFFSx, (*_postDispPop[i]), 0, _postDispPop[i]->size(FEM, OFFSx));
    patch->slice(MAL, OFFSx, (*_postDispPop[i]), 0, _postDispPop[i]->size(MAL, OFFSx));
    _postDispPop[i]->clear(FEM, OFFSx);
    _postDispPop[i]->clear(MAL, OFFSx);
  }
}
// ----------------------------------------------------------------------------------------
// LCE_Disperse_base::reset_counters
// ----------------------------------------------------------------------------------------
void LCE_Disperse_base::reset_counters()
{
  
  
  Patch *patch;
  for(unsigned int i = 0; i < _npatch; i++) {
    
    patch = _popPtr->getPatch(i);
    
    patch->reset_counters();
    
//    patch->flush(PDISPx, _popPtr);
  }  
}
// ----------------------------------------------------------------------------------------
// LCE_Disperse_base::allocateDispMatrix
// ----------------------------------------------------------------------------------------
void LCE_Disperse_base::allocateDispMatrix (sex_t sex, unsigned int dim)
{  
  

  
  if(_DispMatrix[sex] != NULL)
    _DispMatrix[sex]->reset(dim,dim);

  else
    _DispMatrix[sex] = new TMatrix(dim,dim);

}
// ----------------------------------------------------------------------------------------
// LCE_Disperse_base::checkDispMatrix
// ----------------------------------------------------------------------------------------
bool LCE_Disperse_base::checkForwardDispersalMatrix (TMatrix* mat)
{
  unsigned int dim = _popPtr->getPatchNbr();
  double cntr;
  
  if(mat->length() != dim*dim) {
    error("the size of the dispersal matrix is not equal to patch_number X patch_number (%i[%i,%i] != %i)!\n",
          mat->length(),mat->getNbRows(),mat->getNbCols(),dim*dim);
    return false;
  }
  
  for(unsigned int i = 0; i < dim; ++i) {
    cntr = 0;
    for(unsigned int j = 0; j < dim; ++j)
      cntr += mat->get(i,j);
    if(cntr < 0.999999 || cntr > 1.000001) {
      error("the elements of row %i of the dispersal matrix do not sum to 1!\n",i+1);
      return false;
    }
  }
  return true;
}
// ----------------------------------------------------------------------------------------
// LCE_Disperse_base::checkDispMatrix
// ----------------------------------------------------------------------------------------
bool LCE_Disperse_base::checkBackwardDispersalMatrix (TMatrix* mat)
{
  //first check that the matrix is properly set in the forward way:
  //checkForwardDispersalMatrix(mat);
  //normalize the columns:
  unsigned int dim = _popPtr->getPatchNbr();
  double cntr;
  
  if(mat->length() != dim*dim) {
    error("The size of the dispersal matrix is not equal to patch_number X patch_number (%i[%i,%i] != %i)!\n",
          mat->length(),mat->getNbRows(),mat->getNbCols(),dim*dim);
    return false;
  }
  
  for(unsigned int i = 0; i < dim; ++i) {
    cntr = 0;
    for(unsigned int j = 0; j < dim; ++j)
      cntr += mat->get(j,i);
    if(cntr < 0.999999 || cntr > 1.000001) {
      error("The elements of column %i of the dispersal matrix do not sum to 1!\n",i+1);
      return false;
    }
  }
  return true;
}
// ----------------------------------------------------------------------------------------
// LCE_Disperse_base::setPropaguleTargets
// ----------------------------------------------------------------------------------------
void LCE_Disperse_base::setPropaguleTargets ( )
{
  unsigned int nb_patch = _popPtr->getPatchNbr();
  unsigned int tmp_array[nb_patch];
  unsigned int table_index, target_patch;
  unsigned int size, last;
  
  //shuffling algorithm:
  do {
    for(unsigned int i = 0; i < nb_patch; ++i)
      tmp_array[i] = i;
    
    size = nb_patch;
    
    for(unsigned int orig_patch = 0; orig_patch < nb_patch-1; ++orig_patch) {
      
      do{
        
        table_index = RAND::Uniform( size );
        
        target_patch = tmp_array[ table_index ];
        
      }while(target_patch == orig_patch);
      
      size--;
      
      last = tmp_array[size];
      
      tmp_array[table_index] = last;
      
      tmp_array[size] = target_patch;
    }
    //do this until the last element left is not the last patch:
  }while (tmp_array[0] == nb_patch-1);
  
  _PropaguleTargets.assign(nb_patch,0);
  
  unsigned int reverse_i = nb_patch;
  
  //we read the shuffled array in reverse order:
  for(unsigned int i=0; i < _PropaguleTargets.size(); i++) {
    _PropaguleTargets[i] = tmp_array[--reverse_i];
    
#ifdef _DEBUG_
    cout<<" -- Patch "<<i<<" : assigned Patch "<<_PropaguleTargets[i]<<endl;
#endif
    
  }
}
// ----------------------------------------------------------------------------------------
// LCE_Disperse_base::setDispMatrix
// ----------------------------------------------------------------------------------------
bool LCE_Disperse_base::setDispMatrix ()
{  

  switch ( getDispersalModel() ) {
    case 1:
      if(!setIsland_MigrantPool_Matrix()) return false;
      break;
    case 2:
      if( !setIsland_PropagulePool_Matrix()) return false;
      break;
    case 3:
      if (!setSteppingStone1DMatrix()) return false;
      break;
    case 4:
      if(!setLatticeMatrix()) return false;
      break;
    default: {
      error("\nDispersal model '%i' not yet implemented\n",getDispersalModel());
      return false;
    }
  }

  //cout<<"test2"<<endl;

	

  return setReducedDispMatrix();
}
// ----------------------------------------------------------------------------------------
// LCE_Disperse_base::setIsland_MigrantPool_Matrix()  (set the Island dispersal matrix)
// ----------------------------------------------------------------------------------------
bool LCE_Disperse_base::setIsland_MigrantPool_Matrix()
{
#ifdef _DEBUG_
  cout<<"setIsland_MigrantPool_Matrix(_npatch="<<_npatch<<", _mal_rate="
      <<_mal_rate<<", _fem_rate="<<_fem_rate<<")"<<endl;
#endif
  allocateDispMatrix(MAL, _npatch);
  allocateDispMatrix(FEM, _npatch);
  
  TMatrix* mmat = _DispMatrix[MAL];
  TMatrix* fmat = _DispMatrix[FEM];
  double pmal = 1 - _mal_rate;
  double pfem = 1 - _fem_rate;
  double mmal = _mal_rate/(_npatch-1);
  double mfem = _fem_rate/(_npatch-1);
  
  for (unsigned int i=0; i<_npatch; ++i){
    for (unsigned int j=0; j<_npatch; ++j){
      mmat->set(i,j, mmal);
      fmat->set(i,j, mfem);
    }
  }
  
  for (unsigned int i=0; i<_npatch; ++i){
    mmat->set(i,i, pmal);
    fmat->set(i,i, pfem);
  }
  return true;
}
// ----------------------------------------------------------------------------------------
// LCE_Disperse_base::setIsland_PropagulePool_Matrix()  (set the Island dispersal matrix)
// ----------------------------------------------------------------------------------------
bool LCE_Disperse_base::setIsland_PropagulePool_Matrix()
{
  allocateDispMatrix(MAL, _npatch);
  allocateDispMatrix(FEM, _npatch);
  
  if( !_paramSet->isSet(_prefix + "_propagule_prob") ) {
    error("Missing parameter \"dispersal_propagule_prob\" with dispersal model 2!\n");
    return false;
  }
  
  setPropaguleTargets();
  
  double propagulePHI = getPropaguleProb();
  double c1 = (1 - _fem_rate), c2 = (_fem_rate*propagulePHI),
	c3 = (_fem_rate*(1.0 - propagulePHI)/(_npatch-2));
  
  TMatrix* mmat = _DispMatrix[MAL];
  TMatrix* fmat = _DispMatrix[FEM];
  
  for (unsigned int i=0; i < _npatch; ++i){
		
    fmat->set(i, i, c1);
    
    for (unsigned int j=i+1; j < _npatch; ++j){
      fmat->set(i, j, c3);
      fmat->set(j, i, c3);
    }
    fmat->set(i, getPropaguleTarget(i), c2);
  }
  
  c1 = (1 - _mal_rate);
  c2 = (_mal_rate*propagulePHI);
  c3 = (_mal_rate*(1.0 - propagulePHI)/(_npatch-2));
  
  for (unsigned int i=0; i < _npatch; ++i){
    
    mmat->set(i, i, c1);
    
    for (unsigned int j=i+1; j< _npatch; ++j) {
      mmat->set(i, j, c3);
      mmat->set(j, i, c3);
    }
    mmat->set(i, getPropaguleTarget(i), c2);
  }
  
  return true;
}
// ----------------------------------------------------------------------------------------
// LCE_Disperse_base::setSteppingStone1DMatrix
// ----------------------------------------------------------------------------------------
bool LCE_Disperse_base::setSteppingStone1DMatrix()
{
#ifdef _DEBUG_
  message("setSteppingStone1DMatrix()\n");
#endif
  int border_model = (unsigned int)get_parameter_value(_prefix + "_border_model");
  
  //check for the border model, the extra patch is the absorbing patch
  if(border_model == 3) {
    allocateDispMatrix(MAL, _npatch+1);
    allocateDispMatrix(FEM, _npatch+1);
  } else {    
    allocateDispMatrix(MAL, _npatch);
    allocateDispMatrix(FEM, _npatch);
  }
  
  TMatrix* mmat = _DispMatrix[MAL];
  TMatrix* fmat = _DispMatrix[FEM];
  
  //philopatry:
  double pmal = 1 - _mal_rate, pfem = 1 - _fem_rate;
  //migration:
  double mmal = _mal_rate/2, mfem = _fem_rate/2;
  
  fmat->assign(0);
  mmat->assign(0);
  
  //diagonal:
  for (unsigned int i = 0; i < _npatch; ++i){
    fmat->set(i, i, pfem);
    mmat->set(i, i, pmal);
  }
  //around the diagonal
  for (unsigned int i = 0; i < _npatch-1; ++i){
    fmat->set(i, i+1, mfem);
    fmat->set(i+1, i, mfem);
    mmat->set(i, i+1, mmal);
    mmat->set(i+1, i, mmal);
  }
  
  if(border_model == 3) {
    
    //set migration rates into the absorbing patch:
    fmat->set(0, _npatch, mfem);
    mmat->set(0, _npatch, mmal);
    fmat->set(_npatch -1, _npatch, mfem);
    mmat->set(_npatch -1, _npatch, mmal);
    fmat->set(_npatch, _npatch, 1);
    mmat->set(_npatch, _npatch, 1);
    
    if(!_isForward) { fmat->transpose();  mmat->transpose(); }
    
  } else if (border_model == 2) {
    
    //set migration rates of the two border patches:
    fmat->set(0, 1, 2*mfem);
    mmat->set(0, 1, 2*mmal);
    fmat->set(_npatch -1, _npatch -2, 2*mfem);
    mmat->set(_npatch -1, _npatch -2, 2*mmal);
    
    if(!_isForward) { fmat->transpose();  mmat->transpose(); }
    
  } else { //is a torus by default
    //the 2 last elements, this is a ring population!
    fmat->set(0, _npatch -1, mfem);
    mmat->set(0, _npatch -1, mmal);
    fmat->set(_npatch -1, 0, mfem);
    mmat->set(_npatch -1, 0, mmal);
  }
  return true;
}
// ----------------------------------------------------------------------------------------
// LCE_Disperse_base::setLatticeMatrix
// ----------------------------------------------------------------------------------------
/**Sets the dispersal matrices for the Lattice dispersal model.*/
bool LCE_Disperse_base::setLatticeMatrix()
{
#ifdef _DEBUG_
  message("setLatticeMatrix()\n");
#endif
  if(!_paramSet->isSet(_prefix + "_border_model")) {
    error("Missing parameter \"dispersal_border_model\" with dispersal model 4!\n");
    return false;
  }
  
  if(!_paramSet->isSet(_prefix + "_lattice_range") ) {
    error("Missing parameter \"dispersal_lattice_range\" with dispersal model 4!\n");
    return false;
  }
  
  unsigned int side = (unsigned int) sqrt((double)_npatch);
  
  if( side*side != _npatch ) {
    error("The number of patches is not a square number in the lattice dispersal model\n");
    return false;
  }
  
  /**Each matrix has 'patch number' x 'patch number' cells unless the lattice model is the absorbing
   boundaries model where we add the void patch.*/
  if((unsigned int)get_parameter_value(_prefix + "_border_model") == 3) {
    allocateDispMatrix(MAL, _npatch+1);
    allocateDispMatrix(FEM, _npatch+1);
  } else {
    allocateDispMatrix(MAL, _npatch);
    allocateDispMatrix(FEM, _npatch);
  }
  
  /**The "dispersal_lattice_range" parameter defines the number of neighbouring patches to disperse into. 
   Option 1 sets this number to 4 (left and right, up and down patches) whereas option 2 allows to disperse to the 8 neighbouring
   patches, including the patches on the diagonals.*/
  unsigned int range = (unsigned int)get_parameter_value(_prefix + "_lattice_range");
  //philopatry:
  double pmal = 1 - _mal_rate, pfem = 1 - _fem_rate;
  //migration:
  double mmal = _mal_rate/(range == 1 ? 4 : 8), mfem = _fem_rate/(range == 1 ? 4 : 8);
  
  setBasicLatticeMatrix(side, pmal, pfem, mmal, mfem);
  
  switch((unsigned int)get_parameter_value(_prefix + "_border_model")) {
    case 1:
      return setLatticeTorrusMatrix(side, mmal, mfem);
    case 2:
      return setLatticeReflectingMatrix(side);
    case 3:
      return setLatticeAbsorbingMatrix();
  }
  
  return true;
  
  //  TMatrix* fmat = _DispMatrix[FEM];
  //  for (unsigned int i=0; i<_npatch; ++i){
  //    for (unsigned int j=0; j<_npatch; ++j)
  //      cout<<fmat->get(i,j)<<"\t";
  //    cout<<endl;
  //  }
} 
// ----------------------------------------------------------------------------------------
// LCE_Disperse_base::setBasicLatticeMatrix
// ----------------------------------------------------------------------------------------
bool LCE_Disperse_base::setBasicLatticeMatrix(unsigned int side, double phi_mal, double phi_fem,
                                              double disp_mal, double disp_fem) 
{  
  unsigned int I, J;
  TMatrix* mmat = _DispMatrix[MAL];
  TMatrix* fmat = _DispMatrix[FEM];
  //init:
  for (unsigned int i=0; i<_npatch; ++i){
    for (unsigned int j=0; j<_npatch; ++j){
      fmat->set(i, j, 0.);
      mmat->set(i, j, 0.);
    }}
  //diagonal:
  for (unsigned int i = 0; i < _npatch; ++i){
    fmat->set(i, i, phi_fem);
    mmat->set(i, i, phi_mal);
  }
  //around the diagonal:
  for (unsigned int i = 0; i < _npatch-1; ++i){
    //left and right elements:
    fmat->set(i, i+1, disp_fem);
    fmat->set(i+1, i, disp_fem);
    mmat->set(i, i+1, disp_mal);
    mmat->set(i+1, i, disp_mal);
  }
  for (unsigned int i = 0; i < _npatch - side; ++i){
    //one row up and down elements:
    J = i + side;
    fmat->set(i, J, disp_fem);
    fmat->set(J, i, disp_fem);
    mmat->set(i, J, disp_mal);
    mmat->set(J, i, disp_mal);
  }
  //diagonal steps:
  if((unsigned int)get_parameter_value(_prefix + "_lattice_range") == 2) {
    unsigned int step = side + 1;
    for (unsigned int i = 0; i < _npatch - step; ++i){
      //diagonal step elements, to the right:
      J = i + step;
      fmat->set(i, J, disp_fem);
      fmat->set(J, i, disp_fem);
      mmat->set(i, J, disp_mal);
      mmat->set(J, i, disp_mal);
    }
    step = side - 1;
    for (unsigned int i = 0; i < _npatch - step; ++i){
      //diagonal step elements, to the left:
      J = i + step;
      fmat->set(i, J, disp_fem);
      fmat->set(J, i, disp_fem);
      mmat->set(i, J, disp_mal);
      mmat->set(J, i, disp_mal);
    }
    for (unsigned int i = 1; i < side - 1; ++i){
      //inconsistent diagonale steps in any models
      I = i * side - 1;
      J = (i + 1) * side;
      fmat->set(I, J, 0.0);
      fmat->set(J, I, 0.0);
      mmat->set(I, J, 0.0);
      mmat->set(J, I, 0.0);
    }
  }
  
  if((unsigned int)get_parameter_value(_prefix + "_border_model") != 1) {
    //remove corners
    for (unsigned int i = 0; i < side; ++i){
      //inconsistent rigth and left steps (these were set in the case lattice range = 2 and step = side -1 above)
      I = i * side;
      J = ((i + 1) * side) - 1;
      fmat->set(I, J, 0.0);
      fmat->set(J, I, 0.0);
      mmat->set(I, J, 0.0);
      mmat->set(J, I, 0.0);
    }
    for (unsigned int i = 1; i < side; ++i){
      //inconsistent diagonale steps
      I = i * side - 1;
      J = i * side;
      fmat->set(I, J, 0.0);
      fmat->set(J, I, 0.0);
      mmat->set(I, J, 0.0);
      mmat->set(J, I, 0.0);
    }
  }  
  return true;
}
// ----------------------------------------------------------------------------------------
// LCE_Disperse_base::setLatticeTorrusMatrix
// ----------------------------------------------------------------------------------------
bool LCE_Disperse_base::setLatticeTorrusMatrix(unsigned int side, double disp_mal, double disp_fem) 
{
  unsigned int I, J;
  
  TMatrix* mmat = _DispMatrix[MAL];
  TMatrix* fmat = _DispMatrix[FEM];
  
  for (unsigned int i = 0; i < side; ++i){
    //rigth and left steps across the borders
    I = i * side;
    J = ((i + 1) * side) - 1;
    fmat->set(I, J, disp_fem);
    fmat->set(J, I, disp_fem);
    mmat->set(I, J, disp_mal);
    mmat->set(J, I, disp_mal);
  }
  for (unsigned int i = 0; i < side; ++i){
    //up and down steps across the borders
    J = ((side - 1) * side) + i;
    fmat->set(i, J, disp_fem);
    fmat->set(J, i, disp_fem);
    mmat->set(i, J, disp_mal);
    mmat->set(J, i, disp_mal);
  }
  if((unsigned int)get_parameter_value(_prefix + "_lattice_range") != 2) {
    for (unsigned int i = 1; i < side; ++i){
      //remove diagonale steps across the borders
      I = i * side - 1;
      J = i * side;
      fmat->set(I, J, 0.0);
      fmat->set(J, I, 0.0);
      mmat->set(I, J, 0.0);
      mmat->set(J, I, 0.0);
    }
  } else {
    //add diagonale steps across the borders
    for (unsigned int i = 0; i < side - 1; ++i){
      I = i * side;
      J = ((i + 2) * side) - 1;
      fmat->set(I, J, disp_fem);
      fmat->set(J, I, disp_fem);
      mmat->set(I, J, disp_mal);
      mmat->set(J, I, disp_mal);
    }
    for (unsigned int i = 0; i < side-1; ++i){
      J = (side - 1) * side + i + 1;
      fmat->set(i, J, disp_fem);
      fmat->set(J, i, disp_fem);
      mmat->set(i, J, disp_mal);
      mmat->set(J, i, disp_mal);
    }
    for (unsigned int i = 1; i < side; ++i){
      J = (side - 1) * side + i - 1;
      fmat->set(i, J, disp_fem);
      fmat->set(J, i, disp_fem);
      mmat->set(i, J, disp_mal);
      mmat->set(J, i, disp_mal);
    }
    //add the corners:
    unsigned int corners[4] = {0, side -1, side*(side -1), _npatch -1};
    fmat->set(corners[0], corners[3], disp_fem);
    fmat->set(corners[3], corners[0], disp_fem);
    fmat->set(corners[1], corners[2], disp_fem);
    fmat->set(corners[2], corners[1], disp_fem);
    mmat->set(corners[0], corners[3], disp_mal);
    mmat->set(corners[3], corners[0], disp_mal);
    mmat->set(corners[1], corners[2], disp_mal);
    mmat->set(corners[2], corners[1], disp_mal);
  }
  return true;
}
// ----------------------------------------------------------------------------------------
// LCE_Disperse_base::setLatticeReflectingMatrix
// ----------------------------------------------------------------------------------------
bool LCE_Disperse_base::setLatticeReflectingMatrix(unsigned int side)
{
  TMatrix* mmat = _DispMatrix[MAL];
  TMatrix* fmat = _DispMatrix[FEM];
  
  double dm, df;
  unsigned int range = (unsigned int)get_parameter_value(_prefix + "_lattice_range");
  
  //the sum of the corner and border cells disp probs must be one!
  
  //CORNERS
  unsigned int corners[4] = {0, side -1, side*(side -1), _npatch -1};
  
  dm = _mal_rate/(range == 1 ? 2.0 : 3.0);
  df = _fem_rate/(range == 1 ? 2.0 : 3.0);
  
  for(unsigned int i = 0; i < 4; i++) {
    
    fmat->set(corners[i], (corners[i] % side ? corners[i] - 1 : corners[i] + 1), df);
    fmat->set(corners[i], (i < 2 ? corners[i] + side : corners[i] - side), df);
    
    mmat->set(corners[i], (corners[i] % side ? corners[i] - 1 : corners[i] + 1), dm);
    mmat->set(corners[i], (i < 2 ? corners[i] + side : corners[i] - side), dm);
  }
  
  if(range == 2) {
    fmat->set(corners[0], corners[0] + side + 1, df);
    fmat->set(corners[1], corners[1] + side - 1, df);
    fmat->set(corners[2], corners[2] - side + 1, df);
    fmat->set(corners[3], corners[3] - side - 1, df);
    mmat->set(corners[0], corners[0] + side + 1, dm);
    mmat->set(corners[1], corners[1] + side - 1, dm);
    mmat->set(corners[2], corners[2] - side + 1, dm);
    mmat->set(corners[3], corners[3] - side - 1, dm);
  }
  
  //BORDERS
  unsigned int nb_borders = (side - 2) * 4, nb_per_border = side - 2;
  unsigned int borders[nb_borders];
  
  for(unsigned int i = 0; i < nb_per_border; ++i) {
    borders[i] = i + 1; //top row
    borders[i + nb_per_border] = (i+1) * side; //left column
    borders[i + 2*nb_per_border] = (i+2) * side -1; //right column
    borders[i + 3*nb_per_border] = side*(side - 1) + i + 1; //bottom row
  }
  
  dm = _mal_rate/(range == 1 ? 3.0 : 5.0);
  df = _fem_rate/(range == 1 ? 3.0 : 5.0);
  
  for(unsigned int i = 0; i < nb_per_border; ++i) {
    //left and right steps
    //top border:
    fmat->set(borders[i], borders[i] + 1, df);
    fmat->set(borders[i], borders[i] - 1, df);
    //left border, move to the right only
    fmat->set(borders[i + nb_per_border], borders[i + nb_per_border] + 1, df);
    //right border, move to the left only
    fmat->set(borders[i + 2*nb_per_border], borders[i + 2*nb_per_border] - 1, df);
    //bottom border
    fmat->set(borders[i + 3*nb_per_border], borders[i + 3*nb_per_border] + 1, df);
    fmat->set(borders[i + 3*nb_per_border], borders[i + 3*nb_per_border] - 1, df);
    
    mmat->set(borders[i], borders[i] + 1, dm);
    mmat->set(borders[i], borders[i] - 1, dm);
    mmat->set(borders[i + nb_per_border], borders[i+ nb_per_border] + 1, dm);
    mmat->set(borders[i + 2*nb_per_border], borders[i + 2*nb_per_border] - 1, dm);
    mmat->set(borders[i + 3*nb_per_border], borders[i + 3*nb_per_border] + 1, dm);
    mmat->set(borders[i + 3*nb_per_border], borders[i + 3*nb_per_border] - 1, dm);
    //up and down
    fmat->set(borders[i], borders[i] + side, df);
    fmat->set(borders[i + nb_per_border], borders[i + nb_per_border] + side, df);
    fmat->set(borders[i + nb_per_border], borders[i + nb_per_border] - side, df);
    fmat->set(borders[i + 2*nb_per_border], borders[i + 2*nb_per_border] - side, df);
    fmat->set(borders[i + 2*nb_per_border], borders[i + 2*nb_per_border] + side, df);
    fmat->set(borders[i + 3*nb_per_border], borders[i + 3*nb_per_border] - side, df);
    
    mmat->set(borders[i], borders[i] + side, dm);
    mmat->set(borders[i + nb_per_border], borders[i + nb_per_border] + side, dm);
    mmat->set(borders[i + nb_per_border], borders[i + nb_per_border] - side, dm);
    mmat->set(borders[i + 2*nb_per_border], borders[i + 2*nb_per_border] - side, dm);
    mmat->set(borders[i + 2*nb_per_border], borders[i + 2*nb_per_border] + side, dm);
    mmat->set(borders[i + 3*nb_per_border], borders[i + 3*nb_per_border] - side, dm);
  }
  
  if(range == 2) {
    for(unsigned int i = 0; i < nb_per_border; ++i) {
      fmat->set(borders[i], borders[i] + side + 1, df);
      fmat->set(borders[i], borders[i] + side - 1, df);
      fmat->set(borders[i + nb_per_border], borders[i + nb_per_border] - side + 1, df);
      fmat->set(borders[i + nb_per_border], borders[i + nb_per_border] + side + 1, df);
      fmat->set(borders[i + 2*nb_per_border], borders[i + 2*nb_per_border] - side - 1, df);
      fmat->set(borders[i + 2*nb_per_border], borders[i + 2*nb_per_border] + side - 1, df);
      fmat->set(borders[i + 3*nb_per_border], borders[i + 3*nb_per_border] - side + 1, df);
      fmat->set(borders[i + 3*nb_per_border], borders[i + 3*nb_per_border] - side - 1, df);
      
      mmat->set(borders[i], borders[i] + side + 1, dm);
      mmat->set(borders[i], borders[i] + side - 1, dm);
      mmat->set(borders[i + nb_per_border], borders[i + nb_per_border] - side + 1, dm);
      mmat->set(borders[i + nb_per_border], borders[i + nb_per_border] + side + 1, dm);
      mmat->set(borders[i + 2*nb_per_border], borders[i + 2*nb_per_border] - side - 1, dm);
      mmat->set(borders[i + 2*nb_per_border], borders[i + 2*nb_per_border] + side - 1, dm);
      mmat->set(borders[i + 3*nb_per_border], borders[i + 3*nb_per_border] - side + 1, dm);
      mmat->set(borders[i + 3*nb_per_border], borders[i + 3*nb_per_border] - side - 1, dm);
    }    
  }
  
  if(!_isForward) {  fmat->transpose();  mmat->transpose(); }
  
  return true;
}
// ----------------------------------------------------------------------------------------
// LCE_Disperse_base::setLatticeAbsorbingMatrix
// ----------------------------------------------------------------------------------------
bool LCE_Disperse_base::setLatticeAbsorbingMatrix()
{
  TMatrix* mmat = _DispMatrix[MAL];
  TMatrix* fmat = _DispMatrix[FEM];
  
  //set the absorbing patch probs to 1
  for(unsigned int i = 0; i < _npatch+1; ++i) {
    fmat->set(i, _npatch, 1.0);
    mmat->set(i, _npatch, 1.0);
  }
  if(!_isForward) {  fmat->transpose();  mmat->transpose(); }
  
  return true;
}

// ----------------------------------------------------------------------------------------
// LCE_Disperse_base::setLatticeReducedDispMatrix
// ----------------------------------------------------------------------------------------

bool LCE_Disperse_base::setReducedDispMatrix()
{



	if(get_parameter(_prefix + "_patchConnected")->isSet()){ //@TODO remove this, parameter non longer exist

		if(get_parameter(_prefix + "_pDisp")->isSet()){

			if(!_patchConnected) _patchConnected = new TMatrix;

			if(!_pDisp) _pDisp = new TMatrix;


			unsigned int nb_patch_connected;
			unsigned int nb_patch;

			nb_patch_connected=_patchConnected->getNbCols();

			nb_patch=_patchConnected->getNbRows();

			//some check

			if(_npatch != nb_patch){cout<<"matrix dim must equal patch number"<<endl;}

			if(_reducedDispMat[0].size() != 0) _reducedDispMat[0].clear();
			if(_reducedDispMatProba[0].size() != 0) _reducedDispMatProba[0].clear();
			
			if(_reducedDispMat[1].size() != 0) _reducedDispMat[1].clear();
			if(_reducedDispMatProba[1].size() != 0) _reducedDispMatProba[1].clear();

			for(unsigned int i = 0; i < _npatch ;i++){

				_reducedDispMat[0].push_back(vector<double>());

				_reducedDispMat[1].push_back(vector<double>());

				_reducedDispMatProba[0].push_back(vector<double>());

				_reducedDispMatProba[1].push_back(vector<double>());

				for (unsigned int j = 0; j < nb_patch_connected; ++j) {

				      _reducedDispMat[0][i].push_back(_patchConnected->get(i,j));

				      _reducedDispMat[1][i].push_back(_patchConnected->get(i,j));

				      _reducedDispMatProba[0][i].push_back(_pDisp->get(i,j));

				      _reducedDispMatProba[1][i].push_back(_pDisp->get(i,j));
				}
			}
		}

		else {cout<<"**error** both patchConnected and pDisp must be set"<<endl;}
	}

	else {

    _reducedDispMat[0].clear();
    _reducedDispMat[1].clear();
    _reducedDispMatProba[0].clear();
    _reducedDispMatProba[1].clear();
    
  for (unsigned int i = 0; i < _npatch; ++i) {

    _reducedDispMat[0].push_back(vector<double>());
    _reducedDispMat[1].push_back(vector<double>());
    
    _reducedDispMatProba[0].push_back(vector<double>());
    _reducedDispMatProba[1].push_back(vector<double>());
    
    for (unsigned int j = 0; j < _npatch; ++j) {

      if(_DispMatrix[0]->get(i, j) != 0) {
        _reducedDispMat[0][i].push_back(j);
        _reducedDispMatProba[0][i].push_back(_DispMatrix[0]->get(i, j));
      }
      
      if(_DispMatrix[1]->get(i, j) != 0) {
        _reducedDispMat[1][i].push_back(j);
        _reducedDispMatProba[1][i].push_back(_DispMatrix[1]->get(i, j));
      }
    }
  }

}

  return true;
}

// ----------------------------------------------------------------------------------------
// LCE_Disperse_base::Migrate
// ----------------------------------------------------------------------------------------
unsigned int LCE_Disperse_base::getMigrationPatchForward (sex_t SEX, unsigned int LocalPatch)
{
  
  //if(!LCE_Disperse_base::setReducedDispMatrix()) cout<<"setReducedMatrix not set"<<endl;
  
  double sum = 0, random = RAND::Uniform();
  unsigned int AimedPatch = 0;
  
  if(random > 0.999999) random = 0.999999;//this to avoid overflows when random == 1
  
  sum =  _reducedDispMatProba[SEX][LocalPatch][AimedPatch];
  
  //cout<<"sum 1 "<<sum<<endl;
  //cout<<"random "<<random<<endl;
  
  
  while (random > sum) {
  
    AimedPatch++;
    
    sum +=  _reducedDispMatProba[SEX][LocalPatch][AimedPatch];
    
    //cout<<"sum "<<sum<<endl;
    
  }

//cout<< "disperse from"<< LocalPatch <<"to "<<_reducedDispMat[SEX][LocalPatch][AimedPatch]<<endl;

  return _reducedDispMat[SEX][LocalPatch][AimedPatch]-1;
}



// ----------------------------------------------------------------------------------------
// LCE_Disperse_base::Migrate
// ----------------------------------------------------------------------------------------
unsigned int LCE_Disperse_base::getMigrationPatchBackward (sex_t SEX, unsigned int LocalPatch)
{
  double sum = 0, random = RAND::Uniform();
  unsigned int SourcePatch = 0;
  
  if(random > 0.999999) random = 0.999999;//this to avoid overflows when random == 1

  sum = _reducedDispMatProba[SEX][LocalPatch][SourcePatch];
  
  while (random > sum) {
    SourcePatch++;
    sum += _reducedDispMatProba[SEX][LocalPatch][SourcePatch];
  }
    
  return _reducedDispMat[SEX][LocalPatch][SourcePatch];
}
/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

//                             ******** LCE_Disperse_ConstDisp ********/

// ----------------------------------------------------------------------------------------
// LCE_Disperse_ConstDisp
// ----------------------------------------------------------------------------------------
LCE_Disperse_ConstDisp::LCE_Disperse_ConstDisp () 
: LifeCycleEvent ("disperse",""), doMigration(0), doPatchMigration(0)
{
//  cout << "calling LCE_Disperse_ConstDisp()\n";
  
  ParamUpdater< LCE_Disperse_ConstDisp > * updater = 
  new ParamUpdater< LCE_Disperse_ConstDisp > (&LCE_Disperse_ConstDisp::setParameters);
  
  LCE_Disperse_base::addParameters("dispersal", updater);
}
// ----------------------------------------------------------------------------------------
// LCE_Disperse_ConstDisp::setParameters
// ----------------------------------------------------------------------------------------
bool LCE_Disperse_ConstDisp::setParameters(string prefix)
{


  
  if(!LCE_Disperse_base::setBaseParameters(prefix)) return false;
  
  switch ( getDispersalModel() ) {
    case 0: //dispersal matrix given in input
    case 1:
    case 3:
    case 4:
      doMigration = &LCE_Disperse_ConstDisp::Migrate;
      break;
    case 2:
      doMigration = &LCE_Disperse_ConstDisp::Migrate_propagule;
      break;
    default: {
      error("\nDispersal model '%i' not yet implemented\n",getDispersalModel());
      return false;
    }
  }
  
  switch ((int)get_parameter_value(prefix + "_border_model")) {
    case -1:
    case 1:
    case 2:
      doPatchMigration = &LCE_Disperse_ConstDisp::MigratePatch;
      break;
    case 3:
      doPatchMigration = &LCE_Disperse_ConstDisp::MigratePatch_AbsorbingBorder;
      break;
    default:{
      error("\nDispersal border model '%i' not yet implemented\n",
            (int)get_parameter_value(prefix + "_border_model"));
      return false;
    }
  }
  
  cout<<"ConstDisp::setParameter - done"<<endl;
  
  return true;
}
// ----------------------------------------------------------------------------------------
// LCE_Disperse_ConstDisp::execute
// ----------------------------------------------------------------------------------------
void LCE_Disperse_ConstDisp::execute ()
{

	 //check//
	  //cout<<"disperse_starting"<<endl;
	  //


#ifdef _DEBUG_
  message("LCE_Disperse_ConstDisp::execute (Patch nb: %i offsprg nb: %i adlt nb: %i "
          ,_popPtr->getPatchNbr(), _popPtr->size( OFFSPRG ), _popPtr->size( ADULTS ));
#endif

  //check whether the number of patches changed during simulation:
  if(_npatch != _popPtr->getPatchNbr()) {

	  //check
	 // cout<<"TRUE"<<endl;
	  //

    _npatch = _popPtr->getPatchNbr();
    if(!updateDispMatrix()) fatal("bailing out\n");
  }
  
 

  reset_counters();
  
	//check//
 	//cout<<"migrate execute"<<endl;
  	//

  (this->*doMigration)();
  
#ifdef _DEBUG_
  unsigned int c = 0;
  for(unsigned int i = 0; i < _npatch; i++)
    c += _popPtr->getPatch(i)->nbEmigrant;
  message("emigrants nb: %i)\n",c);
#endif
}
// ----------------------------------------------------------------------------------------
// LCE_Disperse_ConstDisp::Migrate
// ----------------------------------------------------------------------------------------
void LCE_Disperse_ConstDisp::Migrate ()
{

//cout<<"check8"<<endl;
  Patch *patch;
  
  for(unsigned int i = 0; i < _npatch; i++) {
  
  	//cout<<"patch to do "<<i<<endl;
    
    patch = _popPtr->getPatch(i);
    
    (this->*doPatchMigration)(FEM, i);
    
    (this->*doPatchMigration)(MAL, i);
    
    //cout<<"patch done "<<i<<endl;
    
    //set coloniser counter
    
    if(patch->get_isExtinct())
     
      patch->nbKolonisers = _postDispPop[i]->size(OFFSx);
      
    else 
    
      patch->nbKolonisers = -1; //value used in stat_demo to check for colonization
        
  }//end for nb_patch
  
  //put back the indviduals into the offspring container
  swapPostDisp();
  
  //cout<<"check9"<<endl;
}
// ----------------------------------------------------------------------------------------
// LCE_Disperse_ConstDisp
// ----------------------------------------------------------------------------------------
void LCE_Disperse_ConstDisp::Migrate_propagule ()
{
  setIsland_PropagulePool_Matrix();
  Migrate();
}
// ----------------------------------------------------------------------------------------
// LCE_Disperse_ConstDisp::MigratePatch
// ----------------------------------------------------------------------------------------
void LCE_Disperse_ConstDisp::MigratePatch (sex_t SEX, unsigned int LocalPatch)
{

//cout<<"check10"<<endl;

  Patch* patch = _popPtr->getPatch(LocalPatch);
  unsigned int AimedPatch;
  unsigned int Limit = _npatch -1;
  Individual* ind;
  
  //cout<<"check11"<<endl;
  
  
  while( patch->size(SEX, OFFSx) != 0 ) {
  
    //cout<<"check12"<<endl;
    
    
    do{
      
      
      AimedPatch = getMigrationPatchForward(SEX, LocalPatch);
      
      //cout<<"local patch (in MigratePatch)"<<LocalPatch<<endl;
      //cout<<"Aimed patch (in MigratePatch)"<<AimedPatch<<endl;
    
    }
    
    while ( AimedPatch > Limit );
    
   
    ind = patch->remove(SEX, OFFSx, 0);
    
    _postDispPop[AimedPatch]->add(SEX, OFFSx, ind);
       
       
    if(LocalPatch != AimedPatch) {
    
      patch->nbEmigrant++;
      
      _popPtr->getPatch(AimedPatch)->nbImigrant++;
      
    } else
    
      patch->nbPhilopat++;
       
  }//end while

}
// ----------------------------------------------------------------------------------------
// MigratePatch_AbsorbingBorder
// ----------------------------------------------------------------------------------------
void LCE_Disperse_ConstDisp::MigratePatch_AbsorbingBorder (sex_t SEX, unsigned int LocalPatch)
{   
  Patch* patch = _popPtr->getPatch(LocalPatch);
  unsigned int AimedPatch; 
  Individual* ind;
  
  while( patch->size(SEX, OFFSx) != 0 ) {
   
    do{
      AimedPatch = getMigrationPatchForward(SEX, LocalPatch);
    }while ( AimedPatch > _npatch );
    
    if(AimedPatch == _npatch) {
      _popPtr->recycle( patch->get(SEX, OFFSx, 0) );
      patch->remove(SEX, OFFSx, 0);
    } else {
      ind = patch->remove(SEX, OFFSx, 0);
      _postDispPop[AimedPatch]->add(SEX, OFFSx, ind);
    }      
//      _popPtr->move(SEX, OFFSx, LocalPatch, PDISPx, AimedPatch, 0);
    
    if(LocalPatch != AimedPatch) {
      patch->nbEmigrant++;
      if(AimedPatch < _npatch)
        _popPtr->getPatch(AimedPatch)->nbImigrant++;
    } else
      patch->nbPhilopat++;
    
  }//end while
}
/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

//                             ******** LCE_SeedDisp ********/

// ----------------------------------------------------------------------------------------
LCE_SeedDisp::LCE_SeedDisp () : LifeCycleEvent ("seed_disp","")
{
  ParamUpdater<LCE_SeedDisp> * updater =
  new ParamUpdater<LCE_SeedDisp> (&LCE_SeedDisp::setParameters);
  
//  cout << "calling LCE_SeedDisp()\n";
  
  set_event_name("seed_disp"); //this resets the paramset, erases the params added by LCE_Disperse_ConstDisp

  addParameters("seed_disp", updater);

//  get_paramset()->show_up();
}
/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

//                             ******** LCE_Disperse_EvolDisp ********/

// ----------------------------------------------------------------------------------------
LCE_Disperse_EvolDisp::LCE_Disperse_EvolDisp () 
: LifeCycleEvent("disperse_evoldisp",""), _fem_cost(-1.0), _mal_cost(-1.0),
_fdisp_trait_link(0), _mdisp_trait_link(0)
{
  ParamUpdater<LCE_Disperse_EvolDisp> * updater =
  new ParamUpdater<LCE_Disperse_EvolDisp> (&LCE_Disperse_EvolDisp::setParameters);

  LCE_Disperse_base::addParameters("dispersal", updater);
  
  add_parameter("dispersal_cost",DBL,false,true,0,1, updater);
  add_parameter("dispersal_cost_fem",DBL,false,true,0,1, updater);
  add_parameter("dispersal_cost_mal",DBL,false,true,0,1, updater);
  add_parameter("dispersal_fixed_trait",STR,false,false,0,0, updater);
  add_parameter("dispersal_fixed_rate",DBL,false,true,0,1, updater);
}
// ----------------------------------------------------------------------------------------
// LCE_Disperse_EvolDisp::setParameters
// ----------------------------------------------------------------------------------------
bool LCE_Disperse_EvolDisp::setParameters ()
{
  //we do not call the LCE_Disperse_base::setParameters here because we don't use
  //dispersal matrices in input
  _npatch = _popPtr->getPatchNbr();
  
  setPostDispPop();
  
  _prefix = "dispersal";
  
  _disp_model = (int)_paramSet->getValue("dispersal_model");
  
  _disp_propagule_prob = _paramSet->getValue("dispersal_propagule_prob");
  
  if(_disp_model == -1) {
    error("dispersal model not specified!\n");
    return false;
  }
  
  if(_disp_model == 2 && _disp_propagule_prob == -1) {
    error("dispersal propagule probability is missing!\n");
    return false;
  }
  
  _fdisp_trait_link = _popPtr->getTraitIndex("fdisp");
  _mdisp_trait_link = _popPtr->getTraitIndex("mdisp");
  
  if(_paramSet->isSet("dispersal_cost")) {
    
    _fem_cost = _mal_cost = _paramSet->getValue("dispersal_cost");
    
  } else if(_paramSet->isSet("dispersal_cost_fem") && _paramSet->isSet("dispersal_cost_mal")) {
    
    _fem_cost = _paramSet->getValue("dispersal_cost_fem");
    
    _mal_cost = _paramSet->getValue("dispersal_cost_mal");
    
  } else {
    error("dispersal cost params are not set !\n");
    return false;
  }
  
  if(_paramSet->isSet("dispersal_fixed_trait")) {
    
    if(_paramSet->getArg("dispersal_fixed_trait").compare("female") == 0) {
      
      exec = &LCE_Disperse_EvolDisp::exec_evolmale;
      
    } else if(_paramSet->getArg("dispersal_fixed_trait").compare("male") == 0) {
      
      exec = &LCE_Disperse_EvolDisp::exec_evolfemale;
      
    } else {
      error("wrong argument value for \"dispersal_fixed_trait\"\n");
      return false;
    }
    
    _fixed_disp_rate = _paramSet->getValue("dispersal_fixed_rate");
    
  } else
    exec = &LCE_Disperse_EvolDisp::exec_evol2sex;
  
  switch ( getDispersalModel() ) {
    case 0:
    case 1:
      getAimedPatch = &LCE_Disperse_EvolDisp::Migrate_Island;
      break;
    case 2:
      getAimedPatch = &LCE_Disperse_EvolDisp::Migrate_Island_Propagule;
      break;
    case 3:
      getAimedPatch = &LCE_Disperse_EvolDisp::Migrate_SteppingStone1D;
      break;
    default: {
      /**@todo implement lattice dispersal models with evolving dispersal.*/
      error("\nDispersal model %i not yet implemented with evolving dispersal!\n",
            getDispersalModel());
      return false;
    }
  }
  return true;
}
// ----------------------------------------------------------------------------------------
// LCE_Disperse_EvolDisp::execute
// ----------------------------------------------------------------------------------------
void LCE_Disperse_EvolDisp::execute ()
{  
#ifdef _DEBUG_
  message("LCE_Disperse_EvolDisp::execute (Patch nb: %i offsprg nbr: %i)\n"
          ,_popPtr->getPatchNbr(),_popPtr->size( OFFSPRG ));
#endif
  if( getDispersalModel() == 2 ) setPropaguleTargets();
  
  reset_counters();
  
  (this->*exec)();
  
  Patch *current_patch;  
  
  _npatch = _popPtr->getPatchNbr();
  
  for(unsigned int i = 0; i < _npatch; i++) {
    
    current_patch = _popPtr->getPatch(i);
    
    //set coloniser counter
    if(current_patch->get_isExtinct()) 
      current_patch->nbKolonisers = _postDispPop[i]->size(OFFSx);
    else 
      current_patch->nbKolonisers = -1;
    
  }//end_for
  
  swapPostDisp();
  
}
// ----------------------------------------------------------------------------------------
// exec_evolfemale
// ----------------------------------------------------------------------------------------
void LCE_Disperse_EvolDisp::exec_evolfemale ()
{
  evoldisp(FEM, _fdisp_trait_link, _fem_cost);
  fixdisp(MAL, _fixed_disp_rate, _mal_cost);
}
// ----------------------------------------------------------------------------------------
// exec_evolmale
// ----------------------------------------------------------------------------------------
void LCE_Disperse_EvolDisp::exec_evolmale ()
{
  evoldisp(MAL, _mdisp_trait_link, _mal_cost);
  fixdisp(FEM, _fixed_disp_rate, _fem_cost);
}
// ----------------------------------------------------------------------------------------
// exec_evol2sex
// ----------------------------------------------------------------------------------------
void LCE_Disperse_EvolDisp::exec_evol2sex ()
{
  evoldisp(FEM, _fdisp_trait_link, _fem_cost);
  evoldisp(MAL, _mdisp_trait_link, _mal_cost);
}
// ----------------------------------------------------------------------------------------
// evoldisp
// ----------------------------------------------------------------------------------------
void LCE_Disperse_EvolDisp::evoldisp (sex_t SEX, int trait_link, double cost)
{
  
  Patch *current_patch;
  unsigned int AimedPatch;
  
  for(unsigned int i = 0; i < _npatch; i++) {
    
    current_patch = _popPtr->getPatch(i);
    
    while( current_patch->size(SEX, OFFSx) != 0 ) {
      
      if(RAND::Uniform() < *(double*)current_patch->get(SEX, OFFSx, 0)->getTraitValue(trait_link)) {
        //this one disperses
        AimedPatch = (this->*getAimedPatch)(i);
        
        current_patch->nbEmigrant++;
        
        if(RAND::Uniform() > cost) {
          //survives
          
          _postDispPop[AimedPatch]->add(SEX, OFFSx, current_patch->remove(SEX, OFFSx, 0) );
          
          _popPtr->getPatch(AimedPatch)->nbImigrant++;
          
        } else {
          
          _popPtr->recycle( current_patch->remove(SEX, OFFSx, 0) );
          
        }
      } else {
        //no dispersal
        _postDispPop[i]->add(SEX, OFFSx, current_patch->remove(SEX, OFFSx, 0) );
//        current_patch->move(SEX,OFFSx,PDISPx,0);
        current_patch->nbPhilopat++;
      }
      
    }//end while
  }//end for
}
// ----------------------------------------------------------------------------------------
// fixdisp
// ----------------------------------------------------------------------------------------
void LCE_Disperse_EvolDisp::fixdisp (sex_t SEX, double rate, double cost)
{
  
  Patch *current_patch;
  unsigned int AimedPatch;
  
  for(unsigned int i = 0; i < _npatch; i++) {
    
    current_patch = _popPtr->getPatch(i);
    
    while( current_patch->size(SEX, OFFSx) != 0 ) {
      
      if(RAND::Uniform() < rate) {
        //this one disperses
        AimedPatch = (this->*getAimedPatch)(i);
        
        current_patch->nbEmigrant++;
        
        if(RAND::Uniform() > cost) {
          //survives
          _postDispPop[AimedPatch]->add(SEX, OFFSx, current_patch->remove(SEX, OFFSx, 0) );
          
          _popPtr->getPatch(AimedPatch)->nbImigrant++;
          
        } else {
          
          _popPtr->recycle( current_patch->remove(SEX, OFFSx, 0) );

        }
      } else {
        //no dispersal
        _postDispPop[i]->add(SEX, OFFSx, current_patch->remove(SEX, OFFSx, 0) );
//        current_patch->move(SEX,OFFSx,PDISPx,0);
        current_patch->nbPhilopat++;
      }
      
    }//end while
  }//end for
}
// ----------------------------------------------------------------------------------------
// Migrate_Island_EvolDisp
// ----------------------------------------------------------------------------------------
unsigned int LCE_Disperse_EvolDisp::Migrate_Island(unsigned int home)
{
  unsigned int AimedPatch;
  //assign a Patch of arrival at random
  do{
    AimedPatch = RAND::Uniform(_npatch);
  }while(AimedPatch == home);
  
  return AimedPatch;
}
// ----------------------------------------------------------------------------------------
// Migrate_Island_Propagule_EvolDisp
// ----------------------------------------------------------------------------------------
unsigned int LCE_Disperse_EvolDisp::Migrate_Island_Propagule(unsigned int home)
{
  unsigned int AimedPatch, PropaguleTarget = getPropaguleTarget(home);
  
  if(!(RAND::Uniform() > getPropaguleProb()) )
    AimedPatch = PropaguleTarget;
  else
    do{
      AimedPatch = RAND::Uniform(_npatch);
    }while(AimedPatch == home || AimedPatch == PropaguleTarget);
  
  return AimedPatch;
}
// ----------------------------------------------------------------------------------------
// Migrate_SteppingStone1D_EvolDisp
// ----------------------------------------------------------------------------------------
unsigned int LCE_Disperse_EvolDisp::Migrate_SteppingStone1D(unsigned int home)
{
  int neighbours[2] = {(int)(home - 1),(int)(home + 1)};
  //if we are at one of the ends of the Patch array, we may migrate to the other end:
  if(neighbours[0] < 0) neighbours[0] = _npatch -1;
  else if(neighbours[1] == (int)_npatch) neighbours[1] = 0;
  
  return(RAND::RandBool() ? neighbours[0] : neighbours[1]);
}
