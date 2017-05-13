/**  $Id: LCEselection.cc,v 1.13.2.9 2016-04-28 13:06:01 fred Exp $

*  @file LCEselection.cc
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
*  created on @date 09.10.2007
*
*  @author fred
*/

#include <sstream>
#include "LCEselection.h"
#include "Uniform.h"
#include "utils.h"
#include "tstring.h"




/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

//                             ******** LCE_Selection_base ********/

// ----------------------------------------------------------------------------------------
// LCE_Selection_base::LCE_Selection_base
// ----------------------------------------------------------------------------------------
LCE_Selection_base::LCE_Selection_base ( ) : LifeCycleEvent("viability_selection", ""),
_selection_matrix(0), _gsl_selection_matrix(0), _diffs(0), _res1(0), _local_optima(0), _phe(0),
_selectTraitDimension(1), _base_fitness(1), _mean_fitness(0), _max_fitness(0), _scaling_factor(1),
_is_local(0), _is_absolute(1), _eVariance(0), _getRawFitness(0), _getFitness(0), _stater(0),_age_classes_selected(0),
_selection_age_variance(0), _getFitnessAge(0), _getRawFitnessAge(0), _rchange(0)
{
  add_parameter("selection_trait",STR,false,false,0,0); //no updaters here, to keep it safe...

  ParamUpdater< LCE_Selection_base > * updater =
	new ParamUpdater<LCE_Selection_base>(&LCE_Selection_base::set_fit_model);

  add_parameter("selection_fitness_model",STR,false,false,0,0, updater);

  updater = new ParamUpdater<LCE_Selection_base>(&LCE_Selection_base::set_sel_model);

  add_parameter("selection_model",STR,false,false,0,0, updater);
  add_parameter("selection_matrix",MAT,false,false,0,0, updater);
  add_parameter("selection_variance",DBL,false,false,0,0, updater);
  add_parameter("selection_correlation",DBL,false,false,0,0, updater);
  add_parameter("selection_trait_dimension",INT,false,false,0,0, updater);
  add_parameter("selection_base_fitness",DBL,false,true,0,1, updater);
  add_parameter("selection_lethal_equivalents",DBL,false,false,0,0, updater);
  add_parameter("selection_pedigree_F",MAT,false,false,0,0, updater);

  //updater = new ParamUpdater<LCE_Selection_base>(&LCE_Selection_base::setSelectionMatrixAge() );
  add_parameter("selection_local_optima",DBL,false,false,0,0, updater);

  add_parameter("selection_randomize",BOOL,false,false,0,0, updater);
  add_parameter("selection_environmental_variance", DBL, false, false, 0, 0, updater);
  add_parameter("rate_environmental_change",MAT,false,false,0,0,updater);
  add_parameter("age_classes_selected",MAT,false,false,0,0,updater);
  add_parameter("selection_age_variance",MAT,false,false,0,0,updater);

}
// ----------------------------------------------------------------------------------------
// LCE_Selection_base::~LCE_Selection_base
// ----------------------------------------------------------------------------------------
LCE_Selection_base::~LCE_Selection_base ( )
{
  for(unsigned int i = 0; i < _selection_matrix.size(); ++i)
    delete _selection_matrix[i];
  for(unsigned int i  = 0; i < _gsl_selection_matrix.size(); ++i)
    gsl_matrix_free(_gsl_selection_matrix[i]);
  if(_local_optima) delete _local_optima;
  if(_diffs) gsl_vector_free(_diffs);
  if(_res1) gsl_vector_free(_res1);
//  if(_phe != 0) delete [] _phe; !!! _phe is never allocated, just a name-holder
  if(_stater) delete _stater;
}

// ----------------------------------------------------------------------------------------
// LCE_Selection_base::loadStatServices
// ----------------------------------------------------------------------------------------
void LCE_Selection_base::loadStatServices ( StatServices* loader )
{
  if(_stater == NULL) _stater = new LCE_SelectionSH(this);
  loader->attach(_stater);
}
// ----------------------------------------------------------------------------------------
// LCE_Selection_base::setParameters
// ----------------------------------------------------------------------------------------
bool LCE_Selection_base::setParameters ( )
{

  if(!attach_trait(get_parameter("selection_trait")->getArg())) return false;


  if(!set_fit_model()) return false;

  if(!set_sel_model()) return false;

  _mean_fitness = _max_fitness = 0;
  _nb_trait = _popPtr->getIndividualProtoype()->getTraitNumber();
  _scaling_factor = 1;
  resetCounters();

  return true;
}
// ----------------------------------------------------------------------------------------
// LCE_Selection_base::set_fit_model
// ----------------------------------------------------------------------------------------
bool LCE_Selection_base::set_fit_model()
{

  if(_paramSet->isSet("selection_fitness_model")) {

    string fit_model = _paramSet->getArg("selection_fitness_model");

    if(fit_model.compare("absolute") == 0) {

      _is_local = false;
      _is_absolute = true;

      if(get_parameter("age_classes_selected")->isSet()){

    	  _getFitnessAge = &LCE_Selection_base::getFitnessAbsoluteAge;
      }
      else {
        _getFitness = &LCE_Selection_base::getFitnessAbsolute;
      }

      _setScalingFactor = &LCE_Selection_base::setScalingFactorAbsolute;

    } else if(fit_model.compare("relative_local") == 0) {

      _is_local = true;
      _is_absolute = false;
      _getFitness = &LCE_Selection_base::getFitnessRelative;
      _setScalingFactor = &LCE_Selection_base::setScalingFactorLocal;

    } else if(fit_model.compare("relative_global") == 0) {

      _is_local = false;
      _is_absolute = false;
      _getFitness = &LCE_Selection_base::getFitnessRelative;
      _setScalingFactor = &LCE_Selection_base::setScalingFactorGlobal;

    } else {
      error("Unknown fitness model \"%s\"", fit_model.c_str());
      return false;
    }
  } //default case:
  else {


    _is_local = false;
    _is_absolute = true;

    if(get_parameter("age_classes_selected")->isSet()){

      _getFitnessAge = &LCE_Selection_base::getFitnessAbsoluteAge;
    }

    else{

      _getFitness = &LCE_Selection_base::getFitnessAbsolute;

    }
    _setScalingFactor = &LCE_Selection_base::setScalingFactorAbsolute;


  }

  //cout<<"2"<<endl;

  return true;
}
// ----------------------------------------------------------------------------------------
// LCE_Selection_base::set_sel_model
// ----------------------------------------------------------------------------------------
bool LCE_Selection_base::set_sel_model()
{

  if(get_parameter("age_classes_selected")->isSet()){

    if(_age_classes_selected) delete _age_classes_selected;
    if(_selection_age_variance) delete _selection_age_variance;

    _age_classes_selected = new TMatrix;
    _selection_age_variance = new TMatrix;

  	get_parameter("age_classes_selected")->getMatrix(_age_classes_selected);
  	get_parameter("selection_age_variance")->getMatrix(_selection_age_variance);

//  	if(_age_classes_selected->getNbCols() != _selection_age_variance->getNbCols()){
//
//  		TMatrix* temp;
//  		temp = new TMatrix(1,_age_classes_selected->getNbCols());
//
//      for(unsigned int i=0; i<_age_classes_selected->getNbCols();i++){
//  			//all variances and identical equal to the first element
//				temp->set(0,i,_selection_age_variance->get(0,0));
//      }
//
//      if(_selection_age_variance) delete _selection_age_variance;
//      _selection_age_variance = temp;
//
//      //_selection_age_variance->show_up();
//
//      //temp->~TMatrix ()   ;
//  	}

//    unsigned int nb_classes_selected = _age_classes_selected->getNbCols();


  }


  ////

  if(get_parameter("rate_environmental_change")->isSet()){

    if(_rchange) delete _rchange;

		_rchange = new TMatrix;

		get_parameter("rate_environmental_change")->getMatrix(_rchange);


  }


  _selectTraitDimension = (int)get_parameter_value("selection_trait_dimension");

  _base_fitness = get_parameter_value("selection_base_fitness");

  if(get_parameter("selection_environmental_variance")->isSet())
    //the variable actually holds the standard dev...
    _eVariance =  sqrt(get_parameter_value("selection_environmental_variance"));
  else
    _eVariance = 0;

  if(!_paramSet->isSet("selection_model")) {


    _getRawFitness = &LCE_Selection_base::getFitnessDirect;

  } else {

    string sel_model = _paramSet->getArg("selection_model");

    if(sel_model == "fix") {

      if(!get_parameter("selection_lethal_equivalents")->isSet()) {
        error("\"selection_lethal_equivalents\" parameter is missing with \"fix\" selection!\n");
        return false;
      } else
        _letheq = get_parameter_value("selection_lethal_equivalents");

      if(!get_parameter("selection_base_fitness")->isSet()) {
        warning("\"selection_base_fitness\" parameter is missing under fix selection model, setting it to 1.\n");
        _base_fitness = 1.0;
      }

      if(!get_parameter("selection_pedigree_F")->isSet()) {
        error("\"selection_pedigree_F\" parameter is missing with \"fix\" selection!\n");
        return false;
      } else {
        TMatrix tmp_mat;
        get_parameter("selection_pedigree_F")->getMatrix(&tmp_mat);
        if(tmp_mat.getNbCols() != 5) {
          error("\"selection_pedigree_F\" must be an array of size 5.\n");
          return false;
        }
        for(unsigned int i = 0; i < 5; i++)
          _Fpedigree[i] = tmp_mat.get(0,i);
      }

      for(unsigned int i = 0; i < 5; i++)
        _FitnessFixModel[i] = _base_fitness * exp( -_letheq * _Fpedigree[i] );

      _getRawFitness = &LCE_Selection_base::getFitnessFixedEffect;

    } else if(sel_model == "direct") {
      //check trait dimensionality!!
      _getRawFitness = &LCE_Selection_base::getFitnessDirect;

    } else if(sel_model == "quadratic") {

      if(_selectTraitDimension == 1){

        if(!setSelectionMatrix()) return false; //this to set the selection variance params
        //now reset the fitness function:
        _getRawFitness = &LCE_Selection_base::getFitnessUnivariateQuadratic;

      } else {
        error("\"quadratic\" fitness model implemented for a single trait only.\n");
        return false;
      }
    } else if(sel_model == "gaussian") {

    	if(get_parameter("age_classes_selected")->isSet()){

        //cout<<"3"<<endl;

      	return (setSelectionMatrixAge());

      }	else cout<<"chevre"<<endl;

      return (setSelectionMatrix());

    } else {
      error("wrong selection model, must be either \"fix\", \"direct\", or \"gaussian\".\n");
      return false;
    }
  } //end_if isSet()

  return true;
}
//----------------------------------------------------------------------------------------
//setSelectionMatrixAge
//----------------------------------------------------------------------------------------
bool LCE_Selection_base::setSelectionMatrixAge(){

//this is called only if age_classes_selected is set

  TMatrix tmp_mat;
  unsigned int patchNbr = _popPtr->getPatchNbr();


  tmp_mat.reset(_selectTraitDimension, _selectTraitDimension);

  //we have to check for spatial variation in variance and covariances
 // _selectTraitDimension = (unsigned int)get_parameter_value("selection_trait_dimension"); (already set above)

  //matrix including the strength of selection for each age-class
  //afterward assumed to be identical in every patches
  //get_parameter("selection_age_variance")->getMatrix(_selection_age_variance);

  int nb_classes_selected = _age_classes_selected->getNbCols();
  int numVarInRow = _selection_age_variance->getNbCols();
  int numRows = _selection_age_variance->getNbRows();
  int numCovariance = (_selectTraitDimension-1)*_selectTraitDimension/2;
  int age_pos, cov_pos;
  
  if(numVarInRow > _selectTraitDimension) {
    if( (numVarInRow - _selectTraitDimension) != numCovariance){
      error("number of covariance terms in \"selection_age_variance\" does not match with the number of traits\n");
      return false;
    }
  }
  
  
  _selection_matrix_age.clear();

// !!! we have to consider that more than one trait is under selection and with different strength !!!

  for(int a = 0; a < nb_classes_selected; a++){

    _selection_matrix_age.push_back( vector<TMatrix*>() ); //reset selection matrix;

    age_pos = a % numRows; //user may pass only one row that will be repeated for all age classes under selection

    //set the selection matrix:
    
    cov_pos = _selectTraitDimension;
    
    for( int i = 0; i < _selectTraitDimension; i++) {
      
      if(numVarInRow < _selectTraitDimension) //repetition of a few values
        tmp_mat.set(i, i, _selection_age_variance->get(age_pos, i % numVarInRow));
      else
        tmp_mat.set(i, i, _selection_age_variance->get(age_pos, i));
      
      
      //check if we have covariance in input, i.e. more than _selectTraitDimension elements in a row
      //we already checked that the number of elements corresponds to the number of covariance terms
      
      if(numVarInRow > _selectTraitDimension) {
        
        //{v1,v2,v3,cov12,cov13,cov23}
        //{v1,v2,v3,v4,cov12,cov13,cov14,cov23,cov24,cov34}
        
        for( int j = i+1; j < _selectTraitDimension; j++) {
          tmp_mat.set(i, j, _selection_age_variance->get(age_pos, cov_pos) );
          tmp_mat.set(j, i, _selection_age_variance->get(age_pos, cov_pos++) );
        }

      } else {
        
        for( int j = i+1; j < _selectTraitDimension; j++) {
          tmp_mat.set(i, j, 0);
          tmp_mat.set(j, i, 0);
        }
        
      }
    }

    for( unsigned int p = 0; p < patchNbr; p++) {

      _selection_matrix_age[a].push_back(new TMatrix(tmp_mat));
    }
  }

  //Create gsl matrices//

  //selection on more than one trait
  if(_selectTraitDimension > 1) {

    _gsl_selection_matrix_age.clear();


    for(int i=0; i< nb_classes_selected; i++){

      //Inverting the selection matrices:

      //cout<<"test selection matrix"<<_selection_matrix_age[i][0]->get(0,0)<<endl;

      _selection_matrix_age[i][0]->inverse();

      _gsl_selection_matrix_age.push_back( vector<gsl_matrix*>() );

      for( unsigned int p = 0; p < patchNbr; p++) {

        _gsl_selection_matrix_age[i].push_back( gsl_matrix_alloc(_selectTraitDimension,
                                                                 _selectTraitDimension) );

        _selection_matrix_age[i][0]->get_gsl_matrix(_gsl_selection_matrix_age[i][p]);
      }


		}

		//allocate the vectors used by the fitness function:
		if(_diffs != NULL) gsl_vector_free(_diffs);
		_diffs = gsl_vector_alloc( _selectTraitDimension );

		if(_res1 != NULL) gsl_vector_free(_res1);
		_res1 = gsl_vector_alloc( _selectTraitDimension );

		if(_eVariance > 0)
			_getRawFitnessAge = &LCE_Selection_base::getFitnessMultivariateGaussianAge_VE;
		else
			_getRawFitnessAge = &LCE_Selection_base::getFitnessMultivariateGaussianAge;



  } else {

    //selection on one trait only



    if(_eVariance > 0){
      _getRawFitnessAge = &LCE_Selection_base::getFitnessUnivariateGaussianAge_VE;
    }
    else
    {
      _getRawFitnessAge = &LCE_Selection_base::getFitnessUnivariateGaussianAge;}

  }

  //cout<<"4"<<endl;

  if(!get_parameter("selection_local_optima")->isSet()) {

	  error("parameter \"selection_local_optima\" must be set to have Gaussian selection.\n");
	  return false;

  } else {

	  //cout<<"current gen "<<_popPtr->getCurrentGeneration()<<endl;

	  if(_popPtr->getCurrentGeneration() == 0){

		  if(_local_optima == 0) _local_optima = new TMatrix();

		  _local_optima->reset(patchNbr, _selectTraitDimension);

		  _paramSet->getMatrix("selection_local_optima", &tmp_mat);

		  return setSpatialMatrix("selection_local_optima", "\"selection_trait_dimension\"", &tmp_mat, _local_optima,
				  _selectTraitDimension, patchNbr, _paramSet->isSet("selection_randomize"));
	  }
  }

  return true;

}

// ----------------------------------------------------------------------------------------
// setSelectionMatrix
// ----------------------------------------------------------------------------------------
bool LCE_Selection_base::setSelectionMatrix()
{
  TMatrix tmp_mat;
  unsigned int patchNbr = _popPtr->getPatchNbr();

  //if(!get_parameter("selection_matrix")->isSet() && !get_parameter("selection_variance")->isSet()) {
    //error("\"selection_matrix\" or \"selection_variance\" must be set with selection model = \"gaussian\".\n");
    //return false;
  //}

  if(get_parameter("selection_variance")->isSet() && !get_parameter("selection_trait_dimension")->isSet()) {
    error("parameter \"selection_trait_dimension\" is missing!\n");
    return false;
  }

  _selection_matrix.clear();


  if(get_parameter("selection_matrix")->isSet()) {

    //selection matrix provided, same selection surface in each patch
    _paramSet->getMatrix("selection_matrix", &tmp_mat);

    if(tmp_mat.getNbCols() != tmp_mat.getNbRows()) {
      error("\"selection_matrix\" must be a square matrix!\n");
      return false;
    }

    _selectTraitDimension = tmp_mat.getNbCols();

    //we have one selection matrix per patch
    for(unsigned int i = 0; i < patchNbr; ++i)
      _selection_matrix.push_back( new TMatrix(tmp_mat) );

  } else {
    //we have to check for spatial variation in variance and covariances
    _selectTraitDimension = (unsigned int)get_parameter_value("selection_trait_dimension");

    TMatrix var_spatmat, corr_spatmat;

    //setting variance spatial matrix:
    var_spatmat.reset(patchNbr, (unsigned)_selectTraitDimension);

    if(get_parameter("selection_variance")->isMatrix()) {

      _paramSet->getMatrix("selection_variance", &tmp_mat);

      if( !setSpatialMatrix("selection_variance","\"selection_trait_dimension\"", &tmp_mat, &var_spatmat,
                            (unsigned)_selectTraitDimension, patchNbr, _paramSet->isSet("selection_randomize") ) )
        return false;

    } else {

      var_spatmat.assign(get_parameter_value("selection_variance"));
    }

    //setting correlation spatial matrix:
    corr_spatmat.reset(patchNbr, (unsigned)_selectTraitDimension*(_selectTraitDimension-1)/2);

    if(get_parameter("selection_correlation")->isMatrix()) {

      _paramSet->getMatrix("selection_correlation", &tmp_mat);

      if( !setSpatialMatrix("selection_correlation","the num of correlation coefficients", &tmp_mat, &corr_spatmat,
                            (unsigned)_selectTraitDimension*(_selectTraitDimension-1)/2, patchNbr, _paramSet->isSet("selection_randomize") ) )
        return false;

    } else {

      corr_spatmat.assign((get_parameter("selection_correlation")->isSet() ?
                           get_parameter_value("selection_correlation") : 0.0 ));
    }

    //set the selection matrix:
    tmp_mat.reset(_selectTraitDimension, _selectTraitDimension);
    double covar;
    unsigned int col;
    for( unsigned int p = 0; p < patchNbr; p++) {
      col = 0;
      for( int i = 0; i < _selectTraitDimension; i++) {
        tmp_mat.set(i, i, var_spatmat.get(p, i));
        for( int j = i+1; j < _selectTraitDimension; j++) {
          covar = corr_spatmat.get(p, col) * sqrt( var_spatmat.get(p, i) * var_spatmat.get(p, j) );
          tmp_mat.set(i, j, covar);
          tmp_mat.set(j, i, covar);
          col++;
        }
      }
      _selection_matrix.push_back(new TMatrix(tmp_mat));
    }
  }

  if(_selectTraitDimension > 1) {

    //selection on more than one trait

    //inversing the selection matrices:
    if(_gsl_selection_matrix.size() != 0)
      for(unsigned int i = 0; i < _gsl_selection_matrix.size(); ++i)
        if(_gsl_selection_matrix[i] != NULL) gsl_matrix_free( _gsl_selection_matrix[i] );
    _gsl_selection_matrix.clear();

    for( unsigned int p = 0; p < patchNbr; p++) {
      _selection_matrix[p]->inverse();

      _gsl_selection_matrix.push_back( gsl_matrix_alloc(_selectTraitDimension, _selectTraitDimension) );

      _selection_matrix[p]->get_gsl_matrix(_gsl_selection_matrix[p]);
    }
    //allocate the vectors used by the fitness function:
    if(_diffs != NULL) gsl_vector_free(_diffs);
    _diffs = gsl_vector_alloc( _selectTraitDimension );

    if(_res1 != NULL) gsl_vector_free(_res1);
    _res1 = gsl_vector_alloc( _selectTraitDimension );

    if(_eVariance > 0)
      _getRawFitness = &LCE_Selection_base::getFitnessMultivariateGaussian_VE;
    else
      _getRawFitness = &LCE_Selection_base::getFitnessMultivariateGaussian;

  } else {
    //selection on one trait only
    _selection_variance.clear();
    for(unsigned int p = 0; p < patchNbr; p++)
      _selection_variance.push_back( _selection_matrix[p]->get(0, 0) );

    if(_eVariance > 0)
      _getRawFitness = &LCE_Selection_base::getFitnessUnivariateGaussian_VE;
    else
      _getRawFitness = &LCE_Selection_base::getFitnessUnivariateGaussian;
  }

  //set the traits' local optima, _selectTraitDimension must be set before that (= nbr of traits to select on)


  if(!get_parameter("selection_local_optima")->isSet()) {

    error("parameter \"selection_local_optima\" must be set to have Gaussian selection.\n");
    return false;

  } else {

	  if(_local_optima == 0) _local_optima = new TMatrix();

    _local_optima->reset(patchNbr, _selectTraitDimension);

    _paramSet->getMatrix("selection_local_optima", &tmp_mat);

    return setSpatialMatrix("selection_local_optima", "\"selection_trait_dimension\"", &tmp_mat, _local_optima,
                            _selectTraitDimension, patchNbr, _paramSet->isSet("selection_randomize"));
  }
  return true;
}
// ----------------------------------------------------------------------------------------
// LCE_Selection_base::getFitnessUnivariateQuadratic
// ----------------------------------------------------------------------------------------
double LCE_Selection_base::getFitnessUnivariateQuadratic ( Individual* ind, unsigned int patch )
{
  register double res2, diff;

  _phe = (double*)ind->getTraitValue(_LCELinkedTraitIndex);

  diff = _phe[0] - _local_optima->get(patch, 0);

  res2 = diff*diff / _selection_variance[patch];

  return 1 - res2;
}
// ----------------------------------------------------------------------------------------
// LCE_Selection_base::getFitnessUnivariateGaussian
// ----------------------------------------------------------------------------------------
double LCE_Selection_base::getFitnessUnivariateGaussian ( Individual* ind, unsigned int patch )
{
  register double res2, diff;

  _phe = (double*)ind->getTraitValue(_LCELinkedTraitIndex);

  diff = _phe[0] - _local_optima->get(patch, 0);

  res2 = diff*diff / _selection_variance[patch];

  return exp( -0.5 * res2 );
}

// ----------------------------------------------------------------------------------------
// LCE_Selection_base::getFitnessUnivariateGaussianAge
// ----------------------------------------------------------------------------------------
double LCE_Selection_base::getFitnessUnivariateGaussianAge ( Individual* ind, unsigned int patch, unsigned int age_s )
{


  register double res2, diff;

  _phe = (double*)ind->getTraitValue(_LCELinkedTraitIndex);

  diff = _phe[0] - _local_optima->get(patch, 0);

  res2 = diff*diff / _selection_matrix_age[age_s][patch]->get(0,0);

  return exp( -0.5 * res2 );
}

// ----------------------------------------------------------------------------------------
// LCE_Selection_base::getFitnessMultivariateGaussianAge
// ----------------------------------------------------------------------------------------
double LCE_Selection_base::getFitnessMultivariateGaussianAge ( Individual* ind, unsigned int patch, unsigned age_s )
{

//cout<<"getmultivariategausssianAge"<<endl;
//cout<<"age_s"<<age_s<<endl;

  register double res2;

  if(!ind) fatal("passing NULL ind ptr to LCE_Selection_base::getFitnessMultivariateGaussian!!!\n");

  ///check///
//  cout<<"fitnessmultivariate_AGE"<<endl;
  //  cout<<"_LCELinkedTraitIndex ="<<_LCELinkedTraitIndex<<endl;

  _phe = (double*)ind->getTraitValue(_LCELinkedTraitIndex);

  for( int i = 0; i < _selectTraitDimension; i++)
    gsl_vector_set(_diffs, i, _phe[i] - _local_optima->get(patch, i));

  //(diff)T * W * diff:
  //right partial product:
  gsl_blas_dsymv(CblasUpper, 1.0, _gsl_selection_matrix_age[age_s][patch], _diffs, 0.0, _res1);
  //left product:
  gsl_blas_ddot(_diffs, _res1, &res2);

	//cout<< gsl_matrix_get (_gsl_selection_matrix_age[age_s][patch], 0, 0)<<endl;



	//cout<<exp( -0.5 * res2 )<<endl;


  return exp( -0.5 * res2 );
}


// ----------------------------------------------------------------------------------------
// LCE_Selection_base::getFitnessMultivariateGaussian
// ----------------------------------------------------------------------------------------
double LCE_Selection_base::getFitnessMultivariateGaussian ( Individual* ind, unsigned int patch )
{
  register double res2;

  if(!ind) fatal("passing NULL ind ptr to LCE_Selection_base::getFitnessMultivariateGaussian!!!\n");

  ///check///
  //cout<<"fitnessmultivariate"<<endl;
  //  cout<<"_LCELinkedTraitIndex ="<<_LCELinkedTraitIndex<<endl;

  _phe = (double*)ind->getTraitValue(_LCELinkedTraitIndex);

  for( int i = 0; i < _selectTraitDimension; i++)
    gsl_vector_set(_diffs, i, _phe[i] - _local_optima->get(patch, i));

  //(diff)T * W * diff:
  //right partial product:
  gsl_blas_dsymv(CblasUpper, 1.0, _gsl_selection_matrix[patch], _diffs, 0.0, _res1);
  //left product:
  gsl_blas_ddot(_diffs, _res1, &res2);

  return exp( -0.5 * res2 );
}
// ----------------------------------------------------------------------------------------
// LCE_Selection_base::getFitnessUnivariateGaussian_VE
// ----------------------------------------------------------------------------------------
double LCE_Selection_base::getFitnessUnivariateGaussian_VE ( Individual* ind, unsigned int patch )
{
  register double res2, diff;

  _phe = (double*)ind->getTraitValue(_LCELinkedTraitIndex);

  //add the environmental variance here:
  diff = _phe[0] + RAND::Gaussian(_eVariance) - _local_optima->get(patch, 0);

  res2 = diff*diff / _selection_variance[patch];

  return exp( -0.5 * res2 );
}
// ----------------------------------------------------------------------------------------
// LCE_Selection_base::getFitnessMultivariateGaussian_VE
// ----------------------------------------------------------------------------------------
double LCE_Selection_base::getFitnessMultivariateGaussian_VE ( Individual* ind, unsigned int patch )
{
  register double res2;

  _phe = (double*)ind->getTraitValue(_LCELinkedTraitIndex);

  //add the environmental variance here:
  for( int i = 0; i < _selectTraitDimension; i++)
    gsl_vector_set(_diffs, i, _phe[i]  + RAND::Gaussian(_eVariance) - _local_optima->get(patch, i));

  //(diff)T * W * diff:
  //right partial product:
  gsl_blas_dsymv(CblasUpper, 1.0, _gsl_selection_matrix[patch], _diffs, 0.0, _res1);
  //left product:
  gsl_blas_ddot(_diffs, _res1, &res2);

  return exp( -0.5 * res2 );
}

// ----------------------------------------------------------------------------------------
// LCE_Selection_base::getFitnessUnivariateGaussianAge_VE
// ----------------------------------------------------------------------------------------
double LCE_Selection_base::getFitnessUnivariateGaussianAge_VE ( Individual* ind, unsigned int patch, unsigned int age_s )
{
  register double res2, diff;


  //get individual age	to get the proper selection variance
//  	int age = ind->getAge(); UNUSED
  	//

  _phe = (double*)ind->getTraitValue(_LCELinkedTraitIndex);

  //add the environmental variance here:
  diff = _phe[0] + RAND::Gaussian(_eVariance) - _local_optima->get(patch, 0);

  res2 = diff*diff / _selection_matrix_age[age_s][patch]->get(0,0);

  return exp( -0.5 * res2 );
}
// ----------------------------------------------------------------------------------------
// LCE_Selection_base::getFitnessMultivariateGaussianAge_VE
// ----------------------------------------------------------------------------------------
double LCE_Selection_base::getFitnessMultivariateGaussianAge_VE ( Individual* ind, unsigned int patch, unsigned int age_s )
{
  register double res2;

  _phe = (double*)ind->getTraitValue(_LCELinkedTraitIndex);

//  cout<<"fitnessMultivariate_AGE-VE"<<endl;

  //add the environmental variance here:
  for( int i = 0; i < _selectTraitDimension; i++)
    gsl_vector_set(_diffs, i, _phe[i]  + RAND::Gaussian(_eVariance) - _local_optima->get(patch, i));

  //(diff)T * W * diff:
  //right partial product:
  gsl_blas_dsymv(CblasUpper, 1.0, _gsl_selection_matrix_age[age_s].at(patch), _diffs, 0.0, _res1);
  //left product:
  gsl_blas_ddot(_diffs, _res1, &res2);

  return exp( -0.5 * res2 );
}
// ----------------------------------------------------------------------------------------
// LCE_Selection_base::getMeanFitness
// ----------------------------------------------------------------------------------------
double LCE_Selection_base::getMeanFitness (age_idx age)
{
  double mean = 0;
  Patch *patch;
//  age_idx age = (AGE == ADULTS ? ADLTx : OFFSx);

  for(unsigned int i = 0, npatch = _popPtr->getPatchNbr(); i < npatch; i++) {
    patch = _popPtr->getPatch(i);
    for(unsigned int j = 0, size = patch->size(FEM, age); j < size; j++)
      mean += getFitness( patch->get(FEM, age, j), i, 0);
    for(unsigned int j = 0, size = patch->size(MAL, age); j < size; j++)
      mean += getFitness( patch->get(MAL, age, j), i, 0);
  }
  return mean/_popPtr->size(age);
}
// ----------------------------------------------------------------------------------------
// LCE_Selection_base::getMeanPatchFitness
// ----------------------------------------------------------------------------------------
double LCE_Selection_base::getMeanPatchFitness (age_idx age, unsigned int p)
{
  double mean = 0;
  Patch *patch = _popPtr->getPatch(p);

  for(unsigned int j = 0, size = patch->size(FEM, age); j < size; j++)
    mean += getFitness( patch->get(FEM, age, j), p, 0);

  for(unsigned int j = 0, size = patch->size(MAL, age); j < size; j++)
    mean += getFitness( patch->get(MAL, age, j), p, 0);

  return mean/patch->size(age);
}
// ----------------------------------------------------------------------------------------
// LCE_Selection_base::setScalingFactorLocal
// ----------------------------------------------------------------------------------------
void LCE_Selection_base::setScalingFactorLocal (age_idx age, unsigned int p)
{
  _scaling_factor = 1; //this to have the raw mean fitness below
  _scaling_factor = 1.0/getMeanPatchFitness(age, p);
}
// ----------------------------------------------------------------------------------------
// LCE_Selection_base::setScalingFactorGlobal
// ----------------------------------------------------------------------------------------
void LCE_Selection_base::setScalingFactorGlobal (age_idx age, unsigned int p)
{
  if(p != 0) return; //we compupte the mean fitness only once
  _scaling_factor = 1;
  _scaling_factor = 1.0/getMeanFitness(age);
}
// ----------------------------------------------------------------------------------------
// LCE_Selection_base::execute
// ----------------------------------------------------------------------------------------
void LCE_Selection_base::execute ()
{
  Patch * patch;
  TMatrix agec;
  unsigned int popSize = _popPtr->size(OFFSPRG);

  resetCounters();

  //if environment changes gradually (for now, each generation)
  //@TODO remove call to parameter interface, use a bool var in the class, faster
  if(get_parameter("rate_environmental_change")->isSet())
  { //if generation 1, reset the selection matrix
    unsigned int gen =_popPtr->getCurrentGeneration();

    if(gen == 1){ setSelectionMatrix(); }

    for(unsigned int p = 0; p < _popPtr->getPatchNbr(); p++) {//change for all patches
      for( int i = 0; i < _selectTraitDimension; i++){
        double change = _rchange->get(p,i);
        _local_optima ->plus(p,i,change);
      }
    }
  }

//    cout << "\nLCE_Selection_base::execute:: \n";

  for(unsigned int p = 0; p < _popPtr->getPatchNbr(); p++) {

    (this->*_setScalingFactor)(OFFSx, p);

//    cout << "--- patch "<<p<<" fitness scaling factor = "<<_scaling_factor<<endl;

    patch = _popPtr->getPatch(p);

    //@TODO remove call to parameter interface, use a bool var in the class, faster
    if(get_parameter("age_classes_selected")->isSet()){


    	for(int i =0 ; i<_age_classes_selected->getNbCols();i++) {

    		//cout<<"num age classes selected "<<_age_classes_selected->getNbCols()<<endl;

    		age_idx age = age_idx(_age_classes_selected->get(0,i));

    		//cout<<"age"<<age<<endl;

    	//@TODO add viability selection for MAL, check compatibility with all mating systems
    		doViabilitySelection(FEM, age, patch, p, i);

    		//cout<<"doViabilitySelection"<<endl;


    	}
    }

    else{

      doViabilitySelection(FEM, OFFSx, patch, p, OFFSx);//selection only on offspring only

      doViabilitySelection(MAL, OFFSx, patch, p, OFFSx);

    }

  }

  setMeans(popSize);
}
// ----------------------------------------------------------------------------------------
// LCE_Selection_base::doViabilitySelection
// ----------------------------------------------------------------------------------------
void LCE_Selection_base::doViabilitySelection (sex_t SEX, age_idx AGE, Patch* patch, unsigned int p ,unsigned int age_s)
{
  register Individual* ind;
  register double fitness;
  register unsigned int cat;

  //check
  //cout<<"selection"<<endl;
  //

//      cout<<"selection"<<endl;
  for(unsigned int i = 0; i < patch->size(SEX, AGE); i++) {

    ind = patch->get(SEX, AGE, i);

//    cout<<"ind "<<ind<<flush;
//    cout<<" (mother: "<<ind->getMother()<<flush<<", father: "<<ind->getFather()<<")\n";

    cat = ind->getPedigreeClass(ind->getMother(), ind->getFather()); //MOD SLIM NEMO 7.11

//    cout<<"pedigree category: "<<cat<<endl;

    fitness = getFitness( ind, p, age_s);

//    cout<<"fitness: "<<fitness<<endl;

  //check
//   cout << "    ind "<<i<<" raw fitness = "<<(this->*_getRawFitness)(ind, p)<<"; scaled = "<<fitness;
   //

    _fitness[cat] += fitness;
    _ind_cntr[cat]++;

    if(RAND::Uniform() > fitness ) {//individual selected or not

      patch->remove(SEX, AGE, i);

      _popPtr->recycle(ind);

      i--;

//      cout << " --> DIED\n";
    } //else; this individual stays in the patch
    else {
      _survival[cat]++;
//      cout << "\n";
    }

  }
}
// ----------------------------------------------------------------------------------------
// LCE_SelectionSH::setStatRecorders
// ----------------------------------------------------------------------------------------
bool LCE_SelectionSH::setStatRecorders (string& token)
{
  string age_tag = token.substr(0,token.find_first_of("."));
  string sub_token;
  age_t AGE = ALL;

  if (age_tag.size() != 0 && age_tag.size() != string::npos) {

    if (age_tag == "adlt") AGE = ADULTS;

    else if (age_tag == "off") AGE = OFFSPRG;

    else age_tag = "";

  } else {
    age_tag = "";
  }

  if (age_tag.size() != 0)
    sub_token = token.substr(token.find_first_of(".") + 1, string::npos);
  else
    sub_token = token;

  if(sub_token.compare("fitness") == 0) {
    add("Mean population fitness","fitness.mean",ALL,0,0,&LCE_SelectionSH::getMeanFitness,0,0,0);
    add("Mean population fitness","fitness.outb",ALL,0,0,0,&LCE_SelectionSH::getFitness,0,0);
    add("Mean population fitness","fitness.outw",ALL,1,0,0,&LCE_SelectionSH::getFitness,0,0);
    add("Mean population fitness","fitness.hsib",ALL,2,0,0,&LCE_SelectionSH::getFitness,0,0);
    add("Mean population fitness","fitness.fsib",ALL,3,0,0,&LCE_SelectionSH::getFitness,0,0);
    add("Mean population fitness","fitness.self",ALL,4,0,0,&LCE_SelectionSH::getFitness,0,0);
  } else if(sub_token.compare("survival") == 0) {
    add("Mean offspring survival","survival.outb",ALL,0,0,0,&LCE_SelectionSH::getSurvival,0,0);
    add("Mean offspring survival","survival.outw",ALL,1,0,0,&LCE_SelectionSH::getSurvival,0,0);
    add("Mean offspring survival","survival.hsib",ALL,2,0,0,&LCE_SelectionSH::getSurvival,0,0);
    add("Mean offspring survival","survival.fsib",ALL,3,0,0,&LCE_SelectionSH::getSurvival,0,0);
    add("Mean offspring survival","survival.self",ALL,4,0,0,&LCE_SelectionSH::getSurvival,0,0);
  } else if(sub_token.compare("fitness.prop") == 0) {
    add("Proportion of b/n demes outbreds","prop.outb",ALL,0,0,0,&LCE_SelectionSH::getPedProp,0,0);
    add("Proportion of w/n demes outbreds","prop.outw",ALL,1,0,0,&LCE_SelectionSH::getPedProp,0,0);
    add("Proportion of half-sib crossings","prop.hsib",ALL,2,0,0,&LCE_SelectionSH::getPedProp,0,0);
    add("Proportion of full-sib crossings","prop.fsib",ALL,3,0,0,&LCE_SelectionSH::getPedProp,0,0);
    add("Proportion of selfed progeny","prop.self",ALL,4,0,0,&LCE_SelectionSH::getPedProp,0,0);

  } else if(sub_token.compare("fitness.patch") == 0) {

    addMeanPerPatch(AGE);

  } else if(sub_token.compare("fitness.var.patch") == 0) {

    addVarPerPatch(AGE);

  } else return false;

  return true;
}
// ----------------------------------------------------------------------------------------
// LCE_SelectionSH::addMeanPerPatch
// ----------------------------------------------------------------------------------------
void LCE_SelectionSH::addMeanPerPatch (age_t AGE)
{
  unsigned int patchNbr = _pop->getPatchNbr();

  if (AGE == ALL) {
    addMeanPerPatch(ADULTS);
    addMeanPerPatch(OFFSPRG);
    return;
  }

  string suffix = (AGE == ADULTS ? "adlt.":"off."); //at this stage, AGE != ALL
  string name = suffix + "W.avg.p";
  string long_name = "Mean fitness in patch ";
  string patch;

  void (LCE_SelectionSH::* setter) (void) = (AGE == ADULTS ?
                                             &LCE_SelectionSH::setAdultTable :
                                             &LCE_SelectionSH::setOffsprgTable);

//  unsigned int int_agex = static_cast<age_idx> ((AGE == ADULTS ? ADLTx : OFFSx));

  //first patch, gets the data table setter:
  add(long_name + "1", name + "1", AGE, 0, AGE,
      0,0,&LCE_SelectionSH::getMeanPatchFitness, setter);

  for(unsigned int p = 1; p < patchNbr; p++) {
    patch = tstring::int2str(p+1);
    add(long_name + patch, name + patch, AGE, p, AGE,
        0,0,&LCE_SelectionSH::getMeanPatchFitness,0);
  }
}
// ----------------------------------------------------------------------------------------
// LCE_SelectionSH::addVarPerPatch
// ----------------------------------------------------------------------------------------
void LCE_SelectionSH::addVarPerPatch (age_t AGE)
{
  unsigned int patchNbr = _pop->getPatchNbr();

  if (AGE == ALL) {
    addVarPerPatch(ADULTS);
    addVarPerPatch(OFFSPRG);
    return;
  }

  string suffix = (AGE == ADULTS ? "adlt.":"off."); //at this stage, AGE != ALL
  string name = suffix + "W.var.p";
  string long_name = "Var fitness in patch ";
  string patch;

  void (LCE_SelectionSH::* setter) (void) = (AGE == ADULTS ?
                                             &LCE_SelectionSH::setAdultTable :
                                             &LCE_SelectionSH::setOffsprgTable);

//  unsigned int int_agex = static_cast<age_idx> ((AGE == ADULTS ? ADLTx : OFFSx));

  //first patch, gets the data table setter:
  add(long_name + "1", name + "1", AGE, 0, AGE,
      0,0,&LCE_SelectionSH::getVarPatchFitness, setter);

  for(unsigned int p = 1; p < patchNbr; p++) {
    patch = tstring::int2str(p+1);
    add(long_name + patch, name + patch, AGE, p, AGE,
        0,0,&LCE_SelectionSH::getVarPatchFitness,0);
  }
}
// ----------------------------------------------------------------------------------------
// setDataTable
// ----------------------------------------------------------------------------------------
void LCE_SelectionSH::setDataTable(age_t AGE)
{
  if(_table_set_age == AGE
     && _table_set_gen == _pop->getCurrentGeneration()
     && _table_set_repl == _pop->getCurrentReplicate()
     ) return;

  unsigned int patchNbr = _pop->getPatchNbr();

  if(_phenoTable.size() != patchNbr) {
    if(_phenoTable.size() < patchNbr) {
      while (_phenoTable.size() < patchNbr)
        _phenoTable.push_back(vector<double>());
    } else {
      while (_phenoTable.size() > patchNbr) {
        _phenoTable.pop_back();
      }
    }
  }

  for(unsigned int i = 0; i < patchNbr; ++i) {
      _phenoTable[i].assign(_pop->size(AGE, i),0);
  }

  Patch* patch;

  age_idx age = (AGE == ADULTS ? ADLTx : OFFSx);

  for(unsigned int i = 0, n; i < patchNbr; i++) {


    patch = _pop->getPatch(i);


    if( !_SHLinkedEvent->_is_absolute ) {
      if(_SHLinkedEvent->_is_local)
        _SHLinkedEvent->setScalingFactorLocal(age, i);
      else
        _SHLinkedEvent->setScalingFactorGlobal(age, i);
    }

//    if ((patch->size(MAL, age)+patch->size(FEM, age)) != _phenoTable[i].size()) {
//      fatal("problem while recording quanti trait values; table size doesn't match patch size.\n");
//    }


    n = 0;

    for(unsigned int j = 0, size = patch->size(FEM, AGE);
        j < size && n < _phenoTable[i].size();
        j++)
    {

      _phenoTable[i][n++] = _SHLinkedEvent->getFitness( patch->get(FEM, AGE, j), i,0);

    }

    for(unsigned int j = 0, size = patch->size(MAL, AGE);
        j < size && n < _phenoTable[i].size();
        j++)
    {

      _phenoTable[i][n++] = _SHLinkedEvent->getFitness( patch->get(MAL, AGE, j), i,0);

    }


    if (n != _phenoTable[i].size()) {
      fatal("problem while recording fitness trait values; size counter doesn't match table size.\n");
    }
  }

  _table_set_age  = AGE;
  _table_set_gen  = _pop->getCurrentGeneration();
  _table_set_repl = _pop->getCurrentReplicate();

}
// ----------------------------------------------------------------------------------------
// LCE_SelectionSH::getMeanPatchFitness
// ----------------------------------------------------------------------------------------
double LCE_SelectionSH::getMeanPatchFitness (unsigned int i, age_t AGE)
{
//  age_idx age = static_cast<age_idx> (int_agex);
  unsigned int patch_size = _pop->getPatchPtr(i)->size(AGE);

  assert(patch_size == _phenoTable[i].size());

  if(patch_size == 0) return (nanf("NULL"));

  return getMeanPatchFitness(i);
}
// ----------------------------------------------------------------------------------------
// LCE_SelectionSH::getMeanPatchFitness
// ----------------------------------------------------------------------------------------
double LCE_SelectionSH::getMeanPatchFitness (unsigned int i)
{
  double mean = 0;

  unsigned int size = _phenoTable[i].size();

  for(unsigned int j = 0; j < size; j++)

    mean += _phenoTable[i][j];

  return mean/size;
}
// ----------------------------------------------------------------------------------------
// LCE_SelectionSH::getVarPatchFitness
// ----------------------------------------------------------------------------------------
double LCE_SelectionSH::getVarPatchFitness (unsigned int i, age_t AGE)
{
//  age_idx age = static_cast<age_idx> (int_agex);
  unsigned int patch_size = _pop->getPatchPtr(i)->size(AGE);

  assert(patch_size == _phenoTable[i].size());

  if(patch_size == 0) return nanf("NULL");

  double mean = getMeanPatchFitness(i);
  double var = 0;

  for(unsigned int j = 0; j < patch_size; j++) {

    var += pow(_phenoTable[i][j] - mean, 2.0);
  }

  return var/patch_size;
}
