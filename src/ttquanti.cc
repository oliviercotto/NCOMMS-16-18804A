/** $Id: ttquanti.cc,v 1.23.2.6 2016-02-09 14:16:01 fred Exp $
 *
 *  @file ttquanti.cc
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
 *  created on @date 14.11.2005
 * 
 *  @author fred
 */

#include <sstream>
#include <fstream>
#include <string.h>
#include <cmath>
#include <algorithm>
#include "ttquanti.h"
#include "filehandler.h"
#include "output.h"
#include "Uniform.h"
#include "tstring.h"
#include "utils.h"

#ifdef HAS_GSL
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#endif

void store_quanti_trait_values (Patch* patch, unsigned int patchID, unsigned int size, unsigned int *cntr,
                                sex_t SEX, age_t AGE, DataTable<double> *ptable, DataTable<double> *gtable,
                                unsigned int nTrait, unsigned int TraitIndex);


// ------------------------------------------------------------------------------

//                             TProtoQuanti

// ----------------------------------------------------------------------------------------
// cstor
// ----------------------------------------------------------------------------------------
TProtoQuanti::TProtoQuanti() :
_nb_locus(0),
_nb_traits(0),
_seq_length(0),
_allele_model(0),
_allele_value(0),
_mutation_matrix(0),
_gsl_mutation_matrix(0),
_evect(0),
_eval(0),
_effects_multivar(0),
_ws(0),
_genomic_mutation_rate(0),
_mutation_correlation(0),
_mutation_sigma(0),
_init_value(0),
_all_chooser(0),
_sizeofLocusType(sizeof(double)),
_eVariance(0),
_mutation_func_ptr(0),
_stats(0),
_writer(0),
_freqExtractor(0)
{
  set_paramset("quantitative_trait", false, this);
  
  add_parameter("quanti_traits",INT,true,false,0,0);
  add_parameter("quanti_loci",INT,true,false,0,0);
  add_parameter("quanti_allele_model",STR,false,false,0,0);
  add_parameter("quanti_allele_value",DBL,false,false,0,0);
  add_parameter("quanti_init_value",MAT,false,false,0,0);
  add_parameter("quanti_init_model",INT,false,true,0,3);
  add_parameter("quanti_environmental_variance",MAT,false,false,0,0);
  
  add_parameter("quanti_mutation_rate",DBL,true,true,0,1, 0);
  add_parameter("quanti_mutation_variance",DBL,false,false,0,0, 0);
  add_parameter("quanti_mutation_correlation",DBL,false,false,0,0, 0);
  add_parameter("quanti_mutation_covariance",DBL,false,false,0,0, 0);
  add_parameter("quanti_mutation_matrix",MAT,false,false,0,0, 0);
  
  //genetic map parameters:
  TTProtoWithMap::addGeneticMapParameters("quanti");
  
  add_parameter("quanti_output",STR,false,false,0,0);
  add_parameter("quanti_logtime",INT,false,false,0,0);
  add_parameter("quanti_dir",STR,false,false,0,0);
  
  add_parameter("quanti_extract_freq",BOOL,false,false,0,0);
  add_parameter("quanti_freq_logtime",INT,false,false,0,0);
  add_parameter("quanti_freq_grain",DBL,false,false,0,0);
}
// ----------------------------------------------------------------------------------------
// copy cstor
// ----------------------------------------------------------------------------------------
TProtoQuanti::TProtoQuanti(const TProtoQuanti& T) : 
_nb_locus(T._nb_locus),
_nb_traits(T._nb_traits),
_seq_length(T._seq_length),
_mutation_matrix(0),
_gsl_mutation_matrix(0),
_evect(0),
_eval(0),
_effects_multivar(0),
_ws(0),
_genomic_mutation_rate(T._genomic_mutation_rate),
_mutation_correlation(T._mutation_correlation),
_mutation_sigma(0),
_init_value(0),
_allele_value(0),
_all_chooser(0),
_eVariance(0),
_mutation_func_ptr(0),
_stats(0),
_writer(0),
_freqExtractor(0)
{ 
  _locusByteSize = T._nb_traits * sizeof(double);
  _sizeofLocusType = sizeof(double);
  _paramSet = new ParamSet( *(T._paramSet) ) ;
}
// ----------------------------------------------------------------------------------------
// dstor
// ----------------------------------------------------------------------------------------
TProtoQuanti::~TProtoQuanti ()
{
  reset_mutation_pointers();
  
  if(_stats != NULL){delete _stats; _stats = NULL;}
  if(_writer != NULL){delete _writer; _writer = NULL;}
  if(_freqExtractor != NULL){delete _freqExtractor; _freqExtractor = NULL;}
  if(_all_chooser) {delete [] _all_chooser; _all_chooser = NULL;}
  if(_init_value) {delete[] _init_value;_init_value = NULL;}
}
// ----------------------------------------------------------------------------------------
// setParameters
// ----------------------------------------------------------------------------------------
bool TProtoQuanti::setParameters()
{  
  if(get_parameter("quanti_mutation_covariance")->isSet()) {
    fatal("\"quanti_mutation_covariance\" is deprecated, use \"quanti_mutation_correlation\" instead.\n");
    return false;
  }
  _nb_traits = (unsigned int)get_parameter_value("quanti_traits");
  _nb_locus = (unsigned int)get_parameter_value("quanti_loci");
  //total nbre of values for both traits in a haploid genome
  _seq_length = _nb_traits * _nb_locus; //mutations are pleiotropic!!!
  _locusByteSize = _nb_traits * sizeof(double);
  _sizeofLocusType = sizeof(double);
  
  TMatrix tmp_evmat;

  if(get_parameter("quanti_environmental_variance")->isSet()){

	  _eVariance.clear();

	  //cout<<"check quanti env variance is set"<<endl;

	  get_parameter("quanti_environmental_variance")->getMatrix(&tmp_evmat);

	  for(unsigned int i = 0; i < _nb_traits; i++) _eVariance.push_back (tmp_evmat.get(0,i));

	   //cout<<_eVariance[0]<<" "<<_eVariance[1]<<" "<<_eVariance[2]<<endl;
  }

  
  TMatrix tmp_matx;
  
  //---------------------------------------------------------------------------------------------
  //initial values:   
  if(_init_value != NULL) {
    delete [] _init_value;
    _init_value = NULL;
  }  
  
  if(get_parameter("quanti_init_value")->isSet()) {
    
    get_parameter("quanti_init_value")->getMatrix(&tmp_matx);
    
    if(tmp_matx.getNbRows() != 1 || tmp_matx.getNbCols() != _nb_traits) {
      fatal("\"quanti_init_value\" must be a vector of length equal to the number of traits!\n");
      return false;
    }
    
    _init_value = new double [_nb_traits];
    
    for(unsigned int i = 0; i < _nb_traits; i++)
      _init_value[i] = tmp_matx.get(0,i);
    
  }
  else {
    
    _init_value = new double [_nb_traits];
    
    for(unsigned int i = 0; i < _nb_traits; i++)
      _init_value[i] = 0.0;
  }
  
  if(get_parameter("quanti_init_model")->isSet()) {
    
    _doInitMutation = (unsigned int)get_parameter_value("quanti_init_model");
    
  } else {
    _doInitMutation = 1;
  }
  
  //---------------------------------------------------------------------------------------
  if( !setMutationParameters() ) return false; //sets _allele_model
  
  if( !setGeneticMapParameters("quanti") ) return false;
  
  //for use in the inherit_free function:
  if(_all_chooser != NULL) delete [] _all_chooser; _all_chooser = new bool[_nb_locus];
  
  //---------------------------------------------------------------------------------------
  if(_nb_traits == 1) {
    
    if (_allele_model == 1 || _allele_model == 2) {
      _mutation_func_ptr = &TProtoQuanti::getMutationEffectUnivariateDiallelic;
    } else {
      _mutation_func_ptr = &TProtoQuanti::getMutationEffectUnivariateGaussian;
    }
    
  } else if (_nb_traits == 2) {
    
    if (_allele_model == 1 || _allele_model == 2) {
      _mutation_func_ptr = &TProtoQuanti::getMutationEffectBivariateDiallelic;
    } else {
      _mutation_func_ptr = &TProtoQuanti::getMutationEffectBivariateGaussian;
    }
    
  } else if (_nb_traits > 2) {
    
    if (_allele_model > 2) {
      _mutation_func_ptr = &TProtoQuanti::getMutationEffectMultivariateGaussian;
    } else {
      fatal("in \"quanti\" trait, the di-allelic model is only allowed for max. 2 quantitative traits.");
    }
    
  }
  
  return true; 
}
// ----------------------------------------------------------------------------------------
// setMutationParameters
// ----------------------------------------------------------------------------------------
bool TProtoQuanti::setMutationParameters ()
{  
  _genomic_mutation_rate = get_parameter_value("quanti_mutation_rate") * 2 * _nb_locus;  
  
  if(get_parameter("quanti_mutation_correlation")->isSet())
    _mutation_correlation = get_parameter_value("quanti_mutation_correlation");
  else
    _mutation_correlation = 0;
  
  reset_mutation_pointers();
  
  //checking allelic model
  if (get_parameter("quanti_allele_model")->isSet()) {
    
    string model = get_parameter("quanti_allele_model")->getArg();
    
    if (model == "diallelic") {
      
      _allele_model = 1;
      
      return setDiallelicMutationModel ();
      
    } else if (model == "diallelic_HC") {
      
      _allele_model = 2;
      
      return setDiallelicMutationModel ();
      
    } else if (model == "continuous") {
      
      _allele_model = 3;
      
      return setContinuousMutationModel ();
      
    } else if (model == "continuous_HC") {
      
      _allele_model = 4;
      
      return setContinuousMutationModel ();
      
    } else {
      error("\"quanti_allele_model\" has options \"diallelic[_HC]\" or \"continuous[_HC]\" only. \n");
      return false;
    }
    
  } 
  else { //default model
    _allele_model = 3;
    return setContinuousMutationModel ();
  }
  
  return true;
}
//---------------------------------------------------------------------------------------------
// setDiallelicMutationModel
//---------------------------------------------------------------------------------------------
bool TProtoQuanti::setDiallelicMutationModel ()
{
  if (!get_parameter("quanti_allele_value")->isSet()) {
    error("in \"quanti\" trait, \"quanti_allele_value\" is not set for the diallelic model.\n");
    return false;
    
  } else {
    
    assert(_allele_value == NULL); 
    //should be the case as reset_mutation_pointers has been called in setMutationParameters
    
    
    _allele_value = new double* [_nb_locus];
    
    for (unsigned i = 0; i < _nb_locus; ++i) _allele_value[i] = new double [2];
    
    if (get_parameter("quanti_allele_value")->isMatrix()) { //locus-specific allelic values
      
      TMatrix tmp;
      
      get_parameter("quanti_allele_value")->getMatrix(&tmp);
      
      if (tmp.ncols() != _nb_locus) {
        fatal("\"quanti_allele_value\" must get a matrix with num. columns = num. loci.\n");
        return false;
      }
      
      if (tmp.nrows() == 1) {
        
        for (unsigned i = 0; i < _nb_locus; ++i) {
          _allele_value[i][0] = tmp.get(0, i);
          _allele_value[i][1] = -_allele_value[i][0];
        }
        
      } else if (tmp.nrows() == 2) {
        
        for (unsigned i = 0; i < _nb_locus; ++i) {
          _allele_value[i][0] = tmp.get(0, i);
          _allele_value[i][1] = tmp.get(1, i);
        }
        
      } else {
        fatal("\"quanti_allele_value\" must get a matrix with a max. of 2 rows (and num. columns = num. loci).\n");
        return false;
      }
      
      
    } else { //no locus-specific allelic values
      
      double val = get_parameter_value("quanti_allele_value");
      for (unsigned i = 0; i < _nb_locus; ++i) {
        _allele_value[i][0] = val;
        _allele_value[i][1] = -_allele_value[i][0];
      }
      
    }
    
  }
  
  return true;
}
//---------------------------------------------------------------------------------------------
// setContinuousMutationModel
//---------------------------------------------------------------------------------------------
bool TProtoQuanti::setContinuousMutationModel ()
{
  unsigned int dims[2];
  //setting the mutation variance-covariance matrix  
  if(get_parameter("quanti_mutation_variance")->isSet()) {
    
    if(get_parameter("quanti_mutation_matrix")->isSet()) {
      warning("both \"quanti_mutation_variance\" and \"quanti_mutation_matrix\" are set, using the matrix only!\n");
    } else {
      
      _mutation_sigma = new double [_nb_traits];
      
      double sigma = sqrt( get_parameter_value("quanti_mutation_variance") );
      
      for(unsigned int i = 0; i < _nb_traits; i++)
        _mutation_sigma[i] = sigma;
      
      if(_nb_traits > 1) {
        //setting the mutation matrix
        _mutation_matrix = new TMatrix();
        
        _gsl_mutation_matrix = gsl_matrix_alloc(_nb_traits, _nb_traits);
        
        double covar, var = get_parameter_value("quanti_mutation_variance");
        
        covar = _mutation_correlation * var;
        
        for(unsigned int i = 0; i < _nb_traits; i++)
          gsl_matrix_set(_gsl_mutation_matrix, i, i, var);
        
        for(unsigned int i = 0; i < _nb_traits - 1; i++)  
          for(unsigned int j = i + 1 ; j < _nb_traits; j++) {
            gsl_matrix_set(_gsl_mutation_matrix, i, j, covar);
            gsl_matrix_set(_gsl_mutation_matrix, j, i, covar);
          }
        _mutation_matrix->set_from_gsl_matrix(_gsl_mutation_matrix);
      }
    }
  }
  
  if(get_parameter("quanti_mutation_matrix")->isSet()) {
    
    _mutation_matrix = new TMatrix();
    
    get_parameter("quanti_mutation_matrix")->getMatrix(_mutation_matrix);
    
    _mutation_matrix->get_dims(&dims[0]);
    
    if( dims[0] != _nb_traits || dims[1] != _nb_traits) {
      error("\"quanti_mutation_matrix\" must be a square matrix of size = \"quanti_traits\" x \"quanti_traits\"\n");
      return false;
    }
    
    _gsl_mutation_matrix = gsl_matrix_alloc(dims[0], dims[1]);
    
    _mutation_matrix->get_gsl_matrix(_gsl_mutation_matrix);
    
    _mutation_sigma = new double [_nb_traits];
    
    for(unsigned int i = 0; i < _nb_traits; i++) 
      _mutation_sigma[i] = sqrt(_mutation_matrix->get(i, i));
    
  }
  else if(!get_parameter("quanti_mutation_variance")->isSet()) {
    error("\"quanti_mutation_matrix\" or \"quanti_mutation_variance\" must be specified!\n");
    return false;
  }
  //some more initializations
  if(_mutation_matrix) {
    
    _evect = gsl_matrix_alloc(_nb_traits, _nb_traits);
    _eval = gsl_vector_alloc(_nb_traits);
    
    set_mutation_matrix_decomposition();
    
    if(_nb_traits == 2)
      _mutation_correlation = _mutation_matrix->get(0,1) / sqrt(_mutation_matrix->get(0,0)*_mutation_matrix->get(1,1));
    
    else if(_nb_traits > 2) {
      _effects_multivar = gsl_vector_alloc(_nb_traits);
      _ws = gsl_vector_alloc(_nb_traits);
    }
#ifdef _DEBUG_
    message("-- Mutation matrix:\n");
    _mutation_matrix->show_up();
    message("-- MMatrix decomposition:\n");
    for(unsigned int i = 0; i < _nb_traits; i++)
      cout<<gsl_vector_get(_eval,i)<<" ";
    cout<<endl;
    if(_nb_traits == 2) message("-- mutation correlation: %f\n",_mutation_correlation);
#endif
  }
  
  return true;
}
// ----------------------------------------------------------------------------------------
// set_mutation_matrix_decomposition
// ----------------------------------------------------------------------------------------
void TProtoQuanti::set_mutation_matrix_decomposition()
{
  gsl_matrix *E = gsl_matrix_alloc(_nb_traits,_nb_traits);
  gsl_matrix_memcpy(E,_gsl_mutation_matrix);
  gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (_nb_traits);
  gsl_eigen_symmv (E, _eval, _evect, w);
  gsl_eigen_symmv_free (w);
  gsl_matrix_free(E);
#ifdef _DEBUG_
  message("-- Mutation matrix eigenvalues:\n");
  for(unsigned int i = 0; i < _nb_traits; i++)
    cout<<gsl_vector_get(_eval,i)<<" ";
  cout<<endl;
#endif
  double eval;
  //take square root of eigenvalues, will be used in Gaussian as stdev
  for(unsigned int i = 0; i < _nb_traits; i++) {
    eval = gsl_vector_get(_eval,i);
    eval = (eval < 0.000001 ? 0 : eval);
    gsl_vector_set( _eval, i, sqrt(eval) );
  }
}
// ----------------------------------------------------------------------------------------
// getMutationEffectMultivariateGaussian
// ----------------------------------------------------------------------------------------
double* TProtoQuanti::getMutationEffectMultivariateGaussian (unsigned int loc)
{
  RAND::MultivariateGaussian(_eval, _evect, _ws, _effects_multivar);
  return _effects_multivar->data;
}
// ----------------------------------------------------------------------------------------
// getMutationEffectBivariateGaussian
// ----------------------------------------------------------------------------------------
double* TProtoQuanti::getMutationEffectBivariateGaussian   (unsigned int loc)
{
  RAND::BivariateGaussian(_mutation_sigma[0], _mutation_sigma[1], _mutation_correlation,
                          &_effects_bivar[0], &_effects_bivar[1]);
  return &_effects_bivar[0];
}
// ----------------------------------------------------------------------------------------
// getMutationEffectUnivariateGaussian
// ----------------------------------------------------------------------------------------
double* TProtoQuanti::getMutationEffectUnivariateGaussian   (unsigned int loc)
{
  _effects_bivar[0] = RAND::Gaussian(_mutation_sigma[0]);
  return &_effects_bivar[0];
}
// ----------------------------------------------------------------------------------------
// getMutationEffectUnivariateDiallelic
// ----------------------------------------------------------------------------------------
double* TProtoQuanti::getMutationEffectUnivariateDiallelic   (unsigned int loc)
{
  _effects_bivar[0] = _allele_value[loc][ RAND::RandBool() ];
  return &_effects_bivar[0];
}
// ----------------------------------------------------------------------------------------
// getMutationEffectBivariateDiallelic
// ----------------------------------------------------------------------------------------
double* TProtoQuanti::getMutationEffectBivariateDiallelic   (unsigned int loc)
{
  bool pos = RAND::RandBool();
  _effects_bivar[0] = _allele_value[loc][pos];
  _effects_bivar[1] = _allele_value[loc][ (RAND::Uniform() < _mutation_correlation ? 
                                           pos : RAND::RandBool()) ];
  return &_effects_bivar[0];
}
// ----------------------------------------------------------------------------------------
// reset_mutation_pointers
// ----------------------------------------------------------------------------------------
void TProtoQuanti::reset_mutation_pointers()
{
  if(_allele_value) {
    for(unsigned int i = 0; i < _nb_locus; ++i)
      delete [] _allele_value[i];
    delete [] _allele_value;
    _allele_value = NULL;
  }
  
  if(_mutation_matrix != NULL) {
    delete _mutation_matrix;
    _mutation_matrix = NULL;
  }
  
  if(_gsl_mutation_matrix != NULL) 
    gsl_matrix_free(_gsl_mutation_matrix);
  
  if(_mutation_sigma != NULL){
    delete [] _mutation_sigma;
    _mutation_sigma = NULL;
  }
  
  if(_evect != NULL) {
    gsl_matrix_free(_evect);
    _evect = NULL;
  }
  
  if(_eval != NULL) {
    gsl_vector_free(_eval);
    _eval = NULL;
  }
  
  if(_effects_multivar != NULL) {
    gsl_vector_free(_effects_multivar);
    _effects_multivar = NULL;
  }
  
  if(_ws) {
    gsl_vector_free(_ws);
    _ws = 0;
  }
  
}
// ----------------------------------------------------------------------------------------
// inherit
// ----------------------------------------------------------------------------------------
inline void TProtoQuanti::inherit_free (sex_t SEX, double* seq, double** parent)
{
  //  bool all_chooser[_nb_locus];
  register unsigned int bloc;
  
  for(unsigned int i = 0; i < _nb_locus; ++i)
    _all_chooser[i] = RAND::RandBool();
  
  for(unsigned int i = 0; i < _nb_locus; ++i) {
    
    bloc = i*_nb_traits;
    
    memcpy(&seq[bloc], &parent[ _all_chooser[i] ][bloc], _locusByteSize);
    
    //    for(unsigned int j = 0; j < _nb_traits; ++j) {
    //      seq[bloc] = parent[ all_chooser[i] ][bloc];
    //      bloc++;
    //    }
  }
}
// ----------------------------------------------------------------------------------------
// inherit
// ----------------------------------------------------------------------------------------
inline void TProtoQuanti::inherit_low (sex_t SEX, double* seq, double** parent)
{
  register unsigned int prevLoc = 0, chrm_bloc, prev_bloc, cpy_bloc;
  register bool flipper;
  
  vector< unsigned int >& recTable = _map.getRecLoci(SEX, _mapIndex);
  vector< bool > & firstRecPos = _map.getFirstRecPosition(SEX);
  
  unsigned int nbRec = recTable.size();
  
  //  cout << "TProtoQuanti::inherit; sex="<<SEX<<"; nb Rec = "<<nbRec;//<<endl;
  
  //  if (!nbRec) return;
  
  for(unsigned int c = 0, stride = 0, rec = 0; c < _numChromosome; ++c) {
    
    flipper = firstRecPos[c];
    
    chrm_bloc = stride + _numLociPerChrmsm[c];
    
    prevLoc = stride;
    
    //    cout<<"chrm "<<c<<" side="<<firstRecPos[c]<<endl;
    
    for(; recTable[rec] < chrm_bloc && rec < nbRec; rec++) {
      
      prev_bloc = prevLoc * _nb_traits;
      
      cpy_bloc = (recTable[rec] - prevLoc) * _nb_traits;
      
      //      cout<<"copy seq from "<<prevLoc<<"("<<prev_bloc<<") to "<<recTable[rec]
      //      <<"("<<(recTable[rec] - prevLoc)<<" loc) ("<<cpy_bloc*_locusByteSize<<"B) side "<<flipper<<endl;
      
      memcpy(&seq[prev_bloc], &parent[flipper][prev_bloc], cpy_bloc * _sizeofLocusType);
      
      prevLoc = recTable[rec];
      
      //      prev_loc = i;
      //      
      //      while(!_recomb_positions[++i] && i < chrm_bloc) ; //!!_recomb_pos should have nb_locus+1 elmnts for the last pass because of the ++i!!
      //      
      //      seq_bloc = (i - prev_loc) * _nb_traits;
      //      
      //      //all this means if rec_pos is true for loc pos 1, recombination occurs between loc 0 and 1...
      
      //      memcpy(&seq[prev_bloc], &parent[flipper][prev_bloc], seq_bloc * sizeof(double));
      
      flipper = !flipper;
    }
    prev_bloc = prevLoc * _nb_traits;
    cpy_bloc = (chrm_bloc - prevLoc) * _nb_traits;
    //    cout << "copy end of chrmsm from "<<prevLoc<<" to "<<chrm_bloc
    //         <<"("<<(chrm_bloc - prevLoc)<<" loc) on side "<<flipper<<endl;
    //copy what's left between the last x-over point and the end of the chrmsme
    memcpy(&seq[prev_bloc], &parent[flipper][prev_bloc], cpy_bloc * _sizeofLocusType);
    
    stride += _numLociPerChrmsm[c];
  }
  //  cout<<" done"<<endl;
}
// ----------------------------------------------------------------------------------------
// hatch
// ----------------------------------------------------------------------------------------
TTQuanti* TProtoQuanti::hatch()
{
  TTQuanti* kid = new TTQuanti();
  
  kid->set_proto(this);
  kid->set_nb_locus(_nb_locus);
  kid->set_nb_traits(_nb_traits);
  kid->set_seq_length(_seq_length);
  kid->set_genomic_mutation_rate(_genomic_mutation_rate);
  kid->set_init_value(_init_value, _doInitMutation);
  kid->set_mutation_fptr(_mutation_func_ptr, (_allele_model != 1 && _allele_model != 3));
  if(_recombRate == 0.5)
    kid->set_inherit_fptr(&TProtoQuanti::inherit_free);
  else
    kid->set_inherit_fptr(&TProtoQuanti::inherit_low);
  
 kid->set_eVariance(_eVariance);
  
  return kid;
}
// ----------------------------------------------------------------------------------------
// loadFileServices
// ----------------------------------------------------------------------------------------
void TProtoQuanti::loadFileServices  (FileServices* loader)
{ 
  int logtime = 0;
  //writer
  if(get_parameter("quanti_output")->isSet()) {
    
    if(_writer == NULL) _writer = new TTQuantiFH(this);
    
    _writer->setOutputOption(get_parameter("quanti_output")->getArg());
    
    Param* param = get_parameter("quanti_logtime");
    
    if(param->isMatrix()) {
      
      TMatrix temp;
      param->getMatrix(&temp);
      _writer->set_multi(true, true, 1, &temp, get_parameter("quanti_dir")->getArg());
      
    } else   //  rpl_per, gen_per, rpl_occ, gen_occ, rank (0), path, self-ref      
      _writer->set(true, true, 1, (param->isSet() ? (int)param->getValue() : 0),
                   0, get_parameter("quanti_dir")->getArg(),this);
    
    loader->attach(_writer);
    
  } else if(_writer != NULL) {
    delete _writer;
    _writer = NULL;
  }
  
  //freq extractor
  if(get_parameter("quanti_extract_freq")->isSet()) {
    
    if(_freqExtractor == NULL) _freqExtractor = new TTQFreqExtractor(this);
    
    Param* param = get_parameter("quanti_freq_logtime");
    
    logtime = (param->isSet() ? (int)param->getValue() : 0);
    
    _freqExtractor->set(true,(logtime != 0),1,logtime,0,get_parameter("quanti_dir")->getArg(),this);
    
    param = get_parameter("quanti_freq_grain");
    
    if( ! param->isSet() ) {
      warning(" parameter \"quanti_freq_grain\" is not set, using 0.1\n");
      _freqExtractor->set_granularity(0.1);
    } else 
      _freqExtractor->set_granularity(param->getValue());
    
    loader->attach(_freqExtractor);
    
  } else if(_freqExtractor != NULL) {
    delete _freqExtractor;
    _freqExtractor = NULL;
  }
  
}
// ----------------------------------------------------------------------------------------
// loadStatServices
// ----------------------------------------------------------------------------------------
void TProtoQuanti::loadStatServices  (StatServices* loader)
{
  //allocate the stat handler
  if(_stats == NULL)
    _stats = new TTQuantiSH(this);
  
  if(_stats != NULL) {
    loader->attach(_stats);
  }
}





// ------------------------------------------------------------------------------

//                             TTQuanti

// ----------------------------------------------------------------------------------------
// operator=
// ----------------------------------------------------------------------------------------
TTQuanti& TTQuanti::operator= (const TTrait& T)
{  
  const TTQuanti& TQ = dynamic_cast<const TTQuanti&>(T);
  
  if(this != &TQ) {
    
    _nb_locus = TQ._nb_locus;
    _nb_traits = TQ._nb_traits;
    _seq_length = TQ._seq_length;
    reset();
    init();
    memcpy(_sequence[0],TQ._sequence[0],_seq_length*sizeof(double));
    memcpy(_sequence[1],TQ._sequence[1],_seq_length*sizeof(double));
    set_value();
  }
  
  return *this;
}
// ----------------------------------------------------------------------------------------
// operator==
// ----------------------------------------------------------------------------------------
bool TTQuanti::operator== (const TTrait& T)
{ 
  if(this->get_type().compare(T.get_type()) != 0) return false;
  
  const TTQuanti& TQ = dynamic_cast<const TTQuanti&>(T);
  
  if(this != &TQ) {
    if(_nb_locus != TQ._nb_locus) return false;
    if(_nb_traits != TQ._nb_traits) return false;
  }
  return true;
}
// ----------------------------------------------------------------------------------------
// operator!=
// ----------------------------------------------------------------------------------------
bool TTQuanti::operator!= (const TTrait& T)
{
  if(!((*this) == T) )
    return true;
  else
    return false;
}
// ----------------------------------------------------------------------------------------
// set_init_value
// ----------------------------------------------------------------------------------------
void TTQuanti::set_init_value             (double* val, unsigned int doInit)
{  
  set_init_value(val);
  
  _doInitMutation = doInit;
}
// ----------------------------------------------------------------------------------------
// set_init_value
// ----------------------------------------------------------------------------------------
void TTQuanti::set_init_value             (double* val) 
{ 
  assert(_nb_traits != 0);
  
  if(_init_value) delete [] _init_value;
  
  _init_value = new double [_nb_traits];
  
  for(unsigned int i = 0; i < _nb_traits; ++i) 
    _init_value[i] = val[i];
}  
// ----------------------------------------------------------------------------------------
// init
// ----------------------------------------------------------------------------------------
inline void TTQuanti::init ()
{  
  _sequence = new double*[2];
  _sequence[0] = new double [_seq_length];
  _sequence[1] = new double [_seq_length];
  if(!_phenotypes) _phenotypes = new double [_nb_traits];
  
  if(!_init_value) {
    _init_value = new double [_nb_traits];
    for(unsigned int i = 0; i < _nb_traits; ++i) 
      _init_value[i] = _myProto->get_init_value(i);
  }
}
// ----------------------------------------------------------------------------------------
// init_sequence
// ----------------------------------------------------------------------------------------
inline void TTQuanti::init_sequence ()
{
  unsigned int pos;
  //options:
  //0: no variation, init value = (trait value)/(2*_nb_locus)
  //1: init value = (trait value)/(2*_nb_locus) + 1 mutation/locus
  //2: init value = (trait value)/(2*_nb_locus) + 1 mutation/locus+make parts ==> mean trait value doesn't change
  //3: init value = (trait value + random deviate N(0,sdVm))/(2*_nb_locus) + 1 mutation/locus+make parts
  
  //decide what initial value to use
  
  //Note: this is kept here as the initial values may have been set individually by LCE_quanti
  //      it wouldn't make sense then to store the init values in the prototype only
  double my_init[_nb_traits];

  if(_doInitMutation == 3) {
    
    double sdVm;
    
    for(unsigned int j = 0; j < _nb_traits; j++) {
      //      sdVm = sqrt(4*_nb_locus*_mut_rate*_myProto->get_trait_var(j)); //trait variance = 2Vm
      sdVm = 0.25;
      my_init[j] = (_init_value[j] + RAND::Gaussian(sdVm)) / (2*_nb_locus);
    }
    
  } else {
    
    for(unsigned int j = 0; j < _nb_traits; j++) my_init[j] = _init_value[j] / (2*_nb_locus);
    
  }
  
  if(_myProto->_allele_model < 3) { //for the di-allelic models

    for(unsigned int i = 0; i < _nb_locus; i++) {
      pos = i * _nb_traits;
      for(unsigned int j = 0; j < _nb_traits; j++) {
        _sequence[0][pos] = _myProto->_allele_value[i][0]; //set with the positive allele value
        _sequence[1][pos] = _myProto->_allele_value[i][0];
        pos++;
      }
    }

  } else {

    //set the allele values from the trait value
    for(unsigned int i = 0; i < _nb_locus; i++) {
      pos = i * _nb_traits;
      for(unsigned int j = 0; j < _nb_traits; j++) {
        _sequence[0][pos] = my_init[j];
        _sequence[1][pos] = my_init[j];
        pos++;
      }
    }
  }
  
  //add random effects to allele values
  if(_doInitMutation != 0) {
    
    double *mut1, *mut2;
    unsigned int L;
    
    for(unsigned int i = 0; i < _nb_locus; i++) {
      
      pos = i * _nb_traits;
      
      mut1 = (_myProto->*_getMutationValues)(i);
      mut2 = (_myProto->*_getMutationValues)(i);
      
      if(_myProto->_allele_model < 3) { //diallelic models

    	  for(unsigned int j = 0; j < _nb_traits; j++) {
    		  _sequence[0][pos] = mut1[j];
    		  _sequence[1][pos] = mut2[j];
    		  pos++;
    	  }

      } else {

    	  for(unsigned int j = 0; j < _nb_traits; j++) {
    		  _sequence[0][pos] += mut1[j];
    		  _sequence[1][pos] += mut2[j];
    		  pos++;
    	  }

    	  if(_doInitMutation > 1) { // the make-parts algorithm

    		  //select a random locus
    		  do{
    			  L = RAND::Uniform(_nb_locus);
    		  }while(L == i);

    		  pos = L * _nb_traits;
    		  //subtract the previous random deviates from that locus
    		  for(unsigned int j = 0; j < _nb_traits; j++) {
    			  _sequence[0][pos] -= mut1[j];
    			  _sequence[1][pos] -= mut2[j];
    			  pos++;
    		  }
    	  }
      }
    }
  } 

}
// ----------------------------------------------------------------------------------------
// reset
// ----------------------------------------------------------------------------------------
inline void TTQuanti::reset ()
{
  if(_sequence != NULL) {
    delete [] _sequence[0]; 
    delete [] _sequence[1];
    delete [] _sequence; 
    _sequence = NULL;
  }
  if(_phenotypes != NULL) delete [] _phenotypes;
  _phenotypes = NULL;
  
  if(_init_value) delete [] _init_value;
  _init_value = NULL;
}
// ----------------------------------------------------------------------------------------
// inherit
// ----------------------------------------------------------------------------------------
inline void TTQuanti::inherit (TTrait* mother, TTrait* father)
{
  double** mother_seq = (double**)mother->get_sequence();
  double** father_seq = (double**)father->get_sequence();
  
  (_myProto->* _inherit) (FEM, _sequence[FEM], mother_seq);
  
  (_myProto->* _inherit) (MAL, _sequence[MAL], father_seq);
}
// ----------------------------------------------------------------------------------------
// mutate_noHC
// ----------------------------------------------------------------------------------------
inline void TTQuanti::mutate_noHC ()
{
  unsigned int NbMut = (unsigned int)RAND::Poisson(_genomic_mutation_rate);
  unsigned int mut_locus, mut_all, pos;
  double *effects;
  
  while(NbMut != 0) {
    mut_locus = RAND::Uniform(_nb_locus);
    //    cout<<"TTQuanti::mutate:: new mutation at <"<<mut_locus<<">"<<flush;
    effects = (_myProto->*_getMutationValues)(mut_locus);
    //    cout<<", new effects "<<effects[0]<<" "<<effects[1]<<flush;
    mut_all = RAND::RandBool();
    //    message(" at all %i\n",mut_all);
    pos = mut_locus*_nb_traits;
    for(unsigned int i = 0; i < _nb_traits; i++)
      _sequence[mut_all][pos + i] += effects[i];
    
    NbMut--; 
  }
  
}
// ----------------------------------------------------------------------------------------
// mutate_HC
// ----------------------------------------------------------------------------------------
inline void TTQuanti::mutate_HC ()
{
  unsigned int NbMut = (unsigned int)RAND::Poisson(_genomic_mutation_rate);
  unsigned int mut_locus, mut_all, pos;
  double *effects;
  
  while(NbMut != 0) {
    mut_locus = RAND::Uniform(_nb_locus);
    //    cout<<"TTQuanti::mutate:: new mutation at <"<<mut_locus<<">"<<flush;
    effects = (_myProto->*_getMutationValues)(mut_locus);
    //    cout<<", new effects "<<effects[0]<<" "<<effects[1]<<flush;
    mut_all = RAND::RandBool();
    //    message(" at all %i\n",mut_all);
    pos = mut_locus*_nb_traits;
    for(unsigned int i = 0; i < _nb_traits; i++)
      _sequence[mut_all][pos + i] = effects[i];
    
    NbMut--; 
  }
} 
// ----------------------------------------------------------------------------------------
// set_value
// ----------------------------------------------------------------------------------------
inline void TTQuanti::set_value ()
{
  register unsigned int loc;
  
  for(unsigned int i = 0; i < _nb_traits; ++i) 
    _phenotypes[i] = 0;
  
  for(unsigned int j = 0; j < _nb_locus; ++j) {
    loc = j * _nb_traits;
    for(unsigned int i = 0; i < _nb_traits; ++i) {
      _phenotypes[i] += _sequence[0][loc] + _sequence[1][loc];
      loc++;
    }
  } 
  
  if(!_eVariance.empty()){

	  //cout<<"check1"<<endl;

	double var;

    for(unsigned int i = 0; i < _nb_traits; ++i) {

    	if(var != 0) var=0;

  	  	var=_eVariance[i];

  	  	//cout<<var<<endl;

      _phenotypes[i] += RAND::Gaussian(sqrt(var));
    }
  }
}

// ----------------------------------------------------------------------------------------
// get_genotype
// ----------------------------------------------------------------------------------------
double TTQuanti::get_genotype (unsigned int trait)
{
  unsigned int loc;
  double genotype = 0;
  
  for(unsigned int j = 0; j < _nb_locus; ++j) {
    loc = j * _nb_traits + trait;
    genotype += _sequence[0][loc] + _sequence[1][loc];
  } 
  return genotype;
}
// ----------------------------------------------------------------------------------------
// show_up
// ----------------------------------------------------------------------------------------
void TTQuanti::show_up  ()
{
  message("\
          Trait's type: QUANTI\n\
          traits: %i\n\
          loci: %i\n",_nb_traits,_nb_locus);
  
  for(unsigned int i = 0; i < _nb_traits; i++)
    message("phenotype %i: %f\n",i+1,_phenotypes[i]);
  
  message("_sequence: \n0:");
  unsigned int loc;
  for(unsigned int i = 0; i < _nb_traits; ++i) {
    message("\nt%i:",i+1);
    for(unsigned int j = 0; (j < _nb_locus); ++j) {
      loc = j * _nb_traits + i;
      message("%.3f,",_sequence[0][loc]);
    }
  }
  message("\n1:");
  for(unsigned int i = 0; i < _nb_traits; ++i) {
    message("\nt%i:",i+1);
    for(unsigned int j = 0; (j < _nb_locus); ++j) {
      loc = j * _nb_traits + i;
      message("%.3f,",_sequence[1][loc]);
    }
  }
  message("\n");
  
}  
/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

//                             ******** TTQuantiSH ********

// ----------------------------------------------------------------------------------------
// 
// ----------------------------------------------------------------------------------------
void TTQuantiSH::resetPtrs ( )
{  
  
  if(_G != NULL) gsl_matrix_free(_G);
  if(_eval != NULL) gsl_vector_free(_eval);
  if(_evec != NULL) gsl_matrix_free(_evec);
  if(_ws != NULL)  gsl_eigen_symmv_free (_ws);
  
  if(_meanP != NULL) delete [] _meanP;
  if(_meanG != NULL) delete [] _meanG;
  if(_Va != NULL) delete [] _Va;
  if(_Vb != NULL) delete [] _Vb;
  if(_Vp != NULL) delete [] _Vp;
  if(_covar != NULL) delete [] _covar;
  if(_eigval != NULL) delete [] _eigval;
  
  
  if(_eigvect) {
    for(unsigned int i=0; i < _nb_trait; ++i) delete [] _eigvect[i];
    delete [] _eigvect;
  }
  if(_pVa) {
    for(unsigned int i=0; i < _nb_trait; ++i) delete [] _pVa[i];
    delete [] _pVa;
  }
  if(_pVp) {
    for(unsigned int i=0; i < _nb_trait; ++i) delete [] _pVp[i];
    delete [] _pVp;
  }
  if(_pmeanP) {
    for(unsigned int i=0; i < _nb_trait; ++i) delete [] _pmeanP[i];
    delete [] _pmeanP;
  }
  if(_pmeanG) {
    for(unsigned int i=0; i < _nb_trait; ++i) delete [] _pmeanG[i];
    delete [] _pmeanG;
  }
  if(_peigval) {
    for(unsigned int i=0; i < _nb_trait; ++i) delete [] _peigval[i];
    delete [] _peigval;
  }
  
  if(_pcovar != NULL) {  
    for(unsigned int i = 0; i < _patchNbr; i++) delete [] _pcovar[i];
    delete [] _pcovar;
  }
  
  if(_peigvect != NULL) {
    for(unsigned int i = 0; i < _patchNbr; i++) delete [] _peigvect[i];
    delete []  _peigvect;
  }
  
}
// ----------------------------------------------------------------------------------------
// init
// ----------------------------------------------------------------------------------------
void TTQuantiSH::init()
{  
  StatHandlerBase::init();
  
  _eVar = (!_SHLinkedTrait->get_env_var().empty());
  
  resetPtrs(); //deallocate anything that was previously allocated
  
  if(_patchNbr != _pop->getPatchNbr())
    _patchNbr = _pop->getPatchNbr();
  
  if(_nb_trait != _SHLinkedTrait->get_nb_traits())
    _nb_trait = _SHLinkedTrait->get_nb_traits();
  
  _G = gsl_matrix_alloc(_nb_trait,_nb_trait);
  _eval = gsl_vector_alloc (_nb_trait);
  _evec = gsl_matrix_alloc (_nb_trait, _nb_trait);
  _ws = gsl_eigen_symmv_alloc (_nb_trait);
  
  _meanP = new double [_nb_trait];
  _meanG = new double [_nb_trait];
  _Va = new double [_nb_trait];
  _Vb = new double [_nb_trait];
  _Vp = new double [_nb_trait];
  _covar = new double [_nb_trait*(_nb_trait -1)/2];
  _eigval = new double [_nb_trait];
  
  _eigvect = new double* [_nb_trait];
  _pVa = new double* [_nb_trait];
  _pVp = new double* [_nb_trait];
  _pmeanP = new double* [_nb_trait];
  _pmeanG = new double* [_nb_trait];
  _peigval = new double* [_nb_trait];
  
  for(unsigned int i=0; i < _nb_trait; ++i) {
    
    _eigvect[i]  = new double [_nb_trait];
    
    _pVa[i] = new double [_patchNbr];
    _pVp[i] = new double [_patchNbr];
    _pmeanP[i] = new double [_patchNbr];
    _pmeanG[i] = new double [_patchNbr];
    _peigval[i] = new double [_patchNbr];
    
  }
  
  _peigvect = new double* [_patchNbr];
  _pcovar = new double* [_patchNbr];
  
  for(unsigned int i = 0; i < _patchNbr; i++) {
    
    _pcovar[i] = new double [_nb_trait*(_nb_trait - 1)/2];
    _peigvect[i] = new double [_nb_trait*_nb_trait];
  }
  
}
// ----------------------------------------------------------------------------------------
// setStatsRecorders
// ----------------------------------------------------------------------------------------
bool TTQuantiSH::setStatRecorders(std::string& token)
{
#ifdef _DEBUG_
  message("-TTQuantiSH::setStatRecorders ");
#endif   
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
  
  //!!! attention, y a un prob ici quand pas de prefix d'age!!!!
  
  if(sub_token == "quanti") {
    addQuanti(AGE);
  } else if(sub_token == "quanti.eigen") {  //based on Vb; among-patch (D) matrix
    addEigen(AGE);
  } else if(sub_token == "quanti.eigenvalues") {  //based on Vb; among-patch (D) matrix
    addEigenValues(AGE);
  } else if(sub_token == "quanti.eigenvect1") {  //based on Vb; among-patch (D) matrix
    addEigenVect1(AGE);
  } else if(sub_token == "quanti.patch") {
    addQuantiPerPatch(AGE);
  } else if(sub_token == "quanti.mean.patch") {
    addAvgPerPatch(AGE);
  } else if(sub_token == "quanti.var.patch") {
    addVarPerPatch(AGE);
  } else if(sub_token == "quanti.covar.patch") {
    addCovarPerPatch(AGE);
  } else if(sub_token == "quanti.eigen.patch") {
    addEigenPerPatch(AGE);
  } else if(sub_token == "quanti.eigenvalues.patch") {
    addEigenValuesPerPatch(AGE);
  } else if(sub_token == "quanti.eigenvect1.patch") {
    addEigenVect1PerPatch(AGE);
  } else if(sub_token == "quanti.skew.patch") {
    addSkewPerPatch(AGE);
  } else
    return false;
  
  return true;
}// ----------------------------------------------------------------------------------------
// addQuanti
// ----------------------------------------------------------------------------------------
void TTQuantiSH::addQuanti (age_t AGE)
{
  if (AGE == ALL) {
    addQuanti(ADULTS);
    addQuanti(OFFSPRG);
    return;
  }
  
  string suffix = (AGE == ADULTS ? "adlt.":"off.");
  string name = suffix + "q";
  string t1, t2;
  
  void (TTQuantiSH::* setter) (void) = (AGE == ADULTS ?
                                        &TTQuantiSH::setAdultStats : &TTQuantiSH::setOffsprgStats);
  
  add("",suffix + "q1",AGE,0,0,0,&TTQuantiSH::getMeanPhenot,0,setter);
  
  for(unsigned int i = 1; i < _nb_trait; i++) {
    t1 = tstring::int2str(i+1);
    add("", name + t1,AGE,i,0,0,&TTQuantiSH::getMeanPhenot,0,0);
  } 
  
  //Va
  for(unsigned int i = 0; i < _nb_trait; i++) {
    t1 = tstring::int2str(i+1);
    add("", name + t1 +".Va",AGE,i,0,0,&TTQuantiSH::getVa,0,0);
  }
  
  //Vb
  for(unsigned int i = 0; i < _nb_trait; i++) {
    t1 = tstring::int2str(i+1);
    add("", name + t1 +".Vb",AGE,i,0,0,&TTQuantiSH::getVb,0,0);
  }
  
  //Vp
  if(_eVar) {

	  //cout<<"check evar print"<<endl;

    for(unsigned int i = 0; i < _nb_trait; i++) {
      t1 = tstring::int2str(i+1);
      add("", name + t1 +".Vp",AGE,i,0,0,&TTQuantiSH::getVp,0,0);
    }
  }
  //Qst
  for(unsigned int i = 0; i < _nb_trait; i++) {
    t1 = tstring::int2str(i+1);
    add("", name + t1 +".Qst",AGE,i,0,0,&TTQuantiSH::getQst,0,0);
  }
  
  if (_nb_trait > 1) {
    unsigned int c = 0;
    for(unsigned int t = 0; t < _nb_trait; ++t) {
      for(unsigned int v = t + 1, cov = (t+1)*10+(v+1) ; v < _nb_trait; ++v) {
        t1 = tstring::int2str(cov);
        add("", name + t1 +".cov",AGE,c++,0,0,&TTQuantiSH::getCovar,0,0);
      }
    }
  }
}
// ----------------------------------------------------------------------------------------
// addEigen
// ----------------------------------------------------------------------------------------
void TTQuantiSH::addEigen (age_t AGE)
{
  if(_nb_trait < 2) {
    warning("not recording G-matrix eigen stats with only one \"quanti\" trait!\n");
    return;
  }
  
  if (AGE == ALL) {
    addEigen(ADULTS);
    addEigen(OFFSPRG);
    return;
  }
  
  string suffix = (AGE == ADULTS ? "adlt.":"off.");
  
  void (TTQuantiSH::* setter) (void) = (AGE == ADULTS ?
                                        &TTQuantiSH::setAdultStats : &TTQuantiSH::setOffsprgStats);
  
  add("", suffix + "q.eval1",AGE,0,0,0,&TTQuantiSH::getEigenValue, 0, setter); //this one calls the setter
  
  for(unsigned int t = 1; t < _nb_trait; ++t)
    add("", suffix + "q.eval" + tstring::int2str(t+1),AGE,t,0,0,&TTQuantiSH::getEigenValue, 0, 0);
  
  for(unsigned int t = 0; t< _nb_trait; ++t)
    for(unsigned int v = 0; v < _nb_trait; ++v)
      add("", suffix + "q.evect" + tstring::int2str((t+1)*10+(v+1)),AGE,t,v,0,0,&TTQuantiSH::getEigenVectorElt,0);
  
}
// ----------------------------------------------------------------------------------------
// addEigenValues
// ----------------------------------------------------------------------------------------
void TTQuantiSH::addEigenValues (age_t AGE)
{
  if(_nb_trait < 2) {
    warning("not recording G-matrix eigen stats with only one \"quanti\" trait!\n");
    return;
  }
  
  if (AGE == ALL) {
    addEigen(ADULTS);
    addEigen(OFFSPRG);
    return;
  }
  
  string suffix = (AGE == ADULTS ? "adlt.":"off.");
  
  void (TTQuantiSH::* setter) (void) = (AGE == ADULTS ?
                                        &TTQuantiSH::setAdultStats : &TTQuantiSH::setOffsprgStats);
  
  add("", suffix + "q.eval1",AGE,0,0,0,&TTQuantiSH::getEigenValue, 0, setter);
  
  for(unsigned int t = 1; t < _nb_trait; ++t)
    add("", suffix + "q.eval" + tstring::int2str(t+1),AGE,t,0,0,&TTQuantiSH::getEigenValue, 0, 0);
  
  
}
// ----------------------------------------------------------------------------------------
// addEigenVect1 : save only the first eigenvector
// ----------------------------------------------------------------------------------------
void TTQuantiSH::addEigenVect1 (age_t AGE)
{
  if(_nb_trait < 2) {
    warning("not recording G-matrix eigen stats with only one \"quanti\" trait!\n");
    return;
  }
  
  if (AGE == ALL) {
    addEigen(ADULTS);
    addEigen(OFFSPRG);
    return;
  }
  
  string suffix = (AGE == ADULTS ? "adlt.":"off.");
  
  void (TTQuantiSH::* setter) (void) = (AGE == ADULTS ?
                                        &TTQuantiSH::setAdultStats : &TTQuantiSH::setOffsprgStats);
  
  add("", suffix + "q.evect11",AGE,0,0,0,0,&TTQuantiSH::getEigenVectorElt,setter);
  
  for(unsigned int v = 1; v < _nb_trait; ++v)
    add("", suffix + "q.evect1" + tstring::int2str(v+1),AGE,0,v,0,0,&TTQuantiSH::getEigenVectorElt,0);
  
}
// ----------------------------------------------------------------------------------------
// addQuantiPerPatch
// ----------------------------------------------------------------------------------------
void TTQuantiSH::addQuantiPerPatch (age_t AGE)
{
  
  if (AGE == ALL) {
    addQuantiPerPatch(ADULTS);
    addQuantiPerPatch(OFFSPRG);
    return;
  }
  
  addAvgPerPatch(AGE);
  addVarPerPatch(AGE);
  addCovarPerPatch(AGE);
  addEigenPerPatch(AGE);
  
}
// ----------------------------------------------------------------------------------------
// addAvgPerPatch
// ----------------------------------------------------------------------------------------
void TTQuantiSH::addAvgPerPatch (age_t AGE)
{
  if (AGE == ALL) {
    addAvgPerPatch(ADULTS);
    addAvgPerPatch(OFFSPRG);
    return;
  }
  
  unsigned int patchNbr = _pop->getPatchNbr();
  
  string suffix = (AGE == ADULTS ? "adlt.":"off.");
  string name;
  string patch;
  string t1;
  
  void (TTQuantiSH::* setter) (void) = (AGE == ADULTS ?
                                        &TTQuantiSH::setAdultStats : &TTQuantiSH::setOffsprgStats);
  
  add("Mean phenotype of trait 1 in patch 1", suffix + "q1.p1",  AGE, 0, 0, 0, 0,
      &TTQuantiSH::getMeanPhenotPerPatch, setter);
  
  for(unsigned int p = 0; p < patchNbr; p++) {
    for(unsigned int i = 0; i < _nb_trait; i++) {
      if(p == 0 && i == 0) continue;
      name = "Mean phenotype of trait " + tstring::int2str(i+1) + " in patch " + tstring::int2str(p+1);
      t1 = "q" + tstring::int2str(i+1);
      patch = ".p" + tstring::int2str(p+1);
      add(name, suffix + t1 + patch, AGE, i, p, 0, 0, &TTQuantiSH::getMeanPhenotPerPatch, 0);
    } 
  }
  
}
// ----------------------------------------------------------------------------------------
// addVarPerPatch
// ----------------------------------------------------------------------------------------
void TTQuantiSH::addVarPerPatch (age_t AGE)
{
  if (AGE == ALL) {
    addVarPerPatch(ADULTS);
    addVarPerPatch(OFFSPRG);
    return;
  }
  
  unsigned int patchNbr = _pop->getPatchNbr();
  
  string suffix = (AGE == ADULTS ? "adlt.":"off.");
  string name;
  string patch;
  string t1;
  
  void (TTQuantiSH::* setter) (void) = (AGE == ADULTS ?
                                        &TTQuantiSH::setAdultStats :
                                        &TTQuantiSH::setOffsprgStats);
  
  add("Genetic variance of trait 1 in patch 1", suffix + "Va.q1.p1",  AGE, 0, 0, 0, 0,
      &TTQuantiSH::getVaPerPatch, setter);
  
  for(unsigned int p = 0; p < patchNbr; p++) {
    for(unsigned int i = 0; i < _nb_trait; i++) {
      if(p == 0 && i == 0) continue;
      name = "Genetic variance of trait " + tstring::int2str(i+1) + " in patch " + tstring::int2str(p+1);
      t1 = "q" + tstring::int2str(i+1);
      patch = ".p" + tstring::int2str(p+1);
      add(name, suffix + "Va." + t1 + patch, AGE, i, p, 0, 0, &TTQuantiSH::getVaPerPatch, 0);
    }
  }
  
  if(_eVar) {
    for(unsigned int p = 0; p < patchNbr; p++) {
      for(unsigned int i = 0; i < _nb_trait; i++) {
        name = "Phenotypic variance of trait " + tstring::int2str(i+1) + " in patch " + tstring::int2str(p+1);
        t1 = "q" + tstring::int2str(i+1);
        patch = ".p" + tstring::int2str(p+1);
        add(name, suffix + "Vp." + t1 + patch, AGE, i, p, 0, 0, &TTQuantiSH::getVpPerPatch, 0);
      }
    }
  }
}
// ----------------------------------------------------------------------------------------
// addCovarPerPatch
// ----------------------------------------------------------------------------------------
void TTQuantiSH::addCovarPerPatch (age_t AGE)
{
  if(_nb_trait < 2) {
    warning("not recording traits covariance with only one \"quanti\" trait!\n");
    return;
  }
  
  if (AGE == ALL) {
    addCovarPerPatch(ADULTS);
    addCovarPerPatch(OFFSPRG);
    return;
  }
  
  string suffix = (AGE == ADULTS ? "adlt.":"off.");
  string patch;
  string cov;
  unsigned int patchNbr = _pop->getPatchNbr();
  
  void (TTQuantiSH::* setter) (void) = (AGE == ADULTS ?
                                        &TTQuantiSH::setAdultStats : &TTQuantiSH::setOffsprgStats);
  
  add("Genetic covariance of trait 1 and trait 2 in patch 1", suffix + "cov.q12.p1",  AGE, 0, 0, 0, 0,
      &TTQuantiSH::getCovarPerPatch, setter);
  
  unsigned int c;
  for(unsigned int p = 0; p < patchNbr; p++) {
    patch = ".p" + tstring::int2str(p+1);
    c = 0;
    for(unsigned int t = 0; t < _nb_trait; ++t) {
      for(unsigned int v = t + 1; v < _nb_trait; ++v){
        if(p==0 && t==0 && v==1) {c++; continue;}
        cov = tstring::int2str((t+1)*10+v+1);
        add("", suffix + "cov.q" + cov + patch,  AGE, p, c++, 0, 0, &TTQuantiSH::getCovarPerPatch, 0);
      }
    }
  }
  
}
// ----------------------------------------------------------------------------------------
// addEigenPerPatch
// ----------------------------------------------------------------------------------------
void TTQuantiSH::addEigenPerPatch (age_t AGE)
{
  if(_nb_trait < 2) {
    warning("not recording G-matrix eigen stats with only one \"quanti\" trait!\n");
    return;
  }
  
  if (AGE == ALL) {
    addEigenPerPatch(ADULTS);
    addEigenPerPatch(OFFSPRG);
    return;
  }
  
  unsigned int patchNbr = _pop->getPatchNbr();
  string suffix = (AGE == ADULTS ? "adlt.":"off.");
  string patch;
  unsigned int pv =0;
  
  void (TTQuantiSH::* setter) (void) = (AGE == ADULTS ?
                                        &TTQuantiSH::setAdultStats : &TTQuantiSH::setOffsprgStats);
  
  
  add("First G-matrix eigenvalue in patch 1", suffix + "qeval1.p1",  AGE, 0, 0, 0, 0,
      &TTQuantiSH::getEigenValuePerPatch, setter);
  
  for(unsigned int p = 0; p < patchNbr; ++p) {
    patch = ".p" + tstring::int2str(p+1);
    for(unsigned int t = 0; t < _nb_trait; ++t) {
      if(p==0 && t==0) continue;
      add("", suffix + "qeval" + tstring::int2str(t+1) + patch,  AGE, t, p, 0, 0, &TTQuantiSH::getEigenValuePerPatch,0);
    }
  }
  for(unsigned int p = 0; p < patchNbr; ++p) {
    patch = ".p" + tstring::int2str(p+1);
    pv = 0;
    for(unsigned int t = 0; t < _nb_trait; ++t)
      for(unsigned int v = 0; v < _nb_trait; ++v)
        add("", suffix + "qevect" + tstring::int2str((t+1)*10+v+1) + patch,  AGE, p, pv++, 0, 0, &TTQuantiSH::getEigenVectorEltPerPatch,0);
  }
  
}
// ----------------------------------------------------------------------------------------
// addEigenValuesPerPatch
// ----------------------------------------------------------------------------------------
void TTQuantiSH::addEigenValuesPerPatch (age_t AGE)
{
  if(_nb_trait < 2) {
    warning("not recording G-matrix eigen stats with only one \"quanti\" trait!\n");
    return;
  }
  
  if (AGE == ALL) {
    addEigenPerPatch(ADULTS);
    addEigenPerPatch(OFFSPRG);
    return;
  }
  
  unsigned int patchNbr = _pop->getPatchNbr();
  string suffix = (AGE == ADULTS ? "adlt.":"off.");
  string patch;
  
  void (TTQuantiSH::* setter) (void) = (AGE == ADULTS ?
                                        &TTQuantiSH::setAdultStats : &TTQuantiSH::setOffsprgStats);
  
  add("First G-matrix eigenvalue in patch 1", suffix + "qeval1.p1",  AGE, 0, 0, 0, 0,
      &TTQuantiSH::getEigenValuePerPatch, setter);
  
  for(unsigned int p = 0; p < patchNbr; ++p) {
    patch = ".p" + tstring::int2str(p+1);
    for(unsigned int t = 0; t < _nb_trait; ++t) {
      if(p==0 && t==0) continue;
      add("", suffix + "qeval" + tstring::int2str(t+1) + patch,  AGE, t, p, 0, 0, &TTQuantiSH::getEigenValuePerPatch,0);
    }
  }  
}
// ----------------------------------------------------------------------------------------
// addEigenVect1PerPatch
// ----------------------------------------------------------------------------------------
void TTQuantiSH::addEigenVect1PerPatch (age_t AGE)
{
  if(_nb_trait < 2) {
    warning("not recording G-matrix eigen stats with only one \"quanti\" trait!\n");
    return;
  }
  
  if (AGE == ALL) {
    addEigenPerPatch(ADULTS);
    addEigenPerPatch(OFFSPRG);
    return;
  }
  
  unsigned int patchNbr = _pop->getPatchNbr();
  string suffix = (AGE == ADULTS ? "adlt.":"off.");
  string patch;
  unsigned int pv =0;
  
  void (TTQuantiSH::* setter) (void) = (AGE == ADULTS ?
                                        &TTQuantiSH::setAdultStats : 
                                        &TTQuantiSH::setOffsprgStats);
  
  
  add("First G-matrix eigenvector in patch 1", suffix + "qevect11.p1",  AGE, 0, 0, 0, 0,
      &TTQuantiSH::getEigenVectorEltPerPatch, setter);
  
  for(unsigned int p = 0; p < patchNbr; ++p) {
    patch = ".p" + tstring::int2str(p+1);
    pv = 0;
    //    for(unsigned int t = 0; t < _nb_trait; ++t)
    for(unsigned int v = 0; v < _nb_trait; ++v){
      if(p==0 && v==0) {pv++; continue;}
      add("", suffix + "qevect1" + tstring::int2str(v+1) + patch,  AGE, p, pv++, 0, 0, &TTQuantiSH::getEigenVectorEltPerPatch,0);
    }
  }
}
// ----------------------------------------------------------------------------------------
// addSkewPerPatch
// ----------------------------------------------------------------------------------------
void TTQuantiSH::addSkewPerPatch(age_t AGE)
{
  if (AGE == ALL) {
    addSkewPerPatch(ADULTS);
    addSkewPerPatch(OFFSPRG);
    return;
  }
  
  unsigned int patchNbr = _pop->getPatchNbr();
  
  string suffix = (AGE == ADULTS ? "adlt.":"off.");
  string name;
  string patch;
  string t1;
  
  void (TTQuantiSH::* setter) (void) = (AGE == ADULTS ?
                                        &TTQuantiSH::setAdultStats :
                                        &TTQuantiSH::setOffsprgStats);
  
  add("Genetic skew of trait 1 in patch 1", suffix + "Sk.q1.p1",  AGE, 0, 0, 0, 0,
      &TTQuantiSH::getSkewPerPatch, setter);
  
  for(unsigned int p = 0; p < patchNbr; p++) {
    for(unsigned int i = 0; i < _nb_trait; i++) {
      if(p == 0 && i == 0) continue;
      name = "Genetic skew of trait " + tstring::int2str(i+1) + " in patch " + tstring::int2str(p+1);
      t1 = "q" + tstring::int2str(i+1);
      patch = ".p" + tstring::int2str(p+1);
      add(name, suffix + "Sk." + t1 + patch, AGE, i, p, 0, 0, &TTQuantiSH::getSkewPerPatch, 0);
    }
  }
  
}
// ----------------------------------------------------------------------------------------
// getSkewPerPatch
// ----------------------------------------------------------------------------------------
double TTQuantiSH::getSkewPerPatch (unsigned int i, unsigned int p)
{
  double skew = 0;
  
  double *phenot = _phenoTable.getClassWithinGroup(i, p);
  unsigned int patch_size = _phenoTable.size(i, p);  
  
  for(unsigned int k = 0; k < patch_size; ++k)
    skew += pow( phenot[k] - _pmeanP[i][p], 3);  //the mean has been set by setStats()
  
  return skew / patch_size;
}
// ----------------------------------------------------------------------------------------
// setDataTables
// ----------------------------------------------------------------------------------------
void TTQuantiSH::setDataTables(age_t AGE) 
{
  unsigned int **sizes;
  unsigned int nb_patch = _pop->getPatchNbr();
  
  sizes = new unsigned int * [_nb_trait];
  
  for(unsigned int i = 0; i < _nb_trait; ++i) {
    sizes[i] = new unsigned int [nb_patch];
    for(unsigned int j = 0; j < nb_patch; ++j)
      sizes[i][j] = _pop->size(AGE, j);
  }
  
  _phenoTable.update(_nb_trait, nb_patch, sizes);
  _genoTable.update(_nb_trait, nb_patch, sizes);
  
  for(unsigned int i = 0; i < _nb_trait; ++i)
    delete [] sizes[i];
  delete [] sizes;
  
  Patch* patch;
  //age_idx age = (AGE == ADULTS ? ADLTx : OFFSx);
  
  for(unsigned int i = 0, n; i < nb_patch; i++) {
    
    patch = _pop->getPatch(i);
    
    n=0;
    
    if ((patch->size(MAL, AGE)+patch->size(FEM, AGE)) != _phenoTable.size(0,i)) {
      fatal("problem while recording quanti trait values; table size doesn't match patch size.\n");
    }
    store_quanti_trait_values(patch, i, patch->size(MAL, AGE), &n, MAL, AGE, &_phenoTable, &_genoTable,
                       _nb_trait, _SHLinkedTraitIndex);
    
    store_quanti_trait_values(patch, i, patch->size(FEM, AGE), &n, FEM, AGE, &_phenoTable, &_genoTable,
                       _nb_trait, _SHLinkedTraitIndex);
    
    if (n != _phenoTable.size(0,i) || n != _genoTable.size(0,i)) {
      fatal("problem while recording quanti trait values; size counter doesn't match table size.\n");
    }
  }
}
// ----------------------------------------------------------------------------------------
// store_trait_values
// ----------------------------------------------------------------------------------------
void store_quanti_trait_values (Patch* patch, unsigned int patchID, unsigned int size, unsigned int *cntr,
                         sex_t SEX, age_t AGE, DataTable<double> *ptable, DataTable<double> *gtable,
                         unsigned int nTrait, unsigned int TraitIndex)
{
  double *phe;
  TTQuanti *trait;
  
  for(unsigned int j = 0; j < size; ++j) {
    
    trait = dynamic_cast<TTQuanti*> (patch->get(SEX, AGE, j)->getTrait( TraitIndex ));
    
    phe = (double*)trait->getValue();
    
    for(unsigned int k = 0; k < nTrait; k++){
      ptable->set( k, patchID, (*cntr), phe[k]);
      gtable->set( k, patchID, (*cntr), trait->get_genotype(k));
    }
    (*cntr)++;
  }
  
}
// ----------------------------------------------------------------------------------------
// setStats
// ----------------------------------------------------------------------------------------
void TTQuantiSH::setStats (age_t AGE)
{  
  if(_table_set_age == AGE 
     && _table_set_gen == _pop->getCurrentGeneration()
     && _table_set_repl == _pop->getCurrentReplicate())
    return;
  
  unsigned int pop_size = _pop->size(AGE);
  unsigned int nb_patch = _pop->getPatchNbr();
  unsigned int patch_size;
  double *phenot1, *genot1, *genot2;
  
  setDataTables(AGE);
  
#ifdef HAS_GSL
  unsigned int c = 0; //covariance position counter
  unsigned int pv = 0; //eigenvector position counter
  //within deme stats:
  for(unsigned int j = 0; j < nb_patch; j++) { 
    
    patch_size = _pop->size(AGE, j);
    
    for(unsigned int t=0; t < _nb_trait; t++) {
      
      phenot1 = _phenoTable.getClassWithinGroup(t,j);
      genot1  = _genoTable.getClassWithinGroup(t,j);
      
      _pmeanP[t][j] = my_mean (phenot1, patch_size );
      _pmeanG[t][j] = my_mean (genot1, patch_size );
      _pVp[t][j] = my_variance_with_fixed_mean (phenot1, patch_size, _pmeanP[t][j]);
      _pVa[t][j] = my_variance_with_fixed_mean (genot1,  patch_size, _pmeanG[t][j]);
           
    }
    
    if(_nb_trait > 1) { 
      c = 0;
      //    calculate the covariances and G, need to adjust dimensions in class declaration
      for(unsigned int t1 = 0; t1 < _nb_trait; t1++) {
        //      set the diagonal elements of G here
        
        gsl_matrix_set(_G, t1, t1, _pVa[t1][j]);
        
        for(unsigned int t2 = t1 + 1; t2 < _nb_trait; t2++) {
          
          genot1  = _genoTable.getClassWithinGroup(t1,j);
          genot2  = _genoTable.getClassWithinGroup(t2,j);
          
          _pcovar[j][c] = gsl_stats_covariance_m (genot1, 1, genot2, 1, patch_size,
                                                  _pmeanG[t1][j], _pmeanG[t2][j]);
          
          gsl_matrix_set(_G, t1, t2, _pcovar[j][c]);
          gsl_matrix_set(_G, t2, t1, _pcovar[j][c++]);
        }
      }
      
      gsl_eigen_symmv (_G, _eval, _evec, _ws);
      gsl_eigen_symmv_sort (_eval, _evec, GSL_EIGEN_SORT_VAL_DESC);
      
      pv = 0;
      
      for(unsigned int t = 0; t < _nb_trait; t++) {
        _peigval[t][j] = gsl_vector_get (_eval, t);
        for(unsigned int v = 0; v < _nb_trait; v++)
          _peigvect[j][pv++] = gsl_matrix_get (_evec, v, t); //read eigenvectors column-wise
      }
    }
  }
  
  double meanGamong1, meanGamong2;
  c = 0; //reset covariance positioner
  //among demes stats:
  for(unsigned int t1 = 0; t1 < _nb_trait; t1++) {
    
    phenot1 = _phenoTable.getGroup(t1);
    genot1  = _genoTable.getGroup(t1);
    
    _meanP[t1] = my_mean (phenot1, pop_size ); //grand mean, pooling all individuals
    _meanG[t1] = my_mean (genot1, pop_size ); //grand mean, pooling all individuals
    
    _Vp[t1] = my_mean_no_nan (_pVp[t1], nb_patch); //mean within patch variance
    _Va[t1] = my_mean_no_nan (_pVa[t1], nb_patch); //mean within patch variance
    
    meanGamong1 = my_mean_no_nan (_pmeanG[t1], nb_patch); //mean of within patch mean genotypic values
    
    _Vb[t1] = my_variance_with_fixed_mean_no_nan (_pmeanG[t1], nb_patch,  meanGamong1); //variance of patch means
    
    gsl_matrix_set(_G, t1, t1, _Vb[t1]); //_G here becomes the D-matrix, the among-deme (Difference) covariance matrix 
    
    for(unsigned int t2 = t1 + 1; t2 < _nb_trait; t2++) {
      
      meanGamong2 = my_mean (_pmeanG[t2], nb_patch);
      
      _covar[c] = gsl_stats_covariance_m (_pmeanG[t1], 1, _pmeanG[t2], 1, nb_patch, meanGamong1, meanGamong2); //covariance of patch means
      
      gsl_matrix_set(_G, t1, t2, _covar[c]);
      gsl_matrix_set(_G, t2, t1, _covar[c++]);
      
    }
  }
  
  if(_nb_trait > 1) { 
    gsl_eigen_symmv (_G, _eval, _evec, _ws);
    gsl_eigen_symmv_sort (_eval, _evec, GSL_EIGEN_SORT_VAL_DESC);
    
    for(unsigned int t1 = 0; t1 < _nb_trait; t1++) {
      _eigval[t1] = gsl_vector_get (_eval, t1);
      for(unsigned int t2 = 0; t2 < _nb_trait; t2++) {
        _eigvect[t2][t1] = gsl_matrix_get (_evec, t2, t1);      
      }
    }
  }
  
#else
  fatal("install the GSL library to get the quanti stats!\n");
#endif
  
  _table_set_age = AGE;
  _table_set_gen = _pop->getCurrentGeneration();  
  _table_set_repl = _pop->getCurrentReplicate();
  
}
// ----------------------------------------------------------------------------------------
// FHwrite
// ----------------------------------------------------------------------------------------
void TTQuantiFH::FHwrite()
{
  
  Metapop* pop = get_pop_ptr();
  
  if (!pop->isAlive()) return;
  
  std::string filename = get_filename();
  
  std::ofstream FILE (filename.c_str(), ios::out);
  
  if(!FILE) fatal("could not open \"%s\" output file!!\n",filename.c_str());
  
  bool print_gene = (_output_option == "genotypes" || _output_option == "genotype");
  bool print_genotype = (!_FHLinkedTrait->get_env_var().empty());
  
  FILE<<"pop ";
  if(print_gene) {
    for(unsigned int k = 0; k < _FHLinkedTrait->get_nb_traits(); k++)
      for(unsigned int l = 0; l < _FHLinkedTrait->get_nb_locus(); l++)
        FILE<<"t"<<k+1<<"l"<<l+1<<"1 "<<"t"<<k+1<<"l"<<l+1<<"2 ";
  }
  
  for(unsigned int k = 0; k < _FHLinkedTrait->get_nb_traits(); k++) {
    FILE<<"P"<<k+1<< " "; 
    if(print_genotype) FILE<<"G"<<k+1<< " ";
  }
  
  FILE<<"age sex home ped isMigrant father mother ID\n";
  
  age_t pop_age = pop->getCurrentAge(); //flag telling which age class should contain individuals
  
  //we print anything that is present in the pop:
  if( (pop_age & OFFSPRG) != 0) print(FILE, OFFSx, print_gene, print_genotype);
  
  unsigned int nb_class = _pop->getLeslieMatrix()->getNbCols();

  if( (pop_age & ADULTS) != 0) {
  
  for(unsigned int i = 0; i < nb_class; i++ ){

		  print(FILE, static_cast<age_idx> (i), print_gene, print_genotype);
  }
  
  }
  FILE.close();  
}
// ----------------------------------------------------------------------------------------
// FHwrite
// ----------------------------------------------------------------------------------------
void TTQuantiFH::print(ofstream& FH, age_idx Ax, bool print_gene, bool print_genotype)
{
  Metapop* pop = get_pop_ptr();
  int patchNbr = pop->getPatchNbr();
  Patch* current_patch;
  Individual* ind;
  TTQuanti* trait;
  double* Tval;
  double **genes;
  unsigned int loc, nb_trait=_FHLinkedTrait->get_nb_traits(), nb_locus = _FHLinkedTrait->get_nb_locus();
  
  for(int i = 0; i < patchNbr; i++) {
    
    current_patch = pop->getPatch(i);
    
    for(unsigned int j = 0, size = current_patch->size(FEM, Ax); j < size; j++) {
      
      ind = current_patch->get(FEM, Ax, j);
      trait = dynamic_cast<TTQuanti*> (ind->getTrait(_FHLinkedTraitIndex));
      
      FH<<i+1<<" ";
      Tval = (double*)trait->getValue();
      
      if(print_gene){
        genes = (double**)trait->get_sequence();
        
        FH.precision(6);
        
        for(unsigned int k = 0; k < nb_trait; k++) {
          for(unsigned int l = 0; l < nb_locus; l++) {
            loc = l * nb_trait + k;
            
            FH<<genes[0][loc]<<" "<<genes[1][loc]<<" ";
          }
        }
      }
      
      FH.precision(4);
      for(unsigned int k = 0; k < nb_trait; k++) {
        FH<<Tval[k]<<" ";
        if(print_genotype) FH << trait->get_genotype(k) << " ";
      }
      
      FH<<Ax<<" "<<ind->getSex()<<" "<<ind->getHome()+1<<" "<<ind->getPedigreeClass(ind->getMother(), ind->getFather())<<" "
      << (ind->getFather() && ind->getMother() ?
          (ind->getFather()->getHome()!=i) + (ind->getMother()->getHome()!=i) : 0)
//      <<" "<<ind->getFatherID()<<" "<<ind->getMotherID()<<" "<<ind->getID()<<std::endl;//MOD SLIM NEMO 7.11
        <<" NA NA "<<ind->getID()<<std::endl;
    }
    
    for(unsigned int j = 0, size = current_patch->size(MAL, Ax); j < size; j++) {
      
      ind = current_patch->get(MAL, Ax, j);
      trait = dynamic_cast<TTQuanti*> (ind->getTrait(_FHLinkedTraitIndex));
      
      FH<<i+1<<" ";
      
      Tval = (double*)trait->getValue();
      
      if(print_gene){
        genes = (double**)trait->get_sequence();
        
        FH.precision(6);
        
        for(unsigned int k = 0; k < nb_trait; k++) {
          for(unsigned int l = 0; l < nb_locus; l++) {
            loc = l * nb_trait + k;
            
            FH<<genes[0][loc]<<" "<<genes[1][loc]<<" ";
          }
        }
      }
      
      FH.precision(4);
      for(unsigned int k = 0; k < nb_trait; k++) {
        FH<<Tval[k]<<" ";
        if(print_genotype) FH << trait->get_genotype(k) << " ";
      }
      
      FH<<Ax<<" "<<ind->getSex()<<" "<<ind->getHome()+1<<" "<<ind->getPedigreeClass(ind->getMother(), ind->getFather())<<" "
      << (ind->getFather() && ind->getMother() ?
          (ind->getFather()->getHome()!=i) + (ind->getMother()->getHome()!=i) : 0)
//	  <<" "<<ind->getFatherID()<<" "<<ind->getMotherID()<<" "<<ind->getID()<<std::endl;//MOD SLIM NEMO 7.11
            <<" NA NA "<<ind->getID()<<std::endl;
    }
  }
}
// ----------------------------------------------------------------------------------------
// FHwrite
// ----------------------------------------------------------------------------------------
void TTQFreqExtractor::FHwrite()
{
  Metapop* pop = get_pop_ptr();
  int patchNbr = pop->getPatchNbr();
  Patch* current_patch;
  Individual* ind;
  double **seq;
  double **val_t1, **val_t2, *dist_t1[2];//, *dist_t2[2];
  unsigned int t1=0, t2=0, bloc, genome_size=_FHLinkedTrait->get_nb_locus(),
  _nb_trait=_FHLinkedTrait->get_nb_traits(), nb_locus = _FHLinkedTrait->get_nb_locus();
  
  std::string filename = get_filename();
  
  std::ofstream FILE (filename.c_str(), ios::out);
  
  if(!FILE) fatal("could not open \"%s\" output file!!\n",filename.c_str());
  
  FILE.close();
  
  age_t pop_age = pop->getCurrentAge();
  
  if( (pop_age & ADULTS) != 0) {
    
    unsigned int val_size = pop->size(ADULTS) * _nb_trait * 2;
    //    cout<<"allocating the val arrays (size "<<val_size<<")"<<endl;
    val_t1 = new double* [ nb_locus ];
    val_t2 = new double* [ nb_locus ];
    
    for(unsigned int i = 0; i < nb_locus; i++) {
      val_t1[i] = new double [ val_size ];
      val_t2[i] = new double [ val_size ];
    }
    //    cout<<"getting the values"<<endl;
    for(int i = 0; i < patchNbr; i++) {
      
      current_patch = pop->getPatch(i);
      
      for(unsigned int j = 0, size = current_patch->size(FEM, ADLTx); j < size; j++) {
        
        ind = current_patch->get(FEM, ADLTx, j);
        
        seq = (double**)ind->getTrait(_FHLinkedTraitIndex)->get_sequence();
        
        for(unsigned int k = 0; k < nb_locus; ++k) {
          bloc = k * _nb_trait;
          t2 = t1 = bloc * 2;
          val_t1[k][ t1++ ] = seq[0][bloc];
          val_t1[k][ t1++ ] = seq[1][bloc];
          val_t1[k][ t1++ ] = seq[0][bloc + genome_size];
          val_t1[k][ t1 ] = seq[1][bloc + genome_size];
          bloc++;
          val_t2[k][ t2++ ] = seq[0][bloc];
          val_t2[k][ t2++ ] = seq[1][bloc];
          val_t2[k][ t2++ ] = seq[0][bloc + genome_size];
          val_t2[k][ t2 ] = seq[1][bloc + genome_size];
        }
        
      }
      
      for(unsigned int j = 0, size = current_patch->size(MAL, ADLTx); j < size; j++) {
        ind = current_patch->get(MAL, ADLTx, j);
        seq = (double**)ind->getTrait(_FHLinkedTraitIndex)->get_sequence();
        
        for(unsigned int k = 0; k < nb_locus; ++k) {
          bloc = k * _nb_trait;
          t2 = t1 = bloc * 2;
          val_t1[k][ t1++ ] = seq[0][bloc];
          val_t1[k][ t1++ ] = seq[1][bloc];
          val_t1[k][ t1++ ] = seq[0][bloc + genome_size];
          val_t1[k][ t1 ] = seq[1][bloc + genome_size];
          bloc++;
          val_t2[k][ t2++ ] = seq[0][bloc];
          val_t2[k][ t2++ ] = seq[1][bloc];
          val_t2[k][ t2++ ] = seq[0][bloc + genome_size];
          val_t2[k][ t2 ] = seq[1][bloc + genome_size];
        }
      }
    }//end for patchNbr
    //dump the values to a text file
    filename += "tot";
    
    FILE.open(filename.c_str(), ios::out);
    
    if(!FILE) fatal("could not open \"%s\" output file!!\n",filename.c_str());
    
    for(unsigned int j = 0; j < nb_locus; j++)
      FILE<<"t1.l"<<j+1<<" ";
    for(unsigned int j = 0; j < nb_locus; j++)
      FILE<<"t2.l"<<j+1<<" ";
    
    FILE<<endl;
    
    for(unsigned int i = 0; i < val_size; i++) {
      
      for(unsigned int j = 0; j < nb_locus; j++)
        FILE<<val_t1[j][i]<<" ";
      
      for(unsigned int j = 0; j < nb_locus; j++)
        FILE<<val_t2[j][i]<<" ";
      
      FILE<<endl;
    }
    FILE.close();
    
    //find min and max:
    double min1=0, max1=0, range1, min2=0, max2=0, range2;
    //    cout<<"find the bounds"<<endl;
    for(unsigned int i = 0; i < nb_locus; i++) 
      for(unsigned int j = 0; j < val_size; i++) {
        min1 = (min1 < val_t1[i][j] ? min1 : val_t1[i][j]);
        max1 = (max1 > val_t1[i][j] ? max1 : val_t1[i][j]);
        min2 = (min2 < val_t2[i][j] ? min2 : val_t2[i][j]);
        max2 = (max2 > val_t2[i][j] ? max2 : val_t2[i][j]);
      }
    
    range1 = max1 - min1;
    range2 = max2 - min2;
    
    //    cout<<"trait1: max "<<max1<<" min "<<min1<<" range "<<range1<<endl;
    //    cout<<"trait2: max "<<max2<<" min "<<min2<<" range "<<range2<<endl;
    max1 = (max1 > max2 ? max1 : max2);
    min1 = (min1 < min2 ? min1 : min2);
    range1 = max1 - min1;
    //    cout<<"total: max "<<max1<<" min "<<min1<<" range "<<range1<<endl;
    
    unsigned int dist_size1 = (unsigned int)ceil(range1 / _granularity);
    //    cout<<"allocating the dist array size "<<dist_size1<<endl;
    dist_t1[0] = new double [dist_size1];
    dist_t1[1] = new double [dist_size1];
    
    for(unsigned int i = 0; i < dist_size1; ++i) 
      dist_t1[0][i] = dist_t1[1][i] = 0;
    
    //    unsigned int dist_size2 = (unsigned int)ceil(range1 / _granularity);
    //    cout<<"allocating the dist array size "<<dist_size<<endl;
    //    dist_t2[0] = new double [dist_size2];
    //    dist_t2[1] = new double [dist_size2];    
    //    
    //    for(unsigned int i = 0; i < dist_size2; ++i) 
    //      dist_t2[0] = dist_t2[1] = 0;
    
    for(unsigned int i = 0; i < nb_locus; i++)
      for(unsigned int j = 0; j < val_size; j++) {
        val_t1[i][j] -= min1;
        dist_t1[1][ (unsigned int)floor(val_t1[i][j] / _granularity) ]++;
        
        val_t2[i][j] -= min1;
        dist_t1[1][ (unsigned int)floor(val_t2[i][j] / _granularity) ]++;
      }
    
    min1 += _granularity/2.0;
    //    min2 += _granularity/2.0;
    filename = get_filename();
    
    FILE.open(filename.c_str(), ios::out);
    
    if(!FILE) fatal("could not open \"%s\" output file!!\n",filename.c_str());
    
    FILE<<"val freq"<<endl;
    
    val_size *= 2 * nb_locus;
    for(unsigned int i = 0; i < dist_size1; ++i) {
      dist_t1[1][i] /= val_size;
      dist_t1[0][i] = min1 + i * _granularity;
      FILE<<dist_t1[0][i]<<" "<<dist_t1[1][i]<<endl;
    }
    
    FILE.close();
    
    //    cout<<"trait2 distribution:"<<endl;
    //    for(unsigned int i = 0; i < dist_size2; ++i) {
    //      dist_t2[1][i] /= val_size;
    //      dist_t2[0][i] = min2 + i * _granularity;
    //      cout<<dist_t2[0][i]<<"\t"<<dist_t2[1][i]<<endl;
    //    }
    
    for(unsigned int i = 0; i < nb_locus; i++){
      delete [] val_t1[i];
      delete [] val_t2[i];
    }
    
    delete [] val_t1;
    delete [] val_t2;
    
  } else
    warning("TTQuantiFH::FHwrite:metapop empty or age flag not set, not writing.\n");
}
