/*
 * ttquantiPhenotypic.cc
 *
 *  Created on: 21.04.2014
 *      Author: cotto_o
 */


#include <sstream>
#include <fstream>
#include "ttquantiphenotypic.h"
#include <string.h>
#include <cmath>
#include <algorithm>
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


//local functions needed for the stats:
void store_pheno_trait_values (Patch* patch, unsigned int patchID, unsigned int size, unsigned int *cntr,
                         sex_t SEX, age_t AGE, DataTable<double> *ptable, DataTable<double> *gtable,
                         unsigned int nTrait, unsigned int TraitIndex);


//CLASS//TProtoQuantiPheno

//Constructor

TProtoQuantiPheno::TProtoQuantiPheno()://(toujours initialiser les pointeurs!)
_nb_traits(0),
_genomic_mutation_rate(0),
_init_value(0),
_init_model(0),
_eVariance(0),
_mutation_variance(0),
_mutation_matrix(0),
_gsl_mutation_matrix(0),
_evect(0), 
_eval(0),
_effects_multivar(0),
_ws(0),
_mutation_sigma(0),
_mutation_correlation(0),
_mutation_func_ptr(0),
_stats(0),
_writer(0)


{

  set_paramset("phenotypic_trait", false, this);
  
  add_parameter("pheno_traits",INT,true,false,0,0);
  add_parameter("pheno_init_value",MAT,false,false,0,0);
  add_parameter("pheno_init_model",INT,false,true,0,1);
  add_parameter("pheno_environmental_variance",DBL,false,false,0,0);
  add_parameter("pheno_mutation_variance",DBL,false,false,0,0,0); // MOD F.G. 28.04.14
  add_parameter("pheno_genomic_mutation_rate",DBL,true,true,0,1,0);// MOD F.G. 02.05.14 - is mandatory
  add_parameter("pheno_mutation_correlation",DBL,false,false,0,0,0);
  add_parameter("pheno_mutation_matrix",MAT,false,false,0,0,0);
  add_parameter("pheno_output",STR,false,false,0,0);
  add_parameter("pheno_logtime",INT,false,false,0,0);
  add_parameter("pheno_dir",STR,false,false,0,0);
  
}

//Copy constructor

TProtoQuantiPheno::TProtoQuantiPheno(const TProtoQuantiPheno& T):
		//(toujours initialiser les pointeurs, meme dans le copy cstor!)

_nb_traits(T._nb_traits),
_genomic_mutation_rate(T._genomic_mutation_rate),
_mutation_variance(T._mutation_variance),
_eVariance(T._eVariance),
_init_value(0),
_init_model(T._init_model),
_mutation_matrix(0),
_gsl_mutation_matrix(0),
_evect(0), 
_eval(0),
_effects_multivar(0),
_ws(0),
_mutation_sigma(0),
_mutation_correlation(T._mutation_correlation),
_mutation_func_ptr(T._mutation_func_ptr),
_stats(0),
_writer(0)
{
  _paramSet = new ParamSet( *(T._paramSet) ) ;
}

//Destructor

TProtoQuantiPheno::~TProtoQuantiPheno ()
{// MOD F.G. 28.04.14
  reset_mutation_pointers();
  if(_stats != NULL){delete _stats; _stats = NULL;}
  if(_writer != NULL){delete _writer; _writer = NULL;}
  if(_init_value) {delete[] _init_value;_init_value = NULL;}
}

//Set Parameters - initial trait values and environmental variance
bool TProtoQuantiPheno::setParameters()
{
  
  
  _nb_traits = (unsigned int)get_parameter_value("pheno_traits");
  
  if(get_parameter("pheno_environmental_variance")->isSet())
    _eVariance = sqrt(get_parameter_value("pheno_environmental_variance"));
  else
    _eVariance = 0;
  
  TMatrix tmp_matx;
  
  if(_init_value != NULL) {
    delete [] _init_value;
    _init_value = NULL;
  }

  if(get_parameter("pheno_init_model")->isSet()) {

	  _init_model = (unsigned int)get_parameter_value("pheno_init_model");

  } else
	  _init_model = 0;

  
  if(get_parameter("pheno_init_value")->isSet()) {
    
    get_parameter("pheno_init_value")->getMatrix(&tmp_matx);
    
    if(tmp_matx.getNbRows() != 1 || tmp_matx.getNbCols() != _nb_traits) {
      fatal("\"pheno_init_value\" must be a vector of length equal to the number of traits!\n");
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
  
  _mutation_variance=get_parameter_value("pheno_mutation_variance") ; // MOD F.G. 28.04.14 (ajout 'pheno_')
  //c'est pas suffisant! les fonctions de mutation gaussiennes ont besoin de plusieurs valeurs
  //qui sont definie dans setMutationModel()
  
  if(_nb_traits == 1) {
    _mutation_func_ptr = &TProtoQuantiPheno::getMutationEffectUnivariateGaussian;
  }
  if (_nb_traits == 2) {
    _mutation_func_ptr = &TProtoQuantiPheno::getMutationEffectBivariateGaussian;
  }
  if (_nb_traits > 2) {
    _mutation_func_ptr = &TProtoQuantiPheno::getMutationEffectMultivariateGaussian;
  }

  // MOD F.G. 28.04.14
  //attention, les modele bivar & multivar on besoin de la matrice de mutation (var+covar)! 
  //il faut appeler setMutationModel() pour initialiser la matrice
  //il faut aussi definir le taux de mutation genomique:


  setMutationParameters () ;
  reset_mutation_pointers(); //doit etre imperativement appelee avant setMutationModel()...
  return setMutationModel();//remontee de setMutationParameters() pour etre plus explicite
}

bool TProtoQuantiPheno::setMutationParameters () 
//fonction n'etait jamais appelee
{
	_genomic_mutation_rate = get_parameter_value("pheno_genomic_mutation_rate");
  
  if(get_parameter("pheno_mutation_correlation")->isSet())
    _mutation_correlation = get_parameter_value("pheno_mutation_correlation");
  else
    _mutation_correlation = 0;
  return true;
}

bool TProtoQuantiPheno::setMutationModel(){  
//doit etre appele dans setParameters() pour que les fonctions de mutations fonctionnent

	unsigned int dims[2];

  //setting the mutation variance-covariance matrix
  if(get_parameter("pheno_mutation_variance")->isSet()) { // MOD F.G. 28.04.14 (quanti_ --> pheno_)
    
    if(get_parameter("pheno_mutation_matrix")->isSet()) {
      warning("both \"pheno_mutation_variance\" and \"pheno_mutation_matrix\" are set, using the matrix only!\n");
    } else {
      
      _mutation_sigma = new double [_nb_traits];
      
      double sigma = sqrt( get_parameter_value("pheno_mutation_variance") );
      
      for(unsigned int i = 0; i < _nb_traits; i++)
        _mutation_sigma[i] = sigma;
      
      if(_nb_traits > 1) {
        //setting the mutation matrix
        _mutation_matrix = new TMatrix();
        _gsl_mutation_matrix = gsl_matrix_alloc(_nb_traits, _nb_traits);
        double covar, var = get_parameter_value("pheno_mutation_variance");
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
  
  if(get_parameter("pheno_mutation_matrix")->isSet()) {
    
    _mutation_matrix = new TMatrix();
    
    get_parameter("pheno_mutation_matrix")->getMatrix(_mutation_matrix);
    
    _mutation_matrix->get_dims(&dims[0]);
    
    if( dims[0] != _nb_traits || dims[1] != _nb_traits) {
      error("\"pheno_mutation_matrix\" must be a square matrix of size = \"pheno_traits\" x \"pheno_traits\"\n");
      return false;
    }
    
    _gsl_mutation_matrix = gsl_matrix_alloc(dims[0], dims[1]);
    
    _mutation_matrix->get_gsl_matrix(_gsl_mutation_matrix);
    
    _mutation_sigma = new double [_nb_traits];
    
    for(unsigned int i = 0; i < _nb_traits; i++)
      _mutation_sigma[i] = sqrt(_mutation_matrix->get(i, i));
    
  }
  else if(!get_parameter("pheno_mutation_variance")->isSet()) {
    error("\"pheno_mutation_matrix\" or \"pheno_mutation_variance\" must be specified!\n");
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

void TProtoQuantiPheno::set_mutation_matrix_decomposition()
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

double* TProtoQuantiPheno::getMutationEffectMultivariateGaussian ()
{
  RAND::MultivariateGaussian(_eval, _evect, _ws, _effects_multivar);
  return _effects_multivar->data;
}

double* TProtoQuantiPheno::getMutationEffectBivariateGaussian ()
{
  RAND::BivariateGaussian(_mutation_sigma[0], _mutation_sigma[1], _mutation_correlation,
                          &_effects_bivar[0], &_effects_bivar[1]);
  return &_effects_bivar[0];
}

double* TProtoQuantiPheno::getMutationEffectUnivariateGaussian   ()
{
  _effects_bivar[0] = RAND::Gaussian(_mutation_sigma[0]);
  return &_effects_bivar[0];
}

void TProtoQuantiPheno::reset_mutation_pointers()
{
  
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

//Hatch

TTQuantiPheno* TProtoQuantiPheno::hatch()
{
  TTQuantiPheno* kid = new TTQuantiPheno();
  
  kid->set_proto(this);
  kid->set_nb_traits(_nb_traits);
  kid->set_init_value(_init_value);
  kid->set_mutation_fptr(_mutation_func_ptr);

  return kid;
}

void TProtoQuantiPheno::loadFileServices  (FileServices* loader)
{ 
  //writer
  if(get_parameter("pheno_output")->isSet()) {
    
    if(_writer == NULL) _writer = new TTQuantiPhenoFH(this);
    
    Param* param = get_parameter("pheno_logtime");
    
    if(param->isMatrix()) {
      
      TMatrix temp;
      param->getMatrix(&temp);
      _writer->set_multi(true, true, 1, &temp, get_parameter("pheno_dir")->getArg());
      
    } else   //  rpl_per, gen_per, rpl_occ, gen_occ, rank (0), path, self-ref      
      _writer->set(true, true, 1, (param->isSet() ? (int)param->getValue() : 0),
                   0, get_parameter("pheno_dir")->getArg(),this);
    
    loader->attach(_writer);
    
  } else if(_writer != NULL) {
    delete _writer;
    _writer = NULL;
  }
}
// ----------------------------------------------------------------------------------------
// loadStatServices
// ----------------------------------------------------------------------------------------
void TProtoQuantiPheno::loadStatServices  (StatServices* loader)
{
  //allocate the stat handler
  if(_stats == NULL)
    _stats = new TTQuantiPhenoSH(this);
  
  if(_stats != NULL) {
    loader->attach(_stats);
  }
}


//CLASS//TTQuantiPheno

//Operator =
TTQuantiPheno& TTQuantiPheno::operator= (const TTrait& T)
{
  const TTQuantiPheno& TQ = dynamic_cast<const TTQuantiPheno&>(T);
  
  if(this != &TQ) {
    _nb_traits = TQ._nb_traits;
    reset();
    init();

    // MOD F.G. 28.04.14
    //copie des valeurs genotypic des traits:
    memcpy(_sequence[0],TQ._sequence[0],_nb_traits*sizeof(double));
    memcpy(_sequence[1],TQ._sequence[1],_nb_traits*sizeof(double));
      
    set_value();
  }
  return *this;
}

//Operator ==
bool TTQuantiPheno::operator== (const TTrait& T)
{
  if(this->get_type().compare(T.get_type()) != 0) return false;
  const TTQuantiPheno& TQ = dynamic_cast<const TTQuantiPheno&>(T);
  if(this != &TQ) {
    if(_nb_traits != TQ._nb_traits) return false;
  }
  return true;
}

//Operator !=
bool TTQuantiPheno::operator!= (const TTrait& T)
{
  if(!((*this) == T) )
    return true;
  else
    return false;
}

//set initial values
void TTQuantiPheno::set_init_value(double* val)
{
	assert(_nb_traits != 0);

  if(_init_value) delete [] _init_value;

  _init_value = new double [_nb_traits];

  for(unsigned int i = 0; i < _nb_traits; ++i)
    _init_value[i] = val[i];
}

//init()
inline void TTQuantiPheno::init ()
{
  _sequence = new double*[2];
  _sequence[0] = new double [_nb_traits];
  _sequence[1] = new double [_nb_traits];
  if(!_phenotypes) _phenotypes = new double [_nb_traits];
  if(!_init_value) _init_value = new double [_nb_traits];
}

inline void TTQuantiPheno::init_sequence (){
  
//  double mut_effect = 0.01; //ne devrait pas plutot etre similaire a Vm?
//
//  mut_effect = _myProto->_mutation_variance;  //attention, peut etre different pour chaque trait
  
	for(unsigned int j = 0; j < _nb_traits; j++) {
    
    // MOD F.G. 28.04.14
    
	if (_myProto->_init_model == 1) {

	    _sequence[0][j] = _init_value[j]/2 + RAND::Gaussian(_myProto->_mutation_sigma[j]);
	    _sequence[1][j] = _init_value[j]/2 + RAND::Gaussian(_myProto->_mutation_sigma[j]);

	} else {

	    _sequence[0][j] = _init_value[j]/2;
	    _sequence[1][j] = _init_value[j]/2;

	}
  }
}

//reset
inline void TTQuantiPheno::reset ()
{
  if(_sequence) {
    if(_sequence[0]) delete [] _sequence[0];
    if(_sequence[1]) delete [] _sequence[1];
    delete [] _sequence;
    _sequence = NULL;
  }
  
  if(_phenotypes) delete [] _phenotypes;
  _phenotypes = NULL;
  
  if(_init_value) delete [] _init_value;
  _init_value = NULL;
}

//inherit
inline void TTQuantiPheno::inherit (TTrait* mother, TTrait* father)
{
  for(unsigned int i = 0; i < _nb_traits; ++i){
    _sequence[1][i] = *(double*)mother->get_allele(i, RAND::RandBool());// MOD F.G. 28.04.14
    _sequence[0][i] = *(double*)father->get_allele(i, RAND::RandBool());// MOD F.G. 28.04.14
  }
}

//set value
inline void TTQuantiPheno::set_value ()
{
	for(unsigned int i = 0; i < _nb_traits; ++i){
	  _phenotypes[i] = _sequence[0][i] + _sequence[1][i];
	}
  
  // MOD F.G. 28.04.14
  
  // on accede aux variables privees du proto directement:
  if(_myProto->_eVariance != 0)
    for(unsigned int i = 0; i < _nb_traits; ++i)
      _phenotypes[i] += RAND::Gaussian(_myProto->_eVariance);
}

//get genotype

double TTQuantiPheno::get_genotype (unsigned int trait)
{
  // MOD F.G. 28.04.14
  
  // attention, ici il faut retourner la valeur d'un seul trait
  
  //double genotype = 0;
  
  //for(unsigned int j = 0; j < _nb_traits; ++j) {
    
  //  genotype = (_sequence[0][trait] + _sequence[1][trait])/2;
  //}
  
  return _sequence[0][trait] + _sequence[1][trait];
}

//Mutate
inline void TTQuantiPheno::mutate_Gaussian()
{
  
  
  //on tire le nombre de mutation
  unsigned int mut_allele;
  double *effects;
  int num_mutation = (int)RAND::Poisson(_myProto->_genomic_mutation_rate);
  
  for(int k=0; k<num_mutation; k++){
    
    mut_allele = RAND::RandBool();
    
    // MOD F.G. 28.04.14
    //on tire les effets une seule fois, doit etre en-dehors de la boucle for
    effects = _myProto->getMutationEffects();
    
    for(unsigned int i = 0; i < _nb_traits; i++){
      _sequence[mut_allele][i]  += effects[i];
    }
  }
}
//check to go back in the header

//SET MUTATION PTR

//MUTATE


void TTQuantiPheno::set_mutation_fptr (double* (TProtoQuantiPheno::*val) (void))
{
  _mutationFuncPtr = &TTQuantiPheno::mutate_Gaussian;
}


void TTQuantiPheno::mutate ()
{	
  //if(RAND::Uniform() <= _myProto->_genomic_mutation_rate)
    (this->*_mutationFuncPtr)();
}



//show up

void TTQuantiPheno::show_up  ()
{
  message("\
          Trait's type: PHENO\n\
          traits: %i\n\
          loci: %i\n",_nb_traits);
  
  for(unsigned int i = 0; i < _nb_traits; i++)
    message("phenotype %i: %f\n",i+1,_phenotypes[i]);
  
  message("_sequence: \n0:");
  
  for(unsigned int i = 0; i < _nb_traits; ++i) {
    message("\nt%i:",i+1);
    
    message("%.3f,",_sequence[0][i]);
  }
  message("\n1:");
  for(unsigned int i = 0; i < _nb_traits; ++i) {
    message("%.3f,",_sequence[1][i]);
  }
  message("\n");
}

// MOD F.G. 28.04.14
//CLASS//TTQuantiPhenoFH


void TTQuantiPhenoFH::FHwrite()
{
  
  Metapop* pop = get_pop_ptr();
  
  if (!pop->isAlive()) return;
  
  std::string filename = get_filename();
  
  std::ofstream FILE (filename.c_str(), ios::out);
  
  if(!FILE) fatal("could not open \"%s\" output file!!\n",filename.c_str());
  
  bool print_genotype = (_FHLinkedTrait->get_env_var() != 0);
  
  FILE<<"pop ";
  for(unsigned int k = 0; k < _FHLinkedTrait->get_nb_traits(); k++)
    FILE<<"t"<<k+1<<"a1 "<<"t"<<k+1<<"a2 ";
  
  for(unsigned int k = 0; k < _FHLinkedTrait->get_nb_traits(); k++) {
    FILE<<"P"<<k+1<< " "; 
    if(print_genotype) FILE<<"G"<<k+1<< " ";
  }
  
  FILE<<"age sex home ped isMigrant father mother ID\n";
  
  for(unsigned int agex = 0; agex < pop->getNumAgeClasses(); ++agex) {
    print(FILE, static_cast<age_idx> (agex), print_genotype); 
  }
  
  FILE.close();  
}
// ----------------------------------------------------------------------------------------
// print
// ----------------------------------------------------------------------------------------
void TTQuantiPhenoFH::print(ofstream& FH, age_idx Ax, bool print_genotype)
{
  Metapop* pop = get_pop_ptr();
  int patchNbr = pop->getPatchNbr();
  Patch* current_patch;
  Individual* ind;
  TTQuantiPheno* trait;
  double* Tval;
  double **genes;
  unsigned int nb_trait = _FHLinkedTrait->get_nb_traits();
  
  for(int i = 0; i < patchNbr; i++) {
    
    current_patch = pop->getPatch(i);
    
    //print females
    for(unsigned int j = 0, size = current_patch->size(FEM, Ax); j < size; j++) {
      
      ind = current_patch->get(FEM, Ax, j);
      trait = dynamic_cast<TTQuantiPheno*> (ind->getTrait(_FHLinkedTraitIndex));
      
      FH<<i+1<<" "; //patch number
      
      Tval = (double*)trait->getValue();

      genes = (double**)trait->get_sequence();
      
      FH.precision(6);
      
      for(unsigned int k = 0; k < nb_trait; k++) {
          FH<<genes[0][k]<<" "<<genes[1][k]<<" ";
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
    
    //print males
    for(unsigned int j = 0, size = current_patch->size(MAL, Ax); j < size; j++) {
      
      ind = current_patch->get(MAL, Ax, j);
      trait = dynamic_cast<TTQuantiPheno*> (ind->getTrait(_FHLinkedTraitIndex));
      
      FH<<i+1<<" ";
      
      Tval = (double*)trait->getValue();
      
      genes = (double**)trait->get_sequence();
      
      FH.precision(6);
      for(unsigned int k = 0; k < nb_trait; k++) {
        FH<<genes[0][k]<<" "<<genes[1][k]<<" ";
      }
      
      FH.precision(4);
      for(unsigned int k = 0; k < nb_trait; k++) {
        FH<<Tval[k]<<" ";
        if(print_genotype) FH << trait->get_genotype(k) << " ";
      }
      
      FH<<Ax<<" "<<ind->getSex()<<" "<<ind->getHome()+1<<" "<<ind->getPedigreeClass(ind->getMother(), ind->getFather())<<" "
      << (ind->getFather() && ind->getMother() ?
          (ind->getFather()->getHome()!=i) + (ind->getMother()->getHome()!=i) : 0)<<" "
//	  <<ind->getFatherID()<<" "<<ind->getMotherID()<<" "<<ind->getID()<<std::endl;//MOD SLIM NEMO 7.11
          <<"NA NA "<<ind->getID()<<std::endl;
    }
  }
}

// MOD F.G. 28.04.14
//CLASS//TTQuantiPhenoSH


void TTQuantiPhenoSH::resetPtrs ( )
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
void TTQuantiPhenoSH::init()
{  
  StatHandlerBase::init();
  
  _eVar = (_SHLinkedTrait->get_env_var() != 0);
  
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
bool TTQuantiPhenoSH::setStatRecorders(std::string& token)
{
#ifdef _DEBUG_
  message("-TTQuantiPhenoSH::setStatRecorders ");
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
  
  if(sub_token == "pheno") {
    addQuanti(AGE);
  } else if(sub_token == "pheno.eigen") {  //based on Vb; among-patch (D) matrix
    addEigen(AGE);
  } else if(sub_token == "pheno.eigenvalues") {  //based on Vb; among-patch (D) matrix
    addEigenValues(AGE);
  } else if(sub_token == "pheno.eigenvect1") {  //based on Vb; among-patch (D) matrix
    addEigenVect1(AGE);
  } else if(sub_token == "pheno.patch") {
    addQuantiPerPatch(AGE);
  } else if(sub_token == "pheno.mean.patch") {
    addAvgPerPatch(AGE);
  } else if(sub_token == "pheno.var.patch") {
    addVarPerPatch(AGE);
  } else if(sub_token == "pheno.covar.patch") {
    addCovarPerPatch(AGE);
  } else if(sub_token == "pheno.eigen.patch") {
    addEigenPerPatch(AGE);
  } else if(sub_token == "pheno.eigenvalues.patch") {
    addEigenValuesPerPatch(AGE);
  } else if(sub_token == "pheno.eigenvect1.patch") {
    addEigenVect1PerPatch(AGE);
  } else if(sub_token == "pheno.skew.patch") {
    addSkewPerPatch(AGE);
  } else
    return false;
  
  return true;
}// ----------------------------------------------------------------------------------------
// addQuanti
// ----------------------------------------------------------------------------------------
void TTQuantiPhenoSH::addQuanti (age_t AGE)
{
  if (AGE == ALL) {
    addQuanti(ADULTS);
    addQuanti(OFFSPRG);
    return;
  }
  
  string suffix = (AGE == ADULTS ? "adlt.":"off.");
  string name = suffix + "T";
  string t1, t2;
  
  void (TTQuantiPhenoSH::* setter) (void) = (AGE == ADULTS ?
                                        &TTQuantiPhenoSH::setAdultStats : &TTQuantiPhenoSH::setOffsprgStats);
  
  add("",suffix + "T1",AGE,0,0,0,&TTQuantiPhenoSH::getMeanPhenot,0,setter);
  
  for(unsigned int i = 1; i < _nb_trait; i++) {
    t1 = tstring::int2str(i+1);
    add("", name + t1,AGE,i,0,0,&TTQuantiPhenoSH::getMeanPhenot,0,0);
  } 
  
  //Va
  for(unsigned int i = 0; i < _nb_trait; i++) {
    t1 = tstring::int2str(i+1);
    add("", name + t1 +".Va",AGE,i,0,0,&TTQuantiPhenoSH::getVa,0,0);
  }
  
  //Vb
  for(unsigned int i = 0; i < _nb_trait; i++) {
    t1 = tstring::int2str(i+1);
    add("", name + t1 +".Vb",AGE,i,0,0,&TTQuantiPhenoSH::getVb,0,0);
  }
  
  //Vp
  if(_eVar) {
    for(unsigned int i = 0; i < _nb_trait; i++) {
      t1 = tstring::int2str(i+1);
      add("", name + t1 +".Vp",AGE,i,0,0,&TTQuantiPhenoSH::getVp,0,0);
    }
  }
  //Qst
  for(unsigned int i = 0; i < _nb_trait; i++) {
    t1 = tstring::int2str(i+1);
    add("", name + t1 +".Qst",AGE,i,0,0,&TTQuantiPhenoSH::getQst,0,0);
  }
  
  if (_nb_trait > 1) {
    unsigned int c = 0;
    for(unsigned int t = 0; t < _nb_trait; ++t) {
      for(unsigned int v = t + 1, cov = (t+1)*10+(v+1) ; v < _nb_trait; ++v) {
        t1 = tstring::int2str(cov);
        add("", name + t1 +".cov",AGE,c++,0,0,&TTQuantiPhenoSH::getCovar,0,0);
      }
    }
  }
}
// ----------------------------------------------------------------------------------------
// addEigen
// ----------------------------------------------------------------------------------------
void TTQuantiPhenoSH::addEigen (age_t AGE)
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
  
  void (TTQuantiPhenoSH::* setter) (void) = (AGE == ADULTS ?
                                        &TTQuantiPhenoSH::setAdultStats : &TTQuantiPhenoSH::setOffsprgStats);
  
  add("", suffix + "T.eval1",AGE,0,0,0,&TTQuantiPhenoSH::getEigenValue, 0, setter); //this one calls the setter
  
  for(unsigned int t = 1; t < _nb_trait; ++t)
    add("", suffix + "T.eval" + tstring::int2str(t+1),AGE,t,0,0,&TTQuantiPhenoSH::getEigenValue, 0, 0);
  
  for(unsigned int t = 0; t< _nb_trait; ++t)
    for(unsigned int v = 0; v < _nb_trait; ++v)
      add("", suffix + "T.evect" + tstring::int2str((t+1)*10+(v+1)),AGE,t,v,0,0,&TTQuantiPhenoSH::getEigenVectorElt,0);
  
}
// ----------------------------------------------------------------------------------------
// addEigenValues
// ----------------------------------------------------------------------------------------
void TTQuantiPhenoSH::addEigenValues (age_t AGE)
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
  
  void (TTQuantiPhenoSH::* setter) (void) = (AGE == ADULTS ?
                                        &TTQuantiPhenoSH::setAdultStats : &TTQuantiPhenoSH::setOffsprgStats);
  
  add("", suffix + "T.eval1",AGE,0,0,0,&TTQuantiPhenoSH::getEigenValue, 0, setter);
  
  for(unsigned int t = 1; t < _nb_trait; ++t)
    add("", suffix + "T.eval" + tstring::int2str(t+1),AGE,t,0,0,&TTQuantiPhenoSH::getEigenValue, 0, 0);
  
  
}
// ----------------------------------------------------------------------------------------
// addEigenVect1 : save only the first eigenvector
// ----------------------------------------------------------------------------------------
void TTQuantiPhenoSH::addEigenVect1 (age_t AGE)
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
  
  void (TTQuantiPhenoSH::* setter) (void) = (AGE == ADULTS ?
                                        &TTQuantiPhenoSH::setAdultStats : &TTQuantiPhenoSH::setOffsprgStats);
  
  add("", suffix + "T.evect11",AGE,0,0,0,0,&TTQuantiPhenoSH::getEigenVectorElt,setter);
  
  for(unsigned int v = 1; v < _nb_trait; ++v)
    add("", suffix + "T.evect1" + tstring::int2str(v+1),AGE,0,v,0,0,&TTQuantiPhenoSH::getEigenVectorElt,0);
  
}
// ----------------------------------------------------------------------------------------
// addQuantiPerPatch
// ----------------------------------------------------------------------------------------
void TTQuantiPhenoSH::addQuantiPerPatch (age_t AGE)
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
void TTQuantiPhenoSH::addAvgPerPatch (age_t AGE)
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
  
  void (TTQuantiPhenoSH::* setter) (void) = (AGE == ADULTS ?
                                        &TTQuantiPhenoSH::setAdultStats : &TTQuantiPhenoSH::setOffsprgStats);
  
  add("Mean phenotype of trait 1 in patch 1", suffix + "T1.p1",  AGE, 0, 0, 0, 0,
      &TTQuantiPhenoSH::getMeanPhenotPerPatch, setter);
  
  for(unsigned int p = 0; p < patchNbr; p++) {
    for(unsigned int i = 0; i < _nb_trait; i++) {
      if(p == 0 && i == 0) continue;
      name = "Mean phenotype of trait " + tstring::int2str(i+1) + " in patch " + tstring::int2str(p+1);
      t1 = "T" + tstring::int2str(i+1);
      patch = ".p" + tstring::int2str(p+1);
      add(name, suffix + t1 + patch, AGE, i, p, 0, 0, &TTQuantiPhenoSH::getMeanPhenotPerPatch, 0);
    } 
  }
  
}
// ----------------------------------------------------------------------------------------
// addVarPerPatch
// ----------------------------------------------------------------------------------------
void TTQuantiPhenoSH::addVarPerPatch (age_t AGE)
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
  
  void (TTQuantiPhenoSH::* setter) (void) = (AGE == ADULTS ?
                                        &TTQuantiPhenoSH::setAdultStats :
                                        &TTQuantiPhenoSH::setOffsprgStats);
  
  add("Genetic variance of trait 1 in patch 1", suffix + "Va.T1.p1",  AGE, 0, 0, 0, 0,
      &TTQuantiPhenoSH::getVaPerPatch, setter);
  
  for(unsigned int p = 0; p < patchNbr; p++) {
    for(unsigned int i = 0; i < _nb_trait; i++) {
      if(p == 0 && i == 0) continue;
      name = "Genetic variance of trait " + tstring::int2str(i+1) + " in patch " + tstring::int2str(p+1);
      t1 = "T" + tstring::int2str(i+1);
      patch = ".p" + tstring::int2str(p+1);
      add(name, suffix + "Va." + t1 + patch, AGE, i, p, 0, 0, &TTQuantiPhenoSH::getVaPerPatch, 0);
    }
  }
  
  if(_eVar) {
    for(unsigned int p = 0; p < patchNbr; p++) {
      for(unsigned int i = 0; i < _nb_trait; i++) {
        name = "Phenotypic variance of trait " + tstring::int2str(i+1) + " in patch " + tstring::int2str(p+1);
        t1 = "T" + tstring::int2str(i+1);
        patch = ".p" + tstring::int2str(p+1);
        add(name, suffix + "Vp." + t1 + patch, AGE, i, p, 0, 0, &TTQuantiPhenoSH::getVpPerPatch, 0);
      }
    }
  }
}
// ----------------------------------------------------------------------------------------
// addCovarPerPatch
// ----------------------------------------------------------------------------------------
void TTQuantiPhenoSH::addCovarPerPatch (age_t AGE)
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
  
  void (TTQuantiPhenoSH::* setter) (void) = (AGE == ADULTS ?
                                        &TTQuantiPhenoSH::setAdultStats : &TTQuantiPhenoSH::setOffsprgStats);
  
  add("Genetic covariance of trait 1 and trait 2 in patch 1", suffix + "cov.T12.p1",  AGE, 0, 0, 0, 0,
      &TTQuantiPhenoSH::getCovarPerPatch, setter);
  
  unsigned int c;
  for(unsigned int p = 0; p < patchNbr; p++) {
    patch = ".p" + tstring::int2str(p+1);
    c = 0;
    for(unsigned int t = 0; t < _nb_trait; ++t) {
      for(unsigned int v = t + 1; v < _nb_trait; ++v){
        if(p==0 && t==0 && v==1) {c++; continue;}
        cov = tstring::int2str((t+1)*10+v+1);
        add("", suffix + "cov.T" + cov + patch,  AGE, p, c++, 0, 0, &TTQuantiPhenoSH::getCovarPerPatch, 0);
      }
    }
  }
  
}
// ----------------------------------------------------------------------------------------
// addEigenPerPatch
// ----------------------------------------------------------------------------------------
void TTQuantiPhenoSH::addEigenPerPatch (age_t AGE)
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
  
  void (TTQuantiPhenoSH::* setter) (void) = (AGE == ADULTS ?
                                        &TTQuantiPhenoSH::setAdultStats : &TTQuantiPhenoSH::setOffsprgStats);
  
  
  add("First G-matrix eigenvalue in patch 1", suffix + "Teval1.p1",  AGE, 0, 0, 0, 0,
      &TTQuantiPhenoSH::getEigenValuePerPatch, setter);
  
  for(unsigned int p = 0; p < patchNbr; ++p) {
    patch = ".p" + tstring::int2str(p+1);
    for(unsigned int t = 0; t < _nb_trait; ++t) {
      if(p==0 && t==0) continue;
      add("", suffix + "Teval" + tstring::int2str(t+1) + patch,  AGE, t, p, 0, 0, &TTQuantiPhenoSH::getEigenValuePerPatch,0);
    }
  }
  for(unsigned int p = 0; p < patchNbr; ++p) {
    patch = ".p" + tstring::int2str(p+1);
    pv = 0;
    for(unsigned int t = 0; t < _nb_trait; ++t)
      for(unsigned int v = 0; v < _nb_trait; ++v)
        add("", suffix + "Tevect" + tstring::int2str((t+1)*10+v+1) + patch,  AGE, p, pv++, 0, 0, &TTQuantiPhenoSH::getEigenVectorEltPerPatch,0);
  }
  
}
// ----------------------------------------------------------------------------------------
// addEigenValuesPerPatch
// ----------------------------------------------------------------------------------------
void TTQuantiPhenoSH::addEigenValuesPerPatch (age_t AGE)
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
  
  void (TTQuantiPhenoSH::* setter) (void) = (AGE == ADULTS ?
                                        &TTQuantiPhenoSH::setAdultStats : &TTQuantiPhenoSH::setOffsprgStats);
  
  add("First G-matrix eigenvalue in patch 1", suffix + "Teval1.p1",  AGE, 0, 0, 0, 0,
      &TTQuantiPhenoSH::getEigenValuePerPatch, setter);
  
  for(unsigned int p = 0; p < patchNbr; ++p) {
    patch = ".p" + tstring::int2str(p+1);
    for(unsigned int t = 0; t < _nb_trait; ++t) {
      if(p==0 && t==0) continue;
      add("", suffix + "Teval" + tstring::int2str(t+1) + patch,  AGE, t, p, 0, 0, &TTQuantiPhenoSH::getEigenValuePerPatch,0);
    }
  }  
}
// ----------------------------------------------------------------------------------------
// addEigenVect1PerPatch
// ----------------------------------------------------------------------------------------
void TTQuantiPhenoSH::addEigenVect1PerPatch (age_t AGE)
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
  
  void (TTQuantiPhenoSH::* setter) (void) = (AGE == ADULTS ?
                                        &TTQuantiPhenoSH::setAdultStats : 
                                        &TTQuantiPhenoSH::setOffsprgStats);
  
  
  add("First G-matrix eigenvector in patch 1", suffix + "Tevect11.p1",  AGE, 0, 0, 0, 0,
      &TTQuantiPhenoSH::getEigenVectorEltPerPatch, setter);
  
  for(unsigned int p = 0; p < patchNbr; ++p) {
    patch = ".p" + tstring::int2str(p+1);
    pv = 0;
    //    for(unsigned int t = 0; t < _nb_trait; ++t)
    for(unsigned int v = 0; v < _nb_trait; ++v){
      if(p==0 && v==0) {pv++; continue;}
      add("", suffix + "Tevect1" + tstring::int2str(v+1) + patch,  AGE, p, pv++, 0, 0, &TTQuantiPhenoSH::getEigenVectorEltPerPatch,0);
    }
  }
}
// ----------------------------------------------------------------------------------------
// addSkewPerPatch
// ----------------------------------------------------------------------------------------
void TTQuantiPhenoSH::addSkewPerPatch(age_t AGE)
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
  
  void (TTQuantiPhenoSH::* setter) (void) = (AGE == ADULTS ?
                                        &TTQuantiPhenoSH::setAdultStats :
                                        &TTQuantiPhenoSH::setOffsprgStats);
  
  add("Genetic skew of trait 1 in patch 1", suffix + "Sk.T1.p1",  AGE, 0, 0, 0, 0,
      &TTQuantiPhenoSH::getSkewPerPatch, setter);
  
  for(unsigned int p = 0; p < patchNbr; p++) {
    for(unsigned int i = 0; i < _nb_trait; i++) {
      if(p == 0 && i == 0) continue;
      name = "Genetic skew of trait " + tstring::int2str(i+1) + " in patch " + tstring::int2str(p+1);
      t1 = "T" + tstring::int2str(i+1);
      patch = ".p" + tstring::int2str(p+1);
      add(name, suffix + "Sk." + t1 + patch, AGE, i, p, 0, 0, &TTQuantiPhenoSH::getSkewPerPatch, 0);
    }
  }
  
}
// ----------------------------------------------------------------------------------------
// getSkewPerPatch
// ----------------------------------------------------------------------------------------
double TTQuantiPhenoSH::getSkewPerPatch (unsigned int i, unsigned int p)
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
void TTQuantiPhenoSH::setDataTables(age_t AGE) 
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
//  age_idx age = (AGE == ADULTS ? ADLTx : OFFSx);
  
  for(unsigned int i = 0, n; i < nb_patch; i++) {
    
    patch = _pop->getPatch(i);
    
    n=0;
    
    if ((patch->size(MAL, AGE)+patch->size(FEM, AGE)) != _phenoTable.size(0,i)) {
      fatal("problem while recording pheno trait values; table size doesn't match patch size.\n");
    }
    store_pheno_trait_values(patch, i, patch->size(MAL, AGE), &n, MAL, AGE, &_phenoTable, &_genoTable,
                       _nb_trait, _SHLinkedTraitIndex);
    
    store_pheno_trait_values(patch, i, patch->size(FEM, AGE), &n, FEM, AGE, &_phenoTable, &_genoTable,
                       _nb_trait, _SHLinkedTraitIndex);
    
    if (n != _phenoTable.size(0,i) || n != _genoTable.size(0,i)) {
      fatal("problem while recording pheno trait values; size counter doesn't match table size.\n");
    }
  }
}
// ----------------------------------------------------------------------------------------
// store_trait_values
// ----------------------------------------------------------------------------------------
void store_pheno_trait_values (Patch* patch, unsigned int patchID, unsigned int size, unsigned int *cntr,
                         sex_t SEX, age_t AGE, DataTable<double> *ptable, DataTable<double> *gtable,
                         unsigned int nTrait, unsigned int TraitIndex)
{
  double *phe;
  TTQuantiPheno *trait;
  
  for(unsigned int j = 0; j < size; ++j) {
    
    trait = dynamic_cast<TTQuantiPheno*> (patch->get(SEX, AGE, j)->getTrait( TraitIndex ));
    
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
void TTQuantiPhenoSH::setStats (age_t AGE)
{  
  if(_table_set_age == AGE 
     && _table_set_gen == _pop->getCurrentGeneration()
     && _table_set_repl == _pop->getCurrentReplicate())
    return;
  
  unsigned int pop_size = _pop->size(AGE);
  unsigned int nb_patch = _pop->getPatchNbr();
  unsigned int nb_extant_patch = 0;
  unsigned int patch_size;
  double *phenot1, *genot1, *genot2;
  
  setDataTables(AGE);
  
#ifdef HAS_GSL
  unsigned int c = 0; //covariance position counter
  unsigned int pv = 0; //eigenvector position counter
  //within deme stats:
  for(unsigned int j = 0; j < nb_patch; j++) { 
    
    patch_size = _pop->size(AGE, j);

    if (patch_size > 0) nb_extant_patch++;

    for(unsigned int t=0; t < _nb_trait; t++) {
      
      phenot1 = _phenoTable.getClassWithinGroup(t,j);
      genot1  = _genoTable.getClassWithinGroup(t,j);
      
      _pmeanP[t][j] = my_mean (phenot1, patch_size ); //this may be NaN if patch_size = 0
      _pmeanG[t][j] = my_mean (genot1, patch_size );  //this may be NaN if patch_size = 0
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
  
  // ------------------------ among demes stats:
  for(unsigned int t1 = 0; t1 < _nb_trait; t1++) {
    
    phenot1 = _phenoTable.getGroup(t1);
    genot1  = _genoTable.getGroup(t1);
    
    _meanP[t1] = my_mean (phenot1, pop_size ); //grand mean, pooling all individuals
    _meanG[t1] = my_mean (genot1, pop_size ); //grand mean, pooling all individuals
    
    _Vp[t1] = my_mean_no_nan (_pVp[t1], nb_patch); //mean within-patch variance (extant demes only)
    _Va[t1] = my_mean_no_nan (_pVa[t1], nb_patch); //mean within-patch variance (extant demes only)
    
    meanGamong1 = my_mean_no_nan (_pmeanG[t1], nb_patch); //mean of within-patch mean genotypic values
    
    _Vb[t1] = my_variance_with_fixed_mean_no_nan (_pmeanG[t1], nb_patch,  meanGamong1); //variance of patch means
    
    gsl_matrix_set(_G, t1, t1, _Vb[t1]); //_G here becomes the D-matrix, the among-deme (Difference) covariance matrix 
    
    for(unsigned int t2 = t1 + 1; t2 < _nb_trait; t2++) {
      
      meanGamong2 = my_mean_no_nan (_pmeanG[t2], nb_patch);
      
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
