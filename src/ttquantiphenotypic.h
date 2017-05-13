/*
 * ttquantiphenotypic.h
 *
 *  Created on: 21.04.2014
 *      Author: cotto_o
 */

#ifndef TTQUANTIPHENOTYPIC_H
#define TTQUANTIPHENOTYPIC_H


#include <cmath>
#include "filehandler.h"
#include "stathandler.h"
#include "metapop.h"
#include "datatable.h"
#include <utility>
#include <string>
#include <vector>
#ifdef HAS_GSL
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#endif


class TProtoQuantiPheno;

class TTQuantiPheno : public TTrait {

public:

//Constructor
TTQuantiPheno():
	_sequence(0),
	_phenotypes(0),
	_nb_traits(0),
	_init_value(0),
	_myProto(0) 
  // MOD F.G. 28.04.14
  //j'ai elimine toutes les variables inutiles qui sont deja dans le proto:
//	_genomic_mutation_rate(0),
//	_eVariance(0),
//	_mutation_variance(0)
{}

//Copy Constructor
TTQuantiPheno(const TTQuantiPheno &T):
	_sequence(0),
  _nb_traits(T._nb_traits),
  _myProto(T._myProto),
  // il manque la copie du pointeur sur la fonction de mutation
  // c'est necessaire car les nouveau individus sont cree en 
  // appelant ce constructeur-ci (suite appel a TTQuantiPheno::clone())
  _mutationFuncPtr(T._mutationFuncPtr)
  // MOD F.G. 28.04.14
  //j'ai elimine toutes les variables inutiles qui sont deja dans le proto:
//	_genomic_mutation_rate(T._genomic_mutation_rate),
//  _eVariance(T._eVariance),
//  _mutation_variance(T._mutation_variance)

{
  _phenotypes = new double [_nb_traits];
  _init_value = new double [_nb_traits];

  for (unsigned int i = 0; i < _nb_traits; ++i) 
    _init_value[i] = T._init_value[i];
}

//Destructor
~TTQuantiPheno(){
  
//  cout << "~TTQuantiPheno" << endl;
//  cout << "_phenotypes " << _phenotypes << endl;
//  cout << "_init_value " << _init_value << endl;
//  cout << "_sequence " << _sequence << endl;
//  if(_sequence) cout << "_sequence[0] " << _sequence[0] << endl;
//  if(_sequence) cout << "_sequence[1] " << _sequence[1] << endl;
  
  reset();
  
//  cout << "_phenotypes " << _phenotypes << endl;
//  cout << "_init_value " << _init_value << endl;
//  cout << "_sequence " << _sequence << endl;
}

//Implement TTrait
virtual   void            init ();
virtual	  void			      init_sequence();
virtual   void            reset ();
virtual   void            inherit (TTrait* mother, TTrait* father);
virtual   void            mutate ();
virtual   void*           set_trait (void* value)                 {return value;}
virtual   void            set_sequence(void **seq) {};
virtual   void            set_value ();
virtual   void*           getValue () const                       {return _phenotypes;}// MOD F.G. 28.04.14
virtual   trait_t         get_type () const                       {return PHENO;}
virtual 	void*			      get_allele(int loc, int allele) const   {return (void*)&_sequence[allele][loc];};
virtual		void**			    get_sequence() const                   	{return (void**)_sequence;}// MOD F.G. 28.04.14
virtual   void            show_up  ();
virtual   TTQuantiPheno*  clone ()                           		  {return new TTQuantiPheno(*this);}
virtual   TTQuantiPheno&  operator= (const TTrait& T);
virtual   bool            operator== (const TTrait& T);
virtual   bool            operator!= (const TTrait& T);

double get_genotype(unsigned int trait);

  //implements StorableComponent:
  virtual void store_data    ( BinaryStorageBuffer* saver  )
  {
    saver->store(_sequence[0],_nb_traits*sizeof(double));
    saver->store(_sequence[1],_nb_traits*sizeof(double));
    saver->store(_phenotypes,_nb_traits*sizeof(double));
  }
  virtual bool retrieve_data ( BinaryStorageBuffer* reader )
  {
    reader->read(_sequence[0],_nb_traits*sizeof(double));
    reader->read(_sequence[1],_nb_traits*sizeof(double));
    reader->read(_phenotypes,_nb_traits*sizeof(double));
    return true;
  }
  
  void mutate_Gaussian   ();
  
  //Accessors
  void set_proto                  (TProtoQuantiPheno* proto) {_myProto = proto;}
  void set_nb_traits              (unsigned int val)    {_nb_traits = val;}
  void set_init_value             (double* val);
  void set_mutation_fptr          (double* (TProtoQuantiPheno::*val) (void));
  // MOD F.G. 28.04.14
  //j'ai elimine toutes les variables inutiles qui sont deja dans le proto:
//void set_genomic_mutation_rate  (double val)          {_genomic_mutation_rate = val;}
//void set_mutation_variance		(double val) {_mutation_variance = val;}
//void set_eVariance (double val)  {_eVariance = val;}

private:

  double **_sequence;
  double *_phenotypes;
  double *_init_value;
  TProtoQuantiPheno* _myProto;
  void (TTQuantiPheno::* _mutationFuncPtr) (void);
  unsigned int _nb_traits;
  
  // MOD F.G. 28.04.14
  //j'ai elimine toutes les variables inutiles qui sont deja dans le proto:
  
  //double _genomic_mutation_rate;
  //unsigned int _doInitMutation;
  //double _eVariance;
  //double _mutation_variance;
  //double* (TProtoQuantiPheno::* _getMutationValues) (unsigned int);
  
};

// MOD F.G. 28.04.14
class TTQuantiPhenoSH;
class TTQuantiPhenoFH;
// juste quelques 'forward declarations' pour les outputs

class TProtoQuantiPheno: public TraitPrototype  {

public:

	//Constructor-destructor
	TProtoQuantiPheno ();//constructor
	TProtoQuantiPheno (const TProtoQuantiPheno& T);//constructor copy
	~TProtoQuantiPheno ();//destructor


	unsigned int get_nb_traits() {return _nb_traits;}
	double       get_env_var () {return _eVariance;}
	double       get_trait_var (unsigned int trait) {return _mutation_matrix->get(trait, trait);}
	double       get_mutation_correlation() {return _mutation_correlation;}


	void reset_mutation_pointers();
	bool setMutationParameters ();
	bool setMutationModel();
	void set_mutation_matrix_decomposition ();

	// MOD F.G. 28.04.14
	double* getMutationEffects () {return (this->*_mutation_func_ptr) ();}
	//ces fonctions n'ont plus besoin d'argument, le taux de mutation est pour le trait
	double* getMutationEffectMultivariateGaussian ();
	double* getMutationEffectBivariateGaussian    ();
	double* getMutationEffectUnivariateGaussian   ();
	// end MOD

	//implementation of trait Prototype
	virtual void init(){};
	virtual void reset(){};
	virtual TTQuantiPheno* hatch();
	virtual TProtoQuantiPheno* clone() {return new TProtoQuantiPheno(*this);}
	virtual trait_t get_type() const {return PHENO;}

	//implementation of StorableComponent:
	virtual void store_data    ( BinaryStorageBuffer* saver  )
	{saver->store(&_nb_traits,sizeof(int));}
	virtual bool retrieve_data ( BinaryStorageBuffer* reader )
	{reader->read(&_nb_traits,sizeof(int));return true;};

	//implementation of SimComponent:
	virtual bool setParameters();
	virtual void loadFileServices ( FileServices* loader );
	virtual void loadStatServices ( StatServices* loader );

private:

  TMatrix *_mutation_matrix;
  gsl_matrix *_gsl_mutation_matrix, *_evect;
  gsl_vector *_eval, *_effects_multivar;
  gsl_vector *_ws;
  double  _mutation_correlation;
  double *_mutation_sigma, *_init_value, _effects_bivar[2];
  
  unsigned int _init_model;

  unsigned int _nb_traits;
  double _genomic_mutation_rate, _mutation_variance; //il faut decider comment fixer Vm
  
  double* (TProtoQuantiPheno::* _mutation_func_ptr) (void);
  
  double _eVariance;
  
  friend class TTQuantiPheno; //ici on permet a TTQuantiPheno l'acces aux variables privees...
  
  // MOD F.G. 28.04.14
  TTQuantiPhenoSH* _stats;
  TTQuantiPhenoFH* _writer;

};

// MOD F.G. 28.04.14
class TTQuantiPhenoFH : public TraitFileHandler<TProtoQuantiPheno> {
  
public:
  TTQuantiPhenoFH(TProtoQuantiPheno* T) : TraitFileHandler<TProtoQuantiPheno>(T,".pheno") {}
  virtual ~TTQuantiPhenoFH(){}
   
  virtual void  FHwrite ();
  void print(ofstream& FH, age_idx Ax, bool print_genotype);
  virtual void FHread (string& filename) {}
};


// MOD F.G. 28.04.14
class TTQuantiPhenoSH : public TraitStatHandler<TProtoQuantiPheno, TTQuantiPhenoSH> {
  
  double *_meanP, *_meanG, *_Va, *_Vb, *_Vp, *_covar,*_eigval,**_eigvect;
  double **_pmeanP, **_pmeanG, **_pVa, **_pVp, **_pcovar, **_peigval, **_peigvect;
  
  unsigned int _nb_trait, _patchNbr;
  bool _eVar;
  
  gsl_matrix *_G, *_evec;
  gsl_vector *_eval;
  gsl_eigen_symmv_workspace *_ws;
  
  DataTable< double > _phenoTable, _genoTable;
  unsigned int _table_set_gen, _table_set_age, _table_set_repl;
  
public:
  
  TTQuantiPhenoSH(TProtoQuantiPheno* TP) 
  : TraitStatHandler<TProtoQuantiPheno, TTQuantiPhenoSH> (TP), 
  _meanP(0), _meanG(0), _Va(0), _Vb(0), _Vp(0), _covar(0), _eigval(0), _eigvect(0),
  _pmeanP(0), _pmeanG(0), _pVa(0), _pVp(0), _pcovar(0), _peigval(0), _peigvect(0),
  _nb_trait(0),_patchNbr(0),_G(0),_evec(0),_eval(0),_ws(0),
  _table_set_gen(999999), _table_set_age(999999), _table_set_repl(999999)
  {}
  
  virtual ~TTQuantiPhenoSH() {resetPtrs();}
  
  void  resetPtrs();
  
  virtual void init ( );
  
  virtual bool      setStatRecorders (std::string& token);
  void addQuanti (age_t AGE);
  void addEigen (age_t AGE);
  void addEigenValues (age_t AGE);
  void addEigenVect1 (age_t AGE);
  void addQuantiPerPatch (age_t AGE);
  void addAvgPerPatch (age_t AGE);
  void addVarPerPatch (age_t AGE);
  void addCovarPerPatch (age_t AGE);
  void addEigenPerPatch (age_t AGE);
  void addEigenValuesPerPatch (age_t AGE);
  void addEigenVect1PerPatch (age_t AGE);
  void addEigenStatsPerPatcg (age_t AGE);
  void addSkewPerPatch(age_t AGE);
  
  void   setDataTables             (age_t AGE);
  void   setAdultStats             ( ) {setStats(ADULTS);}
  void   setOffsprgStats           ( ) {setStats(OFFSPRG);}
  void   setStats                  (age_t AGE);
  double getMeanPhenot             (unsigned int i) {return _meanP[i];}
  double getVa                     (unsigned int i) {return _Va[i];}
  double getVb                     (unsigned int i) {return _Vb[i];}
  double getVp                     (unsigned int i) {return _Vp[i];}
  double getQst                    (unsigned int i) {return _Vb[i]/(_Vb[i]+2*_Va[i]);}
  double getCovar                  (unsigned int i) {return _covar[i];}
  double getEigenValue             (unsigned int i) {return _eigval[i];}
  double getEigenVectorElt         (unsigned int t1, unsigned int t2) {return _eigvect[t2][t1];}//eigenvectors arranged column-wise
  
  double getMeanPhenotPerPatch     (unsigned int i, unsigned int p) {return _pmeanP[i][p];}
  double getVaPerPatch             (unsigned int i, unsigned int p) {return _pVa[i][p];}
  double getVpPerPatch             (unsigned int i, unsigned int p) {return _pVp[i][p];}
  double getEigenValuePerPatch     (unsigned int i, unsigned int p) {return _peigval[i][p];}
  double getCovarPerPatch          (unsigned int p, unsigned int i) {return _pcovar[p][i];}
  double getEigenVectorEltPerPatch (unsigned int p, unsigned int v) {return _peigvect[p][v];}
  double getSkewPerPatch           (unsigned int i, unsigned int p);
};

#endif /* TTQUANTIPHENOTYPIC_H */
