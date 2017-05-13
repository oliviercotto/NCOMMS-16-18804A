/** $Id: ttquanti.h,v 1.11.2.4 2015-08-06 06:46:07 fred Exp $
*
*  @file ttquanti.h
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

#ifndef TTQUANTI_H
#define TTQUANTI_H

#include <cmath>
#include <vector>
#include "ttrait_with_map.h"
#include "filehandler.h"
#include "stathandler.h"
#include "metapop.h"
#include "datatable.h"
#ifdef HAS_GSL
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#endif

class TTQuantiSH;
class TTQuantiFH;
class TTQFreqExtractor;
class TProtoQuanti;

// ------------------------------------------------------------------------------
/**
 *  TTQuanti
 */
// ------------------------------------------------------------------------------
class TTQuanti: public TTraitWithMap {
  
public:
  
  TTQuanti () : _sequence(0),_phenotypes(0),_nb_locus(0),_nb_traits(0),_seq_length(0),
                _genomic_mutation_rate(0),_init_value(0),_myProto(0),_getMutationValues(0),
                _mutationFuncPtr(0),_doInitMutation(1),_inherit(0),_eVariance(0)
  { }
  
  TTQuanti (const TTQuanti& T)  : _sequence(0),_nb_locus(T._nb_locus),_nb_traits(T._nb_traits),
    _seq_length(T._seq_length),_genomic_mutation_rate(T._genomic_mutation_rate),
    _myProto(T._myProto),_getMutationValues(T._getMutationValues),
    _mutationFuncPtr(T._mutationFuncPtr),_doInitMutation(T._doInitMutation),_inherit(T._inherit),
    _eVariance(T._eVariance)
  {
    _phenotypes = new double [_nb_traits];
    _init_value = new double [_nb_traits];
    for (unsigned int i = 0; i < _nb_traits; ++i) _init_value[i] = T._init_value[i];
  }

  ~TTQuanti () {reset();}
  
  //implements TTrait:
  virtual   void            init ();
  virtual   void            init_sequence ();
  virtual   void            reset ();
  virtual   void            inherit (TTrait* mother, TTrait* father);
  virtual   void            mutate ()                               {(this->*_mutationFuncPtr)();}
  virtual   void*           set_trait (void* value)                 {return value;}
  virtual   void            set_sequence (void** seq)               {reset(); _sequence = (double**)seq;}
  virtual   void            set_value ();
  virtual   void*           getValue () const                       {return _phenotypes;}
  virtual   trait_t         get_type () const                       {return QUANT;}
  virtual   void**          get_sequence () const                   {return (void**)_sequence;}
  virtual   void*           get_allele (int loc, int all) const     {return (loc < (int)_seq_length && all < 2 ? (void*)&_sequence[all][loc] : 0);}
  virtual   void            show_up  ();
  virtual   TTQuanti*       clone ()                                {return new TTQuanti(*this);}
  virtual   TTQuanti&       operator= (const TTrait& T);
  virtual   bool            operator== (const TTrait& T);
  virtual   bool            operator!= (const TTrait& T);
  
  double    get_genotype (unsigned int trait);
  
  //implements StorableComponent:
  virtual void store_data    ( BinaryStorageBuffer* saver  ) 
   {saver->store(_sequence[0],_seq_length*sizeof(double));
    saver->store(_sequence[1],_seq_length*sizeof(double));}
  
  virtual bool retrieve_data ( BinaryStorageBuffer* reader ) 
   {reader->read(_sequence[0],_seq_length*sizeof(double));
    reader->read(_sequence[1],_seq_length*sizeof(double));return true;}
  
  void mutate_noHC ();
  void mutate_HC   ();
  //accessors:
  void set_proto                  (TProtoQuanti* proto) {_myProto = proto;}
  void set_nb_locus               (unsigned int val)    {_nb_locus = val;}
  void set_nb_traits              (unsigned int val)    {_nb_traits = val;}
  void set_seq_length             (unsigned int val)    {_seq_length = val;}
  void set_genomic_mutation_rate  (double val)          {_genomic_mutation_rate = val;}
  void set_init_value             (double* val, unsigned int doInit);
  void set_init_value             (double* val);
  void set_mutation_fptr          (double* (TProtoQuanti::* val) (unsigned int) , bool _isHC) 
  {_getMutationValues = val; _mutationFuncPtr = (_isHC ? &TTQuanti::mutate_HC : &TTQuanti::mutate_noHC); }
  void set_inherit_fptr           (void (TProtoQuanti::* val) (sex_t, double*, double**)) {_inherit = val;}
  void set_eVariance              (vector<double> evar)          {_eVariance = evar;}

  void set_allele                 (int locus, int allele, double value) {_sequence[allele][locus] = value;}
  
private:
    
  double **_sequence;
  double *_phenotypes;
  
  unsigned int _nb_locus, _nb_traits, _seq_length;
  double _genomic_mutation_rate, *_init_value;
  unsigned int _doInitMutation;
  vector<double> _eVariance;

  TProtoQuanti* _myProto;
  double* (TProtoQuanti::* _getMutationValues) (unsigned int);
  void (TProtoQuanti::* _inherit) (sex_t, double*, double**);
  void (TTQuanti::* _mutationFuncPtr) (void);
};

// ------------------------------------------------------------------------------
/**
 *  TProtoQuanti
 */
// ------------------------------------------------------------------------------
class TProtoQuanti : public TTProtoWithMap {

public:
  
  TProtoQuanti ();
  TProtoQuanti (const TProtoQuanti& T);
  virtual ~TProtoQuanti ();
  
  unsigned int get_nb_traits() {return _nb_traits;}
  unsigned int get_nb_locus() {return _nb_locus;}
  unsigned int get_seq_length () {return _seq_length;}
  vector<double>       get_env_var () {return _eVariance;}
  double       get_trait_var (unsigned int trait) {return _mutation_matrix->get(trait, trait);}
  double       get_mutation_correlation() {return _mutation_correlation;}
  unsigned int get_allele_model () {return _allele_model;}
  double **    get_allele_values () const {return _allele_value;}

  void reset_mutation_pointers();
  bool setMutationParameters ();
  bool setDiallelicMutationModel();   //@add from ttdouble
  bool setContinuousMutationModel();  //@add from ttdouble
  void set_mutation_matrix_decomposition ();
  
  double  get_init_value (unsigned int i) { return _init_value[i]; }

  double* getMutationEffectMultivariateGaussian (unsigned int loc);
  double* getMutationEffectBivariateGaussian    (unsigned int loc);
  double* getMutationEffectUnivariateGaussian   (unsigned int loc);
  double* getMutationEffectUnivariateDiallelic  (unsigned int loc); //@add from ttdouble
  double* getMutationEffectBivariateDiallelic   (unsigned int loc); //@add from ttdouble
  
  void inherit_free (sex_t SEX, double* seq, double** parent);
  void inherit_low  (sex_t SEX, double* seq, double** parent);
    
  //implements TraitPrototype:
  virtual   void            reset () {TTProtoWithMap::reset(); }
  virtual   TTQuanti*       hatch ();
  virtual   TProtoQuanti*   clone ()          {return new TProtoQuanti(*this);}
  virtual   trait_t         get_type () const {return QUANT;}
  
  //implementation of SimComponent:
  virtual bool setParameters();
  virtual void loadFileServices ( FileServices* loader );
  virtual void loadStatServices ( StatServices* loader );
  
  //implementation of StorableComponent:
  virtual void store_data    ( BinaryStorageBuffer* saver  ) 
    {saver->store(&_seq_length,sizeof(int));}
  
  virtual bool retrieve_data ( BinaryStorageBuffer* reader ) 
    {reader->read(&_seq_length,sizeof(int));return true;}
  
private:    
  
  unsigned int _nb_locus, _nb_traits, _seq_length;
  unsigned int _allele_model; //@add from ttdouble
  double** _allele_value;    //@add from ttdouble
  
  //mutations:
  TMatrix *_mutation_matrix;
  gsl_matrix *_gsl_mutation_matrix, *_evect;
  gsl_vector *_eval, *_effects_multivar;
  gsl_vector *_ws;
  double _genomic_mutation_rate, _mutation_correlation;
  double *_mutation_sigma, *_init_value, _effects_bivar[2];
  unsigned int _doInitMutation;
  
  //recombination:
  bool* _all_chooser;
  size_t _locusByteSize, _sizeofLocusType;
  
  vector<double> _eVariance;

  double* (TProtoQuanti::* _mutation_func_ptr) (unsigned int);  
  
  friend class TTQuanti;
  
  TTQuantiSH* _stats;
  TTQuantiFH* _writer;
  TTQFreqExtractor* _freqExtractor;
};
// ------------------------------------------------------------------------------
/**
 *  TTQuantiSH
 */
// ------------------------------------------------------------------------------
class TTQuantiSH : public TraitStatHandler<TProtoQuanti, TTQuantiSH> {

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

    TTQuantiSH(TProtoQuanti* TP) 
    : TraitStatHandler<TProtoQuanti, TTQuantiSH> (TP), 
    _meanP(0), _meanG(0), _Va(0), _Vb(0), _Vp(0), _covar(0), _eigval(0), _eigvect(0),
    _pmeanP(0), _pmeanG(0), _pVa(0), _pVp(0), _pcovar(0), _peigval(0), _peigvect(0),
    _nb_trait(0),_patchNbr(0),_G(0),_evec(0),_eval(0),_ws(0),
    _table_set_gen(999999), _table_set_age(999999), _table_set_repl(999999)
    {}
  
  virtual ~TTQuantiSH() {resetPtrs();}
  
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


class TTQuantiFH : public TraitFileHandler<TProtoQuanti> {

  string _output_option;
  
public:
  TTQuantiFH(TProtoQuanti* T) : TraitFileHandler<TProtoQuanti>(T,".quanti") {}
  virtual ~TTQuantiFH(){}
  
  void setOutputOption (string opt) {_output_option = opt;}
  
  virtual void  FHwrite ();
  void print(ofstream& FH, age_idx Ax, bool print_gene, bool print_genotype);
  virtual void FHread (string& filename) {}
};

class TTQFreqExtractor : public TraitFileHandler<TProtoQuanti> {
  double _granularity;
public:
  TTQFreqExtractor(TProtoQuanti* T) : TraitFileHandler<TProtoQuanti>(T,".qfreq"), _granularity(0) {}
  void set_granularity (double val) {_granularity = val;}
  virtual ~TTQFreqExtractor () {}
  virtual void FHwrite ();
  virtual void FHread (string& filename) {}
};

#endif
