/** $Id: LCEquanti.cc,v 1.10.2.2 2014-04-29 18:28:31 fred Exp $
 *
 *  @file LCEquanti.cc
 *  Nemo2
 *
 *   Copyright (C) 2011 Frederic Guillaume
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
 *  created on @date 17.08.2011
 * 
 *  @author fred
 */

#include "LCEquanti.h"
#include "ttquanti.h"
//#include "ttneutralgenes.h"
#include "utils.h"
#include "output.h"

/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

//                          ******** LCE_QuantiInit ********

// ----------------------------------------------------------------------------------------
LCE_QuantiInit::LCE_QuantiInit ( ) : LifeCycleEvent("quanti_init","quant"), _doByTraitValue(0),
  _doByAlleleFreq(0)
 {
   add_parameter("quanti_init_trait_values",MAT,0,0,0,0);
   add_parameter("quanti_init_freq",MAT,0,0,0,0);
 }
// ----------------------------------------------------------------------------------------
// LCE_QuantiInit::setParameters
// ----------------------------------------------------------------------------------------
bool LCE_QuantiInit::setParameters ( )
{
  TMatrix pat_mat;
  unsigned int patchNbr = _popPtr->getPatchNbr();
  TProtoQuanti* proto = dynamic_cast<TProtoQuanti*> ( _popPtr->getTraitPrototype("quant") );
  
  _nTraits = proto -> get_nb_traits();
  _nLoci   = proto->get_nb_locus();
  
  if(_paramSet->isSet("quanti_init_trait_values")) {
  
    if(!get_parameter("quanti_init_trait_values")->isMatrix()) {
      fatal("\"quanti_init_trait_values\" accepts only a matrix as argument.\n");
      return false;
    }
    _paramSet->getMatrix("quanti_init_trait_values", &pat_mat);
  
    //TMatrix _pip_show;
   // _pip_show=&pat_mat;

    //check
    //unsigned int number_columns=_pip_show->getNbCols();
    //unsigned int number_rows=*pat_mat->getNbRows();
    //

    _doByTraitValue = setSpatialMatrix("quanti_init_trait_values", "num of quanti traits", &pat_mat, &_init_values, _nTraits, patchNbr);
  
    //unsigned int number_columns=_doByTraitValue->getNbCols();

  }
  
  if(_paramSet->isSet("quanti_init_freq")) {
    
    if (proto->get_allele_model() > 2) {
      error("allelic model of quanti trait incompatible with initialization by allele frequency.");
      error("allelic model must be di-allelic instead of continous.\n");
      return false;
    }
    _paramSet->getMatrix("quanti_init_freq", &pat_mat);
        
    _doByAlleleFreq = setSpatialMatrix("quanti_init_freq", "num of quanti trait loci", &pat_mat, &_init_freq, _nLoci, patchNbr);
    
  }
  
  if (_doByTraitValue && _doByAlleleFreq) {
    warning("LCE_QuantiInit:: trying to init quanti traits both by trait values and allele frequencies.");
    warning("LCE_QuantiInit:: forcing initialization by allele frequencies.\n");
    _doByTraitValue = false;
  }
  
  if (!_doByAlleleFreq && !_doByTraitValue) {
    error("either \"quanti_init_trait_values\" or \"quanti_init_freq\" are missing for LCE quanti_init\n");
  }
    
  return (_doByTraitValue || _doByAlleleFreq);
}
// ----------------------------------------------------------------------------------------
// LCE_QuantiInit::execute
// ----------------------------------------------------------------------------------------
void LCE_QuantiInit::execute ( )
{
  if(!(_popPtr->getCurrentGeneration() == 1)) return;

  unsigned int patchNbr = _popPtr->getPatchNbr();
  
  double *values = NULL; 
  
  if (_doByTraitValue) {
    
    values = new double[_nTraits];
    
    for (unsigned int i = 0; i < patchNbr ; i++) { 
      
      _init_values.getRowView(i, _nTraits, values);  
      
      init_trait_value(FEM, OFFSPRG, _popPtr->size(FEM, OFFSPRG, i), i, values);
      init_trait_value(MAL, OFFSPRG, _popPtr->size(MAL, OFFSPRG, i), i, values);
      init_trait_value(FEM, ADULTS, _popPtr->size(FEM, ADULTS, i), i, values);
      init_trait_value(MAL, ADULTS, _popPtr->size(MAL, ADULTS, i), i, values);
    }
  }
  
  if (_doByAlleleFreq) {
    
    if(values) delete [] values;
    
    values = new double[_nLoci];
    
    double **allele_val = dynamic_cast<TProtoQuanti*> ( _popPtr->getTraitPrototype("quant") )->get_allele_values();
    
    for (unsigned int i = 0; i < patchNbr ; i++) { 
      
      _init_freq.getRowView(i, _nLoci, values);  
      
      init_allele_freq(FEM, OFFSPRG, _popPtr->size(FEM, OFFSPRG, i), i, values, allele_val);
      init_allele_freq(MAL, OFFSPRG, _popPtr->size(MAL, OFFSPRG, i), i, values, allele_val);
      init_allele_freq(FEM, ADULTS, _popPtr->size(FEM, ADULTS, i), i, values, allele_val);
      init_allele_freq(MAL, ADULTS, _popPtr->size(MAL, ADULTS, i), i, values, allele_val);
    }
  }
  
  if(values) delete [] values;
}
// ----------------------------------------------------------------------------------------
// LCE_QuantiInit::init_trait_value
// ----------------------------------------------------------------------------------------
void LCE_QuantiInit::init_trait_value(sex_t SEX, age_t AGE, unsigned int size, 
                                      unsigned int deme, double *values)
{ 
  Individual* ind;
  TTQuanti* trait;
//  for(unsigned h = 0; h < _nTraits; h++)
    
    for (unsigned int i = 0; i < size ; i++) {
      ind = _popPtr->get(SEX, AGE, i, deme);
      trait = dynamic_cast<TTQuanti*> ( ind->getTrait(_LCELinkedTraitIndex) );
      trait->set_init_value(values);
      trait->init_sequence();
      trait->set_value();
    }
  
}
// ----------------------------------------------------------------------------------------
// LCE_QuantiInit::init_allele_freq
// ----------------------------------------------------------------------------------------
void LCE_QuantiInit::init_allele_freq(sex_t SEX, age_t AGE, unsigned int size, 
                                      unsigned int deme, double *values, double **all_val)
{ 
  if(!size) return;
  
  Individual* ind;
  TTQuanti* trait;
  unsigned int num_A, num_a;
  unsigned int *ind_A = 0, *ind_a = 0;
  
  for(unsigned a = 0; a < 2; a++) {
    for(unsigned l = 0; l < _nLoci; l++) {

      num_A = (unsigned int)(size * values[l]);
      num_a = size - num_A;
      
      ind_A = new unsigned int [num_A];
      ind_a = new unsigned int [num_a];
      
      RAND::SampleSeqWithReciprocal(0, size, 1, num_A, (int*)ind_A, num_a, (int*)ind_a); //sample individuals within deme, without replacement
      
      for (unsigned int i = 0; i < num_A ; i++) {
        ind = _popPtr->get(SEX, AGE, ind_A[i], deme);
        trait = dynamic_cast<TTQuanti*> ( ind->getTrait(_LCELinkedTraitIndex) );
        trait->set_allele(l, a, all_val[l][0]);
      }
      for (unsigned int i = 0; i < num_a ; i++) {
        ind = _popPtr->get(SEX, AGE, ind_a[i], deme);
        trait = dynamic_cast<TTQuanti*> ( ind->getTrait(_LCELinkedTraitIndex) );
        trait->set_allele(l, a, all_val[l][1]);
      }
      
      delete [] ind_A;
      delete [] ind_a;
    }
  }
  
  for (unsigned int i = 0; i < size; ++i) {
    ind = _popPtr->get(SEX, AGE, i, deme);
    trait = dynamic_cast<TTQuanti*> ( ind->getTrait(_LCELinkedTraitIndex) );
    trait->set_value();
  }
  
}
/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/
/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/
/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

//                          ******** LCE_NtrlInit ********

// ----------------------------------------------------------------------------------------
//LCE_NtrlInit::LCE_NtrlInit ( ) : LifeCycleEvent("ntrl_init","ntrl"), _nLoci(0)
//{
//  add_parameter("ntrl_init_patch_freq",MAT,1,0,0,0);
//}
//// ----------------------------------------------------------------------------------------
//// LCE_NtrlInit::setParameters
//// ----------------------------------------------------------------------------------------
//bool LCE_NtrlInit::setParameters ( )
//{
//  TMatrix pat_mat;
//  unsigned int patchNbr = _popPtr->getPatchNbr();
//  TProtoNeutralGenes* proto = dynamic_cast<TProtoNeutralGenes*> ( _popPtr->getTraitPrototype("ntrl") );
//
//  if(proto->get_allele_num() > 2) {
//    error("LCE ntrl_init only works for di-allelic neutral loci, sorry.\n\
//            Please consider using an FSTAT input file for more complex cases.\n");
//    return false;
//  }
//
//  _nLoci   = proto->get_locus_num();
//
//  if(_paramSet->isSet("ntrl_init_patch_freq")) {
//
//    _paramSet->getMatrix("ntrl_init_patch_freq", &pat_mat);
//
//    if( !setSpatialMatrix("ntrl_init_patch_freq", "num of neutral loci", &pat_mat, &_init_freq, _nLoci, patchNbr) )
//      return false;
//
//  }
//
//  return true;
//}
//// ----------------------------------------------------------------------------------------
//// LCE_NtrlInit::execute
//// ----------------------------------------------------------------------------------------
//void LCE_NtrlInit::execute ( )
//{
//  if(!(_popPtr->getCurrentGeneration() == 1)) return;
//
//  unsigned int patchNbr = _popPtr->getPatchNbr();
//
//  double *values = NULL;
//
//  if(values) delete [] values;
//
//  values = new double[_nLoci];
//
//  for (unsigned int i = 0; i < patchNbr ; i++) {
//
//    _init_freq.getRowView(i, _nLoci, values);
//
//    init_allele_freq(FEM, OFFSPRG, _popPtr->size(FEM, OFFSPRG, i), i, values);
//    init_allele_freq(MAL, OFFSPRG, _popPtr->size(MAL, OFFSPRG, i), i, values);
//    init_allele_freq(FEM, ADULTS, _popPtr->size(FEM, ADULTS, i), i, values);
//    init_allele_freq(MAL, ADULTS, _popPtr->size(MAL, ADULTS, i), i, values);
//  }
//
//  if(values) delete [] values;
//}
//
//// ----------------------------------------------------------------------------------------
//// LCE_NtrlInit::init_allele_freq
//// ----------------------------------------------------------------------------------------
//void LCE_NtrlInit::init_allele_freq(sex_t SEX, age_t AGE, unsigned int size, unsigned int deme, double *values)
//{
//  if(!size) return;
//
//  Individual* ind;
//  TTNeutralGenes* trait;
//  unsigned int num_A, num_a;
//  unsigned int *ind_A = 0, *ind_a = 0;
//
//  for(unsigned a = 0; a < 2; a++) {
//    for(unsigned l = 0; l < _nLoci; l++) {
//
//      num_A = (unsigned int)(size * values[l]);  //number of '0' allele
//      num_a = size - num_A;                      //number of '1' allele
//
//      ind_A = new unsigned int [num_A];
//      ind_a = new unsigned int [num_a];
//
//      RAND::SampleSeqWithReciprocal(0, size, 1, num_A, (int*)ind_A, num_a, (int*)ind_a); //sample individuals within deme, without replacement
//
//      for (unsigned int i = 0; i < num_A ; i++) {
//        ind = _popPtr->get(SEX, AGE, ind_A[i], deme);
//        trait = dynamic_cast<TTNeutralGenes*> ( ind->getTrait(_LCELinkedTraitIndex) );
//        trait->set_allele(l, a, 0);
//      }
//
//      for (unsigned int i = 0; i < num_a ; i++) {
//        ind = _popPtr->get(SEX, AGE, ind_a[i], deme);
//        trait = dynamic_cast<TTNeutralGenes*> ( ind->getTrait(_LCELinkedTraitIndex) );
//        trait->set_allele(l, a, 1);
//      }
//
//      delete [] ind_A;
//      delete [] ind_a;
//    }
//  }
//}

