/*
 * LCEpheno.cc
 *
 *  Created on: 22.04.2014
 *      Author: cotto_o
 */

#include "LCEpheno.h"
#include "ttquantiphenotypic.h"
#include "utils.h"
#include "output.h"

/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

//                          ******** LCE_QuantiInit ********

// ----------------------------------------------------------------------------------------
LCE_PhenoInit::LCE_PhenoInit ( ) : LifeCycleEvent("pheno_init","pheno"), _doByTraitValue(0)
 {
   add_parameter("pheno_init_trait_values",MAT,0,0,0,0);
 }
// ----------------------------------------------------------------------------------------
// LCE_QuantiInit::setParameters
// ----------------------------------------------------------------------------------------
bool LCE_PhenoInit::setParameters ( )
{
  TMatrix pat_mat;
  unsigned int patchNbr = _popPtr->getPatchNbr();
  TProtoQuantiPheno* proto = dynamic_cast<TProtoQuantiPheno*> ( _popPtr->getTraitPrototype("pheno") );

  _nTraits = proto -> get_nb_traits();

  if(_paramSet->isSet("pheno_init_trait_values")) {

    if(!get_parameter("pheno_init_trait_values")->isMatrix()) {
      fatal("\"pheno_init_trait_values\" accepts only a matrix as argument.\n");
      return false;
    }
    _paramSet->getMatrix("pheno_init_trait_values", &pat_mat);
    _doByTraitValue = setSpatialMatrix("pheno_init_trait_values", "num of quanti traits", &pat_mat, &_init_values, _nTraits, patchNbr);
  }
  return (_doByTraitValue);
}
// ----------------------------------------------------------------------------------------
// LCE_QuantiInit::execute
// ----------------------------------------------------------------------------------------
void LCE_PhenoInit::execute ( )
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

    if(values) delete [] values;
}

// ----------------------------------------------------------------------------------------
// LCE_QuantiInit::init_trait_value
// ----------------------------------------------------------------------------------------
void LCE_PhenoInit::init_trait_value(sex_t SEX, age_t AGE, unsigned int size,
                                      unsigned int deme, double *values)
{
  Individual* ind;
  TTQuantiPheno* trait;
//  for(unsigned h = 0; h < _nTraits; h++)

    for (unsigned int i = 0; i < size ; i++) {
      ind = _popPtr->get(SEX, AGE, i, deme);
      trait = dynamic_cast<TTQuantiPheno*> ( ind->getTrait(_LCELinkedTraitIndex) );
      trait->set_init_value(values);
      trait->init_sequence();
      trait->set_value();
    }

}

