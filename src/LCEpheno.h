/*
 * LCEpheno.h
 *
 *  Created on: 22.04.2014
 *      Author: cotto_o
 */

#ifndef LCEPHENO_H_
#define LCEPHENO_H_

#include "lifecycleevent.h"
#include "Uniform.h"


//LCE_QuantiInit
//
/**Set patch-specifiec initial genotypes values.*/

class LCE_PhenoInit : public virtual LifeCycleEvent
{

  TMatrix _init_values;
  bool _doByTraitValue;
  unsigned int _nTraits;

public:

  LCE_PhenoInit ( );

  virtual ~LCE_PhenoInit ( ) { }

  virtual void execute ();
  void init_trait_value(sex_t SEX, age_t AGE, unsigned int size, unsigned int deme, double *values);

  virtual LifeCycleEvent* clone ( ) {return new LCE_PhenoInit();}

  virtual bool setParameters ();

  //SimComponent implementation:
  virtual void loadFileServices ( FileServices* loader ) {}
  virtual void loadStatServices ( StatServices* loader ) {}
  virtual age_t removeAgeClass ( ) {return 0;}
  virtual age_t addAgeClass ( ) {return 0;}
  virtual age_t requiredAgeClass () {return 0;}

};



#endif
