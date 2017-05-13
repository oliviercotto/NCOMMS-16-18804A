
 /** LCEregulation.h
 *
 *  Created on: 29.04.2014
 *      Author: cotto_o
 */
/*
#ifndef LCEREGULATION_H_
#define LCEREGULATION_H_
#include "lifecycleevent.h"
#include "filehandler.h"
#include "Uniform.h"


class LCE_Carrying_Capacity: public virtual LifeCycleEvent{

private:

unsigned int _carrying_capacity;

public:

  LCE_Carrying_Capacity ( ); //the constructor is empty - constructor defined in .cc

  virtual ~LCE_Carrying_Capacity   ( ) { }

  virtual bool setParameters ();
  virtual void  execute ();

  virtual LifeCycleEvent* clone ( ) {return new LCE_Breed();}

  virtual void loadFileServices ( FileServices* loader ) {}
  virtual void loadStatServices ( StatServices* loader ) {}

  virtual age_t removeAgeClass ( ) {return 0;}
  virtual age_t addAgeClass ( ) {return OFFSPRG;}
  virtual age_t requiredAgeClass () {return ADULTS;}

  ///@}
};





#endif /* LCEREGULATION_H_ */


