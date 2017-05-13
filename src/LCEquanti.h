/** $Id: LCEquanti.h,v 1.5.2.2 2014-04-29 18:28:31 fred Exp $
 *
 *  @file LCEquanti.h
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

#ifndef LCEQUANTI_H
#define LCEQUANTI_H

#include "lifecycleevent.h"
#include "Uniform.h"

//LCE_QuantiInit
//
/**Set patch-specifiec initial genotypes values.*/
class LCE_QuantiInit : public virtual LifeCycleEvent
{
  
  TMatrix _init_values, _init_freq;
  bool _doByTraitValue, _doByAlleleFreq;
  unsigned int _nTraits, _nLoci;
  
public:
  
  LCE_QuantiInit ( );
  
  virtual ~LCE_QuantiInit ( ) { }
  
  virtual void execute ();
  void init_trait_value(sex_t SEX, age_t AGE, unsigned int size, unsigned int deme, double *values);
  void init_allele_freq(sex_t SEX, age_t AGE, unsigned int size, unsigned int deme, double *values, double **all_val);
  
  virtual LifeCycleEvent* clone ( ) {return new LCE_QuantiInit();}
  
  virtual bool setParameters ();
  
  //SimComponent implementation:
  virtual void loadFileServices ( FileServices* loader ) {}
  virtual void loadStatServices ( StatServices* loader ) {}
  virtual age_t removeAgeClass ( ) {return 0;}
  virtual age_t addAgeClass ( ) {return 0;}
  virtual age_t requiredAgeClass () {return 0;}
  
};



//LCE_NtrlInit
//
/**Set patch-specifiec initial genotypes values.*/
class LCE_NtrlInit : public virtual LifeCycleEvent
{
  
  TMatrix _init_freq;
  unsigned int _nLoci;
  
public:
  
  LCE_NtrlInit ( );
  
  virtual ~LCE_NtrlInit ( ) { }
  
  virtual void execute ();

  void init_allele_freq(sex_t SEX, age_t AGE, unsigned int size, unsigned int deme, double *values);
  
  virtual LifeCycleEvent* clone ( ) {return new LCE_NtrlInit();}
  
  virtual bool setParameters ();
  
  //SimComponent implementation:
  virtual void loadFileServices ( FileServices* loader ) {}
  virtual void loadStatServices ( StatServices* loader ) {}
  virtual age_t removeAgeClass ( ) {return 0;}
  virtual age_t addAgeClass ( ) {return 0;}
  virtual age_t requiredAgeClass () {return 0;}
  
};

#endif
