/*
 * LCEregulation.cc
 *
 *  Created on: 29.04.2014
 *      Author: cotto_o
 */
/*

#include <deque>
#include <algorithm>
#include "output.h"
#include "LCEregulation.h"
#include "individual.h"
#include "metapop.h"
#include "Uniform.h"
#include "simenv.h"
#include "output.h"

LCE_Carrying_Capacity::LCE_Carrying_Capacity() {
add_parameter("carrying_capacity",INT,true,false,0,0);
};


bool LCE_Carrying_Capacity::setParameters ( )
{
	_carrying_capacity = get_parameter_value("carrying_capacity");
  return LCE_Carrying_Capacity::setParameters ();
}


void LCE_Carrying_Capacity::execute()
{

double	coeff_survival;

Patch* patch;
Individual* offspring;

unsigned int pop_size_adult;

unsigned int pop_size_offspring;

	for(unsigned int i = 0; i < _popPtr->getPatchNbr(); i++) {

patch = _popPtr->getPatch(i);

pop_size_adult = patch->size(ADULTS);

pop_size_offspring = patch->size(OFFSPRG);


if(_carrying_capacity==0)fatal("Carrying capacity cannot be equal to 0 \n");

coeff_survival=(_carrying_capacity - pop_size_adult)/_carrying_capacity;

for(unsigned int size = pop_size_offspring , indexOfOff = 0;

		indexOfOff < size;

        indexOfOff++)

{
offspring = patch->get(FEM, OFFSx,indexOfOff);

if(RAND::Uniform() > coeff_survival) {

    	    	     _popPtr->recycle(offspring);

    	    	      patch->remove(FEM, OFFSx, indexOfOff);

    	    	          	    	       }
				}

	}


};

*/
