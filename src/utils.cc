/** $Id: utils.cc,v 1.1.2.2 2014-05-05 08:55:02 fred Exp $
 *
 *  @file utils.cc
 *  Nemo2
 *
 *   Copyright (C) 2013 Frederic Guillaume
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
 *  created on @date 02.05.2013
 * 
 *  @author fred
 */

#include "utils.h"
#include "Uniform.h"
#include <string>


// ----------------------------------------------------------------------------------------
// my_mean
// ----------------------------------------------------------------------------------------
double my_mean (double *data, unsigned int size)
{
  double mean = 0;
  for (unsigned int i = 0; i < size; ++i) {
    mean += data[i];
  }
  return mean/size;
}
// ----------------------------------------------------------------------------------------
// my_mean_no_nan
// ----------------------------------------------------------------------------------------
double my_mean_no_nan (double *data, unsigned int size)
{
  double mean = 0, cnt = 0;
  for (unsigned int i = 0; i < size; ++i) {
    if( !isnan(data[i]) ) { //exclude NaN values (happens for empty patches)
      mean += data[i];
      cnt++;
    }
  }
  return mean/cnt;
}
// ----------------------------------------------------------------------------------------
// my_variance_with_fixed_mean
// ----------------------------------------------------------------------------------------
double my_variance_with_fixed_mean (double *data, unsigned int size, double mean)
{
  double var = 0;
  for (unsigned int i = 0; i < size; ++i) {
    var += pow( data[i] - mean, 2 );
  }
  return var/size;
}
// ----------------------------------------------------------------------------------------
// my_variance_with_fixed_mean_no_nan
// ----------------------------------------------------------------------------------------
double my_variance_with_fixed_mean_no_nan (double *data, unsigned int size, double mean)
{
  double var = 0, cnt = 0;
  for (unsigned int i = 0; i < size; ++i) {
    if( !isnan(data[i]) ) {
      var += pow( data[i] - mean, 2 );
      cnt++;
    }
  }
  return var/cnt;
}
// ----------------------------------------------------------------------------------------
// setSpatialMatrix
// ----------------------------------------------------------------------------------------
bool setSpatialMatrix (string param, string numColCondition, TMatrix *inMat, TMatrix *outMat,
                       unsigned int nVal, unsigned int patchNbr, bool doRandomize)
{
  unsigned int ncol = inMat->getNbCols(); //nbr of local optima = nbr of traits = nbr of col
  unsigned int npat = inMat->getNbRows(); //nbr of diff. value / traits, may be <= patchNbr
  
  //we have two possible configurations:
  // 1. ncol == 1; same value for all 'traits' in a patch but varies among patches, with possible repetition of a pattern
  // 2. ncol == no traits && npat <= no patches; trait values change in each patch following a pattern (same if npat == 1)
  if(npat > patchNbr) {
    error("The number of rows in \"%s\" is greater than the number of patches, must be at least equal to this.", param.c_str());
    return false;
  }  
  if(ncol != 1 && ncol != nVal) {
    error("The number of columns in \"%s\" not properly set.\n", param.c_str());
    error("It is expected to be equal to 1 or equal to %s (%i).\n", numColCondition.c_str(), nVal);
    return false;
  }
  
  outMat->reset(patchNbr, nVal);
  
  if(doRandomize) {
    
    for(unsigned int i = 0; i < patchNbr; ++i)
      for(unsigned int j = 0; j < nVal; j++)
        outMat->set(i, j, inMat->get( RAND::Uniform(npat), j));
    
  } else {
    
    for(unsigned int i = 0; i < patchNbr; ++i) {
      
      if(npat < patchNbr) {//repetition of a pattern
        
        if(ncol == 1)
          for(unsigned int j = 0; j < nVal; j++) 
            outMat->set(i, j, inMat->get(i % npat, 0)) ;
        else
          for(unsigned int j = 0; j < nVal; j++)
            outMat->set(i, j, inMat->get(i % npat, j));
        
      } else {//different values for each Patch
        
        if(ncol == 1)
          for(unsigned int j = 0; j < nVal; j++) 
            outMat->set(i, j, inMat->get(i, 0));
        else
          for(unsigned int j = 0; j < nVal; j++)
            outMat->set(i, j, inMat->get(i, j));
      }
    }
  }
  
  

  return true;
}


