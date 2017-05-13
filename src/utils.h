/** $Id: utils.h,v 1.2.2.2 2014-05-05 08:55:01 fred Exp $
 *
 *  @file utils.h
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

#ifndef UTILS_H
#define UTILS_H

#include "types.h"
#include "tmatrix.h"
#include <string>






extern bool setSpatialMatrix (string param, string numColCondition, TMatrix *inMat,
                              TMatrix *outMat, unsigned int nVal, unsigned int patchNbr,
                              bool doRandomize = false);

extern double my_mean (double *data, unsigned int size);

extern double my_mean_no_nan (double *data, unsigned int size);

extern double my_variance_with_fixed_mean (double *data, unsigned int size, double mean);

extern double my_variance_with_fixed_mean_no_nan (double *data, unsigned int size, double mean);

#endif
