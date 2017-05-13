/** $Id: types.h,v 1.5.2.2 2014-04-29 16:22:17 fred Exp $
*  
*  @file types.h
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
*  Created on @date 24.08.2004
*
*  @author fred
*/

#ifndef _TYPES_H
#define _TYPES_H

#include <string>

/**Sex types, males are always 0 and females 1!!**/
typedef enum {
  MAL=0, FEM=1
} sex_t;

/**Array index of the age classes in the patch sizes and containers arrays.**/
typedef enum {
  OFFSx=0, ADLTx=1, A2, A3, A4, A5, A6, A7, A8, A9, A10  //PDISPx=1, 
}age_idx;

/**Age class flags.*/
typedef unsigned int age_t;

/**No age flag.*/
#define NONE 0

/**Offspring age class flag.*/
#define OFFSPRG 1

/**Post-dispersal age class flag (pre-adults in unregulated patches).*/
//#define POSTDISP 2
#ifdef POSTDISP
#undef POSTDISP
#endif

/**Adults age class flag (breeders).
 All classes but the first. Limit is 32 age classes*/
#define ADULTS 0xFFFFFFFE

/**All ages age class flag.*/
#define ALL 0xFFFFFFFF

/**Ordering type used to record statistics in the StatRecorders.**/
typedef enum {
  GEN=2,RPL=4,PATCH=8,FLAT=16
}st_order;

/**Trait types**/
typedef std::string trait_t;
/**Max number of characters in the trait's type descriptor.*/
#define TRAIT_T_MAX 5
#define DELE "delet"
#define DISP "disp"
#define FDISP "fdisp"
#define MDISP "mdisp"
#define NTRL "ntrl"
#define DQUANT "qdisc"
#define QUANT "quant"
#define PHENO "pheno"
#define WOLB "wolb"

/**Param's types**/
typedef enum {
  BOOL,DBL,INT,STR,MAT,DIST
}param_t;

#endif

