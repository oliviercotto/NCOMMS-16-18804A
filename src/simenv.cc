/**  $Id: simenv.cc,v 1.9.2.4 2015-08-06 06:46:07 fred Exp $
 *
 *  @file simenv.cc
 *  Nemo2
 *
 *   Copyright (C) 2008-2011 Frederic Guillaume
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
 *  created on @date 31.12.2008
 * 
 *  @author fred
 */

#include "simenv.h"

#include "binarydatasaver.h"
#include "binarydataloader.h"
#include "LCEmisc.h"
#include "LCEbreed.h"
#include "LCEdisperse.h"
#include "LCEselection.h"
#include "LCEcomposite.h"
#include "LCEquanti.h"
#include "LCEpheno.h"
#include "servicenotifiers.h"
//#include "ttneutralgenes.h"
//#include "ttneutralmito.h"
//#include "ttdeletmutations_bitstring.h"
//#include "ttdispersal.h"
//#include "ttwolbachia.h"
#include "ttquanti.h"
#include "ttquantiphenotypic.h"

#include "LCEplantlifecycle.h"

//#include "ttdbmi.h"

SimRunner* SIMenv::MainSim = NULL;

// ----------------------------------------------------------------------------------------
// loadDefaultTemplates
// ----------------------------------------------------------------------------------------
void SIMenv::loadDefaultComponents( SimRunner* sim ) 
{
  //add TraitPrototype
//  sim->add_trait(new TProtoDispersal(FEM));
//  sim->add_trait(new TProtoDispersal(MAL));
//  sim->add_trait(new TProtoNeutralGenes());
//  sim->add_trait(new TProtoDeletMutations_bitstring());
//  sim->add_trait(new TProtoWolbachia());
//  sim->add_trait(new TProto_DBMI());
  sim->add_trait(new TProtoQuantiPheno());
  sim->add_trait(new TProtoQuanti());
//  sim->add_trait(new TProtoNeutralMito());
  
  //add LifeCycleEvents  
  sim->add_LCE(new BinaryDataSaver());
  sim->add_LCE(new LCE_StatServiceNotifier());
  sim->add_LCE(new LCE_FileServicesNotifier());
  sim->add_LCE(new LCE_Regulation());
  sim->add_LCE(new LCE_Aging());
  sim->add_LCE(new LCE_Aging_Multi());
  sim->add_LCE(new LCE_Patch_Extinction());
  sim->add_LCE(new LCE_Breed());
  sim->add_LCE(new LCE_Breed_Selection());
  sim->add_LCE(new LCE_Breed_Disperse());
  sim->add_LCE(new LCE_Breed_Selection_Disperse());
  sim->add_LCE(new LCE_Breed_Wolbachia());
  sim->add_LCE(new LCE_Disperse_ConstDisp());
  sim->add_LCE(new LCE_SeedDisp());
  sim->add_LCE(new LCE_Disperse_EvolDisp());
  sim->add_LCE(new LCE_Selection_base());
  sim->add_LCE(new LCE_Cross());
  sim->add_LCE(new LCE_Resize());
  sim->add_LCE(new LCE_QuantiInit());
  sim->add_LCE(new LCE_PhenoInit());
//  sim->add_LCE(new LCE_NtrlInit());
  sim->add_LCE(new LCE_clonage());

  sim->add_LCE(new LCE_Breed_Disperse_Clone_Regulate());

  sim->build_allParams();
  
}
