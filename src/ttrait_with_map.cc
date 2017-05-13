 /**  $Id: ttrait_with_map.cc,v 1.11.2.1 2014-04-29 17:53:31 fred Exp $
 *
 *  @file ttrait_with_map.cc
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
 *  created on @date 15.12.2011
 *  @author fred
 */

#include <math.h>
#include "ttrait_with_map.h"
#include "Uniform.h"
#include "tstring.h"
#include <algorithm>
#include <sstream>

GeneticMap TTProtoWithMap::_map;

// ----------------------------------------------------------------------------------

//                                 T R A I T   W I T H   M A P

// ----------------------------------------------------------------------------------
// copy cstor
// ----------------------------------------------------------------------------------
TTProtoWithMap::TTProtoWithMap(const TTProtoWithMap& TP):
_totRecombEventsMean(TP._totRecombEventsMean),_recombRate(TP._recombRate), 
_mapResolution(TP._mapResolution), _mapIndex(TP._mapIndex),
_numChromosome(TP._numChromosome), _numLoci(TP._numLoci)
{
  reset_recombination_pointers();
  
  _recombRatePerChrmsm = new double [_numChromosome];
  for (unsigned int i = 0; i < _numChromosome; ++i)
    _recombRatePerChrmsm[i] = TP._recombRatePerChrmsm[i];
  
  _numLociPerChrmsm = new unsigned int [_numChromosome];
  for (unsigned int i = 0; i < _numChromosome; ++i)
    _numLociPerChrmsm[i] = TP._numLociPerChrmsm[i];
  
  _chrsmLength = new unsigned int [_numChromosome];
  for (unsigned int i = 0; i < _numChromosome; ++i)
    _chrsmLength[i] = TP._chrsmLength[i];
  
  _lociMapPositions = new unsigned int [_numLoci];
  for (unsigned int i = 0; i < _numLoci; ++i)
    _lociMapPositions[i] = TP._lociMapPositions[i];

}

TTProtoWithMap::~TTProtoWithMap ()
{  
  reset_recombination_pointers(); 
}
// ----------------------------------------------------------------------------------
// addGeneticMapParameters
// ----------------------------------------------------------------------------------
void TTProtoWithMap::addGeneticMapParameters (string prefix)
{
  add_parameter(prefix + "_recombination_rate",     DBL, false, true,  0, 0.5, 0);
  add_parameter(prefix + "_genetic_map_resolution", DBL, false, true,  0, 1, 0);
  add_parameter(prefix + "_chromosome_num_locus",   MAT, false, false, 0, 0, 0);
  add_parameter(prefix + "_genetic_map",            MAT, false, false, 0, 0, 0);
  add_parameter(prefix + "_random_genetic_map",     MAT, false, false, 0, 0, 0);  
}
// ----------------------------------------------------------------------------------
// setGeneticMapParameters
// ----------------------------------------------------------------------------------
bool TTProtoWithMap::setGeneticMapParameters (string prefix)
{
  TMatrix tmp_matx;
  bool map_set = false;
  
  _paramPrefix = prefix;
  
  string param_name;
  
  _numLoci = (unsigned int)get_parameter_value(prefix + "_loci");
  
//map resolution:
  param_name = prefix + "_genetic_map_resolution";
  _mapResolution = (get_parameter(param_name)->isSet() ? get_parameter_value(param_name) : 1);
  
  //set map's resolution, depends on current resolution
  //***!! have to update current map if resolution changes !!*** <-----
  //map_resolution = _map.setResolution(map_resolution);
  
  
  reset_recombination_pointers();

// ---------------------------------------------------------------------------
//recombination rates --> "fixed" map  
  param_name = prefix + "_recombination_rate";
  if( get_parameter(param_name)->isSet() ) {
    
    if(!get_parameter(param_name)->isMatrix()) {
      
      _numChromosome = 1;
      
      _recombRate = get_parameter_value(param_name);
      
      _numLociPerChrmsm = new unsigned int[1];
      _numLociPerChrmsm[0] = _numLoci;
      
      _recombRatePerChrmsm = new double[1];
      _recombRatePerChrmsm[0] = _recombRate;
      
    } else { 
      //is a matrix parameter
      
      get_parameter(param_name)->getMatrix(&tmp_matx);
      
      if(tmp_matx.getNbRows() != 1) {
        error("\"%s_recombination_rate\" must be one-dimensional, with chromosome-specific recombination rates as elements.\n", prefix.c_str());
        return false;
      }
      
      _numChromosome = tmp_matx.getNbCols();
      
      if( !setNumLociPerChromosome(prefix) ) return false;

      //store the chromosome recombination rates:
      _recombRatePerChrmsm = new double [_numChromosome];
      
      for (unsigned int i = 0; i < _numChromosome; ++i)
        _recombRatePerChrmsm[i] = tmp_matx.get(0, i);
    }     
    
    if( !setRecombinationMapFixed() ) return false;
    
    map_set = true;
  }
  
//---------------------------------------------------------------------------------------------
//non random map:
  param_name = prefix + "_genetic_map";
  if ( get_parameter(param_name)->isSet() ) {
    
    if(map_set) {
      error("A genetic map is already specified for trait \"%s\" while parameter \"%s\" is set.\n",param_name.c_str());
      return false;
    }
      
    get_parameter(param_name)->getMatrix(&tmp_matx);
    
    _numChromosome = tmp_matx.getNbRows();
     
    //the genetic map gives the position of each locus on the map
    if(tmp_matx.length() != _numLoci) {
      error("Number of loci in \"%s_genetic_map\" different from number of loci for trait \"%s\"\n",prefix.c_str(), prefix.c_str());
      return false;
    }
    
    if( !setNumLociPerChromosome(prefix) ) return false;

    if( !setRecombinationMapNonRandom(tmp_matx) ) return false;
    
    map_set = true;
  }
  
//---------------------------------------------------------------------------------------------
//random map:
  param_name = prefix + "_random_genetic_map";
  if ( get_parameter(param_name)->isSet() ) {
    
    if(map_set) {
      error("A genetic map is already specified for trait \"%s\" while parameter \"%s\" is set.\n",param_name.c_str());
      return false;
    }
    
    get_parameter(param_name)->getMatrix(&tmp_matx);
    
    if(tmp_matx.getNbRows() != 1) {
      error("\"%s_random_genetic_map\" must have one row, with chromosome lengths as elements.\n",get_type().c_str());
      return false;
    }
    
    _numChromosome = tmp_matx.getNbCols();
    
    _chrsmLength = new unsigned int [_numChromosome];
    
    for (unsigned int i = 0; i < _numChromosome; ++i)
      _chrsmLength[i] = tmp_matx.get(0, i);
    
    if( !setNumLociPerChromosome(prefix) ) return false;

    if( !setRecombinationMapRandom() ) return false;
    
    map_set = true;
  }
  
//---------------------------------------------------------------------------------------------
//the map/recombination parameters are missing in the config file, setting default map
  if( !map_set ) {

    _recombRate = 0.5;
    _numChromosome = 1;
    
    _numLociPerChrmsm = new unsigned int[1];
    _numLociPerChrmsm[0] = _numLoci;
    
    _recombRatePerChrmsm = new double[1];
    _recombRatePerChrmsm[0] = _recombRate;

    if( !setRecombinationMapFixed() ) return false;
  }
  
//---------------------------------------------------------------------------------------------
  
  // transmit positions to the map!!
  registerGeneticMap();
  
  return true;
}
//---------------------------------------------------------------------------------------------
// num loci table
//---------------------------------------------------------------------------------------------
bool TTProtoWithMap::setNumLociPerChromosome (string prefix)
{
  assert(_numLociPerChrmsm == NULL);
  
  _numLociPerChrmsm = new unsigned int[_numChromosome];
  
  //fill the table:
  string param_name = prefix + "_chromosome_num_locus";
  if ( !get_parameter(param_name)->isSet()) {
    
    if(_numLoci % _numChromosome != 0) {
      
      error("Loci not evenly distributed on chromosome map for trait %s\n", prefix.c_str());
      return false;
      
    } else {
      
      for (unsigned int i = 0; i < _numChromosome; ++i)
        _numLociPerChrmsm[i] = _numLoci / _numChromosome;
    }
    
  } else {
    
    TMatrix numlocmat;
    
    get_parameter(param_name)->getMatrix(&numlocmat);
    
    if(numlocmat.getNbRows() != 1) {
      error("\"%s_chromosome_num_locus\" must be one-dimensional, with num. loci per chromosome as elements.\n", prefix.c_str());
      return false;
    }
    
    if(numlocmat.getNbCols() > _numChromosome) {
      error("\"%s_chromosome_num_locus\" must have same number of elements (chromosomes) as genetic map parameters.\n", prefix.c_str());
      return false;
    }
    
    unsigned int cntr=0;
    for (unsigned int i = 0; i < _numChromosome; ++i) {
      if(numlocmat.get(0, i) == 0){
        error("the genetic map doesn't accept empty chromosomes with 0 loci\n");
        return false;
      }
      _numLociPerChrmsm[i] = numlocmat.get(0, i);
      cntr += numlocmat.get(0, i);
    }
    if (cntr != _numLoci) {
      error("\"%s_chromosome_num_locus\" has more loci than \"%s_loci\"\n", prefix.c_str(), prefix.c_str());
      return false;
    }
  }
  
  return true;  
}
//---------------------------------------------------------------------------------------------
// recombination map
//---------------------------------------------------------------------------------------------
bool TTProtoWithMap::setRecombinationMapNonRandom (TMatrix& lociPositions)
{
  assert(_lociMapPositions == NULL);

  _lociMapPositions = new unsigned int [_numLoci];
  // double map_resolution = _map.getResolution();
//    cout << "TTProtoWithMap::setRecombinationMapNonRandom\n";

//  unsigned int totLength = 0;
//  
//  //set loci position in absolute size
  for(unsigned int i = 0, stride = 0; i < _numChromosome; i++) {
//    
//    _lociMapPositions[stride] = totLength;
//        
    for(unsigned int j = 0; j < _numLociPerChrmsm[i]; j++){
      _lociMapPositions[stride + j] = lociPositions.get(i, j);
//          cout << _lociMapPositions[stride + j]<<" ";
    }
//    _lociMapPositions[stride + j] = totLength + lociPositions.get(i, j);
//    
//    totLength = _lociMapPositions[stride + _numLociPerChrmsm[i] - 1];
//    
    stride += _numLociPerChrmsm[i];
//    cout << endl;
  }
  
  return true;
}
//---------------------------------------------------------------------------------------------
// recombination map with fix rates
//---------------------------------------------------------------------------------------------
bool TTProtoWithMap::setRecombinationMapFixed ()
{
  assert(_lociMapPositions == NULL);

  _lociMapPositions = new unsigned int [_numLoci];

//  cout << "TTProtoWithMap::setRecombinationMapFixed\n";

  double interval_scale;
  //double map_resolution = _map.getResolution();
  
  //find the minimum recombination rate, this will give the resolution of the map.
  for(unsigned int i = 0; i < _numChromosome; i++)
    _mapResolution = (_recombRatePerChrmsm[i] / (_mapResolution * 0.01) < 1 ? 
                      pow(10.0, floor( log10(_recombRatePerChrmsm[i]) ) + 2 ) :
                      _mapResolution);  
  
//  cout << "   map resolution: "<<_mapResolution<<endl;

  interval_scale = 0.01 * _mapResolution;

//  cout << "   interval scale = "<<interval_scale<<endl;
//  cout << "   loci positions:\n";

  //recombination rates are distance *between* loci!
  for(unsigned int i = 0, step, stride = 0; i < _numChromosome; i++) {

    step = (unsigned int)(_recombRatePerChrmsm[i] / interval_scale);
//    cout << "   chrm "<<i<<" [step="<<step<<"] ";
//    _lociMapPositions[stride] = (i == 0 ? 0 : _lociMapPositions[stride-1]);
    
    _lociMapPositions[stride] = 0;
    
    for(unsigned int j = 1; j < _numLociPerChrmsm[i]; ++j){
      _lociMapPositions[stride + j] = _lociMapPositions[stride + j - 1] + step;
//      cout << _lociMapPositions[stride + j]<<" ";
    }
//    cout << endl;
//      _lociMapPositions[stride + j] = _lociMapPositions[stride + j -1] + step;
    
    stride += _numLociPerChrmsm[i];
  }
  
  return true;
}
//---------------------------------------------------------------------------------------------
// random recombination map
//---------------------------------------------------------------------------------------------
bool TTProtoWithMap::setRecombinationMapRandom ()
{
  assert(_lociMapPositions == NULL);

  _lociMapPositions = new unsigned int [_numLoci];

//  cout << "TTProtoWithMap::setRecombinationMapRandom\n";
  
//  _mapResolution = _map.getResolution();
  
  vector<unsigned int> chr_length(_numChromosome, 0);

  for(unsigned int i = 0; i < _numChromosome; i++)
    chr_length[i] = (unsigned int)(_chrsmLength[i] / _mapResolution);
  
  //draw the positions:
  vector< unsigned int > tmp_map;
  unsigned int **chrm_map;
  
  chrm_map = new unsigned int* [_numChromosome];
  
  for (unsigned int i = 0; i < _numChromosome; ++i) {
    chrm_map[i] = new unsigned int [ _numLociPerChrmsm[i] ];
  }
  
  for(unsigned int i = 0; i < _numChromosome; i++) {
    
    tmp_map.assign(_numLociPerChrmsm[i], 0);
    
    for(unsigned int j = 0; j < _numLociPerChrmsm[i]; ++j)
      tmp_map[j] = RAND::Uniform( chr_length[i] );
    
    vector<unsigned int>::iterator vec_first = tmp_map.begin(), vec_last = tmp_map.end();
    //sort the vector
    std::sort(vec_first, vec_last);
    
    //!!there might be duplicates, we don't deal with those
    
    for(unsigned int j = 0; j < _numLociPerChrmsm[i]; ++j)
      chrm_map[i][j] = tmp_map[j];
  }
  
//  cout << "   loci positions:\n";
  for(unsigned int i = 0, stride = 0; i < _numChromosome; i++) {
    //offset = 0,
    for(unsigned int j = 0; j < _numLociPerChrmsm[i]; j++){
      _lociMapPositions[stride + j] = chrm_map[i][j];// + offset;
//      cout << _lociMapPositions[stride + j]<<" ";
    }
    //total size of chrmsm, to add to pos of next chrmsm loci
//    offset += chr_length[i];
//    cout << endl;
    
    stride += _numLociPerChrmsm[i];
  }
  
//  unsigned int totLength = _lociMapPositions[numLoci -1];
  
  for (unsigned int i = 0; i < _numChromosome; ++i) {
    delete [] chrm_map[i];
  }
  delete [] chrm_map; 
  
  // record positions into genetic_map parameter,
  // will be saved with other params in log file
  ostringstream map;
  
  map<<"{";
  
  for (unsigned int c = 0, l = 0; c < _numChromosome; c++) {
    map<<"{";
    
    for (unsigned int i = 0; i < _numLociPerChrmsm[c]-1 && l < _numLoci; i++) {
      map<<_lociMapPositions[ l++ ]<<", ";
    }
    map<<_lociMapPositions[ l++ ];
    map<<"}";
  }
  map<<"}";
  
  Param* mapParam = get_parameter( _paramPrefix + "_genetic_map" );
  mapParam->setArg(map.str());
  mapParam->setIsSet(true);
  mapParam = get_parameter( _paramPrefix + "_random_genetic_map" );
  mapParam->setIsSet(false);
  mapParam = get_parameter( _paramPrefix + "_genetic_map_resolution" );
  mapParam->setArg(tstring::dble2str(_map.getResolution())); // the local resolution is not updated after registering the trait map, other traits might have changed it
  mapParam->setIsSet(true);
  
  return true;
}
// ----------------------------------------------------------------------------------
// registerGeneticMap
// ----------------------------------------------------------------------------------
void TTProtoWithMap::registerGeneticMap()
{
  _mapIndex = _map.addTrait(this->get_type(), _numChromosome, _numLoci, _numLociPerChrmsm,
                            _mapResolution, _lociMapPositions);
  _isRegistered = true;
}
// ----------------------------------------------------------------------------------
// recombine
// ----------------------------------------------------------------------------------
void TTProtoWithMap::recombine (unsigned long indID)
{
  if(!_map.registerIndForRecombine(indID)) return;
//  cout << "\nTTProtoWithMap::recombine::registered ind "<<indID<<endl;
  _map.recombine(FEM);
  _map.recombine(MAL);
}
// ----------------------------------------------------------------------------------
// reset_recombination_pointers
// ----------------------------------------------------------------------------------
void TTProtoWithMap::reset_recombination_pointers()
{
  if(_recombRatePerChrmsm) {delete [] _recombRatePerChrmsm; _recombRatePerChrmsm = NULL;}
  if(_numLociPerChrmsm) {delete [] _numLociPerChrmsm; _numLociPerChrmsm = NULL;}
  if(_chrsmLength) {delete [] _chrsmLength; _chrsmLength = NULL;}
  if(_lociMapPositions) {delete [] _lociMapPositions; _lociMapPositions = NULL;}
}
// ----------------------------------------------------------------------------------
// reset
// ----------------------------------------------------------------------------------
void TTProtoWithMap::reset ( )
{ 
  _map.unregisterTrait( this->get_type() );
  _isRegistered = false;
}
// ----------------------------------------------------------------------------------

//                                 G E N E T I C    M A P

// ----------------------------------------------------------------------------------
// unregisterTrait
// ----------------------------------------------------------------------------------
void GeneticMap::unregisterTrait (trait_t trait)
{
  map<trait_t, unsigned int>::iterator tIter;
  
  tIter = _traits.find(trait);
    
  if ( tIter != _traits.end() )  {
    
//    cout<<"GeneticMap::unregisterTrait::"<<tIter->first<<" (idx: "<<tIter->second<<")"<<endl;
    
    _traits.erase(tIter);
    
//      cout<<"GeneticMap::unregisterTrait::done\n";
  }
  else fatal("Genetic map::unregisterTrait: trait \"%s\" is not registered\n", trait.c_str());
  
  _nTrait--;
}
// ----------------------------------------------------------------------------------
// unregisterTrait_at
// ----------------------------------------------------------------------------------
//void GeneticMap::unregisterTrait_at (unsigned int traitIdx)
//{
//  
//  if(_lociLookupTable[traitIdx])         delete [] _lociLookupTable[traitIdx];
//   
//  if(_numLociPerChrsmPerTrait[traitIdx]) delete [] _numLociPerChrsmPerTrait[traitIdx];
//  
//  if(_locPositionsPerTrait[traitIdx])    delete [] _locPositionsPerTrait[traitIdx];
//  
//  if(_recPositions[traitIdx])            delete [] _recPositions[traitIdx];
//
//  
//}
// ----------------------------------------------------------------------------------
// clear
// ----------------------------------------------------------------------------------
void GeneticMap::clear()
{
  reset_tables();
  _traits.clear();
  _nTrait = 0;
}
// ----------------------------------------------------------------------------------
// addTrait
// ----------------------------------------------------------------------------------
unsigned int GeneticMap::addTrait (trait_t trait, unsigned int nChrm, 
                                   unsigned int nLoc,
                                   unsigned int* nLocChrm, 
                                   double resolution, 
                                   unsigned int* locPositions)
{  
  map<trait_t, unsigned int>::iterator tIter;
  
  unsigned int traitIdx = 0;
  
  bool do_add = true;
    
  tIter = _traits.find(trait);
  
  if ( tIter != _traits.end() ) {
    //this shouldn't happen... unless we are loading a pop from source!!
    
    traitIdx = tIter->second;
    
    if(_numChrsmPerTrait[traitIdx] != nChrm) 
      fatal("mismatch while loading source population and resetting genetic map for trait \"%s\"\n>>>>> number of chromosomes differ (%i != %i)\n>>>>> please check match between init file and source population\n",trait.c_str(),_numChrsmPerTrait[traitIdx],nChrm);
    
    if(_numLociPerTrait[traitIdx] != nLoc) 
      fatal("mismatch while loading source population and resetting genetic map for trait \"%s\">>>>> number of loci differ (%i != %i)\n>>>>> please check match between init file and source population\n",trait.c_str(),_numLociPerTrait[traitIdx], nLoc);
    
    unsigned int *locTable = _numLociPerChrsmPerTrait[traitIdx];
    
    for (unsigned int i=0; i<nChrm; i++)
      if(locTable[i] != nLocChrm[i]) 
        fatal("mismatch while loading source population and resetting genetic map for trait \"%s\">>>>> number of loci differ (%i != %i) on chromosome %i\n>>>>> please check match between init file and source population\n",trait.c_str(),locTable[i], nLocChrm[i], i+1);
    
//    warning("ignoring locus positions saved in source population while resetting genetic map for trait \"%s\"",trait.c_str());
//    warning("loci positions are set from parameters in the init file only, please check manually for mismatch\n");
    
    do_add = false;
    
  }
  
  if (_nTrait == 0) { //this is the first trait to register its map
 
    _numChromosome = nChrm;
    _totalNumLoci = nLoc;
    _resolution = resolution;
    _totalLength = locPositions[nLoc-1];
    
    reset_tables();
    
    _perChrsmLength = new unsigned int[_numChromosome];
    _chrsmFirstLocusPosition = new unsigned int[_numChromosome];
    
  } else {
    
    if (nChrm != _numChromosome) //check chromosome numbers
      fatal("Traits in init file don't have same number of chromosomes, cannot set the genetic map!\n");

    //number of loci
    _totalNumLoci = max(_totalNumLoci, nLoc);
    
    //setting the resolution:
    if (resolution < _resolution) {
      
      rescaleMap(resolution);
      
    } else if (resolution > _resolution) {
      
      double ratio = resolution/_resolution;
      
      for(unsigned int i = 0; i < nLoc; ++i)
        locPositions[i] = (unsigned int)((double)locPositions[i]*ratio);
    }
    
  }
  
  if(do_add) {
    
    traitIdx = _nTrait;
        
    _traits[trait] = _nTrait++;
  
    
    //---- record num chrmsm
    if ( _numChrsmPerTrait.size() != _nTrait -1 ) {
      fatal("Genetic map::wrong size of table of num chrms per trait (%i != %i)", 
            _numChrsmPerTrait.size(), _nTrait -1);
    } else
      _numChrsmPerTrait.push_back(nChrm);
    
    //---- record num loci
    if ( _numLociPerTrait.size() != _nTrait -1 ) {
      fatal("Genetic map::wrong size of table of num loci per trait (%i != %i)", 
            _numLociPerTrait.size(), _nTrait -1);
    } else
      _numLociPerTrait.push_back(nLoc);
    
    //---- record table of num loci per chrmsm
    if ( _numLociPerChrsmPerTrait.size() != _nTrait -1 )
      fatal("Genetic map::wrong size of table of num loci per chrms per trait (%i != %i)", 
            _numLociPerChrsmPerTrait.size(), _nTrait-1);
    
    unsigned int* locTable = new unsigned int[nChrm];
    
    for (unsigned int i = 0; i < nChrm; i++) locTable[i] = nLocChrm[i];
    
    _numLociPerChrsmPerTrait.push_back(locTable);
        
    //---- record the locus positions 
    if ( _locPositionsPerTrait.size() != _nTrait -1 ) {
      fatal("Genetic map::wrong size of table of loc position per trait (%i != %i)", 
            _locPositionsPerTrait.size(), _nTrait-1);
    }
    
    //copy the locus positions
    unsigned int* posTable = new unsigned int[nLoc];
    
    for (unsigned int i = 0; i < nLoc; ++i) {
      posTable[i] = locPositions[i];
    }
    _locPositionsPerTrait.push_back(posTable);
    
    //---- create table to record positions of x-over
    if ( _recPositions.size() != _nTrait -1 )
      fatal("Genetic map::wrong size of rec positions table (%i != %i)", 
            _recPositions.size(), _nTrait-1);
    _recPositions.push_back(new char[nLoc+1]);//+1 is needed here, see recombine and setLookupTable below
    
    _recPositionsF.push_back(vector<unsigned int>());
    _recPositionsM.push_back(vector<unsigned int>());
    
    //---- set the map positions according to other trait's
    //---- need to correctly set starting position and length of each chromosome
    
    if (_nTrait == 1) { //first trait to register
      
      for (unsigned int c = 0, stride = 0; c < _numChromosome; ++c) {
        
        _chrsmFirstLocusPosition[c] = locPositions[stride];
        _perChrsmLength[c] = locPositions[stride + nLocChrm[c] - 1] - locPositions[stride]; //position of last locus - position of first locus on that chromosome.
        stride += nLocChrm[c];
      }
      
    } else {
      
      for (unsigned int c = 0, lastPos, stride = 0; c < _numChromosome; ++c) {
        
        lastPos = max(locPositions[stride + nLocChrm[c] - 1], 
                      _chrsmFirstLocusPosition[c] + _perChrsmLength[c]);
        
        _chrsmFirstLocusPosition[c] = min(_chrsmFirstLocusPosition[c], locPositions[stride]);
        
        _perChrsmLength[c] = lastPos - _chrsmFirstLocusPosition[c];
        
        stride += nLocChrm[c];
      }
    }
    
    //---- compute total chromosome length
    _totalLength = 0;
    for (unsigned int c = 0; c < _numChromosome; ++c) _totalLength += _perChrsmLength[c];
    
    //---- mean num recombination events, map length of 1M = 1 x-over on average
    _totRecombEventsMean = (double)_totalLength * 0.01 * _resolution;
    
     
    
    //---- create new slot for the lookup table of that trait
    //     will be allocated and filled in setLookupTable
    if ( _lociLookupTable.size() != _nTrait -1 )
      fatal("Genetic map::wrong size of loci lookup tables (%i != %i)", 
            _lociLookupTable.size(), _nTrait-1);
    
    _lociLookupTable.push_back(NULL);
    
    //---- set all lookup tables (all genetic maps):
    //     we reset all tables each time a trait is added because _totalLength changes
    for (unsigned int t = 0; t < _nTrait; ++t) setLookupTable(t);
    
  }
  
  return traitIdx;
}
// ----------------------------------------------------------------------------------
// rescaleMap
// ----------------------------------------------------------------------------------
void GeneticMap::rescaleMap (double val)
{
  assert(val < _resolution);
  
  double ratio = _resolution/val; //it is assumed that val < _resolution
  
  _resolution = val;
  
  //we have to update all positions 
  for (unsigned int i = 0; i < _nTrait; ++i) {
    for (unsigned int j = 0; j < _numLociPerTrait[i]; ++j) {
      _locPositionsPerTrait[i][j] *= ratio;
    }
  }
  
  for (unsigned int i = 0; i < _numChromosome; ++i) {
    _perChrsmLength[i] *= ratio;
    _chrsmFirstLocusPosition[i] *= ratio;
  }
  
  _totalLength *= ratio;
  _totRecombEventsMean = (double)_totalLength * 0.01 * _resolution;

  //reset all lookup tables, the total length has changed!
  for (unsigned int i = 0; i < _nTrait; ++i) setLookupTable(i);
}
// ----------------------------------------------------------------------------------
// setLookupTable
// ----------------------------------------------------------------------------------
void GeneticMap::setLookupTable(unsigned int idx)
{
//  assert(_lociLookupTable[idx] != NULL);
  
  if(_lociLookupTable[idx] != NULL)  delete [] _lociLookupTable[idx];

  _lociLookupTable[idx] = new unsigned int [_totalLength +1];
  //!!!! +1 needed here because map must start with 0 !!!!
  
  if(_lociLookupTable[idx] == NULL)
    fatal("GeneticMap::set lookup table: memory exhausted?\n");
  
  //------------------------------------------------------------------------
  int pre_offset = 0, post_offset = 0, offset;
//  cout<<"\n-- Genetic map #"<<idx<<":"<<endl;
  
  for (unsigned int loc = 0, pos = 0, stride = 0, 
       c = 0; c < _numChromosome; ++c) {

//    cout<<"Chrm "<<c+1<<" ("<<_perChrsmLength[c]<<" ["<<_resolution<<"cM]):\n"<<"{ ";
    
    stride += _numLociPerChrsmPerTrait[idx][c];

    //the map must start with position zero for the first locus with lowest position
    pre_offset = -_chrsmFirstLocusPosition[c];

    //chromosomes are contiguous 
    offset = pre_offset + post_offset;
    
    for(; loc < stride && loc < _numLociPerTrait[idx]; ++loc) 
    {
      
      while (pos <= ( _locPositionsPerTrait[idx][loc] + offset) &&
             ((int)pos - post_offset) <= (int)_perChrsmLength[c]) 
             //table has 1 more elmnt, hence <=
      
            _lociLookupTable[idx][pos++] = loc;

//            cout << loc<<":"<<pos-1<<"("<<(int)pos - post_offset<<"; tab["<<pos-1<<"]="<<_lociLookupTable[idx][pos-1]<<") ";
    
    }
    
    post_offset += _perChrsmLength[c];
    
    while (pos < _totalLength && loc >= _numLociPerTrait[idx]) 
    {
    
      _lociLookupTable[idx][pos++] = _numLociPerTrait[idx]; 
      //!!!! this means mumLocus may be returned, adjust _recPositions accordingly !!!! 
      //(i.e. must contain num loci + 1 elements)!!!!

//            cout << loc<<":"<<pos-1<<"(tab["<<pos-1<<"]="<<_lociLookupTable[idx][pos-1]<<") ";
      
    }
//    cout << "}\n";
  }
   
  
#ifdef _DEBUG_
  cout << "\nGeneticMap::setLookupTable ("<<idx<<")\n";
  
  cout << "lookup table: [";
  for(unsigned int i = 0; i < _totalLength; ++i)
    cout << _lociLookupTable[idx][i] << " ";
  cout << "]\n";
  
  pre_offset = 0; post_offset = 0;
  for(unsigned int c = 0, stride = 0; c < _numChromosome; c++) {
    
    cout<<"+++Chrm "<<c+1<<" ("<<_perChrsmLength[c]<<" ["<<_resolution<<"cM]):\n";
    cout << "Loc positions: [ ";
    
    for(unsigned int i = 0; i < _numLociPerChrsmPerTrait[idx][c]; ++i)
      cout << _locPositionsPerTrait[idx][i+stride] << " " ;
    
    cout << "]\n";
   
    pre_offset = -_chrsmFirstLocusPosition[c];
    offset = pre_offset + post_offset;

    cout << "lookup table: [ ";
    
    for(unsigned int i = _chrsmFirstLocusPosition[c]+offset; i < _chrsmFirstLocusPosition[c]+_perChrsmLength[c]+offset; ++i)
      cout << _lociLookupTable[idx][i] << " " ;
    
    cout << "]\n{";
    
    for(unsigned int i = 0; i < _numLociPerChrsmPerTrait[idx][c]-1; ++i)
      cout<<_locPositionsPerTrait[idx][i+stride]+offset<<":"<<_lociLookupTable[idx][_locPositionsPerTrait[idx][i+stride]+offset]<<", ";    
    
    cout<<_locPositionsPerTrait[idx][stride + _numLociPerChrsmPerTrait[idx][c]-1]+offset
        <<":"<<_lociLookupTable[idx][_locPositionsPerTrait[idx][stride + _numLociPerChrsmPerTrait[idx][c]-1]+offset]<<"}\n";
    
    post_offset += _perChrsmLength[c];
    stride += _numLociPerChrsmPerTrait[idx][c];
  }
  cout<<endl;
  //  cout<<"Loc pos lookup table:"<<endl<<"{";
  //  for(unsigned int i = 0; i < _superChrsmLength; ++i)
  //    cout<<_recLociPositionTable[i]<<", ";
  //  cout<<"}"<<endl;  

  cout<<"   Tot map length: "<<_totalLength<<endl;
  cout<<"   Mean num recombination events: "<<_totRecombEventsMean<<endl;
#endif
}
// ----------------------------------------------------------------------------------
// recombine
// ----------------------------------------------------------------------------------
void GeneticMap::recombine (sex_t SEX)
{
  unsigned int nbRec;
  vector< vector< unsigned int > >* recPosSex = &_recPositionsF;
  vector< bool >* chromFirstRecPos = &_chrsmFirstRecombPositionF;
  
  if (SEX == MAL) {
    recPosSex = &_recPositionsM;
    chromFirstRecPos = &_chrsmFirstRecombPositionM;
  }
    
//  cout << "\nGeneticMap::recombine: sex="<<SEX;

  for (unsigned int t = 0; t < _nTrait; t++) {
    memset( _recPositions[t], 0, (_numLociPerTrait[t]+1)*sizeof(char) );
    recPosSex->at(t).assign(0,0);
  }
  
  //create a gamete:
  nbRec = (unsigned int)RAND::Poisson( _totRecombEventsMean );
 
  if(nbRec > _totalNumLoci - 1) nbRec =  _totalNumLoci - 1; 
    
//  cout<<"\ndrawing "<<nbRec<<" recombination events:"<<endl;
  
  for(unsigned int i = 0, pos, hit; i < nbRec; i++) {
    
    pos = RAND::Uniform()*_totalLength;
//    cout<<i<<": at map pos "<<pos<<endl;
    for (unsigned int t = 0; t < _nTrait; t++) {
      
      hit = _lociLookupTable[t][pos]; //this gives the locus to the right of this position
      //!!!! hit may be == num locus ==> _recPositions must have _numLocus+1 elements!! see line 773
//      cout<<"   <-- locus "<<hit<<flush;
      
      _recPositions[t][hit] ^= 1;  //two crossing-overs at the same location resets recombination flag to 0
      
//      cout<<"("<<(int)_recPositions[t][hit]<<")"<<" on trait "<<t+1<<endl;
    }
  }
 
  if(nbRec > 0) {
    //set up trait recombination spots
    for (unsigned int t = 0; t < _nTrait; t++) {
      
//      cout << "setting recombination spots for trait "<<t+1<<endl;
      
      for(unsigned int l = 0; l < _numLociPerTrait[t];l++) {
        
        while ( !_recPositions[t][l] && l < _numLociPerTrait[t] ) l++;
        
        recPosSex->at(t).push_back(l);
        
//        cout << "  rec at locus "<<l<<endl;
      }
    }
  } 
  else {
    for (unsigned int t = 0; t < _nTrait; t++) {
      recPosSex->at(t).push_back(_numLociPerTrait[t]);
    }
  }

//  cout << "setting chromosome first rec positions:\n";
  //set the starting point for each chromosome
  chromFirstRecPos->assign(_numChromosome, 0);
  for(unsigned int c = 0; c < _numChromosome; ++c) {
    (*chromFirstRecPos)[c] = RAND::RandBool();
#ifdef _DEBUG_
//    cout << "   c"<<c<<":"<<(*chromFirstRecPos)[c]<<endl;
#endif
  }
}
// ----------------------------------------------------------------------------------
// inherit
// ----------------------------------------------------------------------------------
/*inline void GeneticMap::inherit (unsigned int trait, sex_t SEX, char* seq,
                                     char** parent, size_t locByteSize)
{
  register bool flipper;
  
  vector< unsigned int >& recTable = getRecLoci(SEX, trait);
  
  unsigned int nbRec = recTable.size();
  
  cout << "GeneticMap::inherit; nb Rec = "<<nbRec<<endl;

  if (!nbRec) return;
  
  for(unsigned int c = 0, rec = 0, stride = 0, prevLoc = 0, chrm_bloc;
      c < _numChromosome; ++c) {
    
    flipper = _chrsmFirstRecombPosition[c];
    
    chrm_bloc = stride + _numLociPerChrsmPerTrait[trait][c]; //number of loci considered so far
    
    prevLoc = stride; //stride is the first locus of a chromosome
    
    //do recombination chromosome-wise
    cout<<"chrm "<<c<<endl;
    for(; recTable[rec] < chrm_bloc && rec < nbRec; rec++) {
      
      cout << "copy loci from "<<prevLoc<<" to "<<recTable[rec]<<" on chrmsm "<<flipper<<endl;
      
      memcpy(&seq[prevLoc], &parent[flipper][prevLoc], 
             (recTable[rec] - prevLoc) * locByteSize);
      
      prevLoc = recTable[rec];
      
      flipper = !flipper;
    }
    cout << "copy end of chrmsm from "<<prevLoc<<" to "<<recTable[rec]<<" on chrmsm "<<flipper<<endl;
    //copy what's left between the last x-over point and the end of the chrmsme
    memcpy(&seq[prevLoc], &parent[flipper][prevLoc], 
           (chrm_bloc - prevLoc) * locByteSize);
    
    stride += _numLociPerChrsmPerTrait[trait][c];
  }
}
*/
// ----------------------------------------------------------------------------------
// reset_tables
// ----------------------------------------------------------------------------------
void GeneticMap::reset_tables()
{
  if(_perChrsmLength != NULL) delete [] _perChrsmLength;
  _perChrsmLength = NULL;
  
  if(_chrsmFirstLocusPosition != NULL) delete [] _chrsmFirstLocusPosition;
  _chrsmFirstLocusPosition = NULL;
  
  _numChrsmPerTrait.clear();
  _numLociPerTrait.clear();
  
  for (unsigned int i = 0; i < _lociLookupTable.size(); ++i) {
    if(_lociLookupTable[i] != NULL) delete [] _lociLookupTable[i];
    else error("GeneticMap::reset_tables::found null pointer in _lociLookupTable\n");
  }
  _lociLookupTable.clear();
  
  for (unsigned int i = 0; i < _numLociPerChrsmPerTrait.size(); ++i) {
    if(_numLociPerChrsmPerTrait[i] != NULL) delete [] _numLociPerChrsmPerTrait[i];
    else error("GeneticMap::reset_tables::found null pointer in _numLociPerChrsmPerTrait\n");
  }
  _numLociPerChrsmPerTrait.clear();
  
  for (unsigned int i = 0; i < _locPositionsPerTrait.size(); ++i) {
    if(_locPositionsPerTrait[i] != NULL) delete [] _locPositionsPerTrait[i];
    else error("GeneticMap::reset_tables::found null pointer in _locPositionsPerTrait\n");
  }
  _locPositionsPerTrait.clear();
  
  for (unsigned int i = 0; i < _recPositions.size(); ++i) {
    if(_recPositions[i] != NULL) delete [] _recPositions[i];
    else error("GeneticMap::reset_tables::found null pointer in _recPositions\n");
  }
  _recPositions.clear();
  
  //following calls correctly call vector< > d-stor
  _recPositionsF.clear();
  _recPositionsM.clear();
  
  _chrsmFirstRecombPositionF.clear();
  _chrsmFirstRecombPositionM.clear();
  
}
