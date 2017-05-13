/**  $Id: datatable.h,v 1.7 2012-02-17 15:02:41 fred Exp $ 
 *
 *  @file datatable.h
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
 *  created on @date 06.20.2005
 * 
 *  @author fred
 */

#ifndef DATATABLE_H
#define DATATABLE_H

#include "output.h"
/**A class to aggregate structured data in an array.
 * The information is subdivided into groups (e.g., patch) themselves structured in classes (e.g., age or sex classes).
 */
template <class T> class DataTable {
  
private:
  /**length of the table**/
  unsigned int _length;
  /**number of groups in the table**/
  unsigned int _groups;
  /**number of classes within each group**/
  unsigned int _classes;
  /**Stores the indexes of each class of each group present in the table.
   *the aggregated sizes of each classe within each group given by: size(class i!=0) = sum of size( (class 0) to (class i-1) )
   *and size(class i=0) = 0. This is repeated for each group.
   **/
  unsigned int **_cumulClassSizes;
  /**Stores the indexes of each group present in the table.
   *the aggregated sizes of each group with size(group i=0) = 0 and size(group i!=0) = sum of size( (group 0) to (group i-1)) where 
   *the group sizes are the sum of the classe sizes of the groups.
   **/
  unsigned int *_cumulGroupSizes;
  
  unsigned int **_sizes;
  
  T* _table;
  
  void store_sizes (unsigned int nbgroups, unsigned int nbclasses, unsigned int** classSizes) {
    
    if (_sizes != NULL) {
      for(unsigned int i = 0; i < nbgroups; ++i) delete [] _sizes[i];
      delete [] _sizes;
    }
    
    _sizes = new unsigned int * [nbgroups];
    
    for(unsigned int i = 0; i < nbgroups; ++i) {
      _sizes[i] = new unsigned int [nbclasses];
      for(unsigned int j = 0; j < nbclasses; ++j)
        _sizes[i][j] = classSizes[i][j];
    }
  }
  
  bool sizeChange (unsigned int nbgroups, unsigned int nbclasses, unsigned int** classSizes) {
    bool status = 0;
    if( nbgroups != _groups || nbclasses != _classes ) return true;
    else {
      for (unsigned int i = 0; i < _groups; i++)
        for (unsigned int j = 0; j < _classes; j++)
          status |= _sizes[i][j] != classSizes[i][j];
    }
    return status;
  }
  
public:
  /**Default constructor.**/
  DataTable() : _length(0), _groups(0), _classes(0), _cumulClassSizes(0), _cumulGroupSizes(0), _sizes(0), _table(0) { }
  /**Constructor.
   * @param nbgroups the number of rows in the table, typically the number of patches
   * @param nbclasses an array of the sizes of each row, i.e. the Patch sizes
   * @param classSizes a matrix of the sizes of each class within each group
   **/
  DataTable(unsigned int nbgroups, unsigned int nbclasses, unsigned int** classSizes) :
  _length(0), _groups(0), _classes(0), _cumulClassSizes(0), _cumulGroupSizes(0), _table(0) 
  {	allocate(nbgroups, nbclasses, classSizes); }
  
  ~DataTable() {free();}
  
  /**Creates a table of size given by the sum of all classes passed by the 'classSizes' matrix.
   * @param nbgroups the number of rows in the table, typically the number of patches
   * @param nbclasses the number of classes within a group, e.g. the age-classes
   * @param classSizes a matrix (nbgroups X nbclasses) of the sizes of each class within each group
   */
  inline void allocate (unsigned int nbgroups, unsigned int nbclasses, unsigned int** classSizes) { 
    //#ifdef _DEBUG_
    //    message("DataTable::allocate:\n");
    //#endif

    free();
    
    store_sizes(nbgroups, nbclasses, classSizes);
    
    _groups = nbgroups;
    _classes = nbclasses;
    _length = 0;
    
    _cumulClassSizes = new unsigned int*[_groups];

    for(unsigned int i = 0; i < _groups; i++){
      _cumulClassSizes[i] = new unsigned int [_classes];
      _cumulClassSizes[i][0] = 0;
      for(unsigned int j = 1; j < _classes; ++j)
        _cumulClassSizes[i][j] = _cumulClassSizes[i][j-1] + classSizes[i][j-1];
    }
    
    _cumulGroupSizes = new unsigned int[_groups];
    
    _cumulGroupSizes[0] = 0; 
    for(unsigned int i = 1; i < _groups; i++){
      _cumulGroupSizes[i] = _cumulGroupSizes[i-1];
      for(unsigned int j = 0; j < _classes; ++j) {
        _cumulGroupSizes[i] += classSizes[(i-1)][j];//total size of the previous group
        _length += classSizes[i][j];//aggregate to compute total length
      }
    }
    //add the sizes of first group classes to have the total length of the table
    for(unsigned int j = 0; j < _classes; ++j)
      _length += classSizes[0][j];
    
    _table = new T[_length];
    //#ifdef _DEBUG_
//        cout<<"DataTable::allocate:_table="<<_table<<endl;
    //#endif
  }
  
  /**Updates the group and classe sizes and re-allocates the table according to its new length.
   * @param nbgroups the number of rows in the table, typically the number of patches
   * @param nbclasses the number of classes within a group, e.g. the age-class sizes
   * @param classSizes a matrix of the sizes of each class within each group
   */
  inline void update(unsigned int nbgroups, unsigned int nbclasses, unsigned int** classSizes) {
    //#ifdef _DEBUG_
    //    message("DataTable::update:\n");
    //#endif
    if( nbgroups != _groups || nbclasses != _classes )
      
      allocate(nbgroups, nbclasses, classSizes); //free() is called in allocate
    
    else if ( sizeChange(nbgroups, nbclasses, classSizes) ) {
      
      store_sizes(nbgroups, nbclasses, classSizes);
      
      for(unsigned int i = 0; i < _groups; i++){
        for(unsigned int j = 1; j < _classes; ++j)
          _cumulClassSizes[i][j] = _cumulClassSizes[i][j-1] + classSizes[i][j-1];
      }
      _length = 0;
      _cumulGroupSizes[0] = 0; 
      for(unsigned int i = 1; i < _groups; i++){
        _cumulGroupSizes[i] = _cumulGroupSizes[i-1];
        for(unsigned int j = 0; j < _classes; ++j) {
          _cumulGroupSizes[i] += classSizes[(i-1)][j];
          _length += classSizes[i][j];
        }
      }
      
      for(unsigned int j = 0; j < _classes; ++j)
        _length += classSizes[0][j];
      
      if(_table != NULL) delete [] _table;
      
      _table = new T[_length];
    }
  }
  /**Deallocates all the allocated tables.*/
  void free ( ) {
//    cout << "\ncalling DataTable::free()\n";
//    cout <<  "    pointer state: table="<<_table<<" cumulGroupSizes="<<_cumulGroupSizes
//         <<"\n    cumulClassSizes="<<_cumulClassSizes<<" sizes="<<_sizes
//        <<endl;

    if(_table != NULL) delete [] _table;
    
    if(_cumulGroupSizes != NULL) delete [] _cumulGroupSizes;
//    cout << "    delete _cumulClassSizes\n";
    if(_cumulClassSizes != NULL){
      for(unsigned int i = 0; i < _groups; ++i) {
//        cout << "     "<<i<<":"<<_cumulClassSizes[i]<<endl;
        delete [] _cumulClassSizes[i];}
      delete [] _cumulClassSizes;
    }
//    cout << "    delete _sizes\n";
    if(_sizes != NULL){
      for(unsigned int i = 0; i < _groups; ++i) delete [] _sizes[i];
      delete [] _sizes;
    }
    
    _table = NULL;
    _cumulGroupSizes = NULL;
    _cumulClassSizes = NULL;
    _sizes = NULL;
    _length = 0;
//    cout << "DataTable::free() done\n";
}
  /**Returns the length of the table (total number of elements present).*/
  inline unsigned int length ( ) {return _length;}
  /**Accessor to the table array.*/
  inline T* getTable ( ) {return _table;}
  /**Accessor to a group array.*/
  inline T* getGroup (unsigned int group) {return &_table[ _cumulGroupSizes[group] ];}
  /**Accessor to a class array whithin a group.*/
  inline T* getClassWithinGroup (unsigned int group, unsigned int Class) {
    return &_table[ _cumulGroupSizes[group] + _cumulClassSizes[group][Class] ];}
  /**Returns value stored of the element 'elmnt' of the class 'Class' in the group 'group'.*/
  inline T  get (unsigned int group, unsigned int Class, unsigned int elmnt) {
    return _table[ _cumulGroupSizes[group] + _cumulClassSizes[group][Class] + elmnt ];
  } 
  /**Sets the element 'elmnt' of the class 'Class' in the group 'group' to the value 'val'.*/
  inline void set (unsigned int group, unsigned int Class, unsigned int elmnt, T val) {
    _table[ _cumulGroupSizes[group] + _cumulClassSizes[group][Class] + elmnt ] = val;
  }
  /**Increments 'elmnt' of the class 'Class' in the group 'group' by one.*/
  inline void increment (unsigned int group, unsigned int Class, unsigned int elmnt) {
    ++_table[ _cumulGroupSizes[group] + _cumulClassSizes[group][Class] + elmnt ];
  }
  /**Adds 'val' to 'elmnt' in the class 'Class' in the group 'group'.*/
  inline void plus (unsigned int group, unsigned int Class, unsigned int elmnt, T val) {
    _table[ _cumulGroupSizes[group] + _cumulClassSizes[group][Class] + elmnt ] += val;
  }
  /**Substracts 'val' from 'elmnt' in the class 'Class' in the group 'group'.*/
  inline void minus (unsigned int group, unsigned int Class, unsigned int elmnt, T val) {
    _table[ _cumulGroupSizes[group] + _cumulClassSizes[group][Class] + elmnt ] -= val;
  }
  /**Multiplies 'elmnt' of the class 'Class' in the group 'group' by 'val'.*/
  inline void multiply (unsigned int group, unsigned int Class, unsigned int elmnt, T val) {
    _table[ _cumulGroupSizes[group] + _cumulClassSizes[group][Class] + elmnt ] *= val;
  }
  /**Subdivide 'elmnt' of the class 'Class' in the group 'group' by 'val'.*/
  inline void subdivide (unsigned int group, unsigned int Class, unsigned int elmnt, T val) {
    _table[ _cumulGroupSizes[group] + _cumulClassSizes[group][Class] + elmnt ] /= val;
  }
  /**Sets all elements of the table to value 'val'.*/
  inline void init (T val) {for(unsigned int i = 0; i < _length; ++i) _table[i] = val;} 
  
  unsigned int getNumGroups () {return _groups;}
  unsigned int getNumClasses () {return _classes;}
  unsigned int size (unsigned int i, unsigned int j) {return _sizes[i][j];}
  
  void show_up ()
  {
    message("DataTable: got %i groups of %i classes; total length is %i.\n", _groups, _classes, _length);
  }
};

#endif //DATATABLE_H
