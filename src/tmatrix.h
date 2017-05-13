/** $Id: tmatrix.h,v 1.7.2.1 2014-04-29 18:05:32 fred Exp $
 *
 *  @file tmatrix.h
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
 *  created on @date 06.04.2007
 * 
 *  @author fred
 */

#ifndef _T_MATRIX_H
#define _T_MATRIX_H

#ifdef HAS_GSL
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#endif

#include <string.h>
#include "output.h"

using namespace std;

/**A class to handle matrix in params, coerces matrix into a vector of same total size**/
class TMatrix {
private:
  
  unsigned int _rows, _cols, _length;
  
  double* _val;
  
public:
	
  TMatrix  () : _rows(0), _cols(0), _length(0), _val(0) { }
  /**copy constructor.**/
  TMatrix  (const TMatrix& mat)
  {
    _rows = mat._rows;
    _cols = mat._cols;
    _length = _rows*_cols;
    _val = new double [_length];
    memcpy(_val,mat._val,_length*sizeof(double));
  } 
  /**@brief Creates an array of doubles of size = rows*cols**/
  TMatrix ( unsigned int rows, unsigned int cols )
  {
    _length = rows*cols;
    _val = new double [_length];
    _rows = rows; _cols = cols;
  }
  
  ~TMatrix () {if(_val != NULL) delete [] _val;}
  
  /**Sets element at row i and column j to value val**/
  void set (unsigned int i, unsigned int j, double val) {
    if( i*j < _length)
      _val[i*_cols + j] = val;
    else
      error("TMatrix::set overflow!\n"); 
  }
  /**Assigns a value to all element of the matrix.*/
  void assign (double val)
  {
    for(unsigned int i = 0; i < _length; ++i) _val[i] = val;
  }
  /**Re-allocate the existing matrix with assigned rows and cols dimensions**/
  void reset (unsigned int rows, unsigned int cols) {
    if(_val != NULL) delete [] _val;
    _length = rows * cols;
    _val = new double [_length];
    memset(_val, 0, _length*sizeof(double));
    _rows = rows; _cols = cols;
  }
  /**Reset the existing matrix to the new dimensions and copies the array**/
  void reset (unsigned int rows, unsigned int cols, double* array) {
    reset(rows, cols);
    memcpy(_val, array, _length * sizeof(double));
  }  
  /**Reset members to zero state.*/
  void reset ( ) 
  {
    _rows = 0; 
    _cols = 0; 
    _length = 0;
    if(_val != NULL) delete [] _val; 
    _val = NULL;
  }
  
  /**@brief Accessor to element at row i and column j**/
  double get (unsigned int i, unsigned int j) {
    if( !((i+1)*(j+1) > _length) )
      return _val[i*_cols + j];
    else
      fatal("TMatrix::get overflow!\n"); 
    return 0;
  }
  /**Accessor to the whole array.*/
  double* get () const {return _val;}
  double* getValArray () const {return _val;}
  /**Accessor to the matrix dimensions.
   *@param dims an array of at least 2 elements to store the row [0] and column [1] numbers. May be NULL.
   *@return the total size of the matrix
   **/
  unsigned int get_dims (unsigned int* dims) {
    if(dims != NULL) { dims[0] = _rows; dims[1] = _cols; }
    return _length;
  }
  /**Gives the number of rows.*/
  unsigned int getNbRows ( ) {return _rows;}
  unsigned int nrows () {return _rows;}
  /**Gives the number of columns.*/
  unsigned int getNbCols ( ) {return _cols;}
  unsigned int ncols () {return _cols;}
  /**Returns the number of elements in the matrix.*/
  unsigned int length    ( ) {return _length;}
  /**Gives access to a column of the matrix.
   @param col index of col to view
   @param n size of the storing array passed, must be equal to no. of rows
   @param array array where the column values will be stored
   */
  void getColumnView (unsigned int col, unsigned int n, double* array)
  {
    if(col > _cols-1) {
      error("TMatrix::getColumnView: not that many columns in matrix\n");
      return;
    }
    if(n != _rows) {
      error("TMatrix::getColumnView: array size not equal to number of rows in matrix\n");
      return;
    }
    for(unsigned int i = 0; i < _rows; ++i)
      array[i] = _val[i*_cols + col];
  }
  /**Gives access to a row of the matrix.
   @param row index of row to view
   @param n size of the storing array passed, must be equal to no. of columns
   @param array array where the row values will be stored
   */
  void getRowView (unsigned int row, unsigned int n, double* array)
  {
    if(row > _rows-1) {
      error("TMatrix::getRowView: not that many rows in matrix\n");
      return;
    }
    if(n != _cols) {
      error("TMatrix::getRowView: array size not equal to number of columns in matrix\n");
      return;
    }
    for(unsigned int i = 0, stride = row*_cols; i < _cols; ++i)
      array[i] = _val[stride + i];
  }
  /**Adds a value to an element of the matrix.*/
  void plus (unsigned int i, unsigned int j, double value)
  {
    if( i*j < _length)
      _val[i*_cols + j] += value;
    else
      error("TMatrix::set overflow!\n"); 
  }

  // ***add value to all elements in a matrix

  void matrix_increment (double value)

  {for(unsigned int i = 0; i < _rows; ++i){
      for(unsigned int j = 0; j < _cols; ++j){
    	  plus(i,j,value);}}}



  /**Substracts a value from an element of the matrix.*/
  void minus (unsigned int i, unsigned int j, double value)
  {
    if( i*j < _length)
      _val[i*_cols + j] -= value;
    else
      error("TMatrix::set overflow!\n"); 
  }
  /**Multiply an element of the matrix by a value.*/
  void multi (unsigned int i, unsigned int j, double value)
  {
    if( i*j < _length)
      _val[i*_cols + j] *= value;
    else
      error("TMatrix::set overflow!\n"); 
  }
  /**Divide an element of the matrix by a value.*/
  void divide (unsigned int i, unsigned int j, double value)
  {
    if( i*j < _length)
      _val[i*_cols + j] /= value;
    else
      error("TMatrix::set overflow!\n"); 
  }
  
  void transpose ()
  {
    TMatrix tmp(_cols, _rows);
    
    for(unsigned int i = 0; i < _rows; i++)
      for(unsigned int j = 0; j < _cols; j++)
        tmp.set(j, i, get(i, j));
    
    reset(_cols, _rows, tmp.get());
  }
  
  void show_up()
  {
    message("TMatrix dimensions: rows = %i, columns = %i\n",_rows,_cols);
    for(unsigned int i = 0; i < _rows; i++) {
      for(unsigned int j = 0; j < _cols; j++)
        message("%.3f ",_val[i*_cols + j]);
      message("\n");
    }
  }
  
#ifdef HAS_GSL
  /**Copies the current matrix into a GSL matrix.
   @param M a pointer to an already allocated GSL matrix
   */
  void get_gsl_matrix(gsl_matrix* M)
  {
    if(M->size1 != _rows)
      fatal("TMatrix::get_gsl_matrix row size of input matrix doesn't match!\n");
    if(M->size2 != _cols)
      fatal("TMatrix::get_gsl_matrix col size of input matrix doesn't match!\n");
    for(unsigned int i = 0; i < _rows; ++i)
      for(unsigned int j = 0; j < _cols; ++j)
        gsl_matrix_set(M,i,j,_val[i*_cols + j]);
  }
  /**Copies a GSL matrix into the current matrix.
   @param M the original GSL matrix to copy from
   */
  void set_from_gsl_matrix (gsl_matrix* M)
  {
    if(_val != 0) delete [] _val;
    _rows = M->size1;
    _cols = M->size2;
    _length = _rows*_cols;
    _val = new double [_length];
    memcpy(_val, M->data, _length*sizeof(double));
  }
  /**Inverses the matrix using LU decomposition.*/
  void inverse ()
  {
    if(_rows != _cols) {
      error("TMatrix::inverse matrix must be square to invert!\n");
      return;
    }
    
    gsl_matrix *S = gsl_matrix_alloc(_rows,_cols);
    get_gsl_matrix(S);
    gsl_matrix *LU = gsl_matrix_alloc(_rows,_cols);
    gsl_permutation *P = gsl_permutation_alloc(_rows);
    int snum;
    
    gsl_matrix_memcpy(LU,S);
    gsl_linalg_LU_decomp(LU,P,&snum);
    gsl_linalg_LU_invert(LU,P,S);
    
    set_from_gsl_matrix(S);
    
    gsl_matrix_free(S);
    gsl_matrix_free(LU);
    gsl_permutation_free(P);
  }
#endif
};

#endif

