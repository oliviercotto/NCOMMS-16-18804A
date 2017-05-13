/**  $Id: bitstring.h,v 1.9 2013-12-04 11:11:33 fred Exp $
 *
 *  @file bitstring.h
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
 *  @author fred
 */


#ifndef BITSRING_H
#define BITSRING_H

#include <bitset>
#include <assert.h> //for CHAR_BIT
#include <limits.h>
#include <string.h>
#include <iostream>

#ifdef _GLIBCPP_BITSET_BITS_PER_WORD
#define BITS_PER_WORD _GLIBCPP_BITSET_BITS_PER_WORD
#else
#define BITS_PER_WORD (CHAR_BIT*sizeof(unsigned long))
#endif

#ifndef _GLIBCPP_BITSET_WORDS
#define BITSET_WORDS(__n) \
((__n) < 1 ? 1 : ((__n) + BITS_PER_WORD - 1)/BITS_PER_WORD)
#else
#define BITSET_WORDS(__n) _GLIBCPP_BITSET_WORDS(__n)
#endif

#define MASK (BITS_PER_WORD == 64 ? 0xFFFFFFFFFFFFFFFF : 0xFFFFFFFF)

/**Non-template and faster implementation of std::bitset.*/
class bitstring {
  
public:
  
  typedef unsigned long _ul;
  
  class reference
  {
    friend class bitstring;
    
    _ul *_word;
    size_t _bitpos;
    
    // left undefined
    reference();
    
  public:
    reference(bitstring& bs, size_t pos)
    {
      _word = bs.getword_atPos( pos );
      _bitpos = pos % BITS_PER_WORD;
    }
    
    ~reference()
    { }
    
    // For b[i] = x;
    reference& operator=(bool x)
    {
      if (x)
        *_word |= ( 1UL << ( _bitpos % BITS_PER_WORD ) );
      else
        *_word &= ~( 1UL << ( _bitpos % BITS_PER_WORD ) );
      return *this;
    }
    
    // For b[i] = b[j];
    reference& operator=(const reference& j)
    {
      if ((*(j._word) & ( 1UL << ( j._bitpos % BITS_PER_WORD ) ) ))
        *_word |= ( 1UL << ( _bitpos % BITS_PER_WORD ) );
      else
        *_word &= ~( 1UL << ( _bitpos % BITS_PER_WORD ) );
      return *this;
    }
    
    // Flips the bit
    bool operator~() const
    { return (*(_word) & ( 1UL << ( _bitpos % BITS_PER_WORD ) ) ) == 0; }
    
    // For __x = b[i];
    operator bool() const
    { return (*(_word) & ( 1UL << ( _bitpos % BITS_PER_WORD ) ) ) != 0; }
    
    // For b[i].flip();
    reference& flip()
    {
      *_word ^= ( 1UL << ( _bitpos % BITS_PER_WORD ) );
      return *this;
    }
  };
  
  friend class reference;
  
  bitstring (size_t length) 
  : _size(length), _words( BITSET_WORDS(length) ), _data(0)
  {
    _data = new _ul [_words];
    memset(_data, 0, _words * sizeof(_ul));
  }
  
  bitstring (const bitstring& b) 
  : _size(b._size), _words(b._words), _data(0)
  {
    _data = new _ul [_words];
    memcpy(_data, b._data, _words * sizeof(_ul));
  }
  
  ~bitstring () {if(_data != NULL) delete [] _data;}
  
  _ul* getword_atPos (size_t pos)
  { return &_data[ pos / BITS_PER_WORD ]; }
  
  _ul* getword_atIdx (size_t index)
  { return &_data[ index ]; }
  
  size_t size     ( ) {return _size;}
  
  size_t nb_words ( ) {return _words;}
  
  reference operator[] (size_t pos)
  { return reference(*this, pos); }
  
  bool operator[] (size_t n) const
  { return static_cast<bool>( _data[ n / BITS_PER_WORD ] & ( 1UL << ( n % BITS_PER_WORD ) ) ); }
  
  bitstring& operator= (const bitstring& b) 
  {
    _size = b._size;
    _words = b._words;
    if(_data != NULL) delete [] _data;
    _data = new _ul [_words];
    memcpy(_data, b._data, _words * sizeof(_ul));
    return *this;
  }
  
  bitstring& operator&= (const bitstring& x)
  {
    for (size_t i = 0; i < _words; i++)
      _data[i] &= x._data[i];
    return *this;
  }
  
  bitstring& operator|= (const bitstring& x)
  {
    for (size_t i = 0; i < _words; i++)
      _data[i] |= x._data[i];
    return *this;
  }
  
  bitstring& operator^= (const bitstring& x)
  {
    for (size_t i = 0; i < _words; i++)
      _data[i] ^= x._data[i];
    return *this;
  }
  
  bitstring operator~ (void)
  {
    for (size_t i = 0; i < _words; i++)
      _data[i] = ~(_data[i]);
    return *this;
  }
  
  bitstring operator& (const bitstring& x)
  {
    bitstring result(*this);
    result &= x;
    return result;
  }
  
  bitstring operator| (const bitstring& x)
  {
    bitstring result(*this);
    result |= x;
    return result;
  }
  
  bitstring operator^ (const bitstring& x)
  {
    bitstring result(*this);
    result ^= x;
    return result;
  }
  
  /**Counts number of one's in string using a bit table count.*/
  unsigned int local_popcountl (_ul wd)
  {
    unsigned char* c = (unsigned char*)&wd;
    unsigned short cnt = 0;
    
    for(unsigned int i = 0; i < sizeof(_ul); i++)
      cnt += _bit_count[ (unsigned int)c[i] ];
    
    return (unsigned int) cnt;
  }
  
  /**Count number of set bits.*/
  unsigned int count ( ) 
  {
    unsigned int cnt = 0;
    
    for(size_t i = 0; i < _words; i++)
      cnt += local_popcountl(_data[i]);
    
    return cnt;
  }
  
  /**Set a bit to 1.*/
  void set (size_t n) {_data[n / BITS_PER_WORD] |= (1UL << (n % BITS_PER_WORD));}
  
  /**Flip the bit at n.*/
  void flip (size_t n) {_data[n / BITS_PER_WORD] ^= (1UL << (n % BITS_PER_WORD));}
  
  /**Copy bits from an array of unsigned long words.*/
  void set_data (_ul* srce, size_t nbwrd)
  { 
    //    if(nbwrd != _words) { 
    //      std::cerr<<"bitstring::set_data: different sizes in memcpy!!\n";
    //      exit(1);
    //      }
    assert(nbwrd == _words);
    memcpy(_data, srce, nbwrd * sizeof(_ul));
  }
  /**Set all bits to 0.*/
  void reset ( ) { for(unsigned int i = 0; i < _words; i++) _data[i] = 0UL; }
  
  /**Unchecked copy, assumes we have sames sizes.*/
  void copy (const bitstring& b) 
  { memcpy(_data, b._data, _words * sizeof(_ul)); }
  
  /**Copy one word.*/
  void copy (const bitstring& b, size_t word_pos)
  { _data[ word_pos ] = b._data[ word_pos ]; }
  
  void copy (const bitstring& b, size_t from, size_t to)
  {
    assert(to <= _size);
    
    size_t start_w, end_w, start_l, end_l;
    _ul mask, tmpl;
    
    start_w = from / BITS_PER_WORD;
    end_w = to / BITS_PER_WORD;
    
    start_l = from % BITS_PER_WORD;
    end_l = to % BITS_PER_WORD;
    
    if(start_w != end_w) {
      //copy wihtin first word:
      mask =  (MASK << start_l); 
      
      _data[ start_w ]   &= ~mask;
      tmpl = b._data[ start_w ] & mask;
      
      _data[ start_w ] |= tmpl;
      
      //copy words in-between:
      size_t k = start_w + 1;
      
      memcpy(&b._data[k], &_data[k], (end_w - k)*sizeof(_ul));
      
//      for(size_t k = start_w + 1; k < end_w; ++k) {
//        
//        _data[ k ] = b._data[ k ];
//        
//      }
      //copy within last word:
      mask =  (MASK >> (BITS_PER_WORD - end_l) ); 
      
      _data[ end_w ]   &= ~mask;;
      tmpl = b._data[ end_w ] & mask;
      
      _data[ end_w ] |= tmpl;
      
    } else {
      //bits to copy are within a word:
      mask =  (MASK << start_l) & (MASK >> (BITS_PER_WORD - end_l) );
      
      _data[ start_w ] &= ~mask;
      tmpl = b._data[ start_w ] & mask;
      
      _data[ start_w ] |= tmpl;
      
    }
    
  }
  
private:
  
  size_t _size;
  size_t _words;
  
  _ul *_data;
  
  
  static unsigned char _bit_count[256];
  
};

#endif

