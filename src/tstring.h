/** $Id: tstring.h,v 1.10 2013-12-04 11:30:13 fred Exp $
 *
 *  @file tstring.h
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
 *  created on @date 29.12.2008
 * 
 *  @author fred
 */
#ifndef TSTRING_H
#define TSTRING_H

#include <sstream>
#include <string>
#include <vector>
#include "output.h"

using namespace std;

/**A global class to handle string conversions and operations. Allows to parse the parameters
   argument strings into various tokens and to remove encolsing characters from a string.*/
class tstring {

private:
  tstring(){}
  
public:
  // ----------------------------------------------------------------------------------------
  // str2uint
  // ----------------------------------------------------------------------------------------
  /**Converts a string into an unsigned integer.*/
  static unsigned int str2uint (const string& str)
  {
    istringstream IN(str);
    unsigned int i;
    IN >> i;
    return i;
  }
  // ----------------------------------------------------------------------------------------
  // str2int
  // ----------------------------------------------------------------------------------------  
  /**Converts a string into an integer.*/
  static int str2int (const string& str)
  {
    istringstream IN(str);
    int i;
    IN >> i;
    return i;
  }
  // ----------------------------------------------------------------------------------------
  // str2dble
  // ----------------------------------------------------------------------------------------  
  /**Converts a string into a double.*/
  static double str2dble (const string& str) 
  {
    istringstream IN(str);
    double d;
    IN >> d;
    return d;
  }
  // ----------------------------------------------------------------------------------------
  // int2str
  // ----------------------------------------------------------------------------------------  
  /**Writes an integer value into a string.*/
  static string int2str (const int i)
  {
    ostringstream OUT;
    OUT << i;
    return OUT.str();
  }
  // ----------------------------------------------------------------------------------------
  // ulong2str
  // ----------------------------------------------------------------------------------------  
  /**Writes an integer value into a string.*/
  static string ulong2str (const unsigned long i)
  {
    ostringstream OUT;
    OUT << i;
    return OUT.str();
  }  
  // ----------------------------------------------------------------------------------------
  // dble2str
  // ----------------------------------------------------------------------------------------  
  /**Writes a floating-point value into a string.*/
  static string dble2str (const double d) 
  {
    ostringstream OUT;
    OUT << d;
    return OUT.str();
  }  
  // ----------------------------------------------------------------------------------------
  // split
  // ----------------------------------------------------------------------------------------  
  /**splits a string into substrings (tokens) delimited by a single character.
  @param str the input string
  @param delim the token delimiter
  @param splitOnce if true, return the first token only*/
  static vector<string> split(const string& str, const char delim, bool splitOnce = false)
  {
    vector<string> tokens;
    string out(str);
    size_t pos = 0;
    
    while(out.find(delim) != string::npos) {
      pos = out.find_first_of(delim);
      tokens.push_back( removeEnclosingChar(out.substr(0, pos), ' ', ' ', true)  );
      out = out.substr(pos+1);
      if(splitOnce) break;
    }
    
    tokens.push_back( removeEnclosingChar(out, ' ', ' ', true) );
    
    return tokens;
  }
  // ----------------------------------------------------------------------------------------
  // splitExcludeEnclosedDelimiters
  // ----------------------------------------------------------------------------------------  
  /**Splits a string into substrings (tokens) delimited by a single character. The tokens
     may contain the delimiter character if it is enclosed.
     For instance: `(\@g0 {{0,1,0,10}}, \@g10 {{1,0,2,2}})' --> [0] = `\@g0 {{0,1,0,10}}'
     and [1] = `\@g10 {{1,0,2,2}}'
     @param str the input string
     @param delim the token delimiter, is a comma by default
     @param encloser the character that starts an enclosed string */
  static vector<string> splitExcludeEnclosedDelimiters (const string& str, const char delim = ',',
                                                        const string& encloser = "([{\"")
  {
    if(str.find_first_of(encloser) == string::npos) return split(str, delim);
    
    string::size_type open, delim_pos, next_index = 0;
    string pair, block, tail(str);
    vector<string> tokens;
    
    while (tail.find(delim) != string::npos) {

      delim_pos = tail.find_first_of(delim);
      
      //check if delim is within any enclosed block
      open = tail.find_first_of(encloser);
      
      //read and remove blocks
      while (open < delim_pos) {   
        
        block = getBlock(tail.substr(open, string::npos), tail.at(open));
        
        next_index = open + block.size();
        
        pair += tail.substr(0, next_index);
                
        tail = tail.substr(next_index);
                
        delim_pos = tail.find_first_of(delim);
        
        open = tail.find_first_of(encloser);        
      }
      
      pair += tail.substr(0, delim_pos);
      
      pair = removeEnclosingChar(pair, ' ', ' ', true);
      
      tokens.push_back(pair);
      pair.clear();
      
      if(delim_pos != string::npos) tail = tail.substr(delim_pos + 1);
      else break;
    }

    if(tail.size() != 0) tokens.push_back( removeEnclosingChar(tail, ' ', ' ', true) );

    return tokens;
  }
  
  
  //----------------------------------------------------------------------------------------
  // getBlock
  //----------------------------------------------------------------------------------------
  /**Reads a substring delimited by two enclosing character from a string. May or may not retain
     the enclosing characters.*/
  static string getBlock(const string& str, const char start_c, const char end_c = '\0', 
                        bool includeEnclosing = true)
  {
    string out;
    char c, end = 0;
    bool closed = false;
    
    istringstream IN( removeEnclosingChar(str, ' ', ' ', true) );
    
    IN.get(c);
    
    if(end_c == 0) {
      switch (c) {
        case '(':
          end = ')';
          break;
        case '{':
          end = '}';
          break;
        case '[':
          end = ']';
          break;
        case '"':
          end = '"';
          break;
        default:
          error("tstring::getBlock: unknown start of block \'%c\'\n", c);
          break;
      }
    } else end = end_c;
    
    out += c;

    while (IN.get(c) && IN.good() && !IN.eof() ) {
      
      if(c == end) {
        
        out += c;
        closed = true;
        break;
        
      } else if(c == start_c) {//nested blocks
        
        out += __get_block(IN, start_c, end);
        
      } else out += c;
      
    } //__end_while__

    
    if(!closed) fatal("missing closing character '%c' in %s.\n", end, str.c_str());
    
    return out;
  } 
  //----------------------------------------------------------------------------------------
  // getBlock
  //----------------------------------------------------------------------------------------
  /**Internal function used by getBlock to recursively read nested blocks. */
  static string __get_block(istringstream& IN, const char start_c, const char end_c)
  {
    string out;
    char c;
    bool closed = false;
        
    out += start_c;
//    cout << "start block : "<<start_c<<endl;
    while (IN.get(c) && IN.good() && !IN.eof() ) {
      
      if(c == end_c) {
        
        out += c;
        closed = true;
        break;
        
      } else if(c == start_c) {//nested blocks
        
        out += __get_block(IN, start_c, end_c);
        
      } else out += c;
      
    } //__end_while__
//    cout<<out<<endl;
    
    if(!closed) fatal("missing closing character '%c' in %s.\n", end_c, IN.str().c_str());
    
    //    cout << "close block : "<< c << endl;
    
    return out;
  }    
  // ----------------------------------------------------------------------------------------
  // removeChar
  // ---------------------------------------------------------------------------------------- 
  /**Removes a given character from a string.*/
  static string removeChar(const string& str, const char c)
  {
    string s(str);
    
    for (unsigned int i = 0; i < s.size(); i++) {
      if (s[i] == c) {
        s.erase(s.begin() + i);
        i--;
      }
    }
    
    return s;
  }
  // ----------------------------------------------------------------------------------------
  // removeFirstCharOf
  // ----------------------------------------------------------------------------------------  
  /**Removes the first of a character found in a string.*/
  static string removeFirstCharOf(const string& str, const char c)
  {
    string out(str);
    size_t pos;
    
    if( (pos = out.find_first_of(c)) != string::npos )
        out.erase(out.begin() + pos);
    
    return out;
  }
  // ----------------------------------------------------------------------------------------
  // removeLastCharOf
  // ----------------------------------------------------------------------------------------  
  /**Removes the last of a character found in a string.*/
  static string removeLastCharOf(const string& str, const char c)
  {
    string out(str);
    size_t pos;
    
    if( (pos = out.find_last_of(c)) != string::npos )
      out.erase(out.begin() + pos);
    
    return out;
  }
  // ----------------------------------------------------------------------------------------
  // removeEnclosingChar
  // ----------------------------------------------------------------------------------------    
  /**Removes characters enclosing a string. Can be used to remove parenthesis/brackets or trailing
     whitespace characters.
     @param str the string
     @param o the opening (front) character
     @param c the closing (back) character
     @param allowMissing if trus, allows to ignore a missing enclosing character */
  static string removeEnclosingChar (const string& str, const char o, const char c, bool allowMissing = false)
  {
    string s(str);
    size_t first, last;
    
    first = s.find_first_of(o, 0);
    if(first != 0) {
      if(!allowMissing) {
        error("tstring::removeEnclosingChar:: string \"%s\" not starting with \"%c\"\n.", s.c_str(), o);
        return s;
      }
    } else 
      s.erase(s.begin());
    
    last = s.find_last_of(c);
    if(last != s.length()-1 || last == string::npos) {
      if(!allowMissing) {
        error("tstring::removeEnclosingChar:: string \"%s\" not ending with \"%c\"\n.", s.c_str(), c);
        return s;
      }
    } else
      s.erase(s.end()-1);
    
    return s;
  }
  // ----------------------------------------------------------------------------------------
  // isanumber
  // ----------------------------------------------------------------------------------------
  /**Check whether the string is a number.*/
  static bool isanumber(const string& str)
  {
    unsigned int i = 0;
    
    while(i < str.size()) {
      if(!isdigit(str[i]))
        if(str[i] != '.' && str[i] != 'e' && str[i] != '-')
          return false;
      i++;
    }
    
    return true;
  }
  
};

#endif
