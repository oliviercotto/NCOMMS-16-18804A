/** $Id: paramsparser.h,v 1.5.2.1 2016-03-09 15:05:53 fred Exp $
*
*  @file paramsparser.h
*  Nemo2
* 
*   Copyright (C) 2006-2015 Frederic Guillaume
*   frederic.guillaume@ieu.uzh.ch
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
*  created on @date 09.05.2005
* 
*  @author fred
*/

#ifndef PARSER_H
#define PARSER_H

#include <iostream>
#include <vector>
#include <string>
#include <map>

using namespace std;
/**Provides interface to read input parameters from various sources and parses them. */
class ParamsParser {
  
public:
  
  ParamsParser(const char* name) {_sname = name;}
  virtual ~ParamsParser(){}
  
  void setName(const char* name) {_sname = name;}
  
  map< string, string > get_inputParams() {return _inputParams;}
  map< string, vector< string > >& getParsedParameters(const char* stream_name);
  map< string, vector< string > >& getParsedParameters( ) {return _parsedParams;}
  map< string, string >& getParameters(const char* stream_name);
  map< string, string >& getParameters( ) {return _inputParams;}
  
  /**Read/parse params & args from a file or a string or an R object. 
    Params and their args are put in the _inputParams.*/
  virtual bool read(const char* stream) = 0;
  /**Builds the _parsedParams from the _inputParams.
    This defines rules of sequential, matricial params, etc.*/
  void parse();
  static void getBlockArgument (istream& IN, char& c, string& arg);
  static void getArguments (string& arg_str, vector< string >& arg_vect);
  
protected:
  
  void reset_inputParams() {_inputParams.clear();}
  void add_inputParam (string& param, const string& arg) {_inputParams[param] = arg;}
  
private:
  /**Attached file of stream name.*/
  const char* _sname;
  /**The whole, unparsed set of input parameters.*/
  map< string, string > _inputParams;
  /**The parsed set of simulation parameters after sequential parameters have been separated.*/
  map< string, vector< string > > _parsedParams;
    
};

/**Read parameters from a text buffer. */
class StreamParser: public ParamsParser {
  
public:
  StreamParser(const char* stream) : ParamsParser(stream) {}
  virtual ~StreamParser(){}
  
  virtual bool   read                    (const char* stream);
  static  bool   removeComment           (istream& IN, int& l_count);
  static  bool   removeSpaceAndComment   (istream& IN, int& l_count, bool keepLast = false);
  virtual bool   readArguments           (istream& IN, int& l_count, string& args);
  static  string readUntilCharacter      (istream& IN, int& l_count, char& start_c, const char end_c);
  static  void   eatLine                 (istream& IN, int& l_count);
  void           replaceCR                (string& stream, const char rpl = '\n');


private:
  
};
#endif
