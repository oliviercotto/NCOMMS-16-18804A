/** $Id: paramsparser.cc,v 1.10.2.5 2016-03-10 13:24:02 fred Exp $
*
*  @file paramsparser.cc
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

#include <iostream>
#include <sstream>
#include <streambuf>
#include <fstream>
#include <unistd.h>
#include <errno.h>
#include "paramsparser.h"
#include "output.h"

#define EOL '\n'

using namespace std;

//----------------------------------------------------------------------------------------
// getParameters
// ----------------------------------------------------------------------------------------
map< string, string >& ParamsParser::getParameters(const char* stream_name)
{
#ifdef _DEBUG_
  cout<<"ParamsParser::getParameters:"<<endl;
#endif  
  if(stream_name == NULL) {
    if(_sname == NULL)
      fatal("ParamsParser::getParameters::no stream attached to the parser!\n");
    else if( !(read(_sname)) )
      fatal("ParamsParser::getParameters::failed from file \"%s\"\n",_sname);
  } else if( !(read(stream_name)) )
    fatal("ParamsParser::getParameters::failed from file \"%s\"\n",stream_name);
  
  if(stream_name != NULL) _sname = stream_name;
  
  return _inputParams;
}
//----------------------------------------------------------------------------------------
// getParsedParameters
// ----------------------------------------------------------------------------------------
map< string, vector< string > >& ParamsParser::getParsedParameters(const char* stream_name)
{
#ifdef _DEBUG_
  cout<<"ParamsParser::getParsedParameters:"<<endl;
#endif  
  getParameters(stream_name);
  
  parse();
  
  return _parsedParams;
}
//----------------------------------------------------------------------------------------
// parse
// ----------------------------------------------------------------------------------------
void ParamsParser::parse()
{
#ifdef _DEBUG_
  cout<<"ParamsParser::parse:"<<endl;
#endif  
  _parsedParams.clear();
  
  vector< string > argvect;
  string arg;
  
  map<string, string>::iterator param = _inputParams.begin();
  
  while(param != _inputParams.end()) {
	
	argvect.clear();
    
	if((param->first == "stats") || (param->first == "stat") )
    
      argvect.push_back(param->second);  //this one cannot be a sequential parameter

    else { //other params:
	  
      getArguments(param->second, argvect);
      
	}
	
	_parsedParams[param->first] = argvect;
    
#ifdef _DEBUG_
    message(" %s (", param->first.c_str());
    unsigned int nb=argvect.size()   ;
    for(unsigned int i=0 ;i<nb; ++i){
      cout<<argvect[i]<<" ";
    }
    message("| %i args)\n", nb);
#endif
       
	param++;
  }
  
  //scanForModifiers();
}

//----------------------------------------------------------------------------------------
// ParamsParser::getArguments
// ----------------------------------------------------------------------------------------
void ParamsParser::getArguments (string& arg_str, vector< string >& arg_vect)
{
  istringstream IN;
  string arg;
  int dummy_line = 0;
  char c;
  
  IN.str(arg_str);
  
  while(IN.get(c) && IN.good() && !IN.eof()){
    
    if(isspace(c)) {
      StreamParser::removeSpaceAndComment(IN, dummy_line);
      continue;
    }
    
    if(c == '{' || c == '(' || c == '['){
      getBlockArgument(IN, c, arg);
      arg_vect.push_back( arg );
    } 
    else {
      IN.putback(c);
      IN>>arg;
      arg_vect.push_back(arg);
    }
  }
}
//----------------------------------------------------------------------------------------
// ParamsParser::getBlockArgument
// ----------------------------------------------------------------------------------------
void ParamsParser::getBlockArgument (istream& IN, char& start_c, string& arg)
{
  char e = 0;
  int dummy_line = 0;

  switch (start_c) {
    case '{': e = '}';
      break;
    case '(': e = ')';
      break;
    case '[': e = ']';
      break;
    default:
      break;
  }
      
  arg = StreamParser::readUntilCharacter(IN, dummy_line, start_c, (const char)e);
    
//    cout<<"added arg to vector: "<<out<<" next char='"<<IN.peek()<<"'\n";
}
//----------------------------------------------------------------------------------------
// StreamParser::read
// ----------------------------------------------------------------------------------------
bool StreamParser::read(const char* stream)
  {
    int linecnt = 0;    
    //output string to collect parameter arguments
    string args;
    //string to store parameter name
    string key;
    
    ParamsParser::reset_inputParams();
    
    //--------------------------------------------------------------------------------
    char c = 0, eof = '\0';
    //put the char stream into a string
    string input_stream(stream);
    //add the terminating character to the string:
    input_stream += eof;
    
    //check if LF 
    if(input_stream.find_first_of('\r') != string::npos) {
      //guess if NL is CRLF as in DOS files:
      if(input_stream.find_first_of('\n') == input_stream.find_first_of('\r') + 1)
        replaceCR(input_stream, ' '); //remove CR
      else
        replaceCR(input_stream); //replaces CR with LF
    }
    
    //initiate the input stringstream:
    istringstream IN(input_stream);

    //--------------------------------------------------------------------------------
    //read the file line by line
    while(IN.good() && !IN.eof()) {
      
      linecnt++;
      
      // remove the forgoing space of the line, returns false if EOL is reached
      if(!removeSpaceAndComment(IN, linecnt)) continue;
      
      key="";
      //read the parameter name:
      while(IN.get(c) && IN.good() && !IN.eof() && c != EOL && !isspace(c)){
        //whitespace is the split char between a parameter and its argument
        //this basically does the same as IN>>key, but we have to chek for comments:
        if(c == '#'){
          removeComment(IN, linecnt);
          break;
        }
        
        key += c;
      } //__end_while_key__
      
      if(c == EOL || !removeSpaceAndComment(IN, linecnt) ) {
        if(c != eof){
          if(key.size() != 0) {//boolean param
            args = "1";
            ParamsParser::add_inputParam(key, args);
          }
#ifdef _DEBUG_
          message("  line %i: %s (1)\n", linecnt, key.c_str());
#endif 
        }
          continue;
      }
      //remove whitespace before arguments
//      if( !removeSpaceAndComment(IN, linecnt) ) continue;
                      
#ifdef _DEBUG_
      message("  line %i: %s (", linecnt, key.c_str());
#endif
            
      //read the arguments:
      args = "";
      while( readArguments(IN, linecnt, args) );
      
      //remove any trailling space on the argument string
      int i;
      for(i=(int)args.size()-1; i>=0; --i){
        if(!isspace(args[i])) break;
      }
      args[i+1]='\0';
      
      ParamsParser::add_inputParam(key, args.c_str());
      
#ifdef _DEBUG_
      message("%s)\n", args.c_str());
#endif
      //reset the stream state, is changed by operator >>
      IN.clear(ios::goodbit);
            
    }//__END__WHILE__
    
    return true;
  }
//----------------------------------------------------------------------------------------
// StreamParser::removeSpaceAndComment
//----------------------------------------------------------------------------------------
/**Removes whitespace char on a line until a non-ws or EOL is reached.
* Returns false if EOL or EOF is reached or true otherwise.
*/
bool StreamParser::removeSpaceAndComment(istream & IN, int& l_count, bool keepLast)
{
  char c;
  
  while(IN.get(c) && IN.good() && !IN.eof() && c != EOL) {
        
    if(!isspace(c)){
      
      if(c=='#') return removeComment(IN, l_count); // recursively remove comment
      
      IN.putback(c); //this is a parameter: put the character back
      
      if(keepLast) {
        //get back one char:
        IN.unget();
        IN.get(c); //read it again and check if whitespace:
        if(isspace(c)) IN.putback(c);
      }
      
      return true;
    }
  }
  return false;
}
//----------------------------------------------------------------------------------------
// StreamParser::removeComment
//----------------------------------------------------------------------------------------
/**Recusively removes comments until the end of a line/of the file, or of a block comment is reached.
*  
* #: commented line (removed until the end of the line is reached)
* #/ ... /#: a block comment (may span several lines)
* Consecutive lines of comments are also removed, even if not part of a block comment.
* Note: this function always returns false, unless something remains on a line after a 
* block comment (i.e. if removeSpaceAndComment() returns true)
*/
bool StreamParser::removeComment (istream & IN, int& l_count)
{
  char c;
  
  //remember: we enter the loop with c = '#'
  //check if next char is the start of a block comment:
  bool isBlock = (IN.peek() == '/');
  bool prevIsComment = true;
  
  while(IN.get(c) && IN.good() && !IN.eof()){
    
    //break if EOL && next line not a comment, continue otherwise:
    if (c == EOL){
      //check if next line is also a comment, start the loop again:
      if (IN.peek() == '#') {
        
        IN.get(c);
        prevIsComment = true;
        ++l_count;
        continue;
        
      //continue if within a block comment:
      } else if( isBlock ) {
        ++l_count;
        continue;
      //next line not a comment, get out of the loop:
      } else
        return false;
    }
    
    //check if we had the block str '#/' within a commented line
    if (c == '/'){ 
      
      if(IN.peek() == '#') {
        IN.get(c);
        
        if(isBlock)
          //block termination string '/#', remove trailing space
          //note: the string '#/#' is also considered as a terminating string
          return removeSpaceAndComment(IN, l_count);
        
      } else if(prevIsComment) { 
        //we've got '#/', a block comment may start anywhere on a comment line
        isBlock = true;
      }
    }
    
    if(c == '#') prevIsComment = true;
    else prevIsComment = false;
    
  } //__END_WHILE__
  
  return false;
}
//----------------------------------------------------------------------------------------
// StreamParser::readArguments
//----------------------------------------------------------------------------------------
bool StreamParser::readArguments(istream & IN, int& l_count, string& args)
{
  char c;
  
  //read one byte at a time
  while(IN.get(c) && IN.good() && !IN.eof() && c != EOL){

    if(c == '\\')
      eatLine(IN, l_count);
    //block delimiters: the block will the be parsed in Param
    else if(c == '{')
      args += readUntilCharacter(IN, l_count, c, '}');
    else if(c == '(')
      args += readUntilCharacter(IN, l_count, c, ')');
    else if(c == '[')
      args += readUntilCharacter(IN, l_count, c, ']');
    else if(c == '\"')
      args += readUntilCharacter(IN, l_count, c, '\"');
    //comment delimiter:
    else if(c == '#') {
      if(!removeComment(IN, l_count))
        return false; // reached end of line

    //argument in external file
    } else if (c == '&') {
      
      //IN.putback(c);
      //read the file name
      string file;
      int extcnt = 1;
      while(IN.get(c) && IN.good() && !IN.eof() && !isspace(c)){
        if(c=='#'){      // check if a comment is included
          IN.putback(c);                                       // comment
          if(!removeSpaceAndComment(IN, l_count)) break;
        }
        file += c;
      }

      if(c==EOL || isspace(c)) IN.putback(c);

      //open the external file and read the parameter argument from it
      //put the argument in a dummy string, this is only as a check
//      ifstream EXT(file.c_str());
//      string dummy;

//      if(!EXT) fatal("External parameter file '%s' could not be found!\n", file.c_str());

//      while(readArguments(EXT, extcnt, dummy));

//      EXT.close();
      
      //put the file name back into the arg, will be processed later
      args = "&" + file;

    } else
      args += c;
    
  }//end while read args

  if(c == EOL || IN.eof())  return false;
  
  if(!IN.good())  //this may happen if an external file doesn't end with a \n
    fatal("problem reading input file; line %i appears corrupted, skipping line\n",l_count);
  
  return true; 
}
//----------------------------------------------------------------------------------------
// StreamParser::readUntilCharacter
//----------------------------------------------------------------------------------------
string StreamParser::readUntilCharacter(istream &IN, int& l_count, char& start_c, const char end_c)
{
  string out;
  char c;
  bool closed = false;
    
  out += start_c;
  //cout << "start block : "<<start_c<<endl;
  while (IN.get(c) && IN.good() && !IN.eof() ) {
    
    if(c == EOL) //a block can span several lines
    
      ++l_count;
    
    else if(c == end_c) {
      
      out += c;
      //remove trailing spaces but keep last one, it's the arg seperator
      if(!removeSpaceAndComment(IN, l_count, true)) {
        IN.putback(EOL); //we've reached eol, put it back
      }
      
      closed = true;
      break;
    
    } else if(c == start_c) {//nested blocks

      out += readUntilCharacter(IN, l_count, start_c, end_c);
    
    } else if(c == '\\') //line continuation within a block
      
      eatLine(IN, l_count);
      
    else if(c == '#') {
    
      if(!removeComment(IN, l_count)) l_count++;
    
    } else out += c;
  } //__end_while__
  //cout<<out<<endl;
  
  if(!closed) fatal("missing closing character '%c' on line %i.\n", end_c, l_count);
  
  //cout << "close block : "<< c << endl;
  
  return out;
}
//----------------------------------------------------------------------------------------
// StreamParser::read
// ----------------------------------------------------------------------------------------
void StreamParser::eatLine (istream & IN, int& l_count)
{
  char c;
  while(IN.get(c) && IN.good() && c != EOL && !IN.eof());
  ++l_count;
}
//----------------------------------------------------------------------------------------
// StreamParser::removeCR
// ----------------------------------------------------------------------------------------
void StreamParser::replaceCR ( string& stream, const char rpl )
{
  size_t pos = 0;
  
  while ( (pos = stream.find_first_of('\r', pos)) != string::npos) {
    stream[pos] = rpl;
  }
}

