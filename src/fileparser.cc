/** $Id: fileparser.cc,v 1.6 2014-01-13 14:00:48 fred Exp $
*
*  @file fileparser.cc
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
*  Created on @date 05.09.2005
*  @author fred 
*/

#include <fstream>
#include <sstream>
#include <errno.h>
#include <string.h>
#include "fileparser.h"
#include "output.h"

using namespace std;

bool FileParser::read(const char* stream)
{
//  message("FileParser::read\n");
  ifstream FILE(stream);
  ostringstream OUT;
  
  if( !(FILE) ) 
    fatal("\"%s\": %s\n***ERROR*** Nemo couldn't find an input file (default is `Nemo2.ini')\n",
          stream,
          strerror(errno));
    
  message("reading parameters from \"%s\"\n",stream);

  OUT << FILE.rdbuf();
  
  FILE.close();

  return StreamParser::read(OUT.str().c_str());
}

bool BinaryFileParser::read(const char* stream)
{
//  message("reading \"%s\"\n",stream);
  
  ifstream FILE(stream);
  ostringstream OUT;
  char c;
  
  //read by character until the first separator is encountered:
  // '@' is now used in the init file for temporal parameters and external files as well
  while( (FILE.get(c)) ) {
    if(c != '@') OUT << c;
    else if ( FILE.peek() != 'G' ) OUT << c; // '@G' is first delimiter found in binary file
    else break;
  }
  
  if( !(FILE) ) fatal("BinaryFileParser::read::get %s\n",strerror(errno));
  
  FILE.close();
  
  return StreamParser::read(OUT.str().c_str());
}
