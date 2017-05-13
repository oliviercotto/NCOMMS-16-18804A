/** $Id: fileparser.h,v 1.5 2011-06-23 12:56:00 fred Exp $
*
*  @file fileparser.h
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
*  created on @date 09.05.2005
* 
*  @author fred
*/

#ifndef _FILEPARSER_H
#define _FILEPARSER_H

#include "paramsparser.h"
/**Text input parameter file parser.
 * This class provides the StreamParser with the whole content of the input file.
 */
class FileParser: public StreamParser {
  
public:
  FileParser(const char* stream) : StreamParser(stream) {}
  virtual ~FileParser(){}
  
  virtual bool read(const char* stream);
  
private:
  
};
/**Retrieves simulation parameters from a binary data file.
 * The parameters are extracted from the binary file before being passed the StreamParser base class.
 */
class BinaryFileParser: public StreamParser {

public:
  BinaryFileParser(const char* stream) : StreamParser(stream) {}
  virtual ~BinaryFileParser(){}
  
  virtual bool read(const char* stream);
  
private:
  
};

#endif
