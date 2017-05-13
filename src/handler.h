/**  $Id: handler.h,v 1.5 2011-06-23 12:56:00 fred Exp $
*
*  @file handler.h
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
*  created on @date 16.06.2005 
*  @author fred
*/

#ifndef HANDLER_H
#define HANDLER_H

/**Service handler (an observer).
 * Interface for the FileHandler and StatHandler classes.
 * Implements the observer design pattern
 */
class Handler {    

public:
  /**Inits state. */
  virtual void init ( ) = 0;
  /**Updates the handler state.  */
  virtual void  update () = 0;
  
  virtual ~Handler ( ) { }
};
#endif //HANDLER_H

