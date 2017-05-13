/**  $Id: statrecorder.h,v 1.1 2014-01-24 10:53:01 fred Exp $
 *
 *  @file statrecorder.h
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
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA.
 *
 *  Created on @date 01.07.2003
 * 
 *  @author fred
 */
#ifndef __STAT_REC_H__
#define __STAT_REC_H__

#include <string>
#include "types.h"
#include "output.h"

/**Base class for the StatRecorder's, declares the interface to record stat values.
 **/
class StatRecBase {
  
private:
  /**The title of the stat recorder, longer and more explicite than the name.*/
  std::string  _title;
  /**Name of the stat, should be short (20 char) and R compliant (no '-', '+', ' ')*/
  std::string  _name;
  /**A argument to be passed to one of the function variable stored in the StatRecorder structure.*/
  unsigned int _arg1;
  unsigned int _arg2;
  /**The age class for which this stat applies.*/
  age_t        _age;
  
public:
	
  StatRecBase  ( ) : _title(""), _name(""), _arg1(0), _arg2(0), _age(ALL) { }

  
  virtual ~StatRecBase ( ) { }

  /**Sets the recorder attributes.
   * @param T the stat title (not printed in output)
   * @param N the stat name (headers in the output file)
   * @param AGE the age class for which the stat will be recorded
   * @param ARG1 the argument to pass to the S function
   * @param ARG2 the argument to pass to the S function
   **/
  void set     (std::string T, std::string N, age_t AGE, 
                unsigned int ARG1, unsigned int ARG2)
  {
    _title = T;
    _name = N;
    _age = AGE;
    _arg1 = ARG1;
    _arg2 = ARG2;
  }  
  
  ///@name Accessors
  ///@{
  void setName                       (std::string N)            {_name = N;}
  std::string  getTitle              ( )                        {return _title;}
  std::string  getName               ( )                        {return _name;}
  age_t        getAge                ( )                        {return _age;}
  unsigned int getArg1               ( )                        {return _arg1;}
  unsigned int getArg2               ( )                        {return _arg2;}
  
  ///@}
  /**Stores the value in the vector following the ordering option. 
   @param crnt_gen the currently running generation of the simulation
   @param gen_cntr the current generation number devided by the stat logging time (see StatRecorder::setVal)
   @param rpl_cntr the current replicate
   @param value the stat value to store
   */  
  virtual double setVal (age_t AGE) = 0;
};

/**Stores the pointers to the StatHandler's stat functions.*/
template <class S> class StatRecorder : public StatRecBase
{
private:
  /**Pointer to a 'stat getter' function of S using no argument.*/
  double (S::* _getStat)     (void);
  /**Pointer to a 'stat getter' function of S using a single unsigned int argument.**/
  double (S::* _getStatOneArg)   (unsigned int);
  /**Pointer to a 'stat getter' function of S using two unsigned int arguments.**/
  double (S::* _getStatTwoArg)   (unsigned int, unsigned int);
  /**Pointer to a 'stat setter' function of S using no argument.
   * A setter function only sets some inner variables subsequently fetched by getter's
   **/
  void   (S::* _setStat)     (void);
  /**Pointer to the owner of this recorder.*/
  S* _myHandler;
  
public:
  
  StatRecorder() : _getStat(0), _getStatOneArg(0), _getStatTwoArg(0),
  _setStat(0), _myHandler(0) {}

  /**@brief Sets the recorder attributes
   * @param T the stat title
   * @param N the stat name (headers in the output file)
   * @param AGE age on which the stat should be processed
   * @param ARG1 the frist argument to pass to the S function
   * @param ARG2 the second argument to pass to the S function
   * @param getNoArg function ptr to a S getter, taking no arguments
   * @param getOneArg function ptr to a S getter with a single uint argument
   * @param getTwoArg function ptr to a S getter with two uint arguments
   * @param setStat ptr to a setter function in class S
   **/
  void set(std::string T, std::string N, age_t AGE, unsigned int ARG1, unsigned int ARG2,
           double (S::* getNoArg)(void), double(S::* getOneArg)(unsigned int),
           double(S::* getTwoArg)(unsigned int,unsigned int), void(S::* setStat)(void));
  
  /**Sets the pointer to the StatHandler that owns this recorder.*/
  void setHandler (S* theHandler) {_myHandler = theHandler;}
  
  /**@brief Calls the linked stat function and returns the result
   * @return the stat value that will be recorded in the output file
   * @param AGE age class on which the stat should be processed
   **/
  virtual double setVal (age_t AGE);
  
};
/* _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/ */
//
//                               ****** StatRecorder ******

// ----------------------------------------------------------------------------------------
// StatRecorder::set()
// ----------------------------------------------------------------------------------------
template<class S> void StatRecorder<S>::set(std::string T, std::string N, age_t AGE, 
                                            unsigned int ARG1,
                                            unsigned int ARG2,
                                            double (S::* getNoArg) (void),
                                            double (S::* getOneArg)(unsigned int),
                                            double (S::* getTwoArg)(unsigned int,unsigned int),
                                            void   (S::* setStat)  (void) )
{
  StatRecBase::set(T, N, AGE, ARG1, ARG2);
  
  _getStat = getNoArg;
  _getStatOneArg = getOneArg;
  _getStatTwoArg = getTwoArg;
  _setStat = setStat;
}
// ----------------------------------------------------------------------------------------
// StatRecorder::setVal
// ----------------------------------------------------------------------------------------
template<class S> double StatRecorder<S>::setVal(age_t AGE)
{
  double statValue = 0;

  if(getAge() & AGE) {
#ifdef _DEBUG_
    message(" %s",getName().c_str());
#endif
    //run _setStat if exists
    if(_setStat != 0)  (_myHandler->*_setStat)();
    
    //get the value:
    if(_getStat != 0)
      
      statValue = (_myHandler->*_getStat)();
    
    else if(_getStatOneArg != 0)
      
      statValue = (_myHandler->*_getStatOneArg)(getArg1());
    
    else if(_getStatTwoArg != 0)
      
      statValue = (_myHandler->*_getStatTwoArg)(getArg1(), getArg2());
    
    else {

      fatal("StatRecorder \"%s\" has no _getStat funct ptr !!\n", getName().c_str());
      
    }
  }
  
  return statValue;
}

#endif

