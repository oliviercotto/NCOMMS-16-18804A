/** $Id: MPImanager.cc,v 1.7.2.3 2016-02-09 14:03:05 fred Exp $
 *
 *  @file MPImanager.cc
 *  Nemo2
 *
 *  Copyright (C) 2006-2011 Frederic Guillaume
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
 *  Created on @date 07.08.2004
 *
 *  @author jacques
 */
#include <cerrno>
#include <sys/wait.h>
#include <iostream>
#include "simulation.h"
#include "metapop.h"
#include "binarydatasaver.h"

extern MPIenv *_myenv;

/*********************************************************************/

MPIenv::MPIenv( int &argc, char **&argv, MPImanager *&_p )  
{
#ifdef USE_MPI
  MPI::Init( argc, argv );
  size = MPI::COMM_WORLD.Get_size();
  rank = MPI::COMM_WORLD.Get_rank();
  int lhost;
  char *thost = new char[100];
  MPI::Get_processor_name(thost, lhost);
  host = std::string( thost, lhost );
  delete[] thost;
  if ( rank == 0 ) _p = new MPImaster(size-1);
  else             _p = new MPIslave();
#else
  size = 1;
  rank = 0;
  host = std::string( "local" );
#endif
}

void MPIenv::abort( int i ) 
{
#ifdef USE_MPI
  MPI::COMM_WORLD.Abort(i);
#endif
  exit(i);
}

void MPIenv::finish( MPImanager *p )
{
#ifdef USE_MPI
  delete p;
  MPI::Finalize();
#endif
}

/*********************************************************************/

#ifdef USE_MPI
//---------------------------------------------------------------------
// MPImanager::init
//---------------------------------------------------------------------
unsigned int MPImanager::init( StatServices *StatManager )
{
//  buf_stride = 0;
  //--- Count total number of rows (=num of generations) in all stat containers
  //std::list< StatRecBase* > all_stats = SH->get_service()->getAllStats();
//  for ( std::list< StatRecBase* >::iterator I = all_stats.begin();
//       I != all_stats.end(); 
//       I++ )
//    buf_stride += (*I)->getRows(); 
  
  //this is the number of generations recorded x number of stat recorders
  buf_stride = StatManager->getNumOccurrences() * (StatManager->getNumStats() + 2);
  //+2 is to count replicate and generation entries
  
  buf_int = new unsigned int*[2];
  return 0;
}
//---------------------------------------------------------------------
// MPImanager::finish
//---------------------------------------------------------------------
void MPImanager::finish( StatServices *StatManager, unsigned int *_gen, unsigned int _repl )
{ 
  for ( unsigned int i = 0; i < size; i++ )
    delete[] buf_dbl[i];
  delete[] buf_dbl;
  delete[] buf_int[0];
  delete[] buf_int[1];
  delete[] buf_int;
//  if (BinaryDataSaver::PID > 0) waitpid(BinaryDataSaver::PID,NULL,0);
  MPI::COMM_WORLD.Barrier(); 
}
/*********************************************************************/
/*                       M P I _ m a s t e r                         */
/*********************************************************************/
MPImaster::MPImaster( const unsigned int count )
{
  _request = new MPI::Request[count];
}
//---------------------------------------------------------------------
// MPImaster::init
//---------------------------------------------------------------------
unsigned int MPImaster::init( StatServices *StatManager )
{
#ifdef DEBUG_MPI
  message("MPI.Master: init\n");
#endif
  
  MPImanager::init( StatManager ); //sets buf_stride and buf_int
  
  size = _myenv->slaveCount();
  
  buf_dbl = new double*[size];
  
  for ( unsigned int i = 0; i < size; i++ ) {
    buf_dbl[i] = new double[buf_stride];
  
  }
  buf_int[0] = new unsigned int[size];
  buf_int[1] = new unsigned int[size];
  slave_job.resize( size );
  
  //--- Each slave is already doing a replicate
  for ( int i = 1; i <= size; i++ ) {
    slave_job[i-1] = i;
    _request[i-1] = MPI::COMM_WORLD.Irecv( buf_dbl[i-1], buf_stride, 
                                          MPI::DOUBLE, i, 0 );
  }
#ifdef DEBUG_MPI
  message("MPI.Master: init done\n");
#endif
  //--- Master's replicate loop starts at size+1!
  return size+1;
}
//---------------------------------------------------------------------
// MPImaster::finish
//---------------------------------------------------------------------
void MPImaster::finish ( StatServices *StatManager, unsigned int *_gen, unsigned int _repl )
{
#ifdef DEBUG_MPI
  message("MPI.Master: finishing\n");
#endif
  unsigned int *_buf2 = buf_int[1];
  
  //--- First, send last replicate to all slave so they exit the replicate loop
  for ( unsigned int slave = 1; slave <= size; slave++, _buf2++ ) {
    *_buf2 = _repl;
    MPI::COMM_WORLD.Isend( _buf2, 1, MPI::UNSIGNED, slave, 0 );
  }
  
  unsigned int num_stats = StatManager->getNumStats() + 2;
  unsigned int num_rec, done;
  
  //--- Second, wait for all pending communications, and fill the stat containers
  while ( !MPI::Request::Testall( size, _request ) ) {
    
    int reqindx = MPI::Request::Waitany( size, _request );
    
    if ( reqindx == MPI::UNDEFINED ) MPIenv::abort(102); 
    
    MPI::COMM_WORLD.Recv( &buf_int[0][reqindx], 1, MPI::UNSIGNED, reqindx+1, 1 );
    
    num_rec = buf_int[0][reqindx];
    done = slave_job[reqindx];
    
    for(unsigned int i = 0, pos = 0; i < num_rec && pos < buf_stride; ++i) {
      StatManager->copyGenerationStatValues(done, i, &buf_dbl[reqindx][pos], num_stats);
      pos += num_stats;
    }
    
    if(StatManager->getNumOccurrences(done) != num_rec)
      fatal("while recording slave[%i] stats into Master; container size don't match (%i != %i)\n", reqindx+1, StatManager->getNumOccurrences(done), num_rec);
    
    
//    *_gen = buf_int[0][reqindx];
//    unsigned int done = slave_job[reqindx], pos = 0, 
//    occur = SH->getOccurrence(), actualGen;
//    
//    for ( std::list< StatRecBase* >::iterator I = all_stats.begin(); 
//         I != all_stats.end(); 
//         I++ ) {
//      
//	    actualGen = occur;
//      
//	    for ( unsigned int n = 1; n <= (*I)->getRows(); n++ ) {
//        (*I)->setVal( actualGen, n, done, buf_dbl[reqindx][pos] );
//        pos++;
//        actualGen += occur;
//        if ( actualGen > *_gen ) break;
//	    }
//    }
  }
  
  
  
  //--- Free everything
  slave_job.clear();
  
  MPImanager::finish( 0, 0, 0 );

#ifdef DEBUG_MPI
  message("MPI.Master: finished\n");
#endif
}
//---------------------------------------------------------------------
// MPImaster::waitSlave
//---------------------------------------------------------------------
unsigned int MPImaster::waitSlave()
{
#ifdef DEBUG_MPI
  message("MPI.Master.waitSlave: Waiting\n" );
#endif
  int next = MPI::Request::Waitany( size, _request );
  
  if ( next == MPI::UNDEFINED ) MPIenv::abort(101); //no more active request
  
  ++next; //this is the rank, not the index
  
  //receive num of stat records done by slave 'next', recorded in buf_int
  MPI::COMM_WORLD.Recv( &buf_int[0][next-1], 1, MPI::UNSIGNED, next, 1 );
  
#ifdef DEBUG_MPI
  message("MPI.Master.waitSlave: Got slave %i\n", next);
#endif
  
  return next; //return slave number
}
//---------------------------------------------------------------------
// MPImaster::assign
//---------------------------------------------------------------------
unsigned int MPImaster::assign( const unsigned int job, const unsigned int slave )
{
#ifdef DEBUG_MPI
  message("MPI.Master.assign: Assigning replicate %i to slave %i\n", job, slave );
#endif
  unsigned int indx = slave-1;
  
  if ( indx >= slave_job.size() ) return 0;
  
  buf_int[1][indx] = job; //new replicate assigned to slave
  //send new relicate number to slave
  MPI::COMM_WORLD.Send( &buf_int[1][indx], 1, MPI::UNSIGNED, slave, 0 );
  
  unsigned int done = slave_job[indx]; //replicate that just finished
  
  slave_job[indx] = job;
  
#ifdef DEBUG_MPI
  message("MPI.Master.assign: slave %i finished repl %i, assigned repl %i\n", 
          slave, done, job);
#endif
  return done;
}
//---------------------------------------------------------------------
// MPImaster::iterate
//---------------------------------------------------------------------
void MPImaster::iterate( SimRunner *_sim, StatServices *StatManager, 
                        unsigned int *_gen, unsigned int *_repl )
{
#ifdef DEBUG_MPI
  message("MPI.Master.iterate: waiting for a slave to finish (current repl=%i, current gen=%i)\n",
          *_repl, *_gen);
#endif
  //wait for a slave to finish its replicate
  //waiting on pending requests for buf_dbl
  unsigned int slave = waitSlave();
  
  //get the generation at which the slave finished, 
  //and assign it to SimRunner::_current_generation (???is private!!!)
  //*_gen = buf_int[0][slave-1];
  
  //assign a new replicate to the slave
  unsigned int done = assign( *_repl, slave );
  
#ifdef DEBUG_MPI
  message("MPI.Master.iterate: slave %i finished replicate %i, collecting stats\n",
          slave, done);
#endif
  
  //record the stats that the slave just sent to Master in Master's stat recorders  
  unsigned int num_stats = StatManager->getNumStats() + 2;
  unsigned int num_rec   = buf_int[0][slave-1];
  
  for(unsigned int i = 0, pos = 0; i < num_rec && pos < buf_stride; ++i) {
    StatManager->copyGenerationStatValues(done, i, &buf_dbl[slave-1][pos], num_stats);
    pos += num_stats;
  }
  
  if(StatManager->getNumOccurrences(done) != num_rec)
    fatal("while recording slave[%i] stats into Master; container size don't match (%i != %i)\n",
          slave, StatManager->getNumOccurrences(done), num_rec);
  
  //issue a new non-blocking request for next replicate done by this slave
  _request[slave-1] = MPI::COMM_WORLD.Irecv( buf_dbl[slave-1], buf_stride,
                                             MPI::DOUBLE, slave, 0 );
  
  
 #ifdef DEBUG_MPI
  message("MPI.Master.iterate: replicate %i assigned, incrementing counter\n",*_repl);
#endif 
  
  //increment replicate counter, this is the only place where it is done
  ++(*_repl);
}
/*********************************************************************/
/*                         M P I _ s l a v e                         */
/*********************************************************************/
unsigned int MPIslave::init( StatServices *StatManager )
{
#ifdef DEBUG_MPI
  message("MPI.Slave: init\n");
#endif
  MPImanager::init( StatManager );
  size = 1;
  buf_dbl = new double*[1];
  buf_dbl[0] = new double[ buf_stride ];
  buf_int[0] = new unsigned int[1];
  buf_int[1] = new unsigned int[1];
#ifdef DEBUG_MPI
  message("MPI.Slave: init done\n");
#endif
  return _myenv->slaveRank();
}
//---------------------------------------------------------------------
// MPIslave::iterate
//---------------------------------------------------------------------
void MPIslave::iterate( SimRunner *_sim, StatServices *StatManager, 
                       unsigned int *_gen, unsigned int *_repl )
{
#ifdef DEBUG_MPI
  message("MPI.Slave[%i]:starting new replicate %i\n", _myenv->slaveRank(), *_repl);
#endif
  
  //--- Start a cycle
  _sim->setForFirstGeneration();
  _sim->Cycle(NULL);
  
  //check if extinction occurred, if so, call the FileServices to update the FileHandlers
  if( !_sim->get_pop()->isAlive() && _sim->getCurrentGeneration() < _sim->getGenerations() ) {
    _sim->setCurrentGeneration(_sim->getGenerations()); //propagates to metapop
    _sim->_FileServices.notify();
  }
  
#ifdef DEBUG_MPI
  message("MPI.Slave[%i]:done with replicate %i (gen %i)\n",
          _myenv->slaveRank(), *_repl, *_gen);
#endif
//  buf_int[0][0] = *_gen; 
//  //record and send the generation counter, may be less than max num generations
//  MPI::COMM_WORLD.Isend( buf_int[0], 1, MPI::UNSIGNED, 0, 1 );
  
  //--- Collect all stats
  unsigned int num_stats = StatManager->getNumStats() + 2;
  unsigned int num_rec   = StatManager->getNumOccurrences(*_repl);
  
  //record and send the number of stat records saved in the buffer
  buf_int[0][0] = num_rec; 
  MPI::COMM_WORLD.Isend( buf_int[0], 1, MPI::UNSIGNED, 0, 1 );
  
  size_t memsize = num_stats*sizeof(double);
  
  for(unsigned int i = 0, pos = 0; i < num_rec && pos < buf_stride; ++i) {
    memcpy(&buf_dbl[0][pos], StatManager->getGenerationStatValues(*_repl, i), memsize);
    pos += num_stats;
  }
  
  //--- Send stats to Master
  MPI::COMM_WORLD.Send( buf_dbl[0], buf_stride, MPI::DOUBLE, 0, 0 );

#ifdef DEBUG_MPI
  message("MPI.Slave[%i]:finished collecting and sending stats for repl %i\n",
          _myenv->slaveRank(),*_repl);
#endif

  //--- Receive next job
  MPI::COMM_WORLD.Recv( buf_int[1], 1, MPI::UNSIGNED, 0, 0 );

  *_repl = buf_int[1][0];
}

#endif
