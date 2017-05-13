/** $Id: MPImanager.h,v 1.6.2.2 2014-04-29 15:55:56 fred Exp $
*
*  @file MPImanager.h
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
#ifndef __MPIMANAGER_H
#define __MPIMANAGER_H

#include <queue>
#include <string>
#include <map>

#ifdef USE_MPI
#include <mpi.h>
#endif

/**
   The abstract manager class.
   Provides generic function interfaces and communication buffers.
   size is the number of buffers, buf_stride the length of each of them
 **/
class SimRunner;
class StatServices;

class MPImanager {

 public:
    MPImanager() {}
    virtual ~MPImanager() {}
    virtual unsigned int init( StatServices* StatManager ) = 0;
    virtual void finish( StatServices* StatManager, unsigned int *_gen, unsigned int _repl );
    virtual void iterate( SimRunner *_sim, StatServices* StatManager, unsigned int *_gen, unsigned int *_repl ) = 0;

 protected:
    double **buf_dbl; //receives stat records from each slave
    unsigned int **buf_int; //is [2][size]: replicate [0] and generation [1] done by each slave
    unsigned int buf_stride, size; //size = slave count

};

#ifdef USE_MPI
/**
   The master class should be instantiated only once: on the master node.
   Each call to Iterate does:
     1) waitSlave(): waits for a slave to finish and returns its index,
     2) assign( repl, s ): assigns the next replicate to that slave,
     3) collects alls stats sent by that slave and inserts them into the 
        stat containers,
     4) increments the replicate counter.
   The vector slave_job keeps track of running jobs: 
       slave_job[s]=replicate being run by slave s+1
   The array _request contains the pointers to pending MPI communications,
 **/
class MPImaster : public MPImanager {

 public:
    MPImaster( const unsigned int count );
    ~MPImaster()                   { delete[] _request; }
    unsigned int init( StatServices *StatManager );
    void finish( StatServices *StatManager, unsigned int *_gen, unsigned int _repl );
    void iterate( SimRunner *_sim, StatServices* StatManager, unsigned int *_gen, unsigned int *_repl );

 private:
    MPI::Request *_request;
    std::vector< unsigned int > slave_job; //records the current replicate done by slave i
    unsigned int waitSlave();
    unsigned int assign( const unsigned int job, const unsigned int slave );
    
};

/**
   The slave class should be instantiated on every slave node.
   Each call to Iterate does:
     1) setPopulation() and Cycle() in the Metapop caller,
     2) Collects all stats into the buffers and sends them to Master,
     3) Waits for the next assignment from Master.
   Each slave starts a cycle immediately, before receiving any assignment 
     from Master: therefore at least SlaveCount() replicates will be performed.
 **/

class MPIslave : public MPImanager {

 public:
    MPIslave() {}
    unsigned int init( StatServices *StatManager );
    void iterate( SimRunner *_sim, StatServices *StatManager, unsigned int *_gen, unsigned int *_repl );

};
#endif

/**
   MPI environment setup.
   Provides basic node information: rank (0 is master, 1...size-1 are slaves), hostname.
   Constructor creates the MPImanagers (Master or Slave).
   This can be instantiated also in a single cpu situation, will return trivial values.
   The abort() function should always be used instead of exit().
 **/
class MPIenv {

 public:
    MPIenv( int &argc, char **&argv, MPImanager *&_p );
    static void  abort( int i );
    static void  finish( MPImanager *p ); 
    bool         isMaster() const   { return (rank==0); }
    unsigned int slaveCount() const { return size-1; }
    unsigned int slaveRank() const  { return rank; }
    std::string  hostName() const   { return host; }

 private:
    unsigned int size;
    unsigned int rank;
    std::string  host;

};

#endif
