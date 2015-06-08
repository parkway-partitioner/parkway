#ifndef _GLOBAL_COMMUNICATOR_HPP
#define _GLOBAL_COMMUNICATOR_HPP

// ### GlobalCommunicator.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 4/1/2005: Last Modified
//
// ###

#include "mpi.h"
#include "Funct.hpp"
#include "Dyna.hpp"

using namespace std;

class GlobalCommunicator {
protected:
  int myRank;
  int numProcs;

  FastDynaArray<FastDynaArray<int> *> dataOutSets;

  FastDynaArray<int> sendLens;
  FastDynaArray<int> recvLens;
  FastDynaArray<int> sendDispls;
  FastDynaArray<int> recvDispls;
  FastDynaArray<int> sendArray;
  FastDynaArray<int> receiveArray;

public:
  GlobalCommunicator(int rank, int nProcs);
  ~GlobalCommunicator();

  void freeMemory();
  void sendFromDataOutArrays(MPI_Comm comm);

  inline int getMyRank() const { return myRank; }
  inline int getNumProcs() const { return numProcs; }
};

#endif
