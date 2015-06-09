
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
#include "data_structures/dynamic_array.hpp"

using parkway::data_structures::DynamicArray;

class GlobalCommunicator {
 protected:
  const int myRank;
  const int numProcs;

  DynamicArray<DynamicArray<int> *> dataOutSets;

  DynamicArray<int> sendLens;
  DynamicArray<int> recvLens;
  DynamicArray<int> sendDispls;
  DynamicArray<int> recvDispls;
  DynamicArray<int> sendArray;
  DynamicArray<int> receiveArray;

 public:
  GlobalCommunicator(int rank, int nProcs);
  ~GlobalCommunicator();

  void freeMemory();
  void sendFromDataOutArrays(MPI_Comm comm);

  inline int getMyRank() const {
    return myRank;
  }

  inline int getNumProcs() const {
    return numProcs;
  }
};

#endif
