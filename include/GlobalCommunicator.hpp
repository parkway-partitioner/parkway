
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

using parkway::data_structures::dynamic_array;

class GlobalCommunicator {
 protected:
  const int myRank;
  const int numProcs;

  dynamic_array<dynamic_array<int> *> dataOutSets;

  dynamic_array<int> sendLens;
  dynamic_array<int> recvLens;
  dynamic_array<int> sendDispls;
  dynamic_array<int> recvDispls;
  dynamic_array<int> sendArray;
  dynamic_array<int> receiveArray;

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
