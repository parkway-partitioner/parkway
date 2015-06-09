
#ifndef _GLOBAL_COMMUNICATOR_CPP
#define _GLOBAL_COMMUNICATOR_CPP

// ### GlobalCommunicator.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 4/1/2005: Last Modified
//
// ###

#include "GlobalCommunicator.hpp"

GlobalCommunicator::GlobalCommunicator(const int rank, const int nProcs)
    : myRank(rank), numProcs(nProcs) {

  dataOutSets.reserve(numProcs);
  sendLens.reserve(numProcs);
  recvLens.reserve(numProcs);
  sendDispls.reserve(numProcs);
  recvDispls.reserve(numProcs);

  for (int i = 0; i < numProcs; ++i)
    dataOutSets[i] = new dynamic_array<int>(1024);
}

GlobalCommunicator::~GlobalCommunicator() {
  for (int i = 0; i < numProcs; ++i) {
    DynaMem<dynamic_array<int> >::deletePtr(dataOutSets[i]);
  }
}

void GlobalCommunicator::freeMemory() {
  for (int i = 0; i < numProcs; ++i) {
    dataOutSets[i]->reserve(0);
  }

  sendArray.reserve(0);
  receiveArray.reserve(0);
}

void GlobalCommunicator::sendFromDataOutArrays(MPI_Comm comm) {
  int j = 0;
  for (int i = 0; i < numProcs; ++i) {
    sendDispls[i] = j;
#ifdef DEBUG_BASICS
    assert(sendLens[i] >= 0);
#endif
    j += sendLens[i];
  }

  sendArray.reserve(j);

  int ij = 0;
  int *array;
  for (int i = 0; i < numProcs; ++i) {
    array = dataOutSets[i]->data();
    for (j = 0; j < sendLens[i]; ++j) {
      sendArray[ij++] = array[j];
    }
  }

  MPI_Alltoall(sendLens.data(), 1, MPI_INT, recvLens.data(), 1, MPI_INT,
               comm);

  j = 0;
  for (int i = 0; i < numProcs; ++i) {
    recvDispls[i] = j;
    j += recvLens[i];
  }

  receiveArray.reserve(j);

  MPI_Alltoallv(sendArray.data(), sendLens.data(),
                sendDispls.data(), MPI_INT, receiveArray.data(),
                recvLens.data(), recvDispls.data(), MPI_INT, comm);
}

#endif
