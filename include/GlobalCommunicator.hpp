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
 public:
  GlobalCommunicator(int rank, int nProcs);
  ~GlobalCommunicator();

  void freeMemory();
  void sendFromDataOutArrays(MPI_Comm comm);

  inline int rank() const {
    return rank_;
  }

  inline int processors() const {
    return processors_;
  }

protected:
  const int rank_;
  const int processors_;

  dynamic_array<dynamic_array<int> *> data_out_sets_;

  dynamic_array<int> send_lens_;
  dynamic_array<int> receive_lens_;
  dynamic_array<int> send_displs_;
  dynamic_array<int> receive_displs_;
  dynamic_array<int> send_array_;
  dynamic_array<int> receive_array_;
};

#endif
