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

namespace parkway {
namespace ds = data_structures;

class global_communicator {
 public:
  global_communicator(int rank, int processors, int display_option = 0);
  ~global_communicator();

  void free_memory();
  void send_from_data_out(MPI_Comm comm);

  inline int rank() const {
    return rank_;
  }

  inline int processors() const {
    return processors_;
  }

protected:
  const int rank_;
  const int processors_;
  const int display_option_;

  ds::dynamic_array<ds::dynamic_array<int> > data_out_sets_;
  ds::dynamic_array<int> send_lens_;
  ds::dynamic_array<int> receive_lens_;
  ds::dynamic_array<int> send_displs_;
  ds::dynamic_array<int> receive_displs_;
  ds::dynamic_array<int> send_array_;
  ds::dynamic_array<int> receive_array_;
};

}  // namespace parkway

#endif
