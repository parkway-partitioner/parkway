// ### GlobalCommunicator.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 4/1/2005: Last Modified
//
// ###
#include "internal/global_communicator.hpp"

namespace parkway {

global_communicator::global_communicator(const int rank, const int processors,
                                         const int display_option)
    : rank_(rank),
      processors_(processors),
      display_option_(display_option),
      data_out_sets_(processors_),
      send_lens_(processors_),
      receive_lens_(processors_),
      send_displs_(processors_),
      receive_displs_(processors_) {
  for (auto &item : data_out_sets_) {
    item.resize(1024);
  }
}

global_communicator::~global_communicator() {
}

void global_communicator::free_memory() {
  for (auto &out_set : data_out_sets_) {
    out_set.clear_and_shrink();
  }
  send_array_.clear_and_shrink();
  receive_array_.clear_and_shrink();
}

void global_communicator::send_from_data_out(MPI_Comm comm) {
  int capacity = 0;
  for (int i = 0; i < processors_; ++i) {
    send_displs_[i] = capacity;
    capacity += send_lens_[i];
  }
  send_array_.resize(capacity);

  int ij = 0;
  for (int i = 0; i < processors_; ++i) {
    for (int j = 0; j < send_lens_[i]; ++j) {
      send_array_[ij++] = data_out_sets_[i][j];
    }
  }

  // Send data from all to all processes.
  // Send 1 int from send_lens_,
  // Receive 1 into into receive_lens_
  MPI_Alltoall(send_lens_.data(), 1, MPI_INT,
               receive_lens_.data(), 1, MPI_INT,
               comm);

  capacity = 0;
  for (int i = 0; i < processors_; ++i) {
    receive_displs_[i] = capacity;
    capacity += receive_lens_[i];
  }

  receive_array_.resize(capacity);

  MPI_Alltoallv(
      send_array_.data(), send_lens_.data(), send_displs_.data(), MPI_INT,
      receive_array_.data(), receive_lens_.data(), receive_displs_.data(),
      MPI_INT, comm);
}

}  // namespace parkway
