// ### GlobalCommunicator.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 4/1/2005: Last Modified
//
// ###

#include "global_communicator.hpp"

global_communicator::global_communicator(const int rank, const int nProcs)
    : rank_(rank), processors_(nProcs) {

  data_out_sets_.reserve(processors_);
  send_lens_.reserve(processors_);
  receive_lens_.reserve(processors_);
  send_displs_.reserve(processors_);
  receive_displs_.reserve(processors_);

  for (int i = 0; i < processors_; ++i) {
    data_out_sets_[i] = new dynamic_array<int>(1024);
  }
}

global_communicator::~global_communicator() {
  for (int i = 0; i < processors_; ++i) {
    dynamic_memory::delete_pointer<dynamic_array<int> >(data_out_sets_[i]);
  }
}

void global_communicator::free_memory() {
  for (int i = 0; i < processors_; ++i) {
    data_out_sets_[i]->reserve(0);
  }

  send_array_.reserve(0);
  receive_array_.reserve(0);
}

void global_communicator::send_from_data_out(MPI_Comm comm) {
  int j = 0;
  for (int i = 0; i < processors_; ++i) {
    send_displs_[i] = j;
#ifdef DEBUG_BASICS
    assert(send_lens_[i] >= 0);
#endif
    j += send_lens_[i];
  }

  send_array_.reserve(j);

  int ij = 0;
  int *array;
  for (int i = 0; i < processors_; ++i) {
    array = data_out_sets_[i]->data();
    for (j = 0; j < send_lens_[i]; ++j) {
      send_array_[ij++] = array[j];
    }
  }

  MPI_Alltoall(send_lens_.data(), 1, MPI_INT, receive_lens_.data(), 1, MPI_INT,
               comm);

  j = 0;
  for (int i = 0; i < processors_; ++i) {
    receive_displs_[i] = j;
    j += receive_lens_[i];
  }

  receive_array_.reserve(j);

  MPI_Alltoallv(send_array_.data(), send_lens_.data(),
                send_displs_.data(), MPI_INT, receive_array_.data(),
                receive_lens_.data(), receive_displs_.data(), MPI_INT, comm);
}
