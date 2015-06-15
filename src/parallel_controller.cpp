
#ifndef _PARA_CONTROLLER_CPP
#define _PARA_CONTROLLER_CPP

// ### ParaController.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 4/1/2005: Last Modified
//
// ###

#include "parallel_controller.hpp"

parallel_controller::parallel_controller(parallel_coarsener &c, parallel_refiner &r,
                               sequential_controller &con, int rank, int nP,
                               int percentile, int inc, int approxRef,
                               ostream &out)
    : global_communicator(rank, nP), out_stream_(out), coarsener_(c), refiner_(r),
      sequential_controller_(con) {
  shuffled_ = 0;
  number_of_runs_ = -1;
  total_number_of_parts_ = -1;
  display_option_ = -1;
  write_partition_to_file_ = -1;

  start_percentile_ = percentile;
  percentile_increment_ = inc;
  approximate_refine_ = approxRef;

  maximum_part_weight_ = 0;
  number_of_orig_local_vertices_ = 0;
  keep_partitions_within_ = 0;
  reduction_in_keep_threshold_ = 1.0;
  balance_constraint_ = 0;
  best_cutsize_ = 0;
  worst_cutsize_ = 0;
  average_cutsize_ = 0;
  start_time_ = 0;
  total_coarsening_time_ = 0;
  total_sequential_time_ = 0;
  total_refinement_time_ = 0;
  total_time_ = 0;
  accumulator_ = 0;

  hypergraph_ = nullptr;

  best_partition_.reserve(0);
}

parallel_controller::~parallel_controller() {}

/*
void ParaController::setShuffleFile(const char *filename)
{
  shuffleFile.reserve(strlen(filename)+1);
  strcpy(shuffleFile.data_(),filename);
}
*/

void parallel_controller::initialize_map_to_orig_verts() {
  int i;
  int j = hypergraph_->minimum_vertex_index();

  map_to_orig_vertices_.reserve(number_of_orig_local_vertices_);

  for (i = 0; i < number_of_orig_local_vertices_; ++i)
    map_to_orig_vertices_[i] = j + i;

  /*
  if(shuffled)
    {
      int *mapToOrig = hgraph->getToOrigVArray();

      for (i=0;i<numOrigLocVerts;++i)
        mapToOrigVerts[i] = mapToOrig[i];
    }
  else
    {
  */
}

void parallel_controller::set_prescribed_partition(const char *filename,
                                                   MPI_Comm comm) {
#ifdef DEBUG_CONTROLLER
  assert(hgraph);
#endif

  if (shuffled_ == 2) {
    int len;
    int numVPerProc = hypergraph_->total_number_of_vertices() / processors_;
    int myOffset = rank_ * numVPerProc;

    char message[512];
    ifstream in_stream;

    in_stream.open(filename, ifstream::in | ifstream::binary);

    if (!in_stream.is_open()) {
      sprintf(message, "p[%d] could not open partition file %s\n", rank_,
              filename);
      out_stream_ << message;
      MPI_Abort(comm, 0);
    }

    shuffle_partition_.reserve(number_of_orig_local_vertices_);

    in_stream.seekg(myOffset * sizeof(int), ifstream::beg);
    len = number_of_orig_local_vertices_ * sizeof(int);
    in_stream.read((char *)(shuffle_partition_.data()), len);

    if (in_stream.gcount() != len) {
      sprintf(message, "p[%d] could not read in %d elements\n", rank_,
              number_of_orig_local_vertices_);
      out_stream_ << message;
      MPI_Abort(comm, 0);
    }

    in_stream.close();
  }
}

void parallel_controller::store_best_partition(int numV, const int *array,
                                               MPI_Comm comm) {

#ifdef DEBUG_CONTROLLER
  int numLocalVertices = hgraph->getNumLocalVertices();
#endif

  int minLocVertIndex = hypergraph_->minimum_vertex_index();
  int numTotalVertices = hypergraph_->total_number_of_vertices();
  int vertPerProc = numTotalVertices / processors_;
  int totToRecv;
  int totToSend;
  int sendLength;
  int vertex;
  int vPart;

  int i;
  int j;
  int ij;

#ifdef DEBUG_CONTROLLER
  assert(numLocalVertices == numV);
  assert(mapToOrigVerts.getLength() == numV);
  if (rank_ != processors_ - 1)
    assert(numLocalVertices == vertPerProc);
#endif

  int *mapToHgraphVerts = map_to_orig_vertices_.data();
  int *auxArray;

  for (i = 0; i < processors_; ++i)
    send_lens_[i] = 0;

  best_partition_.reserve(numV);

  for (i = 0; i < numV; ++i) {
    vertex = mapToHgraphVerts[i];
    vPart = array[i];

#ifdef DEBUG_CONTROLLER
    assert(vertex >= 0 && vertex < numTotalVertices);
    assert(vPart >= 0 && vPart < numTotalParts);
#endif

    ij = min(vertex / vertPerProc, processors_ - 1);
    assert(ij < processors_);
    if (ij == rank_) {
      best_partition_[vertex - minLocVertIndex] = vPart;
    } else {
      data_out_sets_[ij]->assign(send_lens_[ij]++, vertex);
      data_out_sets_[ij]->assign(send_lens_[ij]++, vPart);
    }
  }

  ij = 0;

  for (i = 0; i < processors_; ++i) {
    send_displs_[i] = ij;
    ij += send_lens_[i];
  }

  send_array_.reserve(ij);
  totToSend = ij;
  ij = 0;

  for (i = 0; i < processors_; ++i) {
    j = 0;
    sendLength = send_lens_[i];
    auxArray = data_out_sets_[i]->data();

    while (j < sendLength) {
      send_array_[ij++] = auxArray[j++];
    }
  }

#ifdef DEBUG_CONTROLLER
  assert(ij == totToSend);
#endif

  // ###
  // get dimension and carry
  // out the communication
  // ###

  MPI_Alltoall(send_lens_.data(), 1, MPI_INT, receive_lens_.data(), 1, MPI_INT,
               comm);

  ij = 0;
  for (i = 0; i < processors_; ++i) {
    receive_displs_[i] = ij;
    ij += receive_lens_[i];
  }

  receive_array_.reserve(ij);
  totToRecv = ij;

  MPI_Alltoallv(send_array_.data(), send_lens_.data(),
                send_displs_.data(), MPI_INT, receive_array_.data(),
                receive_lens_.data(), receive_displs_.data(), MPI_INT, comm);

  // ###
  // now initialise the bestPartition data_
  // using the data_ in receive_array_
  // ###

  for (i = 0; i < totToRecv;) {
    vertex = receive_array_[i++];
    vPart = receive_array_[i++];

#ifdef DEBUG_CONTROLLER
    assert(vertex >= minLocVertIndex && vertex < minLocVertIndex + numV);
    assert(vPart >= 0 && vPart < numTotalParts);
#endif

    best_partition_[vertex - minLocVertIndex] = vPart;
  }

#ifdef DEBUG_CONTROLLER
  for (i = 0; i < numV; ++i)
    assert(bestPartition[i] >= 0 && bestPartition[i] < numTotalParts);
#endif
}

void parallel_controller::partition_to_file(const char *filename,
                                            MPI_Comm comm) const {
  int i;

  char message[512];
  ofstream out;

  if (rank_ == 0)
    remove(filename);

  MPI_Barrier(comm);

  for (i = 0; i < processors_; ++i) {
    if (rank_ == i) {
      out.open(filename, ofstream::out | ofstream::app | ofstream::binary);

      if (!out.is_open()) {
        sprintf(message, "p[%d] cannot open %s\n", rank_, filename);
        out_stream_ << message;
      } else
        out.write((char *)(best_partition_.data()),
                  sizeof(int) * number_of_orig_local_vertices_);

      out.close();
    }
    MPI_Barrier(comm);
  }
}

void parallel_controller::copy_out_partition(int numVertices,
                                             int *pVector) const {
  int i;

  for (i = 0; i < numVertices; ++i)
    pVector[i] = best_partition_[i];
}

void parallel_controller::set_weight_constraints(MPI_Comm comm) {
#ifdef DEBUG_CONTROLLER
  assert(hgraph);
  assert(numTotalParts != -1);
  assert(balConstraint > 0.0 && balConstraint < 1.0);
#endif

  int locGraphWt;
  int totGraphWt;
  int maxVertWt;

  double avePartWt;

  locGraphWt = hypergraph_->vertex_weight();

  MPI_Allreduce(&locGraphWt, &totGraphWt, 1, MPI_INT, MPI_SUM, comm);

  avePartWt = static_cast<double>(totGraphWt) / total_number_of_parts_;
  maximum_part_weight_ = static_cast<int>(floor(avePartWt + avePartWt *
                                                            balance_constraint_));
  maxVertWt = static_cast<int>(floor(avePartWt * balance_constraint_));

  coarsener_.set_maximum_vertex_weight(maxVertWt);
  coarsener_.set_total_hypergraph_weight(totGraphWt);
  sequential_controller_.set_maximum_vertex_weight(maxVertWt);
}

#ifdef DEBUG_TABLES
void ParaController::printHashMemUse() {
  // write_log(rank_, "HedgeIndex diff: %d",
  // HedgeIndexEntry::getNumAllocated()-HedgeIndexEntry::getNumDeleted());
  // write_log(rank_, "MatchEntry diff: %d",
  // MatchEntry::getNumAllocated()-MatchEntry::getNumDeleted());
  write_log(rank_, "MatchRequestEntry diff: %d",
            MatchRequestEntry::getNumAllocated() -
                MatchRequestEntry::getNumDeleted());
  // write_log(rank_, "IndexEntry diff: %d",
  // IndexEntry::getNumAllocated()-IndexEntry::getNumDeleted());
  // write_log(rank_, "ConnVertData diff: %d",
  // ConnVertData::getNumAllocated()-ConnVertData::getNumDeleted());
  // write_log(rank_, "VertexPartEntry diff: %d",
  // VertexPartEntry::getNumAllocated()-VertexPartEntry::getNumDeleted());
  // write_log(rank_, "DuplRemEntry diff: %d",
  // DuplRemEntry::getNumAllocated()-DuplRemEntry::getNumDeleted());
  // write_log(rank_, "PVectorEntry diff: %d",
  // PVectorEntry::getNumAllocated()-PVectorEntry::getNumDeleted());
}
#endif

#endif
