// ### ParaController.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 4/1/2005: Last Modified
//
// ###
#include "internal/parallel_controller.hpp"
#include "utility/logging.hpp"

namespace parkway {
namespace parallel {

controller::controller(coarsener &c, refiner &r,
                       serial::controller &con, int rank, int nP,
                       int percentile, int inc, int approxRef)
    : global_communicator(rank, nP),
      coarsener_(c),
      refiner_(r),
      serial_controller_(con) {
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
  total_serial_time_ = 0;
  total_refinement_time_ = 0;
  total_time_ = 0;
  accumulator_ = 0;

  hypergraph_ = nullptr;

  best_partition_.reserve(0);
}

controller::~controller() {
}

void controller::initialize_map_to_orig_verts() {
  int min_vertex_index = hypergraph_->minimum_vertex_index();
  map_to_orig_vertices_.reserve(number_of_orig_local_vertices_);

  for (int i = 0; i < number_of_orig_local_vertices_; ++i) {
    map_to_orig_vertices_[i] = min_vertex_index + i;
  }
}

void controller::set_prescribed_partition(const char *filename, MPI_Comm comm) {
#ifdef DEBUG_CONTROLLER
  assert(hgraph);
#endif
  if (shuffled_ == 2) {
    int len;
    int numVPerProc = hypergraph_->total_number_of_vertices() / processors_;
    int myOffset = rank_ * numVPerProc;

    std::ifstream in_stream;
    in_stream.open(filename, std::ifstream::in | std::ifstream::binary);

    if (!in_stream.is_open()) {
      error_on_processor("[Processor %i] Could not open partition file %s\n",
                         rank_, filename);
      MPI_Abort(comm, 0);
    }

    shuffle_partition_.reserve(number_of_orig_local_vertices_);

    in_stream.seekg(myOffset * sizeof(int), std::ifstream::beg);
    len = number_of_orig_local_vertices_ * sizeof(int);
    in_stream.read((char *)(shuffle_partition_.data()), len);

    if (in_stream.gcount() != len) {
      error_on_processor("[Processor %i] Could not read in %i elements\n",
                         rank_, number_of_orig_local_vertices_);
      MPI_Abort(comm, 0);
    }

    in_stream.close();
  }
}

void controller::store_best_partition(int numV, const dynamic_array<int> array,
                                      MPI_Comm comm) {

#ifdef DEBUG_CONTROLLER
  int numLocalVertices = hgraph->getNumLocalVertices();
#endif

  int minLocVertIndex = hypergraph_->minimum_vertex_index();
  int numTotalVertices = hypergraph_->total_number_of_vertices();
  int vertPerProc = numTotalVertices / processors_;
  int totToRecv;
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

    ij = std::min(vertex / vertPerProc, processors_ - 1);
    assert(ij < processors_);
    if (ij == rank_) {
      best_partition_[vertex - minLocVertIndex] = vPart;
    } else {
      data_out_sets_[ij][send_lens_[ij]++] = vertex;
      data_out_sets_[ij][send_lens_[ij]++] = vPart;
    }
  }

  ij = 0;

  for (i = 0; i < processors_; ++i) {
    send_displs_[i] = ij;
    ij += send_lens_[i];
  }

  send_array_.reserve(ij);
  ij = 0;

  for (i = 0; i < processors_; ++i) {
    j = 0;
    sendLength = send_lens_[i];
    while (j < sendLength) {
      send_array_[ij++] = data_out_sets_[i][j++];
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

void controller::partition_to_file(const char *filename, MPI_Comm comm) const {
  std::ofstream out;

  if (rank_ == 0)
    remove(filename);

  MPI_Barrier(comm);

  for (int i = 0; i < processors_; ++i) {
    if (rank_ == i) {
      out.open(filename, std::ofstream::out | std::ofstream::app |
               std::ofstream::binary);

      if (!out.is_open()) {
        warning_on_processor(
            "[Processor %i] Cannot open %s\n", rank_, filename);
      } else {
        out.write((char *)(best_partition_.data()),
                  sizeof(int) * number_of_orig_local_vertices_);
      }

      out.close();
    }
    MPI_Barrier(comm);
  }
}

void controller::copy_out_partition(int numVertices, int *pVector) const {
  for (int i = 0; i < numVertices; ++i)
    pVector[i] = best_partition_[i];
}

void controller::set_weight_constraints(MPI_Comm comm) {
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
  serial_controller_.set_maximum_vertex_weight(maxVertWt);
}

}  // namespace parallel
}  // namespace parkway
