#ifndef _SEQ_CONTROLLER_CPP
#define _SEQ_CONTROLLER_CPP

// ### SeqController.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 30/11/2004: Last Modified
//
// ###

#include "sequential_controller.hpp"

sequential_controller::sequential_controller(int rank, int nProcs, int nParts, ostream &out)
    : out_stream_(out) {
  rank_ = rank;
  number_of_processors_ = nProcs;
  number_of_parts_ = nParts;
  number_of_runs_ = 0;
  accept_proportion_ = 0;
  display_option_ = 0;
  maximum_vertex_weight_ = 0;
  accept_proportion_ = 0;

  hypergraph_ = nullptr;

  partition_vector_.reserve(0);
  partition_vector_cuts_.reserve(0);
  partition_vector_offsets_.reserve(0);
}

sequential_controller::~sequential_controller() { dynamic_memory::delete_pointer<serial::hypergraph>(
      hypergraph_); }

void sequential_controller::initialize_coarsest_hypergraph(
    parallel::hypergraph &hgraph,
    MPI_Comm comm) {
  int i;
  int j;
  int ij;

  int recvArrayLen;
  int index;
  int numLocalVertices = hgraph.number_of_vertices();
  int numLocalHedges = hgraph.number_of_hyperedges();
  int numLocalPins = hgraph.number_of_pins();
  int localVertexWt = hgraph.vertex_weight();

  int *localVertWeight = hgraph.vertex_weights();
  int *localHedgeOffsets = hgraph.hyperedge_offsets();
  int *localHedgeWeights = hgraph.hyperedge_weights();
  int *localPins = hgraph.pin_list();

  int numVertices = hgraph.total_number_of_vertices();
  int numHedges;
  int numPins;
  int totVertexWt;

  dynamic_array<int> *vWeights;
  dynamic_array<int> *hEdgeWeights;
  dynamic_array<int> *hEdgeOffsets;
  dynamic_array<int> *pinList;
  dynamic_array<int> recvDispls(number_of_processors_);
  dynamic_array<int> recvLens(number_of_processors_);
  dynamic_array<int> recvArray;

  MPI_Allgather(&numLocalVertices, 1, MPI_INT, recvLens.data(), 1, MPI_INT,
                comm);

  ij = 0;
  for (i = 0; i < number_of_processors_; ++i) {
    recvDispls[i] = ij;
    ij += recvLens[i];
  }

#ifdef DEBUG_CONTROLLER
  assert(numVertices == ij);
#endif

  vWeights = new dynamic_array<int>(numVertices);

  MPI_Allreduce(&localVertexWt, &totVertexWt, 1, MPI_INT, MPI_SUM, comm);
  MPI_Allgatherv(localVertWeight, numLocalVertices, MPI_INT,
                 vWeights->data(), recvLens.data(),
                 recvDispls.data(), MPI_INT, comm);
  MPI_Allgather(&numLocalHedges, 1, MPI_INT, recvLens.data(), 1, MPI_INT,
                comm);

  ij = 0;
  for (i = 0; i < number_of_processors_; ++i) {
    recvDispls[i] = ij;
    ij += recvLens[i];
  }

  numHedges = ij;
  hEdgeWeights = new dynamic_array<int>(numHedges);

  MPI_Allgatherv(localHedgeWeights, numLocalHedges, MPI_INT,
                 hEdgeWeights->data(), recvLens.data(),
                 recvDispls.data(), MPI_INT, comm);

  ij = 0;
  for (i = 0; i < number_of_processors_; ++i) {
    recvDispls[i] = ij;
    ++recvLens[i];
    ij += recvLens[i];
  }

  recvArrayLen = ij;
  recvArray.reserve(recvArrayLen);
  MPI_Allgatherv(localHedgeOffsets, numLocalHedges + 1, MPI_INT,
                 recvArray.data(), recvLens.data(),
                 recvDispls.data(), MPI_INT, comm);
  hEdgeOffsets = new dynamic_array<int>(numHedges + 1);

  j = 1;
  index = 0;
  (*hEdgeOffsets)[0] = 0;

  for (i = 1; i < recvArrayLen; ++i) {
    if (recvArray[i] != 0) {
      ij = recvArray[i] - recvArray[i - 1];
      (*hEdgeOffsets)[j] = (*hEdgeOffsets)[j - 1] + ij;
      ++j;
    } else {
      recvLens[index++] = recvArray[i - 1];
    }
  }

#ifdef DEBUG_CONTROLLER
  assert(index == numProcs - 1);
#endif

  recvLens[index] = recvArray[recvArrayLen - 1];

  ij = 0;
  for (i = 0; i < number_of_processors_; ++i) {
    recvDispls[i] = ij;
    ij += recvLens[i];
  }

  numPins = ij;
  pinList = new dynamic_array<int>(numPins);
  MPI_Allgatherv(localPins, numLocalPins, MPI_INT, pinList->data(),
                 recvLens.data(), recvDispls.data(), MPI_INT, comm);

  hypergraph_ = new serial::hypergraph(vWeights->data(), numVertices);

  hypergraph_->set_number_of_hyperedges(numHedges);
  hypergraph_->set_number_of_pins(numPins);
  hypergraph_->set_total_weight(totVertexWt);
  hypergraph_->set_hyperedge_weights(hEdgeWeights->data(), hEdgeWeights->capacity());
  hypergraph_->set_hyperedge_offsets(hEdgeOffsets->data(), hEdgeOffsets->capacity());
  hypergraph_->set_pin_list(pinList->data(), pinList->capacity());
  hypergraph_->buildVtoHedges();

  if (display_option_ > 0 && rank_ == 0)
    hypergraph_->print_characteristics(out_stream_);
}

void sequential_controller::initialize_sequential_partitions(
    parallel::hypergraph &hgraph, MPI_Comm comm) {
  int i;
  int j;
  int ij;

  int keepMyPartition;
  int proc;
  int numKept;
  int myBestCut = hypergraph_->cut(0);
  int ijk;
  int startOffset;
  int endOffset;
  int totToSend;

  int *hPartitionVector;
  int *hPartVectorOffsets;
  int *hPartCuts;

  int numTotVertices = hypergraph_->number_of_vertices();
  int *pVector = hypergraph_->partition_vector();
  int *pCuts = hypergraph_->partition_cuts();

  dynamic_array<int> numVperProc(number_of_processors_);
  dynamic_array<int> procDispls(number_of_processors_);

  dynamic_array<int> sendLens(number_of_processors_);
  dynamic_array<int> sendDispls(number_of_processors_);
  dynamic_array<int> recvLens(number_of_processors_);
  dynamic_array<int> recvDispls(number_of_processors_);
  dynamic_array<int> sendArray;

  dynamic_array<int> procCuts(number_of_processors_);
  dynamic_array<int> procs(number_of_processors_);
  dynamic_array<int> keepPartitions(number_of_processors_);

  // ###
  // First root processor determines
  // which partitions to keep
  // ###

  MPI_Gather(&myBestCut, 1, MPI_INT, procCuts.data(), 1, MPI_INT, ROOT_PROC,
             comm);

  if (rank_ == ROOT_PROC) {
    numKept = 0;

    for (i = 0; i < number_of_processors_; ++i)
      procs[i] = i;

    Funct::randomPermutation(procs.data(), number_of_processors_);

    for (i = 0; i < number_of_processors_; ++i) {
      proc = procs[i];
      ij = 0;

      for (j = i + 1; j < number_of_processors_; ++j) {
        if (procCuts[proc] == procCuts[procs[j]]) {
          ij = 1;
          break;
        }
      }

      if (ij == 0) {
        keepPartitions[proc] = 1;
        ++numKept;
      } else
        keepPartitions[proc] = 0;
    }
  }

  MPI_Scatter(keepPartitions.data(), 1, MPI_INT, &keepMyPartition, 1,
              MPI_INT, ROOT_PROC, comm);
  MPI_Bcast(&numKept, 1, MPI_INT, ROOT_PROC, comm);

  hgraph.set_number_of_partitions(numKept);

  hPartitionVector = hgraph.partition_vector();
  hPartVectorOffsets = hgraph.partition_offsets();
  hPartCuts = hgraph.partition_cuts();

  // ###
  // communicate partition vector values
  // ###

  j = number_of_processors_ - 1;
  ij = numTotVertices / number_of_processors_;

  for (i = 0; i < j; ++i)
    numVperProc[i] = ij;

  numVperProc[i] = ij + Mod(numTotVertices, number_of_processors_);

  j = 0;
  ij = 0;

  for (i = 0; i < number_of_processors_; ++i) {
    sendDispls[i] = j;
    procDispls[i] = ij;
    sendLens[i] = numVperProc[i] * keepMyPartition;
    j += sendLens[i];
    ij += numVperProc[i];
  }

  sendArray.reserve(j);
  totToSend = j;

  ij = 0;

  for (ijk = 0; ijk < number_of_processors_; ++ijk) {
    for (j = 0; j < keepMyPartition; ++j) {
      startOffset = procDispls[ijk];
      endOffset = startOffset + numVperProc[ijk];

      for (i = startOffset; i < endOffset; ++i) {
        sendArray[ij++] = pVector[i];
      }
    }
  }
#ifdef DEBUG_CONTROLLER
  assert(ij == totToSend);
#endif

  MPI_Alltoall(sendLens.data(), 1, MPI_INT, recvLens.data(), 1, MPI_INT,
               comm);

  ij = 0;

  for (i = 0; i < number_of_processors_; ++i) {
    recvDispls[i] = ij;
    ij += recvLens[i];
  }

#ifdef DEBUG_CONTROLLER
  assert(ij == hPartVectorOffsets[numKept]);
#endif

  MPI_Alltoallv(sendArray.data(), sendLens.data(),
                sendDispls.data(), MPI_INT, hPartitionVector,
                recvLens.data(), recvDispls.data(), MPI_INT, comm);

  // ###
  // communicate partition cuts
  // ###

  MPI_Allgather(&keepMyPartition, 1, MPI_INT, recvLens.data(), 1, MPI_INT,
                comm);

  ij = 0;

  for (i = 0; i < number_of_processors_; ++i) {
    recvDispls[i] = ij;
    ij += recvLens[i];
  }

  MPI_Allgatherv(pCuts, keepMyPartition, MPI_INT, hPartCuts,
                 recvLens.data(), recvDispls.data(), MPI_INT, comm);

  if (display_option_ > 1 && rank_ == 0) {
    for (i = 0; i < numKept; ++i)
      out_stream_ << hPartCuts[i] << " ";

    out_stream_ << endl;
  }
}

int sequential_controller::choose_best_partition() const {
#ifdef DEBUG_CONTROLLER
  assert(numSeqRuns > 0);
  assert(partitionCuts.getLength() > 0);
#endif

  int i = 1;
  int j = 0;

  for (; i < number_of_runs_; ++i)
    if (partition_vector_cuts_[i] < partition_vector_cuts_[j])
      j = i;

  return j;
}

int sequential_controller::accept_cut() const {
  int i;
  int bestCut;

#ifdef DEBUG_CONTROLLER
  assert(h);
  assert(numSeqRuns > 0);
  assert(partitionVector.getLength() == h->getNumVertices() * numSeqRuns);
  assert(partitionCuts.getLength() == numSeqRuns);
  assert(partitionVectorOffsets.getLength() == numSeqRuns + 1);
  for (i = 0; i < numSeqRuns; ++i)
    assert(partitionCuts[i] > 0);
#endif

  bestCut = partition_vector_cuts_[0];

  for (i = 1; i < number_of_runs_; ++i) {
    if (partition_vector_cuts_[i] < bestCut) {
      bestCut = partition_vector_cuts_[i];
    }
  }

  return (static_cast<int>(
      floor(static_cast<double>(bestCut) + bestCut * accept_proportion_)));
}

#endif
