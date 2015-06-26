// ### WebGraphSeqController.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 07/04/2005: Last Modified
//
// ###
#include "web_graph_serial_controller.hpp"
#include "utility/logging.hpp"

web_graph_serial_controller::web_graph_serial_controller(
    int rank, int nProcs, int nParts)
    : parkway::serial::controller(rank, nProcs, nParts) {
}

web_graph_serial_controller::~web_graph_serial_controller() {}

void web_graph_serial_controller::initialize_coarsest_hypergraph(
    parallel::hypergraph &hgraph,
    MPI_Comm comm) {
  /***********************

need to modify this so that we can deal with the
scenario where the number of vertices is very large
compared to the number of hyperedges and hence some
of the vertices are not connected to any others at
all.

We thus need not assemble all the unconnected vertices
into the serial hypergraph, as their presence will not
affect the overall cutsize. However, we need to 'remember'
them, so that we can assign parts of the partition to them on
exit from the serial partitioning algorithm. Need also to
devise an algorithm that will assign these vertices to parts
such that the overall partition is balanced.

Also note that the balance constraints need not be expressed
over the serial hypergraph vertices only. Need to provide a
new 'epsilon' based on the max parts weight as derived from
the weight of the entire graph, including the vertices that
have been omitted because they are not connected to any other
vertices.

  ************************/

  int i;
  int j;
  int ij;

  int recvArrayLen;
  int index;
  int numLocalVertices = hgraph.number_of_vertices();
  int numLocalHedges = hgraph.number_of_hyperedges();
  int numLocalPins = hgraph.number_of_pins();
  int localVertexWt = hgraph.vertex_weight();

  auto localVertWeight = hgraph.vertex_weights();
  auto localHedgeOffsets = hgraph.hyperedge_offsets();
  auto localHedgeWeights = hgraph.hyperedge_weights();
  auto localPins = hgraph.pin_list();

  int numVertices = hgraph.total_number_of_vertices();
  int numHedges;
  int numPins;
  int totVertexWt;

  dynamic_array<int> vWeights;
  dynamic_array<int> hEdgeWeights;
  dynamic_array<int> hEdgeOffsets;
  dynamic_array<int> pinList;
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

  vWeights.resize(numVertices);
  MPI_Allreduce(&localVertexWt, &totVertexWt, 1, MPI_INT, MPI_SUM, comm);
  MPI_Allgatherv(localVertWeight.data(), numLocalVertices, MPI_INT,
                 vWeights.data(), recvLens.data(),
                 recvDispls.data(), MPI_INT, comm);
  MPI_Allgather(&numLocalHedges, 1, MPI_INT, recvLens.data(), 1, MPI_INT,
                comm);

  ij = 0;
  for (i = 0; i < number_of_processors_; ++i) {
    recvDispls[i] = ij;
    ij += recvLens[i];
  }

  numHedges = ij;
  hEdgeWeights.resize(numHedges);
  MPI_Allgatherv(localHedgeWeights.data(), numLocalHedges, MPI_INT,
                 hEdgeWeights.data(), recvLens.data(),
                 recvDispls.data(), MPI_INT, comm);

  ij = 0;
  for (i = 0; i < number_of_processors_; ++i) {
    recvDispls[i] = ij;
    ++recvLens[i];
    ij += recvLens[i];
  }

  recvArrayLen = ij;
  recvArray.reserve(recvArrayLen);
  MPI_Allgatherv(localHedgeOffsets.data(), numLocalHedges + 1, MPI_INT,
                 recvArray.data(), recvLens.data(),
                 recvDispls.data(), MPI_INT, comm);
  hEdgeOffsets.resize(numHedges + 1);

  j = 1;
  index = 0;
  hEdgeOffsets[0] = 0;

  for (i = 1; i < recvArrayLen; ++i) {
    if (recvArray[i] != 0) {
      ij = recvArray[i] - recvArray[i - 1];
      hEdgeOffsets[j] = hEdgeOffsets[j - 1] + ij;
      ++j;
    } else {
      recvLens[index++] = recvArray[i - 1];
    }
  }

  recvLens[index] = recvArray[recvArrayLen - 1];

  ij = 0;
  for (i = 0; i < number_of_processors_; ++i) {
    recvDispls[i] = ij;
    ij += recvLens[i];
  }

  numPins = ij;
  pinList.resize(numPins);
  MPI_Allgatherv(localPins.data(), numLocalPins, MPI_INT, pinList.data(),
                 recvLens.data(), recvDispls.data(), MPI_INT, comm);

  hypergraph_ = new serial::hypergraph(vWeights, numVertices);

  hypergraph_->set_number_of_hyperedges(numHedges);
  hypergraph_->set_number_of_pins(numPins);
  hypergraph_->set_total_weight(totVertexWt);
  hypergraph_->set_hyperedge_weights(hEdgeWeights);
  hypergraph_->set_hyperedge_offsets(hEdgeOffsets);
  hypergraph_->set_pin_list(pinList);
  hypergraph_->buildVtoHedges();
  hypergraph_->print_characteristics();
}

void web_graph_serial_controller::initialize_serial_partitions(
    parallel::hypergraph &hgraph,
    MPI_Comm comm) {
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

  int numTotVertices = hypergraph_->number_of_vertices();
  auto pVector = hypergraph_->partition_vector();
  auto pCuts = hypergraph_->partition_cuts();

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

  auto hPartitionVector = hgraph.partition_vector();
  auto hPartVectorOffsets = hgraph.partition_offsets();
  auto hPartCuts = hgraph.partition_cuts();

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
                sendDispls.data(), MPI_INT, hPartitionVector.data(),
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

  MPI_Allgatherv(pCuts.data(), keepMyPartition, MPI_INT, hPartCuts.data(),
                 recvLens.data(), recvDispls.data(), MPI_INT, comm);

  if (parkway::utility::status::handler::progress_enabled()) {
    for (i = 0; i < numKept; ++i) {
      progress("%i ", hPartCuts[i]);
    }
    progress("\n");
  }
}
