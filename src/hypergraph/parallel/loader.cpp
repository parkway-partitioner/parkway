// ### ParaHypergraphLoader.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 4/1/2005: Last Modified
//
// ###

#include "hypergraph/parallel/loader.hpp"

namespace parkway {
namespace parallel {

loader::loader(int rank, int number_of_processors, int number_of_parts,
               std::ostream &out_stream, int display_option)
    : global_communicator(rank, number_of_processors),
      out_stream(out_stream),
      number_of_parts_(number_of_parts),
      number_of_hyperedges_(0),
      number_of_local_pins_(0),
      number_of_local_vertices_(0),
      number_of_vertices_(0),
      minimum_vertex_index_(0),
      maximum_vertex_index_(0),
      local_vertex_weight_(0),
      number_of_allocated_hyperedges_(0),
      display_options_(display_option),
      percentile_(100),
      vertex_weights_(nullptr),
      match_vector_(nullptr) {
  hyperedge_weights_.reserve(0);
  hyperedge_offsets_.reserve(0);
  local_pin_list_.reserve(0);
  vertex_to_hyperedges_offset_.reserve(0);
  vertex_to_hyperedges_.reserve(0);
  allocated_hyperedges_.reserve(0);
}

loader::~loader() {
}

void loader::compute_hyperedges_to_load(bit_field &toLoad, int numH,
                                        int numLocalPins, int *hEdgeWts,
                                        int *hEdgeOffsets, MPI_Comm comm) {
  int i;
  int j;

  double percentileThreshold;
  double aveLen;

  int maxLen;
  int myMaxLen = 0;
  int myPercentileLen;
  int percentileLen;

  /* 0 = pinNumber, 1 = numH */

  int myData[2];
  int hGraphData[2];

  dynamic_array<int> hEdges(numH);
  dynamic_array<int> hEdgeLens(numH);

  /* compute the hyperedges that will not be communicated */

  for (i = 0; i < numH; ++i) {
    hEdgeLens[i] = hEdgeOffsets[i + 1] - hEdgeOffsets[i];
    hEdges[i] = i;

    if (hEdgeLens[i] > myMaxLen)
      myMaxLen = hEdgeLens[i];
  }

  myData[0] = numLocalPins;
  myData[1] = numH;

  MPI_Allreduce(&myData[0], &hGraphData[0], 2, MPI_INT, MPI_SUM, comm);
  MPI_Allreduce(&myMaxLen, &maxLen, 1, MPI_INT, MPI_MAX, comm);

  aveLen = static_cast<double>(hGraphData[0]) / hGraphData[1];

  j = 0;
  for (i = 0; i < numH; ++i)
    j += hEdgeWts[i];

  percentileThreshold = (static_cast<double>(j) * percentile_) / 100;
  Funct::qsortByAnotherArray(0, numH - 1, hEdges.data(),
                             hEdgeLens.data(), INC);
  toLoad.set();

  j = 0;
  i = 0;

  for (; i < numH && j < percentileThreshold;)
    j += hEdgeWts[hEdges[i++]];

  myPercentileLen = hEdgeLens[hEdges[i]];

  MPI_Allreduce(&myPercentileLen, &percentileLen, 1, MPI_INT, MPI_MAX, comm);

  if (display_options_ > 1 && rank_ == 0) {
    out_stream << " " << percentile_ << " " << maxLen << " " << aveLen << " "
               << percentileLen;
  }

  for (; i < numH; ++i)
    if (hEdgeLens[hEdges[i]] > percentileLen)
      toLoad.unset(hEdges[i]);
}

}  // namespace parallel
}  // namespace parkway
