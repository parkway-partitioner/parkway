// ### ParaApproxCoarsener.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 31/12/2004: Last Modified
//
// ###

#include "coarseners/parallel/approximate_coarsener.hpp"
#include <iostream>

namespace parkway {
namespace parallel {

approximate_coarsener::approximate_coarsener(int _rank, int processors,
                                         int _numParts, int percentile, int inc)
    : coarsener(_rank, processors, _numParts) {
  startPercentile = percentile;
  currPercentile = percentile;
  increment = inc;
}

approximate_coarsener::~approximate_coarsener() {
}

void approximate_coarsener::compute_hyperedges_to_load(
    ds::bit_field &toLoad, int numH, ds::dynamic_array<int> &hEdgeWts,
    ds::dynamic_array<int> &hEdgeOffsets, MPI_Comm comm) {
  int i;
  int j;

  double percentileThreshold;
  double aveLen;

  int totLen;
  int maxLen;
  int numTotHedges;
  int myTotLen = 0;
  int myMaxLen = 0;
  int myPercentileLen;
  int percentileLen;

  dynamic_array<int> hEdges(numH);
  dynamic_array<int> hEdgeLens(numH);

  /* compute the hyperedges that will not be communicated */

  for (i = 0; i < numH; ++i) {
    hEdgeLens[i] = hEdgeOffsets[i + 1] - hEdgeOffsets[i];
    myTotLen += hEdgeLens[i];

    if (hEdgeLens[i] > myMaxLen)
      myMaxLen = hEdgeLens[i];
  }

  MPI_Allreduce(&myTotLen, &totLen, 1, MPI_INT, MPI_SUM, comm);
  MPI_Allreduce(&numH, &numTotHedges, 1, MPI_INT, MPI_SUM, comm);
  MPI_Allreduce(&myMaxLen, &maxLen, 1, MPI_INT, MPI_MAX, comm);

  aveLen = static_cast<double>(totLen) / numTotHedges;

  j = 0;
  for (i = 0; i < numH; ++i) {
    hEdges[i] = i;
    j += hEdgeWts[i];
  }

  percentileThreshold = (static_cast<double>(j) * currPercentile) / 100;
  hEdges.sort_using_another_array(hEdgeLens);
  toLoad.set();

  j = 0;
  i = 0;

  for (; i < numH && j < percentileThreshold;)
    j += hEdgeWts[hEdges[i++]];

  myPercentileLen = hEdgeLens[hEdges[i]];

  MPI_Allreduce(&myPercentileLen, &percentileLen, 1, MPI_INT, MPI_MAX, comm);

  for (; i < numH; ++i)
    if (hEdgeLens[hEdges[i]] > percentileLen)
      toLoad.unset(hEdges[i]);

  /* increment the percentile for subsequent coarsening */

  currPercentile = std::min(currPercentile + increment, 100);
}

}  // namespace parallel
}  // namespace parkway
