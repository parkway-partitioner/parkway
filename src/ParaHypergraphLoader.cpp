
#ifndef _PARA_HYPO_LOADER_CPP
#define _PARA_HYPO_LOADER_CPP

// ### ParaHypergraphLoader.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 4/1/2005: Last Modified
//
// ###

#include "ParaHypergraphLoader.hpp"

ParaHypergraphLoader::ParaHypergraphLoader(int rank, int nProcs, int nParts,
                                           std::ostream &o)
    : global_communicator(rank, nProcs), out_stream(o) {
  numParts = nParts;
  vWeight = nullptr;
  matchVector = nullptr;

  numHedges = 0;
  numLocPins = 0;
  numLocalVertices = 0;
  totalVertices = 0;
  minVertexIndex = 0;
  locVertWt = 0;
  numAllocHedges = 0;

  currPercentile = 100;

  hEdgeWeight.reserve(0);
  hEdgeOffset.reserve(0);
  locPinList.reserve(0);
  vToHedgesOffset.reserve(0);
  vToHedgesList.reserve(0);
  allocHedges.reserve(0);
}

ParaHypergraphLoader::~ParaHypergraphLoader() {}

void ParaHypergraphLoader::computeHedgesToLoad(bit_field &toLoad, int numH,
                                               int numLocalPins, int *hEdgeWts,
                                               int *hEdgeOffsets,
                                               MPI_Comm comm) {
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

  percentileThreshold = (static_cast<double>(j) * currPercentile) / 100;
  Funct::qsortByAnotherArray(0, numH - 1, hEdges.data(),
                             hEdgeLens.data(), INC);
  toLoad.set();

  j = 0;
  i = 0;

  for (; i < numH && j < percentileThreshold;)
    j += hEdgeWts[hEdges[i++]];

  myPercentileLen = hEdgeLens[hEdges[i]];

  MPI_Allreduce(&myPercentileLen, &percentileLen, 1, MPI_INT, MPI_MAX, comm);

  if (dispOption > 1 && rank_ == 0) {
    out_stream << " " << currPercentile << " " << maxLen << " " << aveLen << " "
               << percentileLen;
  }

  for (; i < numH; ++i)
    if (hEdgeLens[hEdges[i]] > percentileLen)
      toLoad.unset(hEdges[i]);
}

#endif
