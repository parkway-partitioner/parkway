
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
                                           ostream &o)
    : GlobalCommunicator(rank, nProcs), out_stream(o) {
  numParts = nParts;
  vWeight = NULL;
  matchVector = NULL;

  numHedges = 0;
  numLocPins = 0;
  numLocalVertices = 0;
  totalVertices = 0;
  minVertexIndex = 0;
  locVertWt = 0;
  numAllocHedges = 0;

  currPercentile = 100;

  hEdgeWeight.setLength(0);
  hEdgeOffset.setLength(0);
  locPinList.setLength(0);
  vToHedgesOffset.setLength(0);
  vToHedgesList.setLength(0);
  allocHedges.setLength(0);
}

ParaHypergraphLoader::~ParaHypergraphLoader() {}

void ParaHypergraphLoader::computeHedgesToLoad(BitField &toLoad, int numH,
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

  FastDynaArray<int> hEdges(numH);
  FastDynaArray<int> hEdgeLens(numH);

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
  Funct::qsortByAnotherArray(0, numH - 1, hEdges.getArray(),
                             hEdgeLens.getArray(), INC);
  toLoad.set1();

  j = 0;
  i = 0;

  for (; i < numH && j < percentileThreshold;)
    j += hEdgeWts[hEdges[i++]];

  myPercentileLen = hEdgeLens[hEdges[i]];

  MPI_Allreduce(&myPercentileLen, &percentileLen, 1, MPI_INT, MPI_MAX, comm);

  if (dispOption > 1 && myRank == 0) {
    out_stream << " " << currPercentile << " " << maxLen << " " << aveLen << " "
               << percentileLen;
  }

  for (; i < numH; ++i)
    if (hEdgeLens[hEdges[i]] > percentileLen)
      toLoad.set0(hEdges[i]);
}

#endif
