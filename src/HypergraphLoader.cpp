#ifndef _LOADER_CPP
#define _LOADER_CPP

// ### HypergraphLoader.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 4/1/2005: Last Modified
//
// ###

#include "HypergraphLoader.hpp"

HypergraphLoader::HypergraphLoader(int disp) {
  dispOption = disp;
  currPercentile = 100;

  numHedges = 0;
  numVertices = 0;
  numPins = 0;
  numPartitions = 0;

  vWeight = nullptr;
  hEdgeWeight = nullptr;
  matchVector = nullptr;
  pinList = nullptr;
  hEdgeOffsets = nullptr;
  vToHedges = nullptr;
  vOffsets = nullptr;
  partitionVectors = nullptr;
  partitionOffsets = nullptr;
  partitionCutsizes = nullptr;
}

HypergraphLoader::~HypergraphLoader() {}

void HypergraphLoader::computeHedgesToLoad(bit_field &toLoad) {
  int i;
  int j = 0;

  int percentileLen;
  double percentileThreshold;

  dynamic_array<int> hEdgeLens(numHedges);
  dynamic_array<int> hEdges(numHedges);

  for (i = 0; i < numHedges; ++i) {
    hEdges[i] = i;
    hEdgeLens[i] = hEdgeOffsets[i + 1] - hEdgeOffsets[i];
    j += hEdgeWeight[i];
  }

  percentileThreshold = (static_cast<double>(j) * currPercentile) / 100;
  Funct::qsortByAnotherArray(0, numHedges - 1, hEdges.getArray(),
                             hEdgeLens.getArray(), INC);

  j = 0;
  i = 0;

  for (; i < numHedges && j < percentileThreshold;)
    j += hEdgeWeight[hEdges[i++]];

  percentileLen = hEdgeLens[hEdges[i]];

  for (; i < numHedges; ++i)
    if (hEdgeLens[hEdges[i]] > percentileLen)
      toLoad.set0(hEdges[i]);
}

#endif
