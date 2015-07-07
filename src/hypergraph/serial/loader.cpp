// ### HypergraphLoader.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 4/1/2005: Last Modified
//
// ###

#include "hypergraph/serial/loader.hpp"

namespace parkway {
namespace serial {

loader::loader() {
  currPercentile = 100;
  numHedges = 0;
  numVertices = 0;
  numPins = 0;
  numPartitions = 0;
}

loader::~loader() {}

void loader::compute_hyperedges_to_load(bit_field &toLoad) {
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
  hEdges.sort_using_another_array(hEdgeLens);

  j = 0;
  i = 0;

  for (; i < numHedges && j < percentileThreshold;)
    j += hEdgeWeight[hEdges[i++]];

  percentileLen = hEdgeLens[hEdges[i]];

  for (; i < numHedges; ++i)
    if (hEdgeLens[hEdges[i]] > percentileLen)
      toLoad.unset(hEdges[i]);
}

}  // serial
}  // parkway
