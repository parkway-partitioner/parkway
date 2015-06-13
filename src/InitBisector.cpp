

#ifndef _INIT_BISECTOR_CPP
#define _INIT_BISECTOR_CPP

// ### InitBisector.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 12/1/2005: Last Modified
//
// ###

#include "InitBisector.hpp"

InitBisector::InitBisector(int nRuns, int insMethod, int ee, int dL)
    : FMRefiner(-1, insMethod, ee, dL) {
  numInitRuns = nRuns;
}

InitBisector::~InitBisector() {}

void InitBisector::dispInitBisectorOptions(std::ostream &out) const {
  switch (dispOption) {
  case SILENT:
    break;

  default:

    out << "|- GIB:"
        << " runs = " << numInitRuns << std::endl
        << "|" << std::endl;
    break;
  }
}

void InitBisector::initBisect(serial::hypergraph &h) {
  int i;

  int gain;
  int balanced;

  h.set_number_of_partitions(numInitRuns);

  load_for_refinement(h);
  removeEEThreshold();
  buildBuckets();

  for (i = 0; i < numInitRuns; ++i) {
    setPartitionVector(i);
    partitionCutsizes[i] = setBaseVertex();
    initPartitionStruct();
    initGains1to0();
    partitionCutsizes[i] -= doGreedyPass();

    balanced = (partWeights[0] <= maxPartWt && partWeights[1] <= maxPartWt);

    if (balanced)
      do {
        gain = doFidMatPass();
#ifdef DEBUG_REFINER
        assert(partWeights[0] <= maxPartWt);
        assert(partWeights[1] <= maxPartWt);
#endif
        partitionCutsizes[i] -= gain;
      } while (gain > 0);

    restoreBuckets();
  }

  destroyBuckets();
}

int InitBisector::setBaseVertex() {
  int i;
  int cut = 0;

  int endOffset;
  int baseVertex = RANDOM(0, numVertices);

  partWeights[1] = 0;

  for (i = 0; i < numVertices; ++i) {
    partitionVector[i] = 1;
    partWeights[1] += vWeight[i];
  }

  partitionVector[baseVertex] = 0;
  partWeights[1] -= vWeight[baseVertex];
  partWeights[0] = vWeight[baseVertex];

  endOffset = vOffsets[baseVertex + 1];

  for (i = vOffsets[baseVertex]; i < endOffset; ++i)
    cut += hEdgeWeight[vToHedges[i]];

  return cut;
}

int InitBisector::chooseBestVertex1to0() {
#ifdef DEBUG_REFINER
  assert((*bucketArrays[1])[maxGainEntries[1]].next);
#endif

  return ((*bucketArrays[1])[maxGainEntries[1]].next->vertex_id);
}

int InitBisector::doGreedyPass() {
  int gain = 0;
  int bestVertex;

  while (partWeights[1] > partWeights[0]) {
    bestVertex = chooseBestVertex1to0();

#ifdef DEBUG_REFINER
    assert(partitionVector[bestVertex] == 1);
#endif

    gain += vertexGains[bestVertex];

    removeFromBucketArray(bestVertex, 1, vertexGains[bestVertex]);
    updateGains1_0(bestVertex);
  }

  removeBucketsFrom1();

  // ###
  // if now unbalanced, do a rebalancing pass
  // ###

  if (partWeights[0] > maxPartWt)
    gain += doRebalancingPass(0);

  return gain;
}

#endif
