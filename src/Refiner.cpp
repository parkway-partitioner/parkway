

#ifndef _REFINER_CPP
#define _REFINER_CPP

// ### Refiner.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 4/1/2005: Last Modified
//
// ###

#include "Refiner.hpp"

Refiner::Refiner(int dL) : HypergraphLoader(dL) {
  maxPartWt = 0;
  numParts = 0;
  acceptProp = 0;
  avePartWt = 0;
  partitionVector = NULL;
  partWeights.setLength(0);
}

Refiner::~Refiner() {}

int Refiner::calcCutsize() const {
#ifdef DEBUG_REFINER
  assert(numParts > 0);
#endif

  int i;
  int j;

  FastDynaArray<int> spanned(numParts);

  int k_1Cut = 0;
  int endOffset;
  int vPart;
  int numSpanned;

  for (i = 0; i < numHedges; ++i) {
    endOffset = hEdgeOffsets[i + 1];
    numSpanned = 0;

    for (j = 0; j < numParts; ++j)
      spanned[j] = 0;

    for (j = hEdgeOffsets[i]; j < endOffset; ++j) {
      vPart = partitionVector[pinList[j]];
#ifdef DEBUG_REFINER
      assert(vPart >= 0 && vPart < numParts);
#endif
      if (!spanned[vPart]) {
        spanned[vPart] = 1;
        ++numSpanned;
      }
    }

    k_1Cut += ((numSpanned - 1) * hEdgeWeight[i]);
  }

  return k_1Cut;
}

#endif
