
#ifndef _RESTR_COARSENER_CPP
#define _RESTR_COARSENER_CPP

// ### RestrCoarsener.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 4/1/2005: Last Modified
//
// ###

#include "RestrCoarsener.hpp"

RestrCoarsener::RestrCoarsener(int min, int maxWt, double ratio, int dispL)
    : HypergraphLoader(dispL) {
  minNodes = min;
  maxVertexWt = maxWt;
  reductionRatio = ratio;
}

RestrCoarsener::~RestrCoarsener() {}

Hypergraph *RestrCoarsener::buildCoarseHypergraph(int *coarseWts,
                                                  int *coarsePartVector,
                                                  int numCoarseVerts,
                                                  int totWt) const {
  int numNewHedges = 0;
  int numNewPins = 0;
  int startOff;
  int newHedgeLen;
  int duplDeg;
  int startV;
  int numDupls;
  int hEdge;
  int duplHedge;
  int startHedgeOffset;
  int endHedgeOffset;
  int hEdgeStart;
  int hEdgeEnd;
  int v;

  int i;
  int j;
  int ij;

  Hypergraph *newHypergraph = new Hypergraph(
      coarseWts, coarsePartVector, numCoarseVerts, partitionCutsizes[0]);

  DynamicArray<int> *newHedgeOffsets = new DynamicArray<int>(1024);
  DynamicArray<int> *newPinList = new DynamicArray<int>(1024);
  DynamicArray<int> *newHedgeWt = new DynamicArray<int>(1024);
  DynamicArray<int> *newVerOffsets =
      new DynamicArray<int>(numCoarseVerts + 1);
  DynamicArray<int> *newVtoHedges = new DynamicArray<int>(1024);

  DynamicArray<int> tempPinList(numPins);
  DynamicArray<int> vDegs(numCoarseVerts);
  DynamicArray<int> duplDegs(numCoarseVerts);
  DynamicArray<int> vHedgOffsets(numCoarseVerts + 1);
  DynamicArray<int> vHedges;

  for (i = 0; i < numCoarseVerts; ++i) {
    vDegs[i] = 0;
    duplDegs[i] = 0;
  }

  // ###
  // here there are two schemes for hyperedge contraction
  // one runs in |E|h^2/2 time
  // other in |E|h(2+log h)
  // pick latter for now
  // ###

  for (i = 0; i < numHedges; ++i) {
    endHedgeOffset = hEdgeOffsets[i + 1];
    startHedgeOffset = hEdgeOffsets[i];

    for (j = startHedgeOffset; j < endHedgeOffset; ++j)
      tempPinList[j] = matchVector[pinList[j]];

    Funct::qsort(startHedgeOffset, endHedgeOffset - 1, tempPinList.getArray());

    ++vDegs[tempPinList[startHedgeOffset]];
  }

  vHedgOffsets[0] = 0;

  for (i = 1; i <= numCoarseVerts; ++i)
    vHedgOffsets[i] = vHedgOffsets[i - 1] + vDegs[i - 1];

  vHedges.setLength(vHedgOffsets[numCoarseVerts]);

  for (i = 0; i < numCoarseVerts; ++i)
    vDegs[i] = 0;

  // ###
  // build the new pin list
  // ###

  newHedgeOffsets->assign(numNewHedges, numNewPins);
  numDupls = 0;

  for (i = 0; i < numHedges; ++i) {
    endHedgeOffset = hEdgeOffsets[i + 1];
    startHedgeOffset = numNewPins;

    for (j = hEdgeOffsets[i]; j < endHedgeOffset; ++j) {
      v = tempPinList[j];

      if ((*newHedgeOffsets)[numNewHedges] == numNewPins ||
          v != (*newPinList)[numNewPins - 1]) {
        newPinList->assign(numNewPins++, v);
      }
    }

    newHedgeLen = numNewPins - (*newHedgeOffsets)[numNewHedges];

    if (newHedgeLen > 1) {
      startV = (*newPinList)[startHedgeOffset];
      duplHedge = -1;
      duplDeg = duplDegs[startV];
      startOff = vHedgOffsets[startV];

      for (j = 0; j < duplDeg; ++j) {
        hEdge = vHedges[startOff + j];
        hEdgeEnd = (*newHedgeOffsets)[hEdge + 1];
        hEdgeStart = (*newHedgeOffsets)[hEdge];

        if (hEdgeEnd - hEdgeStart == newHedgeLen) {
          duplHedge = hEdge;
          for (ij = 1; ij < newHedgeLen; ++ij)
            if ((*newPinList)[hEdgeStart + ij] !=
                (*newPinList)[startHedgeOffset + ij]) {
              duplHedge = -1;
              break;
            }

          if (duplHedge == hEdge)
            break;
        }
      }

      if (duplHedge == -1) {
        vHedges[vHedgOffsets[startV] + duplDeg] = numNewHedges;
        ++duplDegs[startV];

        for (j = (*newHedgeOffsets)[numNewHedges]; j < numNewPins; ++j) {
          v = (*newPinList)[j];
          ++vDegs[v];
        }

        newHedgeWt->assign(numNewHedges++, hEdgeWeight[i]);
        newHedgeOffsets->assign(numNewHedges, numNewPins);
      } else {
        ++numDupls;
        (*newHedgeWt)[duplHedge] += hEdgeWeight[i];
        numNewPins = (*newHedgeOffsets)[numNewHedges];
      }
    } else
      numNewPins = (*newHedgeOffsets)[numNewHedges];
  }

  // ###
  // build the new vToHedges
  // ###

  (*newVerOffsets)[0] = 0;

  for (i = 1; i < numCoarseVerts + 1; ++i) {
    ij = i - 1;
    (*newVerOffsets)[i] = (*newVerOffsets)[ij] + vDegs[ij];
    vDegs[ij] = 0;
  }

#ifdef PRUDENT
  assert((*newVerOffsets)[numCoarseVerts] == numNewPins);
#endif

  newVtoHedges->setLength(numNewPins);

  for (i = 0; i < numNewHedges; ++i) {
    endHedgeOffset = (*newHedgeOffsets)[i + 1];

    for (j = (*newHedgeOffsets)[i]; j < endHedgeOffset; ++j) {
      v = (*newPinList)[j];
      (*newVtoHedges)[(*newVerOffsets)[v] + (vDegs[v]++)] = i;
    }
  }

  // ###
  // init new hypergraph
  // ###

  newHedgeWt->setLength(numNewHedges);
  newPinList->setLength(numNewPins);
  newHedgeOffsets->setLength(numNewHedges + 1);

  newHypergraph->setNumHedges(numNewHedges);
  newHypergraph->setNumPins(numNewPins);
  newHypergraph->setTotWeight(totWt);
  newHypergraph->setHedgeWtArray(newHedgeWt->getArray(), numNewHedges);
  newHypergraph->setPinListArray(newPinList->getArray(), numNewPins);
  newHypergraph->setHedgeOffsetArray(newHedgeOffsets->getArray(),
                                     numNewHedges + 1);
  newHypergraph->setVtoHedgesArray(newVtoHedges->getArray(), numNewPins);
  newHypergraph->setVoffsetsArray(newVerOffsets->getArray(),
                                  numCoarseVerts + 1);

  return newHypergraph;
}

#endif
