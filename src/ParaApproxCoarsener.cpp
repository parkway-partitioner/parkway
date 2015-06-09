
#ifndef _PARA_APPROX_COARSENER_CPP
#define _PARA_APPROX_COARSENER_CPP

// ### ParaApproxCoarsener.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 31/12/2004: Last Modified
//
// ###

#include "ParaApproxCoarsener.hpp"
#include <iostream>

ParaApproxCoarsener::ParaApproxCoarsener(int _rank, int _numProcs,
                                         int _numParts, int percentile, int inc,
                                         std::ostream &out)
    : ParaCoarsener(_rank, _numProcs, _numParts, out) {
  startPercentile = percentile;
  currPercentile = percentile;
  increment = inc;
}

ParaApproxCoarsener::~ParaApproxCoarsener() {}

void ParaApproxCoarsener::loadHyperGraph(const ParaHypergraph &h,
                                         MPI_Comm comm) {
  int i;
  int vertsPerProc;
  int endOffset;
  int startOffset;
  int hEdgeLen;
  int recvLen;
  int numLocalPins;
  int numLocalHedges;

  int *localPins;
  int *localHedgeOffsets;
  int *localHedgeWeights;

  int j;
  int l;
  int proc;
  int locVert;

  DynamicArray<int> sentToProc;
  DynamicArray<int> vDegs;

  BitField toLoad;

  numLocalPins = h.getNumLocalPins();
  numLocalHedges = h.getNumLocalHedges();
  localPins = h.getLocalPinsArray();
  localHedgeOffsets = h.getHedgeOffsetsArray();
  localHedgeWeights = h.getHedgeWeightsArray();

  locVertWt = h.getLocalVertexWt();
  vWeight = h.getWeightArray();
  matchVector = h.getMatchVectorArray();

  numLocalVertices = h.getNumLocalVertices();
  totalVertices = h.getNumTotalVertices();
  minVertexIndex = h.getMinVertexIndex();
  maxVertexIndex = minVertexIndex + numLocalVertices;

  // ###
  // Prepare data structures
  // ###

  numHedges = 0;
  numLocPins = 0;
  vertsPerProc = totalVertices / numProcs;

  vToHedgesOffset.setLength(numLocalVertices + 1);
  sentToProc.setLength(numProcs);
  vDegs.setLength(numLocalVertices);
  toLoad.setLength(numLocalHedges);

  for (i = 0; i < numLocalVertices; ++i) {
    vToHedgesOffset[i] = 0;
    vDegs[i] = 0;
  }

  if (currPercentile < 100)
    computeHedgesToLoad(toLoad, numLocalHedges, localHedgeWeights,
                        localHedgeOffsets, comm);
  else
    toLoad.set1();

  // ###
  // use the request sets to send local hyperedges to other
  // processors and to receive hyperedges from processors
  // ###

  for (i = 0; i < numProcs; ++i) {
    sendLens[i] = 0;
    sentToProc[i] = 0;
  }

  for (i = 0; i < numLocalHedges; ++i) {
    if (toLoad(i) == 1) {
      startOffset = localHedgeOffsets[i];
      endOffset = localHedgeOffsets[i + 1];
      hEdgeLen = endOffset - startOffset;

      for (j = startOffset; j < endOffset; ++j) {
#ifdef DEBUG_COARSENER
        assert(localPins[j] < totalVertices && localPins[j] >= 0);
#endif
        proc = std::min(localPins[j] / vertsPerProc, numProcs - 1);

        if (!sentToProc[proc]) {
          if (proc == myRank) {
            hEdgeWeight.assign(numHedges, localHedgeWeights[i]);
            hEdgeOffset.assign(numHedges++, numLocPins);

            for (l = startOffset; l < endOffset; ++l) {
              locPinList.assign(numLocPins++, localPins[l]);
#ifdef DEBUG_COARSENER
              assert(localPins[l] < totalVertices && localPins[l] >= 0);
#endif
              if (localPins[l] >= minVertexIndex &&
                  localPins[l] < maxVertexIndex)
                ++vToHedgesOffset[localPins[l] - minVertexIndex];
            }
          } else {
            dataOutSets[proc]->assign(sendLens[proc]++, hEdgeLen + 2);
            dataOutSets[proc]->assign(sendLens[proc]++, localHedgeWeights[i]);

            for (l = startOffset; l < endOffset; ++l) {
              dataOutSets[proc]->assign(sendLens[proc]++, localPins[l]);
            }
          }

          sentToProc[proc] = 1;
        }
      }

      for (j = 0; j < numProcs; ++j)
        sentToProc[j] = 0;
    }
  }

  // ###
  // Now exchange the hyperedges
  // ###

  sendFromDataOutArrays(comm);

  // ###
  // Now load the non-local hyperedges
  // ###

  j = 0;
  for (i = 0; i < numProcs; ++i) {
    j += recvLens[i];
  }
  recvLen = j;

  j = 0;
  while (j < recvLen) {
    endOffset = j + receiveArray[j];
    ++j;

    hEdgeWeight.assign(numHedges, receiveArray[j++]);
    hEdgeOffset.assign(numHedges++, numLocPins);

    for (; j < endOffset; ++j) {
      locPinList.assign(numLocPins++, receiveArray[j]);

#ifdef DEBUG_COARSENER
      assert(receiveArray[j] < totalVertices && receiveArray[j] >= 0);
#endif

      locVert = receiveArray[j] - minVertexIndex;

      if (locVert >= 0 && locVert < numLocalVertices)
        ++vToHedgesOffset[locVert];
    }
  }

  hEdgeOffset.assign(numHedges, numLocPins);

#ifdef MEM_OPT
  hEdgeOffset.setLength(numHedges + 1);
  hEdgeWeight.setLength(numHedges);
  locPinList.setLength(numHedges);
#endif

  // ###
  // now initialise the vToHedgesList
  // ###

  j = 0;
  l = 0;

  for (; j < numLocalVertices; ++j) {
#ifdef DEBUG_COARSENER
    assert(vToHedgesOffset[j] >= 0);
#endif

    locVert = vToHedgesOffset[j];
    vToHedgesOffset[j] = l;
    l += locVert;
  }

  vToHedgesOffset[j] = l;
  vToHedgesList.setLength(l);

  for (j = 0; j < numHedges; ++j) {
    endOffset = hEdgeOffset[j + 1];

    for (l = hEdgeOffset[j]; l < endOffset; ++l) {
      if (locPinList[l] >= minVertexIndex && locPinList[l] < maxVertexIndex) {

        locVert = locPinList[l] - minVertexIndex;
        vToHedgesList[vToHedgesOffset[locVert] + (vDegs[locVert]++)] = j;
      }
    }
  }

#ifdef DEBUG_COARSENER
  for (j = 1; j <= numLocalVertices; ++j)
    assert(vDegs[j - 1] == vToHedgesOffset[j] - vToHedgesOffset[j - 1]);
#endif
}

void ParaApproxCoarsener::computeHedgesToLoad(BitField &toLoad, int numH,
                                              int *hEdgeWts, int *hEdgeOffsets,
                                              MPI_Comm comm) {
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

  DynamicArray<int> hEdges(numH);
  DynamicArray<int> hEdgeLens(numH);

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
  Funct::qsortByAnotherArray(0, numH - 1, hEdges.getArray(),
                             hEdgeLens.getArray(), INC);
  toLoad.set1();

  j = 0;
  i = 0;

  for (; i < numH && j < percentileThreshold;)
    j += hEdgeWts[hEdges[i++]];

  myPercentileLen = hEdgeLens[hEdges[i]];

  MPI_Barrier(comm);
  std::cout << "myPercentileLen = " << myPercentileLen << std::endl;
  MPI_Barrier(comm);

  MPI_Allreduce(&myPercentileLen, &percentileLen, 1, MPI_INT, MPI_MAX, comm);

  MPI_Barrier(comm);
  if (myRank == 0) {
    std::cout << "percentileLen = " << percentileLen << std::endl;
    std::cout << "maxHedgeLen = " << maxLen << std::endl;
    std::cout << "aveLen = " << aveLen << std::endl;
    std::cout << "currPercentile = " << currPercentile << std::endl;
  }

  for (; i < numH; ++i)
    if (hEdgeLens[hEdges[i]] > percentileLen)
      toLoad.set0(hEdges[i]);

  /* increment the percentile for subsequent coarsening */

  currPercentile = std::min(currPercentile + increment, 100);
}

#endif
