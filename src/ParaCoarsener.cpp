
#ifndef _PARA_COARSENER_CPP
#define _PARA_COARSENER_CPP

// ### ParaCoarsener.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 15/4/2004:  Optimisation: changed the vertices
//             and hedges structures into pin-list
//
// 25/4/2004:  Introduced base class ParaHypergraphLoader
//
// 4/1/2005: Last Modified
//
// ###

#include "ParaCoarsener.hpp"
#include <iostream>

ParaCoarsener::ParaCoarsener(int rank, int nProcs, int nParts, std::ostream &out)
    : ParaHypergraphLoader(rank, nProcs, nParts, out) {
  totalHypergraphWt = 0;
  maxVertexWt = 0;
  minNodes = 0;
  stopCoarsening = 0;
  clusterIndex = 0;
  totalClusters = 0;
  myMinCluIndex = 0;
  dispOption = 0;
  reductionRatio = 0;
  balConstraint = 0;

  clusterWeights.setLength(0);
}

ParaCoarsener::~ParaCoarsener() {}

void ParaCoarsener::loadHyperGraph(const ParaHypergraph &h, MPI_Comm comm) {
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

  dynamic_array<int> sentToProc;
  dynamic_array<int> vDegs;

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

  for (i = 0; i < numLocalVertices; ++i) {
    vToHedgesOffset[i] = 0;
    vDegs[i] = 0;
  }

  if (dispOption > 1 && myRank == 0)
    out_stream << "[PFCC]";

  // ###
  // use the request sets to send local hyperedges to other
  // processors and to receive hyperedges from processors
  // ###

  for (i = 0; i < numProcs; ++i) {
    sendLens[i] = 0;
    sentToProc[i] = 0;
  }

  if (currPercentile == 100) {
    for (i = 0; i < numLocalHedges; ++i) {
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
  } else if (currPercentile > 0) {
#ifdef DEBUG_COARSENER
    assert(currPercentile > 0 && currPercentile < 100);
#endif

    bit_field toLoad(numLocalHedges);

    computeHedgesToLoad(toLoad, numLocalHedges, numLocalPins, localHedgeWeights,
                        localHedgeOffsets, comm);

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
  } else {
    /* compute a fixed limit on hyperedges to be communicated
       - first try maxLen/2
    */

    int maxHedgeLen = 0;
    int maxTotHedgeLen;
    int limit;

    for (i = 0; i < numLocalHedges; ++i)
      if (localHedgeOffsets[i + 1] - localHedgeOffsets[i] > maxHedgeLen)
        maxHedgeLen = localHedgeOffsets[i + 1] - localHedgeOffsets[i];

    MPI_Allreduce(&maxHedgeLen, &maxTotHedgeLen, 1, MPI_INT, MPI_MAX, comm);

    limit = maxTotHedgeLen / 2;

    for (i = 0; i < numLocalHedges; ++i) {
      startOffset = localHedgeOffsets[i];
      endOffset = localHedgeOffsets[i + 1];
      hEdgeLen = endOffset - startOffset;

      if (hEdgeLen < limit) {
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
  locPinList.setLength(numLocPins);
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

ParaHypergraph *ParaCoarsener::contractHyperedges(ParaHypergraph &h,
                                                  MPI_Comm comm) {
  int i;

  int totalToRecv = 0;
  int numMyClusters;
  int clustersPerProc = totalClusters / numProcs;

  if (myRank != numProcs - 1)
    numMyClusters = clustersPerProc;
  else
    numMyClusters = clustersPerProc + Mod(totalClusters, numProcs);

  dynamic_array<int> minClusterIndex(numProcs);
  dynamic_array<int> maxClusterIndex(numProcs);
  dynamic_array<int> *clusterWts = new dynamic_array<int>(numMyClusters);

  for (i = 0; i < numProcs; ++i) {
    if (i == 0) {
      minClusterIndex[i] = 0;
      maxClusterIndex[i] = clustersPerProc;
    } else {
      minClusterIndex[i] = maxClusterIndex[i - 1];
      if (i == numProcs - 1)
        maxClusterIndex[i] = totalClusters;
      else
        maxClusterIndex[i] = minClusterIndex[i] + clustersPerProc;
    }
  }

  for (i = 0; i < numProcs; ++i) {
    if (i == 0)
      sendDispls[i] = 0;
    else
      sendDispls[i] = sendDispls[i - 1] + sendLens[i - 1];

    sendLens[i] =
        std::max(clusterIndex -
                (std::max(clusterIndex + myMinCluIndex - maxClusterIndex[i], 0) +
                 std::max(minClusterIndex[i] - myMinCluIndex, 0)),
            0);
  }

  MPI_Alltoall(sendLens.getArray(), 1, MPI_INT, recvLens.getArray(), 1, MPI_INT,
               comm);

  for (i = 0; i < numProcs; ++i) {
    recvDispls[i] = totalToRecv;
    totalToRecv += recvLens[i];
  }

#ifdef DEBUG_COARSENER
  assert(totalToRecv == numMyClusters);
#endif

  MPI_Alltoallv(clusterWeights.getArray(), sendLens.getArray(),
                sendDispls.getArray(), MPI_INT, clusterWts->getArray(),
                recvLens.getArray(), recvDispls.getArray(), MPI_INT, comm);

  myMinCluIndex = minClusterIndex[myRank];

  ParaHypergraph *coarseGraph =
      new ParaHypergraph(myRank, numProcs, numMyClusters, totalClusters,
                         myMinCluIndex, stopCoarsening, clusterWts->getArray());

  h.contractHyperedges(*coarseGraph, comm);

  if (dispOption > 1) {
    int numTotCoarseVerts = coarseGraph->getNumTotalVertices();
    int numLocHedges = coarseGraph->getNumLocalHedges();
    int numLocPins = coarseGraph->getNumLocalPins();
    int numTotCoarseHedges;
    int numTotCoarsePins;

    MPI_Reduce(&numLocHedges, &numTotCoarseHedges, 1, MPI_INT, MPI_SUM, 0,
               comm);
    MPI_Reduce(&numLocPins, &numTotCoarsePins, 1, MPI_INT, MPI_SUM, 0, comm);

    if (myRank == 0) {
      out_stream << numTotCoarseVerts << " " << numTotCoarseHedges << " "
                 << numTotCoarsePins << " " << std::endl;
    }
  }

  return coarseGraph;
}

#endif
