
#ifndef _PARA_RESTR_COARSENER_CPP
#define _PARA_RESTR_COARSENER_CPP

// ### ParaRestrCoarsener.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 4/1/2005: Last Modified
//
// ###

#include "ParaRestrCoarsener.hpp"

ParaRestrCoarsener::ParaRestrCoarsener(int rank, int nProcs, int nParts,
                                       ostream &out)
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
  numPartitions = 0;

  partitionVector = NULL;
  partitionVectorOffsets = NULL;
  partitionCuts = NULL;
  clusterWeights = NULL;
  pVector = NULL;
}

ParaRestrCoarsener::~ParaRestrCoarsener() {}

void ParaRestrCoarsener::loadHyperGraph(const ParaHypergraph &h,
                                        MPI_Comm comm) {
  int i;
  int endOffset;
  int startOffset;
  int hEdgeLen;
  int recvLen;
  int numLocalPins;
  int numLocalHedges;

  int *localPins;
  int *localHedgeOffsets;
  int *localHedgeWeights;

  register int j;
  register int l;
  register int proc;
  register int locVertex;

  FastDynaArray<int> procs(numProcs);
  FastDynaArray<int> vDegs;
  FastDynaArray<int> minLocIndices(numProcs);
  FastDynaArray<int> maxLocIndices(numProcs);
  FastDynaArray<int> vPerProc(numProcs);

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

  // here add partition information

  numPartitions = h.getNumPartitions();
  partitionVector = h.getPartitionArray();
  partitionVectorOffsets = h.getPartitionOffsetsArray();
  partitionCuts = h.getCutsizesArray();

#ifdef DEBUG_COARSENER
  assert(numPartitions == 1);
#endif

  // end partition information

  // ###
  // Prepare data structures
  // ###

  MPI_Allgather(&minVertexIndex, 1, MPI_INT, minLocIndices.getArray(), 1,
                MPI_INT, comm);
  MPI_Allgather(&maxVertexIndex, 1, MPI_INT, maxLocIndices.getArray(), 1,
                MPI_INT, comm);

  for (i = 0; i < numProcs; ++i)
    procs[i] = i;

  CompleteBinaryTree<int> vToProc(procs.getArray(), minLocIndices.getArray(),
                                  numProcs);

  numHedges = 0;
  numLocPins = 0;

  vToHedgesOffset.setLength(numLocalVertices + 1);
  vDegs.setLength(numLocalVertices);

  for (i = 0; i < numLocalVertices; ++i) {
    vToHedgesOffset[i] = 0;
    vDegs[i] = 0;
  }

  if (dispOption > 1 && myRank == 0)
    out_stream << "[PRFCC]";

  // ###
  // use the dataOutSets to send local hyperedges to other
  // processors and to receive hyperedges from processors
  // ###

  for (i = 0; i < numProcs; ++i) {
    sendLens[i] = 0;
    vPerProc[i] = 0;
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
        proc = vToProc.getRootVal(localPins[j]);
        ++vPerProc[proc];
      }

      for (j = 0; j < numProcs; ++j) {
        if (vPerProc[j] > 1) {
          if (j == myRank) {
            hEdgeWeight.assign(numHedges, localHedgeWeights[i]);
            hEdgeOffset.assign(numHedges++, numLocPins);

            for (l = startOffset; l < endOffset; ++l) {
#ifdef DEBUG_COARSENER
              assert(localPins[l] < totalVertices && localPins[l] >= 0);
#endif
              locVertex = localPins[l] - minVertexIndex;

              if (locVertex >= 0 && locVertex < numLocalVertices) {
                locPinList.assign(numLocPins++, locVertex);
                ++vToHedgesOffset[locVertex];
              }
            }
          } else {
            dataOutSets[j]->assign(sendLens[j]++, vPerProc[j] + 2);
            dataOutSets[j]->assign(sendLens[j]++, localHedgeWeights[i]);

            for (l = startOffset; l < endOffset; ++l) {
              if (localPins[l] >= minLocIndices[j] &&
                  localPins[l] < maxLocIndices[j])
                dataOutSets[j]->assign(sendLens[j]++,
                                       localPins[l] - minLocIndices[j]);
            }
          }

          vPerProc[j] = 0;
        }

        if (vPerProc[j] == 1)
          vPerProc[j] = 0;
      }
    }
  } else {
#ifdef DEBUG_COARSENER
    assert(currPercentile > 0 && currPercentile < 100);
#endif

    BitField toLoad(numLocalHedges);

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
          proc = vToProc.getRootVal(localPins[j]);
          ++vPerProc[proc];
        }

        for (j = 0; j < numProcs; ++j) {
          if (vPerProc[j] > 1) {
            if (j == myRank) {
              hEdgeWeight.assign(numHedges, localHedgeWeights[i]);
              hEdgeOffset.assign(numHedges++, numLocPins);

              for (l = startOffset; l < endOffset; ++l) {
#ifdef DEBUG_COARSENER
                assert(localPins[l] < totalVertices && localPins[l] >= 0);
#endif
                locVertex = localPins[l] - minVertexIndex;

                if (locVertex >= 0 && locVertex < numLocalVertices) {
                  locPinList.assign(numLocPins++, locVertex);
                  ++vToHedgesOffset[locVertex];
                }
              }
            } else {
              dataOutSets[j]->assign(sendLens[j]++, vPerProc[j] + 2);
              dataOutSets[j]->assign(sendLens[j]++, localHedgeWeights[i]);

              for (l = startOffset; l < endOffset; ++l) {
                if (localPins[l] >= minLocIndices[j] &&
                    localPins[l] < maxLocIndices[j])
                  dataOutSets[j]->assign(sendLens[j]++,
                                         localPins[l] - minLocIndices[j]);
              }
            }

            vPerProc[j] = 0;
          }

          if (vPerProc[j] == 1)
            vPerProc[j] = 0;
        }
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
#ifdef DEBUG_COARSENER
    assert(recvLen >= endOffset);
#endif
    ++j;

    hEdgeWeight.assign(numHedges, receiveArray[j++]);
    hEdgeOffset.assign(numHedges++, numLocPins);

    for (; j < endOffset; ++j) {
      locVertex = receiveArray[j];
#ifdef DEBUG_COARSENER
      assert(locVertex < numLocalVertices && locVertex >= 0);
#endif
      locPinList.assign(numLocPins++, locVertex);
      ++vToHedgesOffset[locVertex];
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

    locVertex = vToHedgesOffset[j];
    vToHedgesOffset[j] = l;
    l += locVertex;
  }

  vToHedgesOffset[j] = l;
  vToHedgesList.setLength(l);

  for (j = 0; j < numHedges; ++j) {
    endOffset = hEdgeOffset[j + 1];

    for (l = hEdgeOffset[j]; l < endOffset; ++l) {
      locVertex = locPinList[l];
      vToHedgesList[vToHedgesOffset[locVertex] + (vDegs[locVertex]++)] = j;
    }
  }

#ifdef DEBUG_COARSENER
  for (j = 1; j <= numLocalVertices; ++j)
    assert(vDegs[j - 1] == vToHedgesOffset[j] - vToHedgesOffset[j - 1]);
#endif
}

ParaHypergraph *ParaRestrCoarsener::contractHyperedges(ParaHypergraph &h,
                                                       MPI_Comm comm) {
  ParaHypergraph *coarseGraph =
      new ParaHypergraph(myRank, numProcs, clusterIndex, totalClusters,
                         myMinCluIndex, stopCoarsening, partitionCuts[0],
                         clusterWeights->getArray(), pVector->getArray());

  clusterWeights = NULL;
  pVector = NULL;

  h.contractRestrHyperedges(*coarseGraph, comm);
  h.setNumberPartitions(0);

  if (dispOption > 1) {
    int numTotCoarseVerts = coarseGraph->getNumTotalVertices();
    int numLocHedges = coarseGraph->getNumLocalHedges();
    int numLocPins = coarseGraph->getNumLocalPins();
    int numTotHedges;
    int numTotPins;

    MPI_Reduce(&numLocHedges, &numTotHedges, 1, MPI_INT, MPI_MAX, 0, comm);
    MPI_Reduce(&numLocPins, &numTotPins, 1, MPI_INT, MPI_MAX, 0, comm);

    if (myRank == 0) {
      out_stream << numTotCoarseVerts << " " << numTotHedges << " "
                 << numTotPins << " " << endl;
    }
  }

  return coarseGraph;
}

#endif
