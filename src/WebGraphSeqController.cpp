#ifndef _WEBGRAPH_SEQ_CONTROLLER_CPP
#define _WEBGRAPH_SEQ_CONTROLLER_CPP

// ### WebGraphSeqController.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 07/04/2005: Last Modified
//
// ###

#include "WebGraphSeqController.hpp"

WebGraphSeqController::WebGraphSeqController(int rank, int nProcs, int nParts,
                                             ostream &out)
    : SeqController(rank, nProcs, nParts, out) {}

WebGraphSeqController::~WebGraphSeqController() {}

void WebGraphSeqController::initCoarsestHypergraph(ParaHypergraph &hgraph,
                                                   MPI_Comm comm) {
  /***********************

need to modify this so that we can deal with the
scenario where the number of vertices is very large
compared to the number of hyperedges and hence some
of the vertices are not connected to any others at
all.

We thus need not assemble all the unconnected vertices
into the serial hypergraph, as their presence will not
affect the overall cutsize. However, we need to 'remember'
them, so that we can assign parts of the partition to them on
exit from the serial partitioning algorithm. Need also to
devise an algorithm that will assign these vertices to parts
such that the overall partition is balanced.

Also note that the balance constraints need not be expressed
over the serial hypergraph vertices only. Need to provide a
new 'epsilon' based on the max parts weight as derived from
the weight of the entire graph, including the vertices that
have been omitted because they are not connected to any other
vertices.

  ************************/

  int i;
  int j;
  int ij;

  int recvArrayLen;
  int index;
  int numLocalVertices = hgraph.getNumLocalVertices();
  int numLocalHedges = hgraph.getNumLocalHedges();
  int numLocalPins = hgraph.getNumLocalPins();
  int localVertexWt = hgraph.getLocalVertexWt();

  int *localVertWeight = hgraph.getWeightArray();
  int *localHedgeOffsets = hgraph.getHedgeOffsetsArray();
  int *localHedgeWeights = hgraph.getHedgeWeightsArray();
  int *localPins = hgraph.getLocalPinsArray();

  int numVertices = hgraph.getNumTotalVertices();
  int numHedges;
  int numPins;
  int totVertexWt;

  DynamicArray<int> *vWeights;
  DynamicArray<int> *hEdgeWeights;
  DynamicArray<int> *hEdgeOffsets;
  DynamicArray<int> *pinList;
  DynamicArray<int> recvDispls(numProcs);
  DynamicArray<int> recvLens(numProcs);
  DynamicArray<int> recvArray;

  MPI_Allgather(&numLocalVertices, 1, MPI_INT, recvLens.getArray(), 1, MPI_INT,
                comm);

  ij = 0;
  for (i = 0; i < numProcs; ++i) {
    recvDispls[i] = ij;
    ij += recvLens[i];
  }

#ifdef DEBUG_CONTROLLER
  assert(numVertices == ij);
#endif

  vWeights = new DynamicArray<int>(numVertices);

  MPI_Allreduce(&localVertexWt, &totVertexWt, 1, MPI_INT, MPI_SUM, comm);
  MPI_Allgatherv(localVertWeight, numLocalVertices, MPI_INT,
                 vWeights->getArray(), recvLens.getArray(),
                 recvDispls.getArray(), MPI_INT, comm);
  MPI_Allgather(&numLocalHedges, 1, MPI_INT, recvLens.getArray(), 1, MPI_INT,
                comm);

  ij = 0;
  for (i = 0; i < numProcs; ++i) {
    recvDispls[i] = ij;
    ij += recvLens[i];
  }

  numHedges = ij;
  hEdgeWeights = new DynamicArray<int>(numHedges);

  MPI_Allgatherv(localHedgeWeights, numLocalHedges, MPI_INT,
                 hEdgeWeights->getArray(), recvLens.getArray(),
                 recvDispls.getArray(), MPI_INT, comm);

  ij = 0;
  for (i = 0; i < numProcs; ++i) {
    recvDispls[i] = ij;
    ++recvLens[i];
    ij += recvLens[i];
  }

  recvArrayLen = ij;
  recvArray.setLength(recvArrayLen);
  MPI_Allgatherv(localHedgeOffsets, numLocalHedges + 1, MPI_INT,
                 recvArray.getArray(), recvLens.getArray(),
                 recvDispls.getArray(), MPI_INT, comm);
  hEdgeOffsets = new DynamicArray<int>(numHedges + 1);

  j = 1;
  index = 0;
  (*hEdgeOffsets)[0] = 0;

  for (i = 1; i < recvArrayLen; ++i) {
    if (recvArray[i] != 0) {
      ij = recvArray[i] - recvArray[i - 1];
      (*hEdgeOffsets)[j] = (*hEdgeOffsets)[j - 1] + ij;
      ++j;
    } else {
      recvLens[index++] = recvArray[i - 1];
    }
  }

#ifdef DEBUG_CONTROLLER
  assert(index == numProcs - 1);
#endif

  recvLens[index] = recvArray[recvArrayLen - 1];

  ij = 0;
  for (i = 0; i < numProcs; ++i) {
    recvDispls[i] = ij;
    ij += recvLens[i];
  }

  numPins = ij;
  pinList = new DynamicArray<int>(numPins);
  MPI_Allgatherv(localPins, numLocalPins, MPI_INT, pinList->getArray(),
                 recvLens.getArray(), recvDispls.getArray(), MPI_INT, comm);

  h = new Hypergraph(vWeights->getArray(), numVertices);

  h->setNumHedges(numHedges);
  h->setNumPins(numPins);
  h->setTotWeight(totVertexWt);
  h->setHedgeWtArray(hEdgeWeights->getArray(), hEdgeWeights->getLength());
  h->setHedgeOffsetArray(hEdgeOffsets->getArray(), hEdgeOffsets->getLength());
  h->setPinListArray(pinList->getArray(), pinList->getLength());
  h->buildVtoHedges();

  if (dispOption > 0 && myRank == 0)
    h->printCharacteristics(out_stream);
}

void WebGraphSeqController::initSeqPartitions(ParaHypergraph &hgraph,
                                              MPI_Comm comm) {
  int i;
  int j;
  int ij;

  int keepMyPartition;
  int proc;
  int numKept;
  int myBestCut = h->getCut(0);
  int ijk;
  int startOffset;
  int endOffset;
  int totToSend;

  int *hPartitionVector;
  int *hPartVectorOffsets;
  int *hPartCuts;

  int numTotVertices = h->getNumVertices();
  int *pVector = h->getPartVectorArray();
  int *pCuts = h->getPartCutArray();

  DynamicArray<int> numVperProc(numProcs);
  DynamicArray<int> procDispls(numProcs);

  DynamicArray<int> sendLens(numProcs);
  DynamicArray<int> sendDispls(numProcs);
  DynamicArray<int> recvLens(numProcs);
  DynamicArray<int> recvDispls(numProcs);
  DynamicArray<int> sendArray;

  DynamicArray<int> procCuts(numProcs);
  DynamicArray<int> procs(numProcs);
  DynamicArray<int> keepPartitions(numProcs);

  // ###
  // First root processor determines
  // which partitions to keep
  // ###

  MPI_Gather(&myBestCut, 1, MPI_INT, procCuts.getArray(), 1, MPI_INT, ROOT_PROC,
             comm);

  if (myRank == ROOT_PROC) {
    numKept = 0;

    for (i = 0; i < numProcs; ++i)
      procs[i] = i;

    Funct::randomPermutation(procs.getArray(), numProcs);

    for (i = 0; i < numProcs; ++i) {
      proc = procs[i];
      ij = 0;

      for (j = i + 1; j < numProcs; ++j) {
        if (procCuts[proc] == procCuts[procs[j]]) {
          ij = 1;
          break;
        }
      }

      if (ij == 0) {
        keepPartitions[proc] = 1;
        ++numKept;
      } else
        keepPartitions[proc] = 0;
    }
  }

  MPI_Scatter(keepPartitions.getArray(), 1, MPI_INT, &keepMyPartition, 1,
              MPI_INT, ROOT_PROC, comm);
  MPI_Bcast(&numKept, 1, MPI_INT, ROOT_PROC, comm);

  hgraph.setNumberPartitions(numKept);

  hPartitionVector = hgraph.getPartitionArray();
  hPartVectorOffsets = hgraph.getPartitionOffsetsArray();
  hPartCuts = hgraph.getCutsizesArray();

  // ###
  // communicate partition vector values
  // ###

  j = numProcs - 1;
  ij = numTotVertices / numProcs;

  for (i = 0; i < j; ++i)
    numVperProc[i] = ij;

  numVperProc[i] = ij + Mod(numTotVertices, numProcs);

  j = 0;
  ij = 0;

  for (i = 0; i < numProcs; ++i) {
    sendDispls[i] = j;
    procDispls[i] = ij;
    sendLens[i] = numVperProc[i] * keepMyPartition;
    j += sendLens[i];
    ij += numVperProc[i];
  }

  sendArray.setLength(j);
  totToSend = j;

  ij = 0;

  for (ijk = 0; ijk < numProcs; ++ijk) {
    for (j = 0; j < keepMyPartition; ++j) {
      startOffset = procDispls[ijk];
      endOffset = startOffset + numVperProc[ijk];

      for (i = startOffset; i < endOffset; ++i) {
        sendArray[ij++] = pVector[i];
      }
    }
  }
#ifdef DEBUG_CONTROLLER
  assert(ij == totToSend);
#endif

  MPI_Alltoall(sendLens.getArray(), 1, MPI_INT, recvLens.getArray(), 1, MPI_INT,
               comm);

  ij = 0;

  for (i = 0; i < numProcs; ++i) {
    recvDispls[i] = ij;
    ij += recvLens[i];
  }

#ifdef DEBUG_CONTROLLER
  assert(ij == hPartVectorOffsets[numKept]);
#endif

  MPI_Alltoallv(sendArray.getArray(), sendLens.getArray(),
                sendDispls.getArray(), MPI_INT, hPartitionVector,
                recvLens.getArray(), recvDispls.getArray(), MPI_INT, comm);

  // ###
  // communicate partition cuts
  // ###

  MPI_Allgather(&keepMyPartition, 1, MPI_INT, recvLens.getArray(), 1, MPI_INT,
                comm);

  ij = 0;

  for (i = 0; i < numProcs; ++i) {
    recvDispls[i] = ij;
    ij += recvLens[i];
  }

  MPI_Allgatherv(pCuts, keepMyPartition, MPI_INT, hPartCuts,
                 recvLens.getArray(), recvDispls.getArray(), MPI_INT, comm);

  if (dispOption > 1 && myRank == 0) {
    for (i = 0; i < numKept; ++i)
      out_stream << hPartCuts[i] << " ";

    out_stream << endl;
  }
}

#endif
