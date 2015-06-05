
#ifndef _SEQ_CONTROLLER_CPP
#define _SEQ_CONTROLLER_CPP

// ### SeqController.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 30/11/2004: Last Modified
//
// ###

#include "SeqController.hpp"

SeqController::SeqController(int rank, int nProcs, int nParts, ostream &out)
    : out_stream(out) {
  myRank = rank;
  numProcs = nProcs;
  numParts = nParts;
  numSeqRuns = 0;
  acceptProp = 0;
  dispOption = 0;
  maxVertexWt = 0;
  acceptProp = 0;

  h = NULL;

  partitionVector.setLength(0);
  partitionCuts.setLength(0);
  partitionVectorOffsets.setLength(0);
}

SeqController::~SeqController() { DynaMem<Hypergraph>::deletePtr(h); }

void SeqController::initCoarsestHypergraph(ParaHypergraph &hgraph,
                                           MPI_Comm comm) {
  register int i;
  register int j;
  register int ij;

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

  FastDynaArray<int> *vWeights;
  FastDynaArray<int> *hEdgeWeights;
  FastDynaArray<int> *hEdgeOffsets;
  FastDynaArray<int> *pinList;
  FastDynaArray<int> recvDispls(numProcs);
  FastDynaArray<int> recvLens(numProcs);
  FastDynaArray<int> recvArray;

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

  vWeights = new FastDynaArray<int>(numVertices);

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
  hEdgeWeights = new FastDynaArray<int>(numHedges);

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
  hEdgeOffsets = new FastDynaArray<int>(numHedges + 1);

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
  pinList = new FastDynaArray<int>(numPins);
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

void SeqController::initSeqPartitions(ParaHypergraph &hgraph, MPI_Comm comm) {
  register int i;
  register int j;
  register int ij;

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

  FastDynaArray<int> numVperProc(numProcs);
  FastDynaArray<int> procDispls(numProcs);

  FastDynaArray<int> sendLens(numProcs);
  FastDynaArray<int> sendDispls(numProcs);
  FastDynaArray<int> recvLens(numProcs);
  FastDynaArray<int> recvDispls(numProcs);
  FastDynaArray<int> sendArray;

  FastDynaArray<int> procCuts(numProcs);
  FastDynaArray<int> procs(numProcs);
  FastDynaArray<int> keepPartitions(numProcs);

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

int SeqController::chooseBestPartition() const {
#ifdef DEBUG_CONTROLLER
  assert(numSeqRuns > 0);
  assert(partitionCuts.getLength() > 0);
#endif

  register int i = 1;
  register int j = 0;

  for (; i < numSeqRuns; ++i)
    if (partitionCuts[i] < partitionCuts[j])
      j = i;

  return j;
}

int SeqController::getAcceptCut() const {
  register int i;
  register int bestCut;

#ifdef DEBUG_CONTROLLER
  assert(h);
  assert(numSeqRuns > 0);
  assert(partitionVector.getLength() == h->getNumVertices() * numSeqRuns);
  assert(partitionCuts.getLength() == numSeqRuns);
  assert(partitionVectorOffsets.getLength() == numSeqRuns + 1);
  for (i = 0; i < numSeqRuns; ++i)
    assert(partitionCuts[i] > 0);
#endif

  bestCut = partitionCuts[0];

  for (i = 1; i < numSeqRuns; ++i) {
    if (partitionCuts[i] < bestCut) {
      bestCut = partitionCuts[i];
    }
  }

  return (static_cast<int>(
      floor(static_cast<double>(bestCut) + bestCut * acceptProp)));
}

#endif
