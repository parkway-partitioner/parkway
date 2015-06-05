
#ifndef _RECUR_BISECT_CONTROLLER_CPP
#define _RECUR_BISECT_CONTROLLER_CPP

// ### RecurBisectController.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 4/1/2005: Last Modified
//
// ###

#include "RecurBisectController.hpp"

RecurBisectController::RecurBisectController(BisectionController *b,
                                             GreedyKwayRefiner *k, int rank,
                                             int nProcs, int nParts,
                                             int nBisectRuns, ostream &out)
    : SeqController(rank, nProcs, nParts, out) {
  bisector = b;
  kWayRefiner = k;
  numBisectRuns = nBisectRuns;
  logK = 0;
  maxPartWt = 0;
  sumOfCuts = 0;
  locVertPartInfoLen = 0;
  numMyPartitions = 0;
  bisectConstraint = 0;
  avePartWt = 0;
  aveInitBisectionWt = 0;

  locVertPartitionInfo.setLength(0);
  allPartitionInfo.setLength(0);

#ifdef DEBUG_CONTROLLER
  assert(bisector);
#endif
}

RecurBisectController::~RecurBisectController() {
  DynaMem<BisectionController>::deletePtr(bisector);
  DynaMem<GreedyKwayRefiner>::deletePtr(kWayRefiner);
}

void RecurBisectController::dispSeqControllerOptions() const {
  switch (dispOption) {
  case SILENT:
    break;

  default:

    out_stream << "|--- SEQ_CON:" << endl
               << "|- RBis:"
               << " seqR = " << numSeqRuns << " bisR = " << numBisectRuns
               << " pkT = " << acceptProp << endl
               << "|" << endl;
#ifdef DEBUG_CONTROLLER
    assert(bisector);
#endif
    bisector->dispBisectionControllerOptions();

    break;
  }
}

void RecurBisectController::convToBisectionConstraints() {
#ifdef DEBUG_CONTROLLER
  assert(h);
#endif

  logK = Funct::log2(numParts);

  register int j;
  register int i;

  double a = (1.0 + kWayConstraint) / numParts;
  double b = 1.0 / logK;

  bisectConstraint = pow(a, b) - 0.5;

  j = h->getTotWeight();

  avePartWt = static_cast<double>(j) / numParts;
  aveInitBisectionWt = static_cast<double>(j) / 2;
  maxPartWt = static_cast<int>(floor(avePartWt + avePartWt * kWayConstraint));

  kWayRefiner->setMaxPartWt(maxPartWt);
  kWayRefiner->setAvePartWt(avePartWt);

// ###
// now initialise the partition
// vector structures
// ###

#ifdef DEBUG_CONTROLLER
  assert(numSeqRuns > 0);
#endif

  // ###
  // now determine how many of the sequential
  // runs' partitions  the processor should
  // k-way refine
  // ###

  if (numSeqRuns <= numProcs) {
    if (myRank < numSeqRuns)
      numMyPartitions = 1;
    else
      numMyPartitions = 0;
  } else {
    i = Mod(numSeqRuns, numProcs);
    j = numSeqRuns / numProcs;

    if (myRank < i)
      numMyPartitions = j + 1;
    else
      numMyPartitions = j;
  }

  partitionCuts.setLength(numMyPartitions);
  partitionVectorOffsets.setLength(numMyPartitions + 1);
  partitionVectorOffsets[0] = 0;

  j = h->getNumVertices();

  for (i = 1; i <= numMyPartitions; ++i) {
    partitionVectorOffsets[i] = partitionVectorOffsets[i - 1] + j;
  }

  partitionVector.setLength(partitionVectorOffsets[numMyPartitions]);
}

void RecurBisectController::runSeqPartitioner(ParaHypergraph &hgraph,
                                              MPI_Comm comm) {
  initCoarsestHypergraph(hgraph, comm);
  convToBisectionConstraints();

  if (dispOption > 1 && myRank == 0) {
    out_stream << "[R-B]: " << numSeqRuns << " | ";
  }

  register int i;
  register int j;
  register int ij;

  int numVertices = h->getNumVertices();
  int *pVector = NULL;
  int destProcessor;
  int myPartitionIdx = 0;
  int v;

  FastDynaArray<int> recvLens(numProcs);
  FastDynaArray<int> recvDispls(numProcs);

  Bisection *b;

  allPartitionInfo.setLength(Shiftl(numVertices, 1));

  for (i = 0; i < numSeqRuns; ++i) {
    destProcessor = Mod(i, numProcs);
    sumOfCuts = 0;
    locVertPartInfoLen = 0;

    if (myRank == destProcessor) {
#ifdef DEBUG_CONTROLLER
      assert(myPartitionIdx < numMyPartitions);
#endif
      pVector = &partitionVector[partitionVectorOffsets[myPartitionIdx]];
    }

    b = new Bisection(h, logK, 0);
    b->initMap();

    recursivelyBisect(*b, comm);

    // ###
    // now recover the partition and
    // partition cutsize
    // ###

    MPI_Reduce(&sumOfCuts, &ij, 1, MPI_INT, MPI_SUM, destProcessor, comm);
    MPI_Gather(&locVertPartInfoLen, 1, MPI_INT, recvLens.getArray(), 1, MPI_INT,
               destProcessor, comm);

    if (myRank == destProcessor) {
      partitionCuts[myPartitionIdx] = ij;
      ij = 0;

      for (j = 0; j < numProcs; ++j) {
        recvDispls[j] = ij;
        ij += recvLens[j];
      }
#ifdef DEBUG_CONTROLLER
      assert(ij == Shiftl(numVertices, 1));
#endif
    }

    MPI_Gatherv(locVertPartitionInfo.getArray(), locVertPartInfoLen, MPI_INT,
                allPartitionInfo.getArray(), recvLens.getArray(),
                recvDispls.getArray(), MPI_INT, destProcessor, comm);

    if (myRank == destProcessor) {
      ij = Shiftl(numVertices, 1);

      for (j = 0; j < ij;) {
        v = allPartitionInfo[j++];
        pVector[v] = allPartitionInfo[j++];
      }

      ++myPartitionIdx;

#ifdef DEBUG_CONTROLLER
      for (j = 0; j < numVertices; ++j)
        assert(pVector[j] >= 0 && pVector[j] < numParts);
#endif
    }

    DynaMem<Bisection>::deletePtr(b);
  }

  // ###
  // k-way refine local partitions
  // ###

  h->setNumPartitions(numMyPartitions);

  if (numMyPartitions > 0) {
    for (i = 0; i < numMyPartitions; ++i) {
      pVector = &partitionVector[partitionVectorOffsets[i]];
      h->copyInPartition(pVector, numVertices, i, partitionCuts[i]);
    }

    kWayRefiner->rebalance(*h);
  }

  // ###
  // project partitions
  // ###

  initSeqPartitions(hgraph, comm);

#ifdef DEBUG_CONTROLLER
  hgraph.checkPartitions(numParts, maxPartWt, comm);
#endif

  DynaMem<Hypergraph>::deletePtr(h);
}

void RecurBisectController::initSeqPartitions(ParaHypergraph &hgraph,
                                              MPI_Comm comm) {
  register int i;
  register int j;
  register int ij;

  int numTotVertices = h->getNumVertices();
  int ijk;
  int startOffset;
  int endOffset;
  int totToSend;

  int *hGraphPartitionVector;
  int *hGraphPartVectorOffsets;
  int *hGraphPartCuts;

  int *hPartitionVector = h->getPartVectorArray();
  int *hPartOffsetsVector = h->getPartOffsetArray();
  int *hPartitionCutsArray = h->getPartCutArray();

  FastDynaArray<int> numVperProc(numProcs);
  FastDynaArray<int> procDispls(numProcs);

  FastDynaArray<int> sendLens(numProcs);
  FastDynaArray<int> sendDispls(numProcs);
  FastDynaArray<int> recvLens(numProcs);
  FastDynaArray<int> recvDispls(numProcs);
  FastDynaArray<int> sendArray;

  hgraph.setNumberPartitions(numSeqRuns);

  hGraphPartitionVector = hgraph.getPartitionArray();
  hGraphPartVectorOffsets = hgraph.getPartitionOffsetsArray();
  hGraphPartCuts = hgraph.getCutsizesArray();

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
    sendLens[i] = numVperProc[i] * numMyPartitions;
    j += sendLens[i];
    ij += numVperProc[i];
  }

  sendArray.setLength(j);
  totToSend = j;

  ij = 0;

  for (ijk = 0; ijk < numProcs; ++ijk) {
    for (j = 0; j < numMyPartitions; ++j) {
      startOffset = hPartOffsetsVector[j] + procDispls[ijk];
      endOffset = startOffset + numVperProc[ijk];

      for (i = startOffset; i < endOffset; ++i) {
        sendArray[ij++] = hPartitionVector[i];
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
  assert(ij == hGraphPartVectorOffsets[numSeqRuns]);
#endif

  MPI_Alltoallv(sendArray.getArray(), sendLens.getArray(),
                sendDispls.getArray(), MPI_INT, hGraphPartitionVector,
                recvLens.getArray(), recvDispls.getArray(), MPI_INT, comm);

  // ###
  // communicate partition cuts
  // ###

  MPI_Allgather(&numMyPartitions, 1, MPI_INT, recvLens.getArray(), 1, MPI_INT,
                comm);

  ij = 0;

  for (i = 0; i < numProcs; ++i) {
    recvDispls[i] = ij;
    ij += recvLens[i];
  }

  MPI_Allgatherv(hPartitionCutsArray, numMyPartitions, MPI_INT, hGraphPartCuts,
                 recvLens.getArray(), recvDispls.getArray(), MPI_INT, comm);

  if (dispOption > 1 && myRank == 0) {
    for (i = 0; i < numSeqRuns; ++i)
      out_stream << hGraphPartCuts[i] << " ";

    out_stream << endl;
  }
}

void RecurBisectController::recursivelyBisect(const Bisection &b,
                                              MPI_Comm comm) {
  int cut;
  int rank;
  int nProcs;
  int bisectAgain = b.getBisectAgain();

  Hypergraph *h = b.getHypergraph();

  if (h->getNumVertices() == 0)
    return;

  MPI_Comm_size(comm, &nProcs);
  MPI_Comm_rank(comm, &rank);

  if (nProcs == 1) {
    Bisection *left = NULL;
    Bisection *right = NULL;

    bisector->setNumRuns(numBisectRuns);

    if (bisectAgain == 1)
      bisector->bisect(h, maxPartWt);
    else
      bisector->bisect(h, computeMaxWt(logK - bisectAgain));

    sumOfCuts += h->keepBestPartition();

    if (bisectAgain == 1) {
      int numVerts = h->getNumVertices();
      int *partV = h->getPartVectorArray();
      int *toOrigVmap = b.getMapArray();
      int bisectionPart = b.getPartID();

      register int i;

      for (i = 0; i < numVerts; ++i) {
        locVertPartitionInfo.assign(locVertPartInfoLen++, toOrigVmap[i]);

        if (partV[i] == 0)
          locVertPartitionInfo.assign(locVertPartInfoLen++, bisectionPart);
        else
          locVertPartitionInfo.assign(locVertPartInfoLen++,
                                      (bisectionPart | (1 << (logK - 1))));
      }
    } else {
      splitBisection(b, left, right);
      recursivelyBisect(*left, comm);
      recursivelyBisect(*right, comm);

      if (left) {
        Hypergraph *hLeft = left->getHypergraph();

        DynaMem<Hypergraph>::deletePtr(hLeft);
        DynaMem<Bisection>::deletePtr(left);
      }

      if (right) {
        Hypergraph *hRight = right->getHypergraph();

        DynaMem<Hypergraph>::deletePtr(hRight);
        DynaMem<Bisection>::deletePtr(right);
      }
    }
  } else {
    Bisection *newB;
    int bestCutProc;

    bisector->setNumRuns(max(1, numBisectRuns / nProcs));

    if (bisectAgain == 1)
      bisector->bisect(h, maxPartWt);
    else
      bisector->bisect(h, computeMaxWt(logK - bisectAgain));

    cut = h->keepBestPartition();
    bestCutProc = getBestPartitionProc(cut, comm);

    if (rank == bestCutProc)
      sumOfCuts += cut;

    if (bisectAgain == 1) {
      if (rank == bestCutProc) {
        int numVerts = h->getNumVertices();
        int *partV = h->getPartVectorArray();
        int *toOrigVmap = b.getMapArray();
        int bisectionPart = b.getPartID();

        register int i;

        for (i = 0; i < numVerts; ++i) {
          locVertPartitionInfo.assign(locVertPartInfoLen++, toOrigVmap[i]);

          if (partV[i] == 0)
            locVertPartitionInfo.assign(locVertPartInfoLen++, bisectionPart);
          else
            locVertPartitionInfo.assign(locVertPartInfoLen++,
                                        (bisectionPart | (1 << (logK - 1))));
        }
      }
    } else {
      MPI_Bcast(h->getPartVectorArray(), h->getNumVertices(), MPI_INT,
                bestCutProc, comm);

      MPI_Comm new_comm;
      MPI_Comm_split(comm, And(rank, 0x1), 0, &new_comm);

      splitBisection(b, newB, comm);

      recursivelyBisect(*newB, new_comm);

      if (newB) {
        Hypergraph *hNew = newB->getHypergraph();

        DynaMem<Hypergraph>::deletePtr(hNew);
        DynaMem<Bisection>::deletePtr(newB);
      }

      MPI_Comm_free(&new_comm);
    }
  }
}

void RecurBisectController::splitBisection(const Bisection &b, Bisection *&newB,
                                           MPI_Comm comm) const {
  int rank;
  int nProcs;

  MPI_Comm_size(comm, &nProcs);
  MPI_Comm_rank(comm, &rank);

  register int i;
  register int j;

  Hypergraph *newH;
  Hypergraph *h = b.getHypergraph();

  // ###
  // h data
  // ###

  int numHVertices = h->getNumVertices();
  int numHHedges = h->getNumHedges();
  int bPartID = b.getPartID();
  int bisectAgain = b.getBisectAgain();

  int *mapToOrig = b.getMapArray();
  int *hPartVector = h->getPartVectorArray();
  int *hVertWt = h->getVerWeightsArray();
  int *hHedgeWt = h->getHedgeWeightsArray();
  int *hHedgeOffsets = h->getHedgeOffsetArray();
  int *hPinList = h->getPinListArray();

  // ###
  // newH data
  // ###

  int numVerts = 0;
  int numHedges = 0;
  int totWt = 0;
  int numPins = 0;

  FastDynaArray<int> *vertWt = new FastDynaArray<int>(64);
  FastDynaArray<int> *mapOrig = new FastDynaArray<int>(64);
  FastDynaArray<int> *hedgeWts = new FastDynaArray<int>(64);
  FastDynaArray<int> *hedgeOffsets = new FastDynaArray<int>(64);
  FastDynaArray<int> *pinList = new FastDynaArray<int>(64);

  // ###
  // auxiliary data
  // ###

  int v;
  int hEdgeLen;
  int endOffset;

  FastDynaArray<int> mapFromHtoNewH(numHVertices);

  if (And(rank, 0x1)) {
    // ###
    // in the odd processor case
    // ###

    for (i = 0; i < numHVertices; ++i) {
#ifdef DEBUG_CONTROLLER
      assert(hPartVector[i] == 1 || hPartVector[i] == 0);
#endif
      if (hPartVector[i] == 1) {
        vertWt->assign(numVerts, hVertWt[i]);
        mapOrig->assign(numVerts, mapToOrig[i]);
        mapFromHtoNewH[i] = numVerts++;
        totWt += hVertWt[i];
      } else {
        mapFromHtoNewH[i] = -1;
      }
    }

    vertWt->setLength(numVerts);
    mapOrig->setLength(numVerts);

    newH = new Hypergraph(vertWt->getArray(), numVerts);

    // ###
    // initialise pin list
    // ###

    hedgeOffsets->assign(numHedges, numPins);

    for (i = 0; i < numHHedges; ++i) {
      endOffset = hHedgeOffsets[i + 1];
      hEdgeLen = 0;

      for (j = hHedgeOffsets[i]; j < endOffset; ++j) {
        v = hPinList[j];

        if (hPartVector[v] == 1) {
#ifdef DEBUG_CONTROLLER
          assert(mapFromHtoNewH[v] != -1);
#endif
          pinList->assign(numPins + (hEdgeLen++), mapFromHtoNewH[v]);
        }
      }

      if (hEdgeLen > 1) {
        numPins += hEdgeLen;
        hedgeWts->assign(numHedges++, hHedgeWt[i]);
        hedgeOffsets->assign(numHedges, numPins);
      }
    }

    hedgeWts->setLength(numHedges);
    hedgeOffsets->setLength(numHedges + 1);
    pinList->setLength(numPins);

    // ###
    // now init the hypergraphs
    // ###

    newH->setNumHedges(numHedges);
    newH->setNumPins(numPins);
    newH->setTotWeight(totWt);
    newH->setHedgeWtArray(hedgeWts->getArray(), numHedges);
    newH->setPinListArray(pinList->getArray(), numPins);
    newH->setHedgeOffsetArray(hedgeOffsets->getArray(), numHedges + 1);

    newH->buildVtoHedges();

    // ###
    // now init the bisections
    // ###

    newB = new Bisection(newH, bisectAgain - 1,
                         Or(bPartID, Shiftl(1, (logK - bisectAgain))));
    newB->setMap(mapOrig->getArray(), numVerts);
  } else {
    // ###
    // in the even processor case
    // ###

    for (i = 0; i < numHVertices; ++i) {
#ifdef DEBUG_CONTROLLER
      assert(hPartVector[i] == 1 || hPartVector[i] == 0);
#endif
      if (hPartVector[i] == 0) {
        vertWt->assign(numVerts, hVertWt[i]);
        mapOrig->assign(numVerts, mapToOrig[i]);
        mapFromHtoNewH[i] = numVerts++;
        totWt += hVertWt[i];
      } else {
        mapFromHtoNewH[i] = -1;
      }
    }

    vertWt->setLength(numVerts);
    mapOrig->setLength(numVerts);

    newH = new Hypergraph(vertWt->getArray(), numVerts);

    // ###
    // initialise pin list
    // ###

    hedgeOffsets->assign(numHedges, numPins);

    for (i = 0; i < numHHedges; ++i) {
      endOffset = hHedgeOffsets[i + 1];
      hEdgeLen = 0;

      for (j = hHedgeOffsets[i]; j < endOffset; ++j) {
        v = hPinList[j];

        if (hPartVector[v] == 0) {
#ifdef DEBUG_CONTROLLER
          assert(mapFromHtoNewH[v] != -1);
#endif
          pinList->assign(numPins + (hEdgeLen++), mapFromHtoNewH[v]);
        }
      }

      if (hEdgeLen > 1) {
        numPins += hEdgeLen;
        hedgeWts->assign(numHedges++, hHedgeWt[i]);
        hedgeOffsets->assign(numHedges, numPins);
      }
    }

    hedgeWts->setLength(numHedges);
    hedgeOffsets->setLength(numHedges + 1);
    pinList->setLength(numPins);

    // ###
    // now init the hypergraph
    // ###

    newH->setNumHedges(numHedges);
    newH->setNumPins(numPins);
    newH->setTotWeight(totWt);
    newH->setHedgeWtArray(hedgeWts->getArray(), numHedges);
    newH->setPinListArray(pinList->getArray(), numPins);
    newH->setHedgeOffsetArray(hedgeOffsets->getArray(), numHedges + 1);

    newH->buildVtoHedges();

    // ###
    // now init the bisections
    // ###

    newB = new Bisection(newH, bisectAgain - 1, bPartID);
    newB->setMap(mapOrig->getArray(), numVerts);
  }
}

void RecurBisectController::splitBisection(const Bisection &b, Bisection *&l,
                                           Bisection *&r) const {
  register int i;
  register int j;

  Hypergraph *leftH;
  Hypergraph *rightH;
  Hypergraph *h = b.getHypergraph();

  // ###
  // h data
  // ###

  int numHVertices = h->getNumVertices();
  int numHHedges = h->getNumHedges();
  int bPartID = b.getPartID();
  int bisectAgain = b.getBisectAgain();

  int *mapToOrig = b.getMapArray();
  int *hPartVector = h->getPartVectorArray();
  int *hVertWt = h->getVerWeightsArray();
  int *hHedgeWt = h->getHedgeWeightsArray();
  int *hHedgeOffsets = h->getHedgeOffsetArray();
  int *hPinList = h->getPinListArray();

  // ###
  // leftH data
  // ###

  int numLeftVerts = 0;
  int numLeftHedges = 0;
  int totLeftWt = 0;
  int numLeftPins = 0;

  FastDynaArray<int> *leftVertWt = new FastDynaArray<int>(64);
  FastDynaArray<int> *leftMapOrig = new FastDynaArray<int>(64);
  FastDynaArray<int> *leftHedgeWts = new FastDynaArray<int>(64);
  FastDynaArray<int> *leftHedgeOffsets = new FastDynaArray<int>(64);
  FastDynaArray<int> *leftPinList = new FastDynaArray<int>(64);

  // ###
  // rightH data
  // ###

  int numRightVerts = 0;
  int numRightHedges = 0;
  int totRightWt = 0;
  int numRightPins = 0;

  FastDynaArray<int> *rightVertWt = new FastDynaArray<int>(64);
  FastDynaArray<int> *rightMapOrig = new FastDynaArray<int>(64);
  FastDynaArray<int> *rightHedgeWts = new FastDynaArray<int>(64);
  FastDynaArray<int> *rightHedgeOffsets = new FastDynaArray<int>(64);
  FastDynaArray<int> *rightPinList = new FastDynaArray<int>(64);

  // ###
  // auxiliary data
  // ###

  int v;
  int endOffset;
  int leftHedgeLen;
  int rightHedgeLen;

  FastDynaArray<int> mapFromHtoNewH(numHVertices);

  for (i = 0; i < numHVertices; ++i) {
#ifdef DEBUG_CONTROLLER
    assert(hPartVector[i] == 1 || hPartVector[i] == 0);
#endif
    if (hPartVector[i] == 0) {
      leftVertWt->assign(numLeftVerts, hVertWt[i]);
      leftMapOrig->assign(numLeftVerts, mapToOrig[i]);
      mapFromHtoNewH[i] = numLeftVerts++;
      totLeftWt += hVertWt[i];
    } else {
      rightVertWt->assign(numRightVerts, hVertWt[i]);
      rightMapOrig->assign(numRightVerts, mapToOrig[i]);
      mapFromHtoNewH[i] = numRightVerts++;
      totRightWt += hVertWt[i];
    }
  }

  leftVertWt->setLength(numLeftVerts);
  leftMapOrig->setLength(numLeftVerts);

  rightVertWt->setLength(numRightVerts);
  rightMapOrig->setLength(numRightVerts);

  leftH = new Hypergraph(leftVertWt->getArray(), numLeftVerts);
  rightH = new Hypergraph(rightVertWt->getArray(), numRightVerts);

  // ###
  // initialise pin list
  // ###

  leftHedgeOffsets->assign(numLeftHedges, numLeftPins);
  rightHedgeOffsets->assign(numRightHedges, numRightPins);

  for (i = 0; i < numHHedges; ++i) {
    endOffset = hHedgeOffsets[i + 1];

    leftHedgeLen = 0;
    rightHedgeLen = 0;

    for (j = hHedgeOffsets[i]; j < endOffset; ++j) {
      v = hPinList[j];

      if (hPartVector[v] == 0) {
        leftPinList->assign(numLeftPins + (leftHedgeLen++), mapFromHtoNewH[v]);
      } else {
        rightPinList->assign(numRightPins + (rightHedgeLen++),
                             mapFromHtoNewH[v]);
      }
    }

    if (leftHedgeLen > 1) {
      numLeftPins += leftHedgeLen;
      leftHedgeWts->assign(numLeftHedges++, hHedgeWt[i]);
      leftHedgeOffsets->assign(numLeftHedges, numLeftPins);
    }

    if (rightHedgeLen > 1) {
      numRightPins += rightHedgeLen;
      rightHedgeWts->assign(numRightHedges++, hHedgeWt[i]);
      rightHedgeOffsets->assign(numRightHedges, numRightPins);
    }
  }

  leftHedgeWts->setLength(numLeftHedges);
  leftHedgeOffsets->setLength(numLeftHedges + 1);
  leftPinList->setLength(numLeftPins);

  rightHedgeWts->setLength(numRightHedges);
  rightHedgeOffsets->setLength(numRightHedges + 1);
  rightPinList->setLength(numRightPins);

  // ###
  // now init the hypergraphs
  // ###

  leftH->setNumHedges(numLeftHedges);
  leftH->setNumPins(numLeftPins);
  leftH->setTotWeight(totLeftWt);
  leftH->setHedgeWtArray(leftHedgeWts->getArray(), numLeftHedges);
  leftH->setPinListArray(leftPinList->getArray(), numLeftPins);
  leftH->setHedgeOffsetArray(leftHedgeOffsets->getArray(), numLeftHedges + 1);

  leftH->buildVtoHedges();

  rightH->setNumHedges(numRightHedges);
  rightH->setNumPins(numRightPins);
  rightH->setTotWeight(totRightWt);
  rightH->setHedgeWtArray(rightHedgeWts->getArray(), numRightHedges);
  rightH->setPinListArray(rightPinList->getArray(), numRightPins);
  rightH->setHedgeOffsetArray(rightHedgeOffsets->getArray(),
                              numRightHedges + 1);

  rightH->buildVtoHedges();

  // ###
  // now init the bisections
  // ###

  l = new Bisection(leftH, bisectAgain - 1, bPartID);
  r = new Bisection(rightH, bisectAgain - 1,
                    Or(bPartID, Shiftl(1, (logK - bisectAgain))));

  l->setMap(leftMapOrig->getArray(), numLeftVerts);
  r->setMap(rightMapOrig->getArray(), numRightVerts);
}

int RecurBisectController::getBestPartitionProc(int cut, MPI_Comm comm) const {
  int rank;
  int nProcs;
  int bestCut;
  int bestProc;

  register int i;

  MPI_Comm_size(comm, &nProcs);
  MPI_Comm_rank(comm, &rank);

  FastDynaArray<int> allCuts(nProcs);
  FastDynaArray<int> procs(nProcs);

  MPI_Allgather(&cut, 1, MPI_INT, allCuts.getArray(), 1, MPI_INT, comm);

  bestCut = allCuts[0];
  procs[0] = 0;

  for (i = 1; i < nProcs; ++i) {
    procs[i] = i;

    if (allCuts[i] < bestCut)
      bestCut = allCuts[i];
  }

  Funct::randomPermutation(procs.getArray(), nProcs);

  if (rank == 0) {
    for (i = 0; i < nProcs; ++i) {
      if (allCuts[procs[i]] == bestCut) {
        bestProc = procs[i];
        break;
      }
    }
  }

  MPI_Bcast(&bestProc, 1, MPI_INT, 0, comm);

  return bestProc;
}

int RecurBisectController::computeMaxWt(int numBs) const {
  double maxPartWeight = recursivelyComputeMax(aveInitBisectionWt, numBs);

  return (static_cast<int>(floor(maxPartWeight)));
}

double RecurBisectController::recursivelyComputeMax(double currAve,
                                                    int depth) const {
  if (depth == 0)
    return (currAve + currAve * bisectConstraint);
  return (recursivelyComputeMax((currAve + currAve * bisectConstraint) / 2,
                                depth - 1));
}

#endif
