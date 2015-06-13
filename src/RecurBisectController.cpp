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
#include "hypergraph/parallel/hypergraph.hpp"
#include "hypergraph/serial/hypergraph.hpp"

namespace parallel = parkway::hypergraph::parallel;
namespace serial = parkway::hypergraph::serial;

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

  locVertPartitionInfo.reserve(0);
  allPartitionInfo.reserve(0);

#ifdef DEBUG_CONTROLLER
  assert(bisector);
#endif
}

RecurBisectController::~RecurBisectController() {
  DynaMem::deletePtr<BisectionController>(bisector);
  DynaMem::deletePtr<GreedyKwayRefiner>(kWayRefiner);
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

  int j;
  int i;

  double a = (1.0 + kWayConstraint) / numParts;
  double b = 1.0 / logK;

  bisectConstraint = pow(a, b) - 0.5;

  j = h->total_weight();

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

  partitionCuts.reserve(numMyPartitions);
  partitionVectorOffsets.reserve(numMyPartitions + 1);
  partitionVectorOffsets[0] = 0;

  j = h->number_of_vertices();

  for (i = 1; i <= numMyPartitions; ++i) {
    partitionVectorOffsets[i] = partitionVectorOffsets[i - 1] + j;
  }

  partitionVector.reserve(partitionVectorOffsets[numMyPartitions]);
}

void RecurBisectController::runSeqPartitioner(parallel::hypergraph &hgraph,
                                              MPI_Comm comm) {
  initCoarsestHypergraph(hgraph, comm);
  convToBisectionConstraints();

  if (dispOption > 1 && myRank == 0) {
    out_stream << "[R-B]: " << numSeqRuns << " | ";
  }

  int i;
  int j;
  int ij;

  int numVertices = h->number_of_vertices();
  int *pVector = nullptr;
  int destProcessor;
  int myPartitionIdx = 0;
  int v;

  dynamic_array<int> recvLens(numProcs);
  dynamic_array<int> recvDispls(numProcs);

  Bisection *b;

  allPartitionInfo.reserve(Shiftl(numVertices, 1));

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
    MPI_Gather(&locVertPartInfoLen, 1, MPI_INT, recvLens.data(), 1, MPI_INT,
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

    MPI_Gatherv(locVertPartitionInfo.data(), locVertPartInfoLen, MPI_INT,
                allPartitionInfo.data(), recvLens.data(),
                recvDispls.data(), MPI_INT, destProcessor, comm);

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

    DynaMem::deletePtr<Bisection>(b);
  }

  // ###
  // k-way refine local partitions
  // ###

  h->set_number_of_partitions(numMyPartitions);

  if (numMyPartitions > 0) {
    for (i = 0; i < numMyPartitions; ++i) {
      pVector = &partitionVector[partitionVectorOffsets[i]];
      h->copy_in_partition(pVector, numVertices, i, partitionCuts[i]);
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

  DynaMem::deletePtr<serial::hypergraph>(h);
}

void RecurBisectController::initSeqPartitions(parallel::hypergraph &hgraph,
                                              MPI_Comm comm) {
  int i;
  int j;
  int ij;

  int numTotVertices = h->number_of_vertices();
  int ijk;
  int startOffset;
  int endOffset;
  int totToSend;

  int *hGraphPartitionVector;
  int *hGraphPartVectorOffsets;
  int *hGraphPartCuts;

  int *hPartitionVector = h->partition_vector();
  int *hPartOffsetsVector = h->partition_offsets();
  int *hPartitionCutsArray = h->partition_cuts();

  dynamic_array<int> numVperProc(numProcs);
  dynamic_array<int> procDispls(numProcs);

  dynamic_array<int> sendLens(numProcs);
  dynamic_array<int> sendDispls(numProcs);
  dynamic_array<int> recvLens(numProcs);
  dynamic_array<int> recvDispls(numProcs);
  dynamic_array<int> sendArray;

  hgraph.set_number_of_partitions(numSeqRuns);

  hGraphPartitionVector = hgraph.partition_vector();
  hGraphPartVectorOffsets = hgraph.partition_offsets();
  hGraphPartCuts = hgraph.partition_cuts();

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

  sendArray.reserve(j);
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

  MPI_Alltoall(sendLens.data(), 1, MPI_INT, recvLens.data(), 1, MPI_INT,
               comm);

  ij = 0;

  for (i = 0; i < numProcs; ++i) {
    recvDispls[i] = ij;
    ij += recvLens[i];
  }

#ifdef DEBUG_CONTROLLER
  assert(ij == hGraphPartVectorOffsets[numSeqRuns]);
#endif

  MPI_Alltoallv(sendArray.data(), sendLens.data(),
                sendDispls.data(), MPI_INT, hGraphPartitionVector,
                recvLens.data(), recvDispls.data(), MPI_INT, comm);

  // ###
  // communicate partition cuts
  // ###

  MPI_Allgather(&numMyPartitions, 1, MPI_INT, recvLens.data(), 1, MPI_INT,
                comm);

  ij = 0;

  for (i = 0; i < numProcs; ++i) {
    recvDispls[i] = ij;
    ij += recvLens[i];
  }

  MPI_Allgatherv(hPartitionCutsArray, numMyPartitions, MPI_INT, hGraphPartCuts,
                 recvLens.data(), recvDispls.data(), MPI_INT, comm);

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

  serial::hypergraph *h = b.getHypergraph();

  if (h->number_of_vertices() == 0)
    return;

  MPI_Comm_size(comm, &nProcs);
  MPI_Comm_rank(comm, &rank);

  if (nProcs == 1) {
    Bisection *left = nullptr;
    Bisection *right = nullptr;

    bisector->setNumRuns(numBisectRuns);

    if (bisectAgain == 1)
      bisector->bisect(h, maxPartWt);
    else
      bisector->bisect(h, computeMaxWt(logK - bisectAgain));

    sumOfCuts += h->keep_best_partition();

    if (bisectAgain == 1) {
      int numVerts = h->number_of_vertices();
      int *partV = h->partition_vector();
      int *toOrigVmap = b.getMapArray();
      int bisectionPart = b.getPartID();

      int i;

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
        serial::hypergraph *hLeft = left->getHypergraph();

        DynaMem::deletePtr<serial::hypergraph>(hLeft);
        DynaMem::deletePtr<Bisection>(left);
      }

      if (right) {
        serial::hypergraph *hRight = right->getHypergraph();

        DynaMem::deletePtr<serial::hypergraph>(hRight);
        DynaMem::deletePtr<Bisection>(right);
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

    cut = h->keep_best_partition();
    bestCutProc = getBestPartitionProc(cut, comm);

    if (rank == bestCutProc)
      sumOfCuts += cut;

    if (bisectAgain == 1) {
      if (rank == bestCutProc) {
        int numVerts = h->number_of_vertices();
        int *partV = h->partition_vector();
        int *toOrigVmap = b.getMapArray();
        int bisectionPart = b.getPartID();

        int i;

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
      MPI_Bcast(h->partition_vector(), h->number_of_vertices(), MPI_INT,
                bestCutProc, comm);

      MPI_Comm new_comm;
      MPI_Comm_split(comm, And(rank, 0x1), 0, &new_comm);

      splitBisection(b, newB, comm);

      recursivelyBisect(*newB, new_comm);

      if (newB) {
        serial::hypergraph *hNew = newB->getHypergraph();

        DynaMem::deletePtr<serial::hypergraph>(hNew);
        DynaMem::deletePtr<Bisection>(newB);
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

  int i;
  int j;

  serial::hypergraph *newH;
  serial::hypergraph *h = b.getHypergraph();

  // ###
  // h data_
  // ###

  int numHVertices = h->number_of_vertices();
  int numHHedges = h->number_of_hyperedges();
  int bPartID = b.getPartID();
  int bisectAgain = b.getBisectAgain();

  int *mapToOrig = b.getMapArray();
  int *hPartVector = h->partition_vector();
  int *hVertWt = h->vertex_weights();
  int *hHedgeWt = h->hyperedge_weights();
  int *hHedgeOffsets = h->hyperedge_offsets();
  int *hPinList = h->pin_list();

  // ###
  // newH data_
  // ###

  int numVerts = 0;
  int numHedges = 0;
  int totWt = 0;
  int numPins = 0;

  dynamic_array<int> *vertWt = new dynamic_array<int>(64);
  dynamic_array<int> *mapOrig = new dynamic_array<int>(64);
  dynamic_array<int> *hedgeWts = new dynamic_array<int>(64);
  dynamic_array<int> *hedgeOffsets = new dynamic_array<int>(64);
  dynamic_array<int> *pinList = new dynamic_array<int>(64);

  // ###
  // auxiliary data_
  // ###

  int v;
  int hEdgeLen;
  int endOffset;

  dynamic_array<int> mapFromHtoNewH(numHVertices);

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

    vertWt->reserve(numVerts);
    mapOrig->reserve(numVerts);

    newH = new serial::hypergraph(vertWt->data(), numVerts);

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

    hedgeWts->reserve(numHedges);
    hedgeOffsets->reserve(numHedges + 1);
    pinList->reserve(numPins);

    // ###
    // now init the hypergraphs
    // ###

    newH->set_number_of_hyperedges(numHedges);
    newH->set_number_of_pins(numPins);
    newH->set_total_weight(totWt);
    newH->set_hyperedge_weights(hedgeWts->data(), numHedges);
    newH->set_pin_list(pinList->data(), numPins);
    newH->set_hyperedge_offsets(hedgeOffsets->data(), numHedges + 1);

    newH->buildVtoHedges();

    // ###
    // now init the bisections
    // ###

    newB = new Bisection(newH, bisectAgain - 1,
                         Or(bPartID, Shiftl(1, (logK - bisectAgain))));
    newB->setMap(mapOrig->data(), numVerts);
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

    vertWt->reserve(numVerts);
    mapOrig->reserve(numVerts);

    newH = new serial::hypergraph(vertWt->data(), numVerts);

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

    hedgeWts->reserve(numHedges);
    hedgeOffsets->reserve(numHedges + 1);
    pinList->reserve(numPins);

    // ###
    // now init the hypergraph
    // ###

    newH->set_number_of_hyperedges(numHedges);
    newH->set_number_of_pins(numPins);
    newH->set_total_weight(totWt);
    newH->set_hyperedge_weights(hedgeWts->data(), numHedges);
    newH->set_pin_list(pinList->data(), numPins);
    newH->set_hyperedge_offsets(hedgeOffsets->data(), numHedges + 1);

    newH->buildVtoHedges();

    // ###
    // now init the bisections
    // ###

    newB = new Bisection(newH, bisectAgain - 1, bPartID);
    newB->setMap(mapOrig->data(), numVerts);
  }
}

void RecurBisectController::splitBisection(const Bisection &b, Bisection *&l,
                                           Bisection *&r) const {
  int i;
  int j;

  serial::hypergraph *leftH;
  serial::hypergraph *rightH;
  serial::hypergraph *h = b.getHypergraph();

  // ###
  // h data_
  // ###

  int numHVertices = h->number_of_vertices();
  int numHHedges = h->number_of_hyperedges();
  int bPartID = b.getPartID();
  int bisectAgain = b.getBisectAgain();

  int *mapToOrig = b.getMapArray();
  int *hPartVector = h->partition_vector();
  int *hVertWt = h->vertex_weights();
  int *hHedgeWt = h->hyperedge_weights();
  int *hHedgeOffsets = h->hyperedge_offsets();
  int *hPinList = h->pin_list();

  // ###
  // leftH data_
  // ###

  int numLeftVerts = 0;
  int numLeftHedges = 0;
  int totLeftWt = 0;
  int numLeftPins = 0;

  dynamic_array<int> *leftVertWt = new dynamic_array<int>(64);
  dynamic_array<int> *leftMapOrig = new dynamic_array<int>(64);
  dynamic_array<int> *leftHedgeWts = new dynamic_array<int>(64);
  dynamic_array<int> *leftHedgeOffsets = new dynamic_array<int>(64);
  dynamic_array<int> *leftPinList = new dynamic_array<int>(64);

  // ###
  // rightH data_
  // ###

  int numRightVerts = 0;
  int numRightHedges = 0;
  int totRightWt = 0;
  int numRightPins = 0;

  dynamic_array<int> *rightVertWt = new dynamic_array<int>(64);
  dynamic_array<int> *rightMapOrig = new dynamic_array<int>(64);
  dynamic_array<int> *rightHedgeWts = new dynamic_array<int>(64);
  dynamic_array<int> *rightHedgeOffsets = new dynamic_array<int>(64);
  dynamic_array<int> *rightPinList = new dynamic_array<int>(64);

  // ###
  // auxiliary data_
  // ###

  int v;
  int endOffset;
  int leftHedgeLen;
  int rightHedgeLen;

  dynamic_array<int> mapFromHtoNewH(numHVertices);

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

  leftVertWt->reserve(numLeftVerts);
  leftMapOrig->reserve(numLeftVerts);

  rightVertWt->reserve(numRightVerts);
  rightMapOrig->reserve(numRightVerts);

  leftH = new serial::hypergraph(leftVertWt->data(), numLeftVerts);
  rightH = new serial::hypergraph(rightVertWt->data(), numRightVerts);

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

  leftHedgeWts->reserve(numLeftHedges);
  leftHedgeOffsets->reserve(numLeftHedges + 1);
  leftPinList->reserve(numLeftPins);

  rightHedgeWts->reserve(numRightHedges);
  rightHedgeOffsets->reserve(numRightHedges + 1);
  rightPinList->reserve(numRightPins);

  // ###
  // now init the hypergraphs
  // ###

  leftH->set_number_of_hyperedges(numLeftHedges);
  leftH->set_number_of_pins(numLeftPins);
  leftH->set_total_weight(totLeftWt);
  leftH->set_hyperedge_weights(leftHedgeWts->data(), numLeftHedges);
  leftH->set_pin_list(leftPinList->data(), numLeftPins);
  leftH->set_hyperedge_offsets(leftHedgeOffsets->data(), numLeftHedges + 1);

  leftH->buildVtoHedges();

  rightH->set_number_of_hyperedges(numRightHedges);
  rightH->set_number_of_pins(numRightPins);
  rightH->set_total_weight(totRightWt);
  rightH->set_hyperedge_weights(rightHedgeWts->data(), numRightHedges);
  rightH->set_pin_list(rightPinList->data(), numRightPins);
  rightH->set_hyperedge_offsets(rightHedgeOffsets->data(),
                                numRightHedges + 1);

  rightH->buildVtoHedges();

  // ###
  // now init the bisections
  // ###

  l = new Bisection(leftH, bisectAgain - 1, bPartID);
  r = new Bisection(rightH, bisectAgain - 1,
                    Or(bPartID, Shiftl(1, (logK - bisectAgain))));

  l->setMap(leftMapOrig->data(), numLeftVerts);
  r->setMap(rightMapOrig->data(), numRightVerts);
}

int RecurBisectController::getBestPartitionProc(int cut, MPI_Comm comm) const {
  int rank;
  int nProcs;
  int bestCut;
  int bestProc;

  int i;

  MPI_Comm_size(comm, &nProcs);
  MPI_Comm_rank(comm, &rank);

  dynamic_array<int> allCuts(nProcs);
  dynamic_array<int> procs(nProcs);

  MPI_Allgather(&cut, 1, MPI_INT, allCuts.data(), 1, MPI_INT, comm);

  bestCut = allCuts[0];
  procs[0] = 0;

  for (i = 1; i < nProcs; ++i) {
    procs[i] = i;

    if (allCuts[i] < bestCut)
      bestCut = allCuts[i];
  }

  Funct::randomPermutation(procs.data(), nProcs);

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
