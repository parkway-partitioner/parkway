
#ifndef _KHMETIS_CONTROLLER_CPP
#define _KHMETIS_CONTROLLER_CPP

// ### KHMetisController.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 30/11/2004: Last Modified
//
// ###

#include "Config.h"
#ifdef LINK_HMETIS
#include "KHMetisController.hpp"

KHMetisController::KHMetisController(GreedyKwayRefiner *k, int rank, int nProcs,
                                     int nParts, const int *options,
                                     ostream &out)
    : SeqController(rank, nProcs, nParts, out) {
  int i;

  kWayRefiner = k;
  lenOfOptions = 9;
  maxPartWt = 0;
  avePartWt = 0;

  khMetisOptions.setLength(9);

  for (i = 0; i < lenOfOptions; ++i)
    khMetisOptions[i] = options[i];
}

KHMetisController::~KHMetisController() {
  DynaMem<GreedyKwayRefiner>::deletePtr(kWayRefiner);
}

void KHMetisController::dispSeqControllerOptions() const {
  switch (dispOption) {
  case SILENT:
    break;

  default:

    out_stream << "|--- SEQ_CON:" << endl
               << "|- KHMETIS:"
               << " c = " << khMetisOptions[2] << " r = " << khMetisOptions[4]
               << endl
               << "|" << endl;
    break;
  }
}

void KHMetisController::runSeqPartitioner(ParaHypergraph &hgraph,
                                          MPI_Comm comm) {
  initCoarsestHypergraph(hgraph, comm);

#ifdef DEBUG_CONTROLLER
  assert(h);
#endif

  h->setNumPartitions(1);

  khMetisOptions[1] = numSeqRuns / numProcs + 1;

  if (dispOption > 1 && myRank == 0) {
    out_stream << "[KHMETIS]: " << numSeqRuns << " " << khMetisOptions[1]
               << " | ";
  }

  int numVertices = h->getNumVertices();
  int numHedges = h->getNumHedges();
  int totWt = h->getTotWeight();

  int numHedgesCut;
  int ubFactor = static_cast<int>(floor(kWayConstraint * 100));

  int *vWeights = h->getVerWeightsArray();
  int *hOffsets = h->getHedgeOffsetArray();
  int *pinList = h->getPinListArray();
  int *hEdgeWts = h->getHedgeWeightsArray();
  int *pArray = h->getPartVectorArray();

  khMetisOptions[7] = RANDOM(1, 10000000);

  HMETIS_PartKway(numVertices, numHedges, vWeights, hOffsets, pinList, hEdgeWts,
                  numParts, ubFactor, khMetisOptions.getArray(), pArray,
                  &numHedgesCut);

  avePartWt = static_cast<double>(totWt) / numParts;
  maxPartWt = static_cast<int>(floor(avePartWt + avePartWt * kWayConstraint));

  h->initCutsizes(numParts);

  kWayRefiner->setMaxPartWt(maxPartWt);
  kWayRefiner->setAvePartWt(avePartWt);
  kWayRefiner->rebalance(*h);

#ifdef DEBUG_CONTROLLER
  MPI_Barrier(comm);
  for (int i = 0; i < numProcs; ++i) {
    if (myRank == i)
      h->checkPartitions(numParts, maxPartWt);

    MPI_Barrier(comm);
  }
#endif

  // ###
  // project partitions
  // ###

  initSeqPartitions(hgraph, comm);

#ifdef DEBUG_CONTROLLER
  hgraph.checkPartitions(numParts, maxPartWt, comm);
#endif

  DynaMem<Hypergraph>::deletePtr(h);
}

#endif
#endif
