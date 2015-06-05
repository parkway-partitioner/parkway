
#ifndef _PATOH_CONTROLLER_CPP
#define _PATOH_CONTROLLER_CPP

// ### PaToHController.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 25/1/2005: Last Modified
//
// ###

#include "Config.h"
#ifdef LINK_PATOH
#include "PaToHController.hpp"

PaToHController::PaToHController(GreedyKwayRefiner *k, int rank, int nProcs,
                                 int nParts, int setting, ostream &out)
    : SeqController(rank, nProcs, nParts, out) {
  kWayRefiner = k;
  maxPartWt = 0;
  avePartWt = 0;
  bisectConstraint = 0;

  switch (setting) {
  case 1:
    patohSettings = PATOH_SUGPARAM_DEFAULT;
    break;

  case 2:
    patohSettings = PATOH_SUGPARAM_SPEED;
    break;

  case 3:
    patohSettings = PATOH_SUGPARAM_QUALITY;
    break;

  default:
    patohSettings = PATOH_SUGPARAM_DEFAULT;
    break;
  }
}

PaToHController::~PaToHController() {
  DynaMem<GreedyKwayRefiner>::deletePtr(kWayRefiner);
}

void PaToHController::dispSeqControllerOptions() const {
  switch (dispOption) {
  case SILENT:
    break;

  default:

    out_stream << "|--- SEQ_CON:" << endl << "|- PaToH:";

    switch (patohSettings) {
    case PATOH_SUGPARAM_DEFAULT:
      out_stream << " setting: DEFAULT ";
      break;
    case PATOH_SUGPARAM_SPEED:
      out_stream << " setting: SPEED ";
      break;
    case PATOH_SUGPARAM_QUALITY:
      out_stream << " setting: QUALITY ";
      break;
    default:
      out_stream << " setting: DEFAULT ";
      break;
    }

    out_stream << endl << "|" << endl;
    break;
  }
}

void PaToHController::runSeqPartitioner(ParaHypergraph &hgraph, MPI_Comm comm) {
  initCoarsestHypergraph(hgraph, comm);

#ifdef DEBUG_CONTROLLER
  assert(h);
#endif

  int j;

  h->setNumPartitions(1);

  /*
  if(myRank == 0)
    h->printPercentiles(out_stream);
  */

  int numVertices = h->getNumVertices();
  int numHedges = h->getNumHedges();
  int totWt = h->getTotWeight();
  int numLocalRuns = numSeqRuns / numProcs + 1;
  int _nconst = 1;
  int bestCut = LARGE_CONSTANT;
  int cut;
  int i;

  int *vWeights = h->getVerWeightsArray();
  int *hOffsets = h->getHedgeOffsetArray();
  int *pinList = h->getPinListArray();
  int *hEdgeWts = h->getHedgeWeightsArray();
  int *pArray = h->getPartVectorArray();

  FastDynaArray<int> bestPartition(numVertices);
  FastDynaArray<int> pWeights(numParts);

  bisectConstraint = Funct::toRecurBal(kWayConstraint, numParts);
#ifdef DEBUG_CONTROLLER
  assert(bisectConstraint > 0 && bisectConstraint < 1);
#endif

  if (dispOption > 1 && myRank == 0) {
    out_stream << "[PaToH]: " << numSeqRuns << " " << bisectConstraint << " | ";
  }

  PaToH_Arguments args;
  PaToH_Initialize_Parameters(&args, PATOH_CONPART, patohSettings);

  args._k = numParts;
  args.seed = RANDOM(1, 10000000);

  args.init_imbal = bisectConstraint;
  args.final_imbal = kWayConstraint;

  args.outputdetail = PATOH_OD_LOW;
  args.cuttype = PATOH_CONPART;

  args.crs_alg = 12; // PATOH_CRS_HCC;
  args.ref_alg = PATOH_REFALG_BFMKL;
  args.bigVcycle = 0;
  args.smallVcycle = 0;
  // args.usesamematchinginVcycles = 1;

  PaToH_Alloc(&args, numVertices, numHedges, _nconst, vWeights, hEdgeWts,
              hOffsets, pinList);

  if (myRank == 0) {
    out_stream << "cut type = " << args.cuttype << endl
               << "crs_alg = " << args.crs_alg << endl
               << "bigVcycle = " << args.bigVcycle << endl
               << "smallVcycle = " << args.smallVcycle << endl
               << "outputdetail = " << args.outputdetail << endl
               << "ref alg = " << args.ref_alg << endl;
  }

  MPI_Barrier(comm);

  double startTime = MPI_Wtime();

  for (i = 0; i < numLocalRuns; ++i) {
    PaToH_Partition(&args, numVertices, numHedges, vWeights, hEdgeWts, hOffsets,
                    pinList, pArray, pWeights.getArray(), &cut);

    if (cut < bestCut) {
      bestCut = cut;

      for (j = 0; j < numVertices; ++j)
        bestPartition[j] = pArray[j];
    }
  }

  PaToH_Free();

  double endTime = MPI_Wtime();
  double totTime = endTime - startTime;
  double overallLongTime;

  MPI_Allreduce(&totTime, &overallLongTime, 1, MPI_DOUBLE, MPI_MAX, comm);

  if (myRank == 0)
    out_stream << "time = " << overallLongTime << endl;

  h->copyInPartition(bestPartition.getArray(), numVertices, 0, bestCut);
  h->initCutsizes(numParts);

  out_stream << "p[" << myRank << "] after pth = " << h->getCut(0) << endl;

  avePartWt = static_cast<double>(totWt) / numParts;
  maxPartWt = static_cast<int>(floor(avePartWt + avePartWt * kWayConstraint));

  kWayRefiner->setMaxPartWt(maxPartWt);
  kWayRefiner->setAvePartWt(avePartWt);
  kWayRefiner->rebalance(*h);

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
