
#ifndef _BISECTION_CONTROLLER_CPP
#define _BISECTION_CONTROLLER_CPP

// ### BisectionController.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 4/1/2005: Last Modified
//
// ###

#include "BisectionController.hpp"

BisectionController::BisectionController(int nRuns, double kT, double redFactor,
                                         int eeP, int percentile, int inc,
                                         int dispL, ostream &out)
    : out_stream(out) {
  numSeqRuns = nRuns;
  eeParam = eeP;
  dispLevel = dispL;
  startPercentile = percentile;
  percentileIncrement = inc;
  keepThreshold = kT;
  reductionFactor = redFactor;
  numOrigVertices = 0;
  maxPartWt = 0;

  coarsener = NULL;
  refiner = NULL;
  initBisector = NULL;
}

BisectionController::~BisectionController() {
  DynaMem<Coarsener>::deletePtr(coarsener);
  DynaMem<Refiner>::deletePtr(refiner);
}

void BisectionController::dispBisectionControllerOptions() const {
  switch (dispLevel) {
  case SILENT:
    break;

  default:

    out_stream << "|- BSECTOR:";
#ifdef DEBUG_CONTROLLER
    assert(coarsener && initBisector && refiner);
#endif
    out_stream << " kT = " << keepThreshold << " rF = " << reductionFactor
               << " %le = " << startPercentile
               << " %inc = " << percentileIncrement << endl
               << "|" << endl;
    coarsener->dispCoarsenerOptions(out_stream);
    initBisector->dispInitBisectorOptions(out_stream);
    refiner->dispRefinerOptions(out_stream);

    break;
  }
}

void BisectionController::buildCoarsener(double redRatio, int cType,
                                         int minNodes) {
  switch (cType) {
  case FCwithFanOutDiv:
    coarsener =
        new FCCoarsener(Shiftl(minNodes, 1), -1, redRatio, 1, 1, dispLevel);
    break;

  case FCwithoutFanOutDiv:
    coarsener =
        new FCCoarsener(Shiftl(minNodes, 1), -1, redRatio, 0, 1, dispLevel);
    break;

  case FCwithFanOut:
    coarsener =
        new FCCoarsener(Shiftl(minNodes, 1), -1, redRatio, 1, 0, dispLevel);
    break;

  case FCwithoutFanOut:
    coarsener =
        new FCCoarsener(Shiftl(minNodes, 1), -1, redRatio, 0, 0, dispLevel);
    break;

  default:
    coarsener =
        new FCCoarsener(Shiftl(minNodes, 1), -1, redRatio, 1, 1, dispLevel);
    break;
  }
}

void BisectionController::buildInitBisector(int numInitRuns) {
  initBisector = new InitBisector(numInitRuns, FIFO, eeParam, dispLevel);
}

void BisectionController::buildRefiner(int queueD) {
  refiner = new FMRefiner(-1, queueD, eeParam, dispLevel);
}

void BisectionController::computeBisection() {
  Hypergraph *origGraph = hGraphs.getTopElem();

  int i;
  int cut;
  int bestCut = LARGE_CONSTANT;
  int hEdgePercentile;

  double accumulator;

  Stack<int> hEdgePercentiles;

  numOrigVertices = origGraph->getNumVertices();
  bestPartition.setLength(numOrigVertices);

  for (i = 0; i < numSeqRuns; ++i) {
    Hypergraph *coarseGraph;
    Hypergraph *finerGraph = origGraph;
    hEdgePercentiles.push(startPercentile);

    do {
      hEdgePercentile = hEdgePercentiles.getTopElem();
      coarsener->setPercentile(hEdgePercentile);
      coarseGraph = coarsener->coarsen(*finerGraph);

      if (coarseGraph) {
        hEdgePercentiles.push(min(hEdgePercentile + percentileIncrement, 100));
        hGraphs.push(coarseGraph);
        finerGraph = coarseGraph;
      }
    } while (coarseGraph);

    coarseGraph = hGraphs.pop();
    initBisector->initBisect(*coarseGraph);
    coarseGraph->removeBadPartitions(keepThreshold);

    accumulator = 1.0;

    while (coarseGraph != origGraph) {
      hEdgePercentiles.pop();
      finerGraph = hGraphs.pop();
      finerGraph->projectPartitions(*coarseGraph);

      refiner->refine(*finerGraph);

      finerGraph->removeBadPartitions(keepThreshold * accumulator);

      accumulator *= reductionFactor;

      DynaMem<Hypergraph>::deletePtr(coarseGraph);

      coarseGraph = finerGraph;
    }

    cut = coarseGraph->keepBestPartition();

    if (cut < bestCut) {
      bestCut = cut;
      coarseGraph->copyOutPartition(bestPartition.getArray(), numOrigVertices,
                                    0);
    }

    coarseGraph->resetVertexMaps();

    hGraphs.push(coarseGraph);
  }

  origGraph->setNumPartitions(1);
  origGraph->copyInPartition(bestPartition.getArray(), numOrigVertices, 0,
                             bestCut);
}

void BisectionController::bisect(Hypergraph *h, int _maxPartWt) {
#ifdef DEBUG_CONTROLLER
  assert(hGraphs.getNumElem() == 0);
#endif

  hGraphs.push(h);

  int totWt = h->getTotWeight();
  int maxVertWt;

  double avePartWt = static_cast<double>(totWt) / 2;

  maxPartWt = _maxPartWt;
  maxVertWt =
      static_cast<int>(floor(static_cast<double>(maxPartWt) - avePartWt));

  coarsener->setMaxVertexWt(maxVertWt);
  refiner->setMaxPartWt(maxPartWt);
  initBisector->setMaxPartWt(maxPartWt);

  computeBisection();

  hGraphs.pop();

#ifdef DEBUG_CONTROLLER
  assert(hGraphs.getNumElem() == 0);
#endif
}

#endif
