
#ifndef _VCYCLEFINAL_BISECTION_CPP
#define _VCYCLEFINAL_BISECTION_CPP

// ### VCycleFinalBisectionController.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 4/1/2005: Last Modified
//
// ###

#include "VCycleFinalBisectionController.hpp"

VCycleFinalBisectionController::VCycleFinalBisectionController(
    int nRuns, double kT, double redFactor, int eeParam, int percentile,
    int inc, int dispL, ostream &out)
    : VCycleBisectionController(nRuns, kT, redFactor, eeParam, percentile, inc,
                                dispL, out) {}

VCycleFinalBisectionController::~VCycleFinalBisectionController() {}

void VCycleFinalBisectionController::printVCycleType() const {
  out_stream << " type = BIG";
}

void VCycleFinalBisectionController::computeBisection() {
  Hypergraph *origGraph = hGraphs.top();
  Hypergraph *coarseGraph;
  Hypergraph *finerGraph;

  int i;
  int bestCut = LARGE_CONSTANT;
  int firstCutsize;
  int secondCutsize;
  int diffInCutsize;
  int vCycleGain;
  int hEdgePercentile;

  double accumulator;

  stack<int> hEdgePercentiles;

  numOrigVertices = origGraph->getNumVertices();
  bestPartition.reserve(numOrigVertices);
  vCyclePartition.reserve(numOrigVertices);

  restrCoarsener->setMaxVertexWt(coarsener->getMaxVertexWt());

  for (i = 0; i < numSeqRuns; ++i) {
    finerGraph = origGraph;
    hEdgePercentiles.push(startPercentile);

    // ###
    // coarsen the hypergraph
    // ###

    do {
      hEdgePercentile = hEdgePercentiles.top();
      coarsener->setPercentile(hEdgePercentile);
      coarseGraph = coarsener->coarsen(*finerGraph);

      if (coarseGraph) {
        hEdgePercentiles.push(min(hEdgePercentile + percentileIncrement, 100));
        hGraphs.push(coarseGraph);
        finerGraph = coarseGraph;
      }
    } while (coarseGraph);

    // ###
    // compute the initial partition
    // ###

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

      DynaMem::deletePtr<Hypergraph>(coarseGraph);

      coarseGraph = finerGraph;
    }
#ifdef DEBUG_CONTROLLER
    assert(coarseGraph == origGraph);
#endif

    // ###
    // select the best partition
    // ###

    firstCutsize = coarseGraph->keepBestPartition();

    // ###
    // call v-cycle

    vCycleGain = 0;
    recordVCyclePartition(coarseGraph->getPartVectorArray(), numOrigVertices);

    do {
      coarseGraph->resetMatchVector();
      hEdgePercentiles.push(startPercentile);
      hGraphs.push(coarseGraph);
      finerGraph = coarseGraph;

      // ###
      // coarsen the hypergraph
      // ###

      do {
        hEdgePercentile = hEdgePercentiles.top();
        restrCoarsener->setPercentile(hEdgePercentile);
        coarseGraph = restrCoarsener->coarsen(*finerGraph);

        if (coarseGraph) {
          hEdgePercentiles.push(
              min(hEdgePercentile + percentileIncrement, 100));
          hGraphs.push(coarseGraph);
          finerGraph = coarseGraph;
        }
      }

      while (coarseGraph);

      // ###
      // compute the initial partition
      // ###

      coarseGraph = hGraphs.pop();

      initBisector->initBisect(*coarseGraph);
      coarseGraph->removeBadPartitions(keepThreshold);

      // ###
      // uncoarsen the initial partition
      // ###

      accumulator = 1.0;

      while (coarseGraph != origGraph) {
        hEdgePercentiles.pop();
        finerGraph = hGraphs.pop();
        finerGraph->projectPartitions(*coarseGraph);

        refiner->refine(*finerGraph);

        finerGraph->removeBadPartitions(keepThreshold * accumulator);

        accumulator *= reductionFactor;

        DynaMem::deletePtr<Hypergraph>(coarseGraph);

        coarseGraph = finerGraph;
      }
#ifdef DEBUG_CONTROLLER
      assert(coarseGraph == origGraph);
#endif

      // ###
      // select the best partition
      // ###

      secondCutsize = coarseGraph->keepBestPartition();
      diffInCutsize = firstCutsize - secondCutsize;

      if (firstCutsize > secondCutsize) {

        recordVCyclePartition(coarseGraph->getPartVectorArray(),
                              numOrigVertices);

        firstCutsize = secondCutsize;
        vCycleGain += diffInCutsize;
      }
    } while (diffInCutsize > 0);

    if (firstCutsize < bestCut) {
      bestCut = firstCutsize;
      storeBestPartition(vCyclePartition.data(), numOrigVertices);
    }

    // ###
    // reset hypergraph
    // ###

    origGraph->resetVertexMaps();
    hGraphs.push(coarseGraph);
  }

  origGraph->setNumPartitions(1);
  origGraph->copyInPartition(bestPartition.data(), numOrigVertices, 0,
                             bestCut);
}

#endif
