
#ifndef _VCYCLEALL_BISECTION_CPP
#define _VCYCLEALL_BISECTION_CPP

// ### VCycleAllBisectionController.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 4/1/2005: Last Modified
//
// ###

#include "VCycleAllBisectionController.hpp"

VCycleAllBisectionController::VCycleAllBisectionController(
    int nRuns, double kT, double redFactor, int eeParam, int percentile,
    int inc, int dispL, ostream &out)
    : VCycleBisectionController(nRuns, kT, redFactor, eeParam, percentile, inc,
                                dispL, out) {}

VCycleAllBisectionController::~VCycleAllBisectionController() {}

void VCycleAllBisectionController::printVCycleType() const {
  out_stream << " type = ALL";
}

void VCycleAllBisectionController::computeBisection() {
  Hypergraph *origGraph = hGraphs.getTopElem();
  Hypergraph *coarseGraph;
  Hypergraph *finerGraph;
  Hypergraph *interMedGraph;

  int i;
  int numInStack;
  int bestCut = LARGE_CONSTANT;
  int firstCutSize;
  int secondCutSize;
  int diffInCutSize;
  int vCycleGain;
  int hEdgePercentile;

  double accumulator;
  double othAccumulator;

  Stack<int> hEdgePercentiles;

  numOrigVertices = origGraph->getNumVertices();
  bestPartition.setLength(numOrigVertices);
  vCyclePartition.setLength(numOrigVertices);

  restrCoarsener->setMaxVertexWt(coarsener->getMaxVertexWt());

  for (i = 0; i < numSeqRuns; ++i) {
    finerGraph = origGraph;
    hEdgePercentiles.push(startPercentile);

    // ###
    // coarsen the hypergraph
    // ###

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

    // ###
    // compute the initial partition
    // ###

    coarseGraph = hGraphs.pop();
    initBisector->initBisect(*coarseGraph);

    accumulator = 1.0;

    while (coarseGraph != origGraph) {
      coarseGraph->removeBadPartitions(keepThreshold * accumulator);

      hEdgePercentiles.pop();
      finerGraph = hGraphs.pop();
      finerGraph->projectPartitions(*coarseGraph);

      refiner->refine(*finerGraph);

      accumulator *= reductionFactor;

      DynaMem<Hypergraph>::deletePtr(coarseGraph);

      // ###
      // prepare to call v-cycle

      vCycleGain = 0;
      firstCutSize = finerGraph->keepBestPartition();
      recordVCyclePartition(finerGraph->getPartVectorArray(),
                            finerGraph->getNumVertices());

      numInStack = hGraphs.getNumElem();
      interMedGraph = finerGraph;

      do {
        finerGraph->resetMatchVector();
        hGraphs.push(finerGraph);

        othAccumulator = 1.0;

        // ###
        // coarsen the hypergraph
        // ###

        do {
          hEdgePercentile = hEdgePercentiles.getTopElem();
          restrCoarsener->setPercentile(hEdgePercentile);
          coarseGraph = restrCoarsener->coarsen(*finerGraph);

          if (coarseGraph) {
            hEdgePercentiles.push(
                min(hEdgePercentile + percentileIncrement, 100));
            hGraphs.push(coarseGraph);
            finerGraph = coarseGraph;
          }
        } while (coarseGraph);

        // ###
        // compute the initial partition
        // ###

        coarseGraph = hGraphs.pop();
        initBisector->initBisect(*coarseGraph);

        while (coarseGraph != interMedGraph) {
          coarseGraph->removeBadPartitions(keepThreshold * othAccumulator);
          othAccumulator *= reductionFactor;

          hEdgePercentiles.pop();
          finerGraph = hGraphs.pop();
          finerGraph->projectPartitions(*coarseGraph);

          refiner->refine(*finerGraph);

          DynaMem<Hypergraph>::deletePtr(coarseGraph);

          coarseGraph = finerGraph;
        }
#ifdef DEBUG_CONTROLLER
        assert(interMedGraph == coarseGraph);
#endif
        secondCutSize = coarseGraph->keepBestPartition();
        diffInCutSize = firstCutSize - secondCutSize;

        if (firstCutSize > secondCutSize) {
          recordVCyclePartition(coarseGraph->getPartVectorArray(),
                                coarseGraph->getNumVertices());
          firstCutSize = secondCutSize;
          vCycleGain += diffInCutSize;
        }
      } while (diffInCutSize > 0);

      // ###
      // end v-cycle

      coarseGraph->copyInPartition(vCyclePartition.getArray(),
                                   coarseGraph->getNumVertices(), 0,
                                   firstCutSize);
    }
#ifdef DEBUG_CONTROLLER
    assert(coarseGraph == origGraph);
#endif

    firstCutSize = coarseGraph->keepBestPartition();

    if (firstCutSize < bestCut) {
      bestCut = firstCutSize;
      storeBestPartition(coarseGraph->getPartVectorArray(), numOrigVertices);
    }

    // ###
    // reset hypergraph
    // ###

    origGraph->copyInPartition(bestPartition.getArray(), numOrigVertices, 0,
                               bestCut);
    origGraph->resetVertexMaps();
    hGraphs.push(coarseGraph);
  }

  origGraph->setNumPartitions(1);
  origGraph->copyInPartition(bestPartition.getArray(), numOrigVertices, 0,
                             bestCut);
}

#endif
