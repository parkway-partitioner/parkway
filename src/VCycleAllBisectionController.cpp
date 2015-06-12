
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
  hypergraph *origGraph = hGraphs.top();
  hypergraph *coarseGraph;
  hypergraph *finerGraph;
  hypergraph *interMedGraph;

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

  stack<int> hEdgePercentiles;

  numOrigVertices = origGraph->number_of_vertices();
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
        coarsener->set_percentile(hEdgePercentile);
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
        coarseGraph->remove_bad_partitions(keepThreshold * accumulator);

      hEdgePercentiles.pop();
      finerGraph = hGraphs.pop();
        finerGraph->project_partitions(*coarseGraph);

      refiner->refine(*finerGraph);

      accumulator *= reductionFactor;

      DynaMem::deletePtr<hypergraph>(coarseGraph);

      // ###
      // prepare to call v-cycle

      vCycleGain = 0;
      firstCutSize = finerGraph->keep_best_partition();
      recordVCyclePartition(finerGraph->partition_vector(),
                            finerGraph->number_of_vertices());

      numInStack = hGraphs.size();
      interMedGraph = finerGraph;

      do {
          finerGraph->reset_match_vector();
        hGraphs.push(finerGraph);

        othAccumulator = 1.0;

        // ###
        // coarsen the hypergraph
        // ###

        do {
          hEdgePercentile = hEdgePercentiles.top();
            restrCoarsener->set_percentile(hEdgePercentile);
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
            coarseGraph->remove_bad_partitions(keepThreshold * othAccumulator);
          othAccumulator *= reductionFactor;

          hEdgePercentiles.pop();
          finerGraph = hGraphs.pop();
            finerGraph->project_partitions(*coarseGraph);

          refiner->refine(*finerGraph);

          DynaMem::deletePtr<hypergraph>(coarseGraph);

          coarseGraph = finerGraph;
        }
#ifdef DEBUG_CONTROLLER
        assert(interMedGraph == coarseGraph);
#endif
        secondCutSize = coarseGraph->keep_best_partition();
        diffInCutSize = firstCutSize - secondCutSize;

        if (firstCutSize > secondCutSize) {
          recordVCyclePartition(coarseGraph->partition_vector(),
                                coarseGraph->number_of_vertices());
          firstCutSize = secondCutSize;
          vCycleGain += diffInCutSize;
        }
      } while (diffInCutSize > 0);

      // ###
      // end v-cycle

        coarseGraph->copy_in_partition(vCyclePartition.data(),
                                       coarseGraph->number_of_vertices(), 0,
                                       firstCutSize);
    }
#ifdef DEBUG_CONTROLLER
    assert(coarseGraph == origGraph);
#endif

    firstCutSize = coarseGraph->keep_best_partition();

    if (firstCutSize < bestCut) {
      bestCut = firstCutSize;
      storeBestPartition(coarseGraph->partition_vector(), numOrigVertices);
    }

    // ###
    // reset hypergraph
    // ###

      origGraph->copy_in_partition(bestPartition.data(), numOrigVertices, 0,
                                   bestCut);
      origGraph->reset_vertex_maps();
    hGraphs.push(coarseGraph);
  }

    origGraph->set_number_of_partitions(1);
    origGraph->copy_in_partition(bestPartition.data(), numOrigVertices, 0,
                                 bestCut);
}

#endif
