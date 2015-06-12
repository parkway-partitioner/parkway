
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
  serial_hypergraph *origGraph = hGraphs.top();
  serial_hypergraph *coarseGraph;
  serial_hypergraph *finerGraph;

  int i;
  int bestCut = LARGE_CONSTANT;
  int firstCutsize;
  int secondCutsize;
  int diffInCutsize;
  int vCycleGain;
  int hEdgePercentile;

  double accumulator;

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
      coarseGraph->remove_bad_partitions(keepThreshold);

    accumulator = 1.0;

    while (coarseGraph != origGraph) {
      hEdgePercentiles.pop();
      finerGraph = hGraphs.pop();
        finerGraph->project_partitions(*coarseGraph);

      refiner->refine(*finerGraph);

        finerGraph->remove_bad_partitions(keepThreshold * accumulator);

      accumulator *= reductionFactor;

      DynaMem::deletePtr<serial_hypergraph>(coarseGraph);

      coarseGraph = finerGraph;
    }
#ifdef DEBUG_CONTROLLER
    assert(coarseGraph == origGraph);
#endif

    // ###
    // select the best partition
    // ###

    firstCutsize = coarseGraph->keep_best_partition();

    // ###
    // call v-cycle

    vCycleGain = 0;
    recordVCyclePartition(coarseGraph->partition_vector(), numOrigVertices);

    do {
        coarseGraph->reset_match_vector();
      hEdgePercentiles.push(startPercentile);
      hGraphs.push(coarseGraph);
      finerGraph = coarseGraph;

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
      }

      while (coarseGraph);

      // ###
      // compute the initial partition
      // ###

      coarseGraph = hGraphs.pop();

      initBisector->initBisect(*coarseGraph);
        coarseGraph->remove_bad_partitions(keepThreshold);

      // ###
      // uncoarsen the initial partition
      // ###

      accumulator = 1.0;

      while (coarseGraph != origGraph) {
        hEdgePercentiles.pop();
        finerGraph = hGraphs.pop();
          finerGraph->project_partitions(*coarseGraph);

        refiner->refine(*finerGraph);

          finerGraph->remove_bad_partitions(keepThreshold * accumulator);

        accumulator *= reductionFactor;

        DynaMem::deletePtr<serial_hypergraph>(coarseGraph);

        coarseGraph = finerGraph;
      }
#ifdef DEBUG_CONTROLLER
      assert(coarseGraph == origGraph);
#endif

      // ###
      // select the best partition
      // ###

      secondCutsize = coarseGraph->keep_best_partition();
      diffInCutsize = firstCutsize - secondCutsize;

      if (firstCutsize > secondCutsize) {

        recordVCyclePartition(coarseGraph->partition_vector(),
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

      origGraph->reset_vertex_maps();
    hGraphs.push(coarseGraph);
  }

    origGraph->set_number_of_partitions(1);
    origGraph->copy_in_partition(bestPartition.data(), numOrigVertices, 0,
                                 bestCut);
}

#endif
