// ### VCycleAllBisectionController.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 4/1/2005: Last Modified
//
// ###
#include "controllers/serial/v_cycle_all.hpp"
#include "utility/logging.hpp"

namespace parkway {
namespace serial {

v_cycle_all::v_cycle_all(
    int nRuns, double kT, double redFactor, int eeParam, int percentile,
    int inc)
    : parkway::serial::v_cycle(nRuns, kT, redFactor, eeParam, percentile, inc) {
}

v_cycle_all::~v_cycle_all() {}

void v_cycle_all::print_type() const {
  info(" type = ALL");
}

void v_cycle_all::compute_bisection() {
  serial::hypergraph *origGraph = hypergraphs_.top();
  serial::hypergraph *coarseGraph;
  serial::hypergraph *finerGraph;
  serial::hypergraph *interMedGraph;

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

  std::stack<int> hEdgePercentiles;

  number_of_orig_vertices_ = origGraph->number_of_vertices();
  best_partition_.reserve(number_of_orig_vertices_);
  v_cycle_partition_.reserve(number_of_orig_vertices_);

  restrictive_coarsener_->set_maximum_vertex_weight(coarsener_->maximum_vertex_weight());

  for (i = 0; i < number_of_serial_runs_; ++i) {
    finerGraph = origGraph;
    hEdgePercentiles.push(start_percentile_);

    // ###
    // coarsen the hypergraph
    // ###

    do {
      hEdgePercentile = hEdgePercentiles.top();
      coarsener_->set_percentile(hEdgePercentile);
      coarseGraph = coarsener_->coarsen(*finerGraph);

      if (coarseGraph) {
        hEdgePercentiles.push(std::min(hEdgePercentile + percentile_increment_, 100));
        hypergraphs_.push(coarseGraph);
        finerGraph = coarseGraph;
      }
    } while (coarseGraph);

    // ###
    // compute the initial partition
    // ###

    coarseGraph = hypergraphs_.top();
    hypergraphs_.pop();
    initial_bisector_->initialize_bisector(*coarseGraph);

    accumulator = 1.0;

    while (coarseGraph != origGraph) {
        coarseGraph->remove_bad_partitions(keep_threshold_ * accumulator);

      hEdgePercentiles.pop();
      finerGraph = hypergraphs_.top();
      hypergraphs_.pop();
      finerGraph->project_partitions(*coarseGraph);

      refiner_->refine(*finerGraph);

      accumulator *= reduction_factor_;

      // ###
      // prepare to call v-cycle

      vCycleGain = 0;
      firstCutSize = finerGraph->keep_best_partition();
      record_v_cycle_partition(finerGraph->partition_vector(),
                               finerGraph->number_of_vertices());

      numInStack = hypergraphs_.size();
      interMedGraph = finerGraph;

      do {
          finerGraph->reset_match_vector();
        hypergraphs_.push(finerGraph);

        othAccumulator = 1.0;

        // ###
        // coarsen the hypergraph
        // ###

        do {
          hEdgePercentile = hEdgePercentiles.top();
            restrictive_coarsener_->set_percentile(hEdgePercentile);
          coarseGraph = restrictive_coarsener_->coarsen(*finerGraph);

          if (coarseGraph) {
            hEdgePercentiles.push(
                std::min(hEdgePercentile + percentile_increment_, 100));
            hypergraphs_.push(coarseGraph);
            finerGraph = coarseGraph;
          }
        } while (coarseGraph);

        // ###
        // compute the initial partition
        // ###

        coarseGraph = hypergraphs_.top();
        hypergraphs_.pop();
        initial_bisector_->initialize_bisector(*coarseGraph);

        while (coarseGraph != interMedGraph) {
            coarseGraph->remove_bad_partitions(keep_threshold_ * othAccumulator);
          othAccumulator *= reduction_factor_;

          hEdgePercentiles.pop();
          finerGraph = hypergraphs_.top();
          hypergraphs_.pop();
          finerGraph->project_partitions(*coarseGraph);

          refiner_->refine(*finerGraph);

          coarseGraph = finerGraph;
        }
#ifdef DEBUG_CONTROLLER
        assert(interMedGraph == coarseGraph);
#endif
        secondCutSize = coarseGraph->keep_best_partition();
        diffInCutSize = firstCutSize - secondCutSize;

        if (firstCutSize > secondCutSize) {
          record_v_cycle_partition(coarseGraph->partition_vector(),
                                   coarseGraph->number_of_vertices());
          firstCutSize = secondCutSize;
          vCycleGain += diffInCutSize;
        }
      } while (diffInCutSize > 0);

      // ###
      // end v-cycle

        coarseGraph->copy_in_partition(v_cycle_partition_,
                                       coarseGraph->number_of_vertices(), 0,
                                       firstCutSize);
    }
#ifdef DEBUG_CONTROLLER
    assert(coarseGraph == origGraph);
#endif

    firstCutSize = coarseGraph->keep_best_partition();

    if (firstCutSize < bestCut) {
      bestCut = firstCutSize;
      store_best_partition(coarseGraph->partition_vector(),
                           number_of_orig_vertices_);
    }

    // ###
    // reset hypergraph
    // ###

      origGraph->copy_in_partition(best_partition_,
                                   number_of_orig_vertices_, 0,
                                   bestCut);
      origGraph->reset_vertex_maps();
    hypergraphs_.push(coarseGraph);
  }

    origGraph->set_number_of_partitions(1);
    origGraph->copy_in_partition(best_partition_, number_of_orig_vertices_, 0,
                                 bestCut);
}

}  // namespace serial
}  // namespace parkway
