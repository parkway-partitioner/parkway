// ### VCycleFinalBisectionController.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 4/1/2005: Last Modified
//
// ###
#include "v_cycle_final_bisection_controller.hpp"
#include "utility/logging.hpp"

namespace parkway {
namespace serial {

v_cycle_final_bisection_controller::v_cycle_final_bisection_controller(
    int nRuns, double kT, double redFactor, int eeParam, int percentile,
    int inc)
    : v_cycle_bisection_controller(nRuns, kT, redFactor, eeParam, percentile, inc) {
}

v_cycle_final_bisection_controller::~v_cycle_final_bisection_controller() {}

void v_cycle_final_bisection_controller::print_type() const {
  info(" type = BIG");
}

void v_cycle_final_bisection_controller::compute_bisection() {
  serial::hypergraph *origGraph = hypergraphs_.top();
  serial::hypergraph *coarseGraph;
  serial::hypergraph *finerGraph;

  int i;
  int bestCut = LARGE_CONSTANT;
  int firstCutsize;
  int secondCutsize;
  int diffInCutsize;
  int vCycleGain;
  int hEdgePercentile;

  double accumulator;

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
      coarseGraph->remove_bad_partitions(keep_threshold_);

    accumulator = 1.0;

    while (coarseGraph != origGraph) {
      hEdgePercentiles.pop();
      finerGraph = hypergraphs_.top();
      hypergraphs_.pop();
      finerGraph->project_partitions(*coarseGraph);

      refiner_->refine(*finerGraph);

        finerGraph->remove_bad_partitions(keep_threshold_ * accumulator);

      accumulator *= reduction_factor_;

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
    record_v_cycle_partition(coarseGraph->partition_vector(),
                             number_of_orig_vertices_);

    do {
        coarseGraph->reset_match_vector();
      hEdgePercentiles.push(start_percentile_);
      hypergraphs_.push(coarseGraph);
      finerGraph = coarseGraph;

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
      }

      while (coarseGraph);

      // ###
      // compute the initial partition
      // ###

      coarseGraph = hypergraphs_.top();
      hypergraphs_.pop();

      initial_bisector_->initialize_bisector(*coarseGraph);
        coarseGraph->remove_bad_partitions(keep_threshold_);

      // ###
      // uncoarsen the initial partition
      // ###

      accumulator = 1.0;

      while (coarseGraph != origGraph) {
        hEdgePercentiles.pop();
        finerGraph = hypergraphs_.top();
        hypergraphs_.pop();
        finerGraph->project_partitions(*coarseGraph);

        refiner_->refine(*finerGraph);

          finerGraph->remove_bad_partitions(keep_threshold_ * accumulator);

        accumulator *= reduction_factor_;

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

        record_v_cycle_partition(coarseGraph->partition_vector(),
                                 number_of_orig_vertices_);

        firstCutsize = secondCutsize;
        vCycleGain += diffInCutsize;
      }
    } while (diffInCutsize > 0);

    if (firstCutsize < bestCut) {
      bestCut = firstCutsize;
      store_best_partition(v_cycle_partition_, number_of_orig_vertices_);
    }

    // ###
    // reset hypergraph
    // ###

      origGraph->reset_vertex_maps();
    hypergraphs_.push(coarseGraph);
  }

    origGraph->set_number_of_partitions(1);
    origGraph->copy_in_partition(best_partition_, number_of_orig_vertices_, 0,
                                 bestCut);
}

}  // namespace serial
}  // namespace parkway
