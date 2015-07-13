// ### BisectionController.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 4/1/2005: Last Modified
//
// ###
#include "controllers/serial/bisection_controller.hpp"
#include "utility/logging.hpp"

namespace parkway {
namespace serial {

bisection_controller::bisection_controller(int nRuns, double kT, double redFactor,
                                         int eeP, int percentile, int inc) {
  number_of_serial_runs_ = nRuns;
  ee_parameter_ = eeP;
  start_percentile_ = percentile;
  percentile_increment_ = inc;
  keep_threshold_ = kT;
  reduction_factor_ = redFactor;
  number_of_orig_vertices_ = 0;
  maximum_part_weight_ = 0;

  coarsener_ = nullptr;
  refiner_ = nullptr;
  initial_bisector_ = nullptr;
}

bisection_controller::~bisection_controller() {
}


void bisection_controller::display_options() const {
  assert(coarsener_ && initial_bisector_ && refiner_);
  info("|- BSECTOR: kT = %.2f rF = %.2f percentile = %i increment = %i\n|\n",
       keep_threshold_, reduction_factor_, start_percentile_,
       percentile_increment_);
  coarsener_->display_options();
  initial_bisector_->display_options();
  refiner_->display_options();
}

void bisection_controller::build_coarsener(const parkway::options &options) {
  int min_nodes = options.get<int>("coarsening.minimum-coarse-vertices") / 2;
  double reduction_ratio = options.get<double>("coarsening.reduction-ratio");
  int max_weight = -1;
  int fan_out = static_cast<int>(
      options.get<bool>("serial-partitioning.util-fan-out"));
  int divide_by_weight = static_cast<int>(
      options.get<bool>("serial-partitioning.divide-by-cluster-weight"));

  coarsener_ = new first_choice_coarsener(min_nodes, max_weight,
                                          reduction_ratio, fan_out,
                                          divide_by_weight);
}

void bisection_controller::build_initial_bisector(int numInitRuns) {
  initial_bisector_ = new initial_bisector(numInitRuns, FIFO, ee_parameter_);
}

void bisection_controller::build_refiner(int queueD) {
  refiner_ = new fm_refiner(-1, queueD, ee_parameter_);
}

void bisection_controller::compute_bisection() {
  serial::hypergraph *origGraph = hypergraphs_.top();

  int i;
  int cut;
  int bestCut = LARGE_CONSTANT;
  int hEdgePercentile;

  double accumulator;

  std::stack<int> hEdgePercentiles;

  number_of_orig_vertices_ = origGraph->number_of_vertices();
  best_partition_.reserve(number_of_orig_vertices_);

  for (i = 0; i < number_of_serial_runs_; ++i) {
    serial::hypergraph *coarseGraph;
    serial::hypergraph *finerGraph = origGraph;
    hEdgePercentiles.push(start_percentile_);

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

    cut = coarseGraph->keep_best_partition();

    if (cut < bestCut) {
      bestCut = cut;
      coarseGraph->copy_out_partition(best_partition_, number_of_orig_vertices_,
                                      0);
    }

    coarseGraph->reset_vertex_maps();

    hypergraphs_.push(coarseGraph);
  }

  origGraph->set_number_of_partitions(1);
  origGraph->copy_in_partition(best_partition_, number_of_orig_vertices_, 0,
                               bestCut);
}

void bisection_controller::bisect(serial::hypergraph *h, int _maxPartWt) {
#ifdef DEBUG_CONTROLLER
  assert(hGraphs.getNumElem() == 0);
#endif

  hypergraphs_.push(h);

  int totWt = h->total_weight();
  int maxVertWt;

  double avePartWt = static_cast<double>(totWt) / 2;

  maximum_part_weight_ = _maxPartWt;
  maxVertWt =
      static_cast<int>(floor(static_cast<double>(maximum_part_weight_) - avePartWt));

  coarsener_->set_maximum_vertex_weight(maxVertWt);
  refiner_->set_maximum_part_weight(maximum_part_weight_);
  initial_bisector_->set_maximum_part_weight(maximum_part_weight_);

  compute_bisection();

  hypergraphs_.pop();

#ifdef DEBUG_CONTROLLER
  assert(hGraphs.getNumElem() == 0);
#endif
}

}  // namespace serial
}  // namespace parkway
