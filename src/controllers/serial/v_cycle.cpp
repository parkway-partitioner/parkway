// ### VCycleBisectionController.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 4/1/2005: Last Modified
//
// ###
#include "controllers/serial/v_cycle.hpp"
#include "utility/logging.hpp"

namespace parkway {
namespace serial {

v_cycle::v_cycle(
    const int nRuns, const double kT, const double redFactor, int eeParam,
    int percentile, int inc)
    : bisection_controller(nRuns, kT, redFactor, eeParam, percentile, inc) {
  restrictive_coarsener_ = nullptr;
  v_cycle_partition_.reserve(0);
}

v_cycle::~v_cycle() {
}

void v_cycle::display_options() const {
  assert(coarsener_ && initial_bisector_ && refiner_ && restrictive_coarsener_);
  info("|- V-CYC BSECTOR: kT = %.2f rF = %.2f percentile = %i increment = %i",
       keep_threshold_, reduction_factor_, start_percentile_,
       percentile_increment_);
  print_type();
  info("\n|\n");
  coarsener_->display_options();
  restrictive_coarsener_->display_options();
  initial_bisector_->display_options();
  refiner_->display_options();
}

void v_cycle::build_restrictive_coarsener(double redRatio,
                                                               int cType,
                                                               int minNodes) {
  switch (cType) {
  case RestrFCwithFanOutDiv:
    restrictive_coarsener_ = new restrictive_first_choice_coarsener(
        (minNodes << 1), -1, redRatio, 1, 1);
    break;

  case RestrFCwithoutFanOutDiv:
    restrictive_coarsener_ = new restrictive_first_choice_coarsener(
        (minNodes << 1), -1, redRatio, 0, 1);
    break;

  case RestrFCwithFanOut:
    restrictive_coarsener_ = new restrictive_first_choice_coarsener(
        (minNodes << 1), -1, redRatio, 1, 0);
    break;

  case RestrFCwithoutFanOut:
    restrictive_coarsener_ = new restrictive_first_choice_coarsener(
        (minNodes << 1), -1, redRatio, 0, 0);
    break;

  default:
    restrictive_coarsener_ = new restrictive_first_choice_coarsener(
        (minNodes << 1), -1, redRatio, 1, 1);
    break;
  }
}

void v_cycle::record_v_cycle_partition(
    ds::dynamic_array<int> pVector, int numV) {
  v_cycle_partition_ = pVector;
  for (int i = 0; i < numV; ++i) {
    v_cycle_partition_[i] = pVector[i];
  }
}

void v_cycle::store_best_partition(
    ds::dynamic_array<int> pVector, int numV) {
  for (int i = 0; i < numV; ++i) {
    best_partition_[i] = pVector[i];
  }
}

}  // namespace serial
}  // namespace parkway
