// ### VCycleBisectionController.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 4/1/2005: Last Modified
//
// ###

#include "v_cycle_bisection_controller.hpp"

namespace parkway {
namespace serial {

v_cycle_bisection_controller::v_cycle_bisection_controller(
    const int nRuns, const double kT, const double redFactor, int eeParam,
    int percentile, int inc, int dispL, std::ostream &out)
    : bisection_controller(nRuns, kT, redFactor, eeParam, percentile, inc, dispL,
                          out) {
  restrictive_coarsener_ = nullptr;

  v_cycle_partition_.reserve(0);
}

v_cycle_bisection_controller::~v_cycle_bisection_controller() {
}

void v_cycle_bisection_controller::display_options() const {
  switch (display_level_) {
  case SILENT:
    break;

  default:

    out_stream << "|- V-CYC BSECTOR:";
#ifdef DEBUG_CONTROLLER
    assert(coarsener && initBisector && refiner && restrCoarsener);
#endif
    out_stream << " kT = " << keep_threshold_ << " rF = " << reduction_factor_
               << " %le = " << start_percentile_
               << " %inc = " << percentile_increment_;
      print_type();
    out_stream << std::endl << "|" << std::endl;
      coarsener_->display_options(out_stream);
      restrictive_coarsener_->display_options(out_stream);
      initial_bisector_->display_options(out_stream);
      refiner_->display_options(out_stream);

    break;
  }
}

void v_cycle_bisection_controller::build_restrictive_coarsener(double redRatio,
                                                               int cType,
                                                               int minNodes) {
  switch (cType) {
  case RestrFCwithFanOutDiv:
    restrictive_coarsener_ = new restrictive_first_choice_coarsener(Shiftl(minNodes, 1), -1, redRatio, 1,
                                          1, display_level_);
    break;

  case RestrFCwithoutFanOutDiv:
    restrictive_coarsener_ = new restrictive_first_choice_coarsener(Shiftl(minNodes, 1), -1, redRatio, 0,
                                          1, display_level_);
    break;

  case RestrFCwithFanOut:
    restrictive_coarsener_ = new restrictive_first_choice_coarsener(Shiftl(minNodes, 1), -1, redRatio, 1,
                                          0, display_level_);
    break;

  case RestrFCwithoutFanOut:
    restrictive_coarsener_ = new restrictive_first_choice_coarsener(Shiftl(minNodes, 1), -1, redRatio, 0,
                                          0, display_level_);
    break;

  default:
    restrictive_coarsener_ = new restrictive_first_choice_coarsener(Shiftl(minNodes, 1), -1, redRatio, 1,
                                          1, display_level_);
    break;
  }
}

void v_cycle_bisection_controller::record_v_cycle_partition(
    ds::dynamic_array<int> pVector, int numV) {
  v_cycle_partition_ = pVector;
  for (int i = 0; i < numV; ++i) {
    v_cycle_partition_[i] = pVector[i];
  }
}

void v_cycle_bisection_controller::store_best_partition(
    ds::dynamic_array<int> pVector, int numV) {
  for (int i = 0; i < numV; ++i) {
    best_partition_[i] = pVector[i];
  }
}

}  // namespace serial
}  // namespace parkway
