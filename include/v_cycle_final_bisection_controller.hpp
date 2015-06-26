#ifndef _VCYCLEFINAL_BISECTION_HPP
#define _VCYCLEFINAL_BISECTION_HPP
// ### VCycleFinalBisectionController.hpp ###
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


class v_cycle_final_bisection_controller : public v_cycle_bisection_controller {
 public:
  v_cycle_final_bisection_controller(int nRuns, double kT, double redFactor,
                                     int eeParam, int percentile, int inc);
  ~v_cycle_final_bisection_controller();

  void compute_bisection();
};

}  // namespace serial
}  // namespace parkway

#endif
