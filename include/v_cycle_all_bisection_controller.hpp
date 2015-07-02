#ifndef _VCYCLEALL_BISECTION_HPP
#define _VCYCLEALL_BISECTION_HPP
// ### VCycleAllBisectionController.hpp ###
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

class v_cycle_all_bisection_controller : public v_cycle_bisection_controller {
 public:
  v_cycle_all_bisection_controller(int nRuns, double kT, double redFactor,
                                   int eeParam, int percentile, int inc,
                                   int dispL, std::ostream &out);
  ~v_cycle_all_bisection_controller();

  void print_type() const;
  void compute_bisection();
};

}  // namespace serial
}  // namespace parkway

#endif
