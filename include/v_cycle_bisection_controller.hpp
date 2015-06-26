#ifndef _VCYCLE_BISECTION_HPP
#define _VCYCLE_BISECTION_HPP
// ### VCycleBisectionController.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 4/1/2005: Last Modified
//
// ###
#include "bisection_controller.hpp"
#include "coarseners/serial/restrictive_first_choice_coarsener.hpp"
#include <string>

namespace parkway {
namespace serial {
namespace ds = parkway::data_structures;

class v_cycle_bisection_controller : public bisection_controller {
 protected:
  ds::dynamic_array<int> v_cycle_partition_;
  restrictive_coarsener *restrictive_coarsener_;
  std::string type_;

 public:
  v_cycle_bisection_controller(const int nRuns, const double kT,
                               const double redFactor, int eeParam,
                               int percentile, int inc, std::string type);
  virtual ~v_cycle_bisection_controller();

  void display_options() const;
  void build_restrictive_coarsener(double rRatio, int cType, int minV);
  void record_v_cycle_partition(ds::dynamic_array<int> pVector,  int n);
  void store_best_partition(ds::dynamic_array<int> pVector, int n);

  virtual void compute_bisection() = 0;
};

}  // namespace serial
}  // namespace parkway

#endif
