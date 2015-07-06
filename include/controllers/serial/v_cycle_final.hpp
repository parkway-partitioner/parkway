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
#include "controllers/serial/v_cycle.hpp"

namespace parkway {
namespace serial {

class v_cycle_final : public parkway::serial::v_cycle {
 public:
  v_cycle_final(int nRuns, double kT, double redFactor, int eeParam,
                int percentile, int inc);
  ~v_cycle_final();

  void print_type() const;
  void compute_bisection();
};

}  // namespace serial
}  // namespace parkway

#endif
