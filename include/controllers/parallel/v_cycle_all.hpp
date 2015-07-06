#ifndef _VCYCLEALL_CONTROLLER_HPP
#define _VCYCLEALL_CONTROLLER_HPP
// ### ParaVCycleController.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 4/1/2005: Last Modified
//
// ###
#include "controllers/parallel/v_cycle.hpp"

namespace parkway {
namespace parallel {

class v_cycle_all : public parkway::parallel::v_cycle {
 public:
  v_cycle_all(parallel::restrictive_coarsening &rc, parallel::coarsener &c,
              parallel::refiner &r, parkway::serial::controller &ref, int rank,
              int nP, int percentile, int inc, int approxRef, int limit,
              double limitAsPercent);
  ~v_cycle_all();

  void run(MPI_Comm comm);
  void print_type() const;
};

}  // namespace parallel
}  // namespace parkway

#endif
