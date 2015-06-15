
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

#include "parallel_v_cycle_controller.hpp"

using namespace std;

class parallel_v_cycle_all_controller : public parallel_v_cycle_controller {

protected:
public:
  parallel_v_cycle_all_controller(parallel_restrictive_coarsening &rc, parallel_coarsener &c,
                          parallel_refiner &r, sequential_controller &ref, int rank, int nP,
                          int percentile, int inc, int approxRef, int limit,
                          double limitAsPercent, ostream &out);
  ~parallel_v_cycle_all_controller();

  void run(MPI_Comm comm);
  void print_type() const;
};

#endif
