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
#include "restrictive_first_choice_coarsener.hpp"

using namespace std;

class v_cycle_bisection_controller : public bisection_controller {
protected:
  dynamic_array<int> v_cycle_partition_;
  restrictive_coarsener *restrictive_coarsener_;

public:
  v_cycle_bisection_controller(const int nRuns, const double kT,
                            const double redFactor, int eeParam, int percentile,
                            int inc, int dispL, ostream &out);
  virtual ~v_cycle_bisection_controller();

  void display_options() const;
  void build_restrictive_coarsener(double rRatio, int cType, int minV);
  void record_v_cycle_partition(const int *pVector, int numV);
  void store_best_partition(const int *pVector, int numV);

  virtual void print_type() const = 0;
  virtual void compute_bisection() = 0;
};

#endif
