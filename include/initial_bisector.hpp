#ifndef _INIT_BISECTOR_HPP
#define _INIT_BISECTOR_HPP

// ### InitBisector.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 30/11/2004: Last Modified
//
// ###

#include <iostream>
#include "refiners/serial/fm_refiner.hpp"
#include "hypergraph/serial/hypergraph.hpp"

namespace parkway {
namespace serial {

class initial_bisector : public fm_refiner {
 protected:
  int number_of_initial_runs_;

 public:
  initial_bisector(int nRuns, int insMethod, int ee, int dL);
  ~initial_bisector();

  void display_options() const;
  void initialize_bisector(serial::hypergraph &h);

  int set_base_vertex();
  int choose_best_vertex_1_to_0();
  int greedy_pass();

  inline void set_number_of_runs(int nR) {
    number_of_initial_runs_ = nR;
  }
};

}  // namespace serial
}  // namespace parkway

#endif
