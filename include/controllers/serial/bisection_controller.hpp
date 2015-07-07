#ifndef _BISECTION_CONTROLLER_HPP
#define _BISECTION_CONTROLLER_HPP
// ### BisectionController.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 4/1/2005: Last Modified
//
// ###
#include <ostream>
#include <stack>
#include "hypergraph/serial/hypergraph.hpp"
#include "coarseners/serial/first_choice_coarsener.hpp"
#include "controllers/serial/initial_bisector.hpp"

namespace parkway {
namespace serial {
namespace ds = parkway::data_structures;

class bisection_controller {
 protected:
  int number_of_serial_runs_;
  int ee_parameter_;

  int start_percentile_;
  int percentile_increment_;

  int number_of_orig_vertices_;
  int maximum_part_weight_;

  double keep_threshold_;
  double reduction_factor_;

  coarsener *coarsener_;
  refiner *refiner_;
  initial_bisector *initial_bisector_;
  std::stack<serial::hypergraph *> hypergraphs_;
  ds::dynamic_array<int> best_partition_;

 public:
  bisection_controller(int nRuns, double kT, double redFactor, int eeParam,
                       int percentile, int inc);

  virtual ~bisection_controller();
  virtual void display_options() const;

  void build_coarsener(double redRatio, int cType, int minNodes);
  void build_initial_bisector(int numInitRuns);
  void build_refiner(int queueD);

  virtual void compute_bisection();

  void bisect(serial::hypergraph *h, int maxPartWt);

  inline void set_number_of_serial_runs(int r) { number_of_serial_runs_ = r; }
};

}  // namespace serial
}  // namespace parkway

#endif
