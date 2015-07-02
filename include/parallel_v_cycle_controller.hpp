#ifndef _VCYCLE_CONTROLLER_HPP
#define _VCYCLE_CONTROLLER_HPP
// ### ParaVCycleController.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 4/1/2005: Last Modified
//
// ###
#include <stack>
#include "internal/parallel_controller.hpp"
#include "hypergraph/parallel/hypergraph.hpp"
#include "coarseners/parallel/restrictive_first_choice_coarsening.hpp"
#include "data_structures/dynamic_array.hpp"

namespace ds = parkway::data_structures;
namespace parallel = parkway::parallel;

class parallel_v_cycle_controller : public parallel::controller {
 protected:
  int limit_on_cycles_;
  int minimum_inter_vertex_index_;

  double limit_as_percent_of_cut_;

  std::stack<int> number_of_local_current_vertices_;
  std::stack<int> minimum_local_current_vertices_;
  std::stack<ds::dynamic_array<int> *> best_v_cycle_partition_;

  dynamic_array<int> map_to_inter_vertices_;

  parallel::restrictive_coarsening &restrictive_coarsening_;

 public:
  parallel_v_cycle_controller(parallel::restrictive_coarsening &rc,
                              parallel::coarsener &c, parallel::refiner &r,
                              parkway::serial::controller &ref, int rank, int nP,
                              int percentile, int inc, int approxRef, int limit,
                              double limitAsPercent, std::ostream &out);

  virtual ~parallel_v_cycle_controller();
  virtual void run(MPI_Comm comm) = 0;
  virtual void print_type() const = 0;

  void set_weight_constraints(MPI_Comm comm);
  void display_options() const;

  void record_v_cycle_partition(const parallel::hypergraph &h, int numIteration);
  void gather_in_v_cycle_partition(parallel::hypergraph &h, int cut, MPI_Comm comm);

  void project_v_cycle_partition(parallel::hypergraph &cG, parallel::hypergraph &fG,
                                 MPI_Comm comm);
  void shuffle_v_cycle_vertices_by_partition(parallel::hypergraph &h, MPI_Comm comm);

  void shift_v_cycle_vertices_to_balance(parallel::hypergraph &h, MPI_Comm comm);
  void update_map_to_orig_vertices(MPI_Comm comm);
  void reset_structures();
};

#endif
