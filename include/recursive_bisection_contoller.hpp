#ifndef _RECUR_BISECT_CONTROLLER_HPP
#define _RECUR_BISECT_CONTROLLER_HPP

// ### RecurBisectController.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 30/11/2004: Last Modified
//
// NOTES:
//
// implement the algorithms such that it is not imperative
// that partition balance is maintained during recursive
// bisection. The greedy k-way refiner can be used to rebalance
// an unbalanced partition
//
// ###

#include "internal/serial_controller.hpp"
#include "v_cycle_final_bisection_controller.hpp"
#include "v_cycle_all_bisection_controller.hpp"
#include "refiners/serial/greedy_k_way_refiner.hpp"
#include "bisection.hpp"
#include "hypergraph/parallel/hypergraph.hpp"

namespace parallel = parkway::parallel;

class recursive_bisection_contoller : public parkway::serial::controller {
 protected:
  int number_of_bisections_;
  int log_k_;
  int maximum_part_weight_;
  int sum_of_cuts_;
  int local_vertex_part_info_length_;

  int number_of_partitions_;

  double bisection_constraint_;
  double average_part_weight_;
  double average_initial_bisection_weight_;

  parkway::serial::bisection_controller *bisector_;
  parkway::serial::greedy_k_way_refiner *refiner_;

  dynamic_array<int> local_vertex_partition_info_;
  dynamic_array<int> all_partition_info_;

 public:
  recursive_bisection_contoller(parkway::serial::bisection_controller *b,
                                parkway::serial::greedy_k_way_refiner *k,
                                int rank, int nProcs, int nParts,
                                int nBisectRuns);
  ~recursive_bisection_contoller();

  void display_options() const;
  void convToBisectionConstraints();

  void run(parallel::hypergraph &hgraph, MPI_Comm comm);
  void initialize_serial_partitions(parallel::hypergraph &hgraph,
                                        MPI_Comm comm);
  void recursively_bisect(const bisection &b, MPI_Comm comm);

  void split_bisection(const bisection &b, bisection *&newB,
                       MPI_Comm comm) const;
  void split_bisection(const bisection &b, bisection *&l, bisection *&r) const;

  int best_partition_processor(int cut, MPI_Comm comm) const;
  int compute_maximum_weight(int numBisections) const;
  double recursively_compute_maximum(double ave, int depth) const;

  inline void set_number_of_bisections(int bRuns) { number_of_bisections_ = bRuns; }
};

#endif
