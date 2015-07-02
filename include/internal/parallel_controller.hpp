#ifndef _PARA_CONTROLLER_HPP
#define _PARA_CONTROLLER_HPP
// ### ParaController.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 4/1/2005: Last Modified
//
// ###
#include <iostream>
#include <stack>
#include "hypergraph/parallel/hypergraph.hpp"
#include "data_structures/dynamic_array.hpp"
#include "coarseners/parallel/first_choice_coarsener.hpp"
#include "coarseners/parallel/model_coarsener_2d.hpp"
#include "coarseners/parallel/approximate_first_choice_coarsener.hpp"
#include "coarseners/parallel/restrictive_first_choice_coarsening.hpp"
#include "refiners/parallel/k_way_greedy_refiner.hpp"
#include "internal/serial_controller.hpp"

namespace parkway {
namespace parallel {

namespace ds = parkway::data_structures;

class controller : public parkway::global_communicator {
 protected:
  /* Parallel partitioning options */
  int random_shuffle_before_refine_;
  int shuffled_;
  int number_of_runs_;
  int total_number_of_parts_;
  int maximum_part_weight_;
  int number_of_orig_local_vertices_;

  double keep_partitions_within_;
  double reduction_in_keep_threshold_;

  /* partition used to assign vertices to processors */
  ds::dynamic_array<int> shuffle_partition_;

  /* Approx coarsening and refinement options */
  int start_percentile_;
  int percentile_increment_;
  int approximate_refine_;

  /* Other options */
  int write_partition_to_file_;
  int display_option_;

  double balance_constraint_;

  /* auxiliary variables */
  int best_cutsize_;
  int worst_cutsize_;
  int total_cutsize_;

  double average_cutsize_;
  double start_time_;
  double total_coarsening_time_;
  double total_serial_time_;
  double total_refinement_time_;
  double total_time_;
  double accumulator_;

  ds::dynamic_array<int> best_partition_;
  ds::dynamic_array<int> map_to_orig_vertices_;

  parallel::hypergraph *hypergraph_;

  coarsener &coarsener_;
  refiner &refiner_;
  serial::controller &serial_controller_;
  std::stack<parallel::hypergraph *> hypergraphs_;

 public:
  controller(coarsener &c, refiner &r,
             serial::controller &ref, int rank, int nP, int percentile,
             int inc, int approxRefine);

  virtual ~controller();
  virtual void display_options() const = 0;
  virtual void reset_structures() = 0;
  virtual void run(MPI_Comm comm) = 0;
  virtual void set_weight_constraints(MPI_Comm comm);

  inline int number_of_runs() const { return number_of_runs_; }
  inline int number_of_parts() const { return total_number_of_parts_; }
  inline int maximum_part_weight() const { return maximum_part_weight_; }
  inline int best_cut_size() const { return best_cutsize_; }

  inline double balance_constraint() const { return balance_constraint_; }

  inline void set_number_of_runs(int n) { number_of_runs_ = n; }
  inline void set_number_of_parts(int n) { total_number_of_parts_ = n; }
  inline void set_balance_constraint(double b) { balance_constraint_ = b; }
  inline void set_display_option(int dO) { display_option_ = dO; }
  inline void set_reduction_in_keep_threshold(double r) { reduction_in_keep_threshold_ = r; }
  inline void set_kt_factor(double kT) { keep_partitions_within_ = kT; }
  inline void set_shuffle_vertices(int s) { shuffled_ = s; }
  inline void set_random_shuffle_before_refine(int s) { random_shuffle_before_refine_ = s; }
  inline void set_hypergraph(parallel::hypergraph *graph) {
    hypergraph_ = graph;
    number_of_orig_local_vertices_ = hypergraph_->number_of_vertices();
  }

  void initialize_map_to_orig_verts();
  void set_prescribed_partition(const char *filename, MPI_Comm comm);
  void store_best_partition(int numV, const dynamic_array<int> array, MPI_Comm comm);
  void partition_to_file(const char *filename, MPI_Comm comm) const;
  void copy_out_partition(int numVertices, int *pVector) const;
  void display_partition_info(MPI_Comm comm) const;
};

}  // namespace parallel
}  // namespace parkway

#endif
