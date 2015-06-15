#ifndef _SEQ_CONTROLLER_HPP
#define _SEQ_CONTROLLER_HPP

// ### SeqController.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 30/11/2004: Last Modified
//
// ###

#include "hypergraph/parallel/hypergraph.hpp"
#include "hypergraph/serial/hypergraph.hpp"

namespace parallel = parkway::parallel;
namespace serial = parkway::serial;

using namespace std;

class sequential_controller {
protected:
  ostream &out_stream_;
  int display_option_;
  int rank_;
  int number_of_processors_;
  int number_of_parts_;
  int number_of_runs_;
  int maximum_vertex_weight_;

  double k_way_constraint_;
  double accept_proportion_;

  serial::hypergraph *hypergraph_;

  dynamic_array<int> partition_vector_;
  dynamic_array<int> partition_vector_cuts_;
  dynamic_array<int> partition_vector_offsets_;

public:
  sequential_controller(int rank, int nProcs, int nParts, ostream &out);

  virtual ~sequential_controller();
  virtual void display_options() const = 0;
  virtual void run(parallel::hypergraph &hgraph, MPI_Comm comm) = 0;
  virtual void initialize_sequential_partitions(parallel::hypergraph &h,
                                                MPI_Comm comm);
  virtual void initialize_coarsest_hypergraph(parallel::hypergraph &hgraph,
                                              MPI_Comm comm);

  int choose_best_partition() const;
  int accept_cut() const;

  inline void set_number_of_runs(int r) { number_of_runs_ = r; }
  inline void set_display_option(int d) { display_option_ = d; }
  inline void set_maximum_vertex_weight(int max) { maximum_vertex_weight_ = max; }
  inline void set_hypergraph(serial::hypergraph *hGraph) { hypergraph_ = hGraph; }
  inline void set_k_way_constraint(double c) { k_way_constraint_ = c; }
  inline void set_accept_proportion(double p) { accept_proportion_ = p; }
};

#endif
