#ifndef _PARA_RESTR_COARSENER_HPP
#define _PARA_RESTR_COARSENER_HPP

// ### ParaRestrCoarsener.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 4/1/2005: Last Modified
//
// ###

#include <iostream>
#include "data_structures/dynamic_array.hpp"
#include "hypergraph/parallel/hypergraph.hpp"
#include "hypergraph/parallel/loader.hpp"

namespace parallel = parkway::parallel;
namespace ds = parkway::data_structures;

class parallel_restrictive_coarsening : public parallel::loader {
 protected:
  /* coarsening auxiliary variables */
  int total_hypergraph_weight_;
  int maximum_vertex_weight_;
  int minimum_nodes_;
  int stop_coarsening_;
  int cluster_index_;
  int total_clusters_;
  int minimum_cluster_index_;

  double reduction_ratio_;
  double balance_constraint_;

  /* partition data_ structures */

  int *partition_vector_;
  int *partition_vector_offsets_;
  int *partition_cuts_;

  ds::dynamic_array<int> *cluster_weights_;
  ds::dynamic_array<int> *oartition_vector_;

 public:
  parallel_restrictive_coarsening(int _rank, int _numProcs, int _numParts, std::ostream &out);

  virtual ~parallel_restrictive_coarsening();
  virtual parallel::hypergraph *coarsen(parallel::hypergraph &h, MPI_Comm comm) = 0;
  virtual void set_cluster_indices(MPI_Comm comm) = 0;
  virtual void release_memory() = 0;
  virtual void display_options() const = 0;
  virtual void build_auxiliary_structures(int numTotPins, double aveVertDeg,
                                          double aveHedgeSize) = 0;

  void load(const parallel::hypergraph &h, MPI_Comm comm);

  parallel::hypergraph *contract_hyperedges(parallel::hypergraph &h,
                                            MPI_Comm comm);

  inline void set_reduction_ratio(double ratio) {
    reduction_ratio_ = ratio;
  }
  inline void set_balance_constraint(double constraint) {
    balance_constraint_ = constraint;
  }
  inline void set_minimum_nodes(int min) { minimum_nodes_ = min; }
  inline void set_maximum_vertex_weight(int m) { maximum_vertex_weight_ = m; }
  inline void set_total_graph_weight(int tot) { total_hypergraph_weight_ = tot; }
};

#endif
