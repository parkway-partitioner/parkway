
#ifndef _PARA_COARSENER_HPP
#define _PARA_COARSENER_HPP

// ### ParaCoarsener.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 4/1/2005: Last Modified
//
// ###

#include "hypergraph/parallel/hypergraph.hpp"
#include "hypergraph/parallel/loader.hpp"
#include "data_structures/dynamic_array.hpp"

using parkway::data_structures::dynamic_array;
namespace parallel = parkway::hypergraph::parallel;

class parallel_coarsener : public parallel::loader {

protected:
  /* coarsening auxiliary variables */

  int total_hypergraph_weight_;
  int maximum_vertex_weight_;
  int minimum_number_of_nodes_;
  int stop_coarsening_;
  int cluster_index_;
  int total_number_of_clusters_;
  int minimum_cluster_index_;

  double reduction_ratio_;
  double balance_constrain_;

  dynamic_array<int> cluster_weights_;

public:
  parallel_coarsener(int _rank, int _numProcs, int _numParts, std::ostream &out);

  virtual ~parallel_coarsener();
  virtual parallel::hypergraph *coarsen(parallel::hypergraph &h, MPI_Comm comm) = 0;
  virtual void set_cluster_indices(MPI_Comm comm) = 0;
  virtual void release_memory() = 0;
  virtual void display_coarsening_options() const = 0;
  virtual void build_auxiliary_structures(int numTotPins, double aveVertDeg,
                                          double aveHedgeSize) = 0;

  void load(const parallel::hypergraph &h, MPI_Comm comm);

  parallel::hypergraph *constract_hyperedges(parallel::hypergraph &h, MPI_Comm comm);

  inline void set_reduction_ratio(double ratio) {
    reduction_ratio_ = ratio;
  }

  inline void set_balance_constraint(double constraint) {
    balance_constrain_ = constraint;
  }

  inline void set_minimum_number_of_nodes(int min) { minimum_number_of_nodes_ = min; }
  inline void set_maximum_vertex_weight(int m) { maximum_vertex_weight_ = m; }
  inline void set_total_hypergraph_weight(int t) { total_hypergraph_weight_ = t; }
};

#endif
