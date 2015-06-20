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
#include "coarseners/parallel/coarsener.hpp"
#include "data_structures/dynamic_array.hpp"
#include "hypergraph/parallel/hypergraph.hpp"

namespace parkway {
namespace parallel {
namespace ds = parkway::data_structures;

class restrictive_coarsening : public coarsener {
 public:
  restrictive_coarsening(int _rank, int _numProcs, int _numParts,
                         std::ostream &out);

  virtual ~restrictive_coarsening();
  // virtual hypergraph *coarsen(hypergraph &h, MPI_Comm comm) = 0;
  // virtual void set_cluster_indices(MPI_Comm comm) = 0;
  // virtual void release_memory() = 0;
  // virtual void display_options() const = 0;
  // virtual void build_auxiliary_structures(int numTotPins, double aveVertDeg,
  //                                         double aveHedgeSize) = 0;

  void load(const hypergraph &h, MPI_Comm comm);

  hypergraph *contract_hyperedges(hypergraph &h, MPI_Comm comm);

  inline void set_minimum_nodes(int min) {
    minimum_nodes_ = min;
  }

  inline void set_total_graph_weight(int tot) {
    total_hypergraph_weight_ = tot;
  }

  void initialize_vertex_to_hyperedges();
  void load_non_local_hyperedges();
  void prepare_data_to_send(int n_local_hyperedges, int n_local_pins,
                            ds::dynamic_array<int> local_hyperedge_weights,
                            ds::dynamic_array<int> local_hyperedge_offsets,
                            ds::dynamic_array<int> local_pins, MPI_Comm comm);

 protected:
  /* coarsening auxiliary variables */
  int minimum_nodes_;
  int total_clusters_;

  /* partition data_ structures */
  ds::dynamic_array<int> partition_vector_;
  ds::dynamic_array<int> partition_vector_offsets_;
  ds::dynamic_array<int> partition_cuts_;
  ds::dynamic_array<int> part_vector_;
};

}  // namespace parallel
}  // namespace parkway

#endif
