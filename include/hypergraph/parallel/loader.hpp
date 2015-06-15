
#ifndef _PARA_HYPO_LOADER_HPP
#define _PARA_HYPO_LOADER_HPP

// ### ParaHypergraphLoader.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 4/1/2005: Last Modified
//
// ###

#include "hypergraph/parallel/hypergraph.hpp"
#include <iostream>
#include "data_structures/bit_field.hpp"
#include "data_structures/dynamic_array.hpp"

using parkway::data_structures::dynamic_array;
using parkway::data_structures::bit_field;

namespace parkway {
namespace parallel {

class loader : public global_communicator {
 protected:
  std::ostream &out_stream;

  /* hypergraph variables */

  int number_of_parts_;
  int number_of_hyperedges_;
  int number_of_local_pins_;
  int number_of_local_vertices_;
  int number_of_vertices_;
  int minimum_vertex_index_;
  int maximum_vertex_index_;
  int local_vertex_weight_;

  int number_of_allocated_hyperedges_;
  int display_options_;

  /* member for approx coarsening and refinement */

  int percentile_;

  int *vertex_weights_;
  int *match_vector_;

  dynamic_array<int> hyperedge_weights_;
  dynamic_array<int> hyperedge_offsets_;
  dynamic_array<int> local_pin_list_;

  dynamic_array<int> vertex_to_hyperedges_offset_;
  dynamic_array<int> vertex_to_hyperedges_;

  dynamic_array<int> allocated_hyperedges_;

 public:
  loader(int rank, int nProcs, int nParts, std::ostream &o);

  virtual ~loader();
  virtual void release_memory() = 0;
  virtual void load(const hypergraph &h, MPI_Comm comm) = 0;

  void compute_hyperedges_to_load(bit_field &toLoad, int numH, int numLocalPins,
                                  int *hEdgeWts, int *hEdgeOffsets,
                                  MPI_Comm comm);

  inline int number_of_parts() const { return number_of_parts_; }
  inline int percentile() const { return percentile_; }

  inline void set_display_option(int d) { display_options_ = d; }
  inline void set_percentile(int p) { percentile_ = p; }
};

}  // namespace parallel
}  // namespace parkway

#endif
