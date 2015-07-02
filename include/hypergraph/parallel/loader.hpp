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

namespace parkway {
namespace parallel {
namespace ds = parkway::data_structures;

class loader : public global_communicator {
 public:
  loader(int rank, int number_of_processors, int number_of_parts,
         std::ostream &out_stream, int display_option = 0);

  virtual ~loader();

  virtual void release_memory() = 0;

  virtual void load(const hypergraph &h, MPI_Comm comm) = 0;

  void compute_hyperedges_to_load(ds::bit_field &to_load,
                                  int number_of_hyperedges,
                                  int num_local_pins,
                                  dynamic_array<int> &hyperedge_weights,
                                  dynamic_array<int> &hyperedge_offsets,
                                  MPI_Comm comm);

  inline int number_of_parts() const {
    return number_of_parts_;
  }

  inline int percentile() const {
    return percentile_;
  }

  inline void set_percentile(int p) {
    percentile_ = p;
  }

  inline void set_display_option(int d) {
    display_options_ = d;
  }

  inline int display_option() const {
    return display_options_;
  }

 protected:
  // Output stream
  std::ostream &out_stream;

  // Hypergraph variables
  int number_of_parts_;
  int number_of_hyperedges_;
  int number_of_local_pins_;
  int number_of_local_vertices_;
  int number_of_vertices_;
  int minimum_vertex_index_;
  int maximum_vertex_index_;
  int local_vertex_weight_;

  // Misc
  int number_of_allocated_hyperedges_;
  int display_options_;

  // Member for approximate coarsening and refinement
  int percentile_;

  ds::dynamic_array<int> vertex_weights_;
  ds::dynamic_array<int> match_vector_;
  ds::dynamic_array<int> hyperedge_weights_;
  ds::dynamic_array<int> hyperedge_offsets_;
  ds::dynamic_array<int> local_pin_list_;
  ds::dynamic_array<int> vertex_to_hyperedges_offset_;
  ds::dynamic_array<int> vertex_to_hyperedges_;
  ds::dynamic_array<int> allocated_hyperedges_;
};

}  // namespace parallel
}  // namespace parkway

#endif
