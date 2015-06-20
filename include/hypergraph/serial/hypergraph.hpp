#ifndef _HYPERGRAPH_HPP
#define _HYPERGRAPH_HPP

// ### Hypergraph.hpp ###
//
// Copyright (C) 2004,  Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 30/11/2004: Last Modified
//
// ###

#include <unistd.h>
#include <cstdio>
#include <cassert>
#include <iostream>

#include "Macros.h"
#include "Funct.hpp"
#include "data_structures/dynamic_array.hpp"
#include "hypergraph/base_hypergraph.hpp"

namespace parkway {
namespace serial {

namespace ds = parkway::data_structures;
namespace hg = parkway::hypergraph;

class hypergraph : public hg::base_hypergraph {
 public:
  hypergraph(ds::dynamic_array<int> vWts, int numV);
  hypergraph(ds::dynamic_array<int> vWts, ds::dynamic_array<int> pVector, int numV, int cut);
  ~hypergraph();

  void load_from_file(const char *filename);
  void buildVtoHedges();

  void reset_match_vector();
  void reset_partition_vector();
  void reset_vertex_maps();

  void project_partitions(const hypergraph &coarseGraph);
  void remove_bad_partitions(double fractionOK);
  void set_number_of_partitions(int nPartitions) override;
  void copy_out_partition(ds::dynamic_array<int> pVector, int numV, int pNo) const;
  void copy_in_partition(const ds::dynamic_array<int> pVector, int numV, int pNo, int cut);
  void print_characteristics(std::ostream &o);
  void print_percentiles(std::ostream &o);

  int keep_best_partition();

  inline int total_weight() const { return total_weight_; }
  inline int cut(int pNo) const { return partition_cuts_[pNo]; }

  inline ds::dynamic_array<int> vertex_to_hyperedges() const {
    return vertex_to_hyperedges_;
  }

  inline ds::dynamic_array<int> vertex_offsets() const {
    return vertex_offsets_;
  }

  inline void set_total_weight(int newWt) {
    total_weight_ = newWt;
  }

  inline void set_vertex_to_hyperedges(ds::dynamic_array<int> &arr) {
   vertex_to_hyperedges_ = arr;
  }

  inline void set_vertex_offsets(ds::dynamic_array<int> &offsets) {
   vertex_offsets_ = offsets;
  }

  int export_hyperedge_weight() const;
  int cut_size(int nP, int partitionNo) const;
  int sum_of_external_degrees(int nP, int partitionNo) const;

  void initialize_cut_sizes(int numParts);
  void check_partitions(int nP, int maxWt) const;
  void check_partition(int partitionNum, int numParts, int maxPartWt) const;

  protected:
  int total_weight_;

  ds::dynamic_array<int> vertex_to_hyperedges_;
  ds::dynamic_array<int> vertex_offsets_;

 void convert_to_DOMACS_graph_file(const char *fName);
 void check_part_weights_are_less_than(ds::dynamic_array<int> &part_weights,
                                       const int number, int maximum) const;
};

}  // serial
}  // parkway

#endif
