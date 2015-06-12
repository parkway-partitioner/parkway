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

using parkway::data_structures::dynamic_array;

namespace parkway {
namespace hypergraph {

class serial_hypergraph : public base_hypergraph {
 public:
  serial_hypergraph(int *vWts, int numV);
  serial_hypergraph(int *vWts, int *pVector, int numV, int cut);
  ~serial_hypergraph();

  void load_from_file(const char *filename);
  void buildVtoHedges();

  void reset_match_vector();
  void reset_partition_vector();
  void reset_vertex_maps();

  void project_partitions(const serial_hypergraph &coarseGraph);
  void remove_bad_partitions(double fractionOK);
  void set_number_of_partitions(int nPartitions) override;
  void copy_out_partition(int *pVector, int numV, int pNo) const;
  void copy_in_partition(const int *pVector, int numV, int pNo, int cut);
  void print_characteristics(std::ostream &o);
  void print_percentiles(std::ostream &o);

  int keep_best_partition();

  inline int total_weight() const { return total_weight_; }
  inline int cut(int pNo) const { return partition_cuts_[pNo]; }

  inline int *vertex_to_hyperedges() const { return vertex_to_hyperedges_.data(); }
  inline int *vertex_offsets() const { return vertex_offsets_.data(); }


  inline void set_total_weight(int newWt) { total_weight_ = newWt; }

  inline void set_vertex_to_hyperedges(int *array, int len) {
   vertex_to_hyperedges_.set_data(array, len);
  }

  inline void set_vertex_offsets(int *array, int len) {
   vertex_offsets_.set_data(array, len);
  }

  int export_hyperedge_weight() const;
  int cut_size(int nP, int partitionNo) const;
  int sum_of_external_degrees(int nP, int partitionNo) const;

  void initialize_cut_sizes(int numParts);
  void check_partitions(int nP, int maxWt) const;
  void check_partition(int partitionNum, int numParts, int maxPartWt) const;

  protected:
  int total_weight_;

  dynamic_array<int> vertex_to_hyperedges_;
  dynamic_array<int> vertex_offsets_;

 void convert_to_DOMACS_graph_file(const char *fName) const;
 void check_part_weights_are_less_than(int *part_weights, const int number,
                                       int maximum) const;
};

}  // hypergraph
}  // parkway

#endif
