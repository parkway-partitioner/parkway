
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

using parkway::data_structures::dynamic_array;

class hypergraph {
 public:
  hypergraph(int *vWts, int numV);
  hypergraph(int *vWts, int *pVector, int numV, int cut);
  ~hypergraph();

  void load_from_file(const char *filename);
  void buildVtoHedges();

  void reset_match_vector();
  void reset_partition_vector();
  void reset_vertex_maps();

  void project_partitions(const hypergraph &coarseGraph);
  void remove_bad_partitions(double fractionOK);
  void set_number_of_partitions(int nPartitions);
  void copy_out_partition(int *pVector, int numV, int pNo) const;
  void copy_in_partition(const int *pVector, int numV, int pNo, int cut);
  void print_characteristics(std::ostream &o);
  void print_percentiles(std::ostream &o);

  int keep_best_partition();

  inline int number_of_vertices() const { return number_of_vertices_; }
  inline int number_of_hyperedges() const { return number_of_hyperedges_; }
  inline int number_of_pins() const { return number_of_pins_; }
  inline int total_weight() const { return total_weight_; }
  inline int number_of_partitions() const { return number_of_partitions_; }
  inline int cut(int pNo) const { return partition_cuts_[pNo]; }

  inline int *vertex_weights() const { return vertex_weights_.data(); }
  inline int *hyperedge_weights() const { return hyperedge_weights_.data(); }
  inline int *match_vector() const { return match_vector_.data(); }
  inline int *pin_list() const { return pin_list_.data(); }
  inline int *hyperedge_offsets() const { return hyperedge_offsets_.data(); }
  inline int *vertex_to_hyperedges() const { return vertex_to_hyperedges_.data(); }
  inline int *vertex_offsets() const { return vertex_offsets_.data(); }

  inline int *partition_vector() const { return partition_vector_.data(); }
  inline int *partition_offsets() const {
    return partition_vector_offsets_.data();
  }
  inline int *partition_cuts() const { return partition_cuts_.data(); }
  inline int *partition_vector(int pNo) const {
    return (&partition_vector_[partition_vector_offsets_[pNo]]);
  }

  inline void set_number_of_hypererges(int newNum) { number_of_hyperedges_ = newNum; }
  inline void set_number_of_pins(int newNum) { number_of_pins_ = newNum; }
  inline void set_number_of_vertices(int newNum) { number_of_vertices_ = newNum; }
  inline void set_total_weight(int newWt) { total_weight_ = newWt; }
  inline void set_weights(int *array, int len) {
   vertex_weights_.set_data(array, len);
  }
  inline void set_hyperedge_weights(int *array, int len) {
   hyperedge_weights_.set_data(array, len);
  }
  inline void set_pin_list(int *array, int len) {
   pin_list_.set_data(array, len);
  }
  inline void set_hyperedge_offsets(int *array, int len) {
   hyperedge_offsets_.set_data(array, len);
  }
  inline void set_vertex_to_hyperedges(int *array, int len) {
   vertex_to_hyperedges_.set_data(array, len);
  }
  inline void set_vertex_offsets(int *array, int len) {
   vertex_offsets_.set_data(array, len);
  }

  inline void set_partition_cuts(int *a, int len) {
   partition_cuts_.set_data(a, len);
  }
  inline void set_partition_vector(int *a, int len) {
   partition_vector_.set_data(a, len);
  }

  int export_hyperedge_weight() const;
  int cut_size(int nP, int partitionNo) const;
  int sum_of_external_degrees(int nP, int partitionNo) const;

  void initialize_cut_sizes(int numParts);
  void check_partitions(int nP, int maxWt) const;
  void check_partition(int partitionNum, int numParts, int maxPartWt) const;

  protected:
  int total_weight_;
  int number_of_vertices_;
  int number_of_hyperedges_;
  int number_of_pins_;
  int number_of_partitions_;

  dynamic_array<int> vertex_weights_;
  dynamic_array<int> hyperedge_weights_;
  dynamic_array<int> match_vector_;

  dynamic_array<int> partition_cuts_;
  dynamic_array<int> partition_vector_;
  dynamic_array<int> partition_vector_offsets_;

  dynamic_array<int> pin_list_;
  dynamic_array<int> hyperedge_offsets_;

  dynamic_array<int> vertex_to_hyperedges_;
  dynamic_array<int> vertex_offsets_;

 void convert_to_DOMACS_graph_file(const char *fName) const;

 void check_part_weights_are_less_than(int *part_weights, const int number,
                                       int maximum) const;
};

#endif
