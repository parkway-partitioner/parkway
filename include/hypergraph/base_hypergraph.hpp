#ifndef _HYPERGRAPH_BASE_HYPERGRAPH_HPP
#define _HYPERGRAPH_BASE_HYPERGRAPH_HPP

#include "data_structures/dynamic_array.hpp"

using parkway::data_structures::dynamic_array;

namespace parkway {
namespace hypergraph {

class base_hypergraph {
 public:
  base_hypergraph(int number_of_vertices = 0, int number_of_partitions = 0);
  virtual ~base_hypergraph();

  inline int number_of_partitions() const {
    return number_of_partitions_;
  }

  virtual inline void set_number_of_partitions(int number) = 0;


  inline int number_of_vertices() const {
    return number_of_vertices_;
  }

  inline void set_number_of_vertices(int number) {
    number_of_vertices_ = number;
  }


  inline int number_of_hyperedges() const {
    return number_of_hyperedges_;
  }

  inline void set_number_of_hyperedges(int number) {
    number_of_hyperedges_ = number;
  }


  inline int number_of_pins() const {
    return number_of_pins_;
  }

  inline void set_number_of_pins(int number) {
    number_of_pins_ = number;
  }


  inline dynamic_array<int> vertex_weights() const {
    return vertex_weights_;
  }

  // inline void set_vertex_weights(int *weights, int length) {
  //   vertex_weights_.set_data(weights, length);
  // }

  inline void set_vertex_weights(dynamic_array<int> weights) {
    vertex_weights_ = weights;
  }


  inline dynamic_array<int> hyperedge_weights() const {
    return hyperedge_weights_;
  }

  // inline void set_hyperedge_weights(int *array, int len) {
  //   hyperedge_weights_.set_data(array, len);
  // }

  inline void set_hyperedge_weights(dynamic_array<int> weights) {
    hyperedge_weights_ = weights;
  }


  inline dynamic_array<int> partition_vector() const {
    return partition_vector_;
  }

  // inline void set_partition_vector(int *array, int length) {
  //   partition_vector_.set_data(array, length);
  // }

  inline void set_partition_vector(dynamic_array<int> partition_vector) {
    partition_vector_ = partition_vector;
  }


  inline dynamic_array<int> partition_offsets() const {
    return partition_vector_offsets_;
  }


  inline dynamic_array<int> partition_cuts() const {
    return partition_cuts_;
  }

  // inline void set_partition_cuts(int *array, int len) {
  //   partition_cuts_.set_data(array, len);
  // }

  inline void set_partition_cuts(dynamic_array<int> cuts) {
    partition_cuts_ = cuts;
  }


  inline dynamic_array<int> match_vector() const {
    return match_vector_;
  }


  inline dynamic_array<int> pin_list() const {
    return pin_list_;
  }

  // inline void set_pin_list(int *array, int len) {
  //   pin_list_.set_data(array, len);
  // }

  inline void set_pin_list(dynamic_array<int> pin_list) {
    pin_list_ = pin_list;
  }


  inline dynamic_array<int> hyperedge_offsets() const {
    return hyperedge_offsets_;
  }

  // inline void set_hyperedge_offsets(int *array, int len) {
  //   hyperedge_offsets_.set_data(array, len);
  // }

  inline void set_hyperedge_offsets(dynamic_array<int> offsets) {
    hyperedge_offsets_ = offsets;
  }

 protected:
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
};

}  // namespace hypergraph
}  // namespace parkway


#endif  // _HYPERGRAPH_BASE_HYPERGRAPH_HPP
