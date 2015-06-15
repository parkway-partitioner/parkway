#ifndef _KWAY_GREEDY_REFINER_HPP
#define _KWAY_GREEDY_REFINER_HPP

// ### GreedyKwayRefiner.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 25/4/2004: Last Modified
//
// ###

#include "refiners/serial/refiner.hpp"
#include "hypergraph/VertexNode.hpp"
#include "hypergraph/serial/hypergraph.hpp"

namespace parkway {
namespace serial {
namespace ds = parkway::data_structures;

class greedy_k_way_refiner : public refiner {
protected:
  int number_of_non_positive_moves_;
  double limit_;

  // ###
  // data_ structures from point of view of vertices
  // ###

  ds::dynamic_array<int> number_of_neighboring_parts_;
  ds::dynamic_array<int> neighbors_of_vertex_;
  ds::dynamic_array<int> neighbors_of_vertex_offsets_;

  // ###
  // data_ structures from point of view of hyperedges
  // ###

  ds::dynamic_array<int> hyperedge_vertices_in_part_;
  ds::dynamic_array<int> hyperedge_vertices_in_part_offsets_;

  // ###
  // auxiliary structures
  // ###

  ds::dynamic_array<int> vertices_;
  ds::dynamic_array<int> vertex_seen_;
  ds::dynamic_array<int> seen_vertices_;
  ds::dynamic_array<int> parts_spanned_;

public:
  greedy_k_way_refiner(int max, int nparts, double ave, double limit, int dL);
  ~greedy_k_way_refiner();

  void display_options(std::ostream &out) const;
  void build_data_structures();
  void destroy_data_structures();
  int initialize_data_structures();
  void update_adjacent_vertex_stats(int v, int sP, int bestDP);

  void refine(hypergraph &h);
  void rebalance(hypergraph &h);

  int greedy_pass();
  int rebalancing_pass();

  inline void set_limit() {
    double a = static_cast<double>(numVertices) * limit_;
    number_of_non_positive_moves_ = static_cast<int>(floor(a));
  }

  inline int find_heaviest_overweight() const {
    int p = -1;
    for (int i = 0; i < number_of_parts_; ++i) {
      if (part_weights_[i] > maximum_part_weight_) {
        if (p == -1)
          p = i;
        else if (part_weights_[i] > part_weights_[p])
          p = i;
      }
    }

    return p;
  }
};

}  // namespace serial
}  // namespace parkway

#endif
