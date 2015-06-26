#ifndef _PARA_FCCOARSENER_HPP
#define _PARA_FCCOARSENER_HPP
// ### ParaFCCoarsener.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 31/12/2004: Last Modified
//
// ###

#include <iostream>
#include "coarseners/parallel/coarsener.hpp"
#include "data_structures/match_request_table.hpp"
#include "hypergraph/parallel/hypergraph.hpp"
#include "types.hpp"

namespace parkway {
namespace parallel {
namespace ds = parkway::data_structures;

class first_choice_coarsener : public coarsener {
 protected:
  parkway::visit_order_t vertex_visit_order_;
  parkway::visit_order_t match_request_visit_order_;
  int divide_by_cluster_weight_;
  int divide_by_hyperedge_length_;
  int limit_on_index_during_coarsening_;

  ds::match_request_table *table_;

 public:
  first_choice_coarsener(int rank, int nProcs, int nParts, int vertVisOrder, int
                         matchReqOrder, int divByWt, int divByLen);
  ~first_choice_coarsener();

  void display_options() const;
  void build_auxiliary_structures(int numPins, double aveVertDeg,
                                  double aveHedgeSize);
  void release_memory();

  hypergraph *coarsen(hypergraph &h, MPI_Comm comm);

  void set_request_arrays(int highToLow);
  void set_reply_arrays(int highToLow, int maxVertexWt);
  void process_request_replies();
  void permute_vertices_arrays(dynamic_array<int> &verts, int numLocVerts);
  void set_cluster_indices(MPI_Comm comm);

  int accept(int _locV, int _nonLocCluWt, int hToLow, int _maxWt);

  inline void set_vertex_visit_order(int vO) {
    vertex_visit_order_ = static_cast<parkway::visit_order_t>(vO - 1);
  }

  inline void set_match_request_visit_order(int mvO) {
    match_request_visit_order_ = static_cast<parkway::visit_order_t>(mvO - 1);
  }

  inline void set_divide_by_cluster_weight(int divBy) {
    divide_by_cluster_weight_ = divBy;
  }
};

}  // namespace parallel
}  // namespace parkway

#endif
