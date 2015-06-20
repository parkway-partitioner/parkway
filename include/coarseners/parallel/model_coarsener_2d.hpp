#ifndef _PARA_2DMODEL_COARSENER_HPP
#define _PARA_2DMODEL_COARSENER_HPP
// ### Para2DModelCoarsener.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 19/04/2005: Last Modified
//
// NOTES: Developed a new algorithm specifically for use
// on 2D hypergraph decomposition of PageRank matrices.
// Algorithm behaves like a HYperedge Coarsening Algorithm
//
// ###

#include "coarseners/parallel/coarsener.hpp"
#include "data_structures/match_request_table.hpp"
#include "hypergraph/parallel/hypergraph.hpp"
#include <iostream>

namespace parkway {
namespace parallel {
namespace ds = parkway::data_structures;

class model_coarsener_2d : public coarsener {
 protected:
  int vertex_visit_order_;
  int match_request_visit_order_;
  int divide_by_cluster_weight_;
  int divide_by_hyperedge_length_;
  int limit_on_index_during_coarsening_;
  ds::match_request_table *table_;

 public:
  model_coarsener_2d(int rank, int nProcs, int nParts, int vertVisOrder,
                     int matchReqOrder, int divByWt, int divByLen,
                     std::ostream &out);
  ~model_coarsener_2d();

  void display_options() const;
  void build_auxiliary_structures(int numPins, double aveVertDeg,
                                  double aveHedgeSize);
  void release_memory();

  hypergraph *coarsen(hypergraph &h, MPI_Comm comm);
  hypergraph *parallel_first_choice_coarsen(hypergraph &h, MPI_Comm comm);
  hypergraph *parallel_hyperedge_coarsen(hypergraph &h, MPI_Comm comm);

  void set_request_arrays(int highToLow);
  void set_reply_arrays(int highToLow, int maxVertexWt);
  void process_request_replies();
  void permute_vertices_arrays(ds::dynamic_array<int> &verts, int numLocVerts);
  void set_cluster_indices(MPI_Comm comm);

  int accept(int _locV, int _nonLocCluWt, int hToLow, int _maxWt);

  void print_visit_order(int variable) const;

  inline void set_vertex_visit_order(int vO) {
    vertex_visit_order_ = vO;
  }

  inline void set_match_request_visit_order(int mvO) {
    match_request_visit_order_ = mvO;
  }

  inline void set_divide_by_cluster_weight(int divBy) {
    divide_by_cluster_weight_ = divBy;
  }
};

}  // namespace parallel
}  // namespace parkway

#endif
