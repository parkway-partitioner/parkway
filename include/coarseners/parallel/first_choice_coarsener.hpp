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

namespace parkway {
namespace parallel {
namespace ds = parkway::data_structures;

class first_choice_coarsener : public coarsener {
 protected:
  int vertex_visit_order_;
  int match_request_visit_order_;
  int divide_by_cluster_weight_;
  int divide_by_hyperedge_length_;
  int limit_on_index_during_coarsening_;

  ds::match_request_table *table_;

 public:
  first_choice_coarsener(int rank, int nProcs, int nParts, int vertVisOrder,
                  int matchReqOrder, int divByWt, int divByLen, std::ostream &out);
  ~first_choice_coarsener();

  void display_options() const;
  void build_auxiliary_structures(int numPins, double aveVertDeg,
                                  double aveHedgeSize);
  void release_memory();

  hypergraph *coarsen(hypergraph &h, MPI_Comm comm);

  void set_request_arrays(int highToLow);
  void set_reply_arrays(int highToLow, int maxVertexWt);
  void process_request_replies();
  void permute_vertices_arrays(int *verts, int numLocVerts);
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
