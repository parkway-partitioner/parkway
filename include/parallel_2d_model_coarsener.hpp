
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

#include "parallel_coarsener.hpp"
#include "data_structures/match_request_table.hpp"
#include "hypergraph/parallel/hypergraph.hpp"
#include <iostream>

namespace parallel = parkway::parallel;
namespace ds = parkway::data_structures;

class parallel_2d_model_coarsener : public parallel_coarsener {
 protected:
  int vertex_visit_order_;
  int match_request_visit_order_;
  int divide_by_cluster_weight_;
  int divide_by_hyperedge_length_;
  int limit_on_index_during_coarsening_;

  ds::match_request_table *table_;

 public:
  parallel_2d_model_coarsener(int rank, int nProcs, int nParts, int vertVisOrder,
                       int matchReqOrder, int divByWt, int divByLen,
                       std::ostream &out);
  ~parallel_2d_model_coarsener();

  void display_options() const;
  void build_auxiliary_structures(int numPins, double aveVertDeg,
                                  double aveHedgeSize);
  void release_memory();

  parallel::hypergraph *coarsen(parallel::hypergraph &h, MPI_Comm comm);
  parallel::hypergraph *parallel_first_choice_coarsen(parallel::hypergraph &h, MPI_Comm comm);
  parallel::hypergraph *parallel_hyperedge_coarsen(parallel::hypergraph &h, MPI_Comm comm);

  void set_request_arrays(int highToLow);
  void set_reply_arrays(int highToLow, int maxVertexWt);
  void process_request_replies();
  void permute_vertices_arrays(int *verts, int numLocVerts);
  void set_cluster_indices(MPI_Comm comm);

  int accept(int _locV, int _nonLocCluWt, int hToLow, int _maxWt);

  void print_visit_order(int variable) const;

  inline void set_vertex_visit_order(int vO) { vertex_visit_order_ = vO; }
  inline void set_match_request_visit_order(int mvO) {
    match_request_visit_order_ = mvO;
  }
  inline void set_divide_by_cluster_weight(int divBy) { divide_by_cluster_weight_ = divBy; }
};

#endif
