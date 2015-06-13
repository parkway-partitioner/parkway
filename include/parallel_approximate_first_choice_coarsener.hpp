
#ifndef _PARA_APPROX_FCCOARSENER_HPP
#define _PARA_APPROX_FCCOARSENER_HPP

// ### ParaApproxFCCoarsener.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 31/12/2004: Last Modified
//
// ###

#include "parallel_approximate_coarsener.hpp"
#include <iostream>
#include "data_structures/match_request_table.hpp"
#include "hypergraph/parallel/hypergraph.hpp"

using parkway::hypergraph::parallel::hypergraph;

class parallel_approximate_first_choice_coarsener : public parallel_approximate_coarsener {
 protected:
  int vertexVisitOrder;
  int matchRequestVisitOrder;
  int divByCluWt;
  int divByHedgeLen;
  int limitOnIndexDuringCoarsening;

  parkway::data_structures::match_request_table *table;

 public:
  parallel_approximate_first_choice_coarsener(int rank, int nProcs, int nParts, int percentile,
                        int inc, int vertVisOrder, int matchReqOrder,
                        int divByWt, int divByLen, std::ostream &out);
  ~parallel_approximate_first_choice_coarsener();

  void dispCoarseningOptions() const;
  void buildAuxiliaryStructs(int numPins, double aveVertDeg,
                             double aveHedgeSize);
  void release_memory();

  hypergraph *coarsen(hypergraph &h, MPI_Comm comm);

  void setRequestArrays(int highToLow);
  void setReplyArrays(int highToLow, int maxVertexWt);
  void processReqReplies();
  void permuteVerticesArray(int *verts, int numLocVerts);
  void setClusterIndices(MPI_Comm comm);

  int accept(int _locV, int _nonLocCluWt, int hToLow, int _maxWt);

  void printVisitOrder(int variable) const;

  inline void setVertexVisitOrder(int vO) { vertexVisitOrder = vO; }
  inline void setMatchRequestVisitOrder(int mvO) {
    matchRequestVisitOrder = mvO;
  }
  inline void setDivByCluWt(int divBy) { divByCluWt = divBy; }
};

#endif
