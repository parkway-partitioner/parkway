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
#include <iostream>
#include "coarseners/parallel/approximate_coarsener.hpp"
#include "data_structures/match_request_table.hpp"
#include "hypergraph/parallel/hypergraph.hpp"

namespace parkway {
namespace parallel {
namespace ds = parkway::data_structures;

class approximate_first_choice_coarsener : public approximate_coarsener {
 protected:
  int vertexVisitOrder;
  int matchRequestVisitOrder;
  int divByCluWt;
  int divByHedgeLen;
  int limitOnIndexDuringCoarsening;
  ds::match_request_table *table;

 public:
  approximate_first_choice_coarsener(int rank, int nProcs, int nParts,
                                     int percentile, int inc, int vertVisOrder,
                                     int matchReqOrder, int divByWt,
                                     int divByLen, std::ostream &out);
  ~approximate_first_choice_coarsener();

  void dispCoarseningOptions() const;
  void buildAuxiliaryStructs(int numPins, double aveVertDeg,
                             double aveHedgeSize);
  void release_memory();

  parallel::hypergraph *coarsen(parallel::hypergraph &h, MPI_Comm comm);

  void setRequestArrays(int highToLow);
  void setReplyArrays(int highToLow, int maxVertexWt);
  void processReqReplies();
  void permuteVerticesArray(int *verts, int numLocVerts);
  void setClusterIndices(MPI_Comm comm);

  int accept(int _locV, int _nonLocCluWt, int hToLow, int _maxWt);

  void printVisitOrder(int variable) const;

  inline void setVertexVisitOrder(int vO) {
    vertexVisitOrder = vO;
  }

  inline void setMatchRequestVisitOrder(int mvO) {
    matchRequestVisitOrder = mvO;
  }

  inline void setDivByCluWt(int divBy) {
    divByCluWt = divBy;
  }
};

}  // namespace parallel
}  // namespace parkway

#endif
