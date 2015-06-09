
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

#include "ParaCoarsener.hpp"
#include "data_structures/match_request_table.hpp"
#include <iostream>

namespace ds = parkway::data_structures;

class Para2DModelCoarsener : public ParaCoarsener {
 protected:
  int vertexVisitOrder;
  int matchRequestVisitOrder;
  int divByCluWt;
  int divByHedgeLen;
  int limitOnIndexDuringCoarsening;

  ds::match_request_table *table;

 public:
  Para2DModelCoarsener(int rank, int nProcs, int nParts, int vertVisOrder,
                       int matchReqOrder, int divByWt, int divByLen,
                       std::ostream &out);
  ~Para2DModelCoarsener();

  void dispCoarseningOptions() const;
  void buildAuxiliaryStructs(int numPins, double aveVertDeg,
                             double aveHedgeSize);
  void releaseMemory();

  ParaHypergraph *coarsen(ParaHypergraph &h, MPI_Comm comm);
  ParaHypergraph *ParaFCCoarsen(ParaHypergraph &h, MPI_Comm comm);
  ParaHypergraph *ParaHedgeCoarsen(ParaHypergraph &h, MPI_Comm comm);

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
