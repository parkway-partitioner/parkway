
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

#include "ParaCoarsener.hpp"
#include <iostream>
#include "data_structures/match_request_table.hpp"

namespace ds = parkway::data_structures;

class ParaFCCoarsener : public ParaCoarsener {
 protected:
  int vertexVisitOrder;
  int matchRequestVisitOrder;
  int divByCluWt;
  int divByHedgeLen;
  int limitOnIndexDuringCoarsening;

  ds::match_request_table *table;

 public:
  ParaFCCoarsener(int rank, int nProcs, int nParts, int vertVisOrder,
                  int matchReqOrder, int divByWt, int divByLen, std::ostream &out);
  ~ParaFCCoarsener();

  void dispCoarseningOptions() const;
  void buildAuxiliaryStructs(int numPins, double aveVertDeg,
                             double aveHedgeSize);
  void releaseMemory();

  ParaHypergraph *coarsen(ParaHypergraph &h, MPI_Comm comm);

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
