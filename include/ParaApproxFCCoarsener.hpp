
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

#include "ParaApproxCoarsener.hpp"

using namespace std;

class ParaApproxFCCoarsener : public ParaApproxCoarsener {

protected:
  int vertexVisitOrder;
  int matchRequestVisitOrder;
  int divByCluWt;
  int divByHedgeLen;
  int limitOnIndexDuringCoarsening;

  MatchRequestTable *table;
  // ConnVertTable *connTable;

public:
  ParaApproxFCCoarsener(int rank, int nProcs, int nParts, int percentile,
                        int inc, int vertVisOrder, int matchReqOrder,
                        int divByWt, int divByLen, ostream &out);
  ~ParaApproxFCCoarsener();

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

  inline void setVertexVisitOrder(register int vO) { vertexVisitOrder = vO; }
  inline void setMatchRequestVisitOrder(register int mvO) {
    matchRequestVisitOrder = mvO;
  }
  inline void setDivByCluWt(register int divBy) { divByCluWt = divBy; }
};

#endif
