
#ifndef _PARA_RESTR_FCC_HPP
#define _PARA_RESTR_FCC_HPP

// ### ParaRestrFCCoarsener.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 31/12/2004: Last Modified
//
// ###

#include "ParaRestrCoarsener.hpp"

using namespace std;

class ParaRestrFCCoarsener : public ParaRestrCoarsener {

protected:
  int vertexVisitOrder;
  int divByCluWt;
  int divByHedgeLen;
  int limitOnIndexDuringCoarsening;

public:
  ParaRestrFCCoarsener(int rank, int nProcs, int nParts, int verVisOrder,
                       int divByWt, int divByLen, ostream &out);
  ~ParaRestrFCCoarsener();

  void dispCoarseningOptions() const;
  void buildAuxiliaryStructs(int numPins, double aveVertDeg,
                             double aveHedgeSize);
  void releaseMemory();

  parallel_hypergraph *coarsen(parallel_hypergraph &h, MPI_Comm comm);

  void permuteVerticesArray(int *verts, int numLocVerts);
  void setClusterIndices(MPI_Comm comm);
  void printVisitOrder(int variable) const;
};

#endif
