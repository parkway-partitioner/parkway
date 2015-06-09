
#ifndef _PARA_COARSENER_HPP
#define _PARA_COARSENER_HPP

// ### ParaCoarsener.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 4/1/2005: Last Modified
//
// ###

#include "ParaHypergraphLoader.hpp"
#include "data_structures/dynamic_array.hpp"

using parkway::data_structures::dynamic_array;

class ParaCoarsener : public ParaHypergraphLoader {

protected:
  /* coarsening auxiliary variables */

  int totalHypergraphWt;
  int maxVertexWt;
  int minNodes;
  int stopCoarsening;
  int clusterIndex;
  int totalClusters;
  int myMinCluIndex;

  double reductionRatio;
  double balConstraint;

  dynamic_array<int> clusterWeights;

public:
  ParaCoarsener(int _rank, int _numProcs, int _numParts, std::ostream &out);

  virtual ~ParaCoarsener();
  virtual ParaHypergraph *coarsen(ParaHypergraph &h, MPI_Comm comm) = 0;
  virtual void setClusterIndices(MPI_Comm comm) = 0;
  virtual void releaseMemory() = 0;
  virtual void dispCoarseningOptions() const = 0;
  virtual void buildAuxiliaryStructs(int numTotPins, double aveVertDeg,
                                     double aveHedgeSize) = 0;

  void loadHyperGraph(const ParaHypergraph &h, MPI_Comm comm);

  ParaHypergraph *contractHyperedges(ParaHypergraph &h, MPI_Comm comm);

  inline void setReductionRatio(double ratio) {
    reductionRatio = ratio;
  }
  inline void setBalConstraint(double constraint) {
    balConstraint = constraint;
  }
  inline void setMinNodes(int min) { minNodes = min; }
  inline void setMaxVertexWt(int m) { maxVertexWt = m; }
  inline void setTotGraphWt(int t) { totalHypergraphWt = t; }
};

#endif
