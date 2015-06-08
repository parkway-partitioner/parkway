
#ifndef _PARA_RESTR_COARSENER_HPP
#define _PARA_RESTR_COARSENER_HPP

// ### ParaRestrCoarsener.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 4/1/2005: Last Modified
//
// ###

#include "ParaHypergraphLoader.hpp"
#include "data_structures/HashTables.hpp"

using namespace std;

class ParaRestrCoarsener : public ParaHypergraphLoader {

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

  /* partition data structures */

  int numPartitions;
  int *partitionVector;
  int *partitionVectorOffsets;
  int *partitionCuts;

  DynamicArray<int> *clusterWeights;
  DynamicArray<int> *pVector;

public:
  ParaRestrCoarsener(int _rank, int _numProcs, int _numParts, ostream &out);

  virtual ~ParaRestrCoarsener();
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
  inline void setTotGraphWt(int tot) { totalHypergraphWt = tot; }
};

#endif
