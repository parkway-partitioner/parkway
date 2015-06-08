
#ifndef _PARA_REFINER_HPP
#define _PARA_REFINER_HPP

// ### ParaRefiner.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 4/1/2005: Last Modified
//
// ###

#include "ParaHypergraphLoader.hpp"

using namespace std;

class ParaRefiner : public ParaHypergraphLoader {
protected:
  int maxPartWt;
  int numPartitions;

  int *partitionVector;
  int *partitionVectorOffsets;
  int *partitionCuts;

  int *currPVector;
  int currPnumber;

  double balConstraint;
  double avePartWt;

  DynamicArray<int> partWeights;

  /* newly added structures */

  int numNonLocVerts;
  int numNonLocVertsHedges;
  int *currNonLocPVector;

  DynamicArray<int> nonLocVerts;
  DynamicArray<int> partIndices;
  DynamicArray<int> indexIntoPartIndices;

  DynamicArray<int> nonLocVToHedges;
  DynamicArray<int> nonLocOffsets;

  MapToPosInt toNonLocVerts;

public:
  ParaRefiner(int rank, int nProcs, int nParts, ostream &out);

  virtual ~ParaRefiner();
  virtual void dispRefinementOptions() const = 0;
  virtual void releaseMemory() = 0;
  virtual void initDataStructs(const ParaHypergraph &h, MPI_Comm comm) = 0;
  virtual void resetDataStructs() = 0;
  virtual void setPartitioningStructs(int pNumber, MPI_Comm comm) = 0;
  virtual void refine(ParaHypergraph &h, MPI_Comm comm) = 0;
  virtual int computeCutsize(MPI_Comm comm) = 0;

  void loadHyperGraph(const ParaHypergraph &h, MPI_Comm comm);
  void initPartitionStructs(const ParaHypergraph &h, MPI_Comm comm);

  inline void setBalConstraint(double b) { balConstraint = b; }
  inline int getMaxPartWt() const { return maxPartWt; }
};

#endif
