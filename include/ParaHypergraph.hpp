#ifndef _PARA_HYPERGRAPH_HPP
#define _PARA_HYPERGRAPH_HPP

// ### ParaHypergraph.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 17/4/2004: Modified hyperedge storage:
//            - removed hyperedges from hash table
//            - instead, added DynamicArray<int>
//              to represent them as a pin list
//              can be indexed via hash table with key
//              or directly via index in pin list
//
// 03/12/2004: Last Modified
//
// ###

#include <unistd.h>
#include <fstream>
#include <iostream>
#include <cstdio>
#include <cassert>

#include "GlobalCommunicator.hpp"
#include "data_structures/complete_binary_tree.hpp"
#include "data_structures/HashTables.hpp"
#include "data_structures/bit_field.hpp"

using namespace parkway::data_structures;

class ParaHypergraph : public GlobalCommunicator {
protected:
  int indexInSequence;
  int doNotCoarsen;
  int numTotalVertices;
  int numLocalVertices;
  int numLocalPins;
  int numLocalHedges;
  int minVertexIndex;
  int localVertexWt;
  int numPartitions;

  DynamicArray<int> vWeight;
  DynamicArray<int> matchVector;

  DynamicArray<int> partitionVector;
  DynamicArray<int> partitionOffsetsVector;
  DynamicArray<int> partitionCutsizesVector;

  DynamicArray<int> localPins;
  DynamicArray<int> hEdgeOffsets;
  DynamicArray<int> hEdgeWeights;

  DynamicArray<int> vToOrigV;

public:
  ParaHypergraph(int myRank, int nProcs, int _numLocVerts, int _totVerts,
                 int _minVertIndex, int coarsen, int *wtArray);
  ParaHypergraph(int myRank, int nProcs, int _numLocVerts, int _totVerts,
                 int _minVertIndex, int coarsen, int cut, int *wtArray,
                 int *partArray);
  ParaHypergraph(int myRank, int nProcs, const char *filename, int dispOption,
                 ostream &out, MPI_Comm comm);
  ParaHypergraph(int myRank, int nProcs, int numLocVerts, int numLocHedges,
                 int maxHedgeLen, const int *vWeights, const int *hEdgeWts,
                 const int *locPinList, const int *hEdgeOffsets, int dispOption,
                 ostream &out, MPI_Comm comm);

  ~ParaHypergraph();

  void hypergraphFromFile(const char *filename, int dispOption, ostream &out,
                          MPI_Comm comm);
  void initPartitionFromFile(const char *filename, int numParts, ostream &out,
                             MPI_Comm comm);

  void allocHedgeMem(int numHedges, int numLocPins);
  void contractHyperedges(ParaHypergraph &coarse, MPI_Comm comm);
  void contractRestrHyperedges(ParaHypergraph &coarse, MPI_Comm comm);
  void projectPartitions(ParaHypergraph &coarse, MPI_Comm comm);
  void resetVectors();

  void removeBadPartitions(double cutThreshold);
  void setNumberPartitions(int nP);
  void computePartitionChars(int pNum, int numParts, double constraint,
                             ostream &out, MPI_Comm comm);
  void copyInPartition(const int *partition, int numV, int nP);
  void copyOutPartition(int *partition, int numV, int nP) const;
  int keepBestPartition();

  void prescribedVertexShuffle(int *prescribedAssignment, int nLocVer,
                               MPI_Comm comm);
  void prescribedVertexShuffle(int *mapToOrigV, int *partitionFile,
                               MPI_Comm comm);
  void shuffleVerticesByPartition(int nParts, MPI_Comm comm);
  void randomVertexShuffle(MPI_Comm comm);
  void randomVertexShuffle(int *mapToOrigV, MPI_Comm comm);
  // void randomVertexShuffle(int *mapToInterV, int *mapToOrigV, MPI_Comm comm);
  void randomVertexShuffle(ParaHypergraph &fineG, MPI_Comm comm);

  void shuffleVertices(int *vToProc, int *locVPerProc, MPI_Comm comm);
  void shuffleVerticesAftRandom(int *vToProc, int *locVPerProc, int *mapToOrigV,
                                MPI_Comm comm);
  void shuffleVerticesAftRandom(int *vToProc, int *locVPerProc,
                                ParaHypergraph &fineG, MPI_Comm comm);
  void shuffleVerticesAftRandom(int *vToProc, int *locVPerProc,
                                int *mapToInterV, int *mapToOrigV,
                                MPI_Comm comm);

  // void shuffleVerticesWithPartVals(int *vToProc, int *locVPerProc, MPI_Comm
  // comm);
  void shiftVerticesToBalance(MPI_Comm comm);

  int calcCutsize(int numParts, int pNum, MPI_Comm comm);
  int checkBalance(int numParts, double balConstraint, int numPartition,
                   MPI_Comm comm);
  int computeTotalNumPins(MPI_Comm comm);

  void checkValidityOfPartitions(int numParts) const;
  void checkPartitions(int numParts, int maxPartWt, MPI_Comm comm);
  void checkPartitions(int numParts, double constraint, ostream &out,
                       MPI_Comm comm);
  void computeBalanceWarning(int numParts, double constraint, ostream &out,
                             MPI_Comm comm);

  int getNumTotPins(MPI_Comm comm);
  int getNumTotHedges(MPI_Comm comm);
  int getExposedHedgeWt(MPI_Comm comm) const;

  double getAveVertDeg(MPI_Comm comm);
  double getAveHedgeSize(MPI_Comm comm);

  int computeNonConnectedVerts(MPI_Comm comm);

  inline int getNumTotalVertices() const { return numTotalVertices; }
  inline int getNumLocalVertices() const { return numLocalVertices; }
  inline int getNumLocalHedges() const { return numLocalHedges; }
  inline int getNumLocalPins() const { return numLocalPins; }
  inline int getMinVertexIndex() const { return minVertexIndex; }
  inline int getLocalVertexWt() const { return localVertexWt; }
  inline int getNumPartitions() const { return numPartitions; }
  inline int getIndexInSeq() const { return indexInSequence; }
  inline int dontCoarsen() const { return doNotCoarsen; }

  inline int *getWeightArray() const { return vWeight.getArray(); }
  inline int *getMatchVectorArray() const { return matchVector.getArray(); }
  inline int *getPartVectorArray() const { return partitionVector.getArray(); }
  inline int *getLocalPinsArray() const { return localPins.getArray(); }
  inline int *getHedgeOffsetsArray() const { return hEdgeOffsets.getArray(); }
  inline int *getHedgeWeightsArray() const { return hEdgeWeights.getArray(); }
  inline int *getPartitionArray() const { return partitionVector.getArray(); }
  inline int *getPartitionOffsetsArray() const {
    return partitionOffsetsVector.getArray();
  }
  inline int *getCutsizesArray() const {
    return partitionCutsizesVector.getArray();
  }
  inline int *getToOrigVArray() const { return vToOrigV.getArray(); }

  // inline HashKey *getHashKeysArray() const { return hashKeys.getArray(); }
  // inline HedgeIndexTable *getHedgeIndexTable() const { return table; }

  inline void setIndexInSeq(int index) { indexInSequence = index; }
  inline void setNumTotalVertices(int v) { numTotalVertices = v; }
  inline void setNumLocalVertices(int v) { numLocalVertices = v; }
  inline void setNumLocalHedges(int h) { numLocalHedges = h; }
  inline void setNumLocalPins(int p) { numLocalPins = p; }
  inline void setMinVertexIndex(int m) { minVertexIndex = m; }
  inline void setLocalVertexWt(int w) { localVertexWt = w; }

  inline void setWeightArray(int *a, int l) {
    vWeight.setArray(a, l);
  }
  inline void setMatchVectorArray(int *a, int l) {
    matchVector.setArray(a, l);
  }
  inline void setPartVectorArray(int *a, int l) {
    partitionVector.setArray(a, l);
  }
  inline void setLocalPinsArray(int *a, int l) {
    localPins.setArray(a, l);
  }
  inline void setHedgeOffsetsArray(int *a, int l) {
    hEdgeOffsets.setArray(a, l);
  }
  inline void setHedgeWeightsArray(int *a, int l) {
    hEdgeWeights.setArray(a, l);
  }
  // inline void setHedgeIndexTable(HedgeIndexTable *t) { table = t; }
  // inline void setHashKeysArray(HashKey *a, int l) {
  // hashKeys.setArray(a,l); }

  inline void setCut(int pNo, int cut) {
#ifdef DEBUG_HYPERGRAPH
    assert(pNo >= 0 && pNo < numPartitions);
#endif
    partitionCutsizesVector[pNo] = cut;
  }

  inline int getCut(int i) const {
#ifdef DEBUG_HYPERGRAPH
    assert(i >= 0 && i < numPartitions);
#endif
    return partitionCutsizesVector[i];
  }
};

#endif
