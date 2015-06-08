
#ifndef _HYPERGRAPH_HPP
#define _HYPERGRAPH_HPP

// ### Hypergraph.hpp ###
//
// Copyright (C) 2004,  Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 30/11/2004: Last Modified
//
// ###

#include <unistd.h>
//#  include <fstream>
#include <cstdio>
#include <cassert>

#include "Macros.h"
#include "Funct.hpp"
#include "data_structures/DynamicArray.h"

using namespace std;

class Hypergraph {

protected:
  int totWeight;
  int numVertices;
  int numHedges;
  int numPins;
  int numPartitions;

  DynamicArray<int> vWeight;
  DynamicArray<int> hEdgeWeight;
  DynamicArray<int> matchVector;

  DynamicArray<int> partitionCuts;
  DynamicArray<int> partitionVector;
  DynamicArray<int> partitionVectorOffsets;

  DynamicArray<int> pinList;
  DynamicArray<int> hEdgeOffsets;

  DynamicArray<int> vToHedges;
  DynamicArray<int> vOffsets;

public:
  Hypergraph(int *vWts, int numV);
  Hypergraph(int *vWts, int *pVector, int numV, int cut);
  ~Hypergraph();

  void hypergraphFromFile(const char *filename);
  void buildVtoHedges();

  void resetMatchVector();
  void resetPartitionVectors();
  void resetVertexMaps();

  void projectPartitions(const Hypergraph &coarseGraph);
  void removeBadPartitions(double fractionOK);
  void setNumPartitions(int nPartitions);
  void copyOutPartition(int *pVector, int numV, int pNo) const;
  void copyInPartition(const int *pVector, int numV, int pNo, int cut);
  void printCharacteristics(ostream &o);
  void printPercentiles(ostream &o);

  int keepBestPartition();

  inline int getNumVertices() const { return numVertices; }
  inline int getNumHedges() const { return numHedges; }
  inline int getNumPins() const { return numPins; }
  inline int getTotWeight() const { return totWeight; }
  inline int getNumPartitions() const { return numPartitions; }
  inline int getCut(int pNo) const { return partitionCuts[pNo]; }

  inline int *getVerWeightsArray() const { return vWeight.getArray(); }
  inline int *getHedgeWeightsArray() const { return hEdgeWeight.getArray(); }
  inline int *getMatchVectorArray() const { return matchVector.getArray(); }
  inline int *getPinListArray() const { return pinList.getArray(); }
  inline int *getHedgeOffsetArray() const { return hEdgeOffsets.getArray(); }
  inline int *getVtoHedgesArray() const { return vToHedges.getArray(); }
  inline int *getVerOffsetsArray() const { return vOffsets.getArray(); }

  inline int *getPartVectorArray() const { return partitionVector.getArray(); }
  inline int *getPartOffsetArray() const {
    return partitionVectorOffsets.getArray();
  }
  inline int *getPartCutArray() const { return partitionCuts.getArray(); }
  inline int *getPartitionVector(int pNo) const {
    return (&partitionVector[partitionVectorOffsets[pNo]]);
  }

  inline void setNumHedges(int newNum) { numHedges = newNum; }
  inline void setNumPins(int newNum) { numPins = newNum; }
  inline void setNumVertices(int newNum) { numVertices = newNum; }
  inline void setTotWeight(int newWt) { totWeight = newWt; }
  inline void setWtsArray(int *array, int len) {
    vWeight.setArray(array, len);
  }
  inline void setHedgeWtArray(int *array, int len) {
    hEdgeWeight.setArray(array, len);
  }
  inline void setPinListArray(int *array, int len) {
    pinList.setArray(array, len);
  }
  inline void setHedgeOffsetArray(int *array, int len) {
    hEdgeOffsets.setArray(array, len);
  }
  inline void setVtoHedgesArray(int *array, int len) {
    vToHedges.setArray(array, len);
  }
  inline void setVoffsetsArray(int *array, int len) {
    vOffsets.setArray(array, len);
  }

  inline void setPartitionCutsArray(int *a, int len) {
    partitionCuts.setArray(a, len);
  }
  inline void setPartitionVectorArray(int *a, int len) {
    partitionVector.setArray(a, len);
  }

  int getExposedHedgeWt() const;
  int calcCutsize(int nP, int partitionNo) const;
  int getSOED(int nP, int partitionNo) const;

  void initCutsizes(int numParts);
  void checkPartitions(int nP, int maxWt) const;
  void checkPartition(int partitionNum, int numParts, int maxPartWt) const;

  void convertToDIMACSGraphFile(const char *fName) const;
};

#endif
