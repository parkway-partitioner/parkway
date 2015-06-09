
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
#include <cstdio>
#include <cassert>
#include <iostream>

#include "Macros.h"
#include "Funct.hpp"
#include "data_structures/dynamic_array.hpp"

using parkway::data_structures::dynamic_array;

class Hypergraph {
 protected:
  int totWeight;
  int numVertices;
  int numHedges;
  int numPins;
  int numPartitions;

  dynamic_array<int> vWeight;
  dynamic_array<int> hEdgeWeight;
  dynamic_array<int> matchVector;

  dynamic_array<int> partitionCuts;
  dynamic_array<int> partitionVector;
  dynamic_array<int> partitionVectorOffsets;

  dynamic_array<int> pinList;
  dynamic_array<int> hEdgeOffsets;

  dynamic_array<int> vToHedges;
  dynamic_array<int> vOffsets;

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
  void printCharacteristics(std::ostream &o);
  void printPercentiles(std::ostream &o);

  int keepBestPartition();

  inline int getNumVertices() const { return numVertices; }
  inline int getNumHedges() const { return numHedges; }
  inline int getNumPins() const { return numPins; }
  inline int getTotWeight() const { return totWeight; }
  inline int getNumPartitions() const { return numPartitions; }
  inline int getCut(int pNo) const { return partitionCuts[pNo]; }

  inline int *getVerWeightsArray() const { return vWeight.data(); }
  inline int *getHedgeWeightsArray() const { return hEdgeWeight.data(); }
  inline int *getMatchVectorArray() const { return matchVector.data(); }
  inline int *getPinListArray() const { return pinList.data(); }
  inline int *getHedgeOffsetArray() const { return hEdgeOffsets.data(); }
  inline int *getVtoHedgesArray() const { return vToHedges.data(); }
  inline int *getVerOffsetsArray() const { return vOffsets.data(); }

  inline int *getPartVectorArray() const { return partitionVector.data(); }
  inline int *getPartOffsetArray() const {
    return partitionVectorOffsets.data();
  }
  inline int *getPartCutArray() const { return partitionCuts.data(); }
  inline int *getPartitionVector(int pNo) const {
    return (&partitionVector[partitionVectorOffsets[pNo]]);
  }

  inline void setNumHedges(int newNum) { numHedges = newNum; }
  inline void setNumPins(int newNum) { numPins = newNum; }
  inline void setNumVertices(int newNum) { numVertices = newNum; }
  inline void setTotWeight(int newWt) { totWeight = newWt; }
  inline void setWtsArray(int *array, int len) {
   vWeight.set_data(array, len);
  }
  inline void setHedgeWtArray(int *array, int len) {
   hEdgeWeight.set_data(array, len);
  }
  inline void setPinListArray(int *array, int len) {
   pinList.set_data(array, len);
  }
  inline void setHedgeOffsetArray(int *array, int len) {
   hEdgeOffsets.set_data(array, len);
  }
  inline void setVtoHedgesArray(int *array, int len) {
   vToHedges.set_data(array, len);
  }
  inline void setVoffsetsArray(int *array, int len) {
   vOffsets.set_data(array, len);
  }

  inline void setPartitionCutsArray(int *a, int len) {
   partitionCuts.set_data(a, len);
  }
  inline void setPartitionVectorArray(int *a, int len) {
   partitionVector.set_data(a, len);
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
