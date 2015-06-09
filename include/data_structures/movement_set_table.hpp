
#ifndef _MOVEMENT_SETS_HPP
#define _MOVEMENT_SETS_HPP

// ### MovementSets.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 30/11/2004: Last Modified
//
// ###

#include "data_structures/dynamic_array.hpp"

using parkway::data_structures::dynamic_array;

struct movement_set {
  int gain;
  int weight;
  int proc;
};

class movement_set_table {
 public:
  movement_set_table(int nParts, int nProcs);
  ~movement_set_table();

  inline void setMaxWt(int max) { maxPartWt = max; }

  inline int *getPartWeightsArray() const { return partWeights.data(); }
  inline int *getRestoringMovesLens() const {
    return restoringMovesLens.data();
  }
  inline dynamic_array<int> **getRestoringMoves() const {
    return restoringMoves.data();
  }

  inline int findHeaviest() const {
    int i = 0;
    int j = -1;

    for (; i < numParts; ++i)
      if (partWeights[i] > maxPartWt)
        if (j == -1 || (j >= 0 && partWeights[i] > partWeights[j]))
          j = i;

#ifdef DEBUG_BASICS
    assert(j >= -1 && j < numParts);
#endif

    return j;
  }

  void initPartWeights(const int *partWts, int nParts);
  void completeProcSets(int proc, int dataLen, const int *data);
  void computeRestoringArray();

 protected:
  int numParts;
  int numProcs;
  int maxPartWt;
  int setArrayLen;

  dynamic_array<int> partWeights;
  dynamic_array<int> restoringMovesLens;

  dynamic_array<dynamic_array<int> *> restoringMoves;
  dynamic_array<dynamic_array<movement_set> *> sets;
};

#endif
