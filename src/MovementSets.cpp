
#ifndef _MOVEMENT_SETS_CPP
#define _MOVEMENT_SETS_CPP

// ### MovementSets.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 30/11/2004: Last Modified
//
// ###

#include "data_structures/MovementSets.hpp"

MovementSetTable::MovementSetTable(int nParts, int nProcs) {
  numParts = nParts;
  numProcs = nProcs;
  setArrayLen = numParts * numParts;
  maxPartWt = 0;

  partWeights.setLength(numParts);
  sets.setLength(setArrayLen);

  int i;
  int j;

  for (i = 0; i < setArrayLen; ++i) {
    if (Mod(i, numParts) != i / numParts) {
      sets[i] = new DynamicArray<MOVE_SET>(numProcs);

      for (j = 0; j < numProcs; ++j)
        (*sets[i])[j].proc = j;
    } else {
      sets[i] = nullptr;
    }
  }

  restoringMoves.setLength(numProcs);
  restoringMovesLens.setLength(numProcs);

  for (i = 0; i < numProcs; ++i)
    restoringMoves[i] = new DynamicArray<int>(256);
}

MovementSetTable::~MovementSetTable() {
  int i;

  for (i = 0; i < setArrayLen; ++i)
    DynaMem<DynamicArray<MOVE_SET> >::deletePtr(sets[i]);

  for (i = 0; i < numProcs; ++i)
    DynaMem<DynamicArray<int> >::deletePtr(restoringMoves[i]);
}

void MovementSetTable::initPartWeights(const int *partWts, int nParts) {
#ifdef DEBUG_BASICS
  assert(nParts == numParts);
#endif

  int i;
  int j;

  for (i = 0; i < numParts; ++i)
    partWeights[i] = partWts[i];

  for (j = 0; j < setArrayLen; ++j) {
    if (Mod(j, numParts) != j / numParts) {
      for (i = 0; i < numProcs; ++i) {
        (*sets[j])[i].gain = -1;
        (*sets[j])[i].weight = 0;
      }
    }
  }
}

void MovementSetTable::completeProcSets(int proc, int dataLen,
                                        const int *data) {
  int i = 0;
  int wt;

  int fromPart;
  int toPart;
  int setArray;

  while (i < dataLen) {
    fromPart = data[i];
    toPart = data[i + 1];
    setArray = fromPart * numParts + toPart;
    wt = data[i + 3];

#ifdef DEBUG_BASICS
    assert((*sets[setArray])[proc].proc == proc);
#endif

    (*sets[setArray])[proc].gain = data[i + 2];
    (*sets[setArray])[proc].weight = wt;

    partWeights[fromPart] -= wt;
    partWeights[toPart] += wt;

    i += 4;
  }
}

void MovementSetTable::computeRestoringArray() {
  int heaviest = findHeaviest();

#ifdef DEBUG_BASICS
  assert(heaviest >= -1);
#endif

  int i;
  int j;

  int minIndex;
  int minProc;
  int weight;
  int gain;
  int minGain;
  int prod;

  for (i = 0; i < numProcs; ++i)
    restoringMovesLens[i] = 0;

  while (heaviest > -1) {
    minIndex = -1;
    minProc = -1;
    minGain = LARGE_CONSTANT;

    for (j = 0; j < numParts; ++j) {
      if (j != heaviest) {
        prod = j * numParts;

        for (i = 0; i < numProcs; ++i) {
          gain = (*sets[prod + heaviest])[i].gain;

          if (gain >= 0 && gain < minGain) {
            minGain = gain;
            minIndex = j;
            minProc = i;
          }
        }
      }
    }

#ifdef DEBUG_BASICS
    assert(minGain != LARGE_CONSTANT);
    assert(minProc != -1);
    assert(minIndex != -1);
#endif

    prod = minIndex * numParts;
    weight = (*sets[prod + heaviest])[minProc].weight;

#ifdef DEBUG_BASICS
    assert(weight > 0);
#endif

    restoringMoves[minProc]->assign(restoringMovesLens[minProc]++, minIndex);
    restoringMoves[minProc]->assign(restoringMovesLens[minProc]++, heaviest);

    partWeights[minIndex] += weight;
    partWeights[heaviest] -= weight;

    (*sets[prod + heaviest])[minProc].weight = 0;
    (*sets[prod + heaviest])[minProc].gain = -1;

    heaviest = findHeaviest();
  }
}

#endif
