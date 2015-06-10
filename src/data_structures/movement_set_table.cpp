
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

#include "data_structures/movement_set_table.hpp"

movement_set_table::movement_set_table(int nParts, int nProcs) {
  number_of_parts_ = nParts;
  number_of_processors_ = nProcs;
  setArrayLen = number_of_parts_ * number_of_parts_;
  max_part_weight_ = 0;

  part_weights_.reserve(number_of_parts_);
  sets_.reserve(setArrayLen);

  int i;
  int j;

  for (i = 0; i < setArrayLen; ++i) {
    if (Mod(i, number_of_parts_) != i / number_of_parts_) {
      sets_[i] = new dynamic_array<movement_set>(number_of_processors_);

      for (j = 0; j < number_of_processors_; ++j)
        (*sets_[i])[j].proc = j;
    } else {
      sets_[i] = nullptr;
    }
  }

  restoring_moves_.reserve(number_of_processors_);
  restoring_move_lens_.reserve(number_of_processors_);

  for (i = 0; i < number_of_processors_; ++i)
    restoring_moves_[i] = new dynamic_array<int>(256);
}

movement_set_table::~movement_set_table() {
  int i;

  for (i = 0; i < setArrayLen; ++i)
    DynaMem::deletePtr<dynamic_array<movement_set> >(sets_[i]);

  for (i = 0; i < number_of_processors_; ++i)
    DynaMem::deletePtr<dynamic_array<int> >(restoring_moves_[i]);
}

void movement_set_table::initialize_part_weights(const int *partWts, int nParts) {
#ifdef DEBUG_BASICS
  assert(nParts == numParts);
#endif

  int i;
  int j;

  for (i = 0; i < number_of_parts_; ++i)
    part_weights_[i] = partWts[i];

  for (j = 0; j < setArrayLen; ++j) {
    if (Mod(j, number_of_parts_) != j / number_of_parts_) {
      for (i = 0; i < number_of_processors_; ++i) {
        (*sets_[j])[i].gain = -1;
        (*sets_[j])[i].weight = 0;
      }
    }
  }
}

void movement_set_table::complete_processor_sets(int proc, int dataLen,
                                                 const int *data) {
  int i = 0;
  int wt;

  int fromPart;
  int toPart;
  int setArray;

  while (i < dataLen) {
    fromPart = data[i];
    toPart = data[i + 1];
    setArray = fromPart * number_of_parts_ + toPart;
    wt = data[i + 3];

#ifdef DEBUG_BASICS
    assert((*sets[setArray])[proc].proc == proc);
#endif

    (*sets_[setArray])[proc].gain = data[i + 2];
    (*sets_[setArray])[proc].weight = wt;

    part_weights_[fromPart] -= wt;
    part_weights_[toPart] += wt;

    i += 4;
  }
}

void movement_set_table::compute_restoring_array() {
  int heaviest = find_heaviest_part_index();

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

  for (i = 0; i < number_of_processors_; ++i)
    restoring_move_lens_[i] = 0;

  while (heaviest > -1) {
    minIndex = -1;
    minProc = -1;
    minGain = LARGE_CONSTANT;

    for (j = 0; j < number_of_parts_; ++j) {
      if (j != heaviest) {
        prod = j * number_of_parts_;

        for (i = 0; i < number_of_processors_; ++i) {
          gain = (*sets_[prod + heaviest])[i].gain;

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

    prod = minIndex * number_of_parts_;
    weight = (*sets_[prod + heaviest])[minProc].weight;

#ifdef DEBUG_BASICS
    assert(weight > 0);
#endif

    restoring_moves_[minProc]->assign(restoring_move_lens_[minProc]++, minIndex);
    restoring_moves_[minProc]->assign(restoring_move_lens_[minProc]++, heaviest);

    part_weights_[minIndex] += weight;
    part_weights_[heaviest] -= weight;

    (*sets_[prod + heaviest])[minProc].weight = 0;
    (*sets_[prod + heaviest])[minProc].gain = -1;

    heaviest = find_heaviest_part_index();
  }
}

#endif
