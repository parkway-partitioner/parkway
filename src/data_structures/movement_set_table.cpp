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

namespace parkway {
namespace data_structures {

movement_set_table::movement_set_table(int number_of_parts,
                                       int number_of_processors)
    : number_of_parts_(number_of_parts),
      number_of_processors_(number_of_processors),
      set_array_len_(number_of_parts_ * number_of_parts_),
      max_part_weight_(0) {

  part_weights_.reserve(number_of_parts_);
  sets_.reserve(set_array_len_);

  for (int i = 0; i < set_array_len_; ++i) {
    if ((i % number_of_parts_) != i / number_of_parts_) {
      sets_[i] = new dynamic_array<movement_set>(number_of_processors_);
      for (int j = 0; j < number_of_processors_; ++j) {
        (*sets_[i])[j].proc = j;
      }
    } else {
      sets_[i] = nullptr;
    }
  }

  restoring_moves_.reserve(number_of_processors_);
  restoring_move_lens_.reserve(number_of_processors_);

  for (int i = 0; i < number_of_processors_; ++i) {
    restoring_moves_[i] = new dynamic_array<int>(256);
  }
}

movement_set_table::~movement_set_table() {
  for (int i = 0; i < set_array_len_; ++i) {
    dynamic_memory::delete_pointer<dynamic_array<movement_set> >(sets_[i]);
  }

  for (int i = 0; i < number_of_processors_; ++i) {
    dynamic_memory::delete_pointer<dynamic_array<int> >(restoring_moves_[i]);
  }
}

void movement_set_table::initialize_part_weights(const int *part_weights,
                                                 int number_of_parts) {
#ifdef DEBUG_BASICS
  assert(number_of_parts == number_of_parts_);
#endif

  for (int i = 0; i < number_of_parts_; ++i) {
    part_weights_[i] = part_weights[i];
  }

  for (int j = 0; j < set_array_len_; ++j) {
    if ((j % number_of_parts_) != j / number_of_parts_) {
      for (int i = 0; i < number_of_processors_; ++i) {
        (*sets_[j])[i].gain = -1;
        (*sets_[j])[i].weight = 0;
      }
    }
  }
}

void movement_set_table::complete_processor_sets(int proc, int data_length,
                                                 const int *data) {
  // Each part of the data is split into 4 parts: from, to, gain and weight.
  for (int i = 0; i < data_length; i += 4) {
    int from_part = data[i];
    int to_part = data[i + 1];
    int set_array = from_part * number_of_parts_ + to_part;
    int weight = data[i + 3];

#ifdef DEBUG_BASICS
    assert((*sets[set_array])[proc].proc == proc);
#endif

    (*sets_[set_array])[proc].gain = data[i + 2];
    (*sets_[set_array])[proc].weight = weight;

    part_weights_[from_part] -= weight;
    part_weights_[to_part] += weight;
  }
}

void movement_set_table::compute_restoring_array() {
  int heaviest = find_heaviest_part_index();

#ifdef DEBUG_BASICS
  assert(heaviest >= -1);
#endif


  for (int i = 0; i < number_of_processors_; ++i) {
    restoring_move_lens_[i] = 0;
  }

  int prod;
  while (heaviest > -1) {
    int minIndex = -1;
    int minProc = -1;
    int minGain = LARGE_CONSTANT;

    for (int j = 0; j < number_of_parts_; ++j) {
      if (j != heaviest) {
        prod = j * number_of_parts_;
        for (int i = 0; i < number_of_processors_; ++i) {
          int gain = (*sets_[prod + heaviest])[i].gain;
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
    int weight = (*sets_[prod + heaviest])[minProc].weight;

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

}  // namespace data_structures
}  // namespace parkway

#endif
