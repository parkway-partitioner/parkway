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

  inline void set_max_part_weight(int weight) {
    max_part_weight_ = weight;
  }

  inline int *part_weights_array() const {
    return part_weights_.data();
  }

  inline int *restoring_move_lens() const {
    return restoring_move_lens_.data();
  }

  inline dynamic_array<int> **restoring_moves() const {
    return restoring_moves_.data();
  }

  inline int find_heaviest_part_index() const {
    int j = -1;
    for (int i = 0; i < number_of_parts_; ++i) {
      if (part_weights_[i] > max_part_weight_) {
        if (j == -1 || (j >= 0 && part_weights_[i] > part_weights_[j])) {
          j = i;
        }
      }
    }

#ifdef DEBUG_BASICS
    assert(j >= -1 && j < numParts);
#endif

    return j;
  }

  void initialize_part_weights(const int *partWts, int nParts);
  void complete_processor_sets(int proc, int dataLen, const int *data);
  void compute_restoring_array();

 protected:
  int number_of_parts_;
  // TODO(gb610): processors or processes?
  int number_of_processors_;
  int max_part_weight_;
  int setArrayLen;

  dynamic_array<int> part_weights_;
  dynamic_array<int> restoring_move_lens_;

  dynamic_array<dynamic_array<int> *> restoring_moves_;
  dynamic_array<dynamic_array<movement_set> *> sets_;
};

#endif
