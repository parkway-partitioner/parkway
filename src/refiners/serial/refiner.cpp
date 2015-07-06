// ### Refiner.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 4/1/2005: Last Modified
//
// ###
#include "refiners/serial/refiner.hpp"

namespace parkway {
namespace serial {

refiner::refiner()
    : number_of_parts_(0),
      accept_proportion_(0),
      partition_vector_(nullptr) {
  part_weights_.resize(0);
}

refiner::~refiner() {}

int refiner::calculate_cut_size() const {
#ifdef DEBUG_REFINER
  assert(numParts > 0);
#endif

  ds::dynamic_array<int> spanned(number_of_parts_);

  int k_1Cut = 0;
  int endOffset;
  int vPart;
  int numSpanned;

  for (int i = 0; i < numHedges; ++i) {
    endOffset = hEdgeOffsets[i + 1];
    numSpanned = 0;

    for (int j = 0; j < number_of_parts_; ++j) {
      spanned[j] = 0;
    }

    for (int j = hEdgeOffsets[i]; j < endOffset; ++j) {
      vPart = partition_vector_[pinList[j]];
#ifdef DEBUG_REFINER
      assert(vPart >= 0 && vPart < numParts);
#endif
      if (!spanned[vPart]) {
        spanned[vPart] = 1;
        ++numSpanned;
      }
    }

    k_1Cut += ((numSpanned - 1) * hEdgeWeight[i]);
  }

  return k_1Cut;
}

}  // namespace serial
}  // namespace parkway
