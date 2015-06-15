#ifndef _REFINER_HPP
#define _REFINER_HPP

// ### Refiner.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 30/11/2004: Last Modified
//
// ###

#include <ostream>
#include "Macros.h"
#include "data_structures/dynamic_array.hpp"
#include "hypergraph/serial/hypergraph.hpp"
#include "hypergraph/serial/loader.hpp"

namespace serial = parkway::serial;
using parkway::data_structures::dynamic_array;

class refiner : public serial::loader {
 protected:
  int maximum_part_weight_;
  int number_of_parts_;
  int *partition_vector_;

  double accept_proportion_;
  double average_part_weight_;

  dynamic_array<int> part_weights_;

 public:
  refiner(int dispL);

  virtual ~refiner();
  virtual void refine(serial::hypergraph &h) = 0;
  virtual void display_options(std::ostream &out) const = 0;

  inline int maximum_part_weight() const { return maximum_part_weight_; }

  inline void set_maximum_part_weight(int max) { maximum_part_weight_ = max; }
  inline void set_average_part_weight(double ave) { average_part_weight_ = ave; }
  inline void set_partition_vector(int nP) {
    partition_vector_ = &partitionVectors[partitionOffsets[nP]];
  }

  int calculate_cut_size() const;
};

#endif
