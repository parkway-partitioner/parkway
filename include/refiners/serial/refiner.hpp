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
#include "internal/base/refiner.hpp"

namespace parkway {
namespace serial {
namespace ds = parkway::data_structures;

class refiner : public loader, public parkway::base::refiner {
 public:
  refiner();

  virtual ~refiner();
  virtual void refine(serial::hypergraph &h) = 0;

  inline void set_partition_vector(int nP) {
    partition_vector_ = &partitionVectors.data()[partitionOffsets[nP]];
  }

  int calculate_cut_size() const;

 protected:
  int number_of_parts_;
  int *partition_vector_;
  double accept_proportion_;
};

}  // namespace serial
}  // namespace parkway

#endif
