#ifndef _FCCOARSENER_HPP
#define _FCCOARSENER_HPP
// ### FCCoarsener.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 30/11/2004: Last Modified
//
// ###

#include <iostream>
#include "coarseners/serial/coarsener.hpp"
#include "hypergraph/serial/hypergraph.hpp"

namespace parkway {
namespace serial {

class first_choice_coarsener : public coarsener {
 public:
  first_choice_coarsener(int _min, int _maxwt, double r, int fanOut, int dbWt, int dL);
  ~first_choice_coarsener();

  hypergraph *coarsen(const hypergraph &h);

  void display_options(std::ostream &out) const;

  inline void set_util_fan_out(int f) {
    util_fan_out_ = f;
  }

  inline int util_fan_out() const {
    return util_fan_out_;
  }

 protected:
  int util_fan_out_;
  int divide_by_cluster_weight_;
};

}  // namespace serial
}  // namespace parkway

#endif
