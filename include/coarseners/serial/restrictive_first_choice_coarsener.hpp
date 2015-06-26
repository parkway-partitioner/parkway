#ifndef _RESTR_FCCOARSENER_HPP
#define _RESTR_FCCOARSENER_HPP
// ### RestrFCCoarsener.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 30/11/2004: Last Modified
//
// ###
#include "coarseners/serial/restrictive_coarsener.hpp"
#include "hypergraph/serial/hypergraph.hpp"

namespace parkway {
namespace serial {

class restrictive_first_choice_coarsener : public restrictive_coarsener {
 public:
  restrictive_first_choice_coarsener(int _min, int _maxwt, double r, int fanOut,
                                     int dbWt);
  ~restrictive_first_choice_coarsener();

  hypergraph *coarsen(const hypergraph &h);

  void display_options() const;

  inline void set_util_fan_out(int f) { util_fan_out_ = f; }
  inline int util_fan_out() const { return util_fan_out_; }

 protected:
  int util_fan_out_;
  int divide_by_cluster_weight_;
};

}  // namespace serial
}  // namespace parkway

#endif
