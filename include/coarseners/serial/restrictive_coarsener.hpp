#ifndef _RESTR_COARSENER_HPP
#define _RESTR_COARSENER_HPP
// ### RestrCoarsener.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 4/1/2005: Last Modified
//
// ###

#include <iostream>
#include "hypergraph/serial/hypergraph.hpp"
#include "hypergraph/serial/loader.hpp"

namespace parkway {
namespace serial {

class restrictive_coarsener : public loader {
 public:
  restrictive_coarsener(int min, int maxwt, double ratio, int dispL);

  virtual ~restrictive_coarsener();
  virtual serial::hypergraph *coarsen(const serial::hypergraph &h) = 0;
  virtual void display_options(std::ostream &out) const = 0;

  serial::hypergraph *build_coarse_hypergraph(int *coarseWts, int *pVector,
                                              int numCoarseVerts, int totWt) const;

  inline void set_maximum_vertex_weight(int maxWt) {
    maximum_vertex_weight_ = maxWt;
  }

 protected:
  int minimum_nodes_;
  int maximum_vertex_weight_;
  double reduction_ratio_;
};

}  // namespace serial
}  // namespace parkway

#endif
