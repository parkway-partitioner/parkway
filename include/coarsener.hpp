
#ifndef _COARSENER_HPP
#define _COARSENER_HPP

// ### Coarsener.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 4/1/2005: Last Modified
//
// ###

#include "hypergraph/serial/hypergraph.hpp"
#include "hypergraph/serial/loader.hpp"

namespace serial = parkway::serial;

class coarsener : public serial::loader {
protected:
  int minimum_number_of_nodes_;
  int maximum_vertex_weight_;

  double reduction_ratio_;

public:
  coarsener(int min, int maxwt, double ratio, int dispL);

  virtual ~coarsener();
  virtual serial::hypergraph *coarsen(const serial::hypergraph &h) = 0;
  virtual void display_options(std::ostream &out) const = 0;

  serial::hypergraph *build_coarse_hypergraph(int *coarseWts,
                                              int numCoarseVerts,
                                              int totWt) const;

  inline void set_maximum_vertex_weight(int maxWt) { maximum_vertex_weight_ = maxWt; }
  inline int maximum_vertex_weight() const { return maximum_vertex_weight_; }
};

#endif
