
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

namespace serial = parkway::hypergraph::serial;

class RestrCoarsener : public serial::loader {

protected:
  int minNodes;
  int maxVertexWt;

  double reductionRatio;

public:
  RestrCoarsener(int min, int maxwt, double ratio, int dispL);

  virtual ~RestrCoarsener();
  virtual serial::hypergraph *coarsen(const serial::hypergraph &h) = 0;
  virtual void dispCoarsenerOptions(std::ostream &out) const = 0;

  serial::hypergraph *buildCoarseHypergraph(int *coarseWts, int *pVector,
                                    int numCoarseVerts, int totWt) const;

  inline void setMaxVertexWt(int maxWt) { maxVertexWt = maxWt; }
};

#endif
