
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

namespace serial = parkway::hypergraph::serial;

class Coarsener : public serial::loader {

protected:
  int minNodes;
  int maxVertexWt;

  double reductionRatio;

public:
  Coarsener(int min, int maxwt, double ratio, int dispL);

  virtual ~Coarsener();
  virtual serial::hypergraph *coarsen(const serial::hypergraph &h) = 0;
  virtual void dispCoarsenerOptions(std::ostream &out) const = 0;

  serial::hypergraph *buildCoarseHypergraph(int *coarseWts, int numCoarseVerts,
                                    int totWt) const;

  inline void setMaxVertexWt(int maxWt) { maxVertexWt = maxWt; }
  inline int getMaxVertexWt() const { return maxVertexWt; }
};

#endif
