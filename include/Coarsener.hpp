
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

#include "hypergraph_loader.hpp"

using namespace std;

class Coarsener : public hypergraph_loader {

protected:
  int minNodes;
  int maxVertexWt;

  double reductionRatio;

public:
  Coarsener(int min, int maxwt, double ratio, int dispL);

  virtual ~Coarsener();
  virtual serial_hypergraph *coarsen(const serial_hypergraph &h) = 0;
  virtual void dispCoarsenerOptions(ostream &out) const = 0;

  serial_hypergraph *buildCoarseHypergraph(int *coarseWts, int numCoarseVerts,
                                    int totWt) const;

  inline void setMaxVertexWt(int maxWt) { maxVertexWt = maxWt; }
  inline int getMaxVertexWt() const { return maxVertexWt; }
};

#endif
