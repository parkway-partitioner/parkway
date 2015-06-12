
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

#include "hypergraph_loader.hpp"

using namespace std;

class RestrCoarsener : public hypergraph_loader {

protected:
  int minNodes;
  int maxVertexWt;

  double reductionRatio;

public:
  RestrCoarsener(int min, int maxwt, double ratio, int dispL);

  virtual ~RestrCoarsener();
  virtual serial_hypergraph *coarsen(const serial_hypergraph &h) = 0;
  virtual void dispCoarsenerOptions(ostream &out) const = 0;

  serial_hypergraph *buildCoarseHypergraph(int *coarseWts, int *pVector,
                                    int numCoarseVerts, int totWt) const;

  inline void setMaxVertexWt(int maxWt) { maxVertexWt = maxWt; }
};

#endif
