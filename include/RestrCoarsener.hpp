
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

#include "HypergraphLoader.hpp"

using namespace std;

class RestrCoarsener : public HypergraphLoader {

protected:
  int minNodes;
  int maxVertexWt;

  double reductionRatio;

public:
  RestrCoarsener(int min, int maxwt, double ratio, int dispL);

  virtual ~RestrCoarsener();
  virtual Hypergraph *coarsen(const Hypergraph &h) = 0;
  virtual void dispCoarsenerOptions(ostream &out) const = 0;

  Hypergraph *buildCoarseHypergraph(int *coarseWts, int *pVector,
                                    int numCoarseVerts, int totWt) const;

  inline void setMaxVertexWt(register int maxWt) { maxVertexWt = maxWt; }
};

#endif
