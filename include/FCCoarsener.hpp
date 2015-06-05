
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

#include "Coarsener.hpp"

using namespace std;

class FCCoarsener : public Coarsener {

protected:
  int utilFanOut;
  int divByCluWt;

public:
  FCCoarsener(int _min, int _maxwt, double r, int fanOut, int dbWt, int dL);
  ~FCCoarsener();

  Hypergraph *coarsen(const Hypergraph &h);

  void dispCoarsenerOptions(ostream &out) const;

  inline void setUtilFanOut(register int f) { utilFanOut = f; }
  inline int getUtilFanOut() const { return utilFanOut; }
};

#endif
