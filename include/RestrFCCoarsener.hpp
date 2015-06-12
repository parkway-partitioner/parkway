
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

#include "RestrCoarsener.hpp"

using namespace std;

class RestrFCCoarsener : public RestrCoarsener {

protected:
  int utilFanOut;
  int divByCluWt;

public:
  RestrFCCoarsener(int _min, int _maxwt, double r, int fanOut, int dbWt,
                   int dL);
  ~RestrFCCoarsener();

  hypergraph *coarsen(const hypergraph &h);

  void dispCoarsenerOptions(ostream &out) const;

  inline void setUtilFanOut(int f) { utilFanOut = f; }
  inline int getUtilFanOut() const { return utilFanOut; }
};

#endif
