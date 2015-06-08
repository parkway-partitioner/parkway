#ifndef _HARWELLBOEING_TOBIN_HPP
#define _HARWELLBOEING_TOBIN_HPP

// ### HarwellBoeingToBin.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 13/1/2005: Last Modified
//
// NOTE: Assume a row-vertex, col-hyperedge conversion from
//       sparse matrix to hypergraph
//
// ###

#include "HarwellBoeingReader.hpp"

class HarwellBoeingToBin : public HarwellBoeingReader {
 protected:
  int hypergraphModelType;

 public:
  HarwellBoeingToBin(int model);
  ~HarwellBoeingToBin();

  void convert(const char *filename);
};

#endif
