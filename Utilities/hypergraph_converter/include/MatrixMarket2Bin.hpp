
#ifndef _MTXMKT_TOBIN_HPP
#define _MTXMKT_TOBIN_HPP

// ### MatrixMarket2Bin.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 12/1/2004: Last Modified
//
// ###

#include "MatrixMarketReader.hpp"

using namespace std;

class MatrixMarket2Bin : public MatrixMarketReader {
protected:
public:
  MatrixMarket2Bin();
  ~MatrixMarket2Bin();

  void convert(const char *filename);
  void writeMatrix(const char *filename);
};

#endif
