
#ifndef _HARWELLBOEING_READER_HPP
#define _HARWELLBOEING_READER_HPP

// ### HarwellBoeingReader.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 30/12/2004: Last Modified
//
// ###

#include "StringUtils.hpp"
#include "TextFileReader.hpp"

using namespace std;

class HarwellBoeingReader : public TextFileReader {

protected:
  int totDataLines;
  int totColPtrLines;
  int totRowIdxLines;
  int totNumValLines;
  int numRHSLines;

  int numRows;
  int numCols;
  int numNonZeros;
  int numElemMatEntries;

  DynamicArray<char> title;
  DynamicArray<char> id;
  DynamicArray<char> mtxType;

public:
  HarwellBoeingReader();
  ~HarwellBoeingReader();

  void readPreamble(ifstream &in);
};

#endif
