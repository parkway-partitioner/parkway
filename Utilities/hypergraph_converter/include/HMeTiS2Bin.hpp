#ifndef _HMETIS2BIN_HPP
#define _HMETIS2BIN_HPP

// ### HMeTiS2Bin.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 01/12/2004: Last Modified
//
// ###

#include "HypergraphTextFileReader.hpp"
#include "TextToBinConverter.hpp"

using namespace std;

class HMeTiS2Bin : public HypergraphTextFileReader, TextToBinConverter {

protected:
public:
  HMeTiS2Bin();
  ~HMeTiS2Bin();

  void convert(const char *filename);
  void readPreamble(ifstream &in);

  int processHedgeLine(register char *line, register int &numP);
  int countPinsInLine(register char *line, register int &nP);
  int processVertLine(register char *line);
  int computeNumPins(ifstream &in_stream);
};

#endif
