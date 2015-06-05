#ifndef _PATOH2BIN_HPP
#define _PATOH2BIN_HPP

// ### PaToH2Bin.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 02/12/2004: Last Modified
//
// NOTE: currently do not deal with constraints
//
// ###

#include "HypergraphTextFileReader.hpp"
#include "TextToBinConverter.hpp"

using namespace std;

class PaToH2Bin : public HypergraphTextFileReader, TextToBinConverter {

protected:
  int numScheme;
  int constraints;

public:
  PaToH2Bin();
  ~PaToH2Bin();

  void convert(const char *filename);
  void readPreamble(ifstream &in);

  int processHedgeLine(char *line, int &numP);
  int processVertLine(char *line);
};

#endif
