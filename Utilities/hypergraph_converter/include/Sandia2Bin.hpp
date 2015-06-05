#ifndef _SANDIA2BIN_HPP
#define _SANDIA2BIN_HPP

// ### Sandia2Bin.hpp ###
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

class Sandia2Bin : public HypergraphTextFileReader, TextToBinConverter {
protected:
public:
  Sandia2Bin();
  ~Sandia2Bin();

  void convert(const char *filename);
  void readPreamble(ifstream &in);

  int processHedgeLine(register char *line, register int &numP);
  int processVertLine(register char *line);
};

#endif
