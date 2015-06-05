#ifndef _FROMBIN_CONVERTER_HPP
#define _FROMBIN_CONVERTER_HPP

// ### FromBinConverter.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 07/12/2004: Last Modified
//
// ###

#include "Converter.hpp"

using namespace std;

class FromBinConverter : public Converter {
protected:
  int numHedges;
  int numPins;

public:
  FromBinConverter();

  virtual ~FromBinConverter();
  virtual void convert(const char *filename) = 0;

  void readPreamble(ifstream &in);
  void readInVertexWts(ifstream &in_stream);
  void readInHedgeData(ifstream &in_stream, int &inStream);
};

#endif
