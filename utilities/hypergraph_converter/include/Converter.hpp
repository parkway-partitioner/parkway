#ifndef _CONVERTER_HPP
#define _CONVERTER_HPP

// ### Converter.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 01/12/2004: Last Modified
//
// ###

#include <fstream>
#include "data_structures/DynamicArray.h"

using namespace std;

class Converter {

protected:
  int locOfLengthParam;
  int dataLength;
  int numVerts;

  DynamicArray<int> hEdgeData;
  DynamicArray<int> vWeights;

public:
  Converter();
  virtual ~Converter();

  inline int getDataLength() { return dataLength; }
  inline int getNumVerts() { return numVerts; }
  inline int *getDataArray() { return hEdgeData.getArray(); }
  inline int *getVertWtsArray() { return vWeights.getArray(); }
  inline void setNumVertices(int n) { numVerts = n; }

  void setLengthParameter();
  void resetConverterParameters();
  void addLengthParameter();

  virtual void convert(const char *filename) = 0;
  virtual void readPreamble(ifstream &in) = 0;
};

#endif
