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
#include "data_structures/dynamic_array.hpp"

using parkway::data_structures::dynamic_array;

class Converter {
 protected:
  int locOfLengthParam;
  int dataLength;
  int numVerts;

  dynamic_array<int> hEdgeData;
  dynamic_array<int> vWeights;

 public:
  Converter();
  virtual ~Converter();

  inline int getDataLength() { return dataLength; }
  inline int getNumVerts() { return numVerts; }
  inline int *getDataArray() { return hEdgeData.data(); }
  inline int *getVertWtsArray() { return vWeights.data(); }
  inline void setNumVertices(int n) { numVerts = n; }

  void setLengthParameter();
  void resetConverterParameters();
  void addLengthParameter();

  virtual void convert(const char *filename) = 0;
  virtual void readPreamble(std::ifstream &in) = 0;
};

#endif
