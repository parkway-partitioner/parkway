
#ifndef _MTXMKT_READER_HPP
#define _MTXMKT_READER_HPP

// ### MatrixMarketReader.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 12/1/2004: Last Modified
//
// ###

#include "TextFileReader.hpp"
#include "StringUtils.hpp"

using namespace std;

class MatrixMarketReader : public TextFileReader {
protected:
  int numVertices;
  int numHyperedges;
  int numPins;

  DynamicArray<int> hEdgeOffsets;
  DynamicArray<int> pinList;

public:
  MatrixMarketReader();
  ~MatrixMarketReader();

  void readMatrix(const char *filename);
  void readPreamble(ifstream &in_stream);
  void readMatrix(ifstream &in_stream);
  void skipComments(ifstream &in_stream);
};

#endif
