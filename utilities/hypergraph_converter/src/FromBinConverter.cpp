#ifndef _FROMBIN_CONVERTER_CPP
#define _FROMBIN_CONVERTER_CPP

// ### FromBinConverter.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 07/12/2004: Last Modified
//
// ###

#include "FromBinConverter.hpp"

FromBinConverter::FromBinConverter() : Converter() {
  numHedges = 0;
  numPins = 0;
}

FromBinConverter::~FromBinConverter() {}

void FromBinConverter::readPreamble(ifstream &in_stream) {
  int inPreamble[3];

  in_stream.seekg(0, ifstream::beg);
  in_stream.read((char *)(&inPreamble[0]), sizeof(int) * 3);

  numVerts = inPreamble[0];
  numHedges = inPreamble[1];
  numPins = inPreamble[2];

  vWeights.reserve(numVerts);
}

void FromBinConverter::readInVertexWts(ifstream &in_stream) {
  int length;
  int offset;
  int vertChunk = numVerts * sizeof(int);

  in_stream.seekg(0, ifstream::end);

  length = in_stream.tellg();
  offset = length - vertChunk;

  in_stream.seekg(offset, ifstream::beg);
  in_stream.read((char *)(vWeights.data()), numVerts * sizeof(int));
  in_stream.seekg(3 * sizeof(int), ifstream::beg);
}

void FromBinConverter::readInHedgeData(ifstream &in_stream, int &inStream) {
  in_stream.read((char *)(&dataLength), sizeof(int));
  hEdgeData.reserve(dataLength);
  in_stream.read((char *)(hEdgeData.data()), sizeof(int) * dataLength);
  inStream = in_stream.tellg();
}

#endif
