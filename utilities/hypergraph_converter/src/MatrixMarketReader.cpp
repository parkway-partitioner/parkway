
#ifndef _MTXMKT_READER_CPP
#define _MTXMKT_READER_CPP

// ### MatrixMarketReader.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 12/1/2004: Last Modified
//
// ###

#include "MatrixMarketReader.hpp"

MatrixMarketReader::MatrixMarketReader() : TextFileReader() {
  numVertices = 0;
  numHyperedges = 0;
  numPins = 0;

  hEdgeOffsets.setLength(0);
  pinList.setLength(0);
}

MatrixMarketReader::~MatrixMarketReader() {}

void MatrixMarketReader::readMatrix(const char *filename) {
  ifstream in_stream;

  in_stream.open(filename, ifstream::in);

  if (!in_stream.is_open()) {
    cout << "error opening " << filename << endl;
    exit(1);
  }

  readPreamble(in_stream);
  readMatrix(in_stream);

  in_stream.close();
}

void MatrixMarketReader::readPreamble(ifstream &in_stream) {
  cout << "read preamble - only interested in non-zero pattern" << endl;

  skipComments(in_stream);

  char *data = buffer.getArray();

  StringUtils::skipNonDigits(data);
  numVertices = StringUtils::stringToDigit(data);

  StringUtils::skipNonDigits(data);
  numHyperedges = StringUtils::stringToDigit(data);

  StringUtils::skipNonDigits(data);
  numPins = StringUtils::stringToDigit(data);
}

void MatrixMarketReader::readMatrix(ifstream &in_stream) {
  int i;
  int j;

  char *data;
  int hEdge;
  int vertex;
  int pointInFile;

  DynamicArray<int> hEdgeLens(numHyperedges);

  pinList.setLength(numPins);
  hEdgeOffsets.setLength(numHyperedges + 1);

  pointInFile = in_stream.tellg();

  for (i = 0; i < numHyperedges; ++i)
    hEdgeLens[i] = 0;

  for (i = 0; i < numPins; ++i) {
    getLine(in_stream);
    data = buffer.getArray();

    StringUtils::skipNonDigits(data);
    vertex = StringUtils::stringToDigit(data) - 1;

    if (vertex < 0 || vertex >= numVertices) {
      cout << "read vertex = " << vertex << ", numVertices = " << numVertices
           << endl;

      in_stream.close();
      exit(1);
    }

    StringUtils::skipNonDigits(data);
    hEdge = StringUtils::stringToDigit(data) - 1;

    if (hEdge < 0 || hEdge >= numHyperedges) {
      cout << "read hyperedge = " << hEdge
           << ", numHyperedges = " << numHyperedges << endl;

      in_stream.close();
      exit(1);
    }

    ++hEdgeLens[hEdge];
  }

  in_stream.seekg(pointInFile, ios::beg);

  i = 0;
  j = 0;

  for (; i < numHyperedges; ++i) {
    hEdgeOffsets[i] = j;
    j += hEdgeLens[i];
    hEdgeLens[i] = 0;
  }

  hEdgeOffsets[i] = j;

  for (i = 0; i < numPins; ++i) {
    getLine(in_stream);
    data = buffer.getArray();

    StringUtils::skipNonDigits(data);
    vertex = StringUtils::stringToDigit(data) - 1;

    StringUtils::skipNonDigits(data);
    hEdge = StringUtils::stringToDigit(data) - 1;

    pinList[hEdgeOffsets[hEdge] + hEdgeLens[hEdge]] = vertex;
    ++hEdgeLens[hEdge];
  }
}

void MatrixMarketReader::skipComments(ifstream &in_stream) {
  char c;

  do {
    getLine(in_stream);
    c = StringUtils::getFirstChar(buffer.getArray());
  } while (c == '%');
}

#endif
