#ifndef _BIN2SANDIA_CPP
#define _BIN2SANDIA_CPP

// ### Bin2Sandia.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 07/12/2004: Last Modified
//
// ###

#include "Bin2Sandia.hpp"

Bin2Sandia::Bin2Sandia() : FromBinConverter() {}

Bin2Sandia::~Bin2Sandia() {}

void Bin2Sandia::convert(const char *filename) {
  ifstream in_stream;
  ofstream out_stream;

  char sandia_file[512];
  char preAmble[512];

  int i;
  int pin;
  int chunkLength;
  int inData;
  int inSourceFile;
  int numReadHedges;

  in_stream.open(filename, ifstream::in | ifstream::binary);

  if (!in_stream.is_open()) {
    cout << "error opening " << filename << endl;
    in_stream.close();
    exit(1);
  }

  sprintf(sandia_file, "%s.hg", filename);
  out_stream.open(sandia_file, ofstream::out);

  if (!out_stream.is_open()) {
    cout << "error opening " << sandia_file << endl;
    in_stream.close();
    exit(1);
  }

  readPreamble(in_stream);
  readInVertexWts(in_stream);

  sprintf(&preAmble[0], "%s Sandia file representation of %s\n", "%", filename);
  cout << &preAmble[0] << endl;
  out_stream << &preAmble[0];

  sprintf(&preAmble[0], "%d %d %d\n", numHedges, numVerts, 11);
  out_stream << &preAmble[0];

  numReadHedges = 0;
  inData = 0;
  inSourceFile = in_stream.tellg();

  for (; numReadHedges < numHedges;) {
    if (inData == 0)
      readInHedgeData(in_stream, inSourceFile);

    for (; inData < dataLength;) {
      chunkLength = hEdgeData[inData];
      out_stream << hEdgeData[inData + 1] << " ";

      for (i = 2; i < chunkLength; ++i) {
        pin = hEdgeData[inData + i];

        if (pin < 0 || pin >= numVerts) {
          cout << "pin = " << pin << ", numVerts = " << numVerts << endl;
          in_stream.close();
          out_stream.close();
          exit(1);
        }

        out_stream << (pin + 1) << " ";
      }

      out_stream << "\n";

      inData += chunkLength;
      ++numReadHedges;

      if (numReadHedges == numHedges)
        break;
    }

    if (inData == dataLength)
      inData = 0;
  }

  for (i = 0; i < numVerts; ++i)
    out_stream << vWeights[i] << "\n";

  out_stream.close();
  in_stream.close();
}

#endif
