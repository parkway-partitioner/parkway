#ifndef _BIN2PARA_CPP
#define _BIN2PARA_CPP

// ### Bin2Para.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 02/12/2004: Last Modified
//
// ###

#include "Bin2Para.hpp"

Bin2Para::Bin2Para(int numP) : FromBinConverter() {
  numProcessors = numP;
  numHedges = 0;
  numLocHedges = 0;
  numLocVertices = 0;
}

Bin2Para::Bin2Para() : FromBinConverter() {
  numProcessors = 0;
  numHedges = 0;
  numLocHedges = 0;
  numLocVertices = 0;
}

Bin2Para::~Bin2Para() {}

void Bin2Para::convert(const char *filename) {
  ifstream in_stream;
  ofstream out_stream;

  char para_file[512];

  int i;
  int inLoc;
  int inData;
  int numHedgesPerProc;
  int numVertsPerProc;
  int minVertexIndex;

  in_stream.open(filename, ifstream::in | ifstream::binary);

  if (!in_stream.is_open()) {
    cout << "error opening " << filename << endl;
    in_stream.close();
    exit(1);
  }

  readPreamble(in_stream);

  numHedgesPerProc = numHedges / numProcessors;
  numVertsPerProc = numVerts / numProcessors;

  readInVertexWts(in_stream);

  inLoc = in_stream.tellg();
  inData = 0;

  for (i = 0; i < numProcessors; ++i) {
    sprintf(para_file, "%s-%d", filename, i);

    minVertexIndex = numVertsPerProc * i;

    if (i == numProcessors - 1) {
      numLocVertices = numVertsPerProc + numVerts % numProcessors;
      numLocHedges = numHedgesPerProc + numHedges % numProcessors;
    } else {
      numLocVertices = numVertsPerProc;
      numLocHedges = numHedgesPerProc;
    }

    buildParaFile(in_stream, out_stream, para_file, inLoc, inData,
                  minVertexIndex);
  }

  in_stream.close();
}

void Bin2Para::buildParaFile(ifstream &in_stream, ofstream &out_stream,
                             const char *p_file, int &inStream, int &inData,
                             int minVertIdx) {
  int i;

  int pin;
  int chunkLength;
  int outLen;
  int numReadHedges;

  dynamic_array<int> outHedgeData;

  out_stream.open(p_file, ofstream::out | ofstream::binary);

  if (!out_stream.is_open()) {
    cout << "error opening " << p_file << endl;
    in_stream.close();
    exit(1);
  }

  out_stream.write((char *)(&numVerts), sizeof(int));
  out_stream.write((char *)(&numLocVertices), sizeof(int));

  numReadHedges = 0;
  outLen = 0;

  for (; numReadHedges < numLocHedges;) {
    if (inData == 0)
      readInHedgeData(in_stream, inStream);

    for (; inData < dataLength;) {
      chunkLength = hEdgeData[inData];

      outHedgeData.assign(outLen++, chunkLength);
      outHedgeData.assign(outLen++, hEdgeData[inData + 1]);

      for (i = 2; i < chunkLength; ++i) {
        pin = hEdgeData[inData + i];

        if (pin < 0 || pin >= numVerts) {
          cout << "pin = " << pin << ", numVerts = " << numVerts << endl;
          exit(1);
        }

        outHedgeData.assign(outLen++, pin);
      }

      inData += chunkLength;
      ++numReadHedges;

      if (numReadHedges == numLocHedges)
        break;
    }

    if (inData == dataLength)
      inData = 0;
  }

  out_stream.write((char *)(&outLen), sizeof(int));
  out_stream.write((char *)(&vWeights[minVertIdx]),
                   sizeof(int) * numLocVertices);
  out_stream.write((char *)(outHedgeData.data()), sizeof(int) * outLen);

  out_stream.close();
}

#endif
