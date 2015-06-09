#ifndef _HMETIS2BIN_CPP
#define _HMETIS2BIN_CPP

// ### HMeTiS2Bin.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 13/12/2004: Last Modified
//
// ###

#include "HMeTiS2Bin.hpp"

HMeTiS2Bin::HMeTiS2Bin() : HypergraphTextFileReader(), TextToBinConverter() {}

HMeTiS2Bin::~HMeTiS2Bin() {}

void HMeTiS2Bin::convert(const char *filename) {
  char bin_file[512];

  int maxPinsRead = TextFileReader::getMaxPinsInChunk();
  int numHedgesRead;
  int numPinsRead;
  int binPreAmble[3];

  ifstream in_stream;
  ofstream out_stream;

  in_stream.open(filename, ifstream::in);

  if (!in_stream.is_open()) {
    cout << "error opening " << filename << endl;
    exit(1);
  }

  sprintf(bin_file, "%s.bin", filename);
  out_stream.open(bin_file, ofstream::out | ofstream::binary);

  if (!out_stream.is_open()) {
    cout << "error opening " << bin_file << endl;
    in_stream.close();
    exit(1);
  }

  readPreamble(in_stream);

  numPins = computeNumPins(in_stream);

  binPreAmble[0] = numVertices;
  binPreAmble[1] = numHedges;
  binPreAmble[2] = numPins;

  out_stream.write((char *)(&binPreAmble), sizeof(int) * 3);

  numPinsRead = 0;
  numHedgesRead = 0;

  while (!in_stream.eof() && numHedgesRead < numHedges) {
    if (numPinsRead == 0) {
      resetConverterParameters();
      addLengthParameter();
    }

    getLine(in_stream);
    numHedgesRead += processHedgeLine(buffer.data(), numPinsRead);

    if (numPinsRead > maxPinsRead || numHedgesRead == numHedges) {
      setLengthParameter();
      out_stream.write((char *)(hEdgeData.data()),
                       sizeof(int) * dataLength);
      numPinsRead = 0;
    }
  }

  if (numHedgesRead < numHedges)
    cout << "warning - could not read " << numHedges << " hyperedges" << endl;

  if (wtsOnVertices) {
    resetConverterParameters();

    while (!in_stream.eof() && numVerts < numVertices) {
      getLine(in_stream);
      processVertLine(buffer.data());
    }

    if (numVerts < numVertices) {
      cout << "warning - do not have all vertex weights" << endl;
      in_stream.close();
      out_stream.close();
      exit(1);
    }
  } else {
    vWeights.reserve(numVertices);
    numVerts = numVertices;

    for (int i = 0; i < numVertices; ++i)
      vWeights[i] = 1;
  }

  out_stream.write((char *)(vWeights.data()), sizeof(int) * numVertices);

  in_stream.close();
  out_stream.close();
}

void HMeTiS2Bin::readPreamble(ifstream &in) {
  char *line;
  int options = 0;

  getLine(in);
  line = buffer.data();
  StringUtils::skipNonDigits(line, '%');

  while (*line == '\0' && !in.eof()) {
    getLine(in);
    line = buffer.data();
    StringUtils::skipNonDigits(line, '%');
  }

  if (in.eof()) {
    cout << "warning - cannot find preamble" << endl;
    exit(1);
  }

  // ###
  // now read in 3 numbers

  StringUtils::skipNonDigits(line, '%');
  numHedges = StringUtils::stringToDigit(line);
  StringUtils::skipNonDigits(line, '%');
  numVertices = StringUtils::stringToDigit(line);
  StringUtils::skipNonDigits(line, '%');

  if (*line)
    options = StringUtils::stringToDigit(line);

  // ###
  // process options

  switch (options) {
  case 0:
    wtsOnVertices = 0;
    wtsOnHedges = 0;
    break;
  case 1:
    wtsOnVertices = 0;
    wtsOnHedges = 1;
    break;
  case 10:
    wtsOnVertices = 1;
    wtsOnHedges = 0;
    break;
  case 11:
    wtsOnVertices = 1;
    wtsOnHedges = 1;
    break;
  default:
    wtsOnVertices = 0;
    wtsOnHedges = 0;
    break;
  }
}

int HMeTiS2Bin::processHedgeLine(char *line, int &numP) {
  int hEdgeLengthSlot = dataLength;
  int pin;

  StringUtils::skipNonDigits(line, '%');

  if (*line) {
    hEdgeData.assign(dataLength++, 0);

    if (wtsOnHedges)
      hEdgeData.assign(dataLength++, StringUtils::stringToDigit(line));
    else
      hEdgeData.assign(dataLength++, 1);

    if (!*line)
      cout << "warning - hyperedge with no vertices" << endl;

    while (*line) {
      StringUtils::skipNonDigits(line, '%');

      if (*line) {
        pin = StringUtils::stringToDigit(line) - 1;

        if (pin < 0 || pin >= numVertices) {
          cout << "pin = " << pin << ", numVertices = " << numVertices << endl;
          exit(1);
        }

        hEdgeData.assign(dataLength++, pin);
        ++numP;
      }
    }

    hEdgeData[hEdgeLengthSlot] = dataLength - hEdgeLengthSlot;

    return 1;
  } else
    return 0;
}

int HMeTiS2Bin::countPinsInLine(char *line, int &nP) {
  StringUtils::skipNonDigits(line, '%');

  if (*line) {
    if (wtsOnHedges) {
      StringUtils::skipDigits(line);
      StringUtils::skipNonDigits(line, '%');
    }

    while (*line) {
      StringUtils::skipDigits(line);
      StringUtils::skipNonDigits(line, '%');
      ++nP;
    }

    return 1;
  } else
    return 0;
}

int HMeTiS2Bin::processVertLine(char *line) {
  StringUtils::skipNonDigits(line, '%');

  if (*line) {
    vWeights.assign(numVerts++, StringUtils::stringToDigit(line));
    return 1;
  } else
    return 0;
}

int HMeTiS2Bin::computeNumPins(ifstream &in_stream) {
  int loc = in_stream.tellg();
  int numPins = 0;
  int numHedgesRead = 0;

  while (!in_stream.eof() && numHedgesRead < numHedges) {
    getLine(in_stream);
    numHedgesRead += countPinsInLine(buffer.data(), numPins);
  }

  in_stream.seekg(loc, ios::beg);

  return numPins;
}

#endif
