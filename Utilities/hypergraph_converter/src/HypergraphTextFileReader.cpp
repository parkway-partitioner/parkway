#ifndef _HGRAPH_TEXTREADER_CPP
#define _HGRAPH_TEXTREADER_CPP

// ### HypergraphTextFileReader.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 01/12/2004: Last Modified
//
// ###

#include "HypergraphTextFileReader.hpp"

HypergraphTextFileReader::HypergraphTextFileReader() : TextFileReader() {
  numVertices = 0;
  numHedges = 0;
  numPins = 0;
  wtsOnVertices = 0;
  wtsOnHedges = 0;
}

HypergraphTextFileReader::~HypergraphTextFileReader() {}

#endif
