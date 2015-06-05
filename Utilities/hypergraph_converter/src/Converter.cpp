#ifndef _CONVERTER_CPP
#define _CONVERTER_CPP

// ### Converter.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 01/12/2004: Last Modified
//
// ###

#include "Converter.hpp"

Converter::Converter() {
  locOfLengthParam = 0;
  dataLength = 0;
  numVerts = 0;
}

Converter::~Converter() {}

void Converter::setLengthParameter() {
  hEdgeData[locOfLengthParam] = dataLength - 1;
}

void Converter::resetConverterParameters() {
  locOfLengthParam = 0;
  numVerts = 0;
  dataLength = 0;
}

void Converter::addLengthParameter() {
  locOfLengthParam = dataLength;
  hEdgeData.assign(dataLength++, -1);
}

#endif
