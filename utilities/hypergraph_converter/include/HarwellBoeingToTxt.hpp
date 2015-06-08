#ifndef _HARWELLBOEING_TOTXT_HPP
#define _HARWELLBOEING_TOTXT_HPP

// ### HarwellBoeingToTxt.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 16/4/2005: Last Modified
//
// NOTE: converts the matrix into text format
//
// ###

#include "HarwellBoeingReader.hpp"
#include "data_structures/Bit.hpp"

using namespace std;

class HarwellBoeingToTxt : public HarwellBoeingReader {

protected:
public:
  HarwellBoeingToTxt();
  ~HarwellBoeingToTxt();

  void convert(const char *filename);
};

#endif
