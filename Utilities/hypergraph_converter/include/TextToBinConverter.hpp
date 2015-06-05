#ifndef _TEXT_TOBIN_CONVERTER_HPP
#define _TEXT_TOBIN_CONVERTER_HPP

// ### TextToBinConverter.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 01/12/2004: Last Modified
//
// ###

#include "Converter.hpp"

using namespace std;

class TextToBinConverter : public Converter {

protected:
public:
  TextToBinConverter();

  virtual ~TextToBinConverter();
  virtual int processHedgeLine(char *line, int &numP) = 0;
  virtual int processVertLine(char *line) = 0;
};

#endif
