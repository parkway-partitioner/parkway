
#ifndef _HARWELLBOEING_READER_CPP
#define _HARWELLBOEING_READER_CPP

// ### HarwellBoeingReader.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 30/12/2004: Last Modified
//
// ###

#include "HarwellBoeingReader.hpp"

HarwellBoeingReader::HarwellBoeingReader() : TextFileReader() {
  title.setLength(73);
  id.setLength(9);
  mtxType.setLength(4);

  totDataLines = 0;
  totColPtrLines = 0;
  totRowIdxLines = 0;
  totNumValLines = 0;
  numRHSLines = 0;
  numRows = 0;
  numCols = 0;
  numNonZeros = 0;
  numElemMatEntries = 0;
}

HarwellBoeingReader::~HarwellBoeingReader() {}

void HarwellBoeingReader::readPreamble(ifstream &in) {
  register int i;
  register int j;

  char *data;

  /* read in 1st line */

  getLine(in);
  // assert(length >= 80);

  i = 0;
  for (; buffer[i] == ' '; ++i)
    ;
  for (j = 0; j < 72; ++j)
    title[j] = buffer[i++];
  title[j] = '\0';

  for (; buffer[i] == ' '; ++i)
    ;
  for (j = 0; j < 8; ++j)
    id[j] = buffer[i++];
  id[j] = '\0';

  /* read in 2nd line */

  getLine(in);
  data = buffer.getArray();

  /* read in num tot data lines */

  StringUtils::skipNonDigits(data);
  totDataLines = StringUtils::stringToDigit(data);

  /* read in number lines for ptrs */

  StringUtils::skipNonDigits(data);
  totColPtrLines = StringUtils::stringToDigit(data);

  /* read in number of lines for row */

  StringUtils::skipNonDigits(data);
  totRowIdxLines = StringUtils::stringToDigit(data);

  /* read in number of lines for numerical vals */

  StringUtils::skipNonDigits(data);
  totNumValLines = StringUtils::stringToDigit(data);

  /* read in number of lines for RHS */

  StringUtils::skipNonDigits(data);
  numRHSLines = StringUtils::stringToDigit(data);

  /* read in 3rd line */

  getLine(in);

  i = 0;
  for (; buffer[i] == ' '; ++i)
    ;
  for (j = 0; j < 3; ++j)
    mtxType[j] = buffer[i++];
  mtxType[j] = '\0';

  data = &buffer[i];

  /* read in number of matrix rows */

  StringUtils::skipNonDigits(data);
  numRows = StringUtils::stringToDigit(data);

  /* read in number of matrix columns */

  StringUtils::skipNonDigits(data);
  numCols = StringUtils::stringToDigit(data);

  /* read in number of non-zeros (in assembled matrix) */

  StringUtils::skipNonDigits(data);
  numNonZeros = StringUtils::stringToDigit(data);

  /* read in number of elemental matrix entries (zero for assembled matrix) */

  StringUtils::skipNonDigits(data);
  numElemMatEntries = StringUtils::stringToDigit(data);

  /* skip remaining pre-amble lines */

  getLine(in);
  if (numRHSLines > 0)
    getLine(in);
}

#endif
