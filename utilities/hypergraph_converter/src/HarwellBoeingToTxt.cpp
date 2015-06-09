#ifndef _HARWELLBOEING_TOTXT_CPP
#define _HARWELLBOEING_TOTXT_CPP

// ### HarwellBoeingToTxt.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 16/4/2004: Last Modified
//
// ###

#include "HarwellBoeingToTxt.hpp"

HarwellBoeingToTxt::HarwellBoeingToTxt() : HarwellBoeingReader() {}

HarwellBoeingToTxt::~HarwellBoeingToTxt() {}

void HarwellBoeingToTxt::convert(const char *filename) {
  /* assume that this is a general, non-symmetric sparse matrix */

  ifstream in_stream;
  ofstream out_stream;

  char txt_file[512];
  char *data;

  dynamic_array<int> hEdgeOffsets;
  dynamic_array<int> pinList;

  int i;
  int j;

  in_stream.open(filename, ifstream::in);

  if (!in_stream.is_open()) {
    cout << "error opening " << filename << endl;
    exit(1);
  }

  sprintf(txt_file, "%s.txt", filename);

  out_stream.open(txt_file, ofstream::out | ofstream::app);

  if (!out_stream.is_open()) {
    cout << "error opening " << txt_file << endl;
    in_stream.close();
    exit(1);
  }

  cout << "opened file " << txt_file << endl;

  readPreamble(in_stream);

  cout << "reading in graph: " << endl
       << "v = " << numRows << endl
       << "h = " << numCols << endl
       << "nz = " << numNonZeros << endl;

  out_stream << numCols << " " << numRows << " " << numNonZeros << endl;

  hEdgeOffsets.reserve(numCols + 1);
  pinList.reserve(numNonZeros);

  j = 0;
  for (i = 0; i < totColPtrLines; ++i) {

    getLine(in_stream);
    data = buffer.data();
#ifdef PRINT_GRAPH
    cout << data_ << endl;
#endif
    while (*data != '\0') {
      StringUtils::skipNonDigits(data);
      hEdgeOffsets[j++] = StringUtils::stringToDigit(data) - 1;
      assert(hEdgeOffsets[j - 1] <= numNonZeros);
    }
  }

  if (j != numCols + 1) {
    cout << "could not find " << numCols + 1 << " col ptr entries" << endl;
    cout << "instead found " << j << " col ptr entries" << endl;

    in_stream.close();
    out_stream.close();
    exit(1);
  }

  assert(hEdgeOffsets[j - 1] == numNonZeros);

  j = 0;
  for (i = 0; i < totRowIdxLines; ++i) {
    getLine(in_stream);
    data = buffer.data();

    while (*data != '\0') {
      StringUtils::skipNonDigits(data);
      if (j >= numNonZeros) {
        cout << "j = " << j << endl << "nz = " << numNonZeros << endl;
        assert(0);
      }

      pinList[j++] = StringUtils::stringToDigit(data) - 1;
    }
  }

  if (j != numNonZeros) {
    cout << "could not find " << numNonZeros << " non-zero entries" << endl;
    cout << "instead found " << j << " non-zero entries" << endl;

    in_stream.close();
    out_stream.close();
    exit(1);
  }

  for (i = 0; i < numCols; ++i) {

    int startOffset = hEdgeOffsets[i];
    int endOffset = hEdgeOffsets[i + 1];

    for (j = startOffset; j < endOffset; ++j) {
      if (j == startOffset)
        out_stream << pinList[j];
      else
        out_stream << " " << pinList[j];
    }

    out_stream << endl;
  }

  in_stream.close();
  out_stream.close();
}

#endif
