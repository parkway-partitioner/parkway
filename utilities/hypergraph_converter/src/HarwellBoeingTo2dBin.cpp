#ifndef _HARWELLBOEING_TO2dBIN_CPP
#define _HARWELLBOEING_TO2dBIN_CPP

// ### HarwellBoeingToBin.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 13/4/2004: Last Modified
//
// ###

#include "HarwellBoeingTo2dBin.hpp"
#include "data_structures/bit_field.hpp"

using namespace parkway::data_structures;

HarwellBoeingTo2dBin::HarwellBoeingTo2dBin() : HarwellBoeingReader() {}

HarwellBoeingTo2dBin::~HarwellBoeingTo2dBin() {}

void HarwellBoeingTo2dBin::convert(const char *filename) {
  /* assume that this is a general, non-symmetric sparse matrix */

  ifstream in_stream;
  ofstream out_stream;

  char bin_file[512];
  char zeroVertsFile[512];
  char rowWtsFile[512];
  char *data;

  int binPreAmble[3];
  int hEdgeLengthSlot;
  int maxIntsWritten = TextFileReader::getMaxPinsInChunk();

  dynamic_array<int> hEdgeData;
  dynamic_array<int> hEdgeOffsets;
  dynamic_array<int> pinList;
  dynamic_array<int> vWeights;

  dynamic_array<int> rowWts;
  dynamic_array<int> transposed_cols;
  dynamic_array<int> transposed_col_offsets;

  int numZeroWeightVerts;
  dynamic_array<int> zeroWeightVerts;

  int i;
  int j;
  int ij;

  in_stream.open(filename, ifstream::in);

  if (!in_stream.is_open()) {
    cout << "error opening " << filename << endl;
    exit(1);
  }

  sprintf(bin_file, "%s.2d.bin", filename);

  out_stream.open(bin_file, ofstream::out | ofstream::app | ofstream::binary);

  if (!out_stream.is_open()) {
    cout << "error opening " << bin_file << endl;
    in_stream.close();
    exit(1);
  }

  readPreamble(in_stream);

  cout << "reading in graph: " << endl
       << "v = " << numRows << endl
       << "h = " << numCols << endl
       << "nz = " << numNonZeros << endl;

  pinList.reserve(numNonZeros);
  hEdgeOffsets.reserve(numCols + 1);

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

  cout << "hEdgeOffsets[" << j - 1 << "] = " << hEdgeOffsets[j - 1] << endl;

  assert(hEdgeOffsets[j - 1] == numNonZeros);

/* Initialise the pin-list */
#ifdef PRINT_GRAPH
  cout << "row idx lines = " << totRowIdxLines << endl;
#endif

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
#ifdef PRINT_GRAPH
      cout << pinList[j - 1] << " ";
#endif
    }
#ifdef PRINT_GRAPH
    cout << endl;
#endif
  }

  if (j != numNonZeros) {
    cout << "could not find " << numNonZeros << " non-zero entries" << endl;
    cout << "instead found " << j << " non-zero entries" << endl;

    in_stream.close();
    out_stream.close();
    exit(1);
  }

  /* check num zeros on main diagonal */

  int numZerosOnMainDiag = 0;
  int hasNonZeroOnDiag;

  bit_field addNonZero(numCols);
  addNonZero.unset();

  rowWts.reserve(numCols);
  for (i = 0; i < numCols; ++i)
    rowWts[i] = 0;

  for (i = 0; i < numCols; ++i) {
    hasNonZeroOnDiag = 0;

    for (j = hEdgeOffsets[i]; j < hEdgeOffsets[i + 1]; ++j) {
      if (pinList[j] == i)
        hasNonZeroOnDiag = 1;
      rowWts[pinList[j]]++;
    }

    if (hasNonZeroOnDiag == 0) {
      ++numZerosOnMainDiag;
      addNonZero.set(i);
    }
  }

  ofstream o_stream;
  sprintf(rowWtsFile, "%s_num_nz_perrow", filename);

  o_stream.open(rowWtsFile, ofstream::out | ofstream::app | ofstream::binary);

  if (!o_stream.is_open()) {
    cout << "error opening " << rowWtsFile << endl;
    o_stream.close();
    exit(1);
  }

  o_stream.write((char *)(&numCols), sizeof(int));
  o_stream.write((char *)(rowWts.data()), sizeof(int) * numCols);
  o_stream.close();

  rowWts.reserve(0);

  cout << "num zeros on main diagonal = " << numZerosOnMainDiag << endl;

  int num2Dvertices = numNonZeros + numZerosOnMainDiag;
  int num2Dhedges = numCols + numRows;
  int num2Dpins = 2 * num2Dvertices;
  int vertIdx2D = 0;
  int pin;

  bit_field hasZeroWeight(num2Dvertices);
  hasZeroWeight.unset();

  dynamic_array<int> convertedCols;
  dynamic_array<int> convertedColOffsets(numCols + 1);

  dynamic_array<int> vertSet;
  dynamic_array<int> vertSetOffsets(numCols + 1);

  convertedColOffsets[0] = 0;
  vertSetOffsets[0] = 0;

  numZeroWeightVerts = 0;
  zeroWeightVerts.reserve(numZerosOnMainDiag);

  for (i = 0; i < numCols; ++i) {
    int startIdx = hEdgeOffsets[i];
    int endIdx = hEdgeOffsets[i + 1];

    qsort(pinList.data(), startIdx, endIdx - 1);
    convertedColOffsets[i + 1] = convertedColOffsets[i];
    vertSetOffsets[i + 1] = vertSetOffsets[i];

    for (j = startIdx; j < endIdx; ++j) {
      pin = pinList[j];
      assert(pin >= 0 && pin < numRows);

      if (addNonZero(i) == 1 && pin > i) {

        assert(vertIdx2D < num2Dvertices);
#ifdef PRINT_GRAPH
        cout << vertIdx2D << " ";
#endif
        convertedCols.assign(convertedColOffsets[i + 1]++, vertIdx2D);
        vertSet.assign(vertSetOffsets[i + 1]++, i);
        vertSet.assign(vertSetOffsets[i + 1]++, vertIdx2D);
        addNonZero.unset(i);
        zeroWeightVerts[numZeroWeightVerts++] = vertIdx2D;
        hasZeroWeight.set(vertIdx2D++);
      }
#ifdef PRINT_GRAPH
      cout << vertIdx2D << " ";
#endif
      convertedCols.assign(convertedColOffsets[i + 1]++, vertIdx2D);
      vertSet.assign(vertSetOffsets[i + 1]++, pin);
      vertSet.assign(vertSetOffsets[i + 1]++, vertIdx2D++);
    }

    if (addNonZero(i) == 1) {

      assert(vertIdx2D < num2Dvertices);
#ifdef PRINT_GRAPH
      cout << vertIdx2D << " ";
#endif
      convertedCols.assign(convertedColOffsets[i + 1]++, vertIdx2D);
      vertSet.assign(vertSetOffsets[i + 1]++, i);
      vertSet.assign(vertSetOffsets[i + 1]++, vertIdx2D);
      zeroWeightVerts[numZeroWeightVerts++] = vertIdx2D;
      hasZeroWeight.set(vertIdx2D++);
      addNonZero.unset(i);
    }
#ifdef PRINT_GRAPH
    cout << endl;
#endif
  }

  assert(vertIdx2D == num2Dvertices);
  assert(numZeroWeightVerts == numZerosOnMainDiag);

  /* now write out first half of 2d pin list */

  binPreAmble[0] = num2Dvertices;
  binPreAmble[1] = num2Dhedges;
  binPreAmble[2] = num2Dpins;

  cout << "numVertices = " << num2Dvertices << ", numHedges = " << num2Dhedges
       << ", numPins = " << binPreAmble[2] << endl;

  out_stream.write((char *)(&binPreAmble), sizeof(int) * 3);

  ij = 0;
  for (i = 0; i < numCols; ++i) {
    if (ij == 0)
      hEdgeData.assign(ij++, -1);

    hEdgeLengthSlot = ij;

    int startIdx = convertedColOffsets[i];
    int endIdx = convertedColOffsets[i + 1];

    hEdgeData.assign(ij++, -1);
    hEdgeData.assign(ij++, 1);

    for (j = startIdx; j < endIdx; ++j)
      hEdgeData.assign(ij++, convertedCols[j]);

    hEdgeData[hEdgeLengthSlot] = ij - hEdgeLengthSlot;

    if (ij > maxIntsWritten || i == numCols - 1) {
      hEdgeData[0] = ij - 1;
      out_stream.write((char *)(hEdgeData.data()), sizeof(int) * ij);
      ij = 0;
    }
  }

  convertedColOffsets.reserve(0);
  convertedCols.reserve(0);
  hEdgeOffsets.reserve(0);
  pinList.reserve(0);

  /* transpose the vertices to get the row hyperedges */

  dynamic_array<dynamic_array<int> *> rows(numRows);
  dynamic_array<int> rowLens(numRows);

  for (i = 0; i < numRows; ++i) {
    rows[i] = new dynamic_array<int>;
    rowLens[i] = 0;
  }

  for (i = 0; i < numCols; ++i) {

    int startIdx = vertSetOffsets[i];
    int endIdx = vertSetOffsets[i + 1];

    assert(((endIdx - startIdx) & 0x1) == 0);

    for (ij = startIdx; ij < endIdx; ij += 2) {

      int row = vertSet[ij];
      int vertexID = vertSet[ij + 1];

      rows[row]->assign(rowLens[row]++, vertexID);
    }
  }

  /* write out second half of 2d pin list */

  ij = 0;
  for (i = 0; i < numRows; ++i) {
    if (ij == 0)
      hEdgeData.assign(ij++, -1);

    hEdgeLengthSlot = ij;

    int endIdx = rowLens[i];

    hEdgeData.assign(ij++, -1);
    hEdgeData.assign(ij++, 1);

    for (j = 0; j < endIdx; ++j) {
#ifdef PRINT_GRAPH
      cout << (*rows[i])[j] << " ";
#endif
      hEdgeData.assign(ij++, (*rows[i])[j]);
    }
#ifdef PRINT_GRAPH
    cout << endl;
#endif
    hEdgeData[hEdgeLengthSlot] = ij - hEdgeLengthSlot;

    if (ij > maxIntsWritten || i == numCols - 1) {
      hEdgeData[0] = ij - 1;
      out_stream.write((char *)(hEdgeData.data()), sizeof(int) * ij);
      ij = 0;
    }
  }

  for (i = 0; i < numRows; ++i)
    if (rows[i])
      delete rows[i];

  vWeights.reserve(num2Dvertices);

  for (i = 0; i < num2Dvertices; ++i) {
    if (hasZeroWeight(i) == 1)
      vWeights[i] = 0;
    else
      vWeights[i] = 1;
#ifdef PRINT_GRAPH
    cout << vWeights[i] << endl;
#endif
  }

  out_stream.write((char *)(vWeights.data()), sizeof(int) * num2Dvertices);

  in_stream.close();
  out_stream.close();

  /* write out zero weight verts */

  sprintf(zeroVertsFile, "%s_zero_wts_verts", bin_file);

  out_stream.open(zeroVertsFile,
                  ofstream::out | ofstream::app | ofstream::binary);

  if (!out_stream.is_open()) {
    cout << "error opening " << zeroVertsFile << endl;
    exit(1);
  }

  out_stream.write((char *)(&numZeroWeightVerts), sizeof(int));
  out_stream.write((char *)(zeroWeightVerts.data()),
                   sizeof(int) * numZeroWeightVerts);

  out_stream.close();
}

void HarwellBoeingTo2dBin::qsort(int *array, int left, int right) {
  int left_arrow = left;
  int right_arrow = right;

  int pivot = array[(left + right) / 2];

  do {
    while (array[right_arrow] > pivot)
      --right_arrow;
    while (array[left_arrow] < pivot)
      ++left_arrow;

    if (left_arrow <= right_arrow) {
      swap(array[left_arrow], array[right_arrow]);
      ++left_arrow;
      --right_arrow;
    }
  } while (right_arrow >= left_arrow);

  if (left < right_arrow)
    qsort(array, left, right_arrow);

  if (left_arrow < right)
    qsort(array, left_arrow, right);
}

#endif
