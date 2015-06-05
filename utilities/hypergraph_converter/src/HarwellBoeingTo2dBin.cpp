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

  FastDynaArray<int> hEdgeData;
  FastDynaArray<int> hEdgeOffsets;
  FastDynaArray<int> pinList;
  FastDynaArray<int> vWeights;

  FastDynaArray<int> rowWts;
  FastDynaArray<int> transposed_cols;
  FastDynaArray<int> transposed_col_offsets;

  int numZeroWeightVerts;
  FastDynaArray<int> zeroWeightVerts;

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

  pinList.setLength(numNonZeros);
  hEdgeOffsets.setLength(numCols + 1);

  j = 0;
  for (i = 0; i < totColPtrLines; ++i) {
    getLine(in_stream);
    data = buffer.getArray();
#ifdef PRINT_GRAPH
    cout << data << endl;
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
    data = buffer.getArray();

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

  BitField addNonZero(numCols);
  addNonZero.clear();

  rowWts.setLength(numCols);
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
      addNonZero.set1(i);
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
  o_stream.write((char *)(rowWts.getArray()), sizeof(int) * numCols);
  o_stream.close();

  rowWts.setLength(0);

  cout << "num zeros on main diagonal = " << numZerosOnMainDiag << endl;

  int num2Dvertices = numNonZeros + numZerosOnMainDiag;
  int num2Dhedges = numCols + numRows;
  int num2Dpins = 2 * num2Dvertices;
  int vertIdx2D = 0;
  int pin;

  BitField hasZeroWeight(num2Dvertices);
  hasZeroWeight.clear();

  FastDynaArray<int> convertedCols;
  FastDynaArray<int> convertedColOffsets(numCols + 1);

  FastDynaArray<int> vertSet;
  FastDynaArray<int> vertSetOffsets(numCols + 1);

  convertedColOffsets[0] = 0;
  vertSetOffsets[0] = 0;

  numZeroWeightVerts = 0;
  zeroWeightVerts.setLength(numZerosOnMainDiag);

  for (i = 0; i < numCols; ++i) {
    int startIdx = hEdgeOffsets[i];
    int endIdx = hEdgeOffsets[i + 1];

    qsort(pinList.getArray(), startIdx, endIdx - 1);
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
        addNonZero.set0(i);
        zeroWeightVerts[numZeroWeightVerts++] = vertIdx2D;
        hasZeroWeight.set1(vertIdx2D++);
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
      hasZeroWeight.set1(vertIdx2D++);
      addNonZero.set0(i);
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
      out_stream.write((char *)(hEdgeData.getArray()), sizeof(int) * ij);
      ij = 0;
    }
  }

  convertedColOffsets.setLength(0);
  convertedCols.setLength(0);
  hEdgeOffsets.setLength(0);
  pinList.setLength(0);

  /* transpose the vertices to get the row hyperedges */

  FastDynaArray<FastDynaArray<int> *> rows(numRows);
  FastDynaArray<int> rowLens(numRows);

  for (i = 0; i < numRows; ++i) {
    rows[i] = new FastDynaArray<int>;
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
      out_stream.write((char *)(hEdgeData.getArray()), sizeof(int) * ij);
      ij = 0;
    }
  }

  for (i = 0; i < numRows; ++i)
    if (rows[i])
      delete rows[i];

  vWeights.setLength(num2Dvertices);

  for (i = 0; i < num2Dvertices; ++i) {
    if (hasZeroWeight(i) == 1)
      vWeights[i] = 0;
    else
      vWeights[i] = 1;
#ifdef PRINT_GRAPH
    cout << vWeights[i] << endl;
#endif
  }

  out_stream.write((char *)(vWeights.getArray()), sizeof(int) * num2Dvertices);

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
  out_stream.write((char *)(zeroWeightVerts.getArray()),
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
