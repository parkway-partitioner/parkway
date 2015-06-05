
#ifndef _PRINT_CHARS_CPP
#define _PRINT_CHARS_CPP

// ### PrintHypergraphChars.cpp ###
//
// Copyright (C) 2005, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 12/1/2005: Last Modified
//
// ###

#include "PrintHypergraphChars.hpp"

int main(int argc, char **argv) {
  if (argc < 2) {
    cout << "\nUSAGE:" << endl
         << endl
         << "% " << argv[0] << " [options] <hypergraph filename>" << endl
         << endl
         << "PARAMETERS:" << endl
         << endl
         << "\t -ptype <type_format>" << endl
         << "\t    - integer specifying how statistics are calculated:" << endl
         << "\t      0 -> calculate the statistics on the hypergraph as on file"
         << endl
         << "\t      1 -> ignore hyperedges of length one" << endl
         << endl
         << "\t -otype <type_format>" << endl
         << "\t    - string specifying the output file. If not specified then"
         << endl
         << "\t      defaults output to screen" << endl
         << "\t -nz <0 or 1>" << endl
         << "\t    - specifies whether to output info about non-zeros on main"
         << endl
         << "\t      diagonal - 1 yes, 0 - no" << endl
         << endl;
    exit(1);
  }

  char *out_file =
      StringUtils::getParameterAsCharPtr(argc, argv, "-otype", NULL);
  int type = StringUtils::getParameterAsInteger(argc, argv, "-ptype", 1);
  int nz = StringUtils::getParameterAsInteger(argc, argv, "-nz", 0);

  ifstream in_stream;
  ostream *output;

  if (out_file) {
    output = new ofstream(out_file, ofstream::app | ofstream::out);

    if (!output->good()) {
      cout << "could not open file " << out_file << " - exit" << endl;
      exit(0);
    }
  } else
    output = &cout;

  HypergraphChars printer(type, nz);

  printer.initStructs(argv[argc - 1]);
  printer.outputChars(*output);

  return 0;
}

HypergraphChars::HypergraphChars(int remove, int nZ) {
  zerosOnDiag = nZ;
  removeSingletons = remove;
  numVertices = 0;
  numHedges = 0;
  numPins = 0;

  graphName = new char[512];
  strcpy(graphName, "");
}

HypergraphChars::~HypergraphChars() {
  if (graphName)
    delete[] graphName;
}

void HypergraphChars::initStructs(const char *filename) {
  ifstream in_stream;

  register int i;
  register int j;
  register int ij;

  int ijk;
  int preamble[3];
  int hEdgeDataLength;
  int numZerosOnMainDiag;
  int zeroOnDiag;

  FastDynaArray<int> hEdgeData;

  in_stream.open(filename, ifstream::in | ifstream::binary);

  if (!in_stream.is_open()) {
    cout << "error opening " << filename << endl;
    in_stream.close();
    exit(1);
  }

  in_stream.read((char *)(&preamble[0]), sizeof(int) * 3);

  numVertices = preamble[0];
  numPins = preamble[2];

  cout << "AS ON FILE: numVertices = " << numVertices
       << ", numHedges = " << preamble[0] << ", numPins = " << numPins << endl;

  strcpy(graphName, filename);

  vWeights.setLength(numVertices);
  numZerosOnMainDiag = 0;

  if (removeSingletons) {
    numHedges = 0;
    numPins = 0;

    hEdgeWeights.setLength(0);
    hEdgeLens.setLength(0);

    for (i = 0; i < preamble[1];) {
      in_stream.read((char *)(&hEdgeDataLength), sizeof(int));
      hEdgeData.setLength(hEdgeDataLength);
      in_stream.read((char *)(hEdgeData.getArray()),
                     sizeof(int) * hEdgeDataLength);

      for (j = 0; j < hEdgeDataLength;) {
        ij = hEdgeData[j];

        if (ij - 2 > 1) {
          hEdgeLens.assign(numHedges, ij - 2);
          hEdgeWeights.assign(numHedges++, hEdgeData[j + 1]);
          numPins += (ij - 2);
        }

        zeroOnDiag = 0;
        for (ijk = j + 2; ijk < j + ij; ++ijk)
          if (hEdgeData[ijk] == i)
            zeroOnDiag = 1;

        if (zeroOnDiag == 0)
          ++numZerosOnMainDiag;

        j += ij;
        ++i;
      }
    }
  } else {
    numHedges = preamble[1];

    hEdgeWeights.setLength(numHedges);
    hEdgeLens.setLength(numHedges);

    for (i = 0; i < numHedges;) {
      in_stream.read((char *)(&hEdgeDataLength), sizeof(int));
      hEdgeData.setLength(hEdgeDataLength);
      in_stream.read((char *)(hEdgeData.getArray()),
                     sizeof(int) * hEdgeDataLength);

      for (j = 0; j < hEdgeDataLength;) {
        ij = hEdgeData[j];

        hEdgeLens[i] = ij - 2;
        hEdgeWeights[i++] = hEdgeData[j + 1];

        zeroOnDiag = 0;
        for (ijk = j + 2; ijk < j + ij; ++ijk)
          if (hEdgeData[ijk] == i)
            zeroOnDiag = 1;

        if (zeroOnDiag == 0)
          ++numZerosOnMainDiag;

        j += ij;
      }
    }
  }

  if (zerosOnDiag)
    cout << "num zeros on main diagonal = " << numZerosOnMainDiag << endl;

  in_stream.seekg(0, ifstream::end);
  i = in_stream.tellg();
  ij = i - (numVertices * sizeof(int));

  in_stream.seekg(ij, ifstream::beg);
  in_stream.read((char *)(vWeights.getArray()), numVertices * sizeof(int));

  in_stream.close();
}

void HypergraphChars::outputChars(ostream &out) {
  register int i;
  register int j;
  register int ij;

  if (removeSingletons)
    out << "AFTER REMOVAL OF SINGLETON HYPEREDGES: " << endl;

  out << "Hypergraph " << graphName << endl
      << " #vertices: " << numVertices << " #hyperedges: " << numHedges
      << " #pins: " << numPins << endl;

  double weighted_ave = 0;
  double percentile_95;
  double percentile_90;
  double percentile_75;
  double percentile_50;
  double percentile_25;

  FastDynaArray<int> indices;

  /* display hyperedge information */

  j = 0;
  indices.setLength(numHedges);

  for (i = 0; i < numHedges; ++i) {
    indices[i] = i;
    j += hEdgeWeights[i];
    weighted_ave += (hEdgeWeights[i] * hEdgeLens[i]);
  }

  out << "hyperedge length percentiles: (weighted ave, 25, 50, 75, 95, "
         "maxLength) " << endl;
  out << "\t" << weighted_ave / j << " ";

  percentile_95 = (static_cast<double>(j) * 95) / 100;
  percentile_90 = (static_cast<double>(j) * 90) / 100;
  percentile_75 = (static_cast<double>(j) * 75) / 100;
  percentile_50 = (static_cast<double>(j) * 50) / 100;
  percentile_25 = (static_cast<double>(j) * 25) / 100;

  qsort(0, numHedges - 1, indices.getArray(), hEdgeLens.getArray());

  j = 0;
  i = 0;
  ij = 0;

  for (; i < numHedges;) {
    j += hEdgeWeights[indices[i++]];

    if (ij == 0 && j > percentile_25) {
      out << hEdgeLens[indices[i]] << " ";
      ++ij;
    }

    if (ij == 1 && j > percentile_50) {
      out << hEdgeLens[indices[i]] << " ";
      ++ij;
    }

    if (ij == 2 && j > percentile_75) {
      out << hEdgeLens[indices[i]] << " ";
      ++ij;
    }

    if (ij == 3 && j > percentile_90) {
      out << hEdgeLens[indices[i]] << " ";
      ++ij;
    }

    if (ij == 4 && j > percentile_95) {
      out << hEdgeLens[indices[i]] << " ";
      ++ij;
    }

    if (i == numHedges - 1) {
      out << hEdgeLens[indices[i]] << endl;
    }
  }

  /* display vertex information */

  j = 0;
  indices.setLength(numVertices);

  for (i = 0; i < numVertices; ++i) {
    indices[i] = i;
    j += vWeights[i];
  }

  percentile_95 = (static_cast<double>(j) * 95) / 100;
  percentile_75 = (static_cast<double>(j) * 75) / 100;
  percentile_50 = (static_cast<double>(j) * 50) / 100;
  percentile_25 = (static_cast<double>(j) * 25) / 100;

  qsort(0, numVertices - 1, indices.getArray(), vWeights.getArray());

  out << "vertex weight percentiles: (ave, 25, 50, 75, 95, maxWeight) " << endl;
  out << "\t" << static_cast<double>(j) / numVertices << " ";

  j = 0;
  i = 0;
  ij = 0;

  for (; i < numVertices;) {
    j += vWeights[indices[i++]];

    if (ij == 0 && j > percentile_25) {
      out << vWeights[indices[i]] << " ";
      ++ij;
    }

    if (ij == 1 && j > percentile_50) {
      out << vWeights[indices[i]] << " ";
      ++ij;
    }

    if (ij == 2 && j > percentile_75) {
      out << vWeights[indices[i]] << " ";
      ++ij;
    }

    if (ij == 3 && j > percentile_95) {
      out << vWeights[indices[i]] << " ";
      ++ij;
    }

    if (i == numVertices - 1) {
      out << vWeights[indices[i]] << endl;
    }
  }
}

void HypergraphChars::qsort(const int left, const int right, int *array,
                            const int *cmp) {
  register int left_arrow = left;
  register int right_arrow = right;

  int pivot = array[(left + right) / 2];

  do {
    while (cmp[array[right_arrow]] > cmp[pivot])
      --right_arrow;
    while (cmp[array[left_arrow]] < cmp[pivot])
      ++left_arrow;

    if (left_arrow <= right_arrow) {
      swap(array[left_arrow++], array[right_arrow--]);
    }
  } while (right_arrow >= left_arrow);

  if (left < right_arrow)
    qsort(left, right_arrow, array, cmp);

  if (left_arrow < right)
    qsort(left_arrow, right, array, cmp);
}

#endif
