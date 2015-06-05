#ifndef _CONVERT_CPP
#define _CONVERT_CPP

// ### Convert.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 12/1/2005: Last Modified
//
// ###

#include "Sandia2Bin.hpp"
#include "HMeTiS2Bin.hpp"
#include "PaToH2Bin.hpp"
#include "Bin2Para.hpp"
#include "Bin2HMeTiS.hpp"
#include "Bin2PaToH.hpp"
#include "Bin2Sandia.hpp"
#include "HarwellBoeingToBin.hpp"
#include "HarwellBoeingTo2dBin.hpp"
#include "HarwellBoeingToTxt.hpp"
#include "MatrixMarket2Bin.hpp"

int main(int argc, char **argv) {
  if (argc < 4) {
    cout
        << "\nUSAGE:" << endl
        << endl
        << "% " << argv[0] << " [convert_options...] <hypergraph filename>"
        << endl
        << endl
        << "PARAMETERS:" << endl
        << endl
        << "\t -form <conversion format>" << endl
        << "\t    - integer specifying the type of file conversion. Takes "
           "following" << endl
        << "\t      possible values:" << endl
        << "\t      1 -> hmetis2bin, converts from hmetis format to binary "
           "format" << endl
        << "\t      2 -> patoh2bin, converts from patoh format to binary format"
        << endl
        << "\t      3 -> sandia2bin, converts from sandia format to binary "
           "format" << endl
        << "\t      4 -> converts from Harwell-Boeing sparse-matrix format to "
           "binary format" << endl
        << "\t      5 -> bin2hmetis, converts from binary to hmetis format"
        << endl
        << "\t      6 -> bin2patoh, converts from binary to patoh format"
        << endl
        << "\t      7 -> bin2sandia, converts from binary to sandia format"
        << endl
        << "\t      8 -> MatrixMarket to Bin" << endl
        << "\t      9 -> bin2para, converts from binary to a number of binary "
           "files" << endl
        << "\t           for use with parkway, where the number of files is "
           "specified" << endl
        << "\t           by the -form option" << endl
        << "\t      10 -> harwell boeing to 2D hypergraph representation"
        << endl
        << "\t      11 -> harwell boeing to text matrix  representation" << endl
        << "\t -np <number of processes>" << endl
        << "\t    - integer specifying the number of files when converting a "
           "hypergraph" << endl
        << "\t      from the binary format to the <number of processes> files "
           "used by" << endl
        << "\t      <number of processes> processes when running parkway"
        << endl
        << "\t -matrix <1 or 0>" << endl
        << "\t    - signifies if the file represents a matrix or not. If it "
           "does not" << endl
        << "\t      reresent a matrix, then it is read in and processed "
           "without modification" << endl
        << "\t      except that singleton hyperedges might be removed. If it "
           "does represent a" << endl
        << "\t      matrix for partitioning, non-zeros are added to main "
           "diagonal elements" << endl
        << "\t      (if not already there) and vertex weights are set to be "
           "the sum of the" << endl
        << "\t      number of non-zero elements in the row (a row-wise "
           "decomposition)" << endl
        << endl;

    exit(1);
  }

  int code;
  int numP;
  int mtx;

  code = StringUtils::getParameterAsInteger(argc, argv, "-form");
  mtx = StringUtils::getParameterAsInteger(argc, argv, "-matrix");

  if (code == -1) {
    cout << code << " incorrect value for option -form" << endl;
    exit(1);
  }

  if (mtx != 0 && mtx != 1) {
    cout << code << " incorrect value for option -matrix (must be 0 or 1)"
         << endl;
    // exit(1);
  }

  numP = StringUtils::getParameterAsInteger(argc, argv, "-np");

  if (code == 1) {
    HMeTiS2Bin converter;
    converter.convert(argv[argc - 1]);
  }

  if (code == 2) {
    PaToH2Bin converter;
    converter.convert(argv[argc - 1]);
  }

  if (code == 3) {
    Sandia2Bin converter;
    converter.convert(argv[argc - 1]);
  }

  if (code == 4) {
    /* for now, by default weight of vertex == 1 */

    HarwellBoeingToBin converter(mtx);
    converter.convert(argv[argc - 1]);
  }

  if (code == 5) {
    Bin2HMeTiS converter;
    converter.convert(argv[argc - 1]);
  }

  if (code == 6) {
    Bin2PaToH converter;
    converter.convert(argv[argc - 1]);
  }

  if (code == 7) {
    Bin2Sandia converter;
    converter.convert(argv[argc - 1]);
  }

  if (code == 8) {
    MatrixMarket2Bin converter;
    converter.convert(argv[argc - 1]);
  }

  if (code == 9) {
    if (numP < 2) {
      cout << "number of processors must be greater than 1" << endl;
      exit(1);
    }

    Bin2Para converter(numP);
    converter.convert(argv[argc - 1]);
  }

  if (code == 10) {
    HarwellBoeingTo2dBin converter;
    converter.convert(argv[argc - 1]);
  }

  if (code == 11) {
    HarwellBoeingToTxt converter;
    converter.convert(argv[argc - 1]);
  }

  return 0;
}

#endif
