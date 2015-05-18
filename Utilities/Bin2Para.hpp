#  ifndef _BIN2PARA_HPP
#  define _BIN2PARA_HPP


// ### Bin2Para.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
// 
// HISTORY: 
// 
// 02/12/2004: Last Modified
//
// ###


#  include <fstream>
#  include <cstdio>
#  include "FromBinConverter.hpp"


using namespace std;


class Bin2Para
  : public FromBinConverter
{

protected:
 
  int numProcessors;
  int numLocHedges;
  int numLocVertices;

public:
  
  Bin2Para(int numP);
  Bin2Para();
  ~Bin2Para();

  void convert(const char* filename);
  void buildParaFile(ifstream &in_stream, ofstream &out_stream, const char *p_file, int &inLoc, int &inData, int minVertIdx); 
};






#  endif
