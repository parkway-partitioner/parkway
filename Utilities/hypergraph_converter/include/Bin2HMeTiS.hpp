#  ifndef _BIN2HMETIS_HPP
#  define _BIN2HMETIS_HPP


// ### Bin2HMeTiS.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
// 
// HISTORY: 
// 
// 07/12/2004: Last Modified
//
// ###


#  include <fstream>
#  include <cstdio>
#  include "TextFileReader.hpp"
#  include "FromBinConverter.hpp"


using namespace std;


class Bin2HMeTiS
  : public FromBinConverter
{

protected:


public:

  Bin2HMeTiS();
  ~Bin2HMeTiS();
  
  void convert(const char *filename);


};


#  endif
