#  ifndef _BIN2PATOH_HPP
#  define _BIN2PATOH_HPP


// ### Bin2PaToH.hpp ###
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


class Bin2PaToH
  : public FromBinConverter
{

protected:


public:

  Bin2PaToH();
  ~Bin2PaToH();

  void convert(const char *filename);

}; 



#  endif
