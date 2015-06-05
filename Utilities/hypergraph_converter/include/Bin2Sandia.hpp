#  ifndef _BIN2SANDIA_HPP
#  define _BIN2SANDIA_HPP


// ### Bin2Sandia.hpp ###
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


class Bin2Sandia
  : public FromBinConverter
{

protected:


public:

  Bin2Sandia();
  ~Bin2Sandia();
  
  void convert(const char *filename);


};


#  endif
