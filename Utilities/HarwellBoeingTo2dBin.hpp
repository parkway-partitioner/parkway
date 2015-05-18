#  ifndef _HARWELLBOEING_TO2DBIN_HPP
#  define _HARWELLBOEING_TO2DBIN_HPP


// ### HarwellBoeingTo2dBin.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
// 
// HISTORY: 
// 
// 13/4/2005: Last Modified
//
// NOTE: converts the matrix into 2d hypergraph format
//
// ###


#  include "HarwellBoeingReader.hpp"
#  include "Bit.hpp"


using namespace std;


class HarwellBoeingTo2dBin
  : public HarwellBoeingReader
{

protected:

public:
  
  HarwellBoeingTo2dBin();
  ~HarwellBoeingTo2dBin();
  
  void convert(const char *filename);
  void qsort(int *array, int left, int right);
  
  inline void swap(int &a, int &b) {
    int temp = a;
    a = b;
    b = temp;
  }
};




#  endif
