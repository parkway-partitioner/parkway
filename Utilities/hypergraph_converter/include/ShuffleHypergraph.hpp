#  ifndef _SHUFFLE_HYPERGRAPH_HPP
#  define _SHUFFLE_HYPERGRAPH_HPP


// ### ShuffleHypergraph.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
// 
// HISTORY: 
// 
// 10/1/2005: Last Modified
//
// ###


#  include "FromBinConverter.hpp"


using namespace std;


class ShuffleHypergraph
  : public FromBinConverter
{
protected:

  int numParts;
 
public:
  
  ShuffleHypergraph(int nP);
  ~ShuffleHypergraph();

  void convert(const char *filename);
  
  
};


#  endif

