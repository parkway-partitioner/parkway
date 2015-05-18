#  ifndef _HGRAPH_TEXTREADER_HPP
#  define _HGRAPH_TEXTREADER_HPP


// ### HypergraphTextFileReader.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
// 
// HISTORY: 
// 
// 01/12/2004: Last Modified
//
// ###


#  include "StringUtils.hpp"
#  include "TextFileReader.hpp"


using namespace std;


class HypergraphTextFileReader
  : public TextFileReader
{

protected:

  int numVertices;
  int numHedges;
  int numPins;
  int wtsOnVertices;
  int wtsOnHedges;
  
public:

  HypergraphTextFileReader();
  ~HypergraphTextFileReader();

};


#  endif
