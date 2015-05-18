#  ifndef _SHUFFLE_HYPERGRAPH_CPP
#  define _SHUFFLE_HYPERGRAPH_CPP


// ### ShuffleHypergraph.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
// 
// HISTORY: 
// 
// 10/1/2005: Last Modified
//
// ###


#  include "ShuffleHypergraph.hpp"


ShuffleHypergraph::ShuffleHypergraph(int nP)
  : FromBinConverter()
{
  numParts = nP;
}

ShuffleHypergraph::~ShuffleHypergraph()
{

}



void ShuffleHypergraph::convert(const char *filename)
{
  ifstream in_stream;
  ofstream out_stream;
  
  in_stream.open(filename, ifstream::in | ifstream::binary);
  
  if(!in_stream.is_open())
    {
      cout << "error opening " << filename << endl;
      in_stream.close();
      exit(1);
    }

  readPreamble(in_stream);
  readInVertexWts(in_stream);

  
  

}




#  endif
