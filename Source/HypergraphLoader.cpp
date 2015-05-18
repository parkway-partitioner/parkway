#  ifndef _LOADER_CPP
#  define _LOADER_CPP


// ### HypergraphLoader.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
// 
// HISTORY: 
// 
// 4/1/2005: Last Modified
//
// ###


#  include "HypergraphLoader.hpp"


HypergraphLoader::HypergraphLoader(int disp)
{
  dispOption = disp;
  currPercentile = 100;

  numHedges = 0;
  numVertices = 0;
  numPins = 0;
  numPartitions = 0;
  
  vWeight = NULL;
  hEdgeWeight = NULL;
  matchVector = NULL;
  pinList = NULL;
  hEdgeOffsets = NULL;
  vToHedges = NULL;
  vOffsets = NULL;
  partitionVectors = NULL;
  partitionOffsets = NULL;
  partitionCutsizes = NULL;
}
  

HypergraphLoader::~HypergraphLoader()
{

}


void HypergraphLoader::computeHedgesToLoad(BitField &toLoad)
{
  register int i;
  register int j = 0;

  int percentileLen;
  double percentileThreshold;

  FastDynaArray<int> hEdgeLens(numHedges);
  FastDynaArray<int> hEdges(numHedges);

  for (i=0;i<numHedges;++i)
    {
      hEdges[i] = i;
      hEdgeLens[i] = hEdgeOffsets[i+1]-hEdgeOffsets[i];
      j += hEdgeWeight[i];
    }

  percentileThreshold = (static_cast<double>(j)*currPercentile)/100;
  Funct::qsortByAnotherArray(0,numHedges-1,hEdges.getArray(),hEdgeLens.getArray(),INC);
  
  j=0;
  i=0;

  for ( ;i<numHedges && j<percentileThreshold; )
    j += hEdgeWeight[hEdges[i++]];
  
  percentileLen = hEdgeLens[hEdges[i]];

  for ( ;i<numHedges;++i)
    if(hEdgeLens[hEdges[i]] > percentileLen)
      toLoad.set0(hEdges[i]);
}




#  endif
