
#  ifndef _HYPERGRAPH_CPP
#  define _HYPERGRAPH_CPP


// ### Hypergraph.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
// 
// HISTORY: 
// 
// 5/1/2005: Last Modified
//
// ###


#  include "Hypergraph.hpp"


Hypergraph::Hypergraph(int *vWts, int numVerts)
{
  numVertices = numVerts;
  
  matchVector.setLength(numVertices);
  vWeight.setArray(vWts,numVertices);

  register int i;

  for (i=0;i<numVertices;++i)     
    matchVector[i] = -1;    
}



Hypergraph::Hypergraph(int *vWts, int *pVector, int numVerts, int cut)
{
  numVertices = numVerts;
  
  matchVector.setLength(numVertices);
  vWeight.setArray(vWts,numVertices);
  partitionVector.setArray(pVector,numVertices);

  partitionCuts.setLength(1);
  partitionVectorOffsets.setLength(2);

  partitionVectorOffsets[0] = 0;
  partitionVectorOffsets[1] = numVertices;
  partitionCuts[0] = cut;

  register int i;

  for (i=0;i<numVertices;++i)     
    matchVector[i] = -1;    
}


Hypergraph::~Hypergraph()
{
  
}



void Hypergraph::buildVtoHedges()
{
#  ifdef DEBUG_HYPERGRAPH
  assert(numVertices >= 0);
  assert(numHedges >= 0);
  assert(numPins >= 0);
  assert(pinList.getLength() >= 0);
  assert(hEdgeOffsets.getLength() >= 0);
#  endif
  
  register int i;
  register int j;
  register int ij;
  
  int endOffset;
  
  FastDynaArray<int> vDegs(numVertices);

  vToHedges.setLength(numPins);
  vOffsets.setLength(numVertices+1);
  
  for (i=0;i<numVertices;++i)    
    vDegs[i] = 0;    

  for (i=0;i<numPins;++i)
    {
#  ifdef DEBUG_HYPERGRAPH
      assert(pinList[i] >= 0 && pinList[i] < numVertices);
#  endif
      ++vDegs[pinList[i]];
    }
  
  j = 0;
  for (i=0;i<numVertices;++i)
    {
      vOffsets[i] = j;
      j += vDegs[i];
      vDegs[i] = 0;
    }
  vOffsets[i] = j;
  
  for (i=0;i<numHedges;++i)
    {
      endOffset = hEdgeOffsets[i+1];
      
      for (j=hEdgeOffsets[i];j<endOffset;++j)
	{
	  ij = pinList[j];
	  vToHedges[vOffsets[ij]+(vDegs[ij]++)] = i;
	}
    }
}



void Hypergraph::resetMatchVector()
{
  register int i;

  for (i=0;i<numVertices;++i)    
    matchVector[i] = -1;    
}


void Hypergraph::resetPartitionVectors()
{
  partitionVector.setLength(0);
  partitionCuts.setLength(0);
  partitionVectorOffsets.setLength(0);

  numPartitions = 0;
}



void Hypergraph::resetVertexMaps()
{
  resetMatchVector();
  resetPartitionVectors();
}




void Hypergraph::projectPartitions(const Hypergraph &coarseGraph)
{
  register int *coarsePartVector = coarseGraph.getPartVectorArray();

  int *coarsePartOffsets = coarseGraph.getPartOffsetArray();
  int *coarsePartCuts = coarseGraph.getPartCutArray();
  int numP = coarseGraph.getNumPartitions();

  int startOffset;
  int v;
  int i;
  
  register int j;
  register int index;

  numPartitions = numP;

  partitionCuts.setLength(numPartitions);
  partitionVectorOffsets.setLength(numPartitions+1);
  partitionVector.setLength(numPartitions*numVertices);
  
  j = 0;
  for (i=0;i<=numPartitions;++i)
    {      
      partitionVectorOffsets[i] = j;
      j += numVertices;
    }
  
  for (i=0;i<numPartitions;++i)
    {
      partitionCuts[i] = coarsePartCuts[i];
      
      index = coarsePartOffsets[i];
      startOffset = partitionVectorOffsets[i];

      for (j=0;j<numVertices;++j)
	{
	  v = matchVector[j];
	  partitionVector[startOffset+j] = coarsePartVector[index+v];
	}
    }
}


void Hypergraph::removeBadPartitions(double fractionOK)
{
  register int i;
  register int indexIntoOld;
  register int indexIntoNew;
 
  int acceptedCut;
  int bestCut = partitionCuts[0];
  int bestPartition = 0;
  int numNewPartitions = 0;
  int diffInCut; 
  int endOffset;
  int pSeenBefore;
  int j;

  for (i=1;i<numPartitions;++i)
    {
      if(partitionCuts[i] < bestCut)
	{
	  bestCut = partitionCuts[i];
	  bestPartition = i;
	}
    }

  diffInCut = static_cast<int>(floor(static_cast<double>(bestCut)*fractionOK));
  acceptedCut = bestCut + diffInCut;
  
  indexIntoNew = 0;
  indexIntoOld = 0;
  
  for (i=0;i<numPartitions;++i)
    {     
      if(partitionCuts[i] <= acceptedCut) 
	{
	  pSeenBefore = 0;

	  for (j=0;j<numNewPartitions;++j)
	    {
	      if(partitionCuts[j] == partitionCuts[i])
		pSeenBefore = 1; 
	    }

	  if(pSeenBefore == 0)
	    {	      
	      if(indexIntoOld > indexIntoNew)
		{
		  endOffset = partitionVectorOffsets[i+1];	      
		  
		  for (; indexIntoOld<endOffset;++indexIntoOld)
		    {
		      partitionVector[indexIntoNew++] = partitionVector[indexIntoOld];		   
		    }	    
		}
	      else
		{
		  indexIntoOld += numVertices;
		  indexIntoNew += numVertices;
		}

	      partitionCuts[numNewPartitions++] = partitionCuts[i];	  
	    }
	  else
	    {	      
	      indexIntoOld += numVertices;
	    }
	}
      else
	{
	  indexIntoOld += numVertices;
	}
    }
  
  numPartitions = numNewPartitions;
}




void Hypergraph::setNumPartitions(int nPartitions)
{
  register int i;
  register int j;
  
  numPartitions = nPartitions;

  partitionCuts.setLength(numPartitions);
  partitionVector.setLength(numPartitions*numVertices);
  partitionVectorOffsets.setLength(numPartitions+1);

  j = 0;

  for (i=0;i<=numPartitions;++i)
    {
      partitionVectorOffsets[i] = j;
      j += numVertices;
    }
}



void Hypergraph::copyOutPartition(register int *pVector, int nV, int pNo) const
{
#  ifdef DEBUG_HYPERGRAPH
  assert(pNo >= 0 && pNo < numPartitions);
  assert(nV == numVertices);
#  endif

  register int i;
  register int *vOffset = &partitionVector[partitionVectorOffsets[pNo]];
  
  for (i=0;i<numVertices;++i)    
    pVector[i] = vOffset[i];    
}




void Hypergraph::copyInPartition(const register int *pVector, int nV, int pNo, int cut)
{
#  ifdef DEBUG_HYPERGRAPH
  assert(pNo >= 0 && pNo < numPartitions);
  assert(nV == numVertices);
#  endif
  
  register int i;
  register int *vOffset = &partitionVector[partitionVectorOffsets[pNo]];

  for (i=0;i<numVertices;++i)          
    vOffset[i] = pVector[i];        
  
  partitionCuts[pNo] = cut;
}



void Hypergraph::printCharacteristics(ostream &o)
{
  register int i;
  register int j;
  register int ij;
  
  o << " |cGraph| " 
    << numVertices << " "
    << numHedges << " " 
    << numPins << " : ";
  
  double weighted_ave = 0;
  double percentile_75;
  double percentile_50;
  double percentile_25;
  double percentile_95;
  
  FastDynaArray<int> hEdgeLens(numHedges);
  FastDynaArray<int> hEdges(numHedges);
  FastDynaArray<int> vertices(numVertices);

  j = 0;
  for (i=0;i<numHedges;++i)
    {
      hEdges[i] = i;
      hEdgeLens[i] = hEdgeOffsets[i+1]-hEdgeOffsets[i];
      weighted_ave += (hEdgeLens[i]*hEdgeWeight[i]);
      j += hEdgeWeight[i];
    }

  percentile_75 = (static_cast<double>(j)*75)/100;
  percentile_50 = (static_cast<double>(j)*50)/100;
  percentile_95 = (static_cast<double>(j)*95)/100;
  percentile_25 = (static_cast<double>(j)*25)/100;

  o << weighted_ave/j << " ";

  Funct::qsortByAnotherArray(0,numHedges-1,hEdges.getArray(),hEdgeLens.getArray(),INC);

  j=0;
  i=0;
  ij=0;

  for ( ;i<numHedges; )
    {
      j += hEdgeWeight[hEdges[i++]];
      
      if(ij == 0 && j > percentile_25)
	{
	  o << hEdgeLens[hEdges[i]] << " ";  
	  ++ij;
	}

      if(ij == 1 && j > percentile_50)
	{
	  o << hEdgeLens[hEdges[i]] << " ";  
	  ++ij;
	}

      if(ij == 2 && j > percentile_75)
	{
      	  o << hEdgeLens[hEdges[i]] << " ";  
	  ++ij;
	}

      if(ij == 3 && j > percentile_95)
	{
      	  o << hEdgeLens[hEdges[i]] << " ";  
	  ++ij;
	}

      if(i == numHedges-1)
	{
	  o << hEdgeLens[hEdges[i]] << " : ";
	}
    }

  j = 0;
  for (i=0;i<numVertices;++i)
    {
      vertices[i] = i;      
      j += vWeight[i];
    }

  percentile_75 = (static_cast<double>(j)*75)/100;
  percentile_50 = (static_cast<double>(j)*50)/100;
  percentile_95 = (static_cast<double>(j)*95)/100;
  percentile_25 = (static_cast<double>(j)*25)/100;

  Funct::qsortByAnotherArray(0,numVertices-1,vertices.getArray(),vWeight.getArray(),INC);

  j=0;
  i=0;
  ij=0;

  for ( ;i<numVertices; )
    {
      j += vWeight[vertices[i++]];
      
      if(ij == 0 && j > percentile_25)
	{
	  o << vWeight[vertices[i]] << " ";  
	  ++ij;
	}

      if(ij == 1 && j > percentile_50)
	{
	  o << vWeight[vertices[i]] << " ";  
	  ++ij;
	}
      
      if(ij == 2 && j > percentile_75)
	{
      	  o << vWeight[vertices[i]] << " ";  
	  ++ij;
	}

      if(ij == 3 && j > percentile_95)
	{
      	  o << vWeight[vertices[i]] << " ";  
	  ++ij;
	}
      
      if(i == numVertices-1)
	{
	  o << vWeight[vertices[i]] << endl;
	}  
    }
}



int Hypergraph::keepBestPartition()
{
  register int i;
  register int bestOffset;
  
  int bestPartition = 0;
  int bestCut = partitionCuts[0];

  for (i=1;i<numPartitions;++i)
    {
      if(partitionCuts[i] < bestCut)
	{
	  bestPartition = i;
	  bestCut = partitionCuts[i];
	}
    }
  
  if(bestPartition != 0)
    {
      bestOffset = partitionVectorOffsets[bestPartition];

      for (i=0;i<numVertices;++i)
	{
	  partitionVector[i] = partitionVector[bestOffset+i];
	}
    }

  partitionCuts[0] = bestCut;
  numPartitions = 1;

  return bestCut;
}



int Hypergraph::getExposedHedgeWt() const
{
  register int i;
  register int ij=0;

  for (i=0;i<numHedges;++i)    
    ij += hEdgeWeight[i];    
  
  return ij;
}


int Hypergraph::calcCutsize(int nP, int partitionNo) const
{
  register int i;
  register int j;
  register int *pVector = &partitionVector[partitionVectorOffsets[partitionNo]];

  FastDynaArray<int> spanned(nP);

  int k_1Cut = 0;
  int endOffset;
  int vPart;
  int numSpanned;

  for (i=0;i<numHedges;++i)
    {
      endOffset = hEdgeOffsets[i+1];
      numSpanned = 0;

      for (j=0;j<nP;++j)
	spanned[j] = 0;

      for (j=hEdgeOffsets[i];j<endOffset;++j)
	{
	  vPart = pVector[pinList[j]];
#  ifdef DEBUG_HYPERGRAPH        
	  assert(vPart >= 0 && vPart < nP);
#  endif
	  if(!spanned[vPart]) 
	    {
	      spanned[vPart] = 1;
	      ++numSpanned;
	    }
	}

      k_1Cut += ((numSpanned-1)*hEdgeWeight[i]);
    }

  return k_1Cut;
}


int Hypergraph::getSOED(int nP, int partitionNo) const
{
  register int i;
  register int j;
  register int *pVector = &partitionVector[partitionVectorOffsets[partitionNo]];

  FastDynaArray<int> spanned(nP);

  int soed = 0;
  int endOffset;
  int vPart;
  int numSpanned;

  for (i=0;i<numHedges;++i)
    {
      endOffset = hEdgeOffsets[i+1];
      numSpanned = 0;

      for (j=0;j<nP;++j)
	spanned[j] = 0;

      for (j=hEdgeOffsets[i];j<endOffset;++j)
	{
	  vPart = pVector[pinList[j]];
#  ifdef DEBUG_HYPERGRAPH        
	  assert(vPart >= 0 && vPart < nP);
#  endif
	  if(!spanned[vPart]) 
	    {
	      spanned[vPart] = 1;
	      ++numSpanned;
	    }
	}

      if(numSpanned > 1)
	soed += (numSpanned*hEdgeWeight[i]);
    }
  
  return soed;
}



void Hypergraph::initCutsizes(int numParts)
{
#  ifdef DEBUG_HYPERGRAPH
  assert(numPartitions >= 1);
  assert(partitionCuts.getLength() >= 1);
  assert(partitionVector.getLength() >= numVertices);
  assert(partitionVectorOffsets.getLength() > 1);
#  endif

  register int i;

  for (i=0;i<numPartitions;++i)    
    partitionCuts[i] = calcCutsize(numParts,i);
}




void Hypergraph::checkPartitions(int nP, int maxWt) const
{
  register int i;
  register int j;
  register int *pVector;

  int cut;

  FastDynaArray<int> partWts(nP);

  for (i=0;i<numPartitions;++i)
    {
      cut = calcCutsize(nP,i);
      
      assert(partitionCuts[i] == cut);
      
      for (j=0;j<nP;++j)	
	partWts[j] = 0;      

      pVector = &partitionVector[partitionVectorOffsets[i]];

      for (j=0;j<numVertices;++j)	
	partWts[pVector[j]] += vWeight[j];    

      for (j=0;j<nP;++j)	
	assert(partWts[j] <= maxWt);	  
    }
}


void Hypergraph::checkPartition(int numPartition, int nP, int maxWt) const
{
  register int i;
  register int *pVector;

  int cut;

  FastDynaArray<int> partWts(nP);

  cut = calcCutsize(nP,numPartition);
  assert(cut == partitionCuts[numPartition]); 

  for (i=0;i<nP;++i)    
    partWts[i] = 0;    
  
  pVector = &partitionVector[partitionVectorOffsets[numPartition]];
  
  for (i=0;i<numVertices;++i)    
    partWts[pVector[i]] += vWeight[i];    
 
  for (i=0;i<nP;++i)	
    assert(partWts[i] <= maxWt);	
}


void Hypergraph::convertToDIMACSGraphFile(const char *fN) const
{
  int i;
  int j;
  int ij;
  int endOffset;
  int startOffset;
  int v1;
  int v2;
  int numEdges = 0;
  int maxVdegree = 0;

  char writeString[2064];

  ofstream out_stream;

  FastDynaArray<int> numVNeighs(numVertices);
  FastDynaArray<int> vDegs(numVertices);

  FastDynaArray<FastDynaArray<int>*> vNeighs(numVertices);

  for (i=0;i<numVertices;++i)
    {
      numVNeighs[i] = 0;
      vNeighs[i] = new FastDynaArray<int>(32);
    }
 
  // ###
  // first qsort the hyperedges 
  // then create the adjacency list
  // ###

  for (i=0;i<numHedges;++i)
    {
      startOffset = hEdgeOffsets[i];
      endOffset = hEdgeOffsets[i+1]-1;
      
      Funct::qsort(startOffset,endOffset,pinList.getArray());

      for (j=startOffset;j<endOffset;++j)
	{
	  v1 = pinList[j];

	  for (ij=j+1;ij<endOffset;++ij)
	    {
	      v2 = pinList[ij];

	      if(Funct::search(vNeighs[v1]->getArray(),numVNeighs[v1],v2)==-1)
		{
		  vNeighs[v1]->assign(numVNeighs[v1]++,v2);
		  ++numEdges;
		}
	    }
	}
    }

  for (i=0;i<numVertices;++i)
    vDegs[i] = 0;

  for (i=0;i<numVertices;++i)
    {
      for (j=0;j<numVNeighs[i];++j)
	{	   
	  v1 = (*vNeighs[i])[j];
	  ++vDegs[i];
	  ++vDegs[v1];
	}
    }

  
  for (i=0;i<numVertices;++i)
    if(vDegs[i] > maxVdegree)
      maxVdegree = vDegs[i];
  
  // ###
  // now translate the adjacency list to file
  // ###
  
  out_stream.open(fN, ofstream::out | ofstream::binary);

  if(!out_stream.is_open())
    {
      cout << "error opening " << fN << endl;     	  
      return;
    }

  strcpy(writeString, "c Hypergraph to Graph representation\n");
  out_stream.write(&writeString[0], strlen(writeString));

  sprintf(writeString, "p edge %d %d\n",numVertices,numEdges);
  out_stream.write(&writeString[0], strlen(writeString));

  for (i=0;i<numVertices;++i)
    {
      endOffset = numVNeighs[i];

      for (j=0;j<endOffset;++j)
	{
	  v2 = (*vNeighs[i])[j];
	  
	  sprintf(writeString, "e %d %d\n",i+1,v2+1);
	  out_stream.write(&writeString[0], strlen(writeString));
	}
    }
  
  out_stream.close();

  for (i=0;i<numVertices;++i)
    DynaMem<FastDynaArray<int> >::deletePtr(vNeighs[i]);
}


void Hypergraph::printPercentiles(ostream &o)
{
  register int i;
  register int j;
  register int ij;
  
  o   << " #vertices: " << numVertices 
      << " #hyperedges: " << numHedges 
      << " #pins: " << numPins << endl;

  double weighted_ave = 0;
  double percentile_95;
  double percentile_75;
  double percentile_50;
  double percentile_25;
  
  FastDynaArray<int> indices;
  FastDynaArray<int> hEdgeLens(numHedges);

  /* display hyperedge information */

  j = 0;
  indices.setLength(numHedges);

  for (i=0;i<numHedges;++i)
    {
      indices[i] = i;     
      j += hEdgeWeight[i];
      hEdgeLens[i] = hEdgeOffsets[i+1]-hEdgeOffsets[i];
      weighted_ave += (hEdgeLens[i]*hEdgeWeight[i]);
    }

  percentile_95 = (static_cast<double>(j)*95)/100;
  percentile_75 = (static_cast<double>(j)*75)/100;
  percentile_50 = (static_cast<double>(j)*50)/100;
  percentile_25 = (static_cast<double>(j)*25)/100;

  o << "hyperedge length percentiles: (weighted ave, 25, 50, 75, 95, maxLength) " << endl;
  o << "\t" << weighted_ave/j << " ";

  Funct::qsortByAnotherArray(0,numHedges-1,indices.getArray(),hEdgeLens.getArray(),INC);

  j=0;
  i=0;
  ij=0;

  for ( ;i<numHedges; )
    {
      j += hEdgeWeight[indices[i++]];
      
      if(ij == 0 && j > percentile_25)
	{ 
	  o << hEdgeLens[indices[i]] << " ";  
	  ++ij;
	}

      if(ij == 1 && j > percentile_50)
	{
	  o << hEdgeLens[indices[i]] << " ";  
	  ++ij;
	}

      if(ij == 2 && j > percentile_75)
	{
      	  o << hEdgeLens[indices[i]] << " ";  
	  ++ij;
	}

      if(ij == 3 && j > percentile_95)
	{
      	  o << hEdgeLens[indices[i]] << " ";  
	  ++ij;
	}

      if(i == numHedges-1)
	{
	  o << hEdgeLens[indices[i]] << endl;
	}
    }

  /* display vertex information */

  j = 0;
  indices.setLength(numVertices);

  for (i=0;i<numVertices;++i)
    {
      indices[i] = i;      
      j += vWeight[i];
    }

  percentile_95 = (static_cast<double>(j)*95)/100;
  percentile_75 = (static_cast<double>(j)*75)/100;
  percentile_50 = (static_cast<double>(j)*50)/100;
  percentile_25 = (static_cast<double>(j)*25)/100;

  Funct::qsortByAnotherArray(0,numVertices-1,indices.getArray(),vWeight.getArray(),INC);  

  o << "vertex weight percentiles: (ave, 25, 50, 75, 95, maxWeight) " << endl;
  o << "\t" << static_cast<double>(j)/numVertices << " ";

  j=0;
  i=0;
  ij=0;

  for ( ;i<numVertices; )
    {
      j += vWeight[indices[i++]];
      
      if(ij == 0 && j > percentile_25)
	{
	  o << vWeight[indices[i]] << " ";  
	  ++ij;
	}

      if(ij == 1 && j > percentile_50)
	{
	  o << vWeight[indices[i]] << " ";  
	  ++ij;
	}
      
      if(ij == 2 && j > percentile_75)
	{
      	  o << vWeight[indices[i]] << " ";  
	  ++ij;
	}

      if(ij == 3 && j > percentile_95)
	{
      	  o << vWeight[indices[i]] << " ";  
	  ++ij;
	}
      
      if(i == numVertices-1)
	{
	  o << vWeight[indices[i]] << endl;
	}  
    }
}



#  endif
