
#  ifndef _PARA_REFINER_CPP
#  define _PARA_REFINER_CPP


// ### ParaRefiner.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
// 
// HISTORY: 
// 
// 3/2/2005: Last Modified
//
// ###


#  include "ParaRefiner.hpp"


ParaRefiner::ParaRefiner(int rank, int nProcs, int nParts, ostream &out)
  : ParaHypergraphLoader(rank,nProcs,nParts,out)
{
  partWeights.setLength(nParts);
}


ParaRefiner::~ParaRefiner()
{

}



void ParaRefiner::loadHyperGraph(const ParaHypergraph &h, MPI_Comm comm)
{
  int i;
  int ij;
  int vertsPerProc;
  int endOffset;
  int startOffset;
  int hEdgeLen;
  int recvLen;
  int nonLocIndex;

  int numLocalPins;
  int numLocalHedges;

  int *localPins;
  int *localHedgeOffsets;
  int *localHedgeWeights;
  
  register int j;
  register int l;
  register int proc;
  register int locVert;

  FastDynaArray<int> sentToProc;
  FastDynaArray<int> vDegs;

  numLocalPins = h.getNumLocalPins();
  numLocalHedges = h.getNumLocalHedges();
  localPins = h.getLocalPinsArray();
  localHedgeOffsets = h.getHedgeOffsetsArray();
  localHedgeWeights = h.getHedgeWeightsArray();

  locVertWt = h.getLocalVertexWt();
  vWeight = h.getWeightArray();
  matchVector = h.getMatchVectorArray();

  numLocalVertices = h.getNumLocalVertices();
  totalVertices = h.getNumTotalVertices();
  minVertexIndex = h.getMinVertexIndex();
  maxVertexIndex = minVertexIndex+numLocalVertices;
  
  // ###
  // Prepare data structures
  // ###

  numAllocHedges = 0;
  numHedges = 0;
  numLocPins = 0;
  vertsPerProc = totalVertices / numProcs;
  
  vToHedgesOffset.setLength(numLocalVertices+1);
  sentToProc.setLength(numProcs);
  vDegs.setLength(numLocalVertices);

  for (i=0;i<numLocalVertices;++i) 
    {    
      vToHedgesOffset[i] = 0;
      vDegs[i] = 0;
    }

  if(dispOption > 1 && myRank == 0)
    out_stream << "[PR]";

  // ###
  // use the request sets to send local hyperedges to other
  // processors and to receive hyperedges from processors
  // ###

  for (i=0;i<numProcs;++i) 
    {
      sendLens[i] = 0;
      sentToProc[i] = 0;
    }

  if(currPercentile == 100)
    {
      int numActiveProcs;
      int activeProc;
      
      FastDynaArray<int> activeProcs(numProcs);

      for (i=0;i<numLocalHedges;++i)
	{      
	  startOffset = localHedgeOffsets[i];
	  endOffset = localHedgeOffsets[i+1];
	  hEdgeLen = endOffset - startOffset;
	  numActiveProcs = 0;
	  
	  for (j=startOffset;j<endOffset;++j) 
	    {	
#  ifdef DEBUG_REFINER
	      assert(localPins[j] < totalVertices && localPins[j] >= 0);
#  endif
	      proc = min(localPins[j]/vertsPerProc, numProcs-1);

	      if(!sentToProc[proc]) 
		{	  
		  if(proc == myRank) 
		    {	    
		      hEdgeWeight.assign(numHedges, localHedgeWeights[i]);
		      hEdgeOffset.assign(numHedges++, numLocPins);
		      
		      for (l=startOffset;l<endOffset;++l) 
			{
			  locPinList.assign(numLocPins++, localPins[l]);	
#  ifdef DEBUG_REFINER    
			  assert(localPins[l]<totalVertices && localPins[l]>=0);
#  endif
			  if(localPins[l]>=minVertexIndex && localPins[l]<maxVertexIndex) 
			    ++vToHedgesOffset[localPins[l]-minVertexIndex];
			}	    
		    }
		  else 
		    {		    
		      dataOutSets[proc]->assign(sendLens[proc]++, hEdgeLen+2);
		      dataOutSets[proc]->assign(sendLens[proc]++, localHedgeWeights[i]);	  
		      
		      for (l=startOffset;l<endOffset;++l)
			{
			  dataOutSets[proc]->assign(sendLens[proc]++, localPins[l]);
			}		
		    }
		  
		  sentToProc[proc] = 1;
		  activeProcs[numActiveProcs++] = proc;
		}  
	    }
	  
	  activeProc = activeProcs[RANDOM(0,numActiveProcs)];
	  
#  ifdef DEBUG_REFINER
	  assert(activeProc >= 0 && activeProc < numProcs);
#  endif
	  
	  if(activeProc == myRank)
	    {
	      allocHedges.assign(numAllocHedges++, numHedges-1);
	    }
	  else
	    {
	      dataOutSets[activeProc]->assign(sendLens[activeProc]++, RESP_FOR_HEDGE);
	    }
	  
	  for (j=0;j<numProcs;++j)
	    {
	      sentToProc[j] = 0;
	    }
	}
    }
  else
    if(currPercentile > 0)
      {
#  ifdef DEBUG_REFINER
	assert(currPercentile > 0 && currPercentile < 100);
#  endif
	
	/* When not loading all the hyperedges, do not need to
	   make processors responsible for hyperedges during
	   the cutsize calculations                            */
	
	BitField toLoad(numLocalHedges);
	
	computeHedgesToLoad(toLoad,numLocalHedges,numLocalPins,localHedgeWeights,localHedgeOffsets, comm);
	
	for (i=0;i<numLocalHedges;++i)
	  {
	    if(toLoad(i) == 1)
	      {
		startOffset = localHedgeOffsets[i];
		endOffset = localHedgeOffsets[i+1];
		hEdgeLen = endOffset - startOffset;
		
		for (j=startOffset;j<endOffset;++j) 
		  {	
#  ifdef DEBUG_REFINER
		    assert(localPins[j] < totalVertices && localPins[j] >= 0);
#  endif
		    proc = min(localPins[j]/vertsPerProc, numProcs-1);
		    
		    if(!sentToProc[proc]) 
		      {	  
			if(proc == myRank) 
			  {	    
			    hEdgeWeight.assign(numHedges, localHedgeWeights[i]);
			    hEdgeOffset.assign(numHedges++, numLocPins);
			    
			    for (l=startOffset;l<endOffset;++l) 
			      {
				locPinList.assign(numLocPins++, localPins[l]);	
#  ifdef DEBUG_REFINER    
				assert(localPins[l]<totalVertices && localPins[l]>=0);
#  endif
				if(localPins[l]>=minVertexIndex && localPins[l]<maxVertexIndex) 
				  ++vToHedgesOffset[localPins[l]-minVertexIndex];
			      }	    
			  }
			else 
			  {		    
			    dataOutSets[proc]->assign(sendLens[proc]++, hEdgeLen+2);
			    dataOutSets[proc]->assign(sendLens[proc]++, localHedgeWeights[i]);	  
			    
			    for (l=startOffset;l<endOffset;++l)
			      {
				dataOutSets[proc]->assign(sendLens[proc]++, localPins[l]);
			      }		
			  }
			
			sentToProc[proc] = 1;
		      }  
		  }	      
		
		for (j=0;j<numProcs;++j)		
		  sentToProc[j] = 0;		
	      }
	  }
      }
    else
      {
	/* compute a fixed limit on hyperedges to be communicated 
	   - first try maxLen/2
	*/

	int maxHedgeLen = 0;
	int maxTotHedgeLen;
	int limit;

	for (i=0;i<numLocalHedges;++i)
	  if(localHedgeOffsets[i+1]-localHedgeOffsets[i]>maxHedgeLen)
	    maxHedgeLen = localHedgeOffsets[i+1]-localHedgeOffsets[i];

	MPI_Allreduce(&maxHedgeLen, &maxTotHedgeLen, 1, MPI_INT, MPI_MAX, comm);
	
	limit = maxTotHedgeLen / 2;

	for (i=0;i<numLocalHedges;++i)
	  {
	    startOffset = localHedgeOffsets[i];
	    endOffset = localHedgeOffsets[i+1];
	    hEdgeLen = endOffset - startOffset;
		
	    if(hEdgeLen < limit)
	      {
		for (j=startOffset;j<endOffset;++j) 
		  {	
#  ifdef DEBUG_REFINER
		    assert(localPins[j] < totalVertices && localPins[j] >= 0);
#  endif
		    proc = min(localPins[j]/vertsPerProc, numProcs-1);
		    
		    if(!sentToProc[proc]) 
		      {	  
			if(proc == myRank) 
			  {	    
			    hEdgeWeight.assign(numHedges, localHedgeWeights[i]);
			    hEdgeOffset.assign(numHedges++, numLocPins);
			    
			    for (l=startOffset;l<endOffset;++l) 
			      {
				locPinList.assign(numLocPins++, localPins[l]);	
#  ifdef DEBUG_REFINER    
				assert(localPins[l]<totalVertices && localPins[l]>=0);
#  endif
				if(localPins[l]>=minVertexIndex && localPins[l]<maxVertexIndex) 
				  ++vToHedgesOffset[localPins[l]-minVertexIndex];
			      }	    
			  }
			else 
			  {		    
			    dataOutSets[proc]->assign(sendLens[proc]++, hEdgeLen+2);
			    dataOutSets[proc]->assign(sendLens[proc]++, localHedgeWeights[i]);	  
			    
			    for (l=startOffset;l<endOffset;++l)
			      {
				dataOutSets[proc]->assign(sendLens[proc]++, localPins[l]);
			      }		
			  }
			
			sentToProc[proc] = 1;
		      }  
		  }	      
		
		for (j=0;j<numProcs;++j)		
		  sentToProc[j] = 0;		
	      }
	  }
      }

  // ###
  // Now exchange the hyperedges
  // ###
  
  sendFromDataOutArrays(comm);

  // ###
  // Now load the non-local hyperedges
  // ###

  j = 0;
  for (i=0;i<numProcs;++i)
    { 
      j += recvLens[i];
    }

  recvLen = j;
  j = 0;

  while(j < recvLen) 
    {
#  ifdef DEBUG_REFINER
      assert(receiveArray[j] > 0);
#  endif
      endOffset = j+receiveArray[j];
      ++j;            
      
      hEdgeWeight.assign(numHedges, receiveArray[j++]);
      hEdgeOffset.assign(numHedges++, numLocPins);
      
      for ( ;j<endOffset;++j) 
	{      
	  locPinList.assign(numLocPins++, receiveArray[j]);	
	  
#  ifdef DEBUG_REFINER
	  assert(receiveArray[j] < totalVertices && receiveArray[j] >= 0);	  
#  endif
	  
	  locVert = receiveArray[j] - minVertexIndex;
	  
	  if(locVert >= 0 && locVert < numLocalVertices) 
	    ++vToHedgesOffset[locVert];
	}
#  ifdef DEBUG_REFINER
      assert(j == endOffset);
#  endif
      
      if(currPercentile == 100 && j < recvLen && receiveArray[j] == RESP_FOR_HEDGE)
	{
	  allocHedges.assign(numAllocHedges++, numHedges-1);
	  ++j;
	}
    }
  
  hEdgeOffset.assign(numHedges, numLocPins);

#  ifdef MEM_OPT
  hEdgeOffset.setLength(numHedges+1);
  hEdgeWeight.setLength(numHedges);
  allocHedges.setLength(numHedges);
  locPinList.setLength(numLocPins);
#  endif
  
  // ###
  // now initialise the vToHedgesList
  // ###

  j=0;
  l=0;
 
  for ( ;j<numLocalVertices;++j) 
    {
#  ifdef DEBUG_REFINER
      assert(vToHedgesOffset[j] >= 0);
#  endif

      locVert = vToHedgesOffset[j];
      vToHedgesOffset[j] = l;
      l += locVert;
    }
  
  vToHedgesOffset[j] = l;
  vToHedgesList.setLength(l);
  
  for (j=0;j<numHedges;++j) 
    {    
      endOffset = hEdgeOffset[j+1];
    
      for (l=hEdgeOffset[j];l<endOffset;++l) 
	{      
	  if(locPinList[l] >= minVertexIndex && locPinList[l] < maxVertexIndex) {
	
	    locVert = locPinList[l]-minVertexIndex;
	    vToHedgesList[vToHedgesOffset[locVert]+(vDegs[locVert]++)] = j;
	  }
	}
    }
  
#  ifdef DEBUG_REFINER
  for (j=1;j<=numLocalVertices;++j)    
    assert(vDegs[j-1] == vToHedgesOffset[j]-vToHedgesOffset[j-1]);    
#  endif

  /* now init the non-local vertex structs */

  numNonLocVerts = 0;
  
  if(numLocalPins < totalVertices/2)
    toNonLocVerts.createTable(numLocalPins,1);
  else
    toNonLocVerts.createTable(totalVertices,0);
  
  for (i=0;i<numHedges;++i) 
    {      
      endOffset = hEdgeOffset[i+1];
      
      for (j=hEdgeOffset[i]; j<endOffset;++j) 
	{
	  ij = locPinList[j];
#  ifdef DEBUG_REFINER
	  assert(ij >= 0 && ij < totalVertices);
#  endif	  
	  if(ij < minVertexIndex || ij >= maxVertexIndex) 
	    {	
	      nonLocIndex = toNonLocVerts.insertIfEmpty(ij,numNonLocVerts);

	      if(nonLocIndex == -1)
		{
		  vDegs.assign(numNonLocVerts, 1);
		  nonLocVerts.assign(numNonLocVerts++, ij);		  
		}
	      else
		{
#  ifdef DEBUG_REFINER
		  assert(nonLocIndex < numNonLocVerts);
#  endif
		  ++vDegs[nonLocIndex];
		}
	    }
	}
    }

  numNonLocVertsHedges = 0;
  for (i=0;i<numNonLocVerts;++i)
    {
#  ifdef DEBUG_REFINER
      assert(vDegs[i] >= 1);
#  endif
      numNonLocVertsHedges += vDegs[i];
    }

  nonLocVerts.setLength(numNonLocVerts);
  nonLocOffsets.setLength(numNonLocVerts+1);
  nonLocVToHedges.setLength(numNonLocVertsHedges);

#  ifdef DEBUG_REFINER
  for (i=0;i<numNonLocVertsHedges;++i)
    nonLocVToHedges[i] = -1;
#  endif

  ij = 0;
  for (i=0;i<numNonLocVerts;++i)
    {
      nonLocOffsets[i] = ij;      
      ij += vDegs[i];
      vDegs[i] = 0;
    }
  nonLocOffsets[i] = ij; 
  
  /* now intialise the vToHedges for non-local vertices */
  
  for (i=0;i<numHedges;++i) 
    {      
      endOffset = hEdgeOffset[i+1];
      
      for (j=hEdgeOffset[i]; j<endOffset; ++j) 
	{
	  ij = locPinList[j];
#  ifdef DEBUG_REFINER
	  assert(ij >= 0 && ij < totalVertices);
#  endif	  
	  if(ij < minVertexIndex || ij >= maxVertexIndex) 
	    {	
	      nonLocIndex = toNonLocVerts.getVal(ij);
#  ifdef DEBUG_REFINER
	      assert(nonLocIndex >= 0 && nonLocIndex < numNonLocVerts);
#  endif
	      l = nonLocOffsets[nonLocIndex]+vDegs[nonLocIndex];
	      nonLocVToHedges[l] = i;
	      ++vDegs[nonLocIndex];	     
	    }
	}
    }
 
#  ifdef DEBUG_REFINER
  for (i=0;i<numNonLocVertsHedges;++i)
    assert(nonLocVToHedges[i] > -1);
#  endif
  
  if(dispOption > 1)
    {
      int numTotHedgesInGraph;
      int numTotPinsInGraph;
      
      MPI_Reduce(&numLocalHedges, &numTotHedgesInGraph, 1, MPI_INT, MPI_SUM, 0, comm);
      MPI_Reduce(&numLocalPins, &numTotPinsInGraph, 1, MPI_INT, MPI_SUM, 0, comm);
      
      if(myRank == 0)
	{
	  out_stream << " " << totalVertices 
		     << " " << numTotHedgesInGraph 
		     << " " << numTotPinsInGraph 
		     << endl;
	}      
    }
}



void ParaRefiner::initPartitionStructs(const ParaHypergraph &h, MPI_Comm comm)
{
  loadHyperGraph(h,comm);
    
  // ###
  // init max part weight
  // ###
  
  int totWt;
  int vPerProc;
  int arraySize;
  int vertex;
  int totalToRecv;
  int endOffset;
  int totalToSend;
  
  int *array;
  
  register int i;
  register int j;
  register int ij;

  FastDynaArray<int> copyOfSendArray;

  MPI_Allreduce(&locVertWt, &totWt, 1, MPI_INT, MPI_SUM, comm);
  
  avePartWt = static_cast<double>(totWt) / numParts;
  maxPartWt = static_cast<int>(floor(avePartWt+avePartWt*balConstraint));

  numPartitions = h.getNumPartitions();
  partitionVector = h.getPartitionArray();
  partitionVectorOffsets = h.getPartitionOffsetsArray();
  partitionCuts = h.getCutsizesArray();

  vPerProc = totalVertices / numProcs;

#  ifdef DEBUG_REFINER 
  for (register int i=0;i<partitionVectorOffsets[numPartitions];++i)
    assert(partitionVector[i] >= 0 && partitionVector[i] < numParts);
#  endif

  partIndices.setLength(numNonLocVerts*numPartitions);
  indexIntoPartIndices.setLength(numPartitions+1);

  j = 0;
  for (i=0;i<numPartitions;++i)
    {
      indexIntoPartIndices[i] = j;
      j += numNonLocVerts;
    }
  
  /*
    now communicate the partition vector
    requests for values of non-local vertices
  */
  
  for (i=0;i<numProcs;++i)    
    sendLens[i] = 0;   

  for (i=0;i<numNonLocVerts;++i)
    {
      j = nonLocVerts[i];
#  ifdef DEBUG_REFINER
      assert(j < minVertexIndex || j >= maxVertexIndex);
#  endif
      ij = min(j/vPerProc,numProcs-1);
#  ifdef DEBUG_REFINER
      assert(ij != myRank);
#  endif            
      dataOutSets[ij]->assign(sendLens[ij]++,j);
    }
  
  ij = 0;
  for (i=0;i<numProcs;++i)
    {
      sendDispls[i] = ij;
      ij += sendLens[i];
    }

#  ifdef DEBUG_REFINER
  assert(ij == numNonLocVerts);
  assert(sendLens[myRank] == 0);
#  endif
  
  sendArray.setLength(ij);
  copyOfSendArray.setLength(ij);
  arraySize = ij;
  
  ij = 0;
  for (i=0;i<numProcs;++i)
    {
      endOffset = sendLens[i];
      array = dataOutSets[i]->getArray();
      
      for (j=0;j<endOffset;++j)
	{
#  ifdef DEBUG_REFINER
	  assert(array[j] < minVertexIndex || array[j] >= maxVertexIndex);
#  endif
	  sendArray[ij] = array[j];
	  copyOfSendArray[ij++] = array[j];
	}
    }

  MPI_Alltoall(sendLens.getArray(), 1, MPI_INT, recvLens.getArray(), 1, MPI_INT, comm);
  
  ij = 0;
  for (i=0;i<numProcs;++i) 
    {    
      recvDispls[i] = ij;
      ij += recvLens[i];
    }

  totalToRecv = ij;
  receiveArray.setLength(ij);
  
  MPI_Alltoallv(sendArray.getArray(), sendLens.getArray(), sendDispls.getArray(), MPI_INT, receiveArray.getArray(), recvLens.getArray(), recvDispls.getArray(), MPI_INT, comm);
  
// ###
  // now communicate the partition vector
  // requests for values of non-local vertices
  // ###  
  
  totalToSend = totalToRecv*numPartitions;
  sendArray.setLength(totalToSend);
  
  for (i=0;i<numProcs;++i)
    {
      ij = sendLens[i];
      sendLens[i] = recvLens[i]*numPartitions;
      recvLens[i] = ij*numPartitions;      
    }
  
  ij = 0;
  j = 0;

  for (i=0;i<numProcs;++i)
    {
      sendDispls[i] = ij;
      recvDispls[i] = j;

      ij += sendLens[i];
      j += recvLens[i];
    }

  // ###
  // now get the partition vector values
  // of the requested local vertices
  // ###  

#  ifdef DEBUG_REFINER
  assert(receiveArray.getLength() == totalToRecv);
  assert(sendArray.getLength() == totalToRecv*numPartitions);
#  endif

  ij = 0;
  for (i=0;i<totalToRecv;++i)
    {
      vertex = receiveArray[i]-minVertexIndex;
      
#  ifdef DEBUG_REFINER
      assert(vertex >= 0 && vertex < maxVertexIndex-minVertexIndex);
#  endif
      
      for (j=0;j<numPartitions;++j)
	{
#  ifdef DEBUG_REFINER
	  int part = partitionVector[partitionVectorOffsets[j]+vertex];
	  assert(part >= 0 && part < numParts);
	  assert(ij < totalToRecv*numPartitions);
#  endif	  
	  sendArray[ij++] = partitionVector[partitionVectorOffsets[j]+vertex];
	}
    }
#  ifdef DEBUG_REFINER
  assert(ij == totalToSend);
#  endif

  totalToRecv = numNonLocVerts*numPartitions;
  receiveArray.setLength(totalToRecv);
  
  MPI_Alltoallv(sendArray.getArray(), sendLens.getArray(), sendDispls.getArray(), MPI_INT, receiveArray.getArray(), recvLens.getArray(), recvDispls.getArray(), MPI_INT, comm);

  ij = 0;
  for (i=0;i<arraySize;++i)
    {     
      vertex = toNonLocVerts.getVal(copyOfSendArray[i]);
#  ifdef DEBUG_REFINER
      assert(vertex >= 0 && vertex < numNonLocVerts);
#  endif            
      for (j=0;j<numPartitions;++j)		  
	partIndices[indexIntoPartIndices[j]+vertex]=receiveArray[ij++];
    }
}








#   endif
