
#  ifndef _WEB_GRAPH_PARTITIONER_CPP
#  define _WEB_GRAPH_PARTITIONER_CPP


// ### WebGraphPartitioner.cpp ###
//
// Copyright (C) 2005, Aleksandar Trifunovic, Imperial College London
// 
// HISTORY: 
// 
// 21/3/2005: Last Modified
//
// ###


#  include "WebGraphPartitioner.hpp"



WebGraphPartitioner::~WebGraphPartitioner()
{


}


//// Need to think about this - take away the bisseling paper and the aykanat paper...


void WebGraphPartitioner::runPartitioner(MPI_Comm comm)
{
  int hEdgePercentile;
  int firstCutSize;
  int secondCutSize;
  int i;
  int diffInCutSize;
  int numIterations;
  int vCycleGain;
  int minIterationGain;

  int percentCoarsening;
  int percentSequential;
  int percentRefinement;
  int percentOther;
  
  double totStartTime;

  Stack<int> hEdgePercentiles;

  ParaHypergraph *coarseGraph;
  ParaHypergraph *finerGraph;

  initMapToOrigVerts();

  bestCutsize = LARGE_CONSTANT;
  worstCutsize = 0;
  totCutsizes = 0;

  totCoaTime = 0;
  totSeqTime = 0;
  totRefTime = 0;

  MPI_Barrier(comm);
  totStartTime = MPI_Wtime();

#  ifdef DEBUG_CONTROLLER
  int checkCutsize;
#  endif

  for (i=0;i<numRuns;++i)
    {
      if(i > 0)
	hgraph->prescribedVertexShuffle(mapToOrigVerts.getArray(),mapToOrigVerts.getArray(),comm); 
      
      /* create zhukov partition and refine using k-way refinement */            
      
      refiner.refineZhukovPartition(*hgraph, comm);
      refiner.releaseMemory();

      /* now use refined partition to guide coarsening + refinement */
      
      firstCutSize = hgraph->keepBestPartition();
      numIterations = 0;
      recordVCyclePartition(*hgraph,numIterations++);
        
      shuffleVCycleVertsByPartition(*hgraph, comm);
      hgraphs.push(hgraph);
      hEdgePercentiles.push(startPercentile);

      finerGraph = hgraph;
      accumulator = 1.0;
	  
      // ### 
      // coarsen the hypergraph
      // ###
      
      MPI_Barrier(comm);
      startTime = MPI_Wtime();
      
      do 
	{
	  hEdgePercentile = hEdgePercentiles.getTopElem();
	  restrCoarsener.setPercentile(hEdgePercentile);
	  coarseGraph = restrCoarsener.coarsen(*finerGraph,comm);	      
	  
	  if(coarseGraph) 
	    {	      
	      hEdgePercentiles.push(min(hEdgePercentile+percentileIncrement,100));
	      hgraphs.push(coarseGraph);
	      finerGraph = coarseGraph;
	    }
	}
	  
      while(coarseGraph);
      
      MPI_Barrier(comm);
      totCoaTime += (MPI_Wtime()-startTime);
      
      restrCoarsener.releaseMemory();
	  
      // ###
      // compute the initial partition
      // ###
      
      coarseGraph = hgraphs.pop();     
      coarseGraph->setNumberPartitions(0);
      coarseGraph->shiftVerticesToBalance(comm);
      
      MPI_Barrier(comm);
      startTime = MPI_Wtime();
      
      seqController.runSeqPartitioner(*coarseGraph,comm);            
      
      MPI_Barrier(comm);
      totSeqTime += (MPI_Wtime()-startTime);

      // ### 
      // uncoarsen the initial partition
      // ###
	  
      MPI_Barrier(comm);
      startTime = MPI_Wtime();
      accumulator = 1.0;
      
      while(hgraphs.getNumElem() > 0)
	{
	  coarseGraph->removeBadPartitions(keepPartitionsWithin*accumulator);
	  accumulator *= reductionInKeepThreshold;
	  
	  hEdgePercentile = hEdgePercentiles.pop();
	  finerGraph = hgraphs.pop();
	  
	  if(finerGraph == hgraph)		
	    shiftVCycleVertsToBalance(*finerGraph,comm);		
	  else		
	    finerGraph->shiftVerticesToBalance(comm);		      		
	  
	  finerGraph->projectPartitions(*coarseGraph,comm);
	  
#  ifdef DEBUG_CONTROLLER
	  finerGraph->checkPartitions(numTotalParts,maxPartWt,comm);
#  endif	  
	  if(randShuffBefRef)
	    {
	      if(hgraphs.getNumElem() == 0)
		finerGraph->randomVertexShuffle(mapToInterVerts.getArray(), comm);
	      else
		finerGraph->randomVertexShuffle(*(hgraphs.getTopElem()), comm);
	    }	  
	  
	  if(approxRefine)
	    refiner.setPercentile(hEdgePercentile);
	  
	  refiner.refine(*finerGraph,comm);
	  
#  ifdef DEBUG_CONTROLLER
	  finerGraph->checkPartitions(numTotalParts,maxPartWt,comm);
#  endif
	  DynaMem<ParaHypergraph>::deletePtr(coarseGraph);
	  
	  coarseGraph = finerGraph;	
	}     
	  
#  ifdef DEBUG_CONTROLLER
      assert(coarseGraph == hgraph);
      checkCutsize = coarseGraph->calcCutsize(numTotalParts,0,comm);	  	
#  endif
      MPI_Barrier(comm);
      totRefTime += (MPI_Wtime()-startTime);
      
      refiner.releaseMemory();
	  
      // ###
      // select the best partition           
      // ###
      
      secondCutSize = coarseGraph->keepBestPartition();
      diffInCutSize = firstCutSize - secondCutSize;	  
      
      if(dispOption > 1 && myRank == 0)
	{
	  out_stream << "\t ------ [" << numIterations << "] " << diffInCutSize << endl;
	}	  
      
#  ifdef DEBUG_CONTROLLER
      checkCutsize = coarseGraph->calcCutsize(numTotalParts,0,comm);	  
      assert(secondCutSize == checkCutsize);
#  endif
      
      if(diffInCutSize > 0 && diffInCutSize < minIterationGain)
	break;
      
      if(numIterations == limitOnCycles)
	break;
      
      if(firstCutSize > secondCutSize)
	{	     
	  recordVCyclePartition(*coarseGraph,numIterations++);
	  firstCutSize = secondCutSize;
	  vCycleGain += diffInCutSize;
	}
    }
  while(diffInCutSize > 0);
  
  gatherInVCyclePartition(*hgraph,firstCutSize,comm);
  updateMapToOrigVerts(comm);
  
  if(dispOption > 1 && myRank == 0)
    {
      out_stream << "\t ------ " << vCycleGain 
		 << " ------" << endl;
    }
      
#  ifdef DEBUG_CONTROLLER
      assert(hgraphs.getNumElem() == 0);
      assert(hgraph->getNumLocalVertices() == numOrigLocVerts);
      hgraph->checkPartitions(numTotalParts, maxPartWt, comm);      
#  endif      
      
      if(myRank == 0 && dispOption > 0)
	{
	  out_stream << "\nPRUN[" << i << "] = "
		     << firstCutSize << endl << endl;
	}

      totCutsizes += firstCutSize;
      
      if(firstCutSize < bestCutsize)
	{
	  bestCutsize = firstCutSize;
	  storeBestPartition(numOrigLocVerts,hgraph->getPartVectorArray(),comm);
	}
      
      if(firstCutSize > worstCutsize)
	worstCutsize = firstCutSize;
      
      // ### 
      // reset hypergraph 
      // ###

      resetStructs();                 





    }




}






#  endif
