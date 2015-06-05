
#  ifndef _VCYCLEALL_CONTROLLER_CPP
#  define _VCYCLEALL_CONTROLLER_CPP


// ### ParaVCycleAllController.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
// 
// HISTORY: 
// 
// 4/1/2005: Last Modified
//
// ###


#  include "ParaVCycleAllController.hpp"


ParaVCycleAllController::ParaVCycleAllController(ParaRestrCoarsener &rc, ParaCoarsener& c, ParaRefiner& r, SeqController& ref, int rank, int nP, int percentile, int inc, int approxRef, int limit, double limitAsPercent, ostream &out)
  : ParaVCycleController(rc,c,r,ref,rank,nP,percentile,inc,approxRef,limit,limitAsPercent,out)
{

}


ParaVCycleAllController::~ParaVCycleAllController()
{


}



void ParaVCycleAllController::runPartitioner(MPI_Comm comm)
{
  int i;
  int firstCutSize;
  int secondCutSize;
  int numInStack;
  int diffInCutSize;
  int vCycleIteration;
  int vCycleGain;
  int minIterationGain;
  
  int percentCoarsening;
  int percentSequential;
  int percentRefinement;
  int percentOther;
  int hEdgePercentile;

  /* experimental limit on number of v-cycle iterations */

  int doVCycle;
  int numOrigTotVerts = hgraph->getNumTotalVertices();

  double othAccumulator;
  double totStartTime;

  Stack<int> hEdgePercentiles;

  ParaHypergraph *coarseGraph;
  ParaHypergraph *interMedGraph;
  ParaHypergraph *finerGraph;

  initMapToOrigVerts();

  bestCutsize = LARGE_CONSTANT;
  worstCutsize = 0;
  totCutsizes = 0;
  
#  ifdef DEBUG_CONTROLLER
  int checkCutsize;
#  endif
  
  MPI_Barrier(comm);
  totStartTime = MPI_Wtime();
  
  for (i=0;i<numParaRuns;++i) 
    {     
      if(shuffled == 1)     
	hgraph->randomVertexShuffle(mapToOrigVerts.getArray(), comm);
      
      if(shuffled == 2)    	
	hgraph->prescribedVertexShuffle(mapToOrigVerts.getArray(), shufflePartition.getArray(), comm);    	

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
	  coarsener.setPercentile(hEdgePercentile);
	  coarseGraph = coarsener.coarsen(*finerGraph,comm);	  
	  
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

      coarsener.releaseMemory();
      
      // ###
      // compute the initial partition
      // ###

      MPI_Barrier(comm);
      startTime = MPI_Wtime();
      
      coarseGraph = hgraphs.pop();     
      seqController.runSeqPartitioner(*coarseGraph,comm);            
      
      MPI_Barrier(comm);
      totSeqTime += (MPI_Wtime()-startTime);

      // ### 
      // uncoarsen the initial partition
      // ###

      while(hgraphs.getNumElem() > 0)
	{
	  coarseGraph->removeBadPartitions(keepPartitionsWithin*accumulator);
	  accumulator *= reductionInKeepThreshold;	  
	  hEdgePercentile = hEdgePercentiles.pop();
	  finerGraph = hgraphs.pop();	    
	  
	  MPI_Barrier(comm);
	  startTime = MPI_Wtime();
  
	  projectVCyclePartition(*coarseGraph,*finerGraph,comm);
	  
#  ifdef DEBUG_CONTROLLER
	  finerGraph->checkPartitions(numTotalParts,maxPartWt,comm);
#  endif		  
	  if(approxRefine)
	    refiner.setPercentile(hEdgePercentile);
	  
	  refiner.refine(*finerGraph,comm);
	  refiner.releaseMemory();

#  ifdef DEBUG_CONTROLLER
	  finerGraph->checkPartitions(numTotalParts,maxPartWt,comm);
#  endif
	  	      
	  MPI_Barrier(comm);
	  totRefTime += (MPI_Wtime()-startTime);
	  
	  DynaMem<ParaHypergraph>::deletePtr(coarseGraph);

	  /* choose/reject v-cycle */
	  
	  if(finerGraph->getNumTotalVertices() > numOrigTotVerts/4)
	    doVCycle = 0;
	  else 
	    doVCycle = 1;

	  if(doVCycle)
	    {
	      // doing v-cycle
	      
	      vCycleIteration = 0;
	      vCycleGain = 0;
	      firstCutSize = finerGraph->keepBestPartition();
	      recordVCyclePartition(*finerGraph,vCycleIteration++);
	      
	      numInStack = hgraphs.getNumElem();
	      interMedGraph = finerGraph;
	      
	      if(dispOption > 1 && myRank == 0)
		{
		  out_stream << "\t ------ PARALLEL V-CYCLE CALL ------" << endl;
		}
	      
	      do
		{      	  	      
		  minIterationGain = static_cast<int>(floor(limAsPercentOfCut*firstCutSize));	      
		  
		  shuffleVCycleVertsByPartition(*finerGraph, comm);
		  hgraphs.push(finerGraph);
		  
		  othAccumulator = 1.0;
		  
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

		  while(hgraphs.getNumElem() > numInStack)
		    {
		      coarseGraph->removeBadPartitions(keepPartitionsWithin*othAccumulator);
		      othAccumulator *= reductionInKeepThreshold;
		      
		      hEdgePercentile = hEdgePercentiles.pop();
		      finerGraph = hgraphs.pop();
		      
		      if(finerGraph == interMedGraph)
			{		   
			  shiftVCycleVertsToBalance(*finerGraph,comm);
			}
		      else
			{
			  finerGraph->shiftVerticesToBalance(comm);		      
			}		  
		      
		      finerGraph->projectPartitions(*coarseGraph,comm);
		 
#  ifdef DEBUG_CONTROLLER
		      finerGraph->checkPartitions(numTotalParts,maxPartWt,comm);
#  endif	  	  	

		      if(randShuffBefRef)
			{
			  if(finerGraph == interMedGraph)
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
	      
		  MPI_Barrier(comm);
		  totRefTime += (MPI_Wtime()-startTime);

#  ifdef DEBUG_CONTROLLER
		  assert(coarseGraph == interMedGraph);
#  endif	      
		  refiner.releaseMemory();
	      
		  // ###
		  // select the best partition           
		  // ###
	      
		  secondCutSize = coarseGraph->keepBestPartition();

#  ifdef DEBUG_CONTROLLER
		  coarseGraph->checkPartitions(numTotalParts,maxPartWt,comm);
#  endif		  	      
		  diffInCutSize = firstCutSize - secondCutSize;	  
	      
		  if(dispOption > 1 && myRank == 0)
		    {
		      out_stream << "\t ------ [" << vCycleIteration << "] " << diffInCutSize << endl;
		    }	      
		  
		  if(diffInCutSize > 0 && diffInCutSize < minIterationGain)
		    break;	 
		  
		  if(vCycleIteration == limitOnCycles)
		    break;
		  
		  if(firstCutSize > secondCutSize)
		    {
		      recordVCyclePartition(*coarseGraph,vCycleIteration++);
		      firstCutSize = secondCutSize;		 
		      vCycleGain += diffInCutSize;
		    }
		}
	      while(diffInCutSize > 0);
	      
#  ifdef DEBUG_CONTROLLER
	      assert(coarseGraph == interMedGraph);
#  endif	  
	      
	      gatherInVCyclePartition(*coarseGraph,firstCutSize,comm);
	      
	      if(dispOption > 1 && myRank == 0)
		{
		  out_stream << "\t ------ " << vCycleGain 
			     << " ------" << endl;
		}
	    }
	  else
	    {
	      // not doing v-cycle
	      coarseGraph = finerGraph;
	    }

#  ifdef DEBUG_CONTROLLER
	  coarseGraph->checkPartitions(numTotalParts,maxPartWt,comm);
#  endif		  
	}     
      
#  ifdef DEBUG_CONTROLLER
      assert(coarseGraph == hgraph);
#  endif

      firstCutSize = coarseGraph->keepBestPartition();

      updateMapToOrigVerts(comm);

#  ifdef DEBUG_CONTROLLER
      checkCutsize = coarseGraph->calcCutsize(numTotalParts,0,comm);      
      assert(firstCutSize == checkCutsize);
#  endif
                  
      if(myRank == 0 && dispOption > 0)
	{
	  out_stream << "\nPRUN[" << i << "] " 
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
    
  MPI_Barrier(comm);
  totalTime = MPI_Wtime()-totStartTime;
  
  percentCoarsening = static_cast<int>(floor((totCoaTime/totalTime)*100));
  percentSequential = static_cast<int>(floor((totSeqTime/totalTime)*100));  
  percentRefinement = static_cast<int>(floor((totRefTime/totalTime)*100));
  percentOther = 100 - (percentCoarsening+percentSequential+percentRefinement);
    
  aveCutSize = static_cast<double>(totCutsizes) / numParaRuns;  
  
  if(myRank == 0 && dispOption > 0)
    {
      out_stream << endl
		 << " --- PARTITIONING SUMMARY ---" << endl
		 << "|" << endl
		 << "|--- Cutsizes statistics:" << endl
		 << "|" << endl
		 << "|-- BEST = " << bestCutsize << endl
		 << "|-- WORST = " << worstCutsize << endl
		 << "|-- AVE = " << aveCutSize << endl
		 << "|" << endl
		 << "|--- Time usage:" << endl
		 << "|" << endl
		 << "|-- TOTAL TIME = " << totalTime << endl
		 << "|-- AVE TIME = " << totalTime/numParaRuns << endl
		 << "|-- PARACOARSENING% = " << percentCoarsening << endl
		 << "|-- SEQPARTITIONING% = " << percentSequential << endl
		 << "|-- PARAREFINEMENT% = " << percentRefinement << endl
		 << "|-- OTHER% = " << percentOther << endl
		 << "|" << endl
		 << " ----------------------------" << endl;      
    }
}


void ParaVCycleAllController::printVCycleType() const
{
  out_stream << " type = ALL";
}
 



#  endif
