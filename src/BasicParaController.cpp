

#  ifndef _BASIC_PARA_CONTROLLER_CPP
#  define _BASIC_PARA_CONTROLLER_CPP


// ### BasicParaController.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
// 
// HISTORY: 
// 
// 18/1/2005: Last Modified
//
// NOTES (18/1/2005) :
//
// - working on random shuffle before refinement.
//   first make work in the simple multilevel case
//   then extend idea to v-cycle refinement
//
// ###



#  include "BasicParaController.hpp"


BasicParaController::BasicParaController(ParaCoarsener& c, ParaRefiner& r, SeqController& ref, int rank, int nP, int percentile, int inc, int approxRef, ostream &out)
  : ParaController(c,r,ref,rank,nP,percentile,inc,approxRef,out)
{

}


BasicParaController::~BasicParaController()
{

}


void BasicParaController::dispParaControllerOptions() const
{
  switch(dispOption)
    {
    case SILENT:
      break;

    default:
      
      out_stream << "|--- PARA_CONTR (# parts = " << numTotalParts << "): " << endl
		 << "|- BASIC:"
		 << " pRuns = " << numParaRuns
		 << " kT = " << keepPartitionsWithin
		 << " rKT = " << reductionInKeepThreshold
		 << " appRef = " << approxRefine
		 << " wTF = " << writePartitionToFile
		 << " start %le " << startPercentile
		 << " %le inc " << percentileIncrement 
		 << endl << "|" << endl;
      break;
    }
}



void BasicParaController::runPartitioner(MPI_Comm comm)
{
  int hEdgePercentile;
  int cutSize;
  int percentCoarsening;
  int percentSequential;
  int percentRefinement;
  int percentOther;
  int i;

  double totStartTime;
 
#  ifdef DEBUG_CONTROLLER
  int checkCutsize;
#  endif

  Stack<int> hEdgePercentiles;
 
  numOrigLocVerts = hgraph->getNumLocalVertices();
 
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

  for (i=0;i<numParaRuns;++i) 
    {
#  ifdef MEM_CHECK
      write_log(myRank, "[begin run]: usage: %f", MemoryTracker::usage());
      Funct::printMemUse(myRank, "[begin run]");
#  endif      
      if(shuffled == 1) {    
	hgraph->randomVertexShuffle(mapToOrigVerts.getArray(), comm);
      }

      if(shuffled == 2) {
	hgraph->prescribedVertexShuffle(mapToOrigVerts.getArray(), shufflePartition.getArray(), comm);    	
      }

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
	  finerGraph->freeMemory();
	  
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

#  ifdef MEM_CHECK
      MPI_Barrier(comm);
      write_log(myRank, "[after coarsening]: usage: %f", MemoryTracker::usage());
      Funct::printMemUse(myRank, "[after coarsening]");
      MPI_Barrier(comm);
#  endif   
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

#  ifdef MEM_CHECK
      MPI_Barrier(comm);
      write_log(myRank, "[after seq partitioning]: usage: %f", MemoryTracker::usage());
      Funct::printMemUse(myRank, "[after seq partitioning]");
      MPI_Barrier(comm);
#  endif   

      MPI_Barrier(comm);
      startTime = MPI_Wtime();

      while(hgraphs.getNumElem() > 0)
	{
	  coarseGraph->removeBadPartitions(keepPartitionsWithin*accumulator);
	  accumulator *= reductionInKeepThreshold;
	  
	  finerGraph = hgraphs.pop();	         	    	     	      
	  hEdgePercentile = hEdgePercentiles.pop();

	  finerGraph->projectPartitions(*coarseGraph,comm);
	  
#  ifdef DEBUG_CONTROLLER
	  finerGraph->checkPartitions(numTotalParts,maxPartWt,comm);
#  endif		
	  if(randShuffBefRef)
	    {
	      if(hgraphs.getNumElem() == 0)
		finerGraph->randomVertexShuffle(mapToOrigVerts.getArray(), comm);
	      else
		finerGraph->randomVertexShuffle(*(hgraphs.getTopElem()), comm);
	    }	  
	  
	  if(approxRefine)
	    refiner.setPercentile(hEdgePercentile);

#  ifdef MEM_CHECK
	  write_log(myRank, "[before refineme]: usage: %f", MemoryTracker::usage());
	  Funct::printMemUse(myRank, "[before refinement]");
#  endif  
	  refiner.refine(*finerGraph,comm);

#  ifdef DEBUG_CONTROLLER
	  finerGraph->checkPartitions(numTotalParts,maxPartWt,comm);
#  endif	  	      	  

	  DynaMem<ParaHypergraph>::deletePtr(coarseGraph);	
	  coarseGraph = finerGraph;	
	}     
      
#  ifdef DEBUG_CONTROLLER
      assert(coarseGraph == hgraph);
#  endif

      MPI_Barrier(comm);
      totRefTime += (MPI_Wtime()-startTime);

      refiner.releaseMemory();

#  ifdef MEM_CHECK
      MPI_Barrier(comm);
      write_log(myRank, "[after refinement]: usage: %f", MemoryTracker::usage());
      Funct::printMemUse(myRank, "[after refinement]");
      MPI_Barrier(comm);
#  endif   

      /* select the best partition */          

      cutSize = coarseGraph->keepBestPartition();

#  ifdef DEBUG_CONTROLLER
      checkCutsize = coarseGraph->calcCutsize(numTotalParts,0,comm);
      assert(cutSize == checkCutsize);
#  endif
      
      if(cutSize < bestCutsize)
	{	  	 
	  bestCutsize = cutSize;
	  storeBestPartition(numOrigLocVerts,hgraph->getPartVectorArray(),comm);
	}

      if(cutSize > worstCutsize)
	worstCutsize = cutSize;
      
      totCutsizes += cutSize;
      
      /* free memory used by the para controller */

      resetStructs();

#  ifdef DEBUG_CONTROLLER
      assert(hgraphs.getNumElem() == 0);
#  endif
      
      if(myRank == 0 && dispOption > 0)
	{
	  out_stream << "\nPRUN[" << i << "] = "
		     << cutSize << endl << endl;
	}
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


void BasicParaController::resetStructs()
{
  hgraph->resetVectors(); 
  freeMemory();
}









#  endif
