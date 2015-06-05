
#  ifndef _PARA_GREEDYKWAY_REFINER_CPP
#  define _PARA_GREEDYKWAY_REFINER_CPP


// ### ParaGreedyKwayRefiner.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
// 
// HISTORY: 
// 
// 3/2/2005: Last Modified
//
// ###


#  include "ParaGreedyKwayRefiner.hpp"


ParaGreedyKwayRefiner::ParaGreedyKwayRefiner(int rank, int nProcs, int nParts, int numVperP, int eExit, double lim, ostream &out)
: ParaRefiner(rank,nProcs,nParts,out)
{
  register int i;
  register int j;
  register int ij;

  // ###
  // create the move set structures
  // ###

  earlyExit = eExit;
  limit = lim;
  numTotVerticesMoved = 0;
  ij = numParts*numParts;

  moveSets.setLength(ij);
  moveSetData.setLength(Shiftl(ij,1));
  indexIntoMoveSetData.setLength(ij);
  numVerticesMoved.setLength(ij);
  
  for (i=0;i<ij;++i)    
    moveSets[i] = new FastDynaArray<int>;    

  for (i=0;i<numParts;++i) 
    {
      ij = i*numParts;

      for (j=0;j<numParts;++j)
	{
	  if(i == j)
	    {
	      indexIntoMoveSetData[ij+j] = -1;
	    }
	  else
	    {
	      indexIntoMoveSetData[ij+j] = Shiftl((ij+j),1);
	    }
	}
    }
  
  // ###
  // build tables
  // ###

  movementSets = new MovementSetTable(numParts, numProcs);

  numNeighParts.setLength(0);
  neighboursOfV.setLength(0);
  neighboursOfVOffsets.setLength(0);
  hEdgeVinPart.setLength(0);
  hEdgeVinPartOffsets.setLength(0);
  vertices.setLength(0);
  movedVertices.setLength(0);
  numPartsSpanned.setLength(0);
  spannedParts.setLength(0);
  locked.setLength(0);
  vertSeen.setLength(0);
}



ParaGreedyKwayRefiner::~ParaGreedyKwayRefiner()
{
  register int i;
  register int j;
  register int ij;

  DynaMem<MovementSetTable>::deletePtr(movementSets);

  for (i=0;i<numParts;++i) 
    {    
      ij = i*numParts;
      
      for (j=0;j<numParts;++j)
	DynaMem<FastDynaArray<int> >::deletePtr(moveSets[ij+j]);
    }
}


void ParaGreedyKwayRefiner::dispRefinementOptions() const
{
  switch(dispOption)
    {
    case SILENT:
      
      break;
      
    default:
      
      out_stream << "|--- PARA_REF: " << endl
		 << "|- PKWAY:"
		 << " eeL = " << limit
		 << " eExit = " << earlyExit
		 << endl << "|" << endl;
      break;
    }
}


void ParaGreedyKwayRefiner::releaseMemory()
{
  hEdgeWeight.setLength(0);
  hEdgeOffset.setLength(0);
  locPinList.setLength(0);

  vToHedgesOffset.setLength(0);
  vToHedgesList.setLength(0);
  allocHedges.setLength(0);

  numNeighParts.setLength(0);  
  neighboursOfV.setLength(0);
  neighboursOfVOffsets.setLength(0);
  hEdgeVinPart.setLength(0);
  hEdgeVinPartOffsets.setLength(0);
  vertices.setLength(0);
  movedVertices.setLength(0);
  seenVertices.setLength(0);
  numPartsSpanned.setLength(0);
  spannedParts.setLength(0);

  nonLocVerts.setLength(0);
  partIndices.setLength(0);
  indexIntoPartIndices.setLength(0);
  nonLocVToHedges.setLength(0);
  nonLocOffsets.setLength(0);

  toNonLocVerts.destroyTable();

  freeMemory();
}



void ParaGreedyKwayRefiner::initDataStructs(const ParaHypergraph &h, MPI_Comm comm)
{
  register int i;
  register int j;
  register int ij;

  initPartitionStructs(h,comm);

  movementSets->setMaxWt(maxPartWt);

  // ###
  // init data structures used in refinement
  // ###

  vertSeen.setLength(numLocalVertices);
  vertSeen.clear();

  locked.setLength(numLocalVertices);
  vertices.setLength(numLocalVertices);
  seenVertices.setLength(numLocalVertices);
  spannedParts.setLength(numParts);
  numPartsSpanned.setLength(numHedges);
  numNeighParts.setLength(numLocalVertices);
  neighboursOfV.setLength(numLocalVertices*numParts);
  neighboursOfVOffsets.setLength(numLocalVertices+1);

  hEdgeVinPart.setLength(numHedges*numParts);
  hEdgeVinPartOffsets.setLength(numHedges+1);

  j = 0;
  ij = numParts;

  for (i=0;i<=numLocalVertices;++i)
    {
      neighboursOfVOffsets[i] = j;
      j += ij; 
    }

  j = 0;
  
  for (i=0;i<=numHedges;++i)
    {
      hEdgeVinPartOffsets[i] = j;
      j += ij;
    }
}



void ParaGreedyKwayRefiner::resetDataStructs()
{
  toNonLocVerts.destroyTable();
  freeMemory();
}



void ParaGreedyKwayRefiner::setPartitioningStructs(int pNo, MPI_Comm comm)
{
#  ifdef DEBUG_REFINER
  assert(pNo >= 0 && pNo < numPartitions);
#  endif
  
  register int i;
  register int j;
  register int ij;

  int neighVertOffset;
  int hEdgeVertOffset;
  int vertOffset;
  int endOffset;
  int numSpanned;
  int vertex;
  int vPart;
  int v;
  int nonLocIndex;

#   ifdef DEBUG_REFINER
  assert(numParts > 1);
  assert(numNeighParts.getLength() > 0);
  assert(neighboursOfVOffsets.getLength() > 0);
  assert(neighboursOfV.getLength() > 0);  
  assert(hEdgeVinPartOffsets.getLength() > 0);
#  endif

  FastDynaArray<int> locPartWts(numParts);

  // ###
  // initialise current partition vector
  // ###

  currPVector = &partitionVector[partitionVectorOffsets[pNo]];
  if(numNonLocVerts == 0) currNonLocPVector = NULL;
  else currNonLocPVector = &partIndices[indexIntoPartIndices[pNo]];
  currPnumber = pNo;
  
  for (i=0;i<numParts;++i)    
    locPartWts[i] = 0;    

  // ###
  // initialise other data structures
  // first reset to zeros...
  // ###

  j = numLocalVertices;

  for (i=0;i<j;++i)
    {
      locPartWts[currPVector[i]] += vWeight[i];
      numNeighParts[i] = 0;
    }

  j = neighboursOfVOffsets[numLocalVertices];

  for (i=0;i<j;++i)    
      neighboursOfV[i] = 0;    
  
  j = numHedges;

  for (i=0;i<j;++i)    
      numPartsSpanned[i] = 0;    

  j = hEdgeVinPartOffsets[numHedges];
  
  for (i=0;i<j;++i)    
    hEdgeVinPart[i] = 0;    
  
  MPI_Allreduce(locPartWts.getArray(), partWeights.getArray(), numParts, MPI_INT, MPI_SUM, comm);

#  ifdef DEBUG_REFINER
  for (i=0;i<numParts;++i)
    assert(partWeights[i] <= maxPartWt);
#  endif

  for (i=0;i<numHedges;++i)
    {
      endOffset = hEdgeOffset[i+1]; 
      hEdgeVertOffset = hEdgeVinPartOffsets[i]; 
      numSpanned = 0;

      for (j=hEdgeOffset[i];j<endOffset;++j)
	{
	  ij = locPinList[j];

	  if(ij >= minVertexIndex && ij < maxVertexIndex)
	    {
	      vPart = currPVector[ij-minVertexIndex];
	    }
	  else
	    {
	      nonLocIndex = toNonLocVerts.getVal(ij); 
#  ifdef DEBUG_REFINER
	      assert(nonLocIndex >= 0 && nonLocIndex < numNonLocVerts);
#  endif
	      vPart = currNonLocPVector[nonLocIndex];      
	    }
#  ifdef DEBUGU_REFINER
	  assert(vPart >= 0 && vPart < numParts);
#  endif
	  ij = hEdgeVertOffset+vPart;
	  
	  if(hEdgeVinPart[ij] == 0)
	    {
	      spannedParts[numSpanned++] = vPart;
	    }
	  
	  ++hEdgeVinPart[ij];
	}

#  ifdef DEBUG_REFINER
      int hEdgeLen = 0;
      for (int ijk=0;ijk<numParts;++ijk)
	hEdgeLen += hEdgeVinPart[hEdgeVertOffset+ijk];
      assert(hEdgeLen == hEdgeOffset[i+1]-hEdgeOffset[i]);
#  endif
      
      numPartsSpanned[i] = numSpanned;

      // ###
      // update the vertex neighbour structs
      // ###
      
      for (j=hEdgeOffset[i];j<endOffset;++j)
	{
	  vertex = locPinList[j];

	  // ###
	  // only consider if vertex a local vertex
	  // ###

	  if(vertex >= minVertexIndex && vertex < maxVertexIndex)
	    {
	      v = vertex-minVertexIndex;
	      vertOffset = neighboursOfVOffsets[v];
 
	      for (ij=0;ij<numSpanned;++ij)
		{	
		  neighVertOffset = vertOffset+spannedParts[ij]; 

		  if(neighboursOfV[neighVertOffset] == 0)
		    ++numNeighParts[v];
	      
		  neighboursOfV[neighVertOffset] = 1;
		}        
	    }
	}      
    }

#  ifdef DEBUG_REFINER
  nonLocVertCheck(); 
  sanityHedgeCheck();
#  endif
}



void ParaGreedyKwayRefiner::refine(ParaHypergraph &h, MPI_Comm comm)
{
  initDataStructs(h,comm);

  register int i;

  int gain;
  int lastGain = 0;
  int totalGain;
  int newCutsize;
  int numPasses;

  for (i=0;i<numPartitions;++i)
    {            
      setPartitioningStructs(i,comm);

      totalGain = 0;
      numPasses = 0;

      if(dispOption > 1 && myRank == 0)
	{
	  out_stream << "\t[" << i << "] ";
	}
      
      do
	{
	  newCutsize = runGreedyKwayRefinement(h,i,comm);
	  gain = partitionCuts[i] - newCutsize;

	  if(dispOption > 1 && myRank == 0)
	    {
	      out_stream << gain << " ";
	    }

	  if(gain < 0)
	    {
	      takeBackPassMoves();	    
	    }
	  else
	    {
	      partitionCuts[i] = newCutsize;
	      totalGain += gain;
	      ++numPasses;
	    }

	  if(earlyExit && numPasses>1 && gain<lastGain)
	    break;

	  lastGain = gain;
	}
      while(gain > 0);      

      if(dispOption > 1 && myRank == 0)
	{	  
	  out_stream << "| " << numPasses << " "
		     << totalGain << " "
		     << partitionCuts[i] << endl;
	}
    }  

  resetDataStructs();
}



int ParaGreedyKwayRefiner::runGreedyKwayRefinement(ParaHypergraph &h, int pNo, MPI_Comm comm)
{
  register int i;
  
  numTotVerticesMoved = 0;
  locked.clear();
  
  for (i=0;i<2;++i)
    {      
      if(myRank == ROOT_PROC)
	{
	  // ###
	  // init movement sets
	  // ###

	  movementSets->initPartWeights(partWeights.getArray(), numParts);	
	}

      doGreedyPass(i,comm);      
      manageBalanceConstraint(comm);
      updateVertexMoveInfo(comm);
    }
  
  if(currPercentile == 100)
    return (computeCutsize(comm));
  else
    {
#  ifdef DEBUG_REFINER
      //assert(currPercentile > 0 && currPercentile < 100);
#  endif
      return (h.calcCutsize(numParts,pNo,comm));
    }
}


int ParaGreedyKwayRefiner::doGreedyPass(int lowToHigh, MPI_Comm comm)
{      
  int i;
  int j;
  int ij;
  int v;
  int sP;
  int tmp;
  int gain;
  int prod;
  int vGain;
  int posGain;
  int hEdge;
  int bestMove;
  int randomNum;
  int vertexWt;
  int vertOffset;
  int vNeighOffset;
  int hEdgeOff;
  int neighOfVOffset;
  int numNonPos = 0;
  int limNonPosMoves = static_cast<int>(ceil(limit*static_cast<double>(numLocalVertices)));;

  double currImbalance = 0;
  double bestImbalance;
  double posImbalance;

  for (i=0;i<numParts;++i) 
    {
      prod = numParts*i;

      for (j=0;j<numParts;++j) 
	{
	  if(i != j) 
	    {	
	      ij = prod+j;
	
	      numVerticesMoved[ij] = 0;
	      v = indexIntoMoveSetData[ij];
	      
	      moveSetData[v] = 0;
	      moveSetData[v+1] = 0;
	    }
	}
    }
  
  for (i=0;i<numLocalVertices;++i)    
      vertices[i] = i;    
  
  for (i=0;i<numParts;++i)    
      currImbalance += fabs(partWeights[i]-avePartWt);    

  i = numLocalVertices;
  gain = 0;

  do 
    {
      randomNum = RANDOM(0,i);
      v = vertices[randomNum];

      if(!locked(v))
	{
	  vertexWt = vWeight[v];
	  sP = currPVector[v];
	  vGain = 0;
	  bestImbalance = currImbalance;
	  bestMove = -1;
	  
	  if(numNeighParts[v] > 1) 
	    {
	      vNeighOffset = neighboursOfVOffsets[v];
	      
	      for (j=0;j<numParts;++j)
		{
		  if(neighboursOfV[vNeighOffset+j] > 0 && ((lowToHigh && j>sP) || (!lowToHigh && j<sP)))
		    {
		      if(partWeights[j]+vertexWt <= maxPartWt)
			{
			  posGain = 0;
			  vertOffset = vToHedgesOffset[v+1];
			  
			  for (ij=vToHedgesOffset[v];ij<vertOffset;++ij)
			    {
			      hEdge = vToHedgesList[ij];
#  ifdef DEBUG_REFINER
			      assert(hEdge >= 0 && hEdge < numHedges);
#  endif			      
			      hEdgeOff = hEdgeVinPartOffsets[hEdge];
			      
			      if(hEdgeVinPart[hEdgeOff+sP] == 1)
				posGain += hEdgeWeight[hEdge];
			      
			      if(hEdgeVinPart[hEdgeOff+j] == 0)
				posGain -= hEdgeWeight[hEdge];
			    }
			  
			  posImbalance = currImbalance+fabs(partWeights[sP]-(vertexWt+avePartWt));
			  posImbalance += fabs((partWeights[j]+vertexWt)-avePartWt);
			  posImbalance -= fabs(partWeights[sP] - avePartWt);
			  posImbalance -= fabs(partWeights[j] - avePartWt);
			  
			  if((posGain > vGain) || (posGain == vGain && posImbalance < bestImbalance))
			    {
			      vGain = posGain;
			      bestMove = j;
			      bestImbalance = posImbalance;
			    }
			}
		    }
		}
	      
	      if(bestMove != -1)
		{	  
#  ifdef DEBUG_REFINER
		  assert(bestMove >= 0 && bestMove < numParts);
#  endif		  
		  vertOffset = vToHedgesOffset[v+1];
		  neighOfVOffset = neighboursOfVOffsets[v];
		  
		  // ###
		  // update the moved vertices' stats
		  // ###
		  
		  if(neighboursOfV[neighOfVOffset+bestMove] == 0)		    
		    ++numNeighParts[v];		    
		  
		  neighboursOfV[neighOfVOffset+bestMove] = 1;
		  neighboursOfV[neighOfVOffset+sP] = 0;	      	     
		  
		  for (j=vToHedgesOffset[v];j<vertOffset;++j)
		    {
		      // ###
		      // update the hyperedge stats: (vInPart etc.)
		      // ###
		      
		      hEdge = vToHedgesList[j];
#  ifdef DEBUG_REFINER
		      assert(hEdge >= 0 && hEdge < numHedges);
#  endif
		      hEdgeOff = hEdgeVinPartOffsets[hEdge];
#  ifdef DEBUG_REFINER
		      int hEdgeLen = 0;
		      for (int ijk=0;ijk<numParts;++ijk)
			hEdgeLen += hEdgeVinPart[hEdgeOff+ijk];
		      assert(hEdgeLen == hEdgeOffset[hEdge+1]-hEdgeOffset[hEdge]);
#  endif		      		      
		      if(hEdgeVinPart[hEdgeOff+bestMove] == 0)
			{
			  ++numPartsSpanned[hEdge];
			}

		      --hEdgeVinPart[hEdgeOff+sP];		
		      ++hEdgeVinPart[hEdgeOff+bestMove];
		      
#  ifdef DEBUG_REFINER
		      assert(hEdgeVinPart[hEdgeOff+sP] >= 0);
		      assert(hEdgeVinPart[hEdgeOff+bestMove] > 0);
#  endif
		      
		      if(hEdgeVinPart[hEdgeOff+sP] > 0)
			{
			  neighboursOfV[neighOfVOffset+sP] = 1;
			}
		      else
			{
			  --numPartsSpanned[hEdge];
			}		  
		    }
		  
		  if(neighboursOfV[neighOfVOffset+sP] == 0)
		    {
		      --numNeighParts[v];
		    }
		  
		  // ###
		  // update the adj vertices stats: 
		  // (num neighbours in part etc.)
		  // ###
		  
		  updateAdjVertStatus(v,sP,bestMove);	      
		  
		  // ###
		  // update other structs
		  // ###
		  
		  locked.set1(v);
		  currPVector[v] = bestMove;
		  partWeights[sP] -= vertexWt;
		  partWeights[bestMove] += vertexWt;
		  currImbalance = bestImbalance;
		  
		  // ###
		  // update the gain...
		  // ###
		  
		  gain += vGain;
		  
		  if(vGain <= 0) ++numNonPos;
		  else numNonPos = 0;
		  
		  // ###
		  // update the movement set structures
		  // ###

		  ij = sP*numParts+bestMove;

#  ifdef DEBUG_REFINER
		  assert(ij < numParts*numParts);
		  assert(indexIntoMoveSetData[ij] != -1);
		  assert(indexIntoMoveSetData[ij] >= 0);
		  assert(indexIntoMoveSetData[ij] < Shiftl(numParts*numParts,1));
#  endif		  		  
		  moveSets[ij]->assign(numVerticesMoved[ij]++, v+minVertexIndex);		  
		  moveSetData[indexIntoMoveSetData[ij]] += vGain;
		  moveSetData[indexIntoMoveSetData[ij]+1] += vertexWt;
		  
		  if(limit < 1.0)
		    {
		      if(numNonPos > limNonPosMoves)
			break;
		    }		  
		}
	    }
	}

      fswap(vertices[randomNum], vertices[i-1], tmp);
      --i;
      
    }
  while(i > 0);
  
#  ifdef DEBUG_REFINER
  sanityHedgeCheck();
#  endif

  return gain;
}



int ParaGreedyKwayRefiner::computeCutsize(MPI_Comm comm)
{  
  register int i;
  register int ij;
  register int locCut = 0;

  int totalCut;
  
  for (i=0;i<numAllocHedges;++i)
    {            	
      ij = allocHedges[i];
#  ifdef DEBUG_REFINER
      assert(numPartsSpanned[ij] > 0 && numPartsSpanned[ij] <= numParts);
      assert(hEdgeWeight[ij] > 0);
#  endif
      locCut += ((numPartsSpanned[ij]-1)*hEdgeWeight[ij]);	
    }
  
  MPI_Allreduce(&locCut, &totalCut, 1, MPI_INT, MPI_SUM, comm);
  
  return totalCut;
}




void ParaGreedyKwayRefiner::manageBalanceConstraint(MPI_Comm comm) 
{

  register int i;
  register int j;
  register int ij; 

  int prod;
  int arrayLen;
  int numToSend;
  int totToRecv;
  int indexIntoMoveSets;

  int *movesLengths;
  int *array;

  FastDynaArray<int> **moves;
  
  numToSend = 0;
    
  for (i=0;i<numParts;++i) 
    {
      prod = i*numParts;

      for (j=0;j<numParts;++j)
	{ 
	  if(i != j) 
	    {	  
	      ij = indexIntoMoveSetData[prod+j];
	  
	      if(moveSetData[ij+1] > 0) 
		{	    
#  ifdef DEBUG_REFINER
		  assert(numVerticesMoved[prod+j] > 0);
#  endif	    
		  sendArray.assign(numToSend++, i);
		  sendArray.assign(numToSend++, j);
		  sendArray.assign(numToSend++, moveSetData[ij]);
		  sendArray.assign(numToSend++, moveSetData[ij+1]);	  
		}
	      else
		{
#  ifdef DEBUG_REFINER
		  if(numVerticesMoved[prod+j] > 0)
		    out_stream << "p[" << myRank << "] moved vertices with weight zero" << endl;
#  endif	    
		}
	    }
	}
    }
  
  MPI_Gather(&numToSend, 1, MPI_INT, recvLens.getArray(), 1, MPI_INT, ROOT_PROC, comm);
  
  if(myRank == ROOT_PROC) 
    {    
      ij = 0;

      for (i=0;i<numProcs;++i) 
	{      
	  recvDispls[i] = ij;
	  ij += recvLens[i];
	}
      
      receiveArray.setLength(ij);
      totToRecv = ij;
    }
  
  MPI_Gatherv(sendArray.getArray(), numToSend, MPI_INT, receiveArray.getArray(), recvLens.getArray(), recvDispls.getArray(), MPI_INT, ROOT_PROC, comm);
  
  if(myRank == ROOT_PROC) 
    {	
      for (i=0;i<numProcs;++i) 
	{      
	  ij = recvDispls[i];
	  
	  if(recvLens[i] > 0)
	    {
	      movementSets->completeProcSets(i,recvLens[i], &(receiveArray[ij])); 
	    }
	}

      movementSets->computeRestoringArray();
    }
  
  moves = movementSets->getRestoringMoves();
  movesLengths = movementSets->getRestoringMovesLens();    
  
  MPI_Scatter(movesLengths, 1, MPI_INT, &totToRecv, 1, MPI_INT, ROOT_PROC, comm);
  
  if(myRank == ROOT_PROC) 
    {      
      ij = 0;
      
      for (i=0;i<numProcs;++i) 
	{      
	  arrayLen = movesLengths[i];
	  sendDispls[i] = ij;
	  array = moves[i]->getArray();	 
	  
	  for (j=0;j<arrayLen;++j) 	    	
	      sendArray.assign(ij++, array[j]);	     	    
	}
    }
  
  receiveArray.setLength(totToRecv);
  
  MPI_Scatterv(sendArray.getArray(), movesLengths, sendDispls.getArray(), MPI_INT, receiveArray.getArray(), totToRecv, MPI_INT, ROOT_PROC, comm);
  
  // ###
  // now everyone has been told which sets of moves to 
  // take back in order to maintain balance condition
  // take back local moves as directed by root processor
  
  // note that can decrease the size of the take-back array
  // since we are not making use of the weight of the set
  // to be taken back
  // ###
  
  ij = 0;
      
  while(ij < totToRecv) 
    {    
      i = receiveArray[ij];
      j = receiveArray[ij+1];

      indexIntoMoveSets = i*numParts+j;      
      unmakeMoves(indexIntoMoveSets, i, j);
    
      ij += 2;    
    }
  
#  ifdef DEBUG_REFINER
  FastDynaArray<int> partWts(numParts);
  
  for (ij=0; ij<numParts; ++ij)
    partWts[ij] = 0;
  
  for (ij=0; ij<numLocalVertices; ++ij)
    partWts[partitionVector[ij]] += vWeight[ij];
  
  MPI_Allreduce(partWts.getArray(), partWeights.getArray(), numParts, MPI_INT, MPI_SUM, comm);
  
  for (ij=0;ij<numParts;++ij)
    {
      assert(partWts[ij] <= maxPartWt);
    }
#  endif
  
  if(myRank == ROOT_PROC) 
    {    
      array = movementSets->getPartWeightsArray();

      for (ij=0;ij<numParts;++ij)
	{
	  partWeights[ij] = array[ij];
	}
    }
  
  MPI_Bcast(partWeights.getArray(), numParts, MPI_INT, ROOT_PROC, comm);

#  ifdef DEBUG_REFINER
  for (i=0;i<numParts;++i)    
    assert(partWeights[i] <= maxPartWt);    
#  endif
}



void ParaGreedyKwayRefiner::takeBackPassMoves()
{
#  ifdef DEBUG_REFINER
  assert(And(numTotVerticesMoved,0x1) == 0);
#  endif

  register int i;
  register int v;
  register int vPart;
  
  for (i=0;i<numTotVerticesMoved;i+=2)
    {
      v = movedVertices[i];
      vPart = movedVertices[i+1];

#  ifdef DEBUG_REFINER
      assert(vPart >= 0 && vPart < numParts);
#  endif

      currPVector[v-minVertexIndex] = vPart;
    }       
}



void ParaGreedyKwayRefiner::updateVertexMoveInfo(MPI_Comm comm)
{
#  ifdef DEBUG_REFINER
  sanityHedgeCheck();
#  endif

  register int i;
  register int j;
  register int ij;
  register int ijk;

  int index;
  int prod;
  int vert;
  int totToSend;
  int totToRecv;
  int minIndex;
  int v;
  int vertexPart;
  int newVertexPart;
  int locVertIndex;
  int hEdge;
  int endOffset;
  int neighOfVOffset;
  int othVOffset;
  int othHedge;
  int hEdgeOff;
  int numVerticesSeen;
  int nonLocIdx;
  int endNonLocOffset;

  int *array;

  totToSend = 0;
  for (i=0;i<numParts;++i) 
    {
      prod = i*numParts;

      for (j=0;j<numParts;++j)
	{ 
	  if(i != j) 
	    {	
	      index = prod + j;
	      endOffset = numVerticesMoved[index];
	      array = moveSets[index]->getArray();

	      for (ij=0;ij<endOffset;++ij)
		{
		  v = array[ij];
#  ifdef DEBUG_REFINER	
		  assert(v >= minVertexIndex && v < maxVertexIndex);
		  assert(currPVector[v-minVertexIndex] == j);
#  endif	  
		  sendArray.assign(totToSend++, v);
		  sendArray.assign(totToSend++, j);
		  
		  movedVertices.assign(numTotVerticesMoved++, v);
		  movedVertices.assign(numTotVerticesMoved++, i);
		}
	    }
	}
    }

  MPI_Allgather(&totToSend, 1, MPI_INT, recvLens.getArray(), 1, MPI_INT, comm);
  
  ij = 0;
  for (i=0;i<numProcs;++i) 
    {
      recvDispls[i] = ij;
      ij += recvLens[i];
    }

  receiveArray.setLength(ij);
  totToRecv = ij;

  MPI_Allgatherv(sendArray.getArray(), totToSend, MPI_INT, receiveArray.getArray(), recvLens.getArray(), recvDispls.getArray(), MPI_INT, comm);
  
  // ###
  // now go through the moved vertices and update local structures
  // ###

  minIndex = recvDispls[myRank];
 
#  ifdef DEBUG_REFINER
  assert(And(totToRecv,0x1) == 0);
#  endif 

  i = 0;

  while(i<minIndex)
    {     
      v = receiveArray[i];

#  ifdef DEBUG_REFINER
      assert(v < minVertexIndex || v >= maxVertexIndex);
#  endif
      nonLocIdx = toNonLocVerts.getCareful(v);

      if(nonLocIdx >= 0)
	{
#  ifdef DEBUG_REFINER
	  assert(nonLocIdx < numNonLocVerts);
#  endif        
	  vertexPart = currNonLocPVector[nonLocIdx];
	  newVertexPart = receiveArray[i+1];
	  endNonLocOffset = nonLocOffsets[nonLocIdx+1];

#  ifdef DEBUG_REFINER	
	  assert(newVertexPart >= 0 && newVertexPart < numParts);
	  assert(vertexPart >= 0 && vertexPart < numParts);
	  assert(vertexPart != newVertexPart);
#  endif
	  
	  // ###
	  // first update the hyperedge vInPart structure
	  // ###	
	  
	  for (j=nonLocOffsets[nonLocIdx];j<endNonLocOffset;++j)
	    {
	      hEdge = nonLocVToHedges[j];

#  ifdef DEBUG_REFINER		   
	      assert(hEdge >= 0 && hEdge < numHedges);
#  endif
	      hEdgeOff = hEdgeVinPartOffsets[hEdge];

#  ifdef DEBUG_REFINER
	      int hEdgeLen = 0;
	      for (int ijkl=0;ijkl<numParts;++ijkl)
		hEdgeLen += hEdgeVinPart[hEdgeOff+ijkl];
	      assert(hEdgeLen == hEdgeOffset[hEdge+1]-hEdgeOffset[hEdge]);
	      assert(hEdgeVinPart[hEdgeOff+vertexPart] > 0);
#  endif
	      if(hEdgeVinPart[hEdgeOff+newVertexPart] == 0)
		{
		  ++numPartsSpanned[hEdge];
		}
	      
	      --hEdgeVinPart[hEdgeOff+vertexPart];		
	      ++hEdgeVinPart[hEdgeOff+newVertexPart];
	      
#  ifdef DEBUG_REFINER
	      assert(hEdgeVinPart[hEdgeOff+newVertexPart] > 0);
	      assert(hEdgeVinPart[hEdgeOff+vertexPart] >= 0);
#  endif
	      
	      if(hEdgeVinPart[hEdgeOff+vertexPart] == 0)
		--numPartsSpanned[hEdge];				      
	    }
	  
	  // ###
	  // update the data for local adjacent vertices
	  // ###

	  numVerticesSeen = 0;
	  
	  for (j=nonLocOffsets[nonLocIdx];j<endNonLocOffset;++j)
	    {      
	      hEdge = nonLocVToHedges[j];      
#  ifdef DEBUG_REFINER
	      assert(hEdge >= 0 && hEdge < numHedges);
#  endif
	      hEdgeOff = hEdgeOffset[hEdge+1];
	      
	      for (ij=hEdgeOffset[hEdge];ij<hEdgeOff;++ij)
		{
		  vert = locPinList[ij];
		  
		  // ###
		  // if the adjacent vertex is a local vertex
		  // ###
		  
		  if(vert >= minVertexIndex && vert < maxVertexIndex)
		    {
		      locVertIndex = vert-minVertexIndex;
		      
		      if(vertSeen(locVertIndex) == 0)
			{	      				  
			  neighOfVOffset = neighboursOfVOffsets[locVertIndex];
			  
			  if(neighboursOfV[neighOfVOffset+newVertexPart] == 0)
			    {
			      ++numNeighParts[locVertIndex];
			    }
			  
			  neighboursOfV[neighOfVOffset+newVertexPart] = 1;
			  
			  if(currPVector[locVertIndex] != vertexPart)
			    {
			      neighboursOfV[neighOfVOffset+vertexPart] = 0;
			      
			      // sorting out this now...go thro each hedge
			      // of neighbouring vertex and check if its
			      // connected to part bestDP to determine
			      // whether to set to 0 or 1!!!
			      
			      othVOffset = vToHedgesOffset[locVertIndex+1];	   
			      
			      for (ijk=vToHedgesOffset[locVertIndex];ijk<othVOffset;++ijk)
				{
				  othHedge = vToHedgesList[ijk];
				  
				  if(hEdgeVinPart[hEdgeVinPartOffsets[othHedge]+vertexPart] > 0)
				    {
				      neighboursOfV[neighOfVOffset+vertexPart] = 1;
				      break;
				    }
				}
			      
			      if(neighboursOfV[neighOfVOffset+vertexPart] == 0)
				{
				  --numNeighParts[locVertIndex];				  
				}
			    }
			  
			  vertSeen.set1(locVertIndex);
			  seenVertices[numVerticesSeen++] = locVertIndex;		  
			}
		    }
		  
		  // ###
		  // if the adjacent vertex is not a local vertex
		  // we don't update anything
		  // ###
		} 
	    }
	  
	  currNonLocPVector[nonLocIdx] = newVertexPart;

	  // ###
	  // restore the 'seen' vertices structure
	  // ###
	  
	  for (j=0;j<numVerticesSeen;++j)		
	    {
	      vertSeen.set0(seenVertices[j]);
	    }	  
	  
	  // ###
	  // end updating structures after remote vertex move
	  // ###
	}
      
      i += 2;
    }      
    
  i += totToSend;
  
  while(i<totToRecv)
    {
      v = receiveArray[i];
            
#  ifdef DEBUG_REFINER
      assert(v < minVertexIndex || v >= maxVertexIndex);
#  endif      
      
      nonLocIdx = toNonLocVerts.getCareful(v);

      if(nonLocIdx >= 0)
	{
#  ifdef DEBUG_REFINER
	  assert(nonLocIdx < numNonLocVerts);
#  endif	
	  vertexPart = currNonLocPVector[nonLocIdx];
	  newVertexPart = receiveArray[i+1];
	  endNonLocOffset = nonLocOffsets[nonLocIdx+1];

#  ifdef DEBUG_REFINER	
	  assert(newVertexPart >= 0 && newVertexPart < numParts);
	  assert(vertexPart >= 0 && vertexPart < numParts);
	  assert(vertexPart != newVertexPart);
#  endif
	  
	  // ###
	  // first update the hyperedge vInPart structure
	  // ###

	  for (j=nonLocOffsets[nonLocIdx];j<endNonLocOffset;++j)
	    {
	      hEdge = nonLocVToHedges[j];

#  ifdef DEBUG_REFINER		   
	      assert(hEdge >= 0 && hEdge < numHedges);
#  endif
	      hEdgeOff = hEdgeVinPartOffsets[hEdge];
#  ifdef DEBUG_REFINER
	      int hEdgeLen = 0;
	      for (int ijkl=0;ijkl<numParts;++ijkl)
		hEdgeLen += hEdgeVinPart[hEdgeOff+ijkl];
	      assert(hEdgeLen == hEdgeOffset[hEdge+1]-hEdgeOffset[hEdge]);
	      assert(hEdgeVinPart[hEdgeOff+vertexPart] > 0);
#  endif
	      if(hEdgeVinPart[hEdgeOff+newVertexPart] == 0)
		{
		  ++numPartsSpanned[hEdge];
		}
	      
	      --hEdgeVinPart[hEdgeOff+vertexPart];		
	      ++hEdgeVinPart[hEdgeOff+newVertexPart];
	      
#  ifdef DEBUG_REFINER
	      assert(hEdgeVinPart[hEdgeOff+newVertexPart] > 0);
	      assert(hEdgeVinPart[hEdgeOff+vertexPart] >= 0);
#  endif	      
	      if(hEdgeVinPart[hEdgeOff+vertexPart] == 0)
		--numPartsSpanned[hEdge];				      
	    }
	  
	  // ###
	  // update the data for local adjacent vertices
	  // ###
	  
	  numVerticesSeen = 0;

	  for (j=nonLocOffsets[nonLocIdx];j<endNonLocOffset;++j)
	    {      
	      hEdge = nonLocVToHedges[j];
#  ifdef DEBUG_REFINER
	      assert(hEdge >= 0 && hEdge < numHedges);
#  endif
	      hEdgeOff = hEdgeOffset[hEdge+1];
	      
	      for (ij=hEdgeOffset[hEdge];ij<hEdgeOff;++ij)
		{
		  vert = locPinList[ij];
		  
		  // ###
		  // if the adjacent vertex is a local vertex
		  // ###
		  
		  if(vert >= minVertexIndex && vert < maxVertexIndex)
		    {
		      locVertIndex = vert-minVertexIndex;
		      
		      if(vertSeen(locVertIndex) == 0)
			{	      				  
			  neighOfVOffset = neighboursOfVOffsets[locVertIndex];
			  
			  if(neighboursOfV[neighOfVOffset+newVertexPart] == 0)
			    {
			      ++numNeighParts[locVertIndex];
			    }
			  
			  neighboursOfV[neighOfVOffset+newVertexPart] = 1;
			  
			  if(currPVector[locVertIndex] != vertexPart)
			    {
			      neighboursOfV[neighOfVOffset+vertexPart] = 0;
			      
			      // sorting out this now...go thro each hedge
			      // of neighbouring vertex and check if its
			      // connected to part bestDP to determine
			      // whether to set to 0 or 1!!!
			      
			      othVOffset = vToHedgesOffset[locVertIndex+1];	   
			      
			      for (ijk=vToHedgesOffset[locVertIndex];ijk<othVOffset;++ijk)
				{
				  othHedge = vToHedgesList[ijk];
				  
				  if(hEdgeVinPart[hEdgeVinPartOffsets[othHedge]+vertexPart] > 0)
				    {
				      neighboursOfV[neighOfVOffset+vertexPart] = 1;
				      break;
				    }
				}
			      
			      if(neighboursOfV[neighOfVOffset+vertexPart] == 0)
				{
				  --numNeighParts[locVertIndex];				  
				}
			    }
			  
			  vertSeen.set1(locVertIndex);
			  seenVertices[numVerticesSeen++] = locVertIndex;		  
			}
		    }
		  
		  // ###
		  // if the adjacent vertex is not a local vertex
		  // we don't update anything
		  // ###
		} 
	    }
	  	
	  currNonLocPVector[nonLocIdx] = newVertexPart;

	  // ###
	  // restore the 'seen' vertices structure
	  // ###
	  
	  for (j=0;j<numVerticesSeen;++j)			    
	      vertSeen.set0(seenVertices[j]);		 
	  
	  // ###
	  // end updating structures after remote vertex move
	  // ###
	}
      
      i += 2;
    }

#  ifdef DEBUG_REFINER
  sanityHedgeCheck();
#  endif
}




void ParaGreedyKwayRefiner::updateAdjVertStatus(int v, int sP, int bestMove)
{
  int i;
  int j;
  int ij;
  
  int numVerticesSeen;
  int vertOffset;
  int vert;
  int locVertIndex;
  int hEdge;
  int hEdgeOff;
  int neighOfVOffset;
  int othVOffset;
  int othHedge;
  
  vertOffset = vToHedgesOffset[v+1];
  numVerticesSeen = 0;
  
  for (i=vToHedgesOffset[v];i<vertOffset;++i)
    {      
      hEdge = vToHedgesList[i];      
      hEdgeOff = hEdgeOffset[hEdge+1];
      
      for (j=hEdgeOffset[hEdge];j<hEdgeOff;++j)
	{
	  vert = locPinList[j];
	  
	  // ###
	  // if the adjacent vertex is a local vertex
	  // ###
	  
	  if(vert >= minVertexIndex && vert < maxVertexIndex)
	    {
	      locVertIndex = vert-minVertexIndex;

	      if(locVertIndex != v && vertSeen(locVertIndex) == 0)
		{	      				  
		  neighOfVOffset = neighboursOfVOffsets[locVertIndex];
		  
		  if(neighboursOfV[neighOfVOffset+bestMove] == 0)
		    {
		      ++numNeighParts[locVertIndex];
		    }
		  
		  neighboursOfV[neighOfVOffset+bestMove] = 1;
		  
		  if(currPVector[locVertIndex] != sP)
		    {
		      neighboursOfV[neighOfVOffset+sP] = 0;
		  
		      // sorting out this now...go thro each hedge
		      // of neighbouring vertex and check if its
		      // connected to part bestDP to determine
		      // whether to set to 0 or 1!!!
		      
		      othVOffset = vToHedgesOffset[locVertIndex+1];	   
		      
		      for (ij=vToHedgesOffset[locVertIndex];ij<othVOffset;++ij)
			{
			  othHedge = vToHedgesList[ij];
			  
			  if(hEdgeVinPart[hEdgeVinPartOffsets[othHedge]+sP] > 0)
			    {
			      neighboursOfV[neighOfVOffset+sP] = 1;
			      break;
			    }
			}
		      
		      if(neighboursOfV[neighOfVOffset+sP] == 0)
			{
			  --numNeighParts[locVertIndex];				  
			}
		    }

		  vertSeen.set1(locVertIndex);
		  seenVertices[numVerticesSeen++] = locVertIndex;		  
		}
	    }

	  // ###
	  // if the adjacent vertex is not a local vertex
	  // we don't update anything
	  // ###
	} 
    }
	
  // ###
  // restore the 'seen' vertices structure
  // ###

  for (i=0;i<numVerticesSeen;++i)		    
    vertSeen.set0(seenVertices[i]);
}	      




void ParaGreedyKwayRefiner::unmakeMoves(int indexIntoMoveSets, int from, int to)
{
#  ifdef DEBUG_REFINER
  sanityHedgeCheck();
#  endif

  register int v;  
  register int j;

  int i;
  int hEdgeOff;
  int neighOfVOffset;
  int vertOffset;
  int hEdge;

  int numVmoved = numVerticesMoved[indexIntoMoveSets];
  int *movedArray = moveSets[indexIntoMoveSets]->getArray();  

  for (i=0;i<numVmoved;++i) 
    {    
      v = movedArray[i]-minVertexIndex;
      
#  ifdef DEBUG_REFINER
      assert(v >= 0 && v < numLocalVertices);
#  endif
      
      // ###
      // make the move and update the vertex's structs
      // ###

      vertOffset = vToHedgesOffset[v+1];
      neighOfVOffset = neighboursOfVOffsets[v];
      
      // ###
      // update the moved vertices' stats
      // ###
      
      if(neighboursOfV[neighOfVOffset+from] == 0)
	{
	  ++numNeighParts[v];
	}
		  
      neighboursOfV[neighOfVOffset+from] = 1;
      neighboursOfV[neighOfVOffset+to] = 0;	      	     
		  
      for (j=vToHedgesOffset[v];j<vertOffset;++j)
	{
	  // ###
	  // update the hyperedge stats: (vInPart etc.)
	  // ###
		      
	  hEdge = vToHedgesList[j];
#  ifdef DEBUG_REFINER
	  assert(hEdge >= 0 && hEdge < numHedges);
#  endif
	  hEdgeOff = hEdgeVinPartOffsets[hEdge];

#  ifdef DEBUG_REFINER
	  int hEdgeLen = 0;
	  for (int ijk=0;ijk<numParts;++ijk)
	    hEdgeLen += hEdgeVinPart[hEdgeOff+ijk];
	  assert(hEdgeLen == hEdgeOffset[hEdge+1]-hEdgeOffset[hEdge]);
#  endif		      
	  if(hEdgeVinPart[hEdgeOff+from] == 0)
	    {
	      ++numPartsSpanned[hEdge];
	    }

	  --hEdgeVinPart[hEdgeOff+to];		
	  ++hEdgeVinPart[hEdgeOff+from];
		      
#  ifdef DEBUG_REFINER
	  assert(hEdgeVinPart[hEdgeOff+to] >= 0);
#  endif
		      
	  if(hEdgeVinPart[hEdgeOff+to] > 0)
	    {
	      neighboursOfV[neighOfVOffset+to] = 1;
	    }
	  else
	    {
	      --numPartsSpanned[hEdge];
	    }		  
	}
		  
      if(neighboursOfV[neighOfVOffset+to] == 0)
	{
	  --numNeighParts[v];
	}
		  
      // ###
      // update the adj vertices stats: 
      // (num neighbours in part etc.)
      // ###
      
      updateAdjVertStatus(v,to,from);	      
		  
      // ###
      // update other structs
      // ###
      
      locked.set0(v);
      currPVector[v] = from;
    }

#  ifdef DEBUG_REFINER
  sanityHedgeCheck();
#  endif
  
  numVerticesMoved[indexIntoMoveSets] = 0;
}




void ParaGreedyKwayRefiner::nonLocVertCheck() const
{
  register int i;
  register int j;
  register int ij;

  int vertPart;
  int endOffset;
  int h;

  for (i=0;i<numNonLocVerts;++i)
    {
      vertPart = currNonLocPVector[i];
      endOffset = nonLocOffsets[i+1];
      
      for (j=nonLocOffsets[i];j<endOffset;++j)
	{
	  assert(j < nonLocVToHedges.getLength());
	  h = nonLocVToHedges[j];
	  assert(h < hEdgeVinPartOffsets.getLength());
	  ij = hEdgeVinPartOffsets[h];
	  assert(hEdgeVinPart[ij+vertPart] > 0);
	}
    }
}




void ParaGreedyKwayRefiner::sanityHedgeCheck() const
{
  register int i;
  register int j;
  register int ij;

  int hEdgeLen;
  int inParts;

  for (i=0;i<numHedges;++i)
    {
      hEdgeLen = hEdgeOffset[i+1]-hEdgeOffset[i]; 
      inParts = 0;
      ij = hEdgeVinPartOffsets[i+1];
      
      for (j=hEdgeVinPartOffsets[i];j<ij;++j)
	{
	  inParts += hEdgeVinPart[j];
	  assert(hEdgeVinPart[j] >= 0);
	}

      assert(inParts == hEdgeLen);
    }
}


#  endif
