
#ifndef _VCYCLEFINAL_CONTROLLER_CPP
#define _VCYCLEFINAL_CONTROLLER_CPP

// ### ParaVCycleFinalController.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 10/1/2005: Last Modified
//
// ###

#include "ParaVCycleFinalController.hpp"

ParaVCycleFinalController::ParaVCycleFinalController(
    ParaRestrCoarsener &rc, ParaCoarsener &c, ParaRefiner &r,
    SeqController &ref, int rank, int nP, int percentile, int inc,
    int approxRef, int limit, double limitAsPercent, ostream &out)
    : ParaVCycleController(rc, c, r, ref, rank, nP, percentile, inc, approxRef,
                           limit, limitAsPercent, out) {}

ParaVCycleFinalController::~ParaVCycleFinalController() {}

void ParaVCycleFinalController::runPartitioner(MPI_Comm comm) {
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

  stack<int> hEdgePercentiles;

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

#ifdef DEBUG_CONTROLLER
  int checkCutsize;
#endif

  for (i = 0; i < numParaRuns; ++i) {
    if (shuffled == 1)
      hgraph->randomVertexShuffle(mapToOrigVerts.data(), comm);

    if (shuffled == 2)
      hgraph->prescribedVertexShuffle(mapToOrigVerts.data(),
                                      shufflePartition.data(), comm);

    hgraphs.push(hgraph);
    hEdgePercentiles.push(startPercentile);

    finerGraph = hgraph;
    accumulator = 1.0;

    // ###
    // coarsen the hypergraph
    // ###

    MPI_Barrier(comm);
    startTime = MPI_Wtime();

    do {
      hEdgePercentile = hEdgePercentiles.top();
      coarsener.setPercentile(hEdgePercentile);
      coarseGraph = coarsener.coarsen(*finerGraph, comm);

      if (coarseGraph) {
        hEdgePercentiles.push(min(hEdgePercentile + percentileIncrement, 100));
        hgraphs.push(coarseGraph);
        finerGraph = coarseGraph;
      }
    }

    while (coarseGraph);

    MPI_Barrier(comm);
    totCoaTime += (MPI_Wtime() - startTime);

    coarsener.releaseMemory();

    // ###
    // compute the initial partition
    // ###

    MPI_Barrier(comm);
    startTime = MPI_Wtime();

    coarseGraph = hgraphs.pop();
    seqController.runSeqPartitioner(*coarseGraph, comm);

    MPI_Barrier(comm);
    totSeqTime += (MPI_Wtime() - startTime);

    // ###
    // uncoarsen the initial partition
    // ###

    MPI_Barrier(comm);
    startTime = MPI_Wtime();

    while (hgraphs.size() > 0) {
      coarseGraph->removeBadPartitions(keepPartitionsWithin * accumulator);
      accumulator *= reductionInKeepThreshold;

      hEdgePercentile = hEdgePercentiles.pop();
      finerGraph = hgraphs.pop();

      if (finerGraph == hgraph)
        projectVCyclePartition(*coarseGraph, *finerGraph, comm);
      else
        finerGraph->projectPartitions(*coarseGraph, comm);

#ifdef DEBUG_CONTROLLER
      finerGraph->checkPartitions(numTotalParts, maxPartWt, comm);
#endif
      if (approxRefine)
        refiner.setPercentile(hEdgePercentile);

      refiner.refine(*finerGraph, comm);

#ifdef DEBUG_CONTROLLER
      finerGraph->checkPartitions(numTotalParts, maxPartWt, comm);
#endif
      DynaMem::deletePtr<ParaHypergraph>(coarseGraph);

      coarseGraph = finerGraph;
    }

#ifdef DEBUG_CONTROLLER
    assert(coarseGraph == hgraph);
#endif

    MPI_Barrier(comm);
    totRefTime += (MPI_Wtime() - startTime);

    refiner.releaseMemory();

    // ###
    // select the best partition
    // ###

    firstCutSize = coarseGraph->keepBestPartition();

#ifdef DEBUG_CONTROLLER
    checkCutsize = coarseGraph->calcCutsize(numTotalParts, 0, comm);
    assert(firstCutSize == checkCutsize);
#endif

    numIterations = 0;
    vCycleGain = 0;
    recordVCyclePartition(*hgraph, numIterations++);

    if (dispOption > 1 && myRank == 0) {
      out_stream << "\t ------ PARALLEL V-CYCLE CALL ------" << endl;
    }

    do {
      minIterationGain =
          static_cast<int>(floor(limAsPercentOfCut * firstCutSize));

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

      do {
        hEdgePercentile = hEdgePercentiles.top();
        restrCoarsener.setPercentile(hEdgePercentile);
        coarseGraph = restrCoarsener.coarsen(*finerGraph, comm);

        if (coarseGraph) {
          hEdgePercentiles.push(
              min(hEdgePercentile + percentileIncrement, 100));
          hgraphs.push(coarseGraph);
          finerGraph = coarseGraph;
        }
      }

      while (coarseGraph);

      MPI_Barrier(comm);
      totCoaTime += (MPI_Wtime() - startTime);

      restrCoarsener.releaseMemory();

      // ###
      // compute the initial partition
      // ###

      coarseGraph = hgraphs.pop();
      coarseGraph->setNumberPartitions(0);
      coarseGraph->shiftVerticesToBalance(comm);

      MPI_Barrier(comm);
      startTime = MPI_Wtime();

      seqController.runSeqPartitioner(*coarseGraph, comm);

      MPI_Barrier(comm);
      totSeqTime += (MPI_Wtime() - startTime);

      // ###
      // uncoarsen the initial partition
      // ###

      MPI_Barrier(comm);
      startTime = MPI_Wtime();
      accumulator = 1.0;

      while (hgraphs.size() > 0) {
        coarseGraph->removeBadPartitions(keepPartitionsWithin * accumulator);
        accumulator *= reductionInKeepThreshold;

        hEdgePercentile = hEdgePercentiles.pop();
        finerGraph = hgraphs.pop();

        if (finerGraph == hgraph)
          shiftVCycleVertsToBalance(*finerGraph, comm);
        else
          finerGraph->shiftVerticesToBalance(comm);

        finerGraph->projectPartitions(*coarseGraph, comm);

#ifdef DEBUG_CONTROLLER
        finerGraph->checkPartitions(numTotalParts, maxPartWt, comm);
#endif
        if (randShuffBefRef) {
          if (hgraphs.size() == 0)
            finerGraph->randomVertexShuffle(mapToInterVerts.data(), comm);
          else
            finerGraph->randomVertexShuffle(*(hgraphs.top()), comm);
        }

        if (approxRefine)
          refiner.setPercentile(hEdgePercentile);

        refiner.refine(*finerGraph, comm);

#ifdef DEBUG_CONTROLLER
        finerGraph->checkPartitions(numTotalParts, maxPartWt, comm);
#endif
        DynaMem::deletePtr<ParaHypergraph>(coarseGraph);

        coarseGraph = finerGraph;
      }

#ifdef DEBUG_CONTROLLER
      assert(coarseGraph == hgraph);
      checkCutsize = coarseGraph->calcCutsize(numTotalParts, 0, comm);
#endif
      MPI_Barrier(comm);
      totRefTime += (MPI_Wtime() - startTime);

      refiner.releaseMemory();

      // ###
      // select the best partition
      // ###

      secondCutSize = coarseGraph->keepBestPartition();
      diffInCutSize = firstCutSize - secondCutSize;

      if (dispOption > 1 && myRank == 0) {
        out_stream << "\t ------ [" << numIterations << "] " << diffInCutSize
                   << endl;
      }

#ifdef DEBUG_CONTROLLER
      checkCutsize = coarseGraph->calcCutsize(numTotalParts, 0, comm);
      assert(secondCutSize == checkCutsize);
#endif

      if (diffInCutSize > 0 && diffInCutSize < minIterationGain)
        break;

      if (numIterations == limitOnCycles)
        break;

      if (firstCutSize > secondCutSize) {
        recordVCyclePartition(*coarseGraph, numIterations++);
        firstCutSize = secondCutSize;
        vCycleGain += diffInCutSize;
      }
    } while (diffInCutSize > 0);

    gatherInVCyclePartition(*hgraph, firstCutSize, comm);
    updateMapToOrigVerts(comm);

    if (dispOption > 1 && myRank == 0) {
      out_stream << "\t ------ " << vCycleGain << " ------" << endl;
    }

#ifdef DEBUG_CONTROLLER
    assert(hgraphs.getNumElem() == 0);
    assert(hgraph->getNumLocalVertices() == numOrigLocVerts);
    hgraph->checkPartitions(numTotalParts, maxPartWt, comm);
#endif

    if (myRank == 0 && dispOption > 0) {
      out_stream << "\nPRUN[" << i << "] = " << firstCutSize << endl << endl;
    }

    totCutsizes += firstCutSize;

    if (firstCutSize < bestCutsize) {
      bestCutsize = firstCutSize;
      storeBestPartition(numOrigLocVerts, hgraph->getPartVectorArray(), comm);
    }

    if (firstCutSize > worstCutsize)
      worstCutsize = firstCutSize;

    // ###
    // reset hypergraph
    // ###

    resetStructs();
  }

  MPI_Barrier(comm);
  totalTime = MPI_Wtime() - totStartTime;

  percentCoarsening = static_cast<int>(floor((totCoaTime / totalTime) * 100));
  percentSequential = static_cast<int>(floor((totSeqTime / totalTime) * 100));
  percentRefinement = static_cast<int>(floor((totRefTime / totalTime) * 100));
  percentOther =
      100 - (percentCoarsening + percentSequential + percentRefinement);

  aveCutSize = static_cast<double>(totCutsizes) / numParaRuns;

  if (myRank == 0 && dispOption > 0) {
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
               << "|-- AVE TIME = " << totalTime / numParaRuns << endl
               << "|-- PARACOARSENING% = " << percentCoarsening << endl
               << "|-- SEQPARTITIONING% = " << percentSequential << endl
               << "|-- PARAREFINEMENT% = " << percentRefinement << endl
               << "|-- OTHER% = " << percentOther << endl
               << "|" << endl
               << " ----------------------------" << endl;
  }
}

void ParaVCycleFinalController::printVCycleType() const {
  out_stream << " type = BIG";
}

#endif
