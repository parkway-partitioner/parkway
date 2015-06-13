#ifndef _BASIC_PARA_CONTROLLER_CPP
#define _BASIC_PARA_CONTROLLER_CPP

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

#include "BasicParaController.hpp"
#include "hypergraph/parallel/hypergraph.hpp"

using parkway::hypergraph::parallel::hypergraph;

BasicParaController::BasicParaController(parallel_coarsener &c, ParaRefiner &r,
                                         SeqController &ref, int rank, int nP,
                                         int percentile, int inc, int approxRef,
                                         ostream &out)
    : ParaController(c, r, ref, rank, nP, percentile, inc, approxRef, out) {}

BasicParaController::~BasicParaController() {}

void BasicParaController::dispParaControllerOptions() const {
  switch (dispOption) {
  case SILENT:
    break;

  default:

    out_stream << "|--- PARA_CONTR (# parts = " << numTotalParts
               << "): " << endl
               << "|- BASIC:"
               << " pRuns = " << numParaRuns << " kT = " << keepPartitionsWithin
               << " rKT = " << reductionInKeepThreshold
               << " appRef = " << approxRefine
               << " wTF = " << writePartitionToFile << " start %le "
               << startPercentile << " %le inc " << percentileIncrement << endl
               << "|" << endl;
    break;
  }
}

void BasicParaController::runPartitioner(MPI_Comm comm) {
  int hEdgePercentile;
  int cutSize;
  int percentCoarsening;
  int percentSequential;
  int percentRefinement;
  int percentOther;
  int i;

  double totStartTime;

#ifdef DEBUG_CONTROLLER
  int checkCutsize;
#endif

  stack<int> hEdgePercentiles;

  numOrigLocVerts = hgraph->number_of_vertices();

  hypergraph *coarseGraph;
  hypergraph *finerGraph;

  initMapToOrigVerts();

  bestCutsize = LARGE_CONSTANT;
  worstCutsize = 0;
  totCutsizes = 0;

  totCoaTime = 0;
  totSeqTime = 0;
  totRefTime = 0;

  MPI_Barrier(comm);
  totStartTime = MPI_Wtime();

  for (i = 0; i < numParaRuns; ++i) {
#ifdef MEM_CHECK
    write_log(rank_, "[begin run]: usage: %f", MemoryTracker::usage());
    Funct::printMemUse(rank_, "[begin run]");
#endif
    if (shuffled == 1) {
        hgraph->shuffle_vertices_randomly(mapToOrigVerts.data(), comm);
    }

    if (shuffled == 2) {
        hgraph->prescribed_vertex_shuffle(mapToOrigVerts.data(),
                                          shufflePartition.data(), comm);
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

    do {
      hEdgePercentile = hEdgePercentiles.top();
      coarsener.set_percentile(hEdgePercentile);

      coarseGraph = coarsener.coarsen(*finerGraph, comm);
        finerGraph->free_memory();

      if (coarseGraph) {
        hEdgePercentiles.push(min(hEdgePercentile + percentileIncrement, 100));
        hgraphs.push(coarseGraph);
        finerGraph = coarseGraph;
      }
    }

    while (coarseGraph);

    MPI_Barrier(comm);
    totCoaTime += (MPI_Wtime() - startTime);

    coarsener.release_memory();

#ifdef MEM_CHECK
    MPI_Barrier(comm);
    write_log(rank_, "[after coarsening]: usage: %f", MemoryTracker::usage());
    Funct::printMemUse(rank_, "[after coarsening]");
    MPI_Barrier(comm);
#endif
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

#ifdef MEM_CHECK
    MPI_Barrier(comm);
    write_log(rank_, "[after seq partitioning]: usage: %f",
              MemoryTracker::usage());
    Funct::printMemUse(rank_, "[after seq partitioning]");
    MPI_Barrier(comm);
#endif

    MPI_Barrier(comm);
    startTime = MPI_Wtime();

    while (hgraphs.size() > 0) {
        coarseGraph->remove_bad_partitions(keepPartitionsWithin * accumulator);
      accumulator *= reductionInKeepThreshold;

      finerGraph = hgraphs.pop();
      hEdgePercentile = hEdgePercentiles.pop();

        finerGraph->project_partitions(*coarseGraph, comm);

#ifdef DEBUG_CONTROLLER
      finerGraph->checkPartitions(numTotalParts, maxPartWt, comm);
#endif
      if (randShuffBefRef) {
        if (hgraphs.size() == 0)
            finerGraph->shuffle_vertices_randomly(mapToOrigVerts.data(), comm);
        else
            finerGraph->shuffle_vertices_randomly(*(hgraphs.top()), comm);
      }

      if (approxRefine)
        refiner.set_percentile(hEdgePercentile);

#ifdef MEM_CHECK
      write_log(rank_, "[before refineme]: usage: %f", MemoryTracker::usage());
      Funct::printMemUse(rank_, "[before refinement]");
#endif
      refiner.refine(*finerGraph, comm);

#ifdef DEBUG_CONTROLLER
      finerGraph->checkPartitions(numTotalParts, maxPartWt, comm);
#endif

      DynaMem::deletePtr<hypergraph>(coarseGraph);
      coarseGraph = finerGraph;
    }

#ifdef DEBUG_CONTROLLER
    assert(coarseGraph == hgraph);
#endif

    MPI_Barrier(comm);
    totRefTime += (MPI_Wtime() - startTime);

    refiner.release_memory();

#ifdef MEM_CHECK
    MPI_Barrier(comm);
    write_log(rank_, "[after refinement]: usage: %f", MemoryTracker::usage());
    Funct::printMemUse(rank_, "[after refinement]");
    MPI_Barrier(comm);
#endif

    /* select the best partition */

    cutSize = coarseGraph->keep_best_partition();

#ifdef DEBUG_CONTROLLER
    checkCutsize = coarseGraph->calcCutsize(numTotalParts, 0, comm);
    assert(cutSize == checkCutsize);
#endif

    if (cutSize < bestCutsize) {
      bestCutsize = cutSize;
      storeBestPartition(numOrigLocVerts, hgraph->partition_vector(), comm);
    }

    if (cutSize > worstCutsize)
      worstCutsize = cutSize;

    totCutsizes += cutSize;

    /* free memory used by the para controller */

    resetStructs();

#ifdef DEBUG_CONTROLLER
    assert(hgraphs.getNumElem() == 0);
#endif

    if (rank_ == 0 && dispOption > 0) {
      out_stream << "\nPRUN[" << i << "] = " << cutSize << endl << endl;
    }
  }

  MPI_Barrier(comm);
  totalTime = MPI_Wtime() - totStartTime;

  percentCoarsening = static_cast<int>(floor((totCoaTime / totalTime) * 100));
  percentSequential = static_cast<int>(floor((totSeqTime / totalTime) * 100));
  percentRefinement = static_cast<int>(floor((totRefTime / totalTime) * 100));
  percentOther =
      100 - (percentCoarsening + percentSequential + percentRefinement);

  aveCutSize = static_cast<double>(totCutsizes) / numParaRuns;

  if (rank_ == 0 && dispOption > 0) {
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

void BasicParaController::resetStructs() {
    hgraph->reset_vectors();
    free_memory();
}

#endif
