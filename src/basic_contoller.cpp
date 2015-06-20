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

#include "basic_contoller.hpp"
#include "hypergraph/parallel/hypergraph.hpp"

namespace parkway {
namespace parallel {

basic_contoller::basic_contoller(coarsener &c, refiner &r,
                                 serial::controller &ref, int rank, int nP, int
                                 percentile, int inc, int approxRef,
                                 std::ostream &out)
    : controller(c, r, ref, rank, nP, percentile, inc, approxRef, out) {
}

basic_contoller::~basic_contoller() {}

void basic_contoller::display_options() const {
  switch (display_option_) {
  case SILENT:
    break;

  default:

    out_stream_ << "|--- PARA_CONTR (# parts = " << total_number_of_parts_
               << "): " << std::endl
               << "|- BASIC:"
               << " pRuns = " << number_of_runs_ << " kT = " <<
                                                    keep_partitions_within_
               << " rKT = " << reduction_in_keep_threshold_
               << " appRef = " << approximate_refine_
               << " wTF = " << write_partition_to_file_ << " start %le "
               << start_percentile_ << " %le inc " << percentile_increment_ << std::endl
               << "|" << std::endl;
    break;
  }
}

void basic_contoller::run(MPI_Comm comm) {
  int hEdgePercentile;
  int cutSize;
  int percentCoarsening;
  int percentserial;
  int percentRefinement;
  int percentOther;
  int i;

  double totStartTime;

#ifdef DEBUG_CONTROLLER
  int checkCutsize;
#endif

  ds::stack<int> hEdgePercentiles;

  number_of_orig_local_vertices_ = hypergraph_->number_of_vertices();

  parallel::hypergraph *coarseGraph;
  parallel::hypergraph *finerGraph;

  initialize_map_to_orig_verts();

  best_cutsize_ = LARGE_CONSTANT;
  worst_cutsize_ = 0;
  total_cutsize_ = 0;

  total_coarsening_time_ = 0;
  total_serial_time_ = 0;
  total_refinement_time_ = 0;

  MPI_Barrier(comm);
  totStartTime = MPI_Wtime();

  for (i = 0; i < number_of_runs_; ++i) {
#ifdef MEM_CHECK
    write_log(rank_, "[begin run]: usage: %f", MemoryTracker::usage());
    Funct::printMemUse(rank_, "[begin run]");
#endif
    if (shuffled_ == 1) {
        hypergraph_->shuffle_vertices_randomly(map_to_orig_vertices_, comm);
    }

    if (shuffled_ == 2) {
        hypergraph_->prescribed_vertex_shuffle(map_to_orig_vertices_,
                                               shuffle_partition_, comm);
    }

    hypergraphs_.push(hypergraph_);
    hEdgePercentiles.push(start_percentile_);

    finerGraph = hypergraph_;
    accumulator_ = 1.0;

    // ###
    // coarsen the hypergraph
    // ###

    MPI_Barrier(comm);
    start_time_ = MPI_Wtime();

    do {
      hEdgePercentile = hEdgePercentiles.top();
      coarsener_.set_percentile(hEdgePercentile);

      coarseGraph = coarsener_.coarsen(*finerGraph, comm);

      // returns hgraph with zero offsets
      // offsets calculated on next coarsen
      finerGraph->free_memory();

      if (coarseGraph) {
        hEdgePercentiles.push(std::min(hEdgePercentile + percentile_increment_, 100));
        hypergraphs_.push(coarseGraph);
        finerGraph = coarseGraph;
      }
    } while (coarseGraph);

    MPI_Barrier(comm);
    total_coarsening_time_ += (MPI_Wtime() - start_time_);

    coarsener_.release_memory();

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
    start_time_ = MPI_Wtime();

    coarseGraph = hypergraphs_.pop();
    serial_controller_.run(*coarseGraph, comm);

    MPI_Barrier(comm);
    total_serial_time_ += (MPI_Wtime() - start_time_);

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
    start_time_ = MPI_Wtime();

    while (hypergraphs_.size() > 0) {
        coarseGraph->remove_bad_partitions(keep_partitions_within_ *
                                           accumulator_);
      accumulator_ *= reduction_in_keep_threshold_;

      finerGraph = hypergraphs_.pop();
      hEdgePercentile = hEdgePercentiles.pop();

        finerGraph->project_partitions(*coarseGraph, comm);

      if (random_shuffle_before_refine_) {
        if (hypergraphs_.size() == 0)
            finerGraph->shuffle_vertices_randomly(map_to_orig_vertices_, comm);
        else
            finerGraph->shuffle_vertices_randomly(*(hypergraphs_.top()), comm);
      }

      if (approximate_refine_)
        refiner_.set_percentile(hEdgePercentile);

      refiner_.refine(*finerGraph, comm);
      coarseGraph = finerGraph;
    }

    MPI_Barrier(comm);
    total_refinement_time_ += (MPI_Wtime() - start_time_);

    refiner_.release_memory();

    /* select the best partition */
    cutSize = coarseGraph->keep_best_partition();

    if (cutSize < best_cutsize_) {
      best_cutsize_ = cutSize;
      store_best_partition(number_of_orig_local_vertices_,
                           hypergraph_->partition_vector(), comm);
    }

    if (cutSize > worst_cutsize_)
      worst_cutsize_ = cutSize;

    total_cutsize_ += cutSize;

    /* free memory used by the para controller */
    reset_structures();
    if (rank_ == 0 && display_option_ > 0) {
      out_stream_ << "\nPRUN[" << i << "] = " << cutSize << std::endl << std::endl;
    }
  }

  MPI_Barrier(comm);
  total_time_ = MPI_Wtime() - totStartTime;

  percentCoarsening = static_cast<int>(floor((total_coarsening_time_ /
                                              total_time_) * 100));
  percentserial = static_cast<int>(floor((total_serial_time_ /
                                              total_time_) * 100));
  percentRefinement = static_cast<int>(floor((total_refinement_time_ /
                                              total_time_) * 100));
  percentOther =
      100 - (percentCoarsening + percentserial + percentRefinement);

  average_cutsize_ = static_cast<double>(total_cutsize_) / number_of_runs_;

  if (rank_ == 0 && display_option_ > 0) {
    out_stream_ << std::endl
               << " --- PARTITIONING SUMMARY ---" << std::endl
               << "|" << std::endl
               << "|--- Cutsizes statistics:" << std::endl
               << "|" << std::endl
               << "|-- BEST = " << best_cutsize_ << std::endl
               << "|-- WORST = " << worst_cutsize_ << std::endl
               << "|-- AVE = " << average_cutsize_ << std::endl
               << "|" << std::endl
               << "|--- Time usage:" << std::endl
               << "|" << std::endl
               << "|-- TOTAL TIME = " << total_time_ << std::endl
               << "|-- AVE TIME = " << total_time_ / number_of_runs_ << std::endl
               << "|-- PARACOARSENING% = " << percentCoarsening << std::endl
               << "|-- SEQPARTITIONING% = " << percentserial << std::endl
               << "|-- PARAREFINEMENT% = " << percentRefinement << std::endl
               << "|-- OTHER% = " << percentOther << std::endl
               << "|" << std::endl
               << " ----------------------------" << std::endl;
  }
}

void basic_contoller::reset_structures() {
    hypergraph_->reset_vectors();
    free_memory();
}

}  // namespace parallel
}  // namespace parkway
