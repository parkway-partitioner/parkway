#ifndef _VCYCLEALL_CONTROLLER_CPP
#define _VCYCLEALL_CONTROLLER_CPP

// ### ParaVCycleAllController.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 4/1/2005: Last Modified
//
// ###

#include "parallel_v_cycle_all_controller.hpp"

parallel_v_cycle_all_controller::parallel_v_cycle_all_controller(
    parallel_restrictive_coarsening &rc, parallel_coarsener &c, parallel_refiner &r,
    sequential_controller &ref, int rank, int nP, int percentile, int inc,
    int approxRef, int limit, double limitAsPercent, ostream &out)
    : parallel_v_cycle_controller(rc, c, r, ref, rank, nP, percentile, inc, approxRef,
                           limit, limitAsPercent, out) {}

parallel_v_cycle_all_controller::~parallel_v_cycle_all_controller() {}

void parallel_v_cycle_all_controller::run(MPI_Comm comm) {
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
  int numOrigTotVerts = hypergraph_->total_number_of_vertices();

  double othAccumulator;
  double totStartTime;

  stack<int> hEdgePercentiles;

  parallel::hypergraph *coarseGraph;
  parallel::hypergraph *interMedGraph;
  parallel::hypergraph *finerGraph;

  initialize_map_to_orig_verts();

  best_cutsize_ = LARGE_CONSTANT;
  worst_cutsize_ = 0;
  total_cutsize_ = 0;

#ifdef DEBUG_CONTROLLER
  int checkCutsize;
#endif

  MPI_Barrier(comm);
  totStartTime = MPI_Wtime();

  for (i = 0; i < number_of_runs_; ++i) {
    if (shuffled_ == 1)
        hypergraph_->shuffle_vertices_randomly(map_to_orig_vertices_.data(), comm);

    if (shuffled_ == 2)
        hypergraph_->prescribed_vertex_shuffle(map_to_orig_vertices_.data(),
                                          shuffle_partition_.data(), comm);

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

      if (coarseGraph) {
        hEdgePercentiles.push(min(hEdgePercentile + percentile_increment_, 100));
        hypergraphs_.push(coarseGraph);
        finerGraph = coarseGraph;
      }
    }

    while (coarseGraph);

    MPI_Barrier(comm);
    total_coarsening_time_ += (MPI_Wtime() - start_time_);

    coarsener_.release_memory();

    // ###
    // compute the initial partition
    // ###

    MPI_Barrier(comm);
    start_time_ = MPI_Wtime();

    coarseGraph = hypergraphs_.pop();
    sequential_controller_.run(*coarseGraph, comm);

    MPI_Barrier(comm);
    total_sequential_time_ += (MPI_Wtime() - start_time_);

    // ###
    // uncoarsen the initial partition
    // ###

    while (hypergraphs_.size() > 0) {
        coarseGraph->remove_bad_partitions(keep_partitions_within_ *
                                           accumulator_);
      accumulator_ *= reduction_in_keep_threshold_;
      hEdgePercentile = hEdgePercentiles.pop();
      finerGraph = hypergraphs_.pop();

      MPI_Barrier(comm);
      start_time_ = MPI_Wtime();

      project_v_cycle_partition(*coarseGraph, *finerGraph, comm);

#ifdef DEBUG_CONTROLLER
      finerGraph->checkPartitions(numTotalParts, maxPartWt, comm);
#endif
      if (approximate_refine_)
        refiner_.set_percentile(hEdgePercentile);

      refiner_.refine(*finerGraph, comm);
      refiner_.release_memory();

#ifdef DEBUG_CONTROLLER
      finerGraph->checkPartitions(numTotalParts, maxPartWt, comm);
#endif

      MPI_Barrier(comm);
      total_refinement_time_ += (MPI_Wtime() - start_time_);

      dynamic_memory::delete_pointer<parallel::hypergraph>(coarseGraph);

      /* choose/reject v-cycle */

      if (finerGraph->total_number_of_vertices() > numOrigTotVerts / 4)
        doVCycle = 0;
      else
        doVCycle = 1;

      if (doVCycle) {
        // doing v-cycle

        vCycleIteration = 0;
        vCycleGain = 0;
        firstCutSize = finerGraph->keep_best_partition();
        record_v_cycle_partition(*finerGraph, vCycleIteration++);

        numInStack = hypergraphs_.size();
        interMedGraph = finerGraph;

        if (display_option_ > 1 && rank_ == 0) {
          out_stream_ << "\t ------ PARALLEL V-CYCLE CALL ------" << endl;
        }

        do {
          minIterationGain =
              static_cast<int>(floor(limit_as_percent_of_cut_ * firstCutSize));

          shuffle_v_cycle_vertices_by_partition(*finerGraph, comm);
          hypergraphs_.push(finerGraph);

          othAccumulator = 1.0;

          // ###
          // coarsen the hypergraph
          // ###

          MPI_Barrier(comm);
          start_time_ = MPI_Wtime();

          do {
            hEdgePercentile = hEdgePercentiles.top();
            restrictive_coarsening_.set_percentile(hEdgePercentile);
            coarseGraph = restrictive_coarsening_.coarsen(*finerGraph, comm);

            if (coarseGraph) {
              hEdgePercentiles.push(
                  min(hEdgePercentile + percentile_increment_, 100));
              hypergraphs_.push(coarseGraph);
              finerGraph = coarseGraph;
            }
          } while (coarseGraph);

          MPI_Barrier(comm);
          total_coarsening_time_ += (MPI_Wtime() - start_time_);

          restrictive_coarsening_.release_memory();

          // ###
          // compute the initial partition
          // ###

          coarseGraph = hypergraphs_.pop();
            coarseGraph->set_number_of_partitions(0);
            coarseGraph->shift_vertices_to_balance(comm);

          MPI_Barrier(comm);
          start_time_ = MPI_Wtime();

          sequential_controller_.run(*coarseGraph, comm);

          MPI_Barrier(comm);
          total_sequential_time_ += (MPI_Wtime() - start_time_);

          // ###
          // uncoarsen the initial partition
          // ###

          MPI_Barrier(comm);
          start_time_ = MPI_Wtime();

          while (hypergraphs_.size() > numInStack) {
              coarseGraph->remove_bad_partitions(keep_partitions_within_ *
                                                 othAccumulator);
            othAccumulator *= reduction_in_keep_threshold_;

            hEdgePercentile = hEdgePercentiles.pop();
            finerGraph = hypergraphs_.pop();

            if (finerGraph == interMedGraph) {
              shift_v_cycle_vertices_to_balance(*finerGraph, comm);
            } else {
                finerGraph->shift_vertices_to_balance(comm);
            }

              finerGraph->project_partitions(*coarseGraph, comm);

#ifdef DEBUG_CONTROLLER
            finerGraph->checkPartitions(numTotalParts, maxPartWt, comm);
#endif

            if (random_shuffle_before_refine_) {
              if (finerGraph == interMedGraph)
                  finerGraph->shuffle_vertices_randomly(map_to_inter_vertices_.data(),
                                                        comm);
              else
                  finerGraph->shuffle_vertices_randomly(*(hypergraphs_.top()), comm);
            }

            if (approximate_refine_)
              refiner_.set_percentile(hEdgePercentile);

            refiner_.refine(*finerGraph, comm);
#ifdef DEBUG_CONTROLLER
            finerGraph->checkPartitions(numTotalParts, maxPartWt, comm);
#endif
            dynamic_memory::delete_pointer<parallel::hypergraph>(coarseGraph);

            coarseGraph = finerGraph;
          }

          MPI_Barrier(comm);
          total_refinement_time_ += (MPI_Wtime() - start_time_);

#ifdef DEBUG_CONTROLLER
          assert(coarseGraph == interMedGraph);
#endif
          refiner_.release_memory();

          // ###
          // select the best partition
          // ###

          secondCutSize = coarseGraph->keep_best_partition();

#ifdef DEBUG_CONTROLLER
          coarseGraph->checkPartitions(numTotalParts, maxPartWt, comm);
#endif
          diffInCutSize = firstCutSize - secondCutSize;

          if (display_option_ > 1 && rank_ == 0) {
            out_stream_ << "\t ------ [" << vCycleIteration << "] "
                       << diffInCutSize << endl;
          }

          if (diffInCutSize > 0 && diffInCutSize < minIterationGain)
            break;

          if (vCycleIteration == limit_on_cycles_)
            break;

          if (firstCutSize > secondCutSize) {
            record_v_cycle_partition(*coarseGraph, vCycleIteration++);
            firstCutSize = secondCutSize;
            vCycleGain += diffInCutSize;
          }
        } while (diffInCutSize > 0);

#ifdef DEBUG_CONTROLLER
        assert(coarseGraph == interMedGraph);
#endif

        gather_in_v_cycle_partition(*coarseGraph, firstCutSize, comm);

        if (display_option_ > 1 && rank_ == 0) {
          out_stream_ << "\t ------ " << vCycleGain << " ------" << endl;
        }
      } else {
        // not doing v-cycle
        coarseGraph = finerGraph;
      }

#ifdef DEBUG_CONTROLLER
      coarseGraph->checkPartitions(numTotalParts, maxPartWt, comm);
#endif
    }

#ifdef DEBUG_CONTROLLER
    assert(coarseGraph == hgraph);
#endif

    firstCutSize = coarseGraph->keep_best_partition();

    update_map_to_orig_vertices(comm);

#ifdef DEBUG_CONTROLLER
    checkCutsize = coarseGraph->calcCutsize(numTotalParts, 0, comm);
    assert(firstCutSize == checkCutsize);
#endif

    if (rank_ == 0 && display_option_ > 0) {
      out_stream_ << "\nPRUN[" << i << "] " << firstCutSize << endl << endl;
    }

    total_cutsize_ += firstCutSize;

    if (firstCutSize < best_cutsize_) {
      best_cutsize_ = firstCutSize;
      store_best_partition(number_of_orig_local_vertices_,
                           hypergraph_->partition_vector(), comm);
    }

    if (firstCutSize > worst_cutsize_)
      worst_cutsize_ = firstCutSize;

    // ###
    // reset hypergraph
    // ###

    reset_structures();
  }

  MPI_Barrier(comm);
  total_time_ = MPI_Wtime() - totStartTime;

  percentCoarsening = static_cast<int>(floor((total_coarsening_time_ /
                                              total_time_) * 100));
  percentSequential = static_cast<int>(floor((total_sequential_time_ /
                                              total_time_) * 100));
  percentRefinement = static_cast<int>(floor((total_refinement_time_ /
                                              total_time_) * 100));
  percentOther =
      100 - (percentCoarsening + percentSequential + percentRefinement);

  average_cutsize_ = static_cast<double>(total_cutsize_) / number_of_runs_;

  if (rank_ == 0 && display_option_ > 0) {
    out_stream_ << endl
               << " --- PARTITIONING SUMMARY ---" << endl
               << "|" << endl
               << "|--- Cutsizes statistics:" << endl
               << "|" << endl
               << "|-- BEST = " << best_cutsize_ << endl
               << "|-- WORST = " << worst_cutsize_ << endl
               << "|-- AVE = " << average_cutsize_ << endl
               << "|" << endl
               << "|--- Time usage:" << endl
               << "|" << endl
               << "|-- TOTAL TIME = " << total_time_ << endl
               << "|-- AVE TIME = " << total_time_ / number_of_runs_ << endl
               << "|-- PARACOARSENING% = " << percentCoarsening << endl
               << "|-- SEQPARTITIONING% = " << percentSequential << endl
               << "|-- PARAREFINEMENT% = " << percentRefinement << endl
               << "|-- OTHER% = " << percentOther << endl
               << "|" << endl
               << " ----------------------------" << endl;
  }
}

void parallel_v_cycle_all_controller::print_type() const {
  out_stream_ << " type = ALL";
}

#endif
