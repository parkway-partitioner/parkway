// ### ParaVCycleFinalController.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 10/1/2005: Last Modified
//
// ###

#include "parallel_v_cycle_final_controller.hpp"

parallel_v_cycle_final_controller::parallel_v_cycle_final_controller(
    parallel::restrictive_coarsening &rc, parallel::coarsener &c, parallel::refiner &r,
    parkway::serial::controller &ref, int rank, int nP, int percentile, int inc,
    int approxRef, int limit, double limitAsPercent, std::ostream &out)
    : parallel_v_cycle_controller(rc, c, r, ref, rank, nP, percentile, inc, approxRef,
                           limit, limitAsPercent, out) {}

parallel_v_cycle_final_controller::~parallel_v_cycle_final_controller() {}

void parallel_v_cycle_final_controller::run(MPI_Comm comm) {
  int hEdgePercentile;
  int firstCutSize;
  int secondCutSize;
  int i;
  int diffInCutSize;
  int numIterations;
  int vCycleGain;
  int minIterationGain;

  int percentCoarsening;
  int percentserial;
  int percentRefinement;
  int percentOther;

  double totStartTime;

  std::stack<int> hEdgePercentiles;

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

#ifdef DEBUG_CONTROLLER
  int checkCutsize;
#endif

  for (i = 0; i < number_of_runs_; ++i) {
    if (shuffled_ == 1)
        hypergraph_->shuffle_vertices_randomly(map_to_orig_vertices_, comm);

    if (shuffled_ == 2)
        hypergraph_->prescribed_vertex_shuffle(map_to_orig_vertices_,
                                          shuffle_partition_, comm);

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
        hEdgePercentiles.push(std::min(hEdgePercentile + percentile_increment_, 100));
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

    coarseGraph = hypergraphs_.top();
    hypergraphs_.pop();
    serial_controller_.run(*coarseGraph, comm);

    MPI_Barrier(comm);
    total_serial_time_ += (MPI_Wtime() - start_time_);

    // ###
    // uncoarsen the initial partition
    // ###

    MPI_Barrier(comm);
    start_time_ = MPI_Wtime();

    while (hypergraphs_.size() > 0) {
        coarseGraph->remove_bad_partitions(keep_partitions_within_ *
                                           accumulator_);
      accumulator_ *= reduction_in_keep_threshold_;

      hEdgePercentile = hEdgePercentiles.top();
      hEdgePercentiles.pop();
      finerGraph = hypergraphs_.top();
      hypergraphs_.pop();

      if (finerGraph == hypergraph_)
        project_v_cycle_partition(*coarseGraph, *finerGraph, comm);
      else
          finerGraph->project_partitions(*coarseGraph, comm);

#ifdef DEBUG_CONTROLLER
      finerGraph->checkPartitions(numTotalParts, maxPartWt, comm);
#endif
      if (approximate_refine_)
        refiner_.set_percentile(hEdgePercentile);

      refiner_.refine(*finerGraph, comm);

#ifdef DEBUG_CONTROLLER
      finerGraph->checkPartitions(numTotalParts, maxPartWt, comm);
#endif
      coarseGraph = finerGraph;
    }

#ifdef DEBUG_CONTROLLER
    assert(coarseGraph == hgraph);
#endif

    MPI_Barrier(comm);
    total_refinement_time_ += (MPI_Wtime() - start_time_);

    refiner_.release_memory();

    // ###
    // select the best partition
    // ###

    firstCutSize = coarseGraph->keep_best_partition();

#ifdef DEBUG_CONTROLLER
    checkCutsize = coarseGraph->calcCutsize(numTotalParts, 0, comm);
    assert(firstCutSize == checkCutsize);
#endif

    numIterations = 0;
    vCycleGain = 0;
    record_v_cycle_partition(*hypergraph_, numIterations++);

    if (display_option_ > 1 && rank_ == 0) {
      out_stream_ << "\t ------ PARALLEL V-CYCLE CALL ------" << std::endl;
    }

    do {
      minIterationGain =
          static_cast<int>(floor(limit_as_percent_of_cut_ * firstCutSize));

      shuffle_v_cycle_vertices_by_partition(*hypergraph_, comm);
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
        restrictive_coarsening_.set_percentile(hEdgePercentile);
        coarseGraph = restrictive_coarsening_.coarsen(*finerGraph, comm);

        if (coarseGraph) {
          hEdgePercentiles.push(
              std::min(hEdgePercentile + percentile_increment_, 100));
          hypergraphs_.push(coarseGraph);
          finerGraph = coarseGraph;
        }
      }

      while (coarseGraph);

      MPI_Barrier(comm);
      total_coarsening_time_ += (MPI_Wtime() - start_time_);

      restrictive_coarsening_.release_memory();

      // ###
      // compute the initial partition
      // ###

      coarseGraph = hypergraphs_.top();
      hypergraphs_.pop();
      coarseGraph->set_number_of_partitions(0);
      coarseGraph->shift_vertices_to_balance(comm);

      MPI_Barrier(comm);
      start_time_ = MPI_Wtime();

      serial_controller_.run(*coarseGraph, comm);

      MPI_Barrier(comm);
      total_serial_time_ += (MPI_Wtime() - start_time_);

      // ###
      // uncoarsen the initial partition
      // ###

      MPI_Barrier(comm);
      start_time_ = MPI_Wtime();
      accumulator_ = 1.0;

      while (hypergraphs_.size() > 0) {
          coarseGraph->remove_bad_partitions(
              keep_partitions_within_ * accumulator_);
        accumulator_ *= reduction_in_keep_threshold_;

        hEdgePercentile = hEdgePercentiles.top();
        hEdgePercentiles.pop();
        finerGraph = hypergraphs_.top();
        hypergraphs_.pop();

        if (finerGraph == hypergraph_)
          shift_v_cycle_vertices_to_balance(*finerGraph, comm);
        else
            finerGraph->shift_vertices_to_balance(comm);

          finerGraph->project_partitions(*coarseGraph, comm);

#ifdef DEBUG_CONTROLLER
        finerGraph->checkPartitions(numTotalParts, maxPartWt, comm);
#endif
        if (random_shuffle_before_refine_) {
          if (hypergraphs_.size() == 0)
              finerGraph->shuffle_vertices_randomly(map_to_inter_vertices_, comm);
          else
              finerGraph->shuffle_vertices_randomly(*(hypergraphs_.top()), comm);
        }

        if (approximate_refine_)
          refiner_.set_percentile(hEdgePercentile);

        refiner_.refine(*finerGraph, comm);

#ifdef DEBUG_CONTROLLER
        finerGraph->checkPartitions(numTotalParts, maxPartWt, comm);
#endif
        coarseGraph = finerGraph;
      }

#ifdef DEBUG_CONTROLLER
      assert(coarseGraph == hgraph);
      checkCutsize = coarseGraph->calcCutsize(numTotalParts, 0, comm);
#endif
      MPI_Barrier(comm);
      total_refinement_time_ += (MPI_Wtime() - start_time_);

      refiner_.release_memory();

      // ###
      // select the best partition
      // ###

      secondCutSize = coarseGraph->keep_best_partition();
      diffInCutSize = firstCutSize - secondCutSize;

      if (display_option_ > 1 && rank_ == 0) {
        out_stream_ << "\t ------ [" << numIterations << "] " << diffInCutSize
                   << std::endl;
      }

#ifdef DEBUG_CONTROLLER
      checkCutsize = coarseGraph->calcCutsize(numTotalParts, 0, comm);
      assert(secondCutSize == checkCutsize);
#endif

      if (diffInCutSize > 0 && diffInCutSize < minIterationGain)
        break;

      if (numIterations == limit_on_cycles_)
        break;

      if (firstCutSize > secondCutSize) {
        record_v_cycle_partition(*coarseGraph, numIterations++);
        firstCutSize = secondCutSize;
        vCycleGain += diffInCutSize;
      }
    } while (diffInCutSize > 0);

    gather_in_v_cycle_partition(*hypergraph_, firstCutSize, comm);
    update_map_to_orig_vertices(comm);

    if (display_option_ > 1 && rank_ == 0) {
      out_stream_ << "\t ------ " << vCycleGain << " ------" << std::endl;
    }

#ifdef DEBUG_CONTROLLER
    assert(hgraphs.getNumElem() == 0);
    assert(hgraph->getNumLocalVertices() == numOrigLocVerts);
    hgraph->checkPartitions(numTotalParts, maxPartWt, comm);
#endif

    if (rank_ == 0 && display_option_ > 0) {
      out_stream_ << "\nPRUN[" << i << "] = " << firstCutSize << std::endl << std::endl;
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

void parallel_v_cycle_final_controller::print_type() const {
  out_stream_ << " type = BIG";
}
