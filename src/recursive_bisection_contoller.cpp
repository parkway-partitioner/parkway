#ifndef _RECUR_BISECT_CONTROLLER_CPP
#define _RECUR_BISECT_CONTROLLER_CPP

// ### RecurBisectController.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 4/1/2005: Last Modified
//
// ###

#include "recursive_bisection_contoller.hpp"
#include "hypergraph/parallel/hypergraph.hpp"
#include "hypergraph/serial/hypergraph.hpp"

namespace parallel = parkway::parallel;
namespace serial = parkway::serial;

recursive_bisection_contoller::recursive_bisection_contoller(bisection_controller *b,
                                             greedy_k_way_refiner *k, int rank,
                                             int nProcs, int nParts,
                                             int nBisectRuns, ostream &out)
    : sequential_controller(rank, nProcs, nParts, out) {
  bisector_ = b;
  refiner_ = k;
  number_of_bisections_ = nBisectRuns;
  log_k_ = 0;
  maximum_part_weight_ = 0;
  sum_of_cuts_ = 0;
  local_vertex_part_info_length_ = 0;
  number_of_partitions_ = 0;
  bisection_constraint_ = 0;
  average_part_weight_ = 0;
  average_initial_bisection_weight_ = 0;

  local_vertex_partition_info_.reserve(0);
  all_partition_info_.reserve(0);

#ifdef DEBUG_CONTROLLER
  assert(bisector);
#endif
}

recursive_bisection_contoller::~recursive_bisection_contoller() {
  dynamic_memory::delete_pointer<bisection_controller>(bisector_);
  dynamic_memory::delete_pointer<greedy_k_way_refiner>(refiner_);
}

void recursive_bisection_contoller::display_options() const {
  switch (display_option_) {
  case SILENT:
    break;

  default:

    out_stream_ << "|--- SEQ_CON:" << endl
               << "|- RBis:"
               << " seqR = " << number_of_runs_ << " bisR = " <<
                                                   number_of_bisections_
               << " pkT = " << accept_proportion_ << endl
               << "|" << endl;
#ifdef DEBUG_CONTROLLER
    assert(bisector);
#endif
      bisector_->display_options();

    break;
  }
}

void recursive_bisection_contoller::convToBisectionConstraints() {
#ifdef DEBUG_CONTROLLER
  assert(h);
#endif

  log_k_ = Funct::log2(number_of_parts_);

  int j;
  int i;

  double a = (1.0 + k_way_constraint_) / number_of_parts_;
  double b = 1.0 / log_k_;

  bisection_constraint_ = pow(a, b) - 0.5;

  j = hypergraph_->total_weight();

  average_part_weight_ = static_cast<double>(j) / number_of_parts_;
  average_initial_bisection_weight_ = static_cast<double>(j) / 2;
  maximum_part_weight_ = static_cast<int>(floor(
      average_part_weight_ + average_part_weight_ * k_way_constraint_));

  refiner_->set_maximum_part_weight(maximum_part_weight_);
  refiner_->set_average_part_weight(average_part_weight_);

// ###
// now initialise the partition
// vector structures
// ###

#ifdef DEBUG_CONTROLLER
  assert(numSeqRuns > 0);
#endif

  // ###
  // now determine how many of the sequential
  // runs' partitions  the processor should
  // k-way refine
  // ###

  if (number_of_runs_ <= number_of_processors_) {
    if (rank_ < number_of_runs_)
      number_of_partitions_ = 1;
    else
      number_of_partitions_ = 0;
  } else {
    i = Mod(number_of_runs_, number_of_processors_);
    j = number_of_runs_ / number_of_processors_;

    if (rank_ < i)
      number_of_partitions_ = j + 1;
    else
      number_of_partitions_ = j;
  }

  partition_vector_cuts_.reserve(number_of_partitions_);
  partition_vector_offsets_.reserve(number_of_partitions_ + 1);
  partition_vector_offsets_[0] = 0;

  j = hypergraph_->number_of_vertices();

  for (i = 1; i <= number_of_partitions_; ++i) {
    partition_vector_offsets_[i] = partition_vector_offsets_[i - 1] + j;
  }

  partition_vector_.reserve(partition_vector_offsets_[number_of_partitions_]);
}

void recursive_bisection_contoller::run(parallel::hypergraph &hgraph,
                                MPI_Comm comm) {
  initialize_coarsest_hypergraph(hgraph, comm);
  convToBisectionConstraints();

  if (display_option_ > 1 && rank_ == 0) {
    out_stream_ << "[R-B]: " << number_of_runs_ << " | ";
  }

  int i;
  int j;
  int ij;

  int numVertices = hypergraph_->number_of_vertices();
  int *pVector = nullptr;
  int destProcessor;
  int myPartitionIdx = 0;
  int v;

  dynamic_array<int> recvLens(number_of_processors_);
  dynamic_array<int> recvDispls(number_of_processors_);

  bisection *b;

  all_partition_info_.reserve(Shiftl(numVertices, 1));

  for (i = 0; i < number_of_runs_; ++i) {
    destProcessor = Mod(i, number_of_processors_);
    sum_of_cuts_ = 0;
    local_vertex_part_info_length_ = 0;

    if (rank_ == destProcessor) {
#ifdef DEBUG_CONTROLLER
      assert(myPartitionIdx < numMyPartitions);
#endif
      pVector = &partition_vector_[partition_vector_offsets_[myPartitionIdx]];
    }

    b = new bisection(hypergraph_, log_k_, 0);
    b->initMap();

    recursively_bisect(*b, comm);

    // ###
    // now recover the partition and
    // partition cutsize
    // ###

    MPI_Reduce(&sum_of_cuts_, &ij, 1, MPI_INT, MPI_SUM, destProcessor, comm);
    MPI_Gather(&local_vertex_part_info_length_, 1, MPI_INT, recvLens.data(), 1, MPI_INT,
               destProcessor, comm);

    if (rank_ == destProcessor) {
      partition_vector_cuts_[myPartitionIdx] = ij;
      ij = 0;

      for (j = 0; j < number_of_processors_; ++j) {
        recvDispls[j] = ij;
        ij += recvLens[j];
      }
#ifdef DEBUG_CONTROLLER
      assert(ij == Shiftl(numVertices, 1));
#endif
    }

    MPI_Gatherv(local_vertex_partition_info_.data(), local_vertex_part_info_length_, MPI_INT,
                all_partition_info_.data(), recvLens.data(),
                recvDispls.data(), MPI_INT, destProcessor, comm);

    if (rank_ == destProcessor) {
      ij = Shiftl(numVertices, 1);

      for (j = 0; j < ij;) {
        v = all_partition_info_[j++];
        pVector[v] = all_partition_info_[j++];
      }

      ++myPartitionIdx;

#ifdef DEBUG_CONTROLLER
      for (j = 0; j < numVertices; ++j)
        assert(pVector[j] >= 0 && pVector[j] < numParts);
#endif
    }

    dynamic_memory::delete_pointer<bisection>(b);
  }

  // ###
  // k-way refine local partitions
  // ###

  hypergraph_->set_number_of_partitions(number_of_partitions_);

  if (number_of_partitions_ > 0) {
    for (i = 0; i < number_of_partitions_; ++i) {
      pVector = &partition_vector_[partition_vector_offsets_[i]];
      hypergraph_->copy_in_partition(pVector, numVertices, i, partition_vector_cuts_[i]);
    }

    refiner_->rebalance(*hypergraph_);
  }

  // ###
  // project partitions
  // ###

  initialize_sequential_partitions(hgraph, comm);

#ifdef DEBUG_CONTROLLER
  hgraph.checkPartitions(numParts, maxPartWt, comm);
#endif

  dynamic_memory::delete_pointer<serial::hypergraph>(hypergraph_);
}

void recursive_bisection_contoller::initialize_sequential_partitions(
    parallel::hypergraph &hgraph,
    MPI_Comm comm) {
  int i;
  int j;
  int ij;

  int numTotVertices = hypergraph_->number_of_vertices();
  int ijk;
  int startOffset;
  int endOffset;
  int totToSend;

  int *hGraphPartitionVector;
  int *hGraphPartVectorOffsets;
  int *hGraphPartCuts;

  int *hPartitionVector = hypergraph_->partition_vector();
  int *hPartOffsetsVector = hypergraph_->partition_offsets();
  int *hPartitionCutsArray = hypergraph_->partition_cuts();

  dynamic_array<int> numVperProc(number_of_processors_);
  dynamic_array<int> procDispls(number_of_processors_);

  dynamic_array<int> sendLens(number_of_processors_);
  dynamic_array<int> sendDispls(number_of_processors_);
  dynamic_array<int> recvLens(number_of_processors_);
  dynamic_array<int> recvDispls(number_of_processors_);
  dynamic_array<int> sendArray;

  hgraph.set_number_of_partitions(number_of_runs_);

  hGraphPartitionVector = hgraph.partition_vector();
  hGraphPartVectorOffsets = hgraph.partition_offsets();
  hGraphPartCuts = hgraph.partition_cuts();

  // ###
  // communicate partition vector values
  // ###

  j = number_of_processors_ - 1;
  ij = numTotVertices / number_of_processors_;

  for (i = 0; i < j; ++i)
    numVperProc[i] = ij;

  numVperProc[i] = ij + Mod(numTotVertices, number_of_processors_);

  j = 0;
  ij = 0;

  for (i = 0; i < number_of_processors_; ++i) {
    sendDispls[i] = j;
    procDispls[i] = ij;
    sendLens[i] = numVperProc[i] * number_of_partitions_;
    j += sendLens[i];
    ij += numVperProc[i];
  }

  sendArray.reserve(j);
  totToSend = j;

  ij = 0;

  for (ijk = 0; ijk < number_of_processors_; ++ijk) {
    for (j = 0; j < number_of_partitions_; ++j) {
      startOffset = hPartOffsetsVector[j] + procDispls[ijk];
      endOffset = startOffset + numVperProc[ijk];

      for (i = startOffset; i < endOffset; ++i) {
        sendArray[ij++] = hPartitionVector[i];
      }
    }
  }
#ifdef DEBUG_CONTROLLER
  assert(ij == totToSend);
#endif

  MPI_Alltoall(sendLens.data(), 1, MPI_INT, recvLens.data(), 1, MPI_INT,
               comm);

  ij = 0;

  for (i = 0; i < number_of_processors_; ++i) {
    recvDispls[i] = ij;
    ij += recvLens[i];
  }

#ifdef DEBUG_CONTROLLER
  assert(ij == hGraphPartVectorOffsets[numSeqRuns]);
#endif

  MPI_Alltoallv(sendArray.data(), sendLens.data(),
                sendDispls.data(), MPI_INT, hGraphPartitionVector,
                recvLens.data(), recvDispls.data(), MPI_INT, comm);

  // ###
  // communicate partition cuts
  // ###

  MPI_Allgather(&number_of_partitions_, 1, MPI_INT, recvLens.data(), 1, MPI_INT,
                comm);

  ij = 0;

  for (i = 0; i < number_of_processors_; ++i) {
    recvDispls[i] = ij;
    ij += recvLens[i];
  }

  MPI_Allgatherv(hPartitionCutsArray, number_of_partitions_, MPI_INT, hGraphPartCuts,
                 recvLens.data(), recvDispls.data(), MPI_INT, comm);

  if (display_option_ > 1 && rank_ == 0) {
    for (i = 0; i < number_of_runs_; ++i)
      out_stream_ << hGraphPartCuts[i] << " ";

    out_stream_ << endl;
  }
}

void recursive_bisection_contoller::recursively_bisect(const bisection &b,
                                                       MPI_Comm comm) {
  int cut;
  int rank;
  int nProcs;
  int bisectAgain = b.bisect_again();

  serial::hypergraph *h = b.hypergraph();

  if (h->number_of_vertices() == 0)
    return;

  MPI_Comm_size(comm, &nProcs);
  MPI_Comm_rank(comm, &rank);

  if (nProcs == 1) {
    bisection *left = nullptr;
    bisection *right = nullptr;

    bisector_->set_number_of_sequential_runs(number_of_bisections_);

    if (bisectAgain == 1)
      bisector_->bisect(h, maximum_part_weight_);
    else
      bisector_->bisect(h, compute_maximum_weight(log_k_ - bisectAgain));

    sum_of_cuts_ += h->keep_best_partition();

    if (bisectAgain == 1) {
      int numVerts = h->number_of_vertices();
      int *partV = h->partition_vector();
      int *toOrigVmap = b.map_to_orig_vertices();
      int bisectionPart = b.part_id();

      int i;

      for (i = 0; i < numVerts; ++i) {
        local_vertex_partition_info_.assign(local_vertex_part_info_length_++, toOrigVmap[i]);

        if (partV[i] == 0)
          local_vertex_partition_info_.assign(local_vertex_part_info_length_++, bisectionPart);
        else
          local_vertex_partition_info_.assign(local_vertex_part_info_length_++,
                                      (bisectionPart | (1 << (log_k_ - 1))));
      }
    } else {
      split_bisection(b, left, right);
      recursively_bisect(*left, comm);
      recursively_bisect(*right, comm);

      if (left) {
        serial::hypergraph *hLeft = left->hypergraph();

        dynamic_memory::delete_pointer<serial::hypergraph>(hLeft);
        dynamic_memory::delete_pointer<bisection>(left);
      }

      if (right) {
        serial::hypergraph *hRight = right->hypergraph();

        dynamic_memory::delete_pointer<serial::hypergraph>(hRight);
        dynamic_memory::delete_pointer<bisection>(right);
      }
    }
  } else {
    bisection *newB;
    int bestCutProc;

    bisector_->set_number_of_sequential_runs(max(1, number_of_bisections_ / nProcs));

    if (bisectAgain == 1)
      bisector_->bisect(h, maximum_part_weight_);
    else
      bisector_->bisect(h, compute_maximum_weight(log_k_ - bisectAgain));

    cut = h->keep_best_partition();
    bestCutProc = best_partition_processor(cut, comm);

    if (rank == bestCutProc)
      sum_of_cuts_ += cut;

    if (bisectAgain == 1) {
      if (rank == bestCutProc) {
        int numVerts = h->number_of_vertices();
        int *partV = h->partition_vector();
        int *toOrigVmap = b.map_to_orig_vertices();
        int bisectionPart = b.part_id();

        int i;

        for (i = 0; i < numVerts; ++i) {
          local_vertex_partition_info_.assign(local_vertex_part_info_length_++, toOrigVmap[i]);

          if (partV[i] == 0)
            local_vertex_partition_info_.assign(local_vertex_part_info_length_++, bisectionPart);
          else
            local_vertex_partition_info_.assign(local_vertex_part_info_length_++,
                                        (bisectionPart | (1 << (log_k_ - 1))));
        }
      }
    } else {
      MPI_Bcast(h->partition_vector(), h->number_of_vertices(), MPI_INT,
                bestCutProc, comm);

      MPI_Comm new_comm;
      MPI_Comm_split(comm, And(rank, 0x1), 0, &new_comm);

      split_bisection(b, newB, comm);

      recursively_bisect(*newB, new_comm);

      if (newB) {
        serial::hypergraph *hNew = newB->hypergraph();

        dynamic_memory::delete_pointer<serial::hypergraph>(hNew);
        dynamic_memory::delete_pointer<bisection>(newB);
      }

      MPI_Comm_free(&new_comm);
    }
  }
}

void recursive_bisection_contoller::split_bisection(const bisection &b,
                                                    bisection *&newB,
                                                    MPI_Comm comm) const {
  int rank;
  int nProcs;

  MPI_Comm_size(comm, &nProcs);
  MPI_Comm_rank(comm, &rank);

  int i;
  int j;

  serial::hypergraph *newH;
  serial::hypergraph *h = b.hypergraph();

  // ###
  // h data_
  // ###

  int numHVertices = h->number_of_vertices();
  int numHHedges = h->number_of_hyperedges();
  int bPartID = b.part_id();
  int bisectAgain = b.bisect_again();

  int *mapToOrig = b.map_to_orig_vertices();
  int *hPartVector = h->partition_vector();
  int *hVertWt = h->vertex_weights();
  int *hHedgeWt = h->hyperedge_weights();
  int *hHedgeOffsets = h->hyperedge_offsets();
  int *hPinList = h->pin_list();

  // ###
  // newH data_
  // ###

  int numVerts = 0;
  int numHedges = 0;
  int totWt = 0;
  int numPins = 0;

  dynamic_array<int> *vertWt = new dynamic_array<int>(64);
  dynamic_array<int> *mapOrig = new dynamic_array<int>(64);
  dynamic_array<int> *hedgeWts = new dynamic_array<int>(64);
  dynamic_array<int> *hedgeOffsets = new dynamic_array<int>(64);
  dynamic_array<int> *pinList = new dynamic_array<int>(64);

  // ###
  // auxiliary data_
  // ###

  int v;
  int hEdgeLen;
  int endOffset;

  dynamic_array<int> mapFromHtoNewH(numHVertices);

  if (And(rank, 0x1)) {
    // ###
    // in the odd processor case
    // ###

    for (i = 0; i < numHVertices; ++i) {
#ifdef DEBUG_CONTROLLER
      assert(hPartVector[i] == 1 || hPartVector[i] == 0);
#endif
      if (hPartVector[i] == 1) {
        vertWt->assign(numVerts, hVertWt[i]);
        mapOrig->assign(numVerts, mapToOrig[i]);
        mapFromHtoNewH[i] = numVerts++;
        totWt += hVertWt[i];
      } else {
        mapFromHtoNewH[i] = -1;
      }
    }

    vertWt->reserve(numVerts);
    mapOrig->reserve(numVerts);

    newH = new serial::hypergraph(vertWt->data(), numVerts);

    // ###
    // initialise pin list
    // ###

    hedgeOffsets->assign(numHedges, numPins);

    for (i = 0; i < numHHedges; ++i) {
      endOffset = hHedgeOffsets[i + 1];
      hEdgeLen = 0;

      for (j = hHedgeOffsets[i]; j < endOffset; ++j) {
        v = hPinList[j];

        if (hPartVector[v] == 1) {
#ifdef DEBUG_CONTROLLER
          assert(mapFromHtoNewH[v] != -1);
#endif
          pinList->assign(numPins + (hEdgeLen++), mapFromHtoNewH[v]);
        }
      }

      if (hEdgeLen > 1) {
        numPins += hEdgeLen;
        hedgeWts->assign(numHedges++, hHedgeWt[i]);
        hedgeOffsets->assign(numHedges, numPins);
      }
    }

    hedgeWts->reserve(numHedges);
    hedgeOffsets->reserve(numHedges + 1);
    pinList->reserve(numPins);

    // ###
    // now init the hypergraphs
    // ###

    newH->set_number_of_hyperedges(numHedges);
    newH->set_number_of_pins(numPins);
    newH->set_total_weight(totWt);
    newH->set_hyperedge_weights(hedgeWts->data(), numHedges);
    newH->set_pin_list(pinList->data(), numPins);
    newH->set_hyperedge_offsets(hedgeOffsets->data(), numHedges + 1);

    newH->buildVtoHedges();

    // ###
    // now init the bisections
    // ###

    newB = new bisection(newH, bisectAgain - 1,
                         Or(bPartID, Shiftl(1, (log_k_ - bisectAgain))));
    newB->setMap(mapOrig->data(), numVerts);
  } else {
    // ###
    // in the even processor case
    // ###

    for (i = 0; i < numHVertices; ++i) {
#ifdef DEBUG_CONTROLLER
      assert(hPartVector[i] == 1 || hPartVector[i] == 0);
#endif
      if (hPartVector[i] == 0) {
        vertWt->assign(numVerts, hVertWt[i]);
        mapOrig->assign(numVerts, mapToOrig[i]);
        mapFromHtoNewH[i] = numVerts++;
        totWt += hVertWt[i];
      } else {
        mapFromHtoNewH[i] = -1;
      }
    }

    vertWt->reserve(numVerts);
    mapOrig->reserve(numVerts);

    newH = new serial::hypergraph(vertWt->data(), numVerts);

    // ###
    // initialise pin list
    // ###

    hedgeOffsets->assign(numHedges, numPins);

    for (i = 0; i < numHHedges; ++i) {
      endOffset = hHedgeOffsets[i + 1];
      hEdgeLen = 0;

      for (j = hHedgeOffsets[i]; j < endOffset; ++j) {
        v = hPinList[j];

        if (hPartVector[v] == 0) {
#ifdef DEBUG_CONTROLLER
          assert(mapFromHtoNewH[v] != -1);
#endif
          pinList->assign(numPins + (hEdgeLen++), mapFromHtoNewH[v]);
        }
      }

      if (hEdgeLen > 1) {
        numPins += hEdgeLen;
        hedgeWts->assign(numHedges++, hHedgeWt[i]);
        hedgeOffsets->assign(numHedges, numPins);
      }
    }

    hedgeWts->reserve(numHedges);
    hedgeOffsets->reserve(numHedges + 1);
    pinList->reserve(numPins);

    // ###
    // now init the hypergraph
    // ###

    newH->set_number_of_hyperedges(numHedges);
    newH->set_number_of_pins(numPins);
    newH->set_total_weight(totWt);
    newH->set_hyperedge_weights(hedgeWts->data(), numHedges);
    newH->set_pin_list(pinList->data(), numPins);
    newH->set_hyperedge_offsets(hedgeOffsets->data(), numHedges + 1);

    newH->buildVtoHedges();

    // ###
    // now init the bisections
    // ###

    newB = new bisection(newH, bisectAgain - 1, bPartID);
    newB->setMap(mapOrig->data(), numVerts);
  }
}

void recursive_bisection_contoller::split_bisection(const bisection &b,
                                                    bisection *&l,
                                                    bisection *&r) const {
  int i;
  int j;

  serial::hypergraph *leftH;
  serial::hypergraph *rightH;
  serial::hypergraph *h = b.hypergraph();

  // ###
  // h data_
  // ###

  int numHVertices = h->number_of_vertices();
  int numHHedges = h->number_of_hyperedges();
  int bPartID = b.part_id();
  int bisectAgain = b.bisect_again();

  int *mapToOrig = b.map_to_orig_vertices();
  int *hPartVector = h->partition_vector();
  int *hVertWt = h->vertex_weights();
  int *hHedgeWt = h->hyperedge_weights();
  int *hHedgeOffsets = h->hyperedge_offsets();
  int *hPinList = h->pin_list();

  // ###
  // leftH data_
  // ###

  int numLeftVerts = 0;
  int numLeftHedges = 0;
  int totLeftWt = 0;
  int numLeftPins = 0;

  dynamic_array<int> *leftVertWt = new dynamic_array<int>(64);
  dynamic_array<int> *leftMapOrig = new dynamic_array<int>(64);
  dynamic_array<int> *leftHedgeWts = new dynamic_array<int>(64);
  dynamic_array<int> *leftHedgeOffsets = new dynamic_array<int>(64);
  dynamic_array<int> *leftPinList = new dynamic_array<int>(64);

  // ###
  // rightH data_
  // ###

  int numRightVerts = 0;
  int numRightHedges = 0;
  int totRightWt = 0;
  int numRightPins = 0;

  dynamic_array<int> *rightVertWt = new dynamic_array<int>(64);
  dynamic_array<int> *rightMapOrig = new dynamic_array<int>(64);
  dynamic_array<int> *rightHedgeWts = new dynamic_array<int>(64);
  dynamic_array<int> *rightHedgeOffsets = new dynamic_array<int>(64);
  dynamic_array<int> *rightPinList = new dynamic_array<int>(64);

  // ###
  // auxiliary data_
  // ###

  int v;
  int endOffset;
  int leftHedgeLen;
  int rightHedgeLen;

  dynamic_array<int> mapFromHtoNewH(numHVertices);

  for (i = 0; i < numHVertices; ++i) {
#ifdef DEBUG_CONTROLLER
    assert(hPartVector[i] == 1 || hPartVector[i] == 0);
#endif
    if (hPartVector[i] == 0) {
      leftVertWt->assign(numLeftVerts, hVertWt[i]);
      leftMapOrig->assign(numLeftVerts, mapToOrig[i]);
      mapFromHtoNewH[i] = numLeftVerts++;
      totLeftWt += hVertWt[i];
    } else {
      rightVertWt->assign(numRightVerts, hVertWt[i]);
      rightMapOrig->assign(numRightVerts, mapToOrig[i]);
      mapFromHtoNewH[i] = numRightVerts++;
      totRightWt += hVertWt[i];
    }
  }

  leftVertWt->reserve(numLeftVerts);
  leftMapOrig->reserve(numLeftVerts);

  rightVertWt->reserve(numRightVerts);
  rightMapOrig->reserve(numRightVerts);

  leftH = new serial::hypergraph(leftVertWt->data(), numLeftVerts);
  rightH = new serial::hypergraph(rightVertWt->data(), numRightVerts);

  // ###
  // initialise pin list
  // ###

  leftHedgeOffsets->assign(numLeftHedges, numLeftPins);
  rightHedgeOffsets->assign(numRightHedges, numRightPins);

  for (i = 0; i < numHHedges; ++i) {
    endOffset = hHedgeOffsets[i + 1];

    leftHedgeLen = 0;
    rightHedgeLen = 0;

    for (j = hHedgeOffsets[i]; j < endOffset; ++j) {
      v = hPinList[j];

      if (hPartVector[v] == 0) {
        leftPinList->assign(numLeftPins + (leftHedgeLen++), mapFromHtoNewH[v]);
      } else {
        rightPinList->assign(numRightPins + (rightHedgeLen++),
                             mapFromHtoNewH[v]);
      }
    }

    if (leftHedgeLen > 1) {
      numLeftPins += leftHedgeLen;
      leftHedgeWts->assign(numLeftHedges++, hHedgeWt[i]);
      leftHedgeOffsets->assign(numLeftHedges, numLeftPins);
    }

    if (rightHedgeLen > 1) {
      numRightPins += rightHedgeLen;
      rightHedgeWts->assign(numRightHedges++, hHedgeWt[i]);
      rightHedgeOffsets->assign(numRightHedges, numRightPins);
    }
  }

  leftHedgeWts->reserve(numLeftHedges);
  leftHedgeOffsets->reserve(numLeftHedges + 1);
  leftPinList->reserve(numLeftPins);

  rightHedgeWts->reserve(numRightHedges);
  rightHedgeOffsets->reserve(numRightHedges + 1);
  rightPinList->reserve(numRightPins);

  // ###
  // now init the hypergraphs
  // ###

  leftH->set_number_of_hyperedges(numLeftHedges);
  leftH->set_number_of_pins(numLeftPins);
  leftH->set_total_weight(totLeftWt);
  leftH->set_hyperedge_weights(leftHedgeWts->data(), numLeftHedges);
  leftH->set_pin_list(leftPinList->data(), numLeftPins);
  leftH->set_hyperedge_offsets(leftHedgeOffsets->data(), numLeftHedges + 1);

  leftH->buildVtoHedges();

  rightH->set_number_of_hyperedges(numRightHedges);
  rightH->set_number_of_pins(numRightPins);
  rightH->set_total_weight(totRightWt);
  rightH->set_hyperedge_weights(rightHedgeWts->data(), numRightHedges);
  rightH->set_pin_list(rightPinList->data(), numRightPins);
  rightH->set_hyperedge_offsets(rightHedgeOffsets->data(),
                                numRightHedges + 1);

  rightH->buildVtoHedges();

  // ###
  // now init the bisections
  // ###

  l = new bisection(leftH, bisectAgain - 1, bPartID);
  r = new bisection(rightH, bisectAgain - 1,
                    Or(bPartID, Shiftl(1, (log_k_ - bisectAgain))));

  l->setMap(leftMapOrig->data(), numLeftVerts);
  r->setMap(rightMapOrig->data(), numRightVerts);
}

int recursive_bisection_contoller::best_partition_processor(int cut,
                                                            MPI_Comm comm) const {
  int rank;
  int nProcs;
  int bestCut;
  int bestProc;

  int i;

  MPI_Comm_size(comm, &nProcs);
  MPI_Comm_rank(comm, &rank);

  dynamic_array<int> allCuts(nProcs);
  dynamic_array<int> procs(nProcs);

  MPI_Allgather(&cut, 1, MPI_INT, allCuts.data(), 1, MPI_INT, comm);

  bestCut = allCuts[0];
  procs[0] = 0;

  for (i = 1; i < nProcs; ++i) {
    procs[i] = i;

    if (allCuts[i] < bestCut)
      bestCut = allCuts[i];
  }

  Funct::randomPermutation(procs.data(), nProcs);

  if (rank == 0) {
    for (i = 0; i < nProcs; ++i) {
      if (allCuts[procs[i]] == bestCut) {
        bestProc = procs[i];
        break;
      }
    }
  }

  MPI_Bcast(&bestProc, 1, MPI_INT, 0, comm);

  return bestProc;
}

int recursive_bisection_contoller::compute_maximum_weight(int numBs) const {
  double maxPartWeight = recursively_compute_maximum(
      average_initial_bisection_weight_, numBs);

  return (static_cast<int>(floor(maxPartWeight)));
}

double recursive_bisection_contoller::recursively_compute_maximum(
    double currAve,
    int depth) const {
  if (depth == 0)
    return (currAve + currAve * bisection_constraint_);
  return (recursively_compute_maximum(
      (currAve + currAve * bisection_constraint_) / 2,
      depth - 1));
}

#endif
