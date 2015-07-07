// ### RecurBisectController.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 4/1/2005: Last Modified
//
// ###
#include "controllers/serial/recursive_bisection_contoller.hpp"
#include "hypergraph/parallel/hypergraph.hpp"
#include "hypergraph/serial/hypergraph.hpp"
#include "utility/logging.hpp"
#include "utility/math.hpp"

namespace parallel = parkway::parallel;
namespace serial = parkway::serial;

recursive_bisection_contoller::recursive_bisection_contoller(
    serial::bisection_controller *b, serial::greedy_k_way_refiner *k, int rank,
    int nProcs, int nParts, int nBisectRuns)
    : parkway::serial::controller(rank, nProcs, nParts) {
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
}

void recursive_bisection_contoller::display_options() const {
  info("|--- SEQ_CON:\n"
       "|- RBis: seqR = %i bisR = %i pkT = %.2f\n|\n", number_of_runs_,
       number_of_bisections_, accept_proportion_);
  bisector_->display_options();
}

void recursive_bisection_contoller::convToBisectionConstraints() {
#ifdef DEBUG_CONTROLLER
  assert(h);
#endif

  log_k_ = parkway::utility::math::log2(number_of_parts_);

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
  // now determine how many of the serial
  // runs' partitions  the processor should
  // k-way refine
  // ###

  if (number_of_runs_ <= number_of_processors_) {
    if (rank_ < number_of_runs_)
      number_of_partitions_ = 1;
    else
      number_of_partitions_ = 0;
  } else {
    i = number_of_runs_ % number_of_processors_;
    j = number_of_runs_ / number_of_processors_;

    if (rank_ < i)
      number_of_partitions_ = j + 1;
    else
      number_of_partitions_ = j;
  }

  partition_vector_cuts_.resize(number_of_partitions_);
  partition_vector_offsets_.resize(number_of_partitions_ + 1);
  partition_vector_offsets_[0] = 0;

  j = hypergraph_->number_of_vertices();

  for (i = 1; i <= number_of_partitions_; ++i) {
    partition_vector_offsets_[i] = partition_vector_offsets_[i - 1] + j;
  }

  partition_vector_.resize(partition_vector_offsets_[number_of_partitions_]);
}

void recursive_bisection_contoller::run(parallel::hypergraph &hgraph,
                                        MPI_Comm comm) {
  initialize_coarsest_hypergraph(hgraph, comm);
  convToBisectionConstraints();

  progress("[R-B]: %i |", number_of_runs_);

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

  all_partition_info_.resize(numVertices << 1);

  for (i = 0; i < number_of_runs_; ++i) {
    destProcessor = i % number_of_processors_;
    sum_of_cuts_ = 0;
    local_vertex_part_info_length_ = 0;

    if (rank_ == destProcessor) {
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
    }

    MPI_Gatherv(local_vertex_partition_info_.data(), local_vertex_part_info_length_, MPI_INT,
                all_partition_info_.data(), recvLens.data(),
                recvDispls.data(), MPI_INT, destProcessor, comm);

    if (rank_ == destProcessor) {
      ij = numVertices << 1;

      for (j = 0; j < ij;) {
        v = all_partition_info_[j++];
        pVector[v] = all_partition_info_[j++];
      }

      ++myPartitionIdx;
    }
  }

  // ###
  // k-way refine local partitions
  // ###

  hypergraph_->set_number_of_partitions(number_of_partitions_);

  if (number_of_partitions_ > 0) {
    for (i = 0; i < number_of_partitions_; ++i) {
      int *start = &(partition_vector_.data()[partition_vector_offsets_[i]]);
      dynamic_array<int> p_vector(numVertices);
      p_vector.set_data(start, numVertices);
      hypergraph_->copy_in_partition(p_vector, numVertices, i, partition_vector_cuts_[i]);
    }

    refiner_->rebalance(*hypergraph_);
  }

  // ###
  // project partitions
  // ###

  initialize_serial_partitions(hgraph, comm);

#ifdef DEBUG_CONTROLLER
  hgraph.checkPartitions(numParts, maxPartWt, comm);
#endif
}

void recursive_bisection_contoller::initialize_serial_partitions(
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

  ds::dynamic_array<int> hGraphPartitionVector;
  ds::dynamic_array<int> hGraphPartVectorOffsets;
  ds::dynamic_array<int> hGraphPartCuts;

  auto hPartitionVector = hypergraph_->partition_vector();
  auto hPartOffsetsVector = hypergraph_->partition_offsets();
  auto hPartitionCutsArray = hypergraph_->partition_cuts();

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

  numVperProc[i] = ij + (numTotVertices % number_of_processors_);

  j = 0;
  ij = 0;

  for (i = 0; i < number_of_processors_; ++i) {
    sendDispls[i] = j;
    procDispls[i] = ij;
    sendLens[i] = numVperProc[i] * number_of_partitions_;
    j += sendLens[i];
    ij += numVperProc[i];
  }

  sendArray.resize(j);
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
                sendDispls.data(), MPI_INT, hGraphPartitionVector.data(),
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

  MPI_Allgatherv(hPartitionCutsArray.data(), number_of_partitions_, MPI_INT,
                 hGraphPartCuts.data(), recvLens.data(), recvDispls.data(),
                 MPI_INT, comm);

  for (i = 0; i < number_of_runs_; ++i) {
    progress("%i ", hGraphPartCuts[i]);
  }
  progress("\n");
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

    bisector_->set_number_of_serial_runs(number_of_bisections_);

    if (bisectAgain == 1)
      bisector_->bisect(h, maximum_part_weight_);
    else
      bisector_->bisect(h, compute_maximum_weight(log_k_ - bisectAgain));

    sum_of_cuts_ += h->keep_best_partition();

    if (bisectAgain == 1) {
      int numVerts = h->number_of_vertices();
      auto partV = h->partition_vector();
      auto toOrigVmap = b.map_to_orig_vertices();
      int bisectionPart = b.part_id();

      int i;

      for (i = 0; i < numVerts; ++i) {
        local_vertex_partition_info_[local_vertex_part_info_length_++] =
            toOrigVmap[i];

        if (partV[i] == 0)
          local_vertex_partition_info_[local_vertex_part_info_length_++] =
              bisectionPart;
        else
          local_vertex_partition_info_[local_vertex_part_info_length_++] =
              (bisectionPart | (1 << (log_k_ - 1)));
      }
    } else {
      split_bisection(b, left, right);
      recursively_bisect(*left, comm);
      recursively_bisect(*right, comm);

      if (left) {
        serial::hypergraph *hLeft = left->hypergraph();
      }

      if (right) {
        serial::hypergraph *hRight = right->hypergraph();
      }
    }
  } else {
    bisection *newB;
    int bestCutProc;

    bisector_->set_number_of_serial_runs(std::max(1, number_of_bisections_ / nProcs));

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
        auto partV = h->partition_vector();
        auto toOrigVmap = b.map_to_orig_vertices();
        int bisectionPart = b.part_id();

        int i;

        for (i = 0; i < numVerts; ++i) {
          local_vertex_partition_info_[local_vertex_part_info_length_++] = toOrigVmap[i];

          if (partV[i] == 0)
            local_vertex_partition_info_[local_vertex_part_info_length_++] = bisectionPart;
          else
            local_vertex_partition_info_[local_vertex_part_info_length_++] = (bisectionPart | (1 << (log_k_ - 1)));
        }
      }
    } else {
      MPI_Bcast(h->partition_vector().data(), h->number_of_vertices(),
                MPI_INT, bestCutProc, comm);

      MPI_Comm new_comm;
      MPI_Comm_split(comm, (rank & 0x1), 0, &new_comm);

      split_bisection(b, newB, comm);

      recursively_bisect(*newB, new_comm);

      if (newB) {
        serial::hypergraph *hNew = newB->hypergraph();
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

  auto mapToOrig = b.map_to_orig_vertices();
  auto hPartVector = h->partition_vector();
  auto hVertWt = h->vertex_weights();
  auto hHedgeWt = h->hyperedge_weights();
  auto hHedgeOffsets = h->hyperedge_offsets();
  auto hPinList = h->pin_list();

  // ###
  // newH data_
  // ###

  int numVerts = 0;
  int numHedges = 0;
  int totWt = 0;
  int numPins = 0;

  dynamic_array<int> vertWt(64);
  dynamic_array<int> mapOrig(64);
  dynamic_array<int> hedgeWts(64);
  dynamic_array<int> hedgeOffsets(64);
  dynamic_array<int> pinList(64);

  // ###
  // auxiliary data_
  // ###

  int v;
  int hEdgeLen;
  int endOffset;

  dynamic_array<int> mapFromHtoNewH(numHVertices);

  if (rank & 0x1) {
    // ###
    // in the odd processor case
    // ###

    for (i = 0; i < numHVertices; ++i) {
#ifdef DEBUG_CONTROLLER
      assert(hPartVector[i] == 1 || hPartVector[i] == 0);
#endif
      if (hPartVector[i] == 1) {
        vertWt[numVerts] = hVertWt[i];
        mapOrig[numVerts] = mapToOrig[i];
        mapFromHtoNewH[i] = numVerts++;
        totWt += hVertWt[i];
      } else {
        mapFromHtoNewH[i] = -1;
      }
    }

    vertWt.resize(numVerts);
    mapOrig.resize(numVerts);

    newH = new serial::hypergraph(vertWt, numVerts);

    // ###
    // initialise pin list
    // ###

    hedgeOffsets[numHedges] = numPins;

    for (i = 0; i < numHHedges; ++i) {
      endOffset = hHedgeOffsets[i + 1];
      hEdgeLen = 0;

      for (j = hHedgeOffsets[i]; j < endOffset; ++j) {
        v = hPinList[j];

        if (hPartVector[v] == 1) {
          pinList[numPins + (hEdgeLen++)] = mapFromHtoNewH[v];
        }
      }

      if (hEdgeLen > 1) {
        numPins += hEdgeLen;
        hedgeWts[numHedges++] = hHedgeWt[i];
        hedgeOffsets[numHedges] = numPins;
      }
    }

    hedgeWts.resize(numHedges);
    hedgeOffsets.resize(numHedges + 1);
    pinList.resize(numPins);

    // ###
    // now init the hypergraphs
    // ###

    newH->set_number_of_hyperedges(numHedges);
    newH->set_number_of_pins(numPins);
    newH->set_total_weight(totWt);
    newH->set_hyperedge_weights(hedgeWts);
    newH->set_pin_list(pinList);
    newH->set_hyperedge_offsets(hedgeOffsets);

    newH->buildVtoHedges();

    // ###
    // now init the bisections
    // ###

    newB = new bisection(newH, bisectAgain - 1,
                         (bPartID | (1 << (log_k_ - bisectAgain))));
    newB->setMap(mapOrig);
  } else {
    // ###
    // in the even processor case
    // ###

    for (i = 0; i < numHVertices; ++i) {
#ifdef DEBUG_CONTROLLER
      assert(hPartVector[i] == 1 || hPartVector[i] == 0);
#endif
      if (hPartVector[i] == 0) {
        vertWt[numVerts] = hVertWt[i];
        mapOrig[numVerts] = mapToOrig[i];
        mapFromHtoNewH[i] = numVerts++;
        totWt += hVertWt[i];
      } else {
        mapFromHtoNewH[i] = -1;
      }
    }

    vertWt.resize(numVerts);
    mapOrig.resize(numVerts);

    newH = new serial::hypergraph(vertWt, numVerts);

    // ###
    // initialise pin list
    // ###

    hedgeOffsets[numHedges] = numPins;

    for (i = 0; i < numHHedges; ++i) {
      endOffset = hHedgeOffsets[i + 1];
      hEdgeLen = 0;

      for (j = hHedgeOffsets[i]; j < endOffset; ++j) {
        v = hPinList[j];

        if (hPartVector[v] == 0) {
#ifdef DEBUG_CONTROLLER
          assert(mapFromHtoNewH[v] != -1);
#endif
          pinList[numPins + (hEdgeLen++)] = mapFromHtoNewH[v];
        }
      }

      if (hEdgeLen > 1) {
        numPins += hEdgeLen;
        hedgeWts[numHedges++] = hHedgeWt[i];
        hedgeOffsets[numHedges] = numPins;
      }
    }

    hedgeWts.resize(numHedges);
    hedgeOffsets.resize(numHedges + 1);
    pinList.resize(numPins);

    // ###
    // now init the hypergraph
    // ###

    newH->set_number_of_hyperedges(numHedges);
    newH->set_number_of_pins(numPins);
    newH->set_total_weight(totWt);
    newH->set_hyperedge_weights(hedgeWts);
    newH->set_pin_list(pinList);
    newH->set_hyperedge_offsets(hedgeOffsets);

    newH->buildVtoHedges();

    // ###
    // now init the bisections
    // ###

    newB = new bisection(newH, bisectAgain - 1, bPartID);
    newB->setMap(mapOrig);
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

  auto mapToOrig = b.map_to_orig_vertices();
  auto hPartVector = h->partition_vector();
  auto hVertWt = h->vertex_weights();
  auto hHedgeWt = h->hyperedge_weights();
  auto hHedgeOffsets = h->hyperedge_offsets();
  auto hPinList = h->pin_list();

  // ###
  // leftH data_
  // ###

  int numLeftVerts = 0;
  int numLeftHedges = 0;
  int totLeftWt = 0;
  int numLeftPins = 0;

  dynamic_array<int> leftVertWt(64);
  dynamic_array<int> leftMapOrig(64);
  dynamic_array<int> leftHedgeWts(64);
  dynamic_array<int> leftHedgeOffsets(64);
  dynamic_array<int> leftPinList(64);

  // ###
  // rightH data_
  // ###

  int numRightVerts = 0;
  int numRightHedges = 0;
  int totRightWt = 0;
  int numRightPins = 0;

  dynamic_array<int> rightVertWt(64);
  dynamic_array<int> rightMapOrig(64);
  dynamic_array<int> rightHedgeWts(64);
  dynamic_array<int> rightHedgeOffsets(64);
  dynamic_array<int> rightPinList(64);

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
      leftVertWt[numLeftVerts] = hVertWt[i];
      leftMapOrig[numLeftVerts] = mapToOrig[i];
      mapFromHtoNewH[i] = numLeftVerts++;
      totLeftWt += hVertWt[i];
    } else {
      rightVertWt[numRightVerts] = hVertWt[i];
      rightMapOrig[numRightVerts] = mapToOrig[i];
      mapFromHtoNewH[i] = numRightVerts++;
      totRightWt += hVertWt[i];
    }
  }

  leftVertWt.resize(numLeftVerts);
  leftMapOrig.resize(numLeftVerts);

  rightVertWt.resize(numRightVerts);
  rightMapOrig.resize(numRightVerts);

  leftH = new serial::hypergraph(leftVertWt, numLeftVerts);
  rightH = new serial::hypergraph(rightVertWt, numRightVerts);

  // ###
  // initialise pin list
  // ###

  leftHedgeOffsets[numLeftHedges] = numLeftPins;
  rightHedgeOffsets[numRightHedges] = numRightPins;

  for (i = 0; i < numHHedges; ++i) {
    endOffset = hHedgeOffsets[i + 1];

    leftHedgeLen = 0;
    rightHedgeLen = 0;

    for (j = hHedgeOffsets[i]; j < endOffset; ++j) {
      v = hPinList[j];

      if (hPartVector[v] == 0) {
        leftPinList[numLeftPins + (leftHedgeLen++)] = mapFromHtoNewH[v];
      } else {
        rightPinList[numRightPins + (rightHedgeLen++)] = mapFromHtoNewH[v];
      }
    }

    if (leftHedgeLen > 1) {
      numLeftPins += leftHedgeLen;
      leftHedgeWts[numLeftHedges++] = hHedgeWt[i];
      leftHedgeOffsets[numLeftHedges] = numLeftPins;
    }

    if (rightHedgeLen > 1) {
      numRightPins += rightHedgeLen;
      rightHedgeWts[numRightHedges++] = hHedgeWt[i];
      rightHedgeOffsets[numRightHedges] = numRightPins;
    }
  }

  leftHedgeWts.resize(numLeftHedges);
  leftHedgeOffsets.resize(numLeftHedges + 1);
  leftPinList.resize(numLeftPins);

  rightHedgeWts.resize(numRightHedges);
  rightHedgeOffsets.resize(numRightHedges + 1);
  rightPinList.resize(numRightPins);

  // ###
  // now init the hypergraphs
  // ###

  leftH->set_number_of_hyperedges(numLeftHedges);
  leftH->set_number_of_pins(numLeftPins);
  leftH->set_total_weight(totLeftWt);
  leftH->set_hyperedge_weights(leftHedgeWts);
  leftH->set_pin_list(leftPinList);
  leftH->set_hyperedge_offsets(leftHedgeOffsets);

  leftH->buildVtoHedges();

  rightH->set_number_of_hyperedges(numRightHedges);
  rightH->set_number_of_pins(numRightPins);
  rightH->set_total_weight(totRightWt);
  rightH->set_hyperedge_weights(rightHedgeWts);
  rightH->set_pin_list(rightPinList);
  rightH->set_hyperedge_offsets(rightHedgeOffsets);

  rightH->buildVtoHedges();

  // ###
  // now init the bisections
  // ###

  l = new bisection(leftH, bisectAgain - 1, bPartID);
  r = new bisection(rightH, bisectAgain - 1,
                    (bPartID | (1 << (log_k_ - bisectAgain))));

  l->setMap(leftMapOrig);
  r->setMap(rightMapOrig);
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

  procs.random_permutation();

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
