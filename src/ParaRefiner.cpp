
#ifndef _PARA_REFINER_CPP
#define _PARA_REFINER_CPP

// ### ParaRefiner.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 3/2/2005: Last Modified
//
// ###

#include "ParaRefiner.hpp"
#include <iostream>

ParaRefiner::ParaRefiner(int rank, int nProcs, int nParts, std::ostream &out)
    : loader(rank, nProcs, nParts, out) {
  partWeights.reserve(nParts);
}

ParaRefiner::~ParaRefiner() {}

void ParaRefiner::load(const parallel::hypergraph &h, MPI_Comm comm) {
  int i;
  int ij;
  int vertsPerProc;
  int endOffset;
  int startOffset;
  int hEdgeLen;
  int recvLen;
  int nonLocIndex;

  int numLocalPins;
  int numLocalHedges;

  int *localPins;
  int *localHedgeOffsets;
  int *localHedgeWeights;

  int j;
  int l;
  int proc;
  int locVert;

  dynamic_array<int> sentToProc;
  dynamic_array<int> vDegs;

  numLocalPins = h.number_of_pins();
  numLocalHedges = h.number_of_hyperedges();
  localPins = h.pin_list();
  localHedgeOffsets = h.hyperedge_offsets();
  localHedgeWeights = h.hyperedge_weights();

  local_vertex_weight_ = h.vertex_weight();
  vertex_weights_ = h.vertex_weights();
  match_vector_ = h.match_vector();

  number_of_local_vertices_ = h.number_of_vertices();
  number_of_vertices_ = h.total_number_of_vertices();
  minimum_vertex_index_ = h.minimum_vertex_index();
  maximum_vertex_index_ = minimum_vertex_index_ + number_of_local_vertices_;

  // ###
  // Prepare data_ structures
  // ###

  number_of_allocated_hyperedges_ = 0;
  number_of_hyperedges_ = 0;
  number_of_local_pins_ = 0;
  vertsPerProc = number_of_vertices_ / processors_;

  vertex_to_hyperedges_offset_.reserve(number_of_local_vertices_ + 1);
  sentToProc.reserve(processors_);
  vDegs.reserve(number_of_local_vertices_);

  for (i = 0; i < number_of_local_vertices_; ++i) {
    vertex_to_hyperedges_offset_[i] = 0;
    vDegs[i] = 0;
  }

  if (display_options_ > 1 && rank_ == 0)
    out_stream << "[PR]";

  // ###
  // use the request sets to send local hyperedges to other
  // processors and to receive hyperedges from processors
  // ###

  for (i = 0; i < processors_; ++i) {
    send_lens_[i] = 0;
    sentToProc[i] = 0;
  }

  if (percentile_ == 100) {
    int numActiveProcs;
    int activeProc;

    dynamic_array<int> activeProcs(processors_);

    for (i = 0; i < numLocalHedges; ++i) {
      startOffset = localHedgeOffsets[i];
      endOffset = localHedgeOffsets[i + 1];
      hEdgeLen = endOffset - startOffset;
      numActiveProcs = 0;

      for (j = startOffset; j < endOffset; ++j) {
#ifdef DEBUG_REFINER
        assert(localPins[j] < totalVertices && localPins[j] >= 0);
#endif
        proc = std::min(localPins[j] / vertsPerProc, processors_ - 1);

        if (!sentToProc[proc]) {
          if (proc == rank_) {
            hyperedge_weights_.assign(number_of_hyperedges_, localHedgeWeights[i]);
            hyperedge_offsets_.assign(number_of_hyperedges_++, number_of_local_pins_);

            for (l = startOffset; l < endOffset; ++l) {
              local_pin_list_.assign(number_of_local_pins_++, localPins[l]);
#ifdef DEBUG_REFINER
              assert(localPins[l] < totalVertices && localPins[l] >= 0);
#endif
              if (localPins[l] >= minimum_vertex_index_ &&
                  localPins[l] < maximum_vertex_index_)
                ++vertex_to_hyperedges_offset_[localPins[l] - minimum_vertex_index_];
            }
          } else {
            data_out_sets_[proc]->assign(send_lens_[proc]++, hEdgeLen + 2);
            data_out_sets_[proc]->assign(send_lens_[proc]++, localHedgeWeights[i]);

            for (l = startOffset; l < endOffset; ++l) {
              data_out_sets_[proc]->assign(send_lens_[proc]++, localPins[l]);
            }
          }

          sentToProc[proc] = 1;
          activeProcs[numActiveProcs++] = proc;
        }
      }

      activeProc = activeProcs[RANDOM(0, numActiveProcs)];

#ifdef DEBUG_REFINER
      assert(activeProc >= 0 && activeProc < processors_);
#endif

      if (activeProc == rank_) {
        allocated_hyperedges_.assign(number_of_allocated_hyperedges_++, number_of_hyperedges_ - 1);
      } else {
        data_out_sets_[activeProc]->assign(send_lens_[activeProc]++, RESP_FOR_HEDGE);
      }

      for (j = 0; j < processors_; ++j) {
        sentToProc[j] = 0;
      }
    }
  } else if (percentile_ > 0) {
#ifdef DEBUG_REFINER
    assert(currPercentile > 0 && currPercentile < 100);
#endif

    /* When not loading all the hyperedges, do not need to
       make processors responsible for hyperedges during
       the cutsize calculations                            */

    bit_field toLoad(numLocalHedges);

    compute_hyperedges_to_load(toLoad, numLocalHedges, numLocalPins,
                               localHedgeWeights,
                               localHedgeOffsets, comm);

    for (i = 0; i < numLocalHedges; ++i) {
      if (toLoad(i) == 1) {
        startOffset = localHedgeOffsets[i];
        endOffset = localHedgeOffsets[i + 1];
        hEdgeLen = endOffset - startOffset;

        for (j = startOffset; j < endOffset; ++j) {
#ifdef DEBUG_REFINER
          assert(localPins[j] < totalVertices && localPins[j] >= 0);
#endif
          proc = std::min(localPins[j] / vertsPerProc, processors_ - 1);

          if (!sentToProc[proc]) {
            if (proc == rank_) {
              hyperedge_weights_.assign(number_of_hyperedges_, localHedgeWeights[i]);
              hyperedge_offsets_.assign(number_of_hyperedges_++, number_of_local_pins_);

              for (l = startOffset; l < endOffset; ++l) {
                local_pin_list_.assign(number_of_local_pins_++, localPins[l]);
#ifdef DEBUG_REFINER
                assert(localPins[l] < totalVertices && localPins[l] >= 0);
#endif
                if (localPins[l] >= minimum_vertex_index_ &&
                    localPins[l] < maximum_vertex_index_)
                  ++vertex_to_hyperedges_offset_[localPins[l] - minimum_vertex_index_];
              }
            } else {
              data_out_sets_[proc]->assign(send_lens_[proc]++, hEdgeLen + 2);
              data_out_sets_[proc]->assign(send_lens_[proc]++, localHedgeWeights[i]);

              for (l = startOffset; l < endOffset; ++l) {
                data_out_sets_[proc]->assign(send_lens_[proc]++, localPins[l]);
              }
            }

            sentToProc[proc] = 1;
          }
        }

        for (j = 0; j < processors_; ++j)
          sentToProc[j] = 0;
      }
    }
  } else {
    /* compute a fixed limit on hyperedges to be communicated
       - first try maxLen/2
    */

    int maxHedgeLen = 0;
    int maxTotHedgeLen;
    int limit;

    for (i = 0; i < numLocalHedges; ++i)
      if (localHedgeOffsets[i + 1] - localHedgeOffsets[i] > maxHedgeLen)
        maxHedgeLen = localHedgeOffsets[i + 1] - localHedgeOffsets[i];

    MPI_Allreduce(&maxHedgeLen, &maxTotHedgeLen, 1, MPI_INT, MPI_MAX, comm);

    limit = maxTotHedgeLen / 2;

    for (i = 0; i < numLocalHedges; ++i) {
      startOffset = localHedgeOffsets[i];
      endOffset = localHedgeOffsets[i + 1];
      hEdgeLen = endOffset - startOffset;

      if (hEdgeLen < limit) {
        for (j = startOffset; j < endOffset; ++j) {
#ifdef DEBUG_REFINER
          assert(localPins[j] < totalVertices && localPins[j] >= 0);
#endif
          proc = std::min(localPins[j] / vertsPerProc, processors_ - 1);

          if (!sentToProc[proc]) {
            if (proc == rank_) {
              hyperedge_weights_.assign(number_of_hyperedges_, localHedgeWeights[i]);
              hyperedge_offsets_.assign(number_of_hyperedges_++, number_of_local_pins_);

              for (l = startOffset; l < endOffset; ++l) {
                local_pin_list_.assign(number_of_local_pins_++, localPins[l]);
#ifdef DEBUG_REFINER
                assert(localPins[l] < totalVertices && localPins[l] >= 0);
#endif
                if (localPins[l] >= minimum_vertex_index_ &&
                    localPins[l] < maximum_vertex_index_)
                  ++vertex_to_hyperedges_offset_[localPins[l] - minimum_vertex_index_];
              }
            } else {
              data_out_sets_[proc]->assign(send_lens_[proc]++, hEdgeLen + 2);
              data_out_sets_[proc]->assign(send_lens_[proc]++, localHedgeWeights[i]);

              for (l = startOffset; l < endOffset; ++l) {
                data_out_sets_[proc]->assign(send_lens_[proc]++, localPins[l]);
              }
            }

            sentToProc[proc] = 1;
          }
        }

        for (j = 0; j < processors_; ++j)
          sentToProc[j] = 0;
      }
    }
  }

  // ###
  // Now exchange the hyperedges
  // ###

    send_from_data_out(comm);

  // ###
  // Now load the non-local hyperedges
  // ###

  j = 0;
  for (i = 0; i < processors_; ++i) {
    j += receive_lens_[i];
  }

  recvLen = j;
  j = 0;

  while (j < recvLen) {
#ifdef DEBUG_REFINER
    assert(receive_array_[j] > 0);
#endif
    endOffset = j + receive_array_[j];
    ++j;

    hyperedge_weights_.assign(number_of_hyperedges_, receive_array_[j++]);
    hyperedge_offsets_.assign(number_of_hyperedges_++, number_of_local_pins_);

    for (; j < endOffset; ++j) {
      local_pin_list_.assign(number_of_local_pins_++, receive_array_[j]);

#ifdef DEBUG_REFINER
      assert(receive_array_[j] < totalVertices && receive_array_[j] >= 0);
#endif

      locVert = receive_array_[j] - minimum_vertex_index_;

      if (locVert >= 0 && locVert < number_of_local_vertices_)
        ++vertex_to_hyperedges_offset_[locVert];
    }
#ifdef DEBUG_REFINER
    assert(j == endOffset);
#endif

    if (percentile_ == 100 && j < recvLen &&
        receive_array_[j] == RESP_FOR_HEDGE) {
      allocated_hyperedges_.assign(number_of_allocated_hyperedges_++, number_of_hyperedges_ - 1);
      ++j;
    }
  }

  hyperedge_offsets_.assign(number_of_hyperedges_, number_of_local_pins_);

#ifdef MEM_OPT
  hyperedge_offsets_.reserve(number_of_hyperedges_ + 1);
  hyperedge_weights_.reserve(number_of_hyperedges_);
  allocated_hyperedges_.reserve(number_of_hyperedges_);
  local_pin_list_.reserve(number_of_local_pins_);
#endif

  // ###
  // now initialise the vToHedgesList
  // ###

  j = 0;
  l = 0;

  for (; j < number_of_local_vertices_; ++j) {
#ifdef DEBUG_REFINER
    assert(vToHedgesOffset[j] >= 0);
#endif

    locVert = vertex_to_hyperedges_offset_[j];
    vertex_to_hyperedges_offset_[j] = l;
    l += locVert;
  }

  vertex_to_hyperedges_offset_[j] = l;
  vertex_to_hyperedges_.reserve(l);

  for (j = 0; j < number_of_hyperedges_; ++j) {
    endOffset = hyperedge_offsets_[j + 1];

    for (l = hyperedge_offsets_[j]; l < endOffset; ++l) {
      if (local_pin_list_[l] >= minimum_vertex_index_ && local_pin_list_[l] <
                                                    maximum_vertex_index_) {

        locVert = local_pin_list_[l] - minimum_vertex_index_;
        vertex_to_hyperedges_[vertex_to_hyperedges_offset_[locVert] + (vDegs[locVert]++)] = j;
      }
    }
  }

#ifdef DEBUG_REFINER
  for (j = 1; j <= numLocalVertices; ++j)
    assert(vDegs[j - 1] == vToHedgesOffset[j] - vToHedgesOffset[j - 1]);
#endif

  /* now init the non-local vertex structs */

  numNonLocVerts = 0;

  if (numLocalPins < number_of_vertices_ / 2)
    toNonLocVerts.create(numLocalPins, 1);
  else
    toNonLocVerts.create(number_of_vertices_, 0);

  for (i = 0; i < number_of_hyperedges_; ++i) {
    endOffset = hyperedge_offsets_[i + 1];

    for (j = hyperedge_offsets_[i]; j < endOffset; ++j) {
      ij = local_pin_list_[j];
#ifdef DEBUG_REFINER
      assert(ij >= 0 && ij < totalVertices);
#endif
      if (ij < minimum_vertex_index_ || ij >= maximum_vertex_index_) {
        nonLocIndex = toNonLocVerts.insert_if_empty(ij, numNonLocVerts);

        if (nonLocIndex == -1) {
          vDegs.assign(numNonLocVerts, 1);
          nonLocVerts.assign(numNonLocVerts++, ij);
        } else {
#ifdef DEBUG_REFINER
          assert(nonLocIndex < numNonLocVerts);
#endif
          ++vDegs[nonLocIndex];
        }
      }
    }
  }

  numNonLocVertsHedges = 0;
  for (i = 0; i < numNonLocVerts; ++i) {
#ifdef DEBUG_REFINER
    assert(vDegs[i] >= 1);
#endif
    numNonLocVertsHedges += vDegs[i];
  }

  nonLocVerts.reserve(numNonLocVerts);
  nonLocOffsets.reserve(numNonLocVerts + 1);
  nonLocVToHedges.reserve(numNonLocVertsHedges);

#ifdef DEBUG_REFINER
  for (i = 0; i < numNonLocVertsHedges; ++i)
    nonLocVToHedges[i] = -1;
#endif

  ij = 0;
  for (i = 0; i < numNonLocVerts; ++i) {
    nonLocOffsets[i] = ij;
    ij += vDegs[i];
    vDegs[i] = 0;
  }
  nonLocOffsets[i] = ij;

  /* now intialise the vToHedges for non-local vertices */

  for (i = 0; i < number_of_hyperedges_; ++i) {
    endOffset = hyperedge_offsets_[i + 1];

    for (j = hyperedge_offsets_[i]; j < endOffset; ++j) {
      ij = local_pin_list_[j];
#ifdef DEBUG_REFINER
      assert(ij >= 0 && ij < totalVertices);
#endif
      if (ij < minimum_vertex_index_ || ij >= maximum_vertex_index_) {
        nonLocIndex = toNonLocVerts.get(ij);
#ifdef DEBUG_REFINER
        assert(nonLocIndex >= 0 && nonLocIndex < numNonLocVerts);
#endif
        l = nonLocOffsets[nonLocIndex] + vDegs[nonLocIndex];
        nonLocVToHedges[l] = i;
        ++vDegs[nonLocIndex];
      }
    }
  }

#ifdef DEBUG_REFINER
  for (i = 0; i < numNonLocVertsHedges; ++i)
    assert(nonLocVToHedges[i] > -1);
#endif

  if (display_options_ > 1) {
    int numTotHedgesInGraph;
    int numTotPinsInGraph;

    MPI_Reduce(&numLocalHedges, &numTotHedgesInGraph, 1, MPI_INT, MPI_SUM, 0,
               comm);
    MPI_Reduce(&numLocalPins, &numTotPinsInGraph, 1, MPI_INT, MPI_SUM, 0, comm);

    if (rank_ == 0) {
      out_stream << " " << number_of_vertices_ << " " << numTotHedgesInGraph << " "
                 << numTotPinsInGraph << std::endl;
    }
  }
}

void ParaRefiner::initPartitionStructs(const parallel::hypergraph &h, MPI_Comm comm) {
  load(h, comm);

  // ###
  // init max part weight
  // ###

  int totWt;
  int vPerProc;
  int arraySize;
  int vertex;
  int totalToRecv;
  int endOffset;
  int totalToSend;

  int *array;

  int i;
  int j;
  int ij;

  dynamic_array<int> copyOfSendArray;

  MPI_Allreduce(&local_vertex_weight_, &totWt, 1, MPI_INT, MPI_SUM, comm);

  avePartWt = static_cast<double>(totWt) / number_of_partitions_;
  maxPartWt = static_cast<int>(floor(avePartWt + avePartWt * balConstraint));

  numPartitions = h.number_of_partitions();
  partitionVector = h.partition_vector();
  partitionVectorOffsets = h.partition_offsets();
  partitionCuts = h.partition_cuts();

  vPerProc = number_of_vertices_ / processors_;

#ifdef DEBUG_REFINER
  for (int i = 0; i < partitionVectorOffsets[numPartitions]; ++i)
    assert(partitionVector[i] >= 0 && partitionVector[i] < numParts);
#endif

  partIndices.reserve(numNonLocVerts * numPartitions);
  indexIntoPartIndices.reserve(numPartitions + 1);

  j = 0;
  for (i = 0; i < numPartitions; ++i) {
    indexIntoPartIndices[i] = j;
    j += numNonLocVerts;
  }

  /*
    now communicate the partition vector
    requests for values of non-local vertices
  */

  for (i = 0; i < processors_; ++i)
    send_lens_[i] = 0;

  for (i = 0; i < numNonLocVerts; ++i) {
    j = nonLocVerts[i];
#ifdef DEBUG_REFINER
    assert(j < minVertexIndex || j >= maxVertexIndex);
#endif
    ij = std::min(j / vPerProc, processors_ - 1);
#ifdef DEBUG_REFINER
    assert(ij != rank_);
#endif
    data_out_sets_[ij]->assign(send_lens_[ij]++, j);
  }

  ij = 0;
  for (i = 0; i < processors_; ++i) {
    send_displs_[i] = ij;
    ij += send_lens_[i];
  }

#ifdef DEBUG_REFINER
  assert(ij == numNonLocVerts);
  assert(send_lens_[rank_] == 0);
#endif

  send_array_.reserve(ij);
  copyOfSendArray.reserve(ij);
  arraySize = ij;

  ij = 0;
  for (i = 0; i < processors_; ++i) {
    endOffset = send_lens_[i];
    array = data_out_sets_[i]->data();

    for (j = 0; j < endOffset; ++j) {
#ifdef DEBUG_REFINER
      assert(data_[j] < minVertexIndex || data_[j] >= maxVertexIndex);
#endif
      send_array_[ij] = array[j];
      copyOfSendArray[ij++] = array[j];
    }
  }

  MPI_Alltoall(send_lens_.data(), 1, MPI_INT, receive_lens_.data(), 1, MPI_INT,
               comm);

  ij = 0;
  for (i = 0; i < processors_; ++i) {
    receive_displs_[i] = ij;
    ij += receive_lens_[i];
  }

  totalToRecv = ij;
  receive_array_.reserve(ij);

  MPI_Alltoallv(send_array_.data(), send_lens_.data(),
                send_displs_.data(), MPI_INT, receive_array_.data(),
                receive_lens_.data(), receive_displs_.data(), MPI_INT, comm);

  // ###
  // now communicate the partition vector
  // requests for values of non-local vertices
  // ###

  totalToSend = totalToRecv * numPartitions;
  send_array_.reserve(totalToSend);

  for (i = 0; i < processors_; ++i) {
    ij = send_lens_[i];
    send_lens_[i] = receive_lens_[i] * numPartitions;
    receive_lens_[i] = ij * numPartitions;
  }

  ij = 0;
  j = 0;

  for (i = 0; i < processors_; ++i) {
    send_displs_[i] = ij;
    receive_displs_[i] = j;

    ij += send_lens_[i];
    j += receive_lens_[i];
  }

// ###
// now get the partition vector values
// of the requested local vertices
// ###

#ifdef DEBUG_REFINER
  assert(receive_array_.getLength() == totalToRecv);
  assert(send_array_.getLength() == totalToRecv * numPartitions);
#endif

  ij = 0;
  for (i = 0; i < totalToRecv; ++i) {
    vertex = receive_array_[i] - minimum_vertex_index_;

#ifdef DEBUG_REFINER
    assert(vertex >= 0 && vertex < maxVertexIndex - minVertexIndex);
#endif

    for (j = 0; j < numPartitions; ++j) {
#ifdef DEBUG_REFINER
      int part = partitionVector[partitionVectorOffsets[j] + vertex];
      assert(part >= 0 && part < numParts);
      assert(ij < totalToRecv * numPartitions);
#endif
      send_array_[ij++] = partitionVector[partitionVectorOffsets[j] + vertex];
    }
  }
#ifdef DEBUG_REFINER
  assert(ij == totalToSend);
#endif

  totalToRecv = numNonLocVerts * numPartitions;
  receive_array_.reserve(totalToRecv);

  MPI_Alltoallv(send_array_.data(), send_lens_.data(),
                send_displs_.data(), MPI_INT, receive_array_.data(),
                receive_lens_.data(), receive_displs_.data(), MPI_INT, comm);

  ij = 0;
  for (i = 0; i < arraySize; ++i) {
    vertex = toNonLocVerts.get(copyOfSendArray[i]);
#ifdef DEBUG_REFINER
    assert(vertex >= 0 && vertex < numNonLocVerts);
#endif
    for (j = 0; j < numPartitions; ++j)
      partIndices[indexIntoPartIndices[j] + vertex] = receive_array_[ij++];
  }
}

#endif
