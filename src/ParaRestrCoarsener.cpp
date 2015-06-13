
#ifndef _PARA_RESTR_COARSENER_CPP
#define _PARA_RESTR_COARSENER_CPP

// ### ParaRestrCoarsener.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 4/1/2005: Last Modified
//
// ###

#include "ParaRestrCoarsener.hpp"
#include "data_structures/complete_binary_tree.hpp"
#include <iostream>

namespace ds = parkway::data_structures;

ParaRestrCoarsener::ParaRestrCoarsener(int rank, int nProcs, int nParts,
                                       std::ostream &out)
    : loader(rank, nProcs, nParts, out) {
  totalHypergraphWt = 0;
  maxVertexWt = 0;
  minNodes = 0;
  stopCoarsening = 0;
  clusterIndex = 0;
  totalClusters = 0;
  myMinCluIndex = 0;
  display_options_ = 0;

  reductionRatio = 0;
  balConstraint = 0;
  numPartitions = 0;

  partitionVector = nullptr;
  partitionVectorOffsets = nullptr;
  partitionCuts = nullptr;
  clusterWeights = nullptr;
  pVector = nullptr;
}

ParaRestrCoarsener::~ParaRestrCoarsener() {}

void ParaRestrCoarsener::load(const parallel::hypergraph &h,
                              MPI_Comm comm) {
  int i;
  int endOffset;
  int startOffset;
  int hEdgeLen;
  int recvLen;
  int numLocalPins;
  int numLocalHedges;

  int *localPins;
  int *localHedgeOffsets;
  int *localHedgeWeights;

  int j;
  int l;
  int proc;
  int locVertex;

  dynamic_array<int> procs(processors_);
  dynamic_array<int> vDegs;
  dynamic_array<int> minLocIndices(processors_);
  dynamic_array<int> maxLocIndices(processors_);
  dynamic_array<int> vPerProc(processors_);

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

  // here add partition information

  numPartitions = h.number_of_partitions();
  partitionVector = h.partition_vector();
  partitionVectorOffsets = h.partition_offsets();
  partitionCuts = h.partition_cuts();

#ifdef DEBUG_COARSENER
  assert(numPartitions == 1);
#endif

  // end partition information

  // ###
  // Prepare data_ structures
  // ###

  MPI_Allgather(&minimum_vertex_index_, 1, MPI_INT, minLocIndices.data(), 1,
                MPI_INT, comm);
  MPI_Allgather(&maximum_vertex_index_, 1, MPI_INT, maxLocIndices.data(), 1,
                MPI_INT, comm);

  for (i = 0; i < processors_; ++i)
    procs[i] = i;

  ds::complete_binary_tree <int> vToProc(procs.data(), minLocIndices.data(),
                                  processors_);

  number_of_hyperedges_ = 0;
  number_of_local_pins_ = 0;

  vertex_to_hyperedges_offset_.reserve(number_of_local_vertices_ + 1);
  vDegs.reserve(number_of_local_vertices_);

  for (i = 0; i < number_of_local_vertices_; ++i) {
    vertex_to_hyperedges_offset_[i] = 0;
    vDegs[i] = 0;
  }

  if (display_options_ > 1 && rank_ == 0)
    out_stream << "[PRFCC]";

  // ###
  // use the data_out_sets_ to send local hyperedges to other
  // processors and to receive hyperedges from processors
  // ###

  for (i = 0; i < processors_; ++i) {
    send_lens_[i] = 0;
    vPerProc[i] = 0;
  }

  if (percentile_ == 100) {
    for (i = 0; i < numLocalHedges; ++i) {
      startOffset = localHedgeOffsets[i];
      endOffset = localHedgeOffsets[i + 1];
      hEdgeLen = endOffset - startOffset;

      for (j = startOffset; j < endOffset; ++j) {
#ifdef DEBUG_COARSENER
        assert(localPins[j] < totalVertices && localPins[j] >= 0);
#endif
        proc = vToProc.root_value(localPins[j]);
        ++vPerProc[proc];
      }

      for (j = 0; j < processors_; ++j) {
        if (vPerProc[j] > 1) {
          if (j == rank_) {
            hyperedge_weights_.assign(number_of_hyperedges_, localHedgeWeights[i]);
            hyperedge_offsets_.assign(number_of_hyperedges_++, number_of_local_pins_);

            for (l = startOffset; l < endOffset; ++l) {
#ifdef DEBUG_COARSENER
              assert(localPins[l] < totalVertices && localPins[l] >= 0);
#endif
              locVertex = localPins[l] - minimum_vertex_index_;

              if (locVertex >= 0 && locVertex < number_of_local_vertices_) {
                local_pin_list_.assign(number_of_local_pins_++, locVertex);
                ++vertex_to_hyperedges_offset_[locVertex];
              }
            }
          } else {
            data_out_sets_[j]->assign(send_lens_[j]++, vPerProc[j] + 2);
            data_out_sets_[j]->assign(send_lens_[j]++, localHedgeWeights[i]);

            for (l = startOffset; l < endOffset; ++l) {
              if (localPins[l] >= minLocIndices[j] &&
                  localPins[l] < maxLocIndices[j])
                data_out_sets_[j]->assign(send_lens_[j]++,
                                       localPins[l] - minLocIndices[j]);
            }
          }

          vPerProc[j] = 0;
        }

        if (vPerProc[j] == 1)
          vPerProc[j] = 0;
      }
    }
  } else {
#ifdef DEBUG_COARSENER
    assert(currPercentile > 0 && currPercentile < 100);
#endif

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
#ifdef DEBUG_COARSENER
          assert(localPins[j] < totalVertices && localPins[j] >= 0);
#endif
          proc = vToProc.root_value(localPins[j]);
          ++vPerProc[proc];
        }

        for (j = 0; j < processors_; ++j) {
          if (vPerProc[j] > 1) {
            if (j == rank_) {
              hyperedge_weights_.assign(number_of_hyperedges_, localHedgeWeights[i]);
              hyperedge_offsets_.assign(number_of_hyperedges_++, number_of_local_pins_);

              for (l = startOffset; l < endOffset; ++l) {
#ifdef DEBUG_COARSENER
                assert(localPins[l] < totalVertices && localPins[l] >= 0);
#endif
                locVertex = localPins[l] - minimum_vertex_index_;

                if (locVertex >= 0 && locVertex < number_of_local_vertices_) {
                  local_pin_list_.assign(number_of_local_pins_++, locVertex);
                  ++vertex_to_hyperedges_offset_[locVertex];
                }
              }
            } else {
              data_out_sets_[j]->assign(send_lens_[j]++, vPerProc[j] + 2);
              data_out_sets_[j]->assign(send_lens_[j]++, localHedgeWeights[i]);

              for (l = startOffset; l < endOffset; ++l) {
                if (localPins[l] >= minLocIndices[j] &&
                    localPins[l] < maxLocIndices[j])
                  data_out_sets_[j]->assign(send_lens_[j]++,
                                         localPins[l] - minLocIndices[j]);
              }
            }

            vPerProc[j] = 0;
          }

          if (vPerProc[j] == 1)
            vPerProc[j] = 0;
        }
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
    endOffset = j + receive_array_[j];
#ifdef DEBUG_COARSENER
    assert(recvLen >= endOffset);
#endif
    ++j;

    hyperedge_weights_.assign(number_of_hyperedges_, receive_array_[j++]);
    hyperedge_offsets_.assign(number_of_hyperedges_++, number_of_local_pins_);

    for (; j < endOffset; ++j) {
      locVertex = receive_array_[j];
#ifdef DEBUG_COARSENER
      assert(locVertex < numLocalVertices && locVertex >= 0);
#endif
      local_pin_list_.assign(number_of_local_pins_++, locVertex);
      ++vertex_to_hyperedges_offset_[locVertex];
    }
  }

  hyperedge_offsets_.assign(number_of_hyperedges_, number_of_local_pins_);

#ifdef MEM_OPT
  hyperedge_offsets_.reserve(number_of_hyperedges_ + 1);
  hyperedge_weights_.reserve(number_of_hyperedges_);
  local_pin_list_.reserve(number_of_local_pins_);
#endif

  // ###
  // now initialise the vToHedgesList
  // ###

  j = 0;
  l = 0;

  for (; j < number_of_local_vertices_; ++j) {
#ifdef DEBUG_COARSENER
    assert(vToHedgesOffset[j] >= 0);
#endif

    locVertex = vertex_to_hyperedges_offset_[j];
    vertex_to_hyperedges_offset_[j] = l;
    l += locVertex;
  }

  vertex_to_hyperedges_offset_[j] = l;
  vertex_to_hyperedges_.reserve(l);

  for (j = 0; j < number_of_hyperedges_; ++j) {
    endOffset = hyperedge_offsets_[j + 1];

    for (l = hyperedge_offsets_[j]; l < endOffset; ++l) {
      locVertex = local_pin_list_[l];
      vertex_to_hyperedges_[vertex_to_hyperedges_offset_[locVertex] + (vDegs[locVertex]++)] = j;
    }
  }

#ifdef DEBUG_COARSENER
  for (j = 1; j <= numLocalVertices; ++j)
    assert(vDegs[j - 1] == vToHedgesOffset[j] - vToHedgesOffset[j - 1]);
#endif
}

parallel::hypergraph *ParaRestrCoarsener::contractHyperedges(parallel::hypergraph &h,
                                                       MPI_Comm comm) {
  parallel::hypergraph *coarseGraph =
      new parallel::hypergraph(rank_, processors_, clusterIndex, totalClusters,
                         myMinCluIndex, stopCoarsening, partitionCuts[0],
                         clusterWeights->data(), pVector->data());

  clusterWeights = nullptr;
  pVector = nullptr;

  h.contractRestrHyperedges(*coarseGraph, comm);
    h.set_number_of_partitions(0);

  if (display_options_ > 1) {
    int numTotCoarseVerts = coarseGraph->total_number_of_vertices();
    int numLocHedges = coarseGraph->number_of_hyperedges();
    int numLocPins = coarseGraph->number_of_pins();
    int numTotHedges;
    int numTotPins;

    MPI_Reduce(&numLocHedges, &numTotHedges, 1, MPI_INT, MPI_MAX, 0, comm);
    MPI_Reduce(&numLocPins, &numTotPins, 1, MPI_INT, MPI_MAX, 0, comm);

    if (rank_ == 0) {
      out_stream << numTotCoarseVerts << " " << numTotHedges << " "
                 << numTotPins << " " << std::endl;
    }
  }

  return coarseGraph;
}

#endif
