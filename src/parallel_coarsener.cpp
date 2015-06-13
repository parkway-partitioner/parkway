
#ifndef _PARA_COARSENER_CPP
#define _PARA_COARSENER_CPP

// ### ParaCoarsener.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 15/4/2004:  Optimisation: changed the vertices
//             and hedges structures into pin-list
//
// 25/4/2004:  Introduced base class ParaHypergraphLoader
//
// 4/1/2005: Last Modified
//
// ###

#include "parallel_coarsener.hpp"
#include <iostream>

parallel_coarsener::parallel_coarsener(int rank, int nProcs, int nParts, std::ostream &out)
    : loader(rank, nProcs, nParts, out) {
  total_hypergraph_weight_ = 0;
  maximum_vertex_weight_ = 0;
  minimum_number_of_nodes_ = 0;
  stop_coarsening_ = 0;
  cluster_index_ = 0;
  total_number_of_clusters_ = 0;
  minimum_cluster_index_ = 0;
  display_options_ = 0;
  reduction_ratio_ = 0;
  balance_constrain_ = 0;

  cluster_weights_.reserve(0);
}

parallel_coarsener::~parallel_coarsener() {}

void parallel_coarsener::load(const parallel::hypergraph &h, MPI_Comm comm) {
  int i;
  int vertsPerProc;
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
    out_stream << "[PFCC]";

  // ###
  // use the request sets to send local hyperedges to other
  // processors and to receive hyperedges from processors
  // ###

  for (i = 0; i < processors_; ++i) {
    send_lens_[i] = 0;
    sentToProc[i] = 0;
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
        proc = std::min(localPins[j] / vertsPerProc, processors_ - 1);

        if (!sentToProc[proc]) {
          if (proc == rank_) {
            hyperedge_weights_.assign(number_of_hyperedges_, localHedgeWeights[i]);
            hyperedge_offsets_.assign(number_of_hyperedges_++, number_of_local_pins_);

            for (l = startOffset; l < endOffset; ++l) {
              local_pin_list_.assign(number_of_local_pins_++, localPins[l]);
#ifdef DEBUG_COARSENER
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
  } else if (percentile_ > 0) {
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
          proc = std::min(localPins[j] / vertsPerProc, processors_ - 1);

          if (!sentToProc[proc]) {
            if (proc == rank_) {
              hyperedge_weights_.assign(number_of_hyperedges_, localHedgeWeights[i]);
              hyperedge_offsets_.assign(number_of_hyperedges_++, number_of_local_pins_);

              for (l = startOffset; l < endOffset; ++l) {
                local_pin_list_.assign(number_of_local_pins_++, localPins[l]);
#ifdef DEBUG_COARSENER
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
#ifdef DEBUG_COARSENER
          assert(localPins[j] < totalVertices && localPins[j] >= 0);
#endif
          proc = std::min(localPins[j] / vertsPerProc, processors_ - 1);

          if (!sentToProc[proc]) {
            if (proc == rank_) {
              hyperedge_weights_.assign(number_of_hyperedges_, localHedgeWeights[i]);
              hyperedge_offsets_.assign(number_of_hyperedges_++, number_of_local_pins_);

              for (l = startOffset; l < endOffset; ++l) {
                local_pin_list_.assign(number_of_local_pins_++, localPins[l]);
#ifdef DEBUG_COARSENER
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
    endOffset = j + receive_array_[j];
    ++j;

    hyperedge_weights_.assign(number_of_hyperedges_, receive_array_[j++]);
    hyperedge_offsets_.assign(number_of_hyperedges_++, number_of_local_pins_);

    for (; j < endOffset; ++j) {
      local_pin_list_.assign(number_of_local_pins_++, receive_array_[j]);

#ifdef DEBUG_COARSENER
      assert(receive_array_[j] < totalVertices && receive_array_[j] >= 0);
#endif

      locVert = receive_array_[j] - minimum_vertex_index_;

      if (locVert >= 0 && locVert < number_of_local_vertices_)
        ++vertex_to_hyperedges_offset_[locVert];
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

#ifdef DEBUG_COARSENER
  for (j = 1; j <= numLocalVertices; ++j)
    assert(vDegs[j - 1] == vToHedgesOffset[j] - vToHedgesOffset[j - 1]);
#endif
}

parallel::hypergraph *parallel_coarsener::constract_hyperedges(
    parallel::hypergraph &h,
    MPI_Comm comm) {
  int i;

  int totalToRecv = 0;
  int numMyClusters;
  int clustersPerProc = total_number_of_clusters_ / processors_;

  if (rank_ != processors_ - 1)
    numMyClusters = clustersPerProc;
  else
    numMyClusters = clustersPerProc + Mod(total_number_of_clusters_, processors_);

  dynamic_array<int> minClusterIndex(processors_);
  dynamic_array<int> maxClusterIndex(processors_);
  dynamic_array<int> *clusterWts = new dynamic_array<int>(numMyClusters);

  for (i = 0; i < processors_; ++i) {
    if (i == 0) {
      minClusterIndex[i] = 0;
      maxClusterIndex[i] = clustersPerProc;
    } else {
      minClusterIndex[i] = maxClusterIndex[i - 1];
      if (i == processors_ - 1)
        maxClusterIndex[i] = total_number_of_clusters_;
      else
        maxClusterIndex[i] = minClusterIndex[i] + clustersPerProc;
    }
  }

  for (i = 0; i < processors_; ++i) {
    if (i == 0)
      send_displs_[i] = 0;
    else
      send_displs_[i] = send_displs_[i - 1] + send_lens_[i - 1];

    send_lens_[i] =
        std::max(cluster_index_ -
                (std::max(cluster_index_ + minimum_cluster_index_ - maxClusterIndex[i], 0) +
                 std::max(minClusterIndex[i] - minimum_cluster_index_, 0)),
            0);
  }

  MPI_Alltoall(send_lens_.data(), 1, MPI_INT, receive_lens_.data(), 1, MPI_INT,
               comm);

  for (i = 0; i < processors_; ++i) {
    receive_displs_[i] = totalToRecv;
    totalToRecv += receive_lens_[i];
  }

#ifdef DEBUG_COARSENER
  assert(totalToRecv == numMyClusters);
#endif

  MPI_Alltoallv(cluster_weights_.data(), send_lens_.data(),
                send_displs_.data(), MPI_INT, clusterWts->data(),
                receive_lens_.data(), receive_displs_.data(), MPI_INT, comm);

  minimum_cluster_index_ = minClusterIndex[rank_];

  parallel::hypergraph *coarseGraph = new parallel::hypergraph(rank_, processors_, numMyClusters,
                                           total_number_of_clusters_,
                                           minimum_cluster_index_,
                                           stop_coarsening_,
                                           clusterWts->data());

    h.contract_hyperedges(*coarseGraph, comm);

  if (display_options_ > 1) {
    int numTotCoarseVerts = coarseGraph->total_number_of_vertices();
    int numLocHedges = coarseGraph->number_of_hyperedges();
    int numLocPins = coarseGraph->number_of_pins();
    int numTotCoarseHedges;
    int numTotCoarsePins;

    MPI_Reduce(&numLocHedges, &numTotCoarseHedges, 1, MPI_INT, MPI_SUM, 0,
               comm);
    MPI_Reduce(&numLocPins, &numTotCoarsePins, 1, MPI_INT, MPI_SUM, 0, comm);

    if (rank_ == 0) {
      out_stream << numTotCoarseVerts << " " << numTotCoarseHedges << " "
                 << numTotCoarsePins << " " << std::endl;
    }
  }

  return coarseGraph;
}

#endif
