#ifndef _PARA_APPROX_COARSENER_CPP
#define _PARA_APPROX_COARSENER_CPP

// ### ParaApproxCoarsener.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 31/12/2004: Last Modified
//
// ###

#include "parallel_approximate_coarsener.hpp"
#include <iostream>

parallel_approximate_coarsener::parallel_approximate_coarsener(int _rank, int processors,
                                         int _numParts, int percentile, int inc,
                                         std::ostream &out)
    : parallel_coarsener(_rank, processors, _numParts, out) {
  startPercentile = percentile;
  currPercentile = percentile;
  increment = inc;
}

parallel_approximate_coarsener::~parallel_approximate_coarsener() {}

void parallel_approximate_coarsener::load(const hypergraph &h, MPI_Comm comm) {
  int i;
  int vertsPerProc;
  int endOffset;
  int startOffset;
  int hEdgeLen;
  int Len;
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

  bit_field toLoad;

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
  toLoad.reserve(numLocalHedges);

  for (i = 0; i < number_of_local_vertices_; ++i) {
    vertex_to_hyperedges_offset_[i] = 0;
    vDegs[i] = 0;
  }

  if (currPercentile < 100)
    compute_hyperedges_to_load(toLoad, numLocalHedges, localHedgeWeights,
                               localHedgeOffsets, comm);
  else
    toLoad.set();

  // ###
  // use the request sets to send local hyperedges to other
  // processors and to receive hyperedges from processors
  // ###

  for (i = 0; i < processors_; ++i) {
    send_lens_[i] = 0;
    sentToProc[i] = 0;
  }

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
  Len = j;

  j = 0;
  while (j < Len) {
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
  local_pin_list_.reserve(number_of_hyperedges_);
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

void parallel_approximate_coarsener::compute_hyperedges_to_load(
    bit_field &toLoad, int numH,
    int *hEdgeWts, int *hEdgeOffsets,
    MPI_Comm comm) {
  int i;
  int j;

  double percentileThreshold;
  double aveLen;

  int totLen;
  int maxLen;
  int numTotHedges;
  int myTotLen = 0;
  int myMaxLen = 0;
  int myPercentileLen;
  int percentileLen;

  dynamic_array<int> hEdges(numH);
  dynamic_array<int> hEdgeLens(numH);

  /* compute the hyperedges that will not be communicated */

  for (i = 0; i < numH; ++i) {
    hEdgeLens[i] = hEdgeOffsets[i + 1] - hEdgeOffsets[i];
    myTotLen += hEdgeLens[i];

    if (hEdgeLens[i] > myMaxLen)
      myMaxLen = hEdgeLens[i];
  }

  MPI_Allreduce(&myTotLen, &totLen, 1, MPI_INT, MPI_SUM, comm);
  MPI_Allreduce(&numH, &numTotHedges, 1, MPI_INT, MPI_SUM, comm);
  MPI_Allreduce(&myMaxLen, &maxLen, 1, MPI_INT, MPI_MAX, comm);

  aveLen = static_cast<double>(totLen) / numTotHedges;

  j = 0;
  for (i = 0; i < numH; ++i) {
    hEdges[i] = i;
    j += hEdgeWts[i];
  }

  percentileThreshold = (static_cast<double>(j) * currPercentile) / 100;
  Funct::qsortByAnotherArray(0, numH - 1, hEdges.data(),
                             hEdgeLens.data(), INC);
  toLoad.set();

  j = 0;
  i = 0;

  for (; i < numH && j < percentileThreshold;)
    j += hEdgeWts[hEdges[i++]];

  myPercentileLen = hEdgeLens[hEdges[i]];

  MPI_Barrier(comm);
  std::cout << "myPercentileLen = " << myPercentileLen << std::endl;
  MPI_Barrier(comm);

  MPI_Allreduce(&myPercentileLen, &percentileLen, 1, MPI_INT, MPI_MAX, comm);

  MPI_Barrier(comm);
  if (rank_ == 0) {
    std::cout << "percentileLen = " << percentileLen << std::endl;
    std::cout << "maxHedgeLen = " << maxLen << std::endl;
    std::cout << "aveLen = " << aveLen << std::endl;
    std::cout << "currPercentile = " << currPercentile << std::endl;
  }

  for (; i < numH; ++i)
    if (hEdgeLens[hEdges[i]] > percentileLen)
      toLoad.unset(hEdges[i]);

  /* increment the percentile for subsequent coarsening */

  currPercentile = std::min(currPercentile + increment, 100);
}

#endif
