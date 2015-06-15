// ### ParaRestrCoarsener.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 4/1/2005: Last Modified
//
// ###
#include "coarseners/parallel/restrictive_coarsening.hpp"
#include "data_structures/complete_binary_tree.hpp"
#include "utility/array.hpp"
#include <iostream>

namespace parkway {
namespace parallel {
namespace ds = parkway::data_structures;
namespace utility = parkway::utility;

restrictive_coarsening::restrictive_coarsening(int rank, int nProcs, int nParts,
                                               std::ostream &out)
    : coarsener(rank, nProcs, nParts, out, "Restrictive Coarsener"),
      minimum_nodes_(0),
      total_clusters_(0),
      partition_vector_(nullptr),
      partition_vector_offsets_(nullptr),
      partition_cuts_(nullptr),
      part_vector_(nullptr) {
}

restrictive_coarsening::~restrictive_coarsening() {}

void restrictive_coarsening::load(const hypergraph &h, MPI_Comm comm) {
  int number_of_local_pins = h.number_of_pins();
  int number_of_local_hedges = h.number_of_hyperedges();
  int *local_pins = h.pin_list();
  int *local_hyperedge_offsets = h.hyperedge_offsets();
  int *local_hyperedge_weights = h.hyperedge_weights();

  update_hypergraph_information(h);

  // Add partition information
  number_of_parts_ = h.number_of_partitions();
  partition_vector_ = h.partition_vector();
  partition_vector_offsets_ = h.partition_offsets();
  partition_cuts_ = h.partition_cuts();

  vertex_to_hyperedges_offset_.reserve(number_of_local_vertices_ + 1);
  utility::set_to_zero<int>(vertex_to_hyperedges_offset_.data(),
                            number_of_local_vertices_ + 1);

  if (display_options_ > 1 && rank_ == 0) {
    print_name(out_stream);
  }

  // ###
  // use the data_out_sets_ to send local hyperedges to other
  // processors and to receive hyperedges from processors
  // ###
  prepare_data_to_send(number_of_local_hedges, number_of_local_pins,
                       local_hyperedge_weights, local_hyperedge_offsets,
                       local_pins, comm);

  // Exchange the hyperedges.
  send_from_data_out(comm);

  // Load the non-local hyperedges.
  load_non_local_hyperedges();

  // Initialise the vertex_to_hyperedges list.
  initialize_vertex_to_hyperedges();
}


hypergraph *restrictive_coarsening::contract_hyperedges(hypergraph &h,
                                                        MPI_Comm comm) {
  hypergraph *coarseGraph =
      new hypergraph(rank_, processors_, cluster_index_, total_clusters_,
                     minimum_cluster_index_, stop_coarsening_,
                     partition_cuts_[0], cluster_weights_.data(),
                     part_vector_->data());

  part_vector_ = nullptr;

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

void restrictive_coarsening::prepare_data_to_send(
    int n_local_hyperedges, int n_local_pins,
    int *local_hyperedge_weights, int *local_hyperedge_offsets,
    int *local_pins, MPI_Comm comm) {

  ds::dynamic_array<int> processors(processors_);
  for (int i = 0; i < processors_; ++i) {
    processors[i] = i;
  }

  ds::dynamic_array<int> min_local_indices(processors_);
  ds::dynamic_array<int> max_local_indices(processors_);

  // Prepare data structures
  MPI_Allgather(&minimum_vertex_index_, 1, MPI_INT, min_local_indices.data(), 1,
                MPI_INT, comm);
  MPI_Allgather(&maximum_vertex_index_, 1, MPI_INT, max_local_indices.data(), 1,
                MPI_INT, comm);

  ds::complete_binary_tree<int> vertex_to_processor(processors.data(),
                                                    min_local_indices.data(),
                                                    processors_);

  ds::dynamic_array<int> vertices_per_processor(processors_);
  utility::set_to_zero<int>(vertices_per_processor.data(), processors_);
  utility::set_to_zero<int>(send_lens_.data(), processors_);

  ds::bit_field to_load(n_local_hyperedges);
  if (percentile_ == 100) {
    to_load.set();
  } else {
    compute_hyperedges_to_load(to_load, n_local_hyperedges,
                               n_local_pins, local_hyperedge_weights,
                               local_hyperedge_offsets, comm);
  }

  ds::dynamic_array<bool> sent_to_processor(processors_);
  utility::set_to<bool>(sent_to_processor.data(), processors_, false);
  utility::set_to_zero<int>(send_lens_.data(), processors_);

  for (int i = 0; i < n_local_hyperedges; ++i) {
    if (to_load.test(i)) {
      int start_offset = local_hyperedge_offsets[i];
      int end_offset = local_hyperedge_offsets[i + 1];
      int hyperedge_length = end_offset - start_offset;

      for (int j = start_offset; j < end_offset; ++j) {
        int proc = vertex_to_processor.root_value(local_pins[j]);
        ++vertices_per_processor[proc];
      }

      for (int j = 0; j < processors_; ++j) {
        if (vertices_per_processor[j] > 1) {
          if (j == rank_) {
            hyperedge_weights_.assign(
                number_of_hyperedges_, local_hyperedge_weights[i]);
            hyperedge_offsets_.assign(
                number_of_hyperedges_, number_of_local_pins_);
            ++number_of_hyperedges_;

            for (int l = start_offset; l < end_offset; ++l) {
              int local_vertex = local_pins[l] - minimum_vertex_index_;
              if (0 <= local_vertex &&
                  local_vertex < number_of_local_vertices_) {
                local_pin_list_.assign(number_of_local_pins_++, local_vertex);
                ++vertex_to_hyperedges_offset_[local_vertex];
              }
            }
          } else {
            data_out_sets_[j]->assign(send_lens_[j]++,
                                      vertices_per_processor[j] + 2);
            data_out_sets_[j]->assign(send_lens_[j]++,
                                      local_hyperedge_weights[i]);
            for (int l = start_offset; l < end_offset; ++l) {
              if (min_local_indices[j] <= local_pins[l] &&
                  local_pins[l] < max_local_indices[j]) {
                data_out_sets_[j]->assign(
                    send_lens_[j]++, local_pins[l] - min_local_indices[j]);
              }
            }
          }
          vertices_per_processor[j] = 0;
        }
        if (vertices_per_processor[j] == 1) {
          vertices_per_processor[j] = 0;
        }
      }
    }
  }
}

void restrictive_coarsening::load_non_local_hyperedges() {
  int receive_length = 0;
  for (int i = 0; i < processors_; ++i) {
    receive_length += receive_lens_[i];
  }

  int i = 0;
  while (i < receive_length) {
    int end_offset = i + receive_array_[i];
    ++i;

    hyperedge_weights_.assign(number_of_hyperedges_, receive_array_[i++]);
    hyperedge_offsets_.assign(number_of_hyperedges_, number_of_local_pins_);
    ++i;
    ++number_of_hyperedges_;

    while (i < end_offset) {
      int local_vertex = receive_array_[i];
      local_pin_list_.assign(number_of_local_pins_, local_vertex);
      ++number_of_local_pins_;
      ++vertex_to_hyperedges_offset_[local_vertex];
      ++i;
    }
  }
  hyperedge_offsets_.assign(number_of_hyperedges_, number_of_local_pins_);
}

void restrictive_coarsening::initialize_vertex_to_hyperedges() {
  ds::dynamic_array<int> vertex_degrees(number_of_local_vertices_);
  utility::set_to_zero<int>(vertex_degrees.data(), number_of_local_vertices_);

  int j = 0;
  for (int i = 0; i < number_of_local_vertices_; ++i) {
    int local_vertex = vertex_to_hyperedges_offset_[i];
    vertex_to_hyperedges_offset_[i] = j;
    j += local_vertex;
  }

  vertex_to_hyperedges_offset_[number_of_local_vertices_] = j;
  vertex_to_hyperedges_.reserve(j);

  for (int k = 0; k < number_of_hyperedges_; ++k) {
    int endOffset = hyperedge_offsets_[k + 1];
    for (int l = hyperedge_offsets_[k]; l < endOffset; ++l) {
      int local_vertex = local_pin_list_[l];
      int offset = vertex_to_hyperedges_offset_[local_vertex] +
                   vertex_degrees[local_vertex]++;
      vertex_to_hyperedges_[offset] = k;
    }
  }
}



}  // namespace parallel
}  // namespace parkway
