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
#include "coarseners/parallel/coarsener.hpp"
#include "utility/array.hpp"
#include <iostream>

namespace parkway {
namespace parallel {
namespace utility = parkway::utility;

coarsener::coarsener(int rank, int number_of_processors, int number_of_parts,
                     std::ostream &out, std::string object_name)
    : loader(rank, number_of_processors, number_of_parts, out, 0),
      parkway::base::coarsener(object_name),
      total_hypergraph_weight_(0),
      stop_coarsening_(0),
      cluster_index_(0),
      total_number_of_clusters_(0),
      minimum_cluster_index_(0),
      balance_constraint_(0) {
}

coarsener::~coarsener() {
}


void coarsener::update_hypergraph_information(const hypergraph &h) {
  local_vertex_weight_ = h.vertex_weight();
  vertex_weights_ = h.vertex_weights();
  match_vector_ = h.match_vector();
  number_of_local_vertices_ = h.number_of_vertices();
  number_of_vertices_ = h.total_number_of_vertices();
  minimum_vertex_index_ = h.minimum_vertex_index();
  maximum_vertex_index_ = minimum_vertex_index_ + number_of_local_vertices_;
  number_of_hyperedges_ = 0;
  number_of_local_pins_ = 0;
}


void coarsener::load(const hypergraph &h, MPI_Comm comm) {
  int number_of_local_pins = h.number_of_pins();
  int number_of_local_hedges = h.number_of_hyperedges();
  ds::dynamic_array<int> local_pins = h.pin_list();
  ds::dynamic_array<int> local_hyperedge_offsets = h.hyperedge_offsets();
  ds::dynamic_array<int> local_hyperedge_weights = h.hyperedge_weights();

  update_hypergraph_information(h);

  // Prepare data structures
  int vertices_per_processor = number_of_vertices_ / processors_;

  vertex_to_hyperedges_offset_.assign(number_of_local_vertices_ + 1, 0);

  if (display_options_ > 1 && rank_ == 0) {
    print_name(out_stream);
  }

  // Use the request sets to send local hyperedges to other processors and to
  // receive hyperedges from processors.
  int limit = INT_MAX;
  bool check_limit = percentile_ <= 0;
  if (check_limit) {
    // Compute a fixed limit on hyperedges to be communicated.
    //  - first try maximum_total_hyperedge_length/2
    int maximum_hyperedge_length = 0;
    for (int i = 0; i < number_of_local_hedges; ++i) {
      int difference = local_hyperedge_offsets[i+1] - local_hyperedge_offsets[i];
      if (difference > maximum_hyperedge_length) {
        maximum_hyperedge_length = difference;
      }
    }

    int maximum_total_hyperedge_length;
    MPI_Allreduce(&maximum_hyperedge_length, &maximum_total_hyperedge_length,
                  1, MPI_INT, MPI_MAX, comm);

    limit = maximum_total_hyperedge_length / 2;
  }

  prepare_data_to_send(number_of_local_hedges, number_of_local_pins,
                       vertices_per_processor, local_hyperedge_weights,
                       local_hyperedge_offsets, local_pins, comm, check_limit,
                       limit);

  // Exchange the hyperedges.
  send_from_data_out(comm);

  // Load the non-local hyperedges.
  load_non_local_hyperedges();

  // Initialise the vertex_to_hyperedges list.
  initialize_vertex_to_hyperedges();
}


hypergraph *coarsener::contract_hyperedges(hypergraph &h, MPI_Comm comm) {
  int clusters_per_proc = total_number_of_clusters_ / processors_;
  int clusters_on_this_rank;

  if (rank_ != processors_ - 1) {
    clusters_on_this_rank = clusters_per_proc;
  } else {
    clusters_on_this_rank = clusters_per_proc +
        (total_number_of_clusters_ % processors_);
  }

  ds::dynamic_array<int> min_cluster_index(processors_);
  ds::dynamic_array<int> max_cluster_index(processors_);
  ds::dynamic_array<int> cluster_weights(clusters_on_this_rank);

  for (int i = 0; i < processors_; ++i) {
    if (i == 0) {
      min_cluster_index[i] = 0;
      max_cluster_index[i] = clusters_per_proc;
    } else {
      min_cluster_index[i] = max_cluster_index[i - 1];
      if (i == processors_ - 1) {
        max_cluster_index[i] = total_number_of_clusters_;
      } else {
        max_cluster_index[i] = min_cluster_index[i] + clusters_per_proc;
      }
    }
  }

  for (int i = 0; i < processors_; ++i) {
    if (i == 0) {
      send_displs_[i] = 0;
    } else {
      send_displs_[i] = send_displs_[i - 1] + send_lens_[i - 1];
    }

    int temp1 = std::max(cluster_index_ + minimum_cluster_index_ -
                         max_cluster_index[i], 0);
    int temp2 = std::max(min_cluster_index[i] - minimum_cluster_index_, 0);
    int temp3 = cluster_index_ - (temp1 + temp2);
    send_lens_[i] = std::max(temp3, 0);
  }

  MPI_Alltoall(send_lens_.data(), 1, MPI_INT, receive_lens_.data(), 1, MPI_INT,
               comm);

  int total_to_recv = 0;
  for (int i = 0; i < processors_; ++i) {
    receive_displs_[i] = total_to_recv;
    total_to_recv += receive_lens_[i];
  }
  MPI_Barrier(comm);

  MPI_Alltoallv(cluster_weights_.data(), send_lens_.data(), send_displs_.data(),
                MPI_INT, cluster_weights.data(), receive_lens_.data(),
                receive_displs_.data(), MPI_INT, comm);

  minimum_cluster_index_ = min_cluster_index[rank_];

  parallel::hypergraph *coarseGraph =
      new parallel::hypergraph(rank_,
                               processors_,
                               clusters_on_this_rank,
                               total_number_of_clusters_,
                               minimum_cluster_index_,
                               stop_coarsening_,
                               cluster_weights,
                               h.display_option());

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

bool coarsener::within_vertex_index_range(int value) const {
  return minimum_vertex_index_ <= value && value < maximum_vertex_index_;
}

void coarsener::initialize_vertex_to_hyperedges() {
  ds::dynamic_array<int> vertex_degrees(number_of_local_vertices_, 0);

  int i = 0;
  int l = 0;
  while (i < number_of_local_vertices_) {
    int local_vertex = vertex_to_hyperedges_offset_[i];
    vertex_to_hyperedges_offset_[i] = l;
    l += local_vertex;
    ++i;
  }
  vertex_to_hyperedges_offset_[i] = l;
  vertex_to_hyperedges_.resize(l);


  for (int j = 0; j < number_of_hyperedges_; ++j) {
    int end_offset = hyperedge_offsets_[j + 1];
    for (int l = hyperedge_offsets_[j]; l < end_offset; ++l) {
      if (within_vertex_index_range(local_pin_list_[l])) {
        int local_vertex = local_pin_list_[l] - minimum_vertex_index_;
        int offset = vertex_to_hyperedges_offset_[local_vertex] +
                     vertex_degrees[local_vertex]++;
        vertex_to_hyperedges_[offset] = j;
      }
    }
  }
}


void coarsener::load_non_local_hyperedges() {
  int receive_length = 0;
  for (int i = 0; i < processors_; ++i) {
   receive_length += receive_lens_[i];
  }

  int i = 0;
  while (i < receive_length) {
    int end_offset = i + receive_array_[i];
    ++i;

    hyperedge_weights_[number_of_hyperedges_] = receive_array_[i++];
    hyperedge_offsets_[number_of_hyperedges_++] = number_of_local_pins_;

    for (; i < end_offset; ++i) {
      local_pin_list_[number_of_local_pins_++] = receive_array_[i];
      int local_vertex = receive_array_[i] - minimum_vertex_index_;
      if (0 <= local_vertex && local_vertex < number_of_local_vertices_) {
        ++vertex_to_hyperedges_offset_[local_vertex];
      }
    }
  }
  hyperedge_offsets_[number_of_hyperedges_] = number_of_local_pins_;
}


void coarsener::prepare_data_to_send(
    int n_local_hyperedges, int n_local_pins, int vertices_per_processor,
    dynamic_array<int> &local_hyperedge_weights,
    dynamic_array<int> &local_hyperedge_offsets,
    dynamic_array<int> &local_pins,
    MPI_Comm comm, bool check_limit, int limit) {
  ds::bit_field to_load(n_local_hyperedges);
  if (percentile_ == 100 || check_limit) {
    to_load.set();
  } else {
    compute_hyperedges_to_load(to_load, n_local_hyperedges,
                               n_local_pins, local_hyperedge_weights,
                               local_hyperedge_offsets, comm);
  }

  ds::bit_field sent_to_processor(processors_);
  sent_to_processor.unset();
  send_lens_.assign(processors_, 0);

  for (int i = 0; i < n_local_hyperedges; ++i) {
    if (to_load.test(i)) {
      int start_offset = local_hyperedge_offsets[i];
      int end_offset = local_hyperedge_offsets[i + 1];
      int hyperedge_length = end_offset - start_offset;

      if (!check_limit || (check_limit && hyperedge_length < limit)) {
        for (int j = start_offset; j < end_offset; ++j) {
          int proc = std::min(local_pins[j] / vertices_per_processor, processors_ - 1);

          if (!sent_to_processor[proc]) {
            if (proc == rank_) {
              hyperedge_weights_[number_of_hyperedges_] = local_hyperedge_weights[i];
              hyperedge_offsets_[number_of_hyperedges_++] = number_of_local_pins_;

              for (int l = start_offset; l < end_offset; ++l) {
                local_pin_list_[number_of_local_pins_++] = local_pins[l];
                if (within_vertex_index_range(local_pins[l])) {
                  int index = local_pins[l] - minimum_vertex_index_;
                  ++vertex_to_hyperedges_offset_[index];
                }
              }
            } else {
              data_out_sets_[proc][send_lens_[proc]++] = hyperedge_length + 2;
              data_out_sets_[proc][send_lens_[proc]++] = local_hyperedge_weights[i];
              for (int l = start_offset; l < end_offset; ++l) {
                data_out_sets_[proc][send_lens_[proc]++] = local_pins[l];
              }
            }
            sent_to_processor.set(proc);
          }
        }
      }
      sent_to_processor.unset();
    }
  }
}

}  // namespace parallel
}  // namespace parkway
