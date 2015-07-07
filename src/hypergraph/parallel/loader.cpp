// ### ParaHypergraphLoader.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 4/1/2005: Last Modified
//
// ###

#include "hypergraph/parallel/loader.hpp"
#include "utility/logging.hpp"

namespace parkway {
namespace parallel {

loader::loader(int rank, int number_of_processors, int number_of_parts)
    : global_communicator(rank, number_of_processors),
      number_of_parts_(number_of_parts),
      number_of_hyperedges_(0),
      number_of_local_pins_(0),
      number_of_local_vertices_(0),
      number_of_vertices_(0),
      minimum_vertex_index_(0),
      maximum_vertex_index_(0),
      local_vertex_weight_(0),
      number_of_allocated_hyperedges_(0),
      percentile_(100) {
}

loader::~loader() {
}

void loader::compute_hyperedges_to_load(ds::bit_field &to_load,
                                        int number_of_hyperedges,
                                        int num_local_pins,
                                        dynamic_array<int> &hyperedge_weights,
                                        dynamic_array<int> &hyperedge_offsets,
                                        MPI_Comm comm) {

  dynamic_array<int> hyperedges(number_of_hyperedges);
  dynamic_array<int> hyperedge_lengths(number_of_hyperedges);

  /* compute the hyperedges that will not be communicated */
  int max_length_on_this_rank = 0;
  for (int i = 0; i < number_of_hyperedges; ++i) {
    hyperedge_lengths[i] = hyperedge_offsets[i + 1] - hyperedge_offsets[i];
    hyperedges[i] = i;
    if (hyperedge_lengths[i] > max_length_on_this_rank)
      max_length_on_this_rank = hyperedge_lengths[i];
  }

  int data_on_this_rank[2] = {number_of_local_pins_, number_of_hyperedges};
  int hypergraph_data[2];
  int max_length;
  MPI_Allreduce(&data_on_this_rank[0], &hypergraph_data[0], 2, MPI_INT,
                MPI_SUM, comm);
  MPI_Allreduce(&max_length_on_this_rank, &max_length, 1, MPI_INT,
                MPI_MAX, comm);

  double average_length =
      static_cast<double>(hypergraph_data[0]) / hypergraph_data[1];

  int total_weigtht = 0;
  for (const int &weight : hyperedge_weights) {
    total_weigtht += weight;
  }

  double percentile_threshold =
      (static_cast<double>(total_weigtht) * percentile_) / 100;

  hyperedges.sort_using_another_array(hyperedge_lengths);
  to_load.set();

  int i = 0;
  int j = 0;
  while (i < number_of_hyperedges && j < percentile_threshold) {
    j += hyperedge_weights[hyperedges[i++]];
  }

  int percentile_length_this_rank = hyperedge_lengths[hyperedges[i]];
  int percentile_length;
  MPI_Allreduce(&percentile_length_this_rank, &percentile_length, 1, MPI_INT,
                MPI_MAX, comm);

  progress(" %i %i %.2f %i ", percentile_, max_length, average_length,
           percentile_length);

  for (; i < number_of_hyperedges; ++i) {
    if (hyperedge_lengths[hyperedges[i]] > percentile_length) {
      to_load.unset(hyperedges[i]);
    }
  }
}

}  // namespace parallel
}  // namespace parkway
