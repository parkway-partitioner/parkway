#ifndef _PARA_GREEDYKWAY_REFINER_HPP
#define _PARA_GREEDYKWAY_REFINER_HPP
// ### ParaGreedyKwayRefiner.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 30/11/2004: Last Modified
//
// ###

#include <iostream>
#include "data_structures/bit_field.hpp"
#include "data_structures/movement_set_table.hpp"
#include "hypergraph/parallel/hypergraph.hpp"
#include "refiners/parallel/refiner.hpp"

namespace parkway {
namespace parallel {
namespace ds = parkway::data_structures;

class k_way_greedy_refiner : public refiner {
 protected:
  int total_number_of_vertices_moved_;
  int early_exit_;

  double limit_;

  // data_ structures from point of view of vertices

  ds::dynamic_array<int> number_of_neighbor_parts_;
  ds::dynamic_array<int> neighbors_of_vertices_;
  ds::dynamic_array<int> neighbors_of_vertices_offsets_;

  // data_ structures from point of view of hyperedges

  ds::dynamic_array<int> hyperedge_vertices_in_part_;
  ds::dynamic_array<int> hyperedge_vertices_in_part_offsets_;

  // auxiliary structures

  ds::dynamic_array<int> vertices_;
  ds::dynamic_array<int> moved_vertices_;
  ds::dynamic_array<int> seen_vertices_;
  ds::dynamic_array<int> number_of_parts_spanned_;
  ds::dynamic_array<int> spanned_parts_;

  ds::bit_field locked_;
  ds::bit_field vertices_seen_;

  // move set structures

  ds::dynamic_array<ds::dynamic_array<int> *> move_sets_;
  ds::dynamic_array<int> move_set_data_;
  ds::dynamic_array<int> index_into_move_set_;
  ds::dynamic_array<int> number_of_vertices_moved_;

  ds::movement_set_table *movement_sets_;

 public:
  k_way_greedy_refiner(int rank, int nProcs, int nParts, int numVperP,
                       int eExit, double lim);
  ~k_way_greedy_refiner();

  void display_options() const;
  void release_memory();
  void initialize_data_structures(const parallel::hypergraph &h, MPI_Comm comm);
  void reset_data_structures();
  void set_partitioning_structures(int pNumber, MPI_Comm comm);
  // void initVertexPartTable(MPI_Comm comm);
  void refine(parallel::hypergraph &h, MPI_Comm comm);

  int greedy_k_way_refinement(parallel::hypergraph &h, int pNo, MPI_Comm comm);
  int greedy_pass(int lowToHigh, MPI_Comm comm);
  int compute_cutsize(MPI_Comm comm);

  void manage_balance_constraint(MPI_Comm comm);
  void undo_pass_moves();

  void update_vertex_move_info(MPI_Comm comm);
  void update_adjacent_vertex_status(int v, int sP, int bestMove);
  void undo_move(int indexIntoMoveSets, int from, int to);

  void sanity_hyperedge_check() const;
  void non_local_vertices_check() const;
};

}  // namespace parallel
}  // namespace parkway

#endif
