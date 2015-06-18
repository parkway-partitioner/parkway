#ifndef _PARA_HYPERGRAPH_HPP
#define _PARA_HYPERGRAPH_HPP
// ### ParaHypergraph.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 17/4/2004: Modified hyperedge storage:
//            - removed hyperedges from hash table
//            - instead, added ds::dynamic_array<int>
//              to represent them as a pin list
//              can be indexed via hash table with key
//              or directly via index in pin list
//
// 03/12/2004: Last Modified
//
// ###

#include <unistd.h>
#include <fstream>
#include <iostream>
#include <cstdio>
#include <cassert>

#include "internal/global_communicator.hpp"
#include "data_structures/dynamic_array.hpp"
#include "data_structures/complete_binary_tree.hpp"
#include "hypergraph/base_hypergraph.hpp"

namespace parkway {
namespace parallel {
namespace ds = parkway::data_structures;
using parkway::hypergraph::base_hypergraph;

class hypergraph : public parkway::global_communicator, public base_hypergraph {
public:
  hypergraph(int rank, int number_of_processors, int number_of_local_vertices,
             int total_vertices, int minimum_vertex_index, int coarsen,
             ds::dynamic_array<int> weight, int display_option);

  hypergraph(int rank, int number_of_processors, int number_of_local_vertices,
             int total_vertices, int minimum_vertex_index, int coarsen,
             int cut, ds::dynamic_array<int> weight,
             ds::dynamic_array<int> part_array, int display_option);

  hypergraph(int rank, int number_of_processors, const char *filename,
             int display_option, std::ostream &out, MPI_Comm comm);

  hypergraph(int rank, int number_of_processors, int number_of_local_vertices,
             int number_of_local_hedges, int max_hyperedge_length,
             ds::dynamic_array<int> vertex_weights,
             ds::dynamic_array<int> hyperedge_weights,
             ds::dynamic_array<int> loc_pin_list,
             ds::dynamic_array<int> hyperedge_offsets,
             int display_option, std::ostream &out, MPI_Comm comm);

  ~hypergraph();

  void load_from_file(const char *filename, std::ostream &out, MPI_Comm comm);

  void initalize_partition_from_file(const char *filename, int numParts,
                                     std::ostream &out, MPI_Comm comm);

  void allocate_hyperedge_memory(int numHedges, int numLocPins);
  void contract_hyperedges(hypergraph &coarse, MPI_Comm comm);
  void contractRestrHyperedges(hypergraph &coarse, MPI_Comm comm);
  void project_partitions(hypergraph &coarse, MPI_Comm comm);
  void reset_vectors();

  void remove_bad_partitions(double cutThreshold);
  void set_number_of_partitions(int nP) override;
  void compute_partition_characteristics(int pNum, int numParts, double constraint,
                                         std::ostream &out, MPI_Comm comm);
  void copy_in_partition(const int *partition, int numV, int nP);
  void copy_out_partition(int *partition, int numV, int nP) const;
  int keep_best_partition();

  void prescribed_vertex_shuffle(int *prescribedAssignment, int nLocVer,
                                 MPI_Comm comm);
  void prescribed_vertex_shuffle(int *mapToOrigV, int *partitionFile,
                                 MPI_Comm comm);
  void shuffle_vertices_by_partition(int nParts, MPI_Comm comm);
  void shuffle_vertices_randomly(MPI_Comm comm);
  void shuffle_vertices_randomly(int *mapToOrigV, MPI_Comm comm);
  void shuffle_vertices_randomly(hypergraph &fineG, MPI_Comm comm);

  void shuffle_vertices(int *vToProc, int *locVPerProc, MPI_Comm comm);
  void shuffleVerticesAftRandom(int *vToProc, int *locVPerProc, int *mapToOrigV,
                                MPI_Comm comm);
  void shuffleVerticesAftRandom(int *vToProc, int *locVPerProc,
                                hypergraph &fineG, MPI_Comm comm);
  void shuffleVerticesAftRandom(int *vToProc, int *locVPerProc,
                                int *mapToInterV, int *mapToOrigV,
                                MPI_Comm comm);

  void shift_vertices_to_balance(MPI_Comm comm);

  int calculate_cut_size(int numParts, int pNum, MPI_Comm comm);
  int check_balance(int numParts, double balConstraint, int numPartition,
                    MPI_Comm comm);

  void check_validity_of_partitions(int numParts) const;
  void check_partitions(int numParts, int maxPartWt, MPI_Comm comm);
  void check_partitions(int numParts, double constraint, std::ostream &out,
                        MPI_Comm comm);
  void compute_balance_warnings(int numParts, double constraint, std::ostream &out,
                                MPI_Comm comm);

  int total_number_of_pins(MPI_Comm comm);
  int total_number_of_hyperedges(MPI_Comm comm);
  int exposed_hyperedge_weight(MPI_Comm comm) const;

  double average_vertex_degree(MPI_Comm comm);
  double average_hyperedge_size(MPI_Comm comm);

  inline int total_number_of_vertices() const { return total_number_of_vertices_; }
  inline int minimum_vertex_index() const { return minimum_vertex_index_; }
  inline int vertex_weight() const { return vertex_weight_; }
  inline int dont_coarsen() const { return do_not_coarsen; }

  inline dynamic_array<int> to_origin_vertex() const {
    return to_origin_vertex_;
  }


  inline void set_total_number_of_vertics(int v) { total_number_of_vertices_ = v; }
  inline void set_minimum_vertex_index(int m) { minimum_vertex_index_ = m; }

  inline void set_cut(int pNo, int cut) {
    partition_cuts_[pNo] = cut;
  }

  inline int cut(int i) const {
    return partition_cuts_[i];
  }

 protected:
  int do_not_coarsen;
  int total_number_of_vertices_;
  int minimum_vertex_index_;
  int vertex_weight_;

  parkway::data_structures::dynamic_array<int> to_origin_vertex_;

  void check_vertex_and_hyperedge_lengths(
      int hyperedge_data_length, const dynamic_array<int> &hyperedge_data,
      MPI_Comm comm);

  void load_data_from_blocks(const int data_length,
                             const dynamic_array<int> &hypergraph_data);

  void check_loaded_vertex_and_hyperedge_lengths() const;

  int compute_number_of_elements_to_send(dynamic_array<int> &copy_of_requests);
  int compute_number_of_elements_to_receive();

  int get_processor(int vertex, int vertices_per_processor = 0,
                    ds::complete_binary_tree<int> *vertex_to_proc = nullptr);

  int compute_requests_for_remote_vertex_matches(
      int vertices_per_processor,
      ds::complete_binary_tree<int> *vertex_to_proc = nullptr);

  void choose_non_local_vertices_format(
      ds::dynamic_array<int> &original_contracted_pin_list);

  void send_coarse_hyperedges(
      ds::dynamic_array<int> &original_contracted_pin_list,
      int &total_to_send, int &total_to_receive);

  void process_new_hyperedges(hypergraph &coarse,
                              ds::new_hyperedge_index_table &table);

};

}  // namespace parallel
}  // namespace parkway

#endif
