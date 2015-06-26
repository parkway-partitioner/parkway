#ifndef _FM_REFINER_HPP
#define _FM_REFINER_HPP
// ### FMRefiner.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 30/11/2004: Last Modified
//
// NOTES:
//
// 21/11/2004:
//
// changed so that does not always expect a balanced bisection
// furthermore, should the initial bisection be unbalanced, it
// does not expect to be able to rebalance it.
//
// ###
#include <ostream>
#include "refiner.hpp"
#include "data_structures/bit_field.hpp"
#include "data_structures/bucket_node.hpp"
#include "hypergraph/serial/hypergraph.hpp"

namespace parkway {
namespace serial {

namespace ds = parkway::data_structures;

class fm_refiner : public refiner {
 protected:
  // ###
  // auxiliary members
  // ###

  int bucket_arrays_length_;
  int maximum_possible_gain_;
  int insert_method_;
  int ee_threshold_;
  int maximum_non_positive_moves_;

  // ###
  // bucket arrays
  // ###

  NodePtrArray buckets_;
  ds::dynamic_array<NodeArray *> bucket_arrays_;

  ds::dynamic_array<int> maximum_gains_;
  ds::dynamic_array<int> maximum_gain_entries_;
  ds::dynamic_array<int> number_of_buckets_in_array_;
  ds::dynamic_array<int> vertex_in_part_;
  ds::dynamic_array<int> vertex_gains_;
  ds::dynamic_array<int> move_list_;

  ds::bit_field loaded_;

 public:
  fm_refiner(int max, int insMethod, int ee, int dL);
  ~fm_refiner();

  void display_options() const;

  void build_buckets();
  void restore_buckets();
  void destroy_buckets();
  void initialise_partition_structure();

  void prepare_vertex_gains();
  void initialise_gains_1_to_0();
  void initialise_gains_0_to_1();

  void remove_buckets_from_1();
  void remove_buckets_from_0();
  void remove_unlocked_from_bucket_arrays();
  void move_to_bucket_array(int vPart, int vGain, int v);
  void remove_from_bucket_array(int v, int vPart, int vGain);

  void adjust_maximum_gain_pointer(int dP);
  void update_gains(int v);
  void update_gains_1_to_0(int v);
  void update_gains_0_to_1(int v);
  void undo_move(int v);

  void refine(serial::hypergraph &h);

  int fm_pass();
  int rebalancing_pass(int largePart);
  int choose_maximum_gain_vertex();
  int choose_legal_move(int sP);

  inline void set_ee_threshold() {
    maximum_non_positive_moves_ = static_cast<int>(
        floor((static_cast<double>(numVertices) / 100) * ee_threshold_));
  }

  inline void remove_ee_threshold() { maximum_non_positive_moves_ = numVertices; }

  inline int checker_high(int array) {
    int i;
    for (i = (maximum_possible_gain_ << 1); i > 0; --i) {
      if ((*bucket_arrays_[array])[i].next)
        return i;
    }
    return i;
  }

  inline void sanity_check() {
    for (int i = 0; i < 2; ++i) {
      assert(maximum_gains_[i] == maximum_gain_entries_[i] - maximum_possible_gain_);
      assert(maximum_gain_entries_[i] == checker_high(i));
    }
  }
};

}  // namespace serial
}  // namespace parkway

#endif
