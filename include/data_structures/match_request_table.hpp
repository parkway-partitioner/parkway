#ifndef _DATA_STRUCTURES_MATCH_REQUEST_TABLE_HPP
#define _DATA_STRUCTURES_MATCH_REQUEST_TABLE_HPP

#include "data_structures/dynamic_array.hpp"

namespace parkway {
namespace data_structures {

class match_request_table {
 public:
  class entry;

  match_request_table(int initial_capacity);
  ~match_request_table();

  int cluster_weight(int _vertex) const;
  int cluster_index(int _vertex) const;
  int local_count(int _vertex) const;

  inline std::size_t size() const {
    return size_;
  }

  inline std::size_t capacity() const {
    return capacity_;
  }

  inline dynamic_array<entry *> get_entries() const {
    return entries_;
  }

  entry *get_entry(int i) const;

  void add_local(int _vertex, int _local, int locWt, int nonLocProc);
  void set_cluster_index(int _vertex, int _index, int _cluWt);
  void remove_local(int _vertex, int _local, int locWt);
  void remove_entry(int _vertex);
  void clear();

 protected:
  std::size_t size_;
  std::size_t capacity_;

  dynamic_array<entry *> table_;
  dynamic_array<entry *> entries_;
};


class match_request_table::entry {
 protected:
  int non_local_vertex_;
  int cluster_weight_;
  int cluster_index_;
  // TODO(gb610): procress or processors?
  int non_local_process_;
  int number_local_;

  dynamic_array<int> local_vertices_;
  entry *next_;

 public:
  inline entry(int _nonlocal, int _local, int _locWt, int _proc, entry *_next)
      : non_local_vertex_(_nonlocal),
        cluster_weight_(_locWt),
        cluster_index_(-1),
        non_local_process_(_proc),
        number_local_(1),
        next_(_next) {
    local_vertices_[0] = _local;
  }

  inline int non_local_vertex() const {
    return non_local_vertex_;
  }

  inline int cluster_weight() const {
    return cluster_weight_;
  }

  inline int cluster_index() const {
    return cluster_index_;
  }

  inline int non_local_process() const {
    return non_local_process_;
  }

  inline int number_local() const {
    return number_local_;
  }

  inline dynamic_array<int> local_vertices_array() const {
    return local_vertices_;
  }

  inline entry *next() const {
    return next_;
  }

  inline void set_next(entry *newNext) {
    next_ = newNext;
  }

  inline void set_cluster_index(int _index) {
    cluster_index_ = _index;
  }

  inline void set_non_local_process(int _proc) {
    non_local_process_ = _proc;
  }

  inline void set_cluster_weight(int _cluWt) {
    cluster_weight_ = _cluWt;
  }

  inline void clear() {
    number_local_ = 0;
    cluster_weight_ = 0;
  }

  inline void add_local(int _loc, int _locWt) {
    local_vertices_[number_local_++] = _loc;
    cluster_weight_ += _locWt;
  }

  inline void remove_local(int loc, int locWt) {
    int i = 0;
    int j = number_local_ - 1;

    for (; i < number_local_; ++i)
      if (local_vertices_[i] == loc)
        break;

    while (i < j) {
      local_vertices_[i] = local_vertices_[i + 1];
      ++i;
    }

    --number_local_;
    cluster_weight_ -= locWt;
  }
};

}  // namespace data_structures
}  // namespace parkway

#endif  // _DATA_STRUCTURES_MATCH_REQUEST_TABLE_HPP
