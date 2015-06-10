#ifndef _DATA_STRUCTURES_MAP_TO_POS_INT_HPP
#define _DATA_STRUCTURES_MAP_TO_POS_INT_HPP

#include "data_structures/dynamic_array.hpp"

namespace parkway {
namespace data_structures {

/* requires positive int keys */
class map_to_pos_int {
 public:
  map_to_pos_int();
  map_to_pos_int(int _size, bool use_hash);
  ~map_to_pos_int() {}

  void create(int _size, bool use_hash);
  void destroy();
  void recover();

  bool insert(int key, int val);
  void clear();
  bool insert_if_empty(int key, int val, int &already_there);

  int get_careful(int key);
  int get(int key);

  inline bool use_hash() {
    return use_hash_;
  }

  inline int size() {
    return size_;
  }

  inline int capacity() {
    return capacity_;
  }

 protected:
  int size_;
  int capacity_;
  bool use_hash_;

  dynamic_array<int> entries;
  dynamic_array<int> table;
  dynamic_array<int> keys;
};


}  // namespace data_structures
}  // namespace parkway

#endif  // _DATA_STRUCTURES_MAP_TO_POS_INT_HPP
