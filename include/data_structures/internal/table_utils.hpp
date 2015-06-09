#ifndef _DATA_SRUCTURES_INTERNAL_TABLE_UTILS_HPP
#define _DATA_SRUCTURES_INTERNAL_TABLE_UTILS_HPP

#include "data_structures/dynamic_array.hpp"
#include "data_structures/complete_binary_tree.hpp"

namespace parkway {
namespace data_structures {
namespace internal {
namespace hashes {

template<typename Type> Type primary(const Type a, const Type b) {
  return a % b;
}

template<typename Type> Type secondary(const Type a, const Type b) {
  return 1 + (a % (b - 1));
}

template<typename Type> Type chained(const int slot, Type a, Type b) {
  return (slot + secondary(a, b)) % b;
}

}  // namespace hashes

class table_utils {
 protected:
  static complete_binary_tree<int> table_size_tree_;
  static dynamic_array<int> scatter_array_;
  static int scatter_size_;

 public:
  table_utils();
  ~table_utils() {}

  static void set_scatter_array(int size);

  static inline int scatter_key(int i) {
    return scatter_array_[i];
  }

  static inline int table_size(int n) {
    return table_size_tree_.root_value(n);
  }
};

}  // namespace internal
}  // namespace data_structures
}  // namespace parkway

#endif  // _DATA_SRUCTURES_INTERNAL_TABLE_UTILS_HPP
