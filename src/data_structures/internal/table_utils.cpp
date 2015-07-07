#include "data_structures/internal/table_utils.hpp"

namespace parkway {
namespace data_structures {
namespace internal {

const int sizes[16] = {
  1091,
  2113,
  4133,
  9067,
  17093,
  37097,
  70099,
  145109,
  300149,
  610217,
  1290151,
  2600177,
  6000109,
  12500197,
  26000111,
  64000147
};

const int mins[16] = {
  0,
  256,
  512,
  1024,
  2056,
  4112,
  8224,
  16448,
  35000,
  70000,
  140000,
  300000,
  700000,
  1800000,
  4000000,
  12500000
};

complete_binary_tree<int> table_utils::table_size_tree_(sizes, mins, 16);
dynamic_array<int> table_utils::scatter_array_;
const int table_utils::SCATTER_ARRAY_NOT_SET = -1;
int table_utils::scatter_size_ = table_utils::SCATTER_ARRAY_NOT_SET;

void table_utils::set_scatter_array(int size) {
  scatter_size_ = size;
  if (size == SCATTER_ARRAY_NOT_SET) {
    scatter_array_.resize(0);
    return;
  }
  scatter_array_.resize(scatter_size_);

  for (int i = 0; i < scatter_size_; ++i) {
    scatter_array_[i] = i;
  }

  scatter_array_.random_permutation();
}

}  // namespace internal
}  // namespace data_structures
}  // namespace parkway
