#include "data_structures/internal/table_utils.hpp"

namespace parkway {
namespace data_structures {
namespace internal {

complete_binary_tree<int> table_utils::table_size_tree_;
dynamic_array<int> table_utils::scatter_array_;
int table_utils::scatter_size_ = 0;

table_utils::table_utils() {
  int sizes[16];
  int mins[16];

  sizes[0] = 1091;
  sizes[1] = 2113;
  sizes[2] = 4133;
  sizes[3] = 9067;
  sizes[4] = 17093;
  sizes[5] = 37097;
  sizes[6] = 70099;
  sizes[7] = 145109;
  sizes[8] = 300149;
  sizes[9] = 610217;
  sizes[10] = 1290151;
  sizes[11] = 2600177;
  sizes[12] = 6000109;
  sizes[13] = 12500197;
  sizes[14] = 26000111;
  sizes[15] = 64000147;

  mins[0] = 0;
  mins[1] = 256;
  mins[2] = 512;
  mins[3] = 1024;
  mins[4] = 2056;
  mins[5] = 4112;
  mins[6] = 8224;
  mins[7] = 16448;
  mins[8] = 35000;
  mins[9] = 70000;
  mins[10] = 140000;
  mins[11] = 300000;
  mins[12] = 700000;
  mins[13] = 1800000;
  mins[14] = 4000000;
  mins[15] = 12500000;

  table_size_tree_.setup(sizes, mins, 16);
}

void table_utils::set_scatter_array(int size) {
  scatter_size_ = size;
  scatter_array_.reserve(scatter_size_);

  for (std::size_t i = 0; i < scatter_size_; ++i) {
    scatter_array_[i] = i;
  }

  Funct::randomPermutation(scatter_array_.data(), scatter_size_);
}

}  // namespace internal
}  // namespace data_structures
}  // namespace parkway
