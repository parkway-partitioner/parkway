#ifndef INTERNAL_BASE_REFINER_HPP_
#define INTERNAL_BASE_REFINER_HPP_

namespace parkway {
namespace base {
namespace ds = parkway::data_structures;

class refiner {
 public:
  refiner() : maximum_part_weight_(0), average_part_weight_(0.0) {
  }

  virtual void display_options() const = 0;

  inline int maximum_part_weight() const {
    return maximum_part_weight_;
  }

  inline void set_maximum_part_weight(int max) {
    maximum_part_weight_ = max;
  }

  inline void set_average_part_weight(double ave) {
    average_part_weight_ = ave;
  }

 protected:
  int maximum_part_weight_;
  double average_part_weight_;
  ds::dynamic_array<int> part_weights_;

};


}  // namespace base
}  // namespace parkway


#endif  // INTERNAL_BASE_REFINER_HPP_
