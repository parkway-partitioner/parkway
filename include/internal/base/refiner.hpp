#ifndef INTERNAL_BASE_REFINER_HPP_
#define INTERNAL_BASE_REFINER_HPP_

namespace parkway {
namespace base {

class refiner {
 public:
  refiner() : maximum_part_weight_(0), average_part_weight_(0.0) {
  }

  virtual ~refiner() {
  }

  inline void set_maximum_part_weight(int max) {
    maximum_part_weight_ = max;
  }

  inline int maximum_part_weight() const {
    return maximum_part_weight_;
  }

  inline void set_average_part_weight(double ave) {
    average_part_weight_ = ave;
  }

  virtual void display_options() const = 0;



 protected:
  int maximum_part_weight_;

  double average_part_weight_;
};

}  // namespace base
}  // namespace parkway

#endif  // INTERNAL_BASE_REFINER_HPP_
