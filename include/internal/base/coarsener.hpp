#ifndef COARSENERS_BASE_COARSENER_HPP_
#define COARSENERS_BASE_COARSENER_HPP_
#include <iostream>
#include <string>

namespace parkway {
namespace base {

class coarsener {
 public:
  coarsener(std::string object_name, int vertex_weight = 0,
            int minimum_nodes = 0, double reduction_ratio = 0.0)
      : maximum_vertex_weight_(vertex_weight),
        minimum_number_of_nodes_(minimum_nodes),
        reduction_ratio_(reduction_ratio),
        object_name_(object_name) {
  }

  virtual ~coarsener() {
  }

  inline int maximum_vertex_weight() const {
    return maximum_vertex_weight_;
  }

  inline void set_maximum_vertex_weight(int max_weight) {
    maximum_vertex_weight_ = max_weight;
  }

  inline int minimum_number_of_nodes() const {
    return minimum_number_of_nodes_;
  }

  inline void set_minimum_number_of_nodes(int nodes) {
    minimum_number_of_nodes_ = nodes;
  }

  inline void set_reduction_ratio(double ratio) {
    reduction_ratio_ = ratio;
  }

  inline int reduction_ratio() {
    return reduction_ratio_;
  }

  inline void print_name(std::ostream &out_stream) {
    out_stream << "[" << object_name_ << "]";
  }

 protected:
  int maximum_vertex_weight_;
  int minimum_number_of_nodes_;
  double reduction_ratio_;

  std::string object_name_;
};

}  // namespace base
}  // namespace parkway

#endif  // COARSENERS_BASE_COARSENER_HPP_
