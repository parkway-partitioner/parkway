#include "hypergraph/base_hypergraph.hpp"

namespace parkway {
namespace hypergraph {

base_hypergraph::base_hypergraph(int number_of_vertices)
    : number_of_vertices_(number_of_vertices),
      number_of_hyperedges_(0),
      number_of_pins_(0),
      number_of_partitions_(0) {
}

base_hypergraph::~base_hypergraph() {
}


}  // namespace hypergraph
}  // namespace parkway
