#ifndef OPTIONS_HPP_
#define OPTIONS_HPP_
#include <unordered_set>
#include "boost/program_options.hpp"

namespace parkway {
namespace po = boost::program_options;

class options {
 public:
  options(int argc, char *argv[]);

  po::variables_map &map() {
    return variables_;
  }

  template<typename Type>
  inline Type get(const char *key) const {
    return variables_[key].as<Type>();
  }

  inline void set_number_of_processors(int number) const {
    number_of_processors_ = number;
  }

  inline int number_of_processors() const {
    return number_of_processors_;
  }


 private:
  mutable int number_of_processors_;

  po::variables_map variables_;

  po::options_description all_options_;
  po::options_description general_;
  po::options_description coarsening_;
  po::options_description serial_partitioning_;
  po::options_description recursive_bisection_;
  po::options_description hmetis_;
  po::options_description patoh_;
  po::options_description refinement_;

  void add_general_options();
  void add_coarsening_options();
  void add_serial_partitioning_options();
  void add_refinement_options();

  void process_variables();
  void validate_variables() const;
  void load_config_file();

  template<typename Type>
  bool check_less_than(const char *key, Type value) const;

  template<typename Type>
  bool check_less_than_equal(const char *key, Type value) const;

  template<typename Type>
  bool check_greater_than(const char *key, Type value) const;

  template<typename Type>
  bool check_greater_than_equal(const char *key, Type value) const;

  template<typename Type>
  bool check_between(const char *key, Type lower, Type upper) const;

  template<typename Type>
  bool check_in_set(const char *key, std::unordered_set<Type> &&values) const;

  bool check_clash(const char *key1, const char *key2) const;

  void generate_options_file() const;
};


}  // namespace parkway

#endif  // OPTIONS_HPP_
