#include "options.hpp"
#include "configuration.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <climits>
#include <sys/ioctl.h>

namespace window {
static inline int get_columns() {
  struct winsize w;
  ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);
  return w.ws_col;
}
}

namespace parkway {

options::options(int argc, char *argv[])
  : number_of_processors_(0),
    all_options_("Parkway options", window::get_columns()),
    general_("General options", window::get_columns()),
    coarsening_("Coarsening options", window::get_columns()),
    serial_partitioning_("Serial partitioning options", window::get_columns()),
    recursive_bisection_("Recursive bisection options [only used if not using "
                         "either '--use-hmetis' or '--use-patoh']",
                         window::get_columns()),
    patoh_("PaToH options", window::get_columns()),
    hmetis_("hMETIS options", window::get_columns()),
    refinement_("Refinement options", window::get_columns()) {
  add_general_options();
  add_coarsening_options();
  add_serial_partitioning_options();
  add_refinement_options();

  all_options_.add(general_)
              .add(coarsening_)
              .add(serial_partitioning_)
              .add(recursive_bisection_)
              .add(patoh_)
              .add(hmetis_)
              .add(refinement_);

  po::store(po::parse_command_line(argc, argv, all_options_), variables_);
  process_variables();
  validate_variables();
}

void options::add_general_options() {
  general_.add_options()
    ("help,h", "Prints this help message.")

    ("version,v", "Prints the version of Parkway.")

    ("generate-config-file,g", "Generates a configuartion file.")

    ("display-all-options,d", "Prints help message for all options.")

    ("hypergraph", po::value<std::string>(), "Hypergraph file to partition.")

    ("config,o", po::value<std::string>(), "Configuration file to load.")

    ("number-of-runs,n", po::value<int>()->default_value(1),
     "Number of parallel partitioning runs.")

    ("number-of-parts,p", po::value<int>()->default_value(4),
     "Number of parts sought in partition.")

    // Add the "0.05" otherwise boost performs a lexical_cast with (far) too
    // many decimal places.
    ("balance-constraint,c", po::value<double>()->default_value(0.05, "0.05"),
     "Balance constraint.")

    ("vertex-to-processor-allocation", po::value<int>()->default_value(0),
     "Vertex to processor allocation. Options:\n"
     "  0: as read in,\n"
     "  1: random,\n"
     "  2: as prescribed in partition file.")

    ("sprng-seed", po::value<int>()->default_value(1),
     "Seed for pseudo-random number generator (if compiled with SPRNG).")

    ("display-info", po::bool_switch()->default_value(false),
     "Display partitioning information.")

    ("display-progress", po::bool_switch()->default_value(false),
     "Display partitioning progress.")

    ("output-to-file", po::bool_switch()->default_value(false),
     "Write output to file (otherwise terminal). Select output file with "
     "'--output-file' option.")

    ("output-file",
     po::value<std::string>()->default_value("parkway_output.log"),
     "Output file for info and progress.")

    ("write-partitions-to-file", po::bool_switch()->default_value(false),
     "Write the computed partition to file.")

    ("random-vertex-shuffle", po::bool_switch()->default_value(false),
     "Randomly shuffle vertices between processes before coarsening.")
  ;
}

void options::add_coarsening_options() {
  coarsening_.add_options()
    ("coarsening.type", po::value<std::string>()->default_value("first-choice"),
     "Type of coarsener: Options:\n"
     "  first-choice,\n"
     "  model-2d.")

    ("coarsening.percentile-cutoff", po::value<int>()->default_value(100),
     "Percentile of hyperedge weight l such that the given percent of the "
     "total percent of the total hyperedge weight in the hypergraph is found "
     "in hyperedges with length less than or equal to l.")

    ("coarsening.percentile-increment", po::value<int>()->default_value(1),
     "Increment for the 'percentile-cutoff' option used in subsequent "
     "coarsening steps assuming the required percentile does not exceed 100.")

    ("coarsening.minimum-coarse-vertices", po::value<int>()->default_value(200),
     "Minimum number of vertices in coarsest hypergraph for the parallel "
     "coarsening algorithm.")

    ("coarsening.reduction-ratio", po::value<double>()->default_value(1.75),
     "Reduction ratio during parallel coarsening. Typically between 1.5 and "
     "2.0.")

    ("coarsening.vertex-visit-order",
     po::value<std::string>()->default_value("random"),
     "Local vertex visit order during parallel coarsening phase. Options:\n"
     "  random: random order,\n"
     "  increasing: increasing order by index,\n"
     "  decreasing: decreasing order by index,\n"
     "  increasing-weight: increasing order by weight,\n"
     "  decreasing-weight: decreasing order by weight.\n")

    ("coarsening.connectivity-metric", po::value<int>()->default_value(0),
     "Connectivity metric. Options:\n"
     "  0: cluster weight and hyperedge length not taken into account,\n"
     "  1: inversely proportional to resulting cluster weight only,\n"
     "  2: inversely proportional to lengths of connecting hyperedges only,\n"
     "  3: inversely proportional to both resulting cluster weight and "
     "connecting hyperedge lenghts.")

    ("coarsening.randomly-process-match-requests",
     po::bool_switch()->default_value(false),
     "Process match requests in a random order from other processes.")
  ;
}

void options::add_serial_partitioning_options() {
  serial_partitioning_.add_options()
    ("serial-partitioning.number-of-runs", po::value<int>()->default_value(1),
     "Number of serial partitioning runs. If hMETIS or PaToH are used, each "
     "process will compute max(number-of-serial-runs/p, 1) runs where p is "
     "the number of processes.")

    ("serial-partitioning.util-fan-out",
     po::bool_switch()->default_value(true),
     "Fan out during serial coarsening.")

    ("serial-partitioning.divide-by-cluster-weight",
     po::bool_switch()->default_value(true),
     "Divide by cluster weight during serial coarsening.")

    ("serial-partitioning.use-hmetis", po::bool_switch()->default_value(false),
     "Use hMETIS to perform serial partitioning (if available).")

    ("serial-partitioning.use-patoh", po::bool_switch()->default_value(false),
     "Use PaToH to perform serial partitioning (if available).")
  ;

  recursive_bisection_.add_options()
    ("recursive-bisection.v-cycles",
     po::value<std::string>()->default_value("final-only"),
     "V-cycle setting. Options:\n"
     "  final-only: only use V-cycle refinement in the final bisection,\n"
     "  all: use V-cycle refinement in each recursive bisection,\n"
     "  off: don't use V-cycle refinement.")

    ("recursive-bisection.use-weight-metric",
     po::bool_switch()->default_value(false),
     "Use a metric which is inversely proportional to the resulting cluster "
     "weight.")

    ("recursive-bisection.number-of-runs", po::value<int>()->default_value(2),
     "Number of runs used during bisection.")

    ("recursive-bisection.number-of-initial-partitioning-runs",
     po::value<int>()->default_value(10),
     "Number of bisection runs used on the coarsest hypergraph.")
  ;

  patoh_.add_options()
    ("patoh.partition-options",
     po::value<std::string>()->default_value("default"),
     "PaToH settings. Options:\n"
     "  default: use PATOH_SUGPARAM_DEFAULT\n"
     "  speed: use PATOH_SUGPARAM_SPEED\n"
     "  quality: use PATOH_SUGPARAM_QUALITY")
  ;

  hmetis_.add_options()
    ("hmetis.coarsening",
     po::value<std::string>()->default_value("first-choice"),
     "Sets the coarsening scheme used by hMETIS. Options:\n"
     "  first-choice\n"
     "  hybrid: hybrid first choice\n"
     "  greedy: greedy first choice\n"
     "  hyperedge\n"
     "  edge")

    ("hmetis.refinement",
     po::value<std::string>()->default_value("best-intermediate"),
     "Sets the V-cycle refinement scheme used by hMETIS. Options:\n"
     "  best-intermediate\n"
     "  each-intermediate\n"
     "  final\n"
     "  off")
  ;
}

void options::add_refinement_options() {
  refinement_.add_options()
    ("refinement.v-cycles",
     po::value<std::string>()->default_value("final-only"),
     "V-cycle setting. Options:\n"
     "  final-only: \tonly iterate from the final partition,\n"
     "  best: \tonly iterate from the best intermediate partition at the end of "
     "each intermediate stage,\n"
     "  off: \tdon't use V-cycle refinement.")

    ("refinement.v-cycle-iteration-limit",
     po::value<int>()->default_value(INT_MAX, "INT_MAX"),
     "Limit on number of V-cycle iterations (if being used).")

    ("refinement.acceptable-gain", po::value<int>()->default_value(0),
     "Minimum acceptable gain of a parallel V-cycle (if being used) iteration "
     "in order to perform another V-cycle iteration. The value is as a "
     "percentage of the cutsize of the partition before the first V-cycle "
     "iteration.")

    ("refinement.acceptance-threshold", po::value<int>()->default_value(70),
     "Percentage threshold used to determine whether to accept or reject "
     "partitions during the parallel uncoarsening phase.")

    ("refinement.threshold-reduction", po::value<int>()->default_value(70),
     "Reduction in 'refinement.acceptance-threshold' appled at each "
     "uncoarsening step. Percentage threshold for iteration i+1 (t[i+1]) "
     "is equal to refinemen.threshold-reduction * t[i]/100.")

    ("refinement.enable-early-exit", po::bool_switch()->default_value(false),
     "Enable early exit criterion.")

    ("refinement.early-exit", po::value<int>()->default_value(100),
     "Percentage of consecutive local vertices that are visited, all of whose "
     "moves do not result in positive gain in the objective function. Must be "
     "between 0 and 100.")

    ("refinement.limit-by-length", po::bool_switch()->default_value(false),
     "Limit the length of hyperedges which may be refined. The length is "
     "determined by 'coarsening.reduction-ratio'.")

    // TODO(gb610): add to config file -- remove previous two options?!
    ("refinement.approximate", po::value<int>()->default_value(0),
     "Approximate refinement. Options:\n"
     "  0: basic,\n"
     "  1: approximate,\n"
     "  2: early exit,\n"
     "  3: approximate and early exit.")
  ;
}

void options::process_variables() {
  if (variables_.count("help")) {
    std::cout << general_ << std::endl;
    std::exit(1);
  } else if (variables_.count("version")) {
    std::cout << "Parkway v" << VERSION_STRING << std::endl;
    std::exit(1);
  } else if (variables_.count("display-all-options")) {
    std::cout << all_options_ << std::endl;
    std::exit(1);
  } else if (variables_.count("generate-config-file")) {
    generate_options_file();
    std::exit(1);
  } else if (variables_.count("config")) {
    load_config_file();
  }
}


void options::load_config_file() {
  std::ifstream config_file(get<std::string>("config"));
  if (config_file.bad() || !config_file.is_open()) {
    std::cerr << "[Error!] Coult not open config file '"
        << get<std::string>("config") << "'" << std::endl;
    std::exit(1);
  }

  po::store(po::parse_config_file(config_file, all_options_, true),
            variables_);
  config_file.close();
}


template<typename Type>
bool options::check_less_than(const char *key, Type limit) const {
  Type value = get<Type>(key);
  if (!(value < limit)) {
    std::cerr << "[Error!] Option '" << key << "' (" << value
        << ") is invalid and must be less than " << limit << std::endl;
    return false;
  }
  return true;
}

template<typename Type>
bool options::check_less_than_equal(const char *key, Type limit) const {
  Type value = get<Type>(key);
  if (!(value <= limit)) {
    std::cerr << "[Error!] Option '" << key << "' (" << value
        << ") is invalid and must be less than or equal to " << limit
        << std::endl;
    return false;
  }
  return true;
}

template<typename Type>
bool options::check_greater_than(const char *key, Type limit) const {
  Type value = get<Type>(key);
  if (!(value > limit)) {
    std::cerr << "[Error!] Option '" << key << "' (" << value
        << ") is invalid and must be greater than " << limit << std::endl;
    return false;
  }
  return true;
}

template<typename Type>
bool options::check_greater_than_equal(const char *key, Type limit) const {
  Type value = get<Type>(key);
  if (!(value >= limit)) {
    std::cerr << "[Error!] Option '" << key << "' (" << value
        << ") is invalid and must be greater than or equal to " << limit
        << std::endl;
    return false;
  }
  return true;
}

template<typename Type>
bool options::check_between(const char *key, Type lower, Type upper) const {
  Type value = get<Type>(key);
  if (!(lower <= value && value <= upper)) {
    std::cerr << "[Error!] Option '" << key << "' (" << value
        << ") is invalid and must be between " << lower << " and " << upper
        << std::endl;
    return false;
  }
  return true;
}

template<typename Type>
bool options::check_in_set(const char *key,
                           std::unordered_set<Type> &&values) const {
  Type value = get<Type>(key);
  if (values.count(value) == 0) {
    std::cerr << "[Error!] Option '" << key << "' (" << value
        << ") is invalid and must be one of: ";
    for (const auto &option : values) {
      std::cerr << "'" << option << "' ";
    }
    std::cerr << std::endl;
    return false;
  }
  return true;
}

bool options::check_clash(const char *key1, const char *key2) const {
  bool option1 = get<bool>(key1);
  bool option2 = get<bool>(key2);
  if (option1 && option2) {
    std::cerr << "[Error!] Option '" << key1 << "' and '" << key2 << "' "
        << "can't both be enabled." << std::endl;
    return false;
  }
  return true;
}

void options::validate_variables() const {
  bool okay = true;
  // General options.
  okay &= check_greater_than<int>("number-of-runs", 0);
  okay &= check_greater_than<int>("number-of-parts", 1);
  okay &= check_greater_than<double>("balance-constraint", 0.0);
  okay &= check_between("vertex-to-processor-allocation", 0, 2);
  okay &= check_greater_than_equal<int>("sprng-seed", 0);

  // Coarsening options.
  okay &= check_in_set<std::string>("coarsening.type", {"first-choice",
                                    "model-2d"});
  okay &= check_between<int>("coarsening.percentile-cutoff", 0, 100);
  okay &= check_greater_than_equal<int>("coarsening.percentile-increment", 0);
  okay &= check_greater_than<int>("coarsening.minimum-coarse-vertices", 0);
  okay &= check_between<double>("coarsening.reduction-ratio", 1.1, 3.0);
  okay &= check_in_set<std::string>("coarsening.vertex-visit-order",
                                    {"random", "increasing", "decreasing",
                                    "increasing-weight", "decreasing-weight"});
  okay &= check_between<int>("coarsening.connectivity-metric", 0, 3);

  // Serial partitioning.
  okay &= check_greater_than<int>("serial-partitioning.number-of-runs", 0);
  okay &= check_clash("serial-partitioning.use-patoh",
                      "serial-partitioning.use-hmetis");
  okay &= check_in_set<std::string>("recursive-bisection.v-cycles",
                            {"final-only", "all", "off"});
  okay &= check_greater_than<int>("recursive-bisection.number-of-runs", 0);
  okay &= check_in_set<std::string>("patoh.partition-options",
                            {"default", "speed", "quality"});
  okay &= check_in_set<std::string>("hmetis.coarsening", {"first-choice", "hybrid",
                            "greedy", "hyperedge", "edge"});
  okay &= check_in_set<std::string>("hmetis.refinement", {"best-intermediate",
                            "each-intermediate", "final", "off"});

  // Refinement options.
  okay &= check_in_set<std::string>("refinement.v-cycles",
                            {"final-only", "best", "off"});
  okay &= check_greater_than<int>("refinement.v-cycle-iteration-limit", 0);
  okay &= check_between<int>("refinement.acceptable-gain", 0, 100);
  okay &= check_between<int>("refinement.acceptance-threshold", 0, 100);
  okay &= check_between<int>("refinement.threshold-reduction", 0, 100);
  okay &= check_between<int>("refinement.early-exit", 0, 100);

  if (!okay) {
    std::exit(1);
  }
}




void options::generate_options_file() const {
  std::string filename("parkway-configuartion.ini");
  std::ofstream options_file(filename);
  if (!options_file.is_open()) {
    std::cerr << "\nCould not open '" << filename << "'" << std::endl;
    std::exit(1);
  }
  options_file <<
    "# --- Parkway Configuration ---------------------------------------------------\n"
    "\n"
    "# General options.\n"
    "\n"
    "# Hypergraph file to partition.\n"
    "hypergraph =\n"
    "# Number of parallel partitioning runs.\n"
    "number-of-runs = 1\n"
    "# Number of parts sought in partition.\n"
    "number-of-parts = 4\n"
    "# Balance constraint.\n"
    "balance-constraint = 0.05\n"
    "# Vertex to processor allocation.\n"
    "# Options:\n"
    "#   0: as read in,\n"
    "#   1: random,\n"
    "#   2: as prescribed in partition file.\n"
    "vertex-to-processor-allocation = 0\n"
    "# Seed for pseudo-random number generator (if compiled with SPRNG).\n"
    "sprng-seed = 1\n"
    "# Display partitioning information (e.g. settings for each component).\n"
    "display-info = false\n"
    "# Display partitioning progress (e.g. number of vertices during coarsening).\n"
    "display-progress = false\n"
    "# Write info and progress to a file (writes to the terminal if this is false).\n"
    "output-to-file = false\n"
    "# If output-to-file is enabled, write the info and/or progress to this file.\n"
    "output-file arg = parkway_output.log\n"
    "# Write the computed partitions to file.\n"
    "write-partitions-to-file = false\n"
    "# Randomly shuffle vertices between processes before coarsening.\n"
    "random-vertex-shuffle = false\n"
    "\n"
    "[coarsening]\n"
    "# Type of coarsener.\n"
    "# Options:\n"
    "#   first-choice,\n"
    "#   model-2d.\n"
    "type = first-choice\n"
    "# Percentile of hyperedge weight l such that the given percent of the total\n"
    "# percent of the total hyperedge weight in the hypergraph is found in hyperedges\n"
    "# with length less than or equal to l.\n"
    "percentile-cutoff = 100\n"
    "# Increment for the 'percentile-cutoff' option used in subsequent coarsening\n"
    "# steps assuming the required percentile does not exceed 100.\n"
    "percentile-increment = 1\n"
    "# Minimum number of vertices in coarsest hypergraph for the parallel\n"
    "# coarsening algorithm.\n"
    "minimum-coarse-vertices = 200\n"
    "# Reduction ratio during parallel coarsening. Typically between 1.5 and 2.0.\n"
    "reduction-ratio = 1.75\n"
    "# Local vertex visit order during parallel coarsening phase.\n"
    "# Options:\n"
    "#   random: random order,\n"
    "#   increasing: increasing order by index,\n"
    "#   decreasing: decreasing order by index,\n"
    "#   increasing-weight: increasing order by weight,\n"
    "#   decreasing-weight: decreasing order by weight.\n"
    "vertex-visit-order = random\n"
    "# Connectivity metric.\n"
    "# Options:\n"
    "#   0: cluster weight and hyperedge length not taken into account,\n"
    "#   1: inversely proportional to resulting cluster weight only,\n"
    "#   2: inversely proportional to lengths of connecting hyperedges only,\n"
    "#   3: inversely proportional to both resulting cluster weight and connecting\n"
    "#      hyperedge lenghts.\n"
    "connectivity-metric = 0\n"
    "# Process match requests in a random order from other processes.\n"
    "randomly-process-match-requests = false\n"
    "\n"
    "[serial-partitioning]\n"
    "# Number of serial partitioning runs. If hMETIS or PaToH are used, each\n"
    "# process will compute max(number-of-serial-runs/p, 1) runs where p is the\n"
    "# number of processes.\n"
    "number-of-runs = 1\n"
    "# Fan out during serial coarsening.\n"
    "util-fan-out = true\n"
    "# Divide by cluster weight during serial coarsening.\n"
    "divide-by-cluster-weight = true\n"
    "# Use hMETIS to perform serial partitioning (if available).\n"
    "use-hmetis = false\n"
    "# Use PaToH to perform serial partitioning (if available).\n"
    "use-patoh = false\n"
    "\n"
    "# Recursive bisection options (only used if both 'use-hmetis' and 'use-patoh'\n"
    "# are false).\n"
    "[recursive-bisection]\n"
    "# V-cycle setting.\n"
    "# Options:\n"
    "#   final-only: only use V-cycle refinement in the final bisection,\n"
    "#   all: use V-cycle refinement in each recursive bisection,\n"
    "#   off: don't use V-cycle refinement.\n"
    "v-cycles = final-only\n"
    "# Use a metric which is inversely proportional to the resulting cluster weight.\n"
    "use-weight-metric = false\n"
    "# Number of runs used during bisection.\n"
    "number-of-runs = 2\n"
    "# Number of bisection runs used on the coarsest hypergraph.\n"
    "number-of-initial-partitioning-runs = 10\n"
    "\n"
    "# PaToH (only if compiled with PaToH and use-patoh = true)\n"
    "[patoh]\n"
    "# PaToH settings.\n"
    "# Options:\n"
    "#   default: use PATOH_SUGPARAM_DEFAULT\n"
    "#   speed: use PATOH_SUGPARAM_SPEED\n"
    "#   quality: use PATOH_SUGPARAM_QUALITY\n"
    "partition-options = default\n"
    "\n"
    "# hMETIS (only if compiled with hMETIS and use-hmetis = true)\n"
    "[hmetis]\n"
    "# Sets the coarsening scheme used by hMETIS.\n"
    "# Options:\n"
    "#   first-choice\n"
    "#   hybrid: hybrid first choice\n"
    "#   greedy: greedy first choice\n"
    "#   hyperedge\n"
    "#   edge\n"
    "coarsening = first-choice\n"
    "# Sets the V-cycle refinement scheme used by hMETIS.\n"
    "# Options:\n"
    "#   best-intermediate\n"
    "#   each-intermediate\n"
    "#   final\n"
    "#   off\n"
    "refinement = best-intermediate\n"
    "\n"
    "[refinement]\n"
    "# V-cycle setting.\n"
    "# Options:\n"
    "#   final-only: only iterate from the final partition,\n"
    "#   best: only iterate from the best intermediate partition at the end of each\n"
    "#         intermediate stage,\n"
    "#   off: don't use V-cycle refinement.\n"
    "v-cycles = final-only\n"
    "# Limit on number of V-cycle iterations (if being used).\n"
    "v-cycle-iteration-limit = 32767\n"
    "# Minimum acceptable gain of a parallel V-cycle (if being used) iteration in\n"
    "# order to perform another V-cycle iteration. The value is as a percentage of\n"
    "# the cutsize of the partition before the first V-cycle iteration.\n"
    "acceptable-gain = 0\n"
    "# Percentage threshold used to determine whether to accept or reject\n"
    "# partitions during the parallel uncoarsening phase.\n"
    "acceptance-threshold = 70\n"
    "# Reduction in 'refinement.acceptance-threshold' appled at each uncoarsening\n"
    "# step. Percentage threshold for iteration i+1 (t[i+1]) is equal to\n"
    "# threshold-reduction * t[i]/100.\n"
    "threshold-reduction = 70\n"
    "# Enable early exit criterion.\n"
    "enable-early-exit = false\n"
    "# Percentage of consecutive local vertices that are visited, all of whose\n"
    "# moves do not result in positive gain in the objective function. Must be\n"
    "# between 0 and 100.\n"
    "early-exit = 100\n"
    "# Limit the length of hyperedges which may be refined. The length is\n"
    "# determined by 'coarsening.reduction-ratio'.\n"
    "limit-by-length = false\n";
  options_file.close();
  std::cout << "Saved configuartion file as '" << filename << "'" << std::endl;
}


}  // namespace parkway
