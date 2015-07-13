#include "utility/component_builders.hpp"

namespace parkway {

parallel::coarsener *build_parallel_coarsener(
    int rank, const parkway::options &options, parallel::hypergraph *h,
    MPI_Comm comm) {
  int numTotPins = h->total_number_of_pins(comm);
  double aveVertDeg = h->average_vertex_degree(comm);
  double aveHedgeSize = h->average_hyperedge_size(comm);

  const std::string &visit_order =
      options.get<std::string>("coarsening.vertex-visit-order");

  int vertexVisitOrder = get_vertex_visit_order(rank, visit_order);
  if (vertexVisitOrder == 0) {
    MPI_Abort(comm, 0);
  }

  int divByCluWt;
  int divByHedgeLen;
  get_connectivity_values(options.get<int>("coarsening.connectivity-metric"),
                          divByCluWt, divByHedgeLen);

  int matchReqVisitOrder;
  if (options.get<bool>("coarsening.randomly-process-match-requests")) {
    matchReqVisitOrder = 3;
  } else {
    matchReqVisitOrder = 2;
  }

  int num_proc = options.number_of_processors();
  int num_parts = options.get<int>("number-of-parts");
  parallel::coarsener *c = nullptr;
  if (options.get<std::string>("coarsening.type") == "first-choice") {
    c = new parallel::first_choice_coarsener(
        rank, num_proc, num_parts, vertexVisitOrder, matchReqVisitOrder,
        divByCluWt, divByHedgeLen);
  } else if (options.get<std::string>("coarsening.type") == "model-2d") {
    c = new parallel::model_coarsener_2d(
        rank, num_proc, num_parts, vertexVisitOrder, matchReqVisitOrder,
        divByCluWt, divByHedgeLen);
  }

  if (c) {
    int min_nodes = options.get<int>("coarsening.minimum-coarse-vertices");
    double reduction_ratio = options.get<double>("coarsening.reduction-ratio");
    double constraint = options.get<double>("balance-constraint");
    c->set_minimum_number_of_nodes(min_nodes * num_parts);
    c->set_balance_constraint(constraint);
    c->set_reduction_ratio(reduction_ratio);
    c->build_auxiliary_structures(numTotPins, aveVertDeg, aveHedgeSize);
    c->display_options();
  }

  return c;
}


parallel::restrictive_coarsening *build_parallel_restrictive_coarsener(
    int rank, const parkway::options &options, parallel::hypergraph *h,
    MPI_Comm comm) {
  const std::string &visit_order =
      options.get<std::string>("coarsening.vertex-visit-order");

  int vertexVisitOrder = get_vertex_visit_order(rank, visit_order);
  if (vertexVisitOrder == 0) {
    MPI_Abort(comm, 0);
  }

  int divByCluWt;
  int divByHedgeLen;
  get_connectivity_values(options.get<int>("coarsening.connectivity-metric"),
                          divByCluWt, divByHedgeLen);

  int num_proc = options.number_of_processors();
  int num_parts = options.get<int>("number-of-parts");
  parallel::restrictive_coarsening *c =
      new parallel::restrictive_first_choice_coarsening(
          rank, num_proc, num_parts, vertexVisitOrder, divByCluWt,
          divByHedgeLen);

  if (c) {
    int min_nodes = options.get<int>("coarsening.minimum-coarse-vertices");
    double reduction_ratio = options.get<double>("coarsening.reduction-ratio");
    double constraint = options.get<double>("balance-constraint");
    c->set_minimum_nodes(min_nodes * num_parts);
    c->set_balance_constraint(constraint);
    c->set_reduction_ratio(reduction_ratio);
    c->display_options();
  }

  return c;
}

parallel::refiner *build_parallel_refiner(
    int rank, const parkway::options &options, parallel::hypergraph *h,
    MPI_Comm comm) {

  int numTotPins = h->total_number_of_pins(comm);

  int ee = options.get<int>("refinement.early-exit");
  double eeLimit = static_cast<double>(ee) / 100;


  int earlyExit = options.get<bool>("refinement.enable-early-exit") ||
                  options.get<bool>("refinement.limit-by-length");

  int num_proc = options.number_of_processors();
  int num_parts = options.get<int>("number-of-parts");
  parallel::refiner *r = new parallel::k_way_greedy_refiner(
      rank, num_proc, num_parts, numTotPins / num_proc, earlyExit, eeLimit);

  if (r) {
    r->set_balance_constraint(options.get<double>("balance-constraint"));
    r->display_options();
  }

  return r;
}

serial::controller *build_sequential_controller(
    int rank, const parkway::options &options) {
  int num_parts = options.get<int>("number-of-parts");
  int numSeqRuns = options.get<int>("serial-partitioning.number-of-runs");
  bool use_patoh = options.get<bool>("serial-partitioning.use-patoh");
  bool use_hmetis = options.get<bool>("serial-partitioning.use-hmetis");
  bool use_recursive_bisection = !(use_patoh || use_hmetis);
  serial::controller *seqC = nullptr;

  if (use_recursive_bisection) {
    const std::string &v_cycle =
        options.get<std::string>("recursive-bisection.v-cycles");

    int numBisectRuns = options.get<int>("recursive-bisection.number-of-runs");
    int eeParam = DEF_EE_PARAM;
    int startPercentile = 95;
    int inc = 0;

    double keepT = DEF_KEEP_THRESHOLD;
    double redFactor = DEF_REDUC_FACTOR;
    int minSeqNodes = MIN_VERT_MULTIPLIER;
    double bRedRatio = DEF_BIS_RATIO;

    serial::bisection_controller *bC = nullptr;
    if (v_cycle == "all") {
      bC = new serial::v_cycle_all(
          numBisectRuns, keepT, redFactor, eeParam, startPercentile, inc);
      bC->build_coarsener(options);
    } else if (v_cycle == "final-only") {
      bC = new serial::v_cycle_final(
          numBisectRuns, keepT, redFactor, eeParam, startPercentile, inc);
      bC->build_coarsener(options);
    } else {
      startPercentile = 90;
      bC = new serial::bisection_controller(
          numBisectRuns, keepT, redFactor, eeParam, startPercentile, inc);
      bC->build_coarsener(options);
    }

    // ###
    // build the seq initial partitioner
    // ###
    int numInitRuns = options.get<int>(
        "recursive-bisection.number-of-initial-partitioning-runs");
    bC->build_initial_bisector(numInitRuns);

    // ###
    // build the seq refiner (FM refiner)
    // ###
    bC->build_refiner(LIFO);

    // ###
    // build the seq GreedyKwayRefiner
    // ###
    double kWayLimit = DEF_SEQ_KWAY_LIM;
    serial::greedy_k_way_refiner *k =
        new serial::greedy_k_way_refiner(-1, num_parts, -1, kWayLimit);

    // ###
    // build the seq controller
    // ###
    double paraKeepT = DEF_KEEP_THRESHOLD;
    int num_proc = options.number_of_processors();
    int num_parts = options.get<int>("number-of-parts");
    seqC = new recursive_bisection_contoller(bC, k, rank, num_proc,
                                             num_parts, numBisectRuns);

    seqC->set_accept_proportion(paraKeepT);
    seqC->set_number_of_runs(numSeqRuns);
    seqC->set_k_way_constraint(options.get<double>("balance-constraint"));
  }

#ifdef PARKWAY_LINK_HMETIS
  if (use_hmetis) {
    // ###
    // build the seq GreedyKwayRefiner
    // ###
    double kWayLimit = DEF_SEQ_KWAY_LIM;

    GreedyKwayRefiner *k = new GreedyKwayRefiner(-1, num_parts, -1, kWayLimit);

    // ###
    // build the seq controller
    // ###
    dynamic_array<int> seq_options(12);

    seq_options[0] = USER_DEFINED;
    seq_options[2] = options[19];
    seq_options[3] = SOED;
    seq_options[4] = options[20];
    // seq_options[5] not used
    // seq_options[6] not used
    // seq_options[7] random seed;
    seq_options[8] = NO_DEBUG;

    double paraKeepT = DEF_PARA_KEEP_THRESHOLD;

    seqC = new KHMetisController(k, my_rank, num_proc, num_parts,
                                 seq_options.getArray();

    seqC->setAcceptProp(paraKeepT);
    seqC->setNumSeqRuns(numSeqRuns);
    seqC->setDispOption();
    seqC->setKwayConstraint(constraint);
  }
#endif
#ifdef PARKWAY_LINK_PATOH
  if (use_patoh) {
    // ###
    // build the seq GreedyKwayRefiner
    // ###
    double kWayLimit = DEF_SEQ_KWAY_LIM;

    GreedyKwayRefiner *k = new GreedyKwayRefiner(-1, num_parts, -1, kWayLimit);

    // ###
    // build the seq controller
    // ###

    double paraKeepT = DEF_PARA_KEEP_THRESHOLD;

    seqC = new PaToHController(k, my_rank, num_proc, num_parts, options[21]);

    seqC->setAcceptProp(paraKeepT);
    seqC->setNumSeqRuns(numSeqRuns);
    seqC->setDispOption();
    seqC->setKwayConstraint(constraint);
  }
#endif

  if (seqC) {
    seqC->display_options();
  }

  return seqC;
}

parallel::controller *build_parallel_controller(
    int rank, int num_tot_verts, parallel::coarsener *c,
    parallel::restrictive_coarsening *rc, parallel::refiner *r,
    serial::controller *s, const parkway::options &options, MPI_Comm comm) {
  const std::string &v_cycle = options.get<std::string>("refinement.v-cycles");

  int numParaRuns = options.get<int>("number-of-runs");
  int shuffleVertices = options.get<int>("vertex-to-processor-allocation");
  int percentile = options.get<int>("coarsening.percentile-cutoff");
  int perCentInc = options.get<int>("coarsening.percentile-increment");
  int approxRef = options.get<int>("refinement.approximate");
  int num_proc = options.number_of_processors();
  int num_parts = options.get<int>("number-of-parts");
  double constraint = options.get<double>("balance-constraint");

  double paraKeepT = static_cast<double>(
      options.get<int>("refinement.acceptance-threshold")) / 100;
  double redFactor = static_cast<double>(
      options.get<int>("refinement.threshold-reduction")) / 100;

  parallel::controller *paraC = nullptr;

  int min_coarse = options.get<int>("coarsening.minimum-coarse-vertices");
  bool noCoarsening = (num_tot_verts <= min_coarse * num_parts);

  if (v_cycle == "off") {
    paraC = new parallel::basic_contoller(*c, *r, *s, rank, num_proc,
                                          percentile, perCentInc, approxRef);
  } else if (!rc || noCoarsening) {
    if (!rc) {
      info("p[%d] no ParaRestrCoarsener - create basic controller\n", rank);
    } else {
      info("p[%d] max coarse hypergraph size > orig hypergraph size - no "
           "v-cycles possible\n", rank);
    }
    paraC = new parallel::basic_contoller(*c, *r, *s, rank, num_proc,
                                          percentile, perCentInc, approxRef);
  } else if (v_cycle == "final-only") {
    int limitOnCycles = options.get<int>("refinement.v-cycle-iteration-limit");
    double limitAsPercent = static_cast<double>(
        options.get<int>("refinement.acceptable-gain")) / 100;

    paraC = new parkway::parallel::v_cycle_final(
        *rc, *c, *r, *s, rank, num_proc, percentile, perCentInc, approxRef,
        limitOnCycles, limitAsPercent);
  } else if (v_cycle == "all") {
    int limitOnCycles = options.get<int>("refinement.v-cycle-iteration-limit");
    double limitAsPercent = static_cast<double>(
        options.get<int>("refinement.acceptable-gain")) / 100;
    paraC = new parkway::parallel::v_cycle_all(
        *rc, *c, *r, *s, rank, num_proc, percentile, perCentInc, approxRef,
        limitOnCycles, limitAsPercent);
  }

  if (paraC) {
    paraC->set_number_of_parts(num_parts);
    paraC->set_number_of_runs(numParaRuns);
    paraC->set_balance_constraint(constraint);
    paraC->set_kt_factor(paraKeepT);
    paraC->set_reduction_in_keep_threshold(redFactor);
    paraC->set_shuffle_vertices(shuffleVertices);

    /* random vertex shuffle before each k-way refinement */
    paraC->set_random_shuffle_before_refine(0);

    paraC->display_options();
  }

  return paraC;
}

void check_parts_and_processors(const parkway::options &options, MPI_Comm comm) {
  bool use_patoh = options.get<bool>("serial-partitioning.use-patoh");
  bool use_hmetis = options.get<bool>("serial-partitioning.use-hmetis");
  bool use_recursive_bisection = !(use_patoh || use_hmetis);
  bool use_parallel_v_cycles = (options.get<std::string>("refinement.v-cycles")
                                != "off");

  int number_of_processors = options.number_of_processors();
  int number_of_parts = options.get<int>("number-of-parts");

  if (use_recursive_bisection &&
      !parkway::utility::math::is_power_of_2(number_of_parts)) {
    error("number of parts in partition must be a power of two when using "
          "generic recursive bisection - abort\n");
    MPI_Abort(comm, 0);
  }

  if (use_parallel_v_cycles) {
    if (!parkway::utility::math::is_power_of_2(number_of_processors)) {
      error("number of processors must be a power of two when using parallel "
            "V-Cycling - abort\n");
      MPI_Abort(comm, 0);
    } else if (number_of_parts < number_of_processors ||
               number_of_parts % number_of_processors > 0) {
      error("number of parts in partition must be divisible by number of "
            "processors when using parallel V-Cycling - abort\n");
      MPI_Abort(comm, 0);
    }
  } else {
    if (number_of_processors < 2) {
      error("number of processors must be greater than one - abort\n");
      MPI_Abort(comm, 0);
    } else if (number_of_parts < 2) {
      error("number of parts must be greater than one - abort\n");
      MPI_Abort(comm, 0);
    }
  }
}

}  // namespace parkway
