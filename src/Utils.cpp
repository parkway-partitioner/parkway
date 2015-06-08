#ifndef _UTILS_CPP
#define _UTILS_CPP

// ### Utils.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 4/1/2005: Last Modified
//
// ###

#include "Utils.h"

using namespace std;

ParaCoarsener *Utils::buildParaCoarsener(int my_rank, int num_proc,
                                         int num_parts, double constraint,
                                         ParaHypergraph *h, ostream &out,
                                         const int *options, MPI_Comm comm) {
  ParaCoarsener *c = NULL;

  int coarsener_type = ParaFCC; // Para2DModel
  int min_nodes = options[8];
  int disp_option = options[2];
  int numTotPins = h->getNumTotPins(comm);

  double r = static_cast<double>(options[9]) / options[10];
  double aveVertDeg = h->getAveVertDeg(comm);
  double aveHedgeSize = h->getAveHedgeSize(comm);

  if (coarsener_type == ParaFCC) {
    int vertexVisitOrder = options[11];
    int divideConnectivity = options[12];
    int matchReqVisitOrder = options[13];
    int divByCluWt;
    int divByHedgeLen;

    switch (divideConnectivity) {
    case 0:
      divByCluWt = 0;
      divByHedgeLen = 0;
      break;
    case 1:
      divByCluWt = 1;
      divByHedgeLen = 0;
      break;
    case 2:
      divByCluWt = 0;
      divByHedgeLen = 1;
      break;
    default:
      divByCluWt = 1;
      divByHedgeLen = 1;
      break;
    }

    c = new ParaFCCoarsener(my_rank, num_proc, num_parts, vertexVisitOrder,
                            matchReqVisitOrder, divByCluWt, divByHedgeLen, out);

    c->setMinNodes(min_nodes * num_parts);
    c->setDispOption(disp_option);
    c->setBalConstraint(constraint);
    c->setReductionRatio(r);
    c->buildAuxiliaryStructs(numTotPins, aveVertDeg, aveHedgeSize);
  }

  if (coarsener_type == Para2DModel) {
    int vertexVisitOrder = options[11];
    int divideConnectivity = options[12];
    int matchReqVisitOrder = options[13];
    int divByCluWt;
    int divByHedgeLen;

    switch (divideConnectivity) {
    case 0:
      divByCluWt = 0;
      divByHedgeLen = 0;
      break;
    case 1:
      divByCluWt = 1;
      divByHedgeLen = 0;
      break;
    case 2:
      divByCluWt = 0;
      divByHedgeLen = 1;
      break;
    default:
      divByCluWt = 1;
      divByHedgeLen = 1;
      break;
    }

    c = new Para2DModelCoarsener(my_rank, num_proc, num_parts, vertexVisitOrder,
                                 matchReqVisitOrder, divByCluWt, divByHedgeLen,
                                 out);

    c->setMinNodes(min_nodes * num_parts);
    c->setDispOption(disp_option);
    c->setBalConstraint(constraint);
    c->setReductionRatio(r);
    c->buildAuxiliaryStructs(numTotPins, aveVertDeg, aveHedgeSize);
  }

  if (c && my_rank == 0)
    c->dispCoarseningOptions();

  return c;
}

ParaRestrCoarsener *Utils::buildParaRestrCoarsener(
    int my_rank, int num_proc, int num_parts, double constraint,
    ParaHypergraph *h, ostream &out, const int *options, MPI_Comm comm) {
  int coarsener_type = ParaRestFCC;
  int min_nodes = options[8];
  int disp_option = options[2];

  double r = static_cast<double>(options[9]) / options[10];

  ParaRestrCoarsener *c = NULL;

  if (coarsener_type == ParaRestFCC) {
    int vertexVisitOrder = options[11];
    int divideConnectivity = options[12];
    int divByCluWt;
    int divByHedgeLen;

    switch (divideConnectivity) {
    case 0:
      divByCluWt = 0;
      divByHedgeLen = 0;
      break;
    case 1:
      divByCluWt = 1;
      divByHedgeLen = 0;
      break;
    case 2:
      divByCluWt = 0;
      divByHedgeLen = 1;
      break;
    default:
      divByCluWt = 1;
      divByHedgeLen = 1;
      break;
    }

    c = new ParaRestrFCCoarsener(my_rank, num_proc, num_parts, vertexVisitOrder,
                                 divByCluWt, divByHedgeLen, out);

    c->setMinNodes(min_nodes * num_parts);
    c->setDispOption(disp_option);
    c->setBalConstraint(constraint);
    c->setReductionRatio(r);
  }

  if (c && my_rank == 0)
    c->dispCoarseningOptions();

  return c;
}

ParaRefiner *Utils::buildParaRefiner(int my_rank, int num_proc, int num_parts,
                                     double constraint, ParaHypergraph *h,
                                     ostream &out, const int *options,
                                     MPI_Comm comm) {
  int refiner_type = ParaGreedyKway;
  int disp_option = options[2];
  int numTotPins = h->getNumTotPins(comm);

  ParaRefiner *r = NULL;

  if (refiner_type == ParaGreedyKway) {
    double eeLimit = static_cast<double>(options[27]) / 100;
    int earlyExit = 0;

    if (options[28] == 2 || options[28] == 3)
      earlyExit = 1;

    r = new ParaGreedyKwayRefiner(my_rank, num_proc, num_parts,
                                  numTotPins / num_proc, earlyExit, eeLimit,
                                  out);

    r->setDispOption(disp_option);
    r->setBalConstraint(constraint);
  }

  if (r && my_rank == 0)
    r->dispRefinementOptions();

  return r;
}

SeqController *Utils::buildSeqController(int my_rank, int num_proc,
                                         int num_parts, double constraint,
                                         ostream &out, const int *options) {
  int numSeqRuns = options[14];
  int seqControllerType = options[15];
  int disp_option = options[2];

  SeqController *seqC = NULL;

  if (seqControllerType == RecurBisect) {
    int numBisectRuns = options[17];
    int eeParam = DEF_EE_PARAM;
    int startPercentile = 90;
    int inc = 0;

    double keepT = DEF_KEEP_THRESHOLD;
    double redFactor = DEF_REDUC_FACTOR;

    BisectionController *bC =
        new BisectionController(numBisectRuns, keepT, redFactor, eeParam,
                                startPercentile, inc, disp_option, out);

    // ###
    // build the sequential coarsener
    // ###

    int seqCoarsener = options[16];
    int minSeqNodes = MIN_VERT_MULTIPLIER;
    double bRedRatio = DEF_BIS_RATIO;

    bC->buildCoarsener(bRedRatio, seqCoarsener, minSeqNodes);

    // ###
    // build the seq initial partitioner
    // ###

    int numInitRuns = options[18];

    bC->buildInitBisector(numInitRuns);

    // ###
    // build the seq refiner (FM refiner)
    // ###

    int queueDiscipline = LIFO;

    bC->buildRefiner(queueDiscipline);

    // ###
    // build the seq GreedyKwayRefiner
    // ###

    double kWayLimit = DEF_SEQ_KWAY_LIM;

    GreedyKwayRefiner *k =
        new GreedyKwayRefiner(-1, num_parts, -1, kWayLimit, disp_option);

    // ###
    // build the seq controller
    // ###

    double paraKeepT = DEF_KEEP_THRESHOLD;

    seqC = new RecurBisectController(bC, k, my_rank, num_proc, num_parts,
                                     numBisectRuns, out);

    seqC->setAcceptProp(paraKeepT);
    seqC->setNumSeqRuns(numSeqRuns);
    seqC->setDispOption(disp_option);
    seqC->setKwayConstraint(constraint);
  }

  if (seqControllerType == VCycleAllRecurBisect) {
    int numBisectRuns = options[17];
    int eeParam = DEF_EE_PARAM;
    int startPercentile = 95;
    int inc = 0;

    double keepT = DEF_KEEP_THRESHOLD;
    double redFactor = DEF_REDUC_FACTOR;

    VCycleAllBisectionController *bC = new VCycleAllBisectionController(
        numBisectRuns, keepT, redFactor, eeParam, startPercentile, inc,
        disp_option, out);

    // ###
    // build the sequential coarsener
    // ###

    int seqCoarsener = options[16];
    int minSeqNodes = MIN_VERT_MULTIPLIER;
    double bRedRatio = DEF_BIS_RATIO;

    bC->buildCoarsener(bRedRatio, seqCoarsener, minSeqNodes);

    // ###
    // build the restr coarsener
    // ###

    int restrSeqCoarsener = options[17];

    bC->buildRestrCoarsener(bRedRatio, restrSeqCoarsener, minSeqNodes);

    // ###
    // build the seq initial partitioner
    // ###

    int numInitRuns = options[18];

    bC->buildInitBisector(numInitRuns);

    // ###
    // build the seq refiner (FM refiner)
    // ###

    int queueDiscipline = LIFO;

    bC->buildRefiner(queueDiscipline);

    // ###
    // build the seq GreedyKwayRefiner
    // ###

    double kWayLimit = DEF_SEQ_KWAY_LIM;

    GreedyKwayRefiner *k =
        new GreedyKwayRefiner(-1, num_parts, -1, kWayLimit, disp_option);

    // ###
    // build the seq controller
    // ###

    double paraKeepT = DEF_KEEP_THRESHOLD;

    seqC = new RecurBisectController(bC, k, my_rank, num_proc, num_parts,
                                     numBisectRuns, out);

    seqC->setAcceptProp(paraKeepT);
    seqC->setNumSeqRuns(numSeqRuns);
    seqC->setDispOption(disp_option);
    seqC->setKwayConstraint(constraint);
  }

  if (seqControllerType == VCycleFinalRecurBisect) {
    int numBisectRuns = options[17];
    int eeParam = DEF_EE_PARAM;
    int startPercentile = 95;
    int inc = 0;

    double keepT = DEF_KEEP_THRESHOLD;
    double redFactor = DEF_REDUC_FACTOR;

    VCycleFinalBisectionController *bC = new VCycleFinalBisectionController(
        numBisectRuns, keepT, redFactor, eeParam, startPercentile, inc,
        disp_option, out);

    // ###
    // build the sequential coarsener
    // ###

    int seqCoarsener = options[16];
    int minSeqNodes = MIN_VERT_MULTIPLIER;
    double bRedRatio = DEF_BIS_RATIO;

    bC->buildCoarsener(bRedRatio, seqCoarsener, minSeqNodes);

    // ###
    // build the restr coarsener
    // ###

    int restrSeqCoarsener = options[16];

    bC->buildRestrCoarsener(bRedRatio, restrSeqCoarsener, minSeqNodes);

    // ###
    // build the seq initial partitioner
    // ###

    int numInitRuns = options[18];

    bC->buildInitBisector(numInitRuns);

    // ###
    // build the seq refiner (FM refiner)
    // ###

    int queueDiscipline = LIFO;

    bC->buildRefiner(queueDiscipline);

    // ###
    // build the seq GreedyKwayRefiner
    // ###

    double kWayLimit = DEF_SEQ_KWAY_LIM;

    GreedyKwayRefiner *k =
        new GreedyKwayRefiner(-1, num_parts, -1, kWayLimit, disp_option);

    // ###
    // build the seq controller
    // ###

    double paraKeepT = DEF_KEEP_THRESHOLD;

    seqC = new RecurBisectController(bC, k, my_rank, num_proc, num_parts,
                                     numBisectRuns, out);

    seqC->setAcceptProp(paraKeepT);
    seqC->setNumSeqRuns(numSeqRuns);
    seqC->setDispOption(disp_option);
    seqC->setKwayConstraint(constraint);
  }

#ifdef LINK_HMETIS
  if (seqControllerType == KhMeTiS) {
    // ###
    // build the seq GreedyKwayRefiner
    // ###

    double kWayLimit = DEF_SEQ_KWAY_LIM;

    GreedyKwayRefiner *k =
        new GreedyKwayRefiner(-1, num_parts, -1, kWayLimit, disp_option);

    // ###
    // build the seq controller
    // ###

    DynamicArray<int> seq_options(12);

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
                                 seq_options.getArray(), out);

    seqC->setAcceptProp(paraKeepT);
    seqC->setNumSeqRuns(numSeqRuns);
    seqC->setDispOption(disp_option);
    seqC->setKwayConstraint(constraint);
  }
#endif
#ifdef LINK_PATOH
  if (seqControllerType == PaToH) {
    // ###
    // build the seq GreedyKwayRefiner
    // ###

    double kWayLimit = DEF_SEQ_KWAY_LIM;

    GreedyKwayRefiner *k =
        new GreedyKwayRefiner(-1, num_parts, -1, kWayLimit, disp_option);

    // ###
    // build the seq controller
    // ###

    double paraKeepT = DEF_PARA_KEEP_THRESHOLD;

    seqC =
        new PaToHController(k, my_rank, num_proc, num_parts, options[21], out);

    seqC->setAcceptProp(paraKeepT);
    seqC->setNumSeqRuns(numSeqRuns);
    seqC->setDispOption(disp_option);
    seqC->setKwayConstraint(constraint);
  }
#endif

  if (seqC && my_rank == 0)
    seqC->dispSeqControllerOptions();

  return seqC;
}

ParaController *Utils::buildParaController(
    int my_rank, int num_proc, int num_parts, int num_tot_verts,
    double constraint, ParaCoarsener *c, ParaRestrCoarsener *rc, ParaRefiner *r,
    SeqController *s, ostream &out, const int *options, MPI_Comm comm) {
  int paraControllerType = options[22];
  int numParaRuns = options[4];
  int disp_option = options[2];
  int shuffleVertices = options[5];
  int percentile = options[6];
  int perCentInc = options[7];
  int approxRef = options[28];

  double paraKeepT = static_cast<double>(options[25]) / 100;
  double redFactor = static_cast<double>(options[26]) / 100;

  ParaController *paraC = NULL;

  if (paraControllerType == BasicParaC) {
    paraC = new BasicParaController(*c, *r, *s, my_rank, num_proc, percentile,
                                    perCentInc, approxRef, out);
  }

  if (paraControllerType == ParaVCycleBig) {
    bool noCoarsening = (num_tot_verts <= options[8] * num_parts);

    if (!rc || noCoarsening) {
      char message[512];

      if (!rc && my_rank == 0) {
        sprintf(message,
                "p[%d] no ParaRestrCoarsener - create basic controller\n",
                my_rank);
        cout << message;
      }

      if (noCoarsening && my_rank == 0) {
        sprintf(message,
                "p[%d] max coarse hypergraph size > orig hypergraph size - no "
                "v-cycles possible\n",
                my_rank);
        cout << message;
      }

      paraC = new BasicParaController(*c, *r, *s, my_rank, num_proc, percentile,
                                      perCentInc, approxRef, out);
    } else {
      int limitOnCycles = options[23];
      double limitAsPercent = static_cast<double>(options[24]) / 100;

      paraC = new ParaVCycleFinalController(*rc, *c, *r, *s, my_rank, num_proc,
                                            percentile, perCentInc, approxRef,
                                            limitOnCycles, limitAsPercent, out);
    }
  }

  if (paraControllerType == ParaVCycleAll) {
    bool noCoarsening = (num_tot_verts <= options[8] * num_parts);

    if (!rc || noCoarsening) {
      char message[512];

      if (!rc && my_rank == 0) {
        sprintf(message,
                "p[%d] no ParaRestrCoarsener - create basic controller\n",
                my_rank);
        out << message;
      }

      if (noCoarsening && my_rank == 0) {
        sprintf(message,
                "p[%d] max coarse hypergraph size > orig hypergraph size - no "
                "v-cycles possible\n",
                my_rank);
        out << message;
      }

      paraC = new BasicParaController(*c, *r, *s, my_rank, num_proc, percentile,
                                      perCentInc, approxRef, out);
    } else {
      int limitOnCycles = options[23];
      double limitAsPercent = static_cast<double>(options[24]) / 100;

      paraC = new ParaVCycleAllController(*rc, *c, *r, *s, my_rank, num_proc,
                                          percentile, perCentInc, approxRef,
                                          limitOnCycles, limitAsPercent, out);
    }
  }

  if (paraC) {
    paraC->setNumParts(num_parts);
    paraC->setNumParaRuns(numParaRuns);
    paraC->setBalConstraint(constraint);
    paraC->setKTFactor(paraKeepT);
    paraC->setReduction(redFactor);
    paraC->setDispOpt(disp_option);
    paraC->setShuffleVertices(shuffleVertices);

    /* random vertex shuffle before each k-way refinement */

    paraC->setRandShuffBefRef(0);

    if (my_rank == 0)
      paraC->dispParaControllerOptions();
  }

  return paraC;
}

void Utils::initDefaultValues(const int *userOptions, int *programOptions) {
  double r;

  programOptions[0] = userOptions[0];

  if (programOptions[0] == 0) {
    // ###
    // set all to defaults

    programOptions[1] = 0;
    programOptions[2] = 0;
    programOptions[3] = 0;
    programOptions[4] = 1;
    programOptions[5] = 0;
    programOptions[6] = 100;
    programOptions[7] = 0;
    programOptions[8] = 200;
    programOptions[9] = 7;
    programOptions[10] = 4;
    programOptions[11] = 3;
    programOptions[12] = 0;
    programOptions[13] = 3;
    programOptions[14] = 1;
    programOptions[15] = 2;
    programOptions[16] = 1;
    programOptions[17] = 2;
    programOptions[18] = 10;
    programOptions[19] = 2;
    programOptions[20] = 2;
    programOptions[21] = 1;
    programOptions[22] = 1;
    programOptions[23] = LARGE_CONSTANT;
    programOptions[24] = 0;
    programOptions[25] = 70;
    programOptions[26] = 70;
    programOptions[27] = 100;
    programOptions[28] = 0;
  } else {
    // ###
    // otherwise set as requested

    if (userOptions[1] >= 0) {
      programOptions[1] = userOptions[1];
    } else {
      programOptions[1] = 0;
    }

    if (userOptions[2] >= 0 && userOptions[2] < 3) {
      programOptions[2] = userOptions[2];
    } else {
      programOptions[2] = 0;
    }

    if (userOptions[3] >= 0 && userOptions[3] < 3) {
      programOptions[3] = userOptions[3];
    } else {
      programOptions[3] = 0;
    }

    if (userOptions[4] > 0 && userOptions[4] < LARGE_CONSTANT) {
      programOptions[4] = userOptions[4];
    } else {
      programOptions[4] = 1;
    }

    if (userOptions[5] < 0 || userOptions[5] > 2) {
      programOptions[5] = 0;
    } else {
      programOptions[5] = userOptions[5];
    }

    if (userOptions[6] <= 100 && userOptions[6] > 0) {
      programOptions[6] = userOptions[6];
    } else {
      programOptions[6] = 100;
    }

    if (userOptions[7] <= (100 - programOptions[6]) && userOptions[7] >= 0) {
      programOptions[7] = userOptions[7];
    } else {
      programOptions[7] = 0;
    }

    if (userOptions[8] > 100 /*&& userOptions[8] < 1000*/) {
      programOptions[8] = userOptions[8];
    } else {
      programOptions[8] = 200;
    }

    r = static_cast<double>(userOptions[9]) / userOptions[10];

    if (r > 1.1 && r < 3.0) {
      programOptions[9] = userOptions[9];
      programOptions[10] = userOptions[10];
    } else {
      programOptions[9] = 7;
      programOptions[10] = 4;
    }

    if (userOptions[11] > 0 && userOptions[11] < 6) {
      programOptions[11] = userOptions[11];
    } else {
      programOptions[11] = 3;
    }

    if (userOptions[12] >= 0 && userOptions[12] < 4) {
      programOptions[12] = userOptions[12];
    } else {
      programOptions[12] = 3;
    }

    if (userOptions[13] == 2 || userOptions[13] == 3) {
      programOptions[13] = userOptions[13];
    } else {
      programOptions[13] = 3;
    }

    if (userOptions[14] > 0 && userOptions[14] < LARGE_CONSTANT) {
      programOptions[14] = userOptions[14];
    } else {
      programOptions[14] = 1;
    }

    if (userOptions[15] > 0 && userOptions[15] < 6) {
#ifndef LINK_HMETIS
      if (userOptions[15] == 4)
        programOptions[15] = 2;
      else
#endif
#ifndef LINK_PATOH
          if (userOptions[15] == 5)
        programOptions[15] = 2;
      else
#endif
        programOptions[15] = userOptions[15];
    } else {
      programOptions[15] = 2;
    }

    if (userOptions[16] == 1 || userOptions[16] == 2) {
      programOptions[16] = userOptions[16];
    } else {
      programOptions[16] = 1;
    }

    if (userOptions[17] > 0 && userOptions[17] < LARGE_CONSTANT) {
      programOptions[17] = userOptions[17];
    } else {
      programOptions[17] = 2;
    }

    if (userOptions[18] > 0 && userOptions[18] < LARGE_CONSTANT) {
      programOptions[18] = userOptions[18];
    } else {
      programOptions[18] = 10;
    }

    if (userOptions[19] > 0 && userOptions[19] < 6) {
      programOptions[19] = userOptions[19];
    } else {
      programOptions[19] = 2;
    }

    if (userOptions[20] >= 0 && userOptions[20] < 4) {
      programOptions[20] = userOptions[20];
    } else {
      programOptions[20] = 2;
    }

    if (userOptions[21] > 0 && userOptions[21] < 4) {
      programOptions[21] = userOptions[21];
    } else {
      programOptions[21] = 1;
    }

    if (userOptions[22] > 0 && userOptions[22] < 4) {
      programOptions[22] = userOptions[22];
    } else {
      programOptions[22] = 1;
    }

    if (userOptions[23] > 0 && userOptions[23] < LARGE_CONSTANT) {
      programOptions[23] = userOptions[23];
    } else {
      programOptions[23] = LARGE_CONSTANT;
    }

    if (userOptions[24] >= 0 && userOptions[24] < 100) {
      programOptions[24] = userOptions[24];
    } else {
      programOptions[24] = 0;
    }

    if (userOptions[25] > 0 && userOptions[25] < 100) {
      programOptions[25] = userOptions[25];
    } else {
      programOptions[25] = 0;
    }

    if (userOptions[26] > 0 && userOptions[26] < 100) {
      programOptions[26] = userOptions[26];
    } else {
      programOptions[26] = 70;
    }

    if (userOptions[27] > 0 && userOptions[27] <= 100) {
      programOptions[27] = userOptions[27];
    } else {
      programOptions[27] = 100;
    }

    if (userOptions[28] >= 0 || userOptions[28] < 4) {
      programOptions[28] = userOptions[28];
    } else {
      programOptions[28] = 0;
    }
  }
}

void Utils::checkPartsAndProcs(int num_parts, int num_procs, int seqOption,
                               int paraOption, ostream &out, MPI_Comm comm) {
  int my_rank;
  MPI_Comm_rank(comm, &my_rank);

  if (paraOption > 1) {
    if (!Funct::isPowerOf2(num_procs)) {
      MPI_Barrier(comm);
      if (my_rank == 0)
        out << "number of processors must be a power of two when using "
               "parallel V-Cycling - abort\n";

      MPI_Abort(comm, 0);
    }

    if (num_parts < num_procs || num_parts % num_procs > 0) {
      MPI_Barrier(comm);
      if (my_rank == 0)
        out << "number of parts in partition must be divisible by number of "
               "processors when using parallel V-Cycling - abort\n";

      MPI_Abort(comm, 0);
    }

    if (seqOption <= 3 && !Funct::isPowerOf2(num_parts)) {
      MPI_Barrier(comm);
      if (my_rank == 0)
        out << "number of parts in partition must be a power of two when using "
               "generic recursive bisection - abort\n";

      MPI_Abort(comm, 0);
    }
  } else {
    if (num_procs < 2) {
      MPI_Barrier(comm);
      if (my_rank == 0)
        out << "number of processors must be greater than one - abort\n";

      MPI_Abort(comm, 0);
    }

    if (num_parts < 2) {
      MPI_Barrier(comm);
      if (my_rank == 0)
        out << "number of parts must be greater than one - abort\n";

      MPI_Abort(comm, 0);
    }

    if (seqOption <= 3 && !Funct::isPowerOf2(num_parts)) {
      MPI_Barrier(comm);
      if (my_rank == 0)
        out << "number of parts must be a power of two when using generic "
               "recursive bisection - abort\n";
      MPI_Abort(comm, 0);
    }
  }
}

#endif
