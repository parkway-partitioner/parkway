#ifndef _UTILS_H
#define _UTILS_H

// ### Utils.h ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 10/1/2005: Last Modified
//
// ###

#include <math.h>
#include <time.h>

#include "mpi.h"
#include "Log.h"
#include "hypergraph/parallel/hypergraph.hpp"
#include "basic_contoller.hpp"
#include "parallel_v_cycle_final_controller.hpp"
#include "parallel_v_cycle_all_controller.hpp"
#include "recursive_bisection_contoller.hpp"
#include "KHMetisController.hpp"
#include "PaToHController.hpp"

namespace parallel = parkway::parallel;
namespace serial = parkway::serial;

namespace Utils {
parallel::coarsener *buildParaCoarsener(int myRank, int numProc, int numParts,
                                  double constraint, parallel::hypergraph *h,
                                  std::ostream &out, const int *options,
                                  MPI_Comm comm);

parallel::restrictive_coarsening *buildParaRestrCoarsener(int myRank, int numProc,
                                            int numParts, double constraint,
                                            parallel::hypergraph *h, std::ostream &out,
                                            const int *options, MPI_Comm comm);

parallel::refiner *buildParaRefiner(int myRank, int numProc, int numParts,
                              double constraint, parallel::hypergraph *h,
                              std::ostream &out, const int *options, MPI_Comm comm);

serial::controller *buildSeqController(int myRank, int numProc, int numParts,
                                  double constraint, std::ostream &out,
                                  const int *options);

parallel::controller *buildParaController(int myRank, int numProcs, int numParts,
                                    int num_tot_verts, double constraint,
                                    parallel::coarsener *c, parallel::restrictive_coarsening *rc,
                                    parallel::refiner *r, serial::controller *s,
                                    std::ostream &out, const int *options,
                                    MPI_Comm comm);

void initDefaultValues(const int *userOptions, int *programOptions);

void checkPartsAndProcs(int num_parts, int num_procs, int seqOption,
                        int paraOption, std::ostream &out, MPI_Comm comm);
}

#endif
