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
#include "controllers/serial/recursive_bisection_contoller.hpp"
#include "controllers/parallel/basic_contoller.hpp"
#include "controllers/parallel/v_cycle_final.hpp"
#include "controllers/parallel/v_cycle_all.hpp"
#include "KHMetisController.hpp"
#include "PaToHController.hpp"

namespace parallel = parkway::parallel;
namespace serial = parkway::serial;

namespace Utils {
parallel::coarsener *
    buildParaCoarsener(int myRank, int numProc, int numParts, double constraint, parallel::hypergraph *h, const int *options, MPI_Comm comm);

parallel::restrictive_coarsening *buildParaRestrCoarsener(int myRank, int numProc, int numParts, double constraint, parallel::hypergraph *h, const int *options, MPI_Comm comm);

parallel::refiner *buildParaRefiner(int myRank, int numProc, int numParts,
                              double constraint, parallel::hypergraph *h,
                              const int *options, MPI_Comm comm);

serial::controller *buildSeqController(int myRank, int numProc, int numParts,
                                  double constraint, const int *options);

parallel::controller *buildParaController(int myRank, int numProcs, int numParts,
                                    int num_tot_verts, double constraint,
                                    parallel::coarsener *c, parallel::restrictive_coarsening *rc,
                                    parallel::refiner *r, serial::controller *s,
                                    const int *options,
                                    MPI_Comm comm);

void initDefaultValues(const int *userOptions, int *programOptions);

void checkPartsAndProcs(int num_parts, int num_procs, int seqOption,
                        int paraOption, MPI_Comm comm);
}

#endif
