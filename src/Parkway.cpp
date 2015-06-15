#ifndef _PARKWAY_CPP
#define _PARKWAY_CPP

// ### Parkway.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 12/12/2004: Last Modified
//
// ###

#include "parkway.h"
#include <fstream>
#include <iostream>
#include "data_structures/internal/table_utils.hpp"
#include "hypergraph/hypergraph.hpp"

namespace parallel = parkway::parallel;
namespace ds = parkway::data_structures;

void k_way_partition(const char *file_name, const char *out_file, int num_parts,
                     double constraint, int &k_1cut, const int *options,
                     MPI_Comm comm) {
#ifdef MEM_CHECK
  MemoryTracker::start();
#endif

  char part_file[512];
  char shuffle_file[512];
  char message[512];

  int init_options[26];
  int disp_option;
  int output_partition_tofile;
  int num_procs;
  int my_rank;

#ifdef USE_SPRNG
  int seed;
#endif

  std::ostream *output;

  parallel::hypergraph *hgraph = nullptr;
  parallel_coarsener *coarsener = nullptr;
  parallel_restrictive_coarsening *restrC = nullptr;
  parallel_refiner *refiner = nullptr;
  sequential_controller *seqController = nullptr;
  parallel_controller *controller = nullptr;

  ds::internal::table_utils tableUtils;

  MPI_Comm_size(comm, &num_procs);
  MPI_Comm_rank(comm, &my_rank);

  if (out_file) {
    output = new std::ofstream(out_file,
                               std::ofstream::app | std::ofstream::out);

    if (!output->good()) {
      sprintf(message, "p[%d] could not open file %s - abort\n", my_rank,
              out_file);
      cout << message;
      MPI_Abort(comm, 0);
    }
  } else {
    output = &std::cout;
  }

  Utils::initDefaultValues(options, init_options);
  Utils::checkPartsAndProcs(num_parts, num_procs, init_options[15],
                            init_options[22], *output, comm);

  disp_option = init_options[2];
  output_partition_tofile = init_options[3];

/* init pseudo-random number generator */

#ifdef USE_SPRNG
  if (init_options[1] == 0) {
    seed = make_sprng_seed();
  } else {
    seed = init_options[1];
  }

  if (!init_sprng(0, seed, SPRNG_DEFAULT)) {
    sprintf(
        message,
        "p[%d] could not initialise sprng random number generator - abort\n",
        my_rank);
    *output << message;
    MPI_Abort(comm, 0);
  }
#else
  if (init_options[1] == 0)
    srand48((my_rank + 1) * RAND_SEED);
  else
    srand48(init_options[1]);
#endif

  sprintf(part_file, "%s.part.%d", file_name, num_parts);
  sprintf(shuffle_file, "%s.part.%d", file_name, num_procs);

  if (my_rank == 0 && disp_option > 0) {
    Funct::printIntro(*output);
  }

  hgraph = new parallel::hypergraph(my_rank, num_procs, file_name, disp_option,
                              *output, comm);

  if (!hgraph) {
    sprintf(message,
            "p[%d] not able to build local hypergraph from %s - abort\n",
            my_rank, file_name);
    *output << message;
    MPI_Abort(comm, 0);
  }

  ds::internal::table_utils::set_scatter_array(hgraph->total_number_of_vertices());

  coarsener =
      Utils::buildParaCoarsener(my_rank, num_procs, num_parts, constraint,
                                hgraph, *output, init_options, comm);
  restrC =
      Utils::buildParaRestrCoarsener(my_rank, num_procs, num_parts, constraint,
                                     hgraph, *output, init_options, comm);
  refiner = Utils::buildParaRefiner(my_rank, num_procs, num_parts, constraint,
                                    hgraph, *output, init_options, comm);
  seqController = Utils::buildSeqController(my_rank, num_procs, num_parts,
                                            constraint, *output, init_options);

  if (!coarsener) {
    sprintf(message, "p[%d] not able to build ParaCoarsener - abort\n",
            my_rank);
    *output << message;
    MPI_Abort(comm, 0);
  }

  if (!refiner) {
    sprintf(message, "p[%d] not able to build ParaRefiner - abort\n", my_rank);
    *output << message;
    MPI_Abort(comm, 0);
  }

  if (!seqController) {
    sprintf(message, "p[%d] not able to build SeqController - abort\n",
            my_rank);
    *output << message;
    MPI_Abort(comm, 0);
  }

  controller = Utils::buildParaController(
      my_rank, num_procs, num_parts, hgraph->total_number_of_vertices(), constraint,
      coarsener, restrC, refiner, seqController, *output, init_options, comm);

  if (!controller) {
    sprintf(message, "p[%d] not able to build ParaController - abort\n",
            my_rank);
    *output << message;
    MPI_Abort(comm, 0);
  }

  if (my_rank == 0 && disp_option > 0) {
    Funct::printEnd(*output);
  }

  hgraph->compute_balance_warnings(num_parts, constraint, *output, comm);

  controller->set_hypergraph(hgraph);
  controller->set_prescribed_partition(shuffle_file, comm);
  controller->set_weight_constraints(comm);
  controller->run(comm);

  if (output_partition_tofile) {
    controller->partition_to_file(part_file, comm);
  }

  k_1cut = controller->best_cut_size();

  dynamic_memory::delete_pointer<parallel_controller>(controller);
  dynamic_memory::delete_pointer<sequential_controller>(seqController);
  dynamic_memory::delete_pointer<parallel_restrictive_coarsening>(restrC);
  dynamic_memory::delete_pointer<parallel_coarsener>(coarsener);
  dynamic_memory::delete_pointer<parallel_refiner>(refiner);
  dynamic_memory::delete_pointer<parallel::hypergraph>(hgraph);

  if (out_file) {
    dynamic_memory::delete_pointer<std::ostream>(output);
  }

#ifdef MEM_CHECK
  MemoryTracker::stop();
#endif
}

void k_way_partition(int numVertices, int numHedges, const int *vWeights,
                     const int *hEdgeWts, const int *offsets,
                     const int *pinList,
                     int numParts, double constraint, int &k_1cut,
                     const int *options, int *pVector, const char *out_file,
                     MPI_Comm comm) {
  int num_procs;
  int my_rank;
  int maxHedgeLen = 0;
  int globMaxHedgeLen;
  int init_options[26];
  int disp_option;

#ifdef USE_SPRNG
  int seed;
#endif

  char message[512];
  std::ostream *output;

  int i;
  int j;

  parallel::hypergraph *hgraph = nullptr;
  parallel_coarsener *coarsener = nullptr;
  parallel_restrictive_coarsening *restrC = nullptr;
  parallel_refiner *refiner = nullptr;
  sequential_controller *seqController = nullptr;
  parallel_controller *controller = nullptr;

  ds::internal::table_utils tableUtils;

  MPI_Comm_size(comm, &num_procs);
  MPI_Comm_rank(comm, &my_rank);

  if (out_file) {
    output = new std::ofstream(out_file,
                               std::ofstream::app | std::ofstream::out);

    if (!output->good()) {
      sprintf(message, "p[%d] could not open file %s - abort\n", my_rank,
              out_file);
      std::cout << message;
      MPI_Abort(comm, 0);
    }
  } else
    output = &std::cout;

  Utils::initDefaultValues(options, init_options);
  Utils::checkPartsAndProcs(numParts, num_procs, init_options[15],
                            init_options[22], *output, comm);

  disp_option = init_options[2];

/* init pseudo-random number generator */

#ifdef USE_SPRNG
  if (init_options[1] == 0)
    seed = make_sprng_seed();
  else
    seed = init_options[1];

  if (!init_sprng(0, seed, SPRNG_DEFAULT)) {
    sprintf(
        message,
        "p[%d] could not initialise sprng random number generator - abort\n",
        my_rank);
    *output << message;
    MPI_Abort(comm, 0);
  }
#else
  if (init_options[1] == 0)
    srand48((my_rank + 1) * RAND_SEED);
  else
    srand48(init_options[1]);
#endif

  for (i = 0; i < numHedges; ++i) {
    j = offsets[i + 1] - offsets[i];

    if (j > maxHedgeLen)
      maxHedgeLen = j;
  }

  MPI_Allreduce(&maxHedgeLen, &globMaxHedgeLen, 1, MPI_INT, MPI_MAX, comm);

  if (my_rank == 0 && disp_option > 0)
    Funct::printIntro(*output);

  hgraph = new parallel::hypergraph(my_rank, num_procs, numVertices, numHedges,
                              globMaxHedgeLen, vWeights, hEdgeWts, pinList,
                              offsets, disp_option, *output, comm);

  if (!hgraph) {
    sprintf(message, "p[%d] could not initialise hypergraph - abort\n",
            my_rank);
    *output << message;
    MPI_Abort(comm, 0);
  }

  ds::internal::table_utils::set_scatter_array(hgraph->total_number_of_vertices());

  coarsener =
      Utils::buildParaCoarsener(my_rank, num_procs, numParts, constraint,
                                hgraph, *output, init_options, comm);
  restrC =
      Utils::buildParaRestrCoarsener(my_rank, num_procs, numParts, constraint,
                                     hgraph, *output, init_options, comm);
  refiner = Utils::buildParaRefiner(my_rank, num_procs, numParts, constraint,
                                    hgraph, *output, init_options, comm);
  seqController = Utils::buildSeqController(my_rank, num_procs, numParts,
                                            constraint, *output, init_options);

  if (!coarsener) {
    sprintf(message, "p[%d] not able to build ParaCoarsener - abort\n",
            my_rank);
    *output << message;
    MPI_Abort(comm, 0);
  }

  if (!refiner) {
    sprintf(message, "p[%d] not able to build ParaRefiner - abort\n", my_rank);
    *output << message;
    MPI_Abort(comm, 0);
  }

  if (!seqController) {
    sprintf(message, "p[%d] not able to build SeqController - abort\n",
            my_rank);
    *output << message;
    MPI_Abort(comm, 0);
  }

  controller = Utils::buildParaController(
      my_rank, num_procs, numParts, hgraph->total_number_of_vertices(), constraint,
      coarsener, restrC, refiner, seqController, *output, init_options, comm);

  if (!controller) {
    sprintf(message, "p[%d] not able to build ParaController - abort\n",
            my_rank);
    *output << message;
    MPI_Abort(comm, 0);
  }

  if (my_rank == 0 && disp_option > 0)
    Funct::printEnd(*output);

  hgraph->compute_balance_warnings(numParts, constraint, *output, comm);

  controller->set_hypergraph(hgraph);
  controller->set_weight_constraints(comm);
  controller->run(comm);
  controller->copy_out_partition(numVertices, pVector);

  k_1cut = controller->best_cut_size();

  dynamic_memory::delete_pointer<parallel_controller>(controller);
  dynamic_memory::delete_pointer<sequential_controller>(seqController);
  dynamic_memory::delete_pointer<parallel_restrictive_coarsening>(restrC);
  dynamic_memory::delete_pointer<parallel_coarsener>(coarsener);
  dynamic_memory::delete_pointer<parallel_refiner>(refiner);
  dynamic_memory::delete_pointer<parallel::hypergraph>(hgraph);

  if (out_file)
    dynamic_memory::delete_pointer<std::ostream>(output);
}

#endif
