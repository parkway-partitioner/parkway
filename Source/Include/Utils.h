

#  ifndef _UTILS_H
#  define _UTILS_H


// ### Utils.h ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
// 
// HISTORY: 
// 
// 10/1/2005: Last Modified
//
// ###


#  include <math.h>
#  include <time.h>

#  include "mpi.h"
#  include "Log.h"
#  include "BasicParaController.hpp"
#  include "ParaVCycleFinalController.hpp"
#  include "ParaVCycleAllController.hpp"
#  include "RecurBisectController.hpp"
#  include "KHMetisController.hpp"
#  include "PaToHController.hpp"


namespace Utils
{
  ParaCoarsener *buildParaCoarsener(int myRank, int numProc, int numParts, double constraint, ParaHypergraph *h, ostream &out, const int *options, MPI_Comm comm);
  ParaRestrCoarsener *buildParaRestrCoarsener(int myRank, int numProc, int numParts, double constraint, ParaHypergraph *h, ostream &out, const int *options, MPI_Comm comm);
  ParaRefiner *buildParaRefiner(int myRank, int numProc, int numParts, double constraint, ParaHypergraph *h, ostream &out, const int *options, MPI_Comm comm);
  SeqController *buildSeqController(int myRank, int numProc, int numParts, double constraint, ostream &out, const int *options);
  ParaController *buildParaController(int myRank, int numProcs, int numParts, int num_tot_verts, double constraint, ParaCoarsener *c, ParaRestrCoarsener *rc, ParaRefiner *r, SeqController *s, ostream &out, const int *options, MPI_Comm comm);
  
  void initDefaultValues(const int *userOptions, int *programOptions);
  void checkPartsAndProcs(int num_parts, int num_procs, int seqOption, int paraOption, ostream &out, MPI_Comm comm);
}




#  endif
