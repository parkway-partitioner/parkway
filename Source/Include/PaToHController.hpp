
#  ifndef _PATOH_CONTROLLER_HPP
#  define _PATOH_CONTROLLER_HPP


// ### PaToHController.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
// 
// HISTORY: 
// 
// 10/1/2005: Last Modified
//
// ###


#  include "Config.h"
#  ifdef LINK_PATOH
#  include "patoh.h"
#  include "SeqController.hpp"
#  include "GreedyKwayRefiner.hpp"


using namespace std;



class PaToHController
  : public SeqController
{
protected:

  int maxPartWt;
  int patohSettings;

  double avePartWt;
  double bisectConstraint;

  GreedyKwayRefiner *kWayRefiner;

public:

  PaToHController(GreedyKwayRefiner* k, int rank, int nProcs, int nParts, int setting, ostream &out);
  ~PaToHController();

  void dispSeqControllerOptions() const;
  void runSeqPartitioner(ParaHypergraph &h, MPI_Comm comm);


};




#  endif
#  endif
