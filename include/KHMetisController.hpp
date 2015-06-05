
#ifndef _KHMETIS_CONTROLLER_HPP
#define _KHMETIS_CONTROLLER_HPP

// ### KHMetisController.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 30/11/2004: Last Modified
//
// ###

#include "Config.h"
#ifdef LINK_HMETIS

// ## hMeTiS library function declaration ##
extern "C" {
void HMETIS_PartKway(int nv, int nh, int *vwt, int *hptr, int *hv, int *hwt,
                     int p, int ub, int *opt, int *pvector, int *cut);
}

#include "SeqController.hpp"
#include "GreedyKwayRefiner.hpp"

using namespace std;

class KHMetisController : public SeqController {
protected:
  int lenOfOptions;
  int maxPartWt;

  double avePartWt;

  FastDynaArray<int> khMetisOptions;

  GreedyKwayRefiner *kWayRefiner;

public:
  KHMetisController(GreedyKwayRefiner *k, int rank, int nProcs, int nParts,
                    const int *options, ostream &out);
  ~KHMetisController();

  void dispSeqControllerOptions() const;
  void runSeqPartitioner(ParaHypergraph &h, MPI_Comm comm);
};

#endif
#endif
