
#  ifndef _PARA_APPROX_COARSENER_HPP
#  define _PARA_APPROX_COARSENER_HPP


// ### ParaApproxCoarsener.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
// 
// HISTORY: 
//
// 31/12/2004: Last Modified
//
// ###


#  include "Bit.hpp"
#  include "ParaCoarsener.hpp"


using namespace std;


class ParaApproxCoarsener  
  : public ParaCoarsener 
{

protected:

  int startPercentile;
  int currPercentile;
  int increment;
  
public:

  ParaApproxCoarsener(int _rank, int _numProcs, int _numParts, int percentile, int inc, ostream &out);
  
  virtual ~ParaApproxCoarsener();
  virtual ParaHypergraph *coarsen(ParaHypergraph &h, MPI_Comm comm)=0;
  virtual void setClusterIndices(MPI_Comm comm)=0;
  virtual void releaseMemory()=0;
  virtual void dispCoarseningOptions() const=0;
  virtual void buildAuxiliaryStructs(int numTotPins, double aveVertDeg, double aveHedgeSize)=0;
  
  void loadHyperGraph(const ParaHypergraph &h, MPI_Comm comm);
  void computeHedgesToLoad(BitField &toLoad, int numH, int *hEdgeWts, int *hEdgeOffsets, MPI_Comm comm);

};



#  endif