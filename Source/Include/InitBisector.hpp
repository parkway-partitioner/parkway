

#  ifndef _INIT_BISECTOR_HPP
#  define _INIT_BISECTOR_HPP


// ### InitBisector.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
// 
// HISTORY: 
// 
// 30/11/2004: Last Modified
//
// ###


#  include "FMRefiner.hpp"


using namespace std;


class InitBisector
  : public FMRefiner
{

protected:

  int numInitRuns;
  
public:

  InitBisector(int nRuns, int insMethod, int ee, int dL);
  ~InitBisector();

  void dispInitBisectorOptions(ostream &out) const;
  void initBisect(Hypergraph &h);  

  int setBaseVertex();
  int chooseBestVertex1to0();
  int doGreedyPass();  

  inline void setNumRuns(register int nR) { numInitRuns = nR; }
};




#  endif
