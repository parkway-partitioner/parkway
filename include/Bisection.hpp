
#ifndef _BISECTION_HPP
#define _BISECTION_HPP

// ### Bisection.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 30/11/2004: Last Modified
//
// ###

#include "Hypergraph.hpp"

using namespace std;

class Bisection {

protected:
  int bisectAgain;
  int numVertices;
  int partID;

  Hypergraph *hypergraph;

  dynamic_array<int> mapToOrigVerts;

public:
  Bisection(Hypergraph *h, int bAgain, int pID) {
    hypergraph = h;

    numVertices = h->getNumVertices();
    bisectAgain = bAgain;
    partID = pID;
  }

  inline void initMap() {
    int i;

    mapToOrigVerts.setLength(numVertices);

    for (i = 0; i < numVertices; ++i)
      mapToOrigVerts[i] = i;
  }

  inline void setMap(int *vMap, int nV) {
    mapToOrigVerts.setArray(vMap, nV);
  }

  inline int getPartID() const { return partID; }
  inline int getBisectAgain() const { return bisectAgain; }
  inline int getNumVertices() const { return numVertices; }
  inline int *getMapArray() const { return mapToOrigVerts.getArray(); }

  inline Hypergraph *getHypergraph() const { return hypergraph; }
};

#endif
