
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

#include "hypergraph.hpp"

using namespace std;

class Bisection {

protected:
  int bisect_again_;
  int number_of_vertices_;
  int part_id_;

  hypergraph *hypergraph_;

  dynamic_array<int> mapToOrigVerts;

public:
  Bisection(hypergraph *h, int bAgain, int pID) {
    hypergraph_ = h;

    number_of_vertices_ = h->number_of_vertices();
    bisect_again_ = bAgain;
    part_id_ = pID;
  }

  inline void initMap() {
    int i;

    mapToOrigVerts.reserve(number_of_vertices_);

    for (i = 0; i < number_of_vertices_; ++i)
      mapToOrigVerts[i] = i;
  }

  inline void setMap(int *vMap, int nV) {
    mapToOrigVerts.set_data(vMap, nV);
  }

  inline int getPartID() const { return part_id_; }
  inline int getBisectAgain() const { return bisect_again_; }
  inline int getNumVertices() const { return number_of_vertices_; }
  inline int *getMapArray() const { return mapToOrigVerts.data(); }

  inline hypergraph *getHypergraph() const { return hypergraph_; }
};

#endif
