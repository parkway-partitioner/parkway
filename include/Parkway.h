#ifndef _PARKWAY_H
#define _PARKWAY_H

// ### parkway.h ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 03/12/2004: Last Modified
//
// ###

#include "Utils.h"

namespace parkway {

// Partition hypergraph stored on disk.
void ParaPartKway(
    const char *file_name, // base name of input hypergraph
    const char *out_file,  // name of output file (null for cout)
    int num_parts,         // number of parts in partition
    double constraint,     // balance constraint
    int &k_1cut,           // cutsize of best partition
    const int *options,    // user defined options
    MPI_Comm comm);

// Partition hypergraph stored in memory.
void ParaPartKway(
    int numVertices,      // number of vertices stored by process
    int numHedges,        // number of hyperedges
    const int *vWeights,  // weights of locally stored vertices
    const int *hEdgeWts,  // weights of locally stored hyperedges
    const int *offsets,   // offsets of hyperedges stored in pinlist; hyperedge
                          // i is stored between index offsets[i] and
                          // index offsets[i+1] in the pinList
    const int *pinList,   // list of vertex indices of hyperedges
    int numParts,         // number of parts in partition
    double constraint,    // balance constraint
    int &k_1cut,          // cutsize of best partition
    const int *options,   // user defined options
    int *pVector,         // computed part values of local vertices
    const char *outfile,  // name of output file (null for cout)
    MPI_Comm comm);
}  // namesapce parkway

#endif
