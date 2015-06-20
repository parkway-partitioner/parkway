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

void k_way_partition(const char *file_name, const char *out_file, int num_parts,
                     double constraint, int &k_1cut, const int *options,
                     MPI_Comm comm);
void k_way_partition(int numVertices, int numHedges, const int *vWeights,
                     const int *hEdgeWts, const int *offsets,
                     const int *pinList,
                     int numParts, double constraint, int &k_1cut,
                     const int *options, int *pVector, const char *outfile,
                     MPI_Comm comm);

#endif
