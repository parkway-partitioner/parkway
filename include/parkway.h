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
#include "options.hpp"
#include "mpi.h"

int k_way_partition(const parkway::options &options, MPI_Comm comm);

#endif
