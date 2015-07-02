#ifndef _MACROS_H
#define _MACROS_H
// ### Macros.h ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 12/1/2005: Last Modified
//
// ###

#include "configurtion.h"

/* debug options */

#ifdef DEBUG_ALL
#define DEBUG_BASICS
#define DEBUG_COARSENER
#define DEBUG_REFINER
#define DEBUG_FM_REFINER
#define DEBUG_TABLES
#define DEBUG_HYPERGRAPH
#define DEBUG_CONTROLLER
#define DEBUG_LOADER
#endif

/* hash key parameters */

#define KEY_UINT
#ifdef KEY_UINT
#define HashKey unsigned int
#define SLIDE1 0
#define SLIDE2 16
#else
#define HashKey unsigned long long
#define SLIDE1 16
#define SLIDE2 48
#endif

// ### display options

#define SILENT 0
#define OPT_ONLY 1
#define DISP_ALL 2

// # random seed

#define RAND_SEED 117

// # vertex visit order

#define INCREASING_ORDER 1
#define DECREASING_ORDER 2
#define RANDOM_ORDER 3
#define INCREASING_WEIGHT_ORDER 4
#define DECREASING_WEIGHT_ORDER 5

// # type of FCCoarsening

#define FCwithFanOut 1
#define FCwithoutFanOut 2
#define FCwithFanOutDiv 3
#define FCwithoutFanOutDiv 4

// # type of RestrFCCoarsening

#define RestrFCwithFanOut 1
#define RestrFCwithoutFanOut 2
#define RestrFCwithFanOutDiv 3
#define RestrFCwithoutFanOutDiv 4

// # ParaFCCoarsening

#define ParaFCC 3
#define ParaRestFCC 4
#define ParaApproxFCC 5
#define Para2DModel 6

// # min nodes before init partitioning phase

#define MIN_VERT_MULTIPLIER 80
#define MIN_PARA_VERT_MULTIPLIER 200

// # default reduction ratios

#define DEF_PARA_RATIO 1.75
#define DEF_BIS_RATIO 1.75
#define DEF_SEQ_RATIO 1.75

#define INC 0
#define DEC 1
#define ROOT_PROC 0

// predefined threshold values
// for anticipated number of entries

#define VSMALL 2000
#define SMALL 20000
#define MEDIUM 125000
#define LARGE 500000
#define VLARGE 1250000

// predefined hash table sizes
// based on anticipated number of entries

#define VSMALL_TABLE 971
#define SMALL_TABLE 5081
#define MEDIUM_TABLE 25057
#define LARGE_TABLE 120097
#define VLARGE_TABLE 555109
#define HUGE_TABLE 1294031

// auxiliary during coarsening/refinement

#define LARGE_CONSTANT 0xFFFFFFF
#define MATCHED_LOCALLY 0xFFFFFFF
#define NON_LOCAL_MATCH 50000000
#define NO_MATCH 50000001
#define MIN_ALLOWED_REDUCTION_RATIO 1.0005
#define RESP_FOR_HEDGE -1
#define HEDGE_TABLE_REDUC_FACTOR 0.75

// fm bucket instert order

#define LIFO 0
#define FIFO 1

// reduction in acceptable cut during refinement

#define DEF_REDUC_FACTOR 0.5
#define DEF_PARA_REDUC_FACTOR 0.7

// default balance constraint

#define DEF_BAL_CONSTRAINT 0.05

// default numbers of partitioning runs

#define DEF_BISECT_RUNS 20
#define DEF_SEQ_INIT_RUNS 20
#define DEF_PARA_RUNS 1
#define DEF_SEQ_RUNS 20

// early exit fm parameter

#define DEF_EE_PARAM 5
#define DEF_PARA_EELIMIT 1.0
#define DEF_SEQ_KWAY_LIM 1.0

// default for keep threshold percentage

#define DEF_KEEP_THRESHOLD 0.05
#define DEF_PARA_KEEP_THRESHOLD 0.05

// Para refiner types

#define ParaGreedyKway 1

// Para controller types

#define BasicParaC 1
#define ParaVCycleBig 2
#define ParaVCycleAll 3

// Seq Controller types

#define RecurBisect 1
#define VCycleFinalRecurBisect 2
#define VCycleAllRecurBisect 3
#define KhMeTiS 4
#define PaToH 5
//#  define Kway 6

// khmetis options

// - coarsening
#define HFC_COARSE 1
#define FC_COARSE 2
#define GFC_COARSE 3
#define HEDGE_COARSE 4
#define EDGE_COARSE 5

#define NO_DEBUG 0
#define RANDOM_SEED -1
// - uncoarsening options
#define NO_VCYCLE 0
#define VCYCLE_FIN 1
#define VCYCLE_BEST_INTER 2
#define VCYCLE_EACH_INTER 3

#define USER_DEFINED 1
#define SOED 2

#endif
