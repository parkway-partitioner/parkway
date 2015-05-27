
#  ifndef _LOG_CPP
#  define _LOG_CPP


// ### Log.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, William Knottenbelt, Imperial College London
// 
// HISTORY: 
// 
// 30/11/2004: Last Modified
//
// ###


#  include <cstdio>
#  include <cstdarg>
#  include <ctime>
#  include <unistd.h>
#  include <iostream>
#  include <cassert>
#  include "mpi.h"
#  include "Funct.hpp"


using namespace std;


void write_log(int rank,const char *fmt, ...) 
{
  char buffer[2048];
  char output[4096];
  char logfile[1024];
  char hostname[1024];
  char curr_path[512];
  char *t;

  long seconds;
  double start = MPI_Wtime();

  va_list argp;
  va_start(argp, fmt);
  vsprintf(buffer, fmt, argp);
  va_end(argp);
  
  seconds = time(NULL);
  
  t = ctime(&seconds);
  t[19]='\0';
  t += 11;

  gethostname(hostname, 1024);
  sprintf(curr_path,".");
  sprintf(logfile, "%s/log-%02d",curr_path,rank);

  FILE *fp = fopen(logfile, "a");
  assert(fp);

  sprintf(output, "-%s- %s [%d] %s [%5.3g]", hostname, t, rank, buffer, MPI_Wtime() - start);

  cout << output << endl;
  cout.flush();
  
  fprintf(fp, "%s\n", output);
  fclose(fp);
}      
  
                                                                             


#  endif

