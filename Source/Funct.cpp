
#  ifndef _FUNCT_CPP
#  define _FUNCT_CPP


// ### Funct.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY: 
//
// 12/1/2005: Last Modified
//
// ###


#  include "Funct.hpp"


char *Funct::startHeap = 0;
int Funct::shift1 = 7;
int Funct::shift2 = 13;
int Funct::maxHedgeLen = 0;


Funct::Funct()
{
  
}


Funct::~Funct()
{

}



int Funct::setTableSize(register int approxNumElem) 
{
  if(approxNumElem < VSMALL)    
    return VSMALL_TABLE;    
  
  if(approxNumElem < SMALL)    
    return SMALL_TABLE;    
  
  if(approxNumElem < MEDIUM)    
    return MEDIUM_TABLE;    
  
  if(approxNumElem < LARGE)    
    return LARGE_TABLE;    
  
  if(approxNumElem < VLARGE)    
      return VLARGE_TABLE;    
  
  else    
    return HUGE_TABLE;    
}


int Funct::getParameterAsInteger(int argc, char **argv, const char *cmpr, int def) 
{
  register int i;
  register int j = argc-1;

  for(i=0;i<j;++i) 
    {
      if(strcmp(cmpr,argv[i]) == 0) 
	return atoi(argv[i+1]);
    }

  return def;
}



double Funct::getParameterAsDouble(int argc, char **argv, const char *cmpr, double def)
{
  register int i;
  register int j = argc-1;
  
  for (i=0;i<j;++i)
    {
      if(strcmp(cmpr,argv[i]) == 0)
	return atof(argv[i+1]);
    }

  return def;
}


char* Funct::getParameterAsCharPtr(int argc, char **argv, const char *cmpr, char *def)
{
  register int i;
  register int j = argc-1;
  
  for (i=0;i<j;++i)
    {
      if(strcmp(cmpr,argv[i]) == 0)
	return argv[i+1];
    }

  return def;
}



double Funct::toRecurBal(double e, int nP) 
{
  double a;
  double b;
  
  if(isPowerOf2(nP))
    {
      a = (1.0+e) / nP;
      b = 1.0 / log2(nP);            
    }
  else
    {
      /* give the balance for nearest 
	 power of 2 above nP          */
      
      int pow2 = 1;
      int index = 1;
      
      while(index < nP)
	{
	  index = Shiftl(index,1);
	  ++pow2;
	}
      
      a = (1.0+e) / index;
      b = 1.0 / pow2;
    }
  
  return (pow(a,b)-0.5);
}



HashKey Funct::computeHash(const int *vs, int len) 
{
  HashKey key = 0;

  register unsigned int slide1 = SLIDE1;
  register unsigned int slide2 = SLIDE2;

  int sum = 0;
  int i;

  for (i=0;i<maxHedgeLen;++i) 
    {
      if(i < len) 
	{	
	  sum += vs[i];
	  key = Xor(key, RotateLeft(vs[i], slide1));
	}
      else 
	{	
	  sum += 1;
	  key = Xor(key, RotateLeft(1, slide1));
	}
      
      key = Xor(key, RotateLeft(sum, slide2));
      
      slide1 = Mod((slide1 + shift1), sizeof(HashKey));
      slide2 = Mod((slide2 + shift2), sizeof(HashKey));
    }
  
  return key;
}


void Funct::randomPermutation(int *array, int size) 
{
  register int i;
  register int j;
  register int ij;

  for (i=0;i<size;++i)           
    {
      ij = RANDOM(i,size);
      fswap(array[i],array[ij],j);
    }    
}


void Funct::qsortByAnotherArray(const int left, const int right, int *array, const int *valArray, const int order) 
{
  register int left_arrow = left;
  register int right_arrow = right;

  int pivot = array[(left+right)/2];  

  if(order == INC)
    {      
      do 
	{
	  while(valArray[array[right_arrow]] > valArray[pivot])
	    --right_arrow;
	  while(valArray[array[left_arrow]] < valArray[pivot])
	    ++left_arrow;
	  
	  if(left_arrow <= right_arrow) 
	    {
	      swap(array[left_arrow++],array[right_arrow--]);	      
	    }
	} 
      while (right_arrow >= left_arrow);
      
      if(left < right_arrow)
	qsortByAnotherArray(left, right_arrow, array, valArray, order);

      if(left_arrow < right)
	qsortByAnotherArray(left_arrow, right, array, valArray, order);  
    }
  else
    {
      do 
	{
	  while(valArray[array[right_arrow]] < valArray[pivot])
	    --right_arrow;
	  while(valArray[array[left_arrow]] > valArray[pivot])
	    ++left_arrow;
	
	  if(left_arrow <= right_arrow) 
	    {
	      swap(array[left_arrow++],array[right_arrow--]);
	    }
	} 
      while (right_arrow >= left_arrow);
      
      if(left < right_arrow)
	qsortByAnotherArray(left, right_arrow, array, valArray, order);

      if(left_arrow < right)
	qsortByAnotherArray(left_arrow, right, array, valArray, order);  
    }
}



void Funct::printIntro(ostream &out)
{
  out << endl
      << " ------- PARKWAY2.0 -------" << endl
      << "|" << endl;
}



void Funct::printEnd(ostream &out) 
{
  out << " --------------------------" << endl << endl;
}


void Funct::printMemUse(int myRank, const char *info)
{  
  write_log(myRank, "Mem use at point: %s", info);

  char filename[512];
  char word[1024];
  int pid = getpid();
  int next = 0;

  sprintf(filename, "/proc/%d/status", pid);
  
  FILE *fp = fopen(filename, "r");
  assert(fp);

  fscanf(fp, " %s ", word);
  
  while (!feof(fp))
    {
      if(next) 
	{
	  write_log(myRank, "using %s kB", word);
	  break;
	}

      if(strstr( word, "VmRSS"))	 
	next = 1;
      
      fscanf(fp, " %s ", word);
    }

  fclose(fp);
}


void Funct::initStartMem()
{
  startHeap = (char*)(sbrk(0));
}


/*
void Funct::write_to_memory() 
{
  int x;
  int y;
  int z;
  int **block;
  
  register int i;

  block = new (int*)[1024];

  for(i=0; i<1024; ++i) 
    {
      block[i] = new int[2];      
    }

  for(i=0; i<1024 * 2; ++i) 
    {
      x = static_cast<int>(drand48() * 1024);
      y = static_cast<int>(drand48() * 2);
      z = static_cast<int>(drand48() * 10);
      
      block[x][y] = z;  
    }
  
  for(i=0; i<1024; ++i)
    delete [] block[i];
  delete [] block;
}

*/


#endif
