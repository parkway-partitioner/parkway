
#  ifndef _STRING_UTILS_CPP
#  define _STRING_UTILS_CPP


// ### StringUtils.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
// 
// HISTORY: 
// 
// 12/1/2005: Last Modified
//
// ###


#  include "StringUtils.hpp"


StringUtils::StringUtils()
{


}


StringUtils::~StringUtils()
{


}


void StringUtils::skipNonDigits(register char* &strPtr, register char comm) 
{ 
  while(*strPtr != '\0' && *strPtr != comm && (*strPtr < '0' || *strPtr > '9'))      
    ++strPtr;    
  
  if(*strPtr == comm)
    *strPtr = '\0';
}



void StringUtils::skipNonDigits(register char * &strPtr)  
{
  while(*strPtr != '\0' && (*strPtr < '0' || *strPtr > '9'))  
    ++strPtr;
}



void StringUtils::skipDigits(register char* &strPtr) 
{ 
  while(*strPtr != '\0' && *strPtr >= '0' && *strPtr <= '9')  
    ++strPtr;
}



int StringUtils::stringToDigit(register char* &strPtr)
{
  register int number = 0;
  
  while(*strPtr >= '0' && *strPtr <= '9')
    {
      number *= 10;
      number += (*strPtr++ - '0');
    }
  
  return number;
}


int StringUtils::getParameterAsInteger(int argc, char **argv, const char *cmpr)
{
  register int i;
  register int j = argc-1;

  for(i=0;i<j;++i) 
    {
      if(strcmp(cmpr,argv[i]) == 0) 
	return atoi(argv[i+1]);
    }

  return -1;
}



int StringUtils::getParameterAsInteger(int argc, char **argv, const char *cmpr, int def)
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



char *StringUtils::getParameterAsCharPtr(int argc, char **argv, const char *cmpr, char *def)
{
  register int i;
  register int j = argc-1;

  for(i=0;i<j;++i) 
    {
      if(strcmp(cmpr,argv[i]) == 0) 
	return argv[i+1];
    }

  return def;
}


char StringUtils::getFirstChar(char *c)
{
  while(*c == ' ' && *c != '\0')
    ++c;
  
  return *c;
}



#  endif
