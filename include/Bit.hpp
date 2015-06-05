
#  ifndef _BITFIELD_HPP
#  define _BITFIELD_HPP


// ### BitField.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
// Modified from code by William Knottenbelt
// 
// HISTORY: 
// 
// 30/11/2004: Last Modified
//
// NOTES:
//
//  - Assume that sizeof(unsigned int) == 4 bytes
//  
// ###



#  include <iostream>
#  include "Dyna.hpp"


using namespace std;


class BitField 
{
  int bitLength;
  int length;
  
  FastDynaArray<unsigned int> data;
  
public:

  BitField() { }
  BitField(register int bits) 
  {
    setLength(bits);
  }

  ~BitField() { }  

  inline void check(int bit) 
  {
    int old = data.getLength();
    int bits = bit + 1;
    int target = (bits >> 5) + ( (bits & 31) ? 1 : 0 );
    
    data.check(target);
    length = data.getLength();
    
    if (length > old) 
      {
	for (register int n = old; n<length; ++n)	  
	  data[n] = 0;	  
      }
    bitLength = (length << 5);
  }
  
  inline int getNumBits() const { return bitLength; }
  inline int getLength() const { return length; }
  
  inline unsigned int *getData() const { return data.getArray(); }
  inline unsigned int getChunk(register int i) const { return data[i]; }
  
  inline void setLength(register int bits) 
  {
    bitLength = bits;
    length = (bits >> 5) + ( (bits & 31) ? 1 : 0 );
    data.setLength(length);
  }

  inline void clear() 
  {
    for (register int l=0; l<length; ++l)      
      data[l] = 0;    
  }

  inline void set1() 
  {
    for (register int l=0;l<length;++l)       
      data[l] = 0xFFFFFFFF;      
  } 

  inline int allSet() 
  {
    register int allset=1;
    for (register int i=0;i<length;++i)      
      if(data[i] != 0xFFFFFFFF) 	  
	return 0;           
    return allset;
  }
  
  inline int operator()(register int index) const 
  {
    return (data(index>>5) & (1 << (index & 31))) ? 1 : 0;
  }
  
  inline void set1(register int index) 
  {
    data[index >> 5] |= (1 << (index & 31));
  }

  inline void set0(register int index) 
  {
    data[index >> 5] &= (~(1 << (index & 31)));
  }

};




#  endif
