#  ifndef _HARWELLBOEING_TOBIN_CPP
#  define _HARWELLBOEING_TOBIN_CPP


// ### HarwellBoeingToBin.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
// 
// HISTORY: 
// 
// 30/12/2004: Last Modified
//
// ###


#  include "HarwellBoeingToBin.hpp"


HarwellBoeingToBin::HarwellBoeingToBin(int model)
  : HarwellBoeingReader()
{
  hypergraphModelType = model;
}


HarwellBoeingToBin::~HarwellBoeingToBin()
{

}



void HarwellBoeingToBin::convert(const char *filename)
{
  /* assume that this is a general, non-symmetric sparse matrix */

  ifstream in_stream;
  ofstream out_stream;
  
  char bin_file[512];
  char *data;

  int binPreAmble[3];
  int endOffset;
  int hEdgeLengthSlot;
  int maxIntsWritten = TextFileReader::getMaxPinsInChunk();
  
  FastDynaArray<int> hEdgeData;
  FastDynaArray<int> hEdgeOffsets;
  FastDynaArray<int> hEdgeWeights;
  FastDynaArray<int> pinList;
  FastDynaArray<int> vWeights;

  register int i;
  register int j;
  register int ij;

  in_stream.open(filename, ifstream::in);
  
  if(!in_stream.is_open())
    {
      cout << "error opening " << filename << endl;
      exit(1);
    }

  if(hypergraphModelType == 1)
    sprintf(bin_file, "%s.mat.bin", filename);
  else
    sprintf(bin_file, "%s.bin", filename);

  out_stream.open(bin_file, ofstream::out | ofstream::app | ofstream::binary);
  
  if(!out_stream.is_open())
    {
      cout << "error opening " << bin_file << endl;
      in_stream.close();
      exit(1);
    }
    
  readPreamble(in_stream);

  pinList.setLength(numNonZeros);
  hEdgeOffsets.setLength(numCols+1);
  hEdgeWeights.setLength(numCols);
  vWeights.setLength(numRows);

  for (i=0;i<numCols;++i)
    hEdgeWeights[i] = 1;

  /* Initialise the hyperedge offsets within pin-list */

  j=0;
  for (i=0;i<totColPtrLines;++i)
    {
      getLine(in_stream);
      data = buffer.getArray();
     
      while(*data != '\0')
	{
	  StringUtils::skipNonDigits(data);
	  hEdgeOffsets[j++] = StringUtils::stringToDigit(data)-1;
	}
    }
  
  if(j != numCols+1)
    {
      cout << "could not find " << numCols+1 << " col ptr entries" << endl;
      cout << "instead found " << j << " col ptr entries" << endl;

      in_stream.close();
      out_stream.close();
      exit(1);
    }

  /* Initialise the pin-list */

  j=0;
  for (i=0;i<totRowIdxLines;++i)
    {
      getLine(in_stream);
      data = buffer.getArray();
      
      while(*data != '\0')
	{
	  StringUtils::skipNonDigits(data);
	  pinList[j++] = StringUtils::stringToDigit(data)-1;
	}   
    }

  if(j != numNonZeros)
    {
      cout << "could not find " << numNonZeros << " non-zero entries" << endl;
      cout << "instead found " << j << " non-zero entries" << endl;

      in_stream.close();
      out_stream.close();
      exit(1);
    }

  /* check num zeros on main diagonal */ 
  
  int numZerosOnMainDiag = 0;
  int hasNonZeroOnDiag;

  BitField addNonZero(numCols);
  addNonZero.clear();

  for (i=0;i<numCols;++i)
    {    
      hasNonZeroOnDiag = 0;
      
      for (j=hEdgeOffsets[i];j<hEdgeOffsets[i+1];++j)		  
	if(pinList[j] == i)
	  hasNonZeroOnDiag = 1;
      
      if(hasNonZeroOnDiag == 0)
	{
	  ++numZerosOnMainDiag;
	  addNonZero.set1(i);
	}
    }

  cout << "num zeros on main diagonal = " << numZerosOnMainDiag << endl;

  /* compute vertex weights */
  
  if(hypergraphModelType == 1)
    {
      for (i=0;i<numRows;++i)
	vWeights[i] = 0;
      
      /* set weight to sum non-zeros in row */
      
      for (i=0;i<numNonZeros;++i)
	++vWeights[pinList[i]];

      /* NB - NO NEED TO INC WEIGHT AS NON-ZEROS ONLY ADDED IN AS 
	 MODELLING NECESSITY AND WILL NOT CONTRIBUTE TO ACTUAL 
	 COMPUTATIONAL LOAD */
      /* increment weight to take into account added non-zeros on diagonal */
      /*
      for (i=0;i<numRows;++i)
	if(addNonZero(i) == 1)
	  ++vWeights[i];
      */
    }
  else
    {    
      /* set the weight to 1 */
      
      for (i=0;i<numRows;++i)
	vWeights[i] = 1;
    }  
  
  /* output hypergraph on file */

  if(hypergraphModelType == 1)
    cout << "will add non-zeros on main diagonal" << endl;

  binPreAmble[0] = numRows;
  binPreAmble[1] = numCols;
 
  if(hypergraphModelType == 1)
    binPreAmble[2] = numNonZeros+numZerosOnMainDiag;
  else
    binPreAmble[2] = numNonZeros;

  cout << "numVertices = " << numRows 
       << ", numHedges = " << numCols
       << ", numPins = " << binPreAmble[2] << endl; 

  out_stream.write((char*)(&binPreAmble), sizeof(int)*3);

  ij = 0;
  for (i=0;i<numCols;++i)
    {
      if(ij == 0)	
	hEdgeData.assign(ij++, -1);	
      
      hEdgeLengthSlot = ij;
      endOffset = hEdgeOffsets[i+1];

      hEdgeData.assign(ij++, -1);
      hEdgeData.assign(ij++, hEdgeWeights[i]);

      for (j=hEdgeOffsets[i];j<endOffset;++j)	
	hEdgeData.assign(ij++, pinList[j]);

      if(hypergraphModelType == 1 && addNonZero(i))
	hEdgeData.assign(ij++, i);

      hEdgeData[hEdgeLengthSlot] = ij-hEdgeLengthSlot;
      
      if(ij > maxIntsWritten || i == numCols-1)
	{
	  hEdgeData[0] = ij-1;
	  out_stream.write((char*)(hEdgeData.getArray()), sizeof(int)*ij);
	  ij = 0;
	}
    }

  out_stream.write((char*)(vWeights.getArray()), sizeof(int)*numRows);
  
  in_stream.close();
  out_stream.close();
}

  



#  endif