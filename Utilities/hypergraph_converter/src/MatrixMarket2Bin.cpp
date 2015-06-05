
#  ifndef _MTXMKT_TOBIN_CPP
#  define _MTXMKT_TOBIN_CPP


// ### MatrixMarket2Bin.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
// 
// HISTORY: 
// 
// 12/1/2004: Last Modified
//
// ###


#  include "MatrixMarket2Bin.hpp"


MatrixMarket2Bin::MatrixMarket2Bin()
  : MatrixMarketReader()
{

}



MatrixMarket2Bin::~MatrixMarket2Bin()
{


}


void MatrixMarket2Bin::convert(const char *filename)
{
  char bin_file[512];
  sprintf(bin_file, "%s.bin", filename);

  readMatrix(filename);
  writeMatrix(bin_file);
}


void MatrixMarket2Bin::writeMatrix(const char *filename)
{
  ofstream out_stream;

  int binPreamble[3];
  int dataLength;
  int startOffset;
  int endOffset;
  int hEdgeChunkEntry;
  int i;
  int j;
  
  FastDynaArray<int> writeData;

  out_stream.open(filename, ofstream::out | ofstream::app | ofstream::binary);
  
  if(!out_stream.is_open())
    {
      cout << "error opening " << filename << endl;  
      exit(1);
    }

  binPreamble[0] = numVertices;
  binPreamble[1] = numHyperedges;
  binPreamble[2] = numPins;

  out_stream.write((char*)(&binPreamble[0]), sizeof(int)*3);
  
  dataLength = 0;
  writeData.assign(dataLength++, -1);
  
  for (i=0;i<numHyperedges;++i)
    {
      hEdgeChunkEntry = dataLength;
      writeData.assign(dataLength++, -1);
      writeData.assign(dataLength++, 1);
      
      startOffset = hEdgeOffsets[i];
      endOffset = hEdgeOffsets[i+1];
      
      for (j=startOffset;j<endOffset;++j)		
	writeData.assign(dataLength++, pinList[j]);	
      
      writeData[hEdgeChunkEntry] = (endOffset-startOffset)+2;
    }
  
  writeData[0] = dataLength-1;
  out_stream.write((char*)(writeData.getArray()), sizeof(int)*dataLength);	 

  writeData.setLength(numVertices);
  for (i=0;i<numVertices;++i)
    writeData[i] = 1;
  
  out_stream.write((char*)(writeData.getArray()), sizeof(int)*numVertices);	 
  out_stream.close();
}



#  endif
