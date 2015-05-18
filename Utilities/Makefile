
include ../Makefile.in

FLAGS = $(OPTS) $(DIRS) 
CONOBJ = StringUtils.o TextFileReader.o HypergraphTextFileReader.o Converter.o FromBinConverter.o TextToBinConverter.o HMeTiS2Bin.o PaToH2Bin.o Sandia2Bin.o Bin2Para.o Bin2HMeTiS.o Bin2Sandia.o Bin2PaToH.o HarwellBoeingReader.o HarwellBoeingToBin.o HarwellBoeingTo2dBin.o HarwellBoeingToTxt.o MatrixMarketReader.o MatrixMarket2Bin.o Convert.o 
DISPOBJ = PrintHypergraphChars.o StringUtils.o
OTHPROGS = printgraph

###
# Make convert
###
convert: $(CONOBJ) $(OTHPROGS)
	$(CC) $(FLAGS) $(CONOBJ) -o convert 

printgraph: $(DISPOBJ)
	$(CC) $(FLAGS) $(DISPOBJ) -o printgraph 

###
# Make objects
###
.cpp.o:
	$(CC) $(FLAGS) -c $*.cpp

###
# Clean
###
clean:
	rm -f *.o convert printgraph