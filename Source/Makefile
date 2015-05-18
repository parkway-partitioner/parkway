
include ../Makefile.in

###
# Compiler flags
###
CFLAGS = $(OPTS) -fPIC $(DIRS)

###
# General Objects
###
BASIC_OBJ = MemoryTracker.o Funct.o CompleteBinaryTree.o Log.o HashTables.o MovementSets.o Hypergraph.o HypergraphLoader.o GlobalCommunicator.o ParaHypergraph.o ParaHypergraphLoader.o

###
# Coarsening Objects
###
COARSE_OBJ = Coarsener.o FCCoarsener.o RestrCoarsener.o RestrFCCoarsener.o ParaCoarsener.o ParaRestrCoarsener.o ParaFCCoarsener.o Para2DModelCoarsener.o ParaRestrFCCoarsener.o

###
# Refiner Objects
###
REFINER_OBJ = Refiner.o FMRefiner.o GreedyKwayRefiner.o InitBisector.o ParaRefiner.o ParaGreedyKwayRefiner.o

###
# Controller Objects
###
CONTROL_OBJ = BisectionController.o VCycleBisectionController.o VCycleFinalBisectionController.o VCycleAllBisectionController.o SeqController.o RecurBisectController.o KHMetisController.o PaToHController.o ParaController.o BasicParaController.o ParaVCycleController.o ParaVCycleFinalController.o ParaVCycleAllController.o

###
# Auxiliary Objects
###
UTIL = Utils.o Parkway.o

###
# All Objects
###
OBJECTS = $(BASIC_OBJ) $(COARSE_OBJ) $(REFINER_OBJ) $(CONTROL_OBJ) $(UTIL) 

###
# Make Objects
###
.cpp.o:
	$(CC) $(CFLAGS) -c $*.cpp

###
# Make lib
###
../libparkway.a: $(OBJECTS)
	ar rcs ../libparkway.a $(OBJECTS)

###
# Delete Objects
###
clean:
	rm -f *.o 

