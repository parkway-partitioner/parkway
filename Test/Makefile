
include ../Makefile.in

FLAGS = $(OPTS) $(DIRS) 
LIBDIR = -L ../
PARKWAYLIB = ../libparkway.a
TESTOBJ = driver.o reader.o

###
# Make test binary
###
test: $(TESTOBJ) $(PARKWAYLIB)
	 $(LINKER) $(FLAGS) $(TESTOBJ) $(LIBDIR) -o parkway $(LIBS) 

###
# Make Test Objects
###
.cpp.o:
	$(CC) $(FLAGS) -c $*.cpp

###
# Clean
###
clean:
	rm -f *.o parkway