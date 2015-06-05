# Convenience makefile -- passes through to the Makefile in the build
# directory. Runs the configure script otherwise.
SHELL := /bin/bash

all: ./build/Makefile
	@ $(MAKE) -C build

./build/Makefile:
	@ echo "-- Running configure script..."
	@ ./configure

$(MAKECMDGOALS): ./build/Makefile
	@ $(MAKE) -C build $(MAKECMDGOALS)
