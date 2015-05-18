###
# Make Parkway
###
parkway:
	(cd Source ; make)
	(cd Test ; make)
	(cd Utilities ; make)

###
# Clean
###
clean:
	(rm -f libparkway.a)
	(cd Source ; make clean)
	(cd Test ; make clean)
	(cd Utilities ; make clean)
