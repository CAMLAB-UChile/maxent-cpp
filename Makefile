# makefile
# Purpose: Create the EXECUTABLE for maxent-cpp

SRCDIR = src
EXEC = bin/maxent-cpp
MAKE = make
RM   = rm -rf

# To create the executable file in the bin directory, execute "make" or 
# "make maxent-cpp" on the command line

maxent:	
	    cd $(SRCDIR); $(MAKE) maxent-cpp

clean :	
	    cd $(SRCDIR); $(MAKE) clean
