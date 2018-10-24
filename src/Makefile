# makefile
# Purpose: Create the executable maxent-cpp for the code provided in main.cpp

# Includes
INC_EIGEN = -I ../include/eigen 
INC_ALL = $(INC_EIGEN)

# Objects
ALL_OBJECTS = main.o 

# Compiler options
CC = g++
CC_LINK = -ansi -Wall
CC_COMPILE = -c

# Executable
EXE = ../bin/maxent-cpp

# Heads
H_SHARED =  Nodes.hpp  Maxent.hpp
ALL_HEADS = $(H_SHARED) 

main.o: main.cpp $(ALL_HEADS)
	$(CC) $(CC_COMPILE) main.cpp $(INC_ALL)  

# Compile and build
maxent-cpp: $(ALL_OBJECTS)
	$(CC) -o $(EXE) $(ALL_OBJECTS)
  
# Remove objects and executable
clean:
	rm *.o $(EXE)


