INFILE=main
OUTFILE=Simulate
CC=g++-4.9
CFLAGS=-O3 -Wall -Werror -pedantic -ansi -lm -lconfig++ -I ~/.include/ -fopenmp -Wno-long-long 
COMPILE_COMMAND=$(CC) $(CFLAGS)

all: ${INFILE}.cpp
	$(COMPILE_COMMAND) -o $(OUTFILE) $(INFILE).cpp HPR.cpp

clean: 
	find . -name '*~' -delete
	find . -name '*.o' -delete
