INFILE=main
OUTFILE=Simulate
CC=g++
CFLAGS=-O3 -Wall -Werror -pedantic -ansi -I ~/.include/ -fopenmp -Wno-long-long 
LDFLAGS=-lm -lconfig++

all: ${INFILE}.cpp
	$(CC) $(CFLAGS) -o $(OUTFILE) $(INFILE).cpp HPR.cpp $(LDFLAGS)

clean: 
	find . -name '*~' -delete
	find . -name '*.o' -delete
