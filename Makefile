FILE=Boundaries
CC=g++-5
CFLAGS=-O3 -Wall -Werror -pedantic -ansi -lm -I ~/.include/ -fopenmp
COMPILE_COMMAND=$(CC) $(CFLAGS)

all: setup ${FILE}.cpp
	$(COMPILE_COMMAND) -o $(FILE) $(FILE).cpp HPR.cpp

setup: 
	export OMP_NUM_THREADS=4
	export OMP_SCHEDULE=static

clean: 
	find . -name '*~' -delete
	find . -name '*.o' -delete
