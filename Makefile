FILE=test
CC=g++-5
CFLAGS=-O3 -Wall -Werror -pedantic -ansi -lm -fopenmp -I ~/.include/
COMPILE_COMMAND=$(CC) $(CFLAGS)
OUTPUT=test

all: setup ${FILE}.cpp
	$(COMPILE_COMMAND) -o $(OUTPUT) $(FILE).cpp Materials.cpp

setup: 
	export OMP_NUM_THREADS=4
	export OMP_SCHEDULE=static

clean: 
	find . -name '*~' -delete
	find . -name '*.o' -delete
