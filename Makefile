FILE=CylindricalShock
CC=g++-5
CFLAGS=-O3 -Wall -Werror -pedantic -ansi -lm -I ~/.include/ -fopenmp
COMPILE_COMMAND=$(CC) $(CFLAGS)

all: ${FILE}.cpp
	$(COMPILE_COMMAND) -o $(FILE) $(FILE).cpp HPR.cpp

clean: 
	find . -name '*~' -delete
	find . -name '*.o' -delete
