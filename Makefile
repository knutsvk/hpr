FILE=test
CXX=g++
XXFLAGS=-O3 -Wall -Werror -pedantic -lm -fopenmp
OFLAGS=-ansi -I ~/.include/

all: setup ${FILE} 

setup:
	OMP_NUM_THREADS=4
	OMP_SCHEDULE=static

%.o: %.cpp
	${CXX} -x c++ -c $< -o $@ ${OFLAGS}

${FILE}: ${FILE}.o Materials.o
	${CXX} $^ -o $@ ${XXFLAGS}

clean: 
	find . -name '*~' -delete
	find . -name '*.o' -delete
