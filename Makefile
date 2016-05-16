FILE=hello
CXX=g++
XXFLAGS=-O3 -Wall -Werror -pedantic -lm -fopenmp
OFLAGS=-ansi -I ~/.include/

all: setup ${FILE} 

setup: 
	export OMP_NUM_THREADS=4
	export OMP_SCHEDULE=static

%.o: %.cpp
	${CXX} -x c++ -c $< -o $@ ${OFLAGS}

${FILE}: ${FILE}.o Materials.o
	${CXX} $^ -o $@ ${XXFLAGS}

clean: 
	find . -name '*~' -delete
	find . -name '*.o' -delete
