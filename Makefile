FILE=test
CXX=g++
XXFLAGS=-O3 -Wall -Werror -pedantic -lm 
OFLAGS=-ansi -I ~/.include/

all: ${FILE} 

%.o: %.cpp
	${CXX} -x c++ -c $< -o $@ ${OFLAGS}

${FILE}: ${FILE}.o Materials.o
	${CXX} $^ -o $@ ${XXFLAGS}

clean: 
	find . -name '*~' -delete
	find . -name '*.o' -delete
