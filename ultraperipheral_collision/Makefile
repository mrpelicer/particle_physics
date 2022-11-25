SOURCE= main.cpp
DEPS= sigma.cpp MV_nucleus.cpp MV_proton.cpp S_funcs.cpp FF_Charm.cpp

WALL = -std=c++14 -O3
GSLLIBS = -lgsl -lgslcblas -lm
LHAPDFLIBS = -lLHAPDF

.PHONY: all run

all: cmp run

cmp:
	g++ -g $(SOURCE) $(DEPS) $(WALL) $(GSLLIBS) $(LHAPDFLIBS) -o run
run:
	./run
