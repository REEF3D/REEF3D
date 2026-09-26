# Standalone tests for the REEF3D rod-tree solver (no MPI/hypre needed)
CXX      ?= g++
CXXFLAGS ?= -std=c++20 -O2 -march=native
EIGEN    ?= ../../ThirdParty/eigen-5.0.0

test_rodtree: test_rodtree.cpp ../../src/rodtree.cpp ../../src/rodtree.h
	$(CXX) $(CXXFLAGS) -I $(EIGEN) -DEIGEN_MPL2_ONLY -o $@ test_rodtree.cpp ../../src/rodtree.cpp

run: test_rodtree
	./test_rodtree

clean:
	rm -f test_rodtree
