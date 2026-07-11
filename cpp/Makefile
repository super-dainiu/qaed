CXX      ?= g++
CXXFLAGS ?= -O3 -march=native -funroll-loops -std=c++17 -fopenmp -Wall -Wextra
# FlexiBLAS provides BLAS + LAPACK on RHEL9 (only the .so.3 is installed).
LDLIBS   ?= -l:libflexiblas.so.3

qaed_bench: bench.cpp quat.hpp qmat.hpp qaed.hpp
	$(CXX) $(CXXFLAGS) bench.cpp -o $@ $(LDLIBS)

test: qaed_bench
	./qaed_bench --selftest

clean:
	rm -f qaed_bench

.PHONY: test clean
