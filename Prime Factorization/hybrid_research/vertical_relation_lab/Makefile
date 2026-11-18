CXX ?= g++
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -Wpedantic
LDLIBS = -lgmpxx -lgmp
.PHONY: all clean test
all: vertical_lattice relation_bench
vertical_lattice: src/vertical_lattice.cpp
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDLIBS)
relation_bench: src/relation_bench.cpp
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDLIBS)
test: all
	python3 tests/test_vertical.py
	python3 tests/test_newton.py
	python3 tests/test_relations.py
	python3 tests/test_one_sided.py
	python3 tests/test_prime_power_region.py
clean:
	rm -f vertical_lattice relation_bench
