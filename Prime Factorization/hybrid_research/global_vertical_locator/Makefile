CXX ?= g++
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -Wpedantic
ENGINE ?= ../vf-factor-lab/cpp
all: global_vertical y_fermat_wheel
global_vertical: src/global_vertical.cpp
	$(CXX) $(CXXFLAGS) -I"$(ENGINE)" $< -lgmpxx -lgmp -o $@
y_fermat_wheel: src/y_fermat_wheel.cpp
	$(CXX) $(CXXFLAGS) $< -lgmpxx -lgmp -o $@
test: all
	python3 tests/test_global_vertical.py
	python3 tests/test_y_fermat_wheel.py
clean:
	rm -f global_vertical y_fermat_wheel
.PHONY: all test clean
