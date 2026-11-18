CXX ?= g++
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -Wpedantic
ng_gmp: ng_gmp.cpp
	$(CXX) $(CXXFLAGS) $< -lgmpxx -lgmp -o $@
check: ng_gmp
	python3 test_ng.py
clean:
	rm -f ng_gmp
