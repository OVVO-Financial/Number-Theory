CXX ?= g++
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -Wpedantic
LDLIBS = -lgmpxx -lgmp
all: dyadic_scan
dyadic_scan: dyadic_scan.cpp
	$(CXX) $(CXXFLAGS) $< $(LDLIBS) -o $@
test: dyadic_scan
	python3 test_dyadic_vertical.py
	./dyadic_scan 8051 3 2 2
	./dyadic_scan 978508015703 6 50746 256
	./dyadic_scan 978508015703 8 12686 256
	./dyadic_scan 187 4 5 4
clean:
	rm -f dyadic_scan
.PHONY: all test clean
