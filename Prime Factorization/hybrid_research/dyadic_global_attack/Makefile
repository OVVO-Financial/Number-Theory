CXX ?= g++
CXXFLAGS ?= -std=c++17 -O3 -Wall -Wextra -Wpedantic
LDLIBS = -lgmpxx -lgmp
TOOLS=gaussian_norm_one smooth_order newton2_primary dual_residual hybrid_native
all: $(TOOLS)

smooth_order: smooth_order.cpp
	$(CXX) $(CXXFLAGS) $< $(LDLIBS) -o $@
gaussian_norm_one: gaussian_norm_one.cpp
	$(CXX) $(CXXFLAGS) $< $(LDLIBS) -o $@
newton2_primary: newton2_primary.cpp
	$(CXX) $(CXXFLAGS) $< $(LDLIBS) -o $@

dual_residual: dual_residual.cpp
	$(CXX) $(CXXFLAGS) $< $(LDLIBS) -o $@

smooth_order.o: smooth_order.cpp
	$(CXX) $(CXXFLAGS) -DDYADIC_LIB -c $< -o $@
gaussian_norm_one.o: gaussian_norm_one.cpp
	$(CXX) $(CXXFLAGS) -DDYADIC_LIB -c $< -o $@
newton2_primary.o: newton2_primary.cpp
	$(CXX) $(CXXFLAGS) -DDYADIC_LIB -c $< -o $@

hybrid_native: hybrid_native.cpp smooth_order.o gaussian_norm_one.o newton2_primary.o
	$(CXX) $(CXXFLAGS) $^ $(LDLIBS) -o $@

test: all
	python3 lane_cover.py --output lane_census.json
	python3 test_newton_conjugacy.py
	python3 test_gaussian.py
	python3 test_vertical_group_bridge.py
	python3 test_lane_selector.py
	python3 test_norm_channel.py
	python3 test_native_commands.py

clean:
	rm -f $(TOOLS) *.o
.PHONY: all clean test
