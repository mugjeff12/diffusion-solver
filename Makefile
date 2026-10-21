CXX      = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -O2

test_thomas: test_thomas.o thomas.o
	$(CXX) $(CXXFLAGS) test_thomas.o thomas.o -o test_thomas

thomas.o: thomas.cpp thomas.hpp
	$(CXX) $(CXXFLAGS) -c thomas.cpp

test_thomas.o: test_thomas.cpp thomas.hpp
	$(CXX) $(CXXFLAGS) -c test_thomas.cpp

clean:
	rm -f *.o test_thomas

.PHONY: clean
