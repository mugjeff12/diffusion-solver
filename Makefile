CXX      = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -O2

all: diffusion test_thomas

diffusion: main.o thomas.o
	$(CXX) $(CXXFLAGS) main.o thomas.o -o diffusion

test_thomas: test_thomas.o thomas.o
	$(CXX) $(CXXFLAGS) test_thomas.o thomas.o -o test_thomas

main.o: main.cpp thomas.hpp
	$(CXX) $(CXXFLAGS) -c main.cpp

thomas.o: thomas.cpp thomas.hpp
	$(CXX) $(CXXFLAGS) -c thomas.cpp

test_thomas.o: test_thomas.cpp thomas.hpp
	$(CXX) $(CXXFLAGS) -c test_thomas.cpp

clean:
	rm -f *.o test_thomas diffusion

.PHONY: all clean
