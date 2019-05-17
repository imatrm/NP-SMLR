CXX = g++
CXXFLAGS = -O3

all :
	$(CXX) $(CXXFLAGS) ./src/Likelihood.cpp -o ./bin/Likelihood
	$(CXX) $(CXXFLAGS) ./src/Detection.cpp -o ./bin/Detection
	$(CXX) $(CXXFLAGS) ./src/NclsPos.cpp -o ./bin/NclsPos

