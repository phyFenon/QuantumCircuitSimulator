CC = g++-13
CXXFLAGS = -std=c++17 -fopenmp
TARGET = out
SRCS = main.cpp quantum_circuit/quantum_simulator.cpp quantum_circuit/singlegate_mat.cpp
OBJS = $(SRCS:.cpp=.o)

$(TARGET): $(OBJS)
	$(CC) $(CXXFLAGS) $(OBJS) -o $(TARGET)

%.o: %.cpp
	$(CC) $(CXXFLAGS) -c $< -o $@

.PHONY: clean
clean:
	rm -f $(OBJS) $(TARGET)
