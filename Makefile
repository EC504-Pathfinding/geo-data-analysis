# Makefile for compiling the C++ program

# Compiler
CXX = g++

# Compiler flags
CXXFLAGS = -Wall -Wextra -std=c++11

# Target executable name
TARGET = GeoDataAnalysis

# Source files
SOURCES = construct_graph.cpp

# Build the program
all: $(TARGET)

$(TARGET): $(SOURCES)
	$(CXX) $(CXXFLAGS) $(SOURCES) -o $(TARGET)

# Clean up
clean:
	rm -f $(TARGET)
