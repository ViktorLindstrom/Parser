# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++20 -Wall -O2 -march=native

SRCS = parser.cpp
OBJS = $(SRCS:.cpp=.o)
TARGET = parser_executable

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)

.PHONY: all clean
