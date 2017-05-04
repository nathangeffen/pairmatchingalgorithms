# the compiler: gcc for C program, define as g++ for C++
CXX = g++

# compiler flags:
#  -g    adds debugging information to the executable file
#  -Wall turns on most, but not all, compiler warnings
CXXFLAGS  = -g -Wall -std=c++11
RELFLAGS = -O3 -std=c++11
EXTRAFLAGS =

# the build target executable:
TARGET = partnermatch

all: $(TARGET)

$(TARGET): $(TARGET).cc
	$(CXX) $(CXXFLAGS) $(EXTRAFLAGS) -o $(TARGET) $(TARGET).cc

dev_attract: clean
	$(CXX) $(CXXFLAGS) $(EXTRAFLAGS) -D ATTRACT_REJECT -o $(TARGET) $(TARGET).cc

release: clean
	$(CXX) $(RELFLAGS) $(EXTRAFLAGS) -o $(TARGET) $(TARGET).cc

release_attract: clean
	$(CXX) $(RELFLAGS) -D ATTRACT_REJECT -o $(TARGET) $(TARGET).cc

clean:
	$(RM) $(TARGET)
