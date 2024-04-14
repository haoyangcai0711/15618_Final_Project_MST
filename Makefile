# Define the compiler
CC=g++

# Define any compile-time flags
CFLAGS=-Wall -g -std=c++17 -fopenmp

# Define the target executable
TARGET=mst

# Rule to link the program
$(TARGET): mst.o
	$(CC) $(CFLAGS) mst.o -o $(TARGET)

# Rule to compile the source files
mst.o: mst.cpp
	$(CC) $(CFLAGS) -c mst.cpp

# Rule for cleaning up
clean:
	rm -f $(TARGET) *.o
