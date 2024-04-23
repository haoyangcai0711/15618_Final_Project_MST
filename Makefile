# Define the compiler
CC=g++
# Define any compile-time flags
CFLAGS=-Wall -g -std=c++17 -fopenmp
# Define the target executables
TARGET1=mst
TARGET2=dynamic

# Default target
all: $(TARGET1) $(TARGET2)

# Rule to compile and link the first program
$(TARGET1): mst.cpp
	$(CC) $(CFLAGS) mst.cpp -o $(TARGET1)

# Rule to compile and link the second program
$(TARGET2): mst_dynamic.cpp
	$(CC) $(CFLAGS) mst_dynamic.cpp -o $(TARGET2)

# Rule for cleaning up
clean:
	rm -f $(TARGET1) $(TARGET2)
