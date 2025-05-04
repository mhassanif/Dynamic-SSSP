#!/bin/bash

# Shell script to compile and run the make_subgraphs program

# Function to display usage instructions
usage() {
  echo "Usage: $0 <num_partitions>"
  echo "Example: $0 3"
  exit 1
}

# Check if the correct number of arguments is provided
if [ "$#" -ne 1 ]; then
  usage
fi

# Assign the number of partitions to a variable
NUM_PARTITIONS=$1

# Define the executable name
EXECUTABLE_NAME="make_subgraphs"

# Compile the make_subgraphs.cpp program
g++ -o "$EXECUTABLE_NAME" make_subgraphs.cpp -I/usr/include -L/usr/lib/x86_64-linux-gnu
if [ $? -ne 0 ]; then
  echo "Compilation failed."
  exit 1
fi

# Run the compiled program with the specified number of partitions
./"$EXECUTABLE_NAME" "$NUM_PARTITIONS"
if [ $? -ne 0 ]; then
  echo "Execution failed."
  exit 1
fi
