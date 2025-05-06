#!/bin/bash

# Shell script to compile and run the partitioner program

# Function to display usage
usage() {
  echo "Usage: $0 <num_partitions>"
  echo "Example: $0 2"
  exit 1
}

# Check if the correct number of arguments is provided
if [ "$#" -ne 1 ]; then
  usage
fi

# Assign the number of partitions to a variable
NUM_PARTITIONS=$1

# Define the executable name
EXECUTABLE_NAME="partitioner"

# Compile the program
g++ -o "$EXECUTABLE_NAME" main.cpp partitioner.cpp -I/usr/include -L/usr/lib/x86_64-linux-gnu -lmetis
if [ $? -ne 0 ]; then
  echo "Compilation failed."
  exit 1
fi

# Run the program with the specified number of partitions
./"$EXECUTABLE_NAME" "$NUM_PARTITIONS"
if [ $? -ne 0 ]; then
  echo "Execution failed."
  exit 1
fi
