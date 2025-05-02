#!/bin/bash

# Shell script to compile and run the partitioner program

# Function to display usage
usage() {
  echo "Usage: $0 <graph_file> <partition_file> <num_partitions>"
  echo "Example: $0 4elt.graph 4elt.graph.part.2 2"
  exit 1
}

# Check if the correct number of arguments is provided
if [ "$#" -ne 3 ]; then
  usage
fi

# Assign arguments to variables
GRAPH_FILE=$1
PARTITION_FILE=$2
NUM_PARTITIONS=$3

# Define a standard executable name
EXECUTABLE_NAME="partitioner"

# Compile the program
g++ -o "$EXECUTABLE_NAME" main.cpp partitioner.cpp -I/usr/include -L/usr/lib/x86_64-linux-gnu -lmetis
if [ $? -ne 0 ]; then
  echo "Compilation failed."
  exit 1
fi

# Run the program
./"$EXECUTABLE_NAME" "$GRAPH_FILE" "$PARTITION_FILE" "$NUM_PARTITIONS"
if [ $? -ne 0 ]; then
  echo "Execution failed."
  exit 1
fi
