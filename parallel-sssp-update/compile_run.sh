mpic++ -fopenmp -o sssp_update main/run_driver.cpp mpi/comm_utils.cpp core/sssp_update.cpp core/graph_structs.h
mpirun -np 2 ./sssp_update insertions.txt deletions.txt