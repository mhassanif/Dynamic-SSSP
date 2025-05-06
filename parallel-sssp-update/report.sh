# Compile with profiling
mpic++ -pg -fopenmp -o sssp_update main/run_driver.cpp mpi/comm_utils.cpp core/sssp_update.cpp core/graph_structs.h

# Run with one rank to profile rank 0
mpirun -np 2./sssp_update insertions.txt deletions.txt

# Generate report
gprof ./sssp_update gmon.out > profile_report.txt

# View the report
cat profile_report.txt