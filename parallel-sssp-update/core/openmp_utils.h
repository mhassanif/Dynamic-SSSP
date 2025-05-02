#ifndef OPENMP_UTILS_H
#define OPENMP_UTILS_H

#ifdef _OPENMP
#include <omp.h>
#endif

// Set OpenMP thread count
inline void set_omp_threads(int n) {
#ifdef _OPENMP
    omp_set_num_threads(n);
#endif
}

#endif // OPENMP_UTILS_H