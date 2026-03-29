#pragma once

#ifdef _OPENMP
#include <omp.h>
#else
inline int omp_get_max_threads() { return 1; }
inline int omp_get_thread_num() { return 0; }
inline void omp_set_num_threads(int) {}
inline void omp_set_dynamic(int) {}
inline void omp_set_max_active_levels(int) {}
inline void omp_set_nested(int) {}
#endif
