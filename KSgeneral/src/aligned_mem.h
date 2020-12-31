#ifndef __aligned_mem__
#define __aligned_mem__

#include <complex>
#include <cstdlib>
#include "mm_malloc.h"

#define ALIGNMENT (32)

inline double* allocate_aligned_doubles(int n)
{
    return static_cast<double*>(_mm_malloc(n*sizeof(double), ALIGNMENT));
}

inline std::complex<double>* allocate_aligned_complexes(int n)
{
    return static_cast<std::complex<double>*>(_mm_malloc(n*sizeof(std::complex<double>), ALIGNMENT));
}


inline void free_aligned_mem(void* p)
{
    _mm_free(p);
}

#endif
