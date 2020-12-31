#ifndef __fftw_wrappers_h__
#define __fftw_wrappers_h__

// C++ wrappers for the FFTW3 library.
// Currently only one-dimensional transforms of real data to complex and back are supported.
//
// Typical use of FFTW3 entails the following steps:
// 1. Allocate input and output buffers.
// 2. Compute a "plan" struct. This tells FFTW which algorithms to use when actually computing the FFT.
// 3. Execute the FFT/IFFT operation on the buffers from step 1.

#include <complex>
#include <fftw3.h>
//#include "fftw3.h"

// Usage:
// 1. Fill input_buffer with input containing n_real_samples double numbers
//    (note, set_input_zeropadded will copy your buffer with optional zero padding)
// 2. Run execute().
// 3. Read output from output_buffer[0], ..., output_buffer[output_size-1].
//    Note that the output is composed of n_real_samples/2 + 1 complex numbers.
//
// These 3 steps can be repeated many times.

class FFTW_R2C_1D_Executor {
public:
    FFTW_R2C_1D_Executor(int n_real_samples);
    ~FFTW_R2C_1D_Executor();
    void set_input_zeropadded(const double* buffer, int size);
    void execute();

    const int input_size;
    double* input_buffer;

    const int output_size;
    std::complex<double>* output_buffer;

private:
    fftw_plan plan;
};

// Usage of this class is similar to that of FFTW_R2C_1D_Executor, only the input is n_real_samples/2+1 complex samples.
class FFTW_C2R_1D_Executor {
public:
    FFTW_C2R_1D_Executor(int n_real_samples);
    ~FFTW_C2R_1D_Executor();
//    void set_input_zeropadded(const std::complex<double>* buffer, int size);
    void execute();

    const int input_size;
    std::complex<double>* input_buffer;

    const int output_size;
    double* output_buffer;

private:
    fftw_plan plan;
};

#endif

