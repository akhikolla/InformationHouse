#ifndef __fftwconvolver_h__
#define __fftwconvolver_h__

#include <vector>
#include <complex>

class FFTW_R2C_1D_Executor;
class FFTW_C2R_1D_Executor;

class FFTWConvolver {
public:
    FFTWConvolver(int maximum_input_size);
    ~FFTWConvolver();
    void convolve_same_size(int input_buffers_size, const double* input_buffer_0, const double* input_buffer_1, double* output_buffer);
private:
    std::vector<FFTW_R2C_1D_Executor*> r2c_executors;
    std::vector<FFTW_C2R_1D_Executor*> c2r_executors;

    std::complex<double>* tmp_complex;

    int maximum_input_size;
};

#endif
