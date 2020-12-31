#include <iostream>
#include <cassert>
#include <cstring>
#include "fftw_wrappers.h"
#include "aligned_mem.h"
#include <Rcpp.h>

using namespace std;
using namespace Rcpp;

FFTW_R2C_1D_Executor::FFTW_R2C_1D_Executor(int n_real_samples) :
    input_size(n_real_samples),
    input_buffer(allocate_aligned_doubles(n_real_samples)),
    output_size(n_real_samples/2 + 1),
    output_buffer(allocate_aligned_complexes(n_real_samples/2 + 1))
{
    plan = fftw_plan_dft_r2c_1d(n_real_samples, input_buffer, reinterpret_cast<fftw_complex*>(output_buffer), FFTW_ESTIMATE);
}

FFTW_R2C_1D_Executor::~FFTW_R2C_1D_Executor()
{
    fftw_destroy_plan(plan);
    free_aligned_mem(input_buffer);
    free_aligned_mem(output_buffer);
}

void FFTW_R2C_1D_Executor::set_input_zeropadded(const double* buffer, int size)
{
    if (size > input_size) {
        //std::cerr << "size: " << size << "input_size: " << input_size << std::endl;
        Rcpp::Rcerr << "size: " << size << "input_size: " << input_size << std::endl;
    }
    assert(size <= input_size);
    memcpy(input_buffer, buffer, sizeof(double)*size);
    memset(&input_buffer[size], 0, sizeof(double)*(input_size - size));
}

void FFTW_R2C_1D_Executor::execute()
{
    fftw_execute(plan);
}

FFTW_C2R_1D_Executor::FFTW_C2R_1D_Executor(int n_real_samples) :
    input_size(n_real_samples/2 + 1),
    input_buffer(allocate_aligned_complexes(n_real_samples/2 + 1)),
    output_size(n_real_samples),
    output_buffer(allocate_aligned_doubles(n_real_samples))
{
    plan = fftw_plan_dft_c2r_1d(n_real_samples, reinterpret_cast<fftw_complex*>(input_buffer), output_buffer, FFTW_ESTIMATE);
}

FFTW_C2R_1D_Executor::~FFTW_C2R_1D_Executor()
{
    fftw_destroy_plan(plan);
    free_aligned_mem(input_buffer);
    free_aligned_mem(output_buffer);
}

//void FFTW_C2R_1D_Executor::set_input_zeropadded(const complex<double>* buffer, int size)
//{
//    assert(size == input_size);
//    memcpy(input_buffer, buffer, sizeof(complex<double>)*size);
//    memset(&input_buffer[size], 0, sizeof(complex<double>)*(input_size - size));
//}

void FFTW_C2R_1D_Executor::execute()
{
    fftw_execute(plan);
}
