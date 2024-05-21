#include <cmath>
#include <cstdlib>
#include <complex>
#include <limits.h>

std::complex<double>* dft(double* wave, int num_of_inputs);

double* hamming(int num_of_inputs);

double* hamming_data(double* ft_coe, int num_of_inputs);

std::complex<double>* fft(std::complex<double> *wave, int size);
