#include <cmath>
#include <cstdlib>
#include <complex>
#include <limits.h>

std::complex<double>* dft(double* wave, int num_of_inputs);

std::complex<double>* dft_rfft(double* wave, int num_of_inputs);

double* hamming(int num_of_inputs);

double* hamming_data(double* ft_coe, int num_of_inputs);

double** stft(double *wave, int num_of_inputs, int window_size, int hop_length);

double** mel_spec(double *wave, int num_of_inputs, int window_size, int hop_length, int n_mels, int sample_rate);

double** filterbank_gen(double *filterbin, int n_mels, int n_ffts);

double** power_to_db(double** mel_spec, int num_of_ffts, int n_mels);

std::complex<double>* fft(std::complex<double> *wave, int size);
