/*
Reference: https://www.youtube.com/watch?v=nl9TZanwbBk
           https://www.youtube.com/watch?v=Xw4voABxU5c
*/
#include "dft.h"
#include <iostream>

std::complex<double>* dft(double* wave, int num_of_inputs){
        int n = num_of_inputs;

	std::complex<double>* result = (std::complex<double>*)malloc(n*sizeof(std::complex<double>));
        std::complex<double>* wave_to_complex = (std::complex<double>*)malloc(n*sizeof(std::complex<double>));
        for(int i=0; i<n; i++){
                std::complex<double> c_wave(*(wave+i), 0);
                *(wave_to_complex+i) = c_wave;
	}

	result = fft(wave_to_complex, n);

        return result;
}

std::complex<double>* dft_rfft(double* wave, int num_of_inputs){
        int n = num_of_inputs;

        std::complex<double>* result_fft = (std::complex<double>*)malloc(n*sizeof(std::complex<double>));
	std::complex<double>* result_rfft = (std::complex<double>*)malloc((n/2)*sizeof(std::complex<double>));
        std::complex<double>* wave_to_complex = (std::complex<double>*)malloc(n*sizeof(std::complex<double>));
        for(int i=0; i<n; i++){
                std::complex<double> c_wave(*(wave+i), 0);
                *(wave_to_complex+i) = c_wave;
        }

        result_fft = fft(wave_to_complex, n);

	for(int i=0; i<n/2; i++){
		*(result_rfft+i) = *(result_fft+i);
	}

        return result_rfft;
}


double* hamming(int num_of_inputs){
	double *weight = (double*)malloc(num_of_inputs*sizeof(double));

	// Define min and max values for normalising the data
	double max;
	double min;
	for(int  i=0; i<num_of_inputs; i++){
		*(weight+i) = 0.54 - 0.46*cos((2.*M_PI*double(i))/(double(num_of_inputs)-1.));
		if(i==0){
			max = *(weight+i);
			min = *(weight+i);
		}
		if(i>0){
			if(*(weight+i)>=max){
				max = *(weight+i);
			}
			if(*(weight+i)<=min){
				min = *(weight+i);
			}
		}
	}
	// Normalise the data
	for(int i=0; i<num_of_inputs; i++){
		*(weight+i) = (*(weight+i)-min)/(max-min);
	}

	return weight;
}

double* hamming_data(double* wave, int num_of_inputs){
        double *weight = (double*)malloc(num_of_inputs*sizeof(double));
	double *result = (double*)malloc(num_of_inputs*sizeof(double));

        // Define min and max values for normalising the data
        double max;
        double min;

        for(int  i=0; i<num_of_inputs; i++){
                *(weight+i) = 0.54 - 0.46*cos((2.*M_PI*double(i))/(double(num_of_inputs)-1.));
                if(i==0){
                        max = *(weight+i);
                        min = *(weight+i);
                }
                if(i>0){
                        if(*(weight+i)>=max){
                                max = *(weight+i);
                        }
                        if(*(weight+i)<=min){
                                min = *(weight+i);
                        }
                }

        }
        // Normalise the data

        for(int i=0; i<num_of_inputs; i++){
                *(weight+i) = (*(weight+i)-min)/(max-min);
        }

	// weight the data
	for(int i=0; i<num_of_inputs; i++){
		*(result+i) = *(wave+i)*weight[i];
	}

        return result;
}

double** stft(double *wave, int num_of_inputs, int window_size, int step){
	int num_of_ffts = (num_of_inputs - window_size)/step;
//	std::cout<<num_of_ffts<<std::endl;

	double** fft_array = (double**)malloc(num_of_ffts*sizeof(double*));
	for(int i=0; i<num_of_ffts; i++){
		*(fft_array+i) = (double*)malloc(window_size*sizeof(double));
	}

	double* window = (double*)malloc(window_size*sizeof(double));
	std::complex<double>* window_dft = (std::complex<double>*)malloc((window_size/2)*sizeof(std::complex<double>));
	int pointer = 0;

	for(int i=0; i<num_of_ffts; i++){
		for(int j=0; j<window_size; j++){
			*(window+j) = *(wave+pointer+j);
		}
		window = hamming_data(window, window_size);
		pointer += step;
		window_dft = dft_rfft(window, window_size);
		for(int k=0; k<window_size/2; k++){
			*(*(fft_array+i)+k) = abs(*(window_dft+k));
		}
	}

	return fft_array;
}

std::complex<double>* fft(std::complex<double> *wave, int size){
	int n = size;
	if(log2(n) != (float)(int)log2(n)){
		std::cout<<"input size is not 2^n"<<std::endl;
		return wave;
	}
	if(n==1){
		return wave;
	}

	int m = n/2;
	std::complex<double>* feven = (std::complex<double>*)malloc(m*sizeof(std::complex<double>));
	std::complex<double>* fodd = (std::complex<double>*)malloc(m*sizeof(std::complex<double>));

	for(int i=0; i<m; i++){
		*(feven+i) = *(wave+2*i);
		*(fodd+i) = *(wave+2*i+1);
	}

        std::complex<double>* Feven = (std::complex<double>*)malloc(m*sizeof(std::complex<double>));
	Feven = fft(feven, m);
        std::complex<double>* Fodd = (std::complex<double>*)malloc(m*sizeof(std::complex<double>));
	Fodd = fft(fodd, m);

	std::complex<double>* freqbins = (std::complex<double>*)malloc(n*sizeof(std::complex<double>));
	for(int k=0; k<m; k++){
	        std::complex<double> c_e = std::polar(1.0, -2*M_PI*k/n)*Fodd[k];

		freqbins[k] = Feven[k] + c_e;
		freqbins[k+m] = Feven[k] - c_e;
	}

	return freqbins;
}

