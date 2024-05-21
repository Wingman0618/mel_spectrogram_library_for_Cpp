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

/*
        // Define the omega w
        // C++ complex implementation tutrial:
        //  https://www.geeksforgeeks.org/complex-numbers-c-set-1/
        std::complex<float> c(0., ((-2.*M_PI)/n));
        std::complex<float> w = exp(c);

        for(int k=0; k<n; k++){
                std::complex<float> coe(0, 0);
                for(int j=0; j<n; j++){
                         coe += wave_to_complex[j]*pow(w,j*k);
                }
                *(result+k) = coe;
        }
*/
	result = fft(wave_to_complex, n);

        return result;
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
