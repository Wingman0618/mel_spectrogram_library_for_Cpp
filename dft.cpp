/*
Reference: https://www.youtube.com/watch?v=nl9TZanwbBk
           https://www.youtube.com/watch?v=Xw4voABxU5c
*/
#include "dft.h"

std::complex<float>* dft(float* wave, int num_of_inputs){
        int n = num_of_inputs;

	std:: complex<float>* result = (std::complex<float>*)malloc(n*sizeof(std::complex<float>));

        // Define the omega w
        // C++ complex implementation tutrial:
        //  https://www.geeksforgeeks.org/complex-numbers-c-set-1/
        std::complex<float> i(0., ((-2.*M_PI)/n));
        std::complex<float> w = exp(i);

        for(int k=0; k<n; k++){
                std::complex<float> coe = (0, 0);
                for(int j=0; j<n; j++){
                         coe += wave[j]*pow(w,j*k);
                }
                *(result+k) = coe;
        }

        return result;
}

float* hamming(int num_of_inputs){
	float *weight = (float*)malloc(num_of_inputs*sizeof(float));

	// Define min and max values for normalising the data
	float max;
	float min;
	for(int  i=0; i<num_of_inputs; i++){
		*(weight+i) = 0.54 - 0.46*cos((2.*M_PI*float(i))/(float(num_of_inputs)-1.));
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
