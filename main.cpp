/*
Reference: https://www.youtube.com/watch?v=nl9TZanwbBk
 	   https://www.youtube.com/watch?v=Xw4voABxU5c
*/

#include <iostream>
#include "dft.h"

using namespace std;

int main(){

	int size = 16;

	float x_input [size] = {0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15.};

	float* y_output = (float*)malloc(size*sizeof(float));

        for(int i=0; i<size; i++){
                y_output[i] = sin(x_input[i]);
        }
	cout<<"basic inputs and outputs: "<<endl;
	for(int i=0; i<size; i++){
		cout<<x_input[i]<<" ";
	}

	cout<<endl;

	for(int i=0; i<size; i++){
		cout<<y_output[i]<<" ";
	}

	cout<<endl;

	complex<float>* dft_complex = (complex<float>*)malloc(size*sizeof(complex<float>));
	dft_complex = dft(y_output, size);

	cout<<"DFT without hamming: "<<endl;
	for(int i=0; i<size; i++){
		cout<<abs(dft_complex[i])<<" ";
	}
	cout<<endl;
	cout<<"Hamming: "<<endl;
	// add hamming
	float *weight = (float*)malloc(size*sizeof(float));
	weight = hamming(size);
	for(int i=0; i<size; i++){
		y_output[i] = y_output[i]*weight[i];
	}

        dft_complex = dft(y_output, size);
        for(int i=0; i<size; i++){
                cout<<abs(dft_complex[i])<<" ";
        }
        cout<<endl;

	free(y_output);
	free(dft_complex);
	free(weight);

	return 0;
}
