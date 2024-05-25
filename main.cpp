/*
Reference: https://www.youtube.com/watch?v=nl9TZanwbBk
 	   https://www.youtube.com/watch?v=Xw4voABxU5c
*/
#include <iostream>
#include "dft.h"
#include "test.h"

using namespace std;

int main(){

	int size = 16;

	float x_input [size] = {0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15.};

	double* y_output = (double*)malloc(size*sizeof(double));

        for(int i=0; i<size; i++){
                y_output[i] = sin(x_input[i]);
        }
	cout<<"basic inputs and outputs---------------------"<<endl;
	cout<<"inputs: ";
	for(int i=0; i<size; i++){
		cout<<x_input[i]<<" ";
	}

	cout<<endl;
	cout<<"outputs: ";
	for(int i=0; i<size; i++){
		cout<<y_output[i]<<" ";
	}

	cout<<endl;

	complex<double>* dft_complex = (complex<double>*)malloc(size*sizeof(complex<double>));
	dft_complex = dft(y_output, size);

	cout<<"DFT without hamming---------------------"<<endl;
	for(int i=0; i<size; i++){
		cout<<abs(dft_complex[i])<<" ";
	}
	cout<<endl;

	// add hamming
	cout<<"Hamming---------------------"<<endl;
	double *weight = (double*)malloc(size*sizeof(double));
	weight = hamming(size);
	for(int i=0; i<size; i++){
		y_output[i] = y_output[i]*weight[i];
	}

        dft_complex = dft(y_output, size);
        for(int i=0; i<size; i++){
                cout<<abs(dft_complex[i])<<" ";
        }
        cout<<endl;

	// test audio input
	cout<<"audio input test---------------------"<<endl;
	int n = 512;
	cout<<"size: "<<n<<endl;

	double *data = (double*)realloc(y_output, n*sizeof(double));

	for(int i=0; i<n; i++){
		*(data+i) = double(test_input[i]);
	}

	complex<double>* dft_audio = (complex<double>*)realloc(dft_complex, n*sizeof(complex<double>));
	data = hamming_data(data, n);

	dft_audio = dft(data, n);
	for(int i=0; i<5; i++){
		cout<<abs(dft_audio[i])<<" ";
	}
	cout<<endl;

	// test stdf
	cout<<"audio input test (stft)---------------------"<<endl;
	n = 16000-128;
	cout<<"size: "<<n<<endl;

	double* new_data = (double*)realloc(data, n*sizeof(double));
	for(int i=0; i<n; i++){
		*(new_data+i) = double(test_input[i]);
        }

	for(int i=0; i<8; i++){
		cout<<new_data[i]<<" ";
	}
	cout<<endl;

	int window_size = 512;
	int step = 256;
	int num_of_fft = 60;
	// double** stft_arr = (double**)malloc(num_of_fft*sizeof(double*));

	double** stft_arr =  stft(new_data, n, window_size, step);

	// print first 5 elements of each fft
	for(int i=0; i<num_of_fft; i++){
        	for(int j=0; j<5; j++){
        	        cout<<*(*(stft_arr+i)+j)<<" ";
	        }
		cout<<endl;
		cout<<"-------------------------"<<endl;
	}
	// Check the size
	cout<<stft_arr[59][256]<<endl;
	cout<<stft_arr[59][257]<<endl;
	cout<<"-------------------------"<<endl;

	// Check mel spectrogram
	cout<<"-------------------------"<<endl;
	cout<<"mel spectrogram check: "<<endl;

	double** mel_spectrogram = mel_spec(new_data, n, window_size, step, 40, 16000);
	for(int i=0; i<num_of_fft; i++){
		for(int j=0; j<5; j++){
                	cout<<*(*(mel_spectrogram+i)+j)<<" ";
		}
		cout<<endl;
		cout<<"-------------------------"<<endl;
        }

        // Check the size
        cout<<mel_spectrogram[59][39]<<endl;
        cout<<mel_spectrogram[59][40]<<endl;
	cout<<"-------------------------"<<endl;

	// Check power to db
	cout<<"-------------------------"<<endl;
	cout<<"Power to db: "<<endl;

	mel_spectrogram = power_to_db(mel_spectrogram, 60, 40);
        for(int i=0; i<num_of_fft; i++){
                for(int j=0; j<5; j++){
                        cout<<*(*(mel_spectrogram+i)+j)<<" ";
                }
                cout<<endl;
                cout<<"-------------------------"<<endl;
        }

        // Check amplitude to db
        cout<<"-------------------------"<<endl;
        cout<<"Amplitude to db: "<<endl;
	double** mel_spectrogram_new = mel_spec(new_data, n, window_size, step, 40, 16000);

        mel_spectrogram_new = amplitude_to_db(mel_spectrogram_new, 60, 40, 80);
        for(int i=0; i<num_of_fft; i++){
                for(int j=0; j<5; j++){
                        cout<<*(*(mel_spectrogram_new+i)+j)<<" ";
                }
                cout<<endl;
                cout<<"-------------------------"<<endl;
	}

	free(new_data);
	free(dft_audio);
	free(weight);
	free(stft_arr);
	free(mel_spectrogram);

	return 0;
}
