/*
Reference: https://www.youtube.com/watch?v=nl9TZanwbBk
           https://www.youtube.com/watch?v=Xw4voABxU5c
	   http://practicalcryptography.com/miscellaneous/machine-learning/guide-mel-frequency-cepstral-coefficients-mfccs/
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
	int n_fft = (n/2)+1;
        std::complex<double>* result_fft = (std::complex<double>*)malloc(n*sizeof(std::complex<double>));
	std::complex<double>* result_rfft = (std::complex<double>*)malloc(n_fft*sizeof(std::complex<double>));
        std::complex<double>* wave_to_complex = (std::complex<double>*)malloc(n*sizeof(std::complex<double>));
        for(int i=0; i<n; i++){
                std::complex<double> c_wave(*(wave+i), 0);
                *(wave_to_complex+i) = c_wave;
        }

        result_fft = fft(wave_to_complex, n);

	for(int i=0; i<n_fft; i++){
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

double** stft(double *wave, int num_of_inputs, int window_size, int hop_length){
	int num_of_ffts = (num_of_inputs - window_size)/hop_length;
	int n_fft = (window_size/2)+1;
//	std::cout<<num_of_ffts<<std::endl;

	double** fft_array = (double**)malloc(num_of_ffts*sizeof(double*));
	for(int i=0; i<num_of_ffts; i++){
		*(fft_array+i) = (double*)malloc(window_size*sizeof(double));
	}

	double* window = (double*)malloc(window_size*sizeof(double));
	std::complex<double>* window_dft = (std::complex<double>*)malloc(n_fft*sizeof(std::complex<double>));
	int pointer = 0;

	for(int i=0; i<num_of_ffts; i++){
		for(int j=0; j<window_size; j++){
			*(window+j) = *(wave+pointer+j);
		}
		window = hamming_data(window, window_size);
		pointer += hop_length;
		window_dft = dft_rfft(window, window_size);
		for(int k=0; k<n_fft; k++){
			*(*(fft_array+i)+k) = abs(*(window_dft+k));
		}
	}

	return fft_array;
}

double** mel_spec(double *wave, int num_of_inputs, int window_size, int hop_length, int n_mels, int sample_rate){
	int num_of_ffts = (num_of_inputs - window_size)/hop_length;
	int n_fft = (window_size/2)+1;
	// Define the output array
	double** mel_spec_arr = (double**)malloc(num_of_ffts*sizeof(double*));
	for(int i=0; i<num_of_ffts; i++){
		*(mel_spec_arr+i) = (double*)malloc(n_mels*sizeof(double));
	}

	// Define the upper and lower frequency (20 hz --- samplerate)
	double lower_freq = 20;
	double upper_freq = sample_rate;
	// Calculate the mel scale
	double mel_lower_freq = 1125.*log(1.+lower_freq/700.);
	double mel_upper_freq = 1125.*log(1.+upper_freq/700.);
	// Calculate the mal gap
	double mel_gap = (mel_upper_freq - mel_lower_freq)/(n_mels+1);
/*
	// visualise the values
	std::cout<<mel_lower_freq<<std::endl;
	std::cout<<mel_upper_freq<<std::endl;
	std::cout<<mel_gap<<std::endl;
	std::cout<<"number of n_mel: "<<n_mels<<std::endl;
	std::cout<<"number of n_fft: "<<n_fft<<std::endl;
	std::cout<<"number of num_of_ffts: "<<num_of_ffts<<std::endl;
*/
	//Define filter bank
	double mel_value = mel_lower_freq;
	double* filterbank = (double*)malloc((n_mels+2)*sizeof(double));
	for(int i=0; i<n_mels+2; i++){
		// Revert mel scale to hz
		double mel_to_hz = 700*(exp(mel_value/1125)-1);
		// Round the frequency
		double rounding_freq = floor((n_fft*mel_to_hz)/sample_rate);
		*(filterbank+i) = rounding_freq;
		mel_value += mel_gap;
	}
/*
	// Visualise the  filterbank
	for(int i=0; i<n_mels+2; i++){
		std::cout<<filterbank[i]<<" ";
	}
	std::cout<<std::endl;
*/
	// Get the stft array
	double** stft_arr = stft(wave, num_of_inputs, window_size, hop_length);

	// Compute the mel spectrogram
	for(int i=0; i<num_of_ffts; i++){
		for(int m=1; m<n_mels+1; m++){
	                double mel_scale = 0.;
        	        double mel_weight = 0.;
                	for(int j=0; j<n_fft; j++){
                        	if(*(*(stft_arr+i)+j) < *(filterbank+m-1)){
                                	mel_weight = 0;
                        	}else if(*(*(stft_arr+i)+j)>=*(filterbank+m-1) && *(*(stft_arr+i)+j)<=*(filterbank+m)){
                                	mel_weight = (*(*(stft_arr+i)+j) - *(filterbank+m-1))/(*(filterbank+m) - *(filterbank+m-1));
                        	}else if(*(*(stft_arr+i)+j)>=*(filterbank+m) && *(*(stft_arr+i)+j)<=*(filterbank+m+1)){
                                	mel_weight = (*(filterbank+m+1) - *(*(stft_arr+i)+j))/(*(filterbank+m+1) - *(filterbank+m));
                        	}else if(*(*(stft_arr+i)+j) > *(filterbank+m+1)){
                                	mel_weight = 0;
                        	}
                        	mel_scale += *(*(stft_arr+i)+j)*mel_weight;
                	}
			*(*(mel_spec_arr+i)+m-1) = mel_scale;
		}
	}

	return mel_spec_arr;

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

