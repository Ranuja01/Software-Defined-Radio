/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Copyright by Nicola Nicolici
Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>

#include "PLL.cpp"
//#include "MonoBlock.cpp"

#define PI 3.14159265358979323846

void makeSubList (std::vector<float> &subList, const std::vector<float> &list, int first, int last);

// function for computing the impulse response (reuse from previous experiment)
void impulseResponseLPF(float Fs, float Fc, unsigned short int num_taps, std::vector<float> &h)
{
	// allocate memory for the impulse response
	h.clear(); h.resize(num_taps, 0.0);

	// the rest of the code in this function is to be completed by you
	// based on your understanding and the Python code from the first lab

	float norm_cutoff = Fc/(Fs/2);

	for (int i = 0; i < num_taps; i++){
    //  ****TRY SPLITTING IF STATEMENT INTO TWO FOR LOOPS****
		if (i == (num_taps - 1)/2){
			h[i] = norm_cutoff;
		} else {
			h[i] = norm_cutoff * sin (PI * norm_cutoff * (i - (num_taps - 1)/2)) / (PI * norm_cutoff * (i - (num_taps - 1)/2));
		}
		h[i] *= pow(sin(PI * i / num_taps),2);

	}
}

void impulseResponseBPF(float Fs, float Fb, float Fe, unsigned short int num_taps, std::vector<float> &h)
{
	// allocate memory for the impulse response
	h.clear(); h.resize(num_taps, 0.0);

	// the rest of the code in this function is to be completed by you
	// based on your understanding and the Python code from the first lab

	float norm_center = ((Fe + Fb)/2)/(Fs/2);
  float norm_pass = (Fe - Fb)/(Fs/2);



	for (int i = 0; i < num_taps; i++){
    //  ****TRY SPLITTING IF STATEMENT INTO TWO FOR LOOPS****
		if (i == (num_taps - 1)/2){
			h[i] = norm_pass;
		} else {
			h[i] = norm_pass * sin (PI * norm_pass * (i - (num_taps - 1)/2)) / (PI * norm_pass * (i - (num_taps - 1)/2));
		}
    h[i] *= cos(i*PI*norm_center);
		h[i] *= pow(sin(PI * i / num_taps),2);

	}
}

// function for computing the impulse response (reuse from previous experiment)
void convolveFIR(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h)
{

	// allocate memory for the output (filtered) data
	y.clear(); y.resize(x.size()+h.size()-1, 0.0);
	for (int n = 0; n < y.size(); n++){
		for (int k = 0; k < h.size(); k++){
			if (n - k >= 0 && n - k < x.size()){
				y[n] += h[k] * x[n-k];
			}
		}
	}
	// the rest of the code in this function is to be completed by you
	// based on your understanding and the Python code from the first lab

}

// function to read audio data from a binary file that contains raw samples
// represented as 32-bit floats; we also assume two audio channels
// note: check the Python script that can prepare this type of files
// directly from .wav files
void read_audio_data(const std::string in_fname, std::vector<uint8_t> &audio_data)
{
	// file descriptor for the input to be read
	std::ifstream fdin(in_fname, std::ios::binary);
	if(!fdin) {
		std::cout << "File " << in_fname << " not found ... exiting\n";
		exit(1);
	} else {
		//std::cout << "Reading raw audio from \"" << in_fname << "\"\n";
    std::cout << "Read raw RF data from \"" << in_fname <<"\" in unsigned 8-bit format" << std::endl;
	}
	// search for end of file to count the number of samples to be read
	fdin.seekg(0, std::ios::end);
	// we assume the Python script has written data in 8-bit unsigned integer
	const unsigned int num_samples = fdin.tellg() / sizeof(uint8_t);

	// allocate memory space to store all the samples
	audio_data.clear(); audio_data.resize(num_samples);
	// back to the beginning of the file to read all samples at once
	fdin.seekg(0, std::ios::beg);
	// do a single read for audio data from the input file stream
	fdin.read(reinterpret_cast<char*>(&audio_data[0]), \
						num_samples*sizeof(uint8_t));

	// close the input file
	fdin.close();
}

// function to split an audio data where the left channel is in even samples
// and the right channel is in odd samples
void split_audio_into_channels(const std::vector<float> &audio_data, std::vector<float> &audio_left, std::vector<float> &audio_right)
{
	for (int i=0; i<(int)audio_data.size(); i++) {
		if (i%2 == 0)
			audio_left.push_back(audio_data[i]);
		else
			audio_right.push_back(audio_data[i]);
	}
}

// function to write audio data to a binary file that contains raw samples
// represented as 32-bit floats; we also assume two audio channels
// note: check the python script that can read this type of files
// and then reformat them to .wav files to be run on third-party players
void write_audio_data(const std::string out_fname, std::vector<float> &audio)
{
	// file descriptor for the output to be written
	if (audio.size() == 0) {
		std::cout << "Something got messed up with audio channels\n";
		std::cout << "They must have the same size ... exiting\n";
		exit(1);
	} else {
		std::cout << "Writing raw audio to \"" << out_fname << "\"\n";
	}
	std::ofstream fdout(out_fname, std::ios::binary);
	for (int i=0; i<(int)audio.size(); i++) {
		// we assume we have handled a stereo audio file
		// hence, we must interleave the two channels
		// (change as needed if testing with mono files)
    audio[i] = (uint16_t)(audio[i]/2 * 32767);
		fdout.write(reinterpret_cast<const char*>(&audio[i]),\
								sizeof(audio[i]));
		//fdout.write(reinterpret_cast<const char*>(&audio_right[i]),\
								sizeof(audio_right[i]));
	}
	fdout.close();
}

void blockConvolve(std::vector<float> &h, const std::vector<float> &block, std::vector<float> &state, int num_taps, std::vector<float> &filtered_block){

  std::cout << filtered_block.size() << " " << h.size() << "  " << block.size() << std::endl;
  for (int n = 0; n < block.size(); n++){
		for (int k = 0; k < h.size(); k++){
			if (n - k >=0){
				if (n - k < block.size()){
					filtered_block[n] += h[k] * block[n-k];
        //  std::cout<< "if   " << filtered_block[n] << std::endl;

				}
			}else {
		  	if (n - k + num_taps - 1 < state.size()){
					filtered_block[n] += h[k] * state[(n - k) + num_taps - 1];
          //          std::cout<< "else   " << filtered_block[n] << std::endl;
				}
			}
		}
	}

  std::cout << "blockconv" << std::endl;
	makeSubList(state,block,block.size() - num_taps + 1, block.size());
}

//blockProcessing(rf_coeff, block, state_i, rf_taps, filtered_i,startIndex_i,rf_decim);
void blockProcessing(std::vector<float> &h, const std::vector<float> &block, std::vector<float> &state, int num_taps, std::vector<float> &filtered_block, unsigned short int &startIndex, int dRate){
	int endIndex = 0;
	for (int n = startIndex; n < block.size(); n+=dRate){
		//std::cout << "bbbb: " << n << std::endl;
		for (int k = 0; k < h.size(); k++){
			if (n - k >= 0){
				if (n - k < block.size()){
					filtered_block[(n - startIndex)/dRate] += h[k] * block[n-k];

				}
			}else {
		  	if (n - k + num_taps - 1 < state.size()){
					filtered_block[(n - startIndex)/dRate] += h[k] * state[(n - k) + num_taps - 1];
				}
			}





		}
		endIndex = n;
	}
	//startIndex = 1 + (endIndex + block.size() % dRate) - block.size();
	startIndex = (endIndex - block.size()) + dRate;
/*
	for (int i = filtered_block.size() - 5140; i < filtered_block.size() - 5100;i++){
		//std::cout << "aaaa: " << filtered_block[i] << std::endl;
		//std::cout << "bbbb: " <<  block[i] << std::endl;
	//	std::cout << "cccc: " <<  h[i] << std::endl;

}
	/*
	for (int i = 51000; i < block.size();i++){
	//	std::cout << "aaaa: " << filtered_block[i] << std::endl;
		//std::cout << "bbbb: " <<  block[i] << std::endl;

			std::cout << i << " dddd: " <<  block[i] << std::endl;
			std::cout << (i)/dRate << " eeee: " <<  filtered_block[(i)/dRate] << std::endl;
	}*/
	//std::cout << "asdfghjkl: " << filtered_block[filtered_block.size() - 2] << std::endl;
	//std::cout << startIndex << std::endl;
	//startIndex = 0;
	makeSubList(state,block,block.size() - num_taps + 1, block.size());
}

void blockResample(std::vector<float> &h, const std::vector<float> &block, std::vector<float> &state, int num_taps, std::vector<float> &filtered_block, int dRate, int uRate){
	int p = 0;
  int j = 0;

	for (int n = state.size(); n < block.size() + state.size(); n++){
		//std::cout << "bbbb: " << n << std::endl;
    p = (n * dRate)%uRate;
		for (int k = p; k < h.size(); k+= uRate){

      j = (int)((n*dRate - k)/uRate);

			if (j >= 0) {
        if (j < state.size()){
          filtered_block[n - state.size()] += h[k] * state[j - state.size()];
        } else if (j < block.size()){
          filtered_block[n - state.size()] += h[k] * block[j - state.size()];
        }
      }

		}

	}
	//startIndex = 1 + (endIndex + block.size() % dRate) - block.size();
/*
	for (int i = filtered_block.size() - 5140; i < filtered_block.size() - 5100;i++){
		//std::cout << "aaaa: " << filtered_block[i] << std::endl;
		//std::cout << "bbbb: " <<  block[i] << std::endl;
	//	std::cout << "cccc: " <<  h[i] << std::endl;

}
	/*
	for (int i = 51000; i < block.size();i++){
	//	std::cout << "aaaa: " << filtered_block[i] << std::endl;
		//std::cout << "bbbb: " <<  block[i] << std::endl;

			std::cout << i << " dddd: " <<  block[i] << std::endl;
			std::cout << (i)/dRate << " eeee: " <<  filtered_block[(i)/dRate] << std::endl;
	}*/
	//std::cout << "asdfghjkl: " << filtered_block[filtered_block.size() - 2] << std::endl;
	//std::cout << startIndex << std::endl;
	//startIndex = 0;
	makeSubList(state,block,block.size() - num_taps + 1, block.size());
}

void makeSubList (std::vector<float> &subList, const std::vector<float> &list, int first, int last){
	subList.clear();
	for (int i = first; i < last; i++){


		//	std::cout << "uytrew: " << list[i] << std::endl;


		subList.push_back(list[i]);
	}

}

void makeOddEvenSubList (std::vector<float> &subList, const std::vector<float> &list, int first, int last){
	subList.clear();
	for (int i = first; i < last; i+=2){
		subList.push_back(list[i]);
	}

}

void fmDemod (std::vector<float> &demodulatedSignal, const std::vector<float> &I, const std::vector<float> &Q, float &prevI, float &prevQ){
	//demodulatedSignal.clear();
	std::fill (demodulatedSignal.begin(),demodulatedSignal.end(),0);
	demodulatedSignal[0] = 0;
	std::cout <<"I "<< I[0] <<std::endl;
	std::cout <<"Q "<<Q[0] <<std::endl;
	std::cout <<"prevI "<<prevI <<std::endl;
	std::cout <<"prevQ "<<prevQ <<std::endl;

	for (int k = 0; k < I.size(); k++){
		if (!( (I[k] * I[k]) + (Q[k] * Q[k]) == 0)){
			demodulatedSignal[k] = (1.0/( (I[k] * I[k]) + (Q[k] * Q[k]) ) ) * (I[k] * (Q[k] - prevQ) - Q[k] * (I[k] - prevI));


		}else {
			demodulatedSignal[k] = 0;
		}
		prevI = I[k];
		prevQ = Q[k];
	}

	std::cout <<demodulatedSignal[0] <<std::endl;
}


int main()
{
  int mode = 0;
	// assume the wavio.py script was run beforehand to produce a binary file
	const std::string in_fname = "iq_samples.raw";
	// declare vector where the audio data will be stored
	std::vector<uint8_t> iq_data;

	// note: we allocate memory for audio_data from within this read function
	read_audio_data(in_fname, iq_data);

	std::vector<float> audio_data(iq_data.size(),0);
	//float a = iq_data[0]<<24;
	/*
	for (int i = 0; i < 50; i++){
		std::cout << (float)iq_data[i] <<std::endl;
	}*/

  for (int i = 0; i < iq_data.size(); i++){
    audio_data[i] = ((float)iq_data[i] - 128.0)/128.0;
  }
	std::cout << iq_data.size() <<std::endl;
	std::cout << audio_data.size() <<std::endl;

  // RF variables

  float rf_Fs = 2400000.0;
  float rf_Fc = 100000.0;

  float audio_Fs = 48000.0;	// sample rate for our "assumed" audio (change as needed for 48 ksamples/sec audio files)
  float audio_Fb = 23000.0;
  float audio_Fe = 53000.0;

  float RDS_Fb = 54000.0;
  float RDS_Fe = 60000.0;

  float RDS_carrier_Fb = 113500.0;
  float RDS_carrier_Fe = 114500.0;

	//variables for PLL block process
	std::vector<float> pll_variables(5,0);

	pll_variables[0] = 0;
	pll_variables[1] = 0;
	pll_variables[2] = 1;
	pll_variables[3] = 0;
	pll_variables[4] = 0;

/*
  float integrator = 0;
  float phaseEst = 0;
  float feedbackI = 1.0;
  float feedbackQ = 0.0;
  float trigOffset = 0.0;
  float ncoScale = 2.0;
*/

	float freq = 114000;
	float fs = 48000;
	float phaseadjust = 0.0;

	float normBandwidth = 0.01;
	float trigOffset = 0.0;
	float ncoScale = 0.5;

  unsigned short int rf_upSample = 1;
  unsigned short int rf_taps = 151;
  unsigned short int rf_decim = 10;

  unsigned short int audio_upSample = 1;
  unsigned short int audio_taps = 151;
  unsigned short int audio_decim = 5;

  if(mode == 0){
    rf_Fs = 2400000.0;
  	//rf_Fc = 100000.0;
  	rf_taps = 151;
    rf_decim = 10;

  	// audio variables
  	audio_Fs = 48000.0;	// sample rate for our "assumed" audio (change as needed for 48 ksamples/sec audio files)
  	//audio_Fc = 16000.0;	// cutoff frequency (explore ... but up-to Nyquist only!)
  	audio_taps = 151;
    audio_decim = 5;
  }else if (mode == 1){
    rf_Fs = 1440000.0;
  	//rf_Fc = 100000.0;
  	rf_taps = 151;
    rf_decim = 5;

  	// audio variables
  	audio_Fs = 48000.0;	// sample rate for our "assumed" audio (change as needed for 48 ksamples/sec audio files)
  	//audio_Fc = 16000.0;	// cutoff frequency (explore ... but up-to Nyquist only!)
  	audio_taps = 151;
    audio_decim = 6;
  }else if (mode == 2){

    rf_Fs = 2400000.0;
  	//rf_Fc = 100000.0;
  	rf_taps = 151;
    rf_decim = 10;
  	// audio variables
  	audio_Fs = 44100.0;	// sample rate for our "assumed" audio (change as needed for 48 ksamples/sec audio files)
  	//audio_Fc = 16000.0;	// cutoff frequency (explore ... but up-to Nyquist only!)
  	audio_taps = 151;
    audio_upSample = 147;
    audio_decim = 800;

  }else if (mode == 3){
    rf_Fs = 1920000.0;
    //rf_Fc = 100000.0;
    rf_taps = 151;

    rf_decim = 5;

    // audio variables
    audio_Fs = 44100.0;	// sample rate for our "assumed" audio (change as needed for 48 ksamples/sec audio files)
    //audio_Fc = 16000.0;	// cutoff frequency (explore ... but up-to Nyquist only!)
    audio_taps = 151;
    audio_upSample = 147;
    audio_decim = 1280;
  }

  int blockSize = 1024 * rf_decim * audio_decim * 2;
  int blockCount = 0;
  int input = 2;

	// impulse response (reuse code from the previous experiment)
	std::vector<float> rf_coeff;
	impulseResponseLPF(rf_Fs, rf_Fc, rf_taps, rf_coeff);

  std::vector<float> audio_coeff;
	impulseResponseBPF(audio_Fs, RDS_Fb, RDS_Fe, audio_taps, audio_coeff);

  std::vector<float> pll_coeff;
	impulseResponseBPF(audio_Fs, RDS_carrier_Fb, RDS_carrier_Fe, audio_taps, pll_coeff);

  std::vector<float> allPass_coeff;
	impulseResponseBPF(audio_Fs, 0, 100000, audio_taps, allPass_coeff);

	std::vector<float> RDS_coeff;
	impulseResponseLPF(48000, 3000, audio_taps, RDS_coeff);

// MAYBE WE NEED HANN WINDOW

  std::vector<float> state_i(rf_taps - 1,0);
  std::vector<float> state_q(rf_taps - 1,0);
  std::vector<float> state_audio(audio_taps - 1,0);
  std::vector<float> state_allPass(audio_taps - 1,0);
  std::vector<float> state_pll(audio_taps - 1,0);
	std::vector<float> RDS_state(audio_taps - 1, 0);
  //std::vector<float> filtered_block((blockSize/rf_decim), 0);

  float prevI = 0;
  float prevQ = 0;

	//std::vector<float> audio_state;
	std::vector<float> RDS_data_final;

  std::vector<float> filtered_i((blockSize/(2 * rf_decim)), 0);
	std::vector<float> filtered_q((blockSize/(2 * rf_decim)), 0);
  std::vector<float> block;

	unsigned short int startIndex_i = 0;
	unsigned short int startIndex_q = 0;
	unsigned short int startIndex_audio = 0;
  unsigned short int startIndex_pll = 0;



  //std::vector<float> RDS_filtered(RDS_block.size(), 0);
  //std::vector<float> pll_block(RDS_block.size(), 0);

//  std::vector<float> pll;


  while ((blockCount + 1) * blockSize < audio_data.size()){

    std::cout <<"Processing block: " << blockCount << std::endl;

    std::fill (filtered_i.begin(),filtered_i.end(),0);
    makeOddEvenSubList(block,audio_data,blockCount*blockSize,(blockCount + 1)*blockSize);
		//std::cout << "aaaaaa: " << filtered_i.size() << std::endl;
	//	std::cout << "bbbbbb: " << block.size() << std::endl;
		//std::cout << "bbbbbb: " << blockSize << std::endl;

    blockProcessing(rf_coeff, block, state_i, rf_taps, filtered_i,startIndex_i,rf_decim);
/*

			//std::cout << "filtered block: " << filtered_i[n] << std::endl;
			std::cout << "h[k]: " << rf_coeff[rf_coeff.size() - 1] << std::endl;
			std::cout << "block[n-k]: " << block[block.size() - 1] << std::endl;
			std::cout << "state[n-k + num_taps]: " <<  state_i[rf_taps - 5] << std::endl;
*/
    //filtered_i.insert( filtered_i.end(), filtered_block.begin(), filtered_block.end() );
    std::fill (filtered_q.begin(),filtered_q.end(),0);
    makeOddEvenSubList(block,audio_data,blockCount*blockSize + 1,(blockCount + 1)*blockSize);
    blockProcessing(rf_coeff, block, state_q, rf_taps, filtered_q,startIndex_q,rf_decim);
    //filtered_q.insert(filtered_q.end(), filtered_block.begin(), filtered_block.end() );



	//	std::cout << "blocksize: " << blockSize<< std::endl;
	//	std::cout << "filteredSize: " << filtered_i.size()<< std::endl;
	//	std::cout << "bfd: " << filtered_i[filtered_i.size() - 2] << std::endl;
		std::vector<float> demodulatedSignal ((int)(filtered_i.size()),0);
		fmDemod (demodulatedSignal, filtered_i, filtered_q,prevI,prevQ);
		for (int i = 0; i < 5; i++){

		//	std::cout << "bfd: " << demodulatedSignal[i] << std::endl;

		}
		//std::cout << "qwertyu: " << block[1] << std::endl;


		std::vector<float> audio_block((demodulatedSignal.size()/audio_decim)*audio_upSample, 0);
		std::vector<float> pll_block(audio_block.size(), 0);
/*
		if (blockCount == 0){
			for (int i = 0; i < audio_taps - 1; i++){
				audio_state.push_back(0);
			}
		}*/
		for (int i = 0; i < 5; i++){

			std::cout << "abc: " << demodulatedSignal[i] << std::endl;

		}
		//std::fill (filtered_block.begin(),filtered_block.end(),0);
		// blockProcessing(rf_coeff, block, state_q, rf_taps, filtered_q,startIndex_q,rf_decim);
    if (audio_upSample == 1){
      blockProcessing(audio_coeff, demodulatedSignal, state_audio, audio_taps, audio_block,startIndex_audio, audio_decim);
      //blockConvolve(allPass_coeff, RDS_block, state_allPass, audio_taps, RDS_filtered);
      for (int i = 0; i < pll_block.size(); i++){
        pll_block[i] *= pll_block[i];
      }
    }else{
      blockResample(audio_coeff, demodulatedSignal, state_audio, audio_taps, audio_block,audio_decim,audio_upSample);
      //blockConvolve(allPass_coeff, RDS_block, state_allPass, audio_taps, RDS_filtered);

      for (int i = 0; i < pll_block.size(); i++){
        pll_block[i] *= pll_block[i];
      }
    }
		std::vector<float> pll_filtered (pll_block.size(), 0);
// ALL pass filter
    blockConvolve(allPass_coeff, audio_block, state_allPass, audio_taps, audio_block);
// Bandpass filter
    blockConvolve(pll_coeff, pll_block, state_pll, audio_taps, pll_filtered);

		std::vector<float> pll_processed (pll_block.size(), 0);
    fmPLL(pll_processed,pll_filtered, freq, fs, ncoScale, phaseadjust, normBandwidth, pll_variables);
		//fmPLL(pll_processed, pll_block, freq, fs, ncoScale, phaseadjust, normBandwidth,pll_variables);
		std::vector<float> RDS_block ((pll_processed.size()),0);
    std::vector<float> filtered_RDS ((pll_processed.size()),0);

		for (int i = 0; i < pll_block.size(); i++){
	    RDS_block [i] = audio_block[i] * pll_processed[i];
	  }
		blockConvolve(RDS_coeff, RDS_block, RDS_state, audio_taps, filtered_RDS);

		RDS_data_final.insert(RDS_data_final.end(), audio_block.begin(), audio_block.end() );
  //  pll.insert(pll.end(), pll_block.begin(), pll_block.end());

		blockCount += 1;
  }


  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	return 0;
}

// FOR RATIONAL RESAMPLER, TURN SIGNAL INTO 57KS
// FOR DEMODULATION FIND PYTHON FILE FOR RRCF
