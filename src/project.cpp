#include "PLL.h"
//#include "dy4.h"
#include "filter.h"
#include "iofunc.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <vector>
#include <complex>
#include <cmath>

//#include "MonoBlock.cpp"

#define PI 3.14159265358979323846

//makefile
//project uses filter.cpp iofunc.cpp PLL.cpp project.cpp and filter.h iofunc.h PLL.h

//to run the file do
//make -f makefile
// ./project | aplay -c 2 -f S16_LE -r 48000
//cat file | ./project mode | aplay -c 2 -f S16_LE -r 48000

//rtl_sdr -f 107.1M –s 2.4M - | ./project mode | aplay -c 2 -f S16_LE -r 48000

//make this file pipeline live radio
// radio input | ./project | aplay

void readStdinBlockData(unsigned int num_samples, unsigned int block_id, std::vector<float> &block_data){
	std::vector<char> raw_data(num_samples);
	std::cin.read(reinterpret_cast<char*>(&raw_data[0]), num_samples*sizeof(char));
	for (int k=0; k<(int)num_samples; k++){
		block_data[k] = float(((unsigned char)raw_data[k]-128)/128.0);
	}
	}

int main(int argc, char* argv[])
{
  
  
  int mode = 0;
  // Mode Selection
  if (argc<2){
    std::cerr << "Operating in default mode 0 and mono" << std::endl;
  } else if (argc==2){
    mode=atoi(argv[1]);
    if (mode>3){
      std::cerr << "Wrong mode " << mode << std::endl;
      exit(1);}
  } else {
    std::cerr << "Usage: " << argv[0] << std::endl;
    std::cerr << "or " << std::endl;
    std::cerr << "Usage: " << argv[0] << " <mode>" << std::endl;
    std::cerr << "\t\t <mode> is a value from 0-3" << argv[0] << std::endl;
    exit(1);}

  std::cerr << "Operating in mode " << mode << std::endl;



	// RF variables

	float rf_Fs = 2400000.0;
	float rf_Fc = 100000.0;

	float audio_Fs = 48000.0;	// sample rate for our "assumed" audio (change as needed for 48 ksamples/sec audio files)
	float audio_Fb = 23000.0;
	float audio_Fe = 53000.0;

	float carrier_Fb = 18500.0;
	float carrier_Fe = 19500.0;

	//variables for PLL block process
	std::vector<float> pll_variables(5,0);

	pll_variables[0] = 0;
	pll_variables[1] = 0;
	pll_variables[2] = 1;
	pll_variables[3] = 0;
	pll_variables[4] = 0;

	float freq = 19000;
	float fs = 48000;
	float phaseadjust = 0.0;

	float normBandwidth = 0.01;
	float trigOffset = 0.0;
	float ncoScale = 2.0;

	unsigned short int rf_upSample = 1;
	unsigned short int rf_taps = 151;
	unsigned short int rf_decim = 10;

	unsigned short int audio_upSample = 1;
	unsigned short int audio_taps = 151;
	unsigned short int audio_decim = 5;
	float audio_Fc = 16000.0;
	int blockSize = 1024 * rf_decim * audio_decim * 2;

  if(mode == 0){
    rf_Fs = 2400000.0;
    rf_taps = 151;
    rf_decim = 10;
    // audio variables
    audio_Fs = 48000.0;
    audio_taps = 151;
    audio_decim = 5;
    
  }else if (mode == 1){
    rf_Fs = 1440000.0;
    rf_taps = 151;
    rf_decim = 5;
    blockSize = 1024 * rf_decim * audio_decim * 2;
    // audio variables
    audio_Fs = 48000.0;
    audio_taps = 151;
    audio_decim = 6;
    
  }else if (mode == 2){

    rf_Fs = 2400000.0;
    rf_taps = 151;
    rf_decim = 10;
    blockSize = 1024 * rf_decim * audio_decim * 2;

    // audio variables
    audio_Fs = 44100.0;
    audio_taps = 151;
    audio_upSample = 147;
    audio_decim = 800;

  }else if (mode == 3){
    rf_Fs = 1920000.0;
    rf_taps = 151;
    rf_decim = 5;
    blockSize = 1024 * rf_decim * audio_decim * 2;

    // audio variables
    audio_Fs = 44100.0;
    audio_taps = 151;
    audio_upSample = 147;
    audio_decim = 1280;
  }

      // impulse response (reuse code from the previous experiment)
	std::vector<float> rf_coeff;
	impulseResponseLPF(rf_Fs, rf_Fc, rf_taps, rf_coeff);
	std::vector<float> state_i(rf_taps - 1,0);
	std::vector<float> state_q(rf_taps - 1,0);
	float prevI = 0;
	float prevQ = 0;

	std::vector<float> filtered_i((blockSize/(2 * rf_decim)), 0);
	std::vector<float> filtered_q((blockSize/(2 * rf_decim)), 0);
	std::vector<float> block;

	unsigned short int startIndex_i = 0;
	unsigned short int startIndex_q = 0;


	// Stereo Specific Variables
	std::vector<float> stereo_extraction_state(audio_taps - 1,0);
	std::vector<float> stereo_extraction_coeff;
	impulseResponseBPF(audio_Fs, audio_Fb, audio_Fe, audio_taps, stereo_extraction_coeff);
	std::vector<float> carrier_coeff;
	impulseResponseBPF(audio_Fs, carrier_Fb, carrier_Fe, audio_taps, carrier_coeff);
	std::vector<float> stereo_coeff;
	impulseResponseLPF(48000, 19000, audio_taps, stereo_coeff);

	std::vector<float> stereo_state(audio_taps - 1, 0);
	std::vector<float> state_pll(audio_taps - 1,0);
	std::vector<float> stereo_data_final;

	unsigned short int startIndex_stereo = 0;
	unsigned short int startIndex_pll = 0;

	std::vector<float> pll_block(filtered_i.size(), 0);
	std::vector<float> pll;

  // End stereo specific

  // Start Mono Specific variables

      std::vector<float> mono_extraction_state(audio_taps - 1,0);
      unsigned short int allPass_taps = 75;
      std::vector<float> mono_extraction_coeff;
      impulseResponseLPF(audio_Fs, audio_Fc, audio_taps, mono_extraction_coeff);  //Match the delay
      std::vector<float> allPass_coeff;
      impulseResponseLPF(48000,20000, allPass_taps, allPass_coeff);
      unsigned short int startIndex_mono = 0;
      std::vector<float> allPass_State(allPass_taps - 1,0);
      std::vector<float> audio_data(blockSize);
  // End Mono Specific Variables
  
  for (unsigned int blockCount =0; ;blockCount++){
    
    //reading
    readStdinBlockData(blockSize,blockCount,audio_data);
    if (std::cin.rdstate()!=0){
      exit(1);
      }
    
    // Processing IQ samples to acquire them within 100khz range
    std::fill (filtered_i.begin(),filtered_i.end(),0);
    makeOddEvenSubList(block,audio_data,0,blockSize);
    blockProcessing(rf_coeff, block, state_i, rf_taps, filtered_i,startIndex_i,rf_decim);
    std::fill (filtered_q.begin(),filtered_q.end(),0);
    makeOddEvenSubList(block,audio_data,1,blockSize);
    blockProcessing(rf_coeff, block, state_q, rf_taps, filtered_q,startIndex_q,rf_decim);

    // Demodulating to merge IQ samples
    std::vector<float> demodulatedSignal ((int)(filtered_i.size()),0);
    fmDemod (demodulatedSignal, filtered_i, filtered_q,prevI,prevQ);
    std::vector<float> stereo_extraction_block((demodulatedSignal.size()/audio_decim)*audio_upSample, 0);
    std::vector<float> pll_block(stereo_extraction_block.size(), 0);

    if (audio_upSample == 1){
      blockProcessing(stereo_extraction_coeff, demodulatedSignal, stereo_extraction_state, audio_taps, stereo_extraction_block,startIndex_stereo, audio_decim);
      blockProcessing(carrier_coeff, demodulatedSignal, state_pll, audio_taps, pll_block,startIndex_pll, audio_decim);
    }else{
      blockResample(stereo_extraction_coeff, demodulatedSignal, stereo_extraction_state, audio_taps, stereo_extraction_block,audio_decim,audio_upSample);
      blockResample(carrier_coeff, demodulatedSignal, state_pll, audio_taps, pll_block, audio_decim,audio_upSample);
    }


    std::vector<float> mono_extraction_block((demodulatedSignal.size()/audio_decim)*audio_upSample, 0);


    if (audio_upSample == 1){
      blockProcessing(mono_extraction_coeff, demodulatedSignal, mono_extraction_state, audio_taps, mono_extraction_block,startIndex_mono, audio_decim);}
      else{blockResample(mono_extraction_coeff, demodulatedSignal, mono_extraction_state, audio_taps, mono_extraction_block,audio_decim,audio_upSample);}
      
    std::vector<float> mono_block_filtered((demodulatedSignal.size()/audio_decim), 0);
    blockConvolve(allPass_coeff, mono_extraction_block, allPass_State, audio_taps, mono_block_filtered);
    std::vector<float> pll_processed ((pll_block.size()),0);
    

    fmPLL(pll_processed, pll_block, freq, fs, ncoScale, phaseadjust, normBandwidth,pll_variables);
    std::vector<float> stereo_block ((pll_processed.size()),0);
    std::vector<float> filtered_stereo ((pll_processed.size()),0);
    for (int i = 0; i < pll_block.size(); i++){
	    stereo_block [i] = stereo_extraction_block[i] * pll_processed[i];
	  }
    blockConvolve(stereo_coeff, stereo_block, stereo_state, audio_taps, filtered_stereo);

    std::vector<float> stereoLeft ((filtered_stereo.size()),0);
    std::vector<float> stereoRight ((filtered_stereo.size()),0);
    std::vector<float> stereo ((filtered_stereo.size()*2),0);
  
    for (int i = 0; i < stereoLeft.size(); i++){
      
      stereoLeft [i] = (filtered_stereo[i] + mono_block_filtered[i])/2;
      stereoRight [i] = (mono_block_filtered[i] - filtered_stereo[i])/2;

  		stereo[2*i] = stereoLeft[i];
  		stereo[2*i+1] = stereoRight[i];
  	}


  
      std::vector<short int> final_data(stereo.size());
       for (unsigned int k=0;k<stereo.size();k++){
         if(std::isnan(stereo[k])) final_data[k]=0;
         else final_data[k] = static_cast<short int>(stereo[k]*16384);

       }
       
       fwrite(&final_data[0], sizeof(short int), final_data.size(), stdout);}
	
  


	return 0;
}
