/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Copyright by Nicola Nicolici
Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#include "PLL.h"
#include "dy4.h"
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

//make this file pipeline live radio
// radio input | ./project | aplay

void mono(std::vector<float> &audio_data_final)
{

  int mode = 0;


  // assume the wavio.py script was run beforehand to produce a binary file
  const std::string in_fname = "stereo.raw";

  // declare vector where the audio data will be stored
  std::vector<uint8_t> iq_data;
  std::vector <short int> play;

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
  unsigned short int rf_taps = 151;
  unsigned short int rf_decim = 10;

  // audio variables
  float audio_Fs = 48000.0;	// sample rate for our "assumed" audio (change as needed for 48 ksamples/sec audio files)
  float audio_Fc = 16000.0;	// cutoff frequency (explore ... but up-to Nyquist only!)

  unsigned short int rf_upSample = 1;
    unsigned short int audio_upSample = 1;
  unsigned short int audio_taps = 151;
  unsigned short int audio_decim = 5;

  int blockSize = 1024 * rf_decim * audio_decim * 2;
  int blockCount = 0;
  int input = 2;

  if(mode == 0){
     rf_Fs = 2400000.0;
   	rf_Fc = 100000.0;
   	rf_taps = 151;
     rf_decim = 10;

   	// audio variables
   	audio_Fs = 48000.0;	// sample rate for our "assumed" audio (change as needed for 48 ksamples/sec audio files)
   	audio_Fc = 16000.0;	// cutoff frequency (explore ... but up-to Nyquist only!)
   	audio_taps = 151;
     audio_decim = 5;
   }else if (mode == 1){
     rf_Fs = 1440000.0;
   	rf_Fc = 100000.0;
   	rf_taps = 151;
     rf_decim = 5;

   	// audio variables
   	audio_Fs = 48000.0;	// sample rate for our "assumed" audio (change as needed for 48 ksamples/sec audio files)
   	audio_Fc = 16000.0;	// cutoff frequency (explore ... but up-to Nyquist only!)
   	audio_taps = 151;
     audio_decim = 6;
 		blockSize = 1024 * rf_decim  * audio_decim * 4;
   }else if (mode == 2){
 		//std::cout << "uytghgtghytg" << std::endl;
     rf_Fs = 2400000.0;
   	rf_Fc = 100000.0;
   	rf_taps = 151;
     rf_decim = 10;
   	// audio variables
   	audio_Fs = 44100.0;	// sample rate for our "assumed" audio (change as needed for 48 ksamples/sec audio files)
   	audio_Fc = 16000.0;	// cutoff frequency (explore ... but up-to Nyquist only!)
   	audio_taps = 151;

     audio_upSample = 147;
     audio_decim = 800;
 		blockSize = 8 * rf_decim  * audio_decim * 2;

   }else if (mode == 3){
     rf_Fs = 1920000.0;
     rf_Fc = 100000.0;
     rf_taps = 151;

     rf_decim = 5;

     // audio variables
     audio_Fs = 44100.0;	// sample rate for our "assumed" audio (change as needed for 48 ksamples/sec audio files)
     audio_Fc = 16000.0;	// cutoff frequency (explore ... but up-to Nyquist only!)
     audio_taps = 151;
     audio_upSample = 147;
     audio_decim = 1280;
   }



  // impulse response (reuse code from the previous experiment)
  std::vector<float> rf_coeff;
  impulseResponseLPF(rf_Fs, rf_Fc, rf_taps, rf_coeff);

  std::vector<float> audio_coeff;
  impulseResponseLPF(audio_Fs, audio_Fc, audio_taps, audio_coeff);

  std::vector<float> allPass_coeff;

	impulseResponseLPF(48000,20000, 75, allPass_coeff);

// MAYBE WE NEED HANN WINDOW

  std::vector<float> state_i(rf_taps - 1,0);
  std::vector<float> state_q(rf_taps - 1,0);
  std::vector<float> state_audio(audio_taps - 1,0);
  std::vector<float> filtered_block((blockSize/rf_decim), 0);
  std::vector<float> allPass_State(75 - 1,0);

  float prevI = 0;
  float prevQ = 0;

//	std::vector<float> audio_state;
  //std::vector<float> audio_data_final;

  std::vector<float> filtered_i((blockSize/(2 * rf_decim)), 0);
  std::vector<float> filtered_q((blockSize/(2 * rf_decim)), 0);
  std::vector<float> block;

  unsigned short int startIndex_i = 0;
  unsigned short int startIndex_q = 0;
  unsigned short int startIndex_audio = 0;

  while ((blockCount + 1) * blockSize < audio_data.size()){

    std::cout <<"Processing block: " << blockCount << std::endl;

    std::fill (filtered_i.begin(),filtered_i.end(),0);

    makeOddEvenSubList(block,audio_data,blockCount*blockSize,(blockCount + 1)*blockSize);
    //std::cout << "aaaaaa: " << filtered_i.size() << std::endl;
  //	std::cout << "bbbbbb: " << block.size() << std::endl;
    //std::cout << "bbbbbb: " << blockSize << std::endl;
    blockProcessing(rf_coeff, block, state_i, rf_taps, filtered_i,startIndex_i,rf_decim);
    for (int i = block.size() - rf_taps + 1 ; i < block.size(); i++){

  //		std::cout << "bfd: " << block[i] << std::endl;

    }
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


    std::vector<float> demodulatedSignal ((int)(filtered_i.size()),0);
  //	std::cout << "blocksize: " << blockSize<< std::endl;
  //	std::cout << "filteredSize: " << filtered_i.size()<< std::endl;
  //	std::cout << "bfd: " << filtered_i[filtered_i.size() - 2] << std::endl;
    fmDemod (demodulatedSignal, filtered_i, filtered_q,prevI,prevQ);

    //std::cout << "qwertyu: " << block[1] << std::endl;


    std::vector<float> audio_block((demodulatedSignal.size()/audio_decim)*audio_upSample, 0);


    if (audio_upSample == 1){
      blockProcessing(audio_coeff, demodulatedSignal, state_audio, audio_taps, audio_block,startIndex_audio, audio_decim);
			//std::cout << "FFF: " << audio_block.size() << std::endl;
		}else{
      blockResample(audio_coeff, demodulatedSignal, state_audio, audio_taps, audio_block,audio_decim,audio_upSample);
		}
    //write_audio_data(audio_block, audio_Fs/2, play);
    std::vector<float> audio_block_filtered((demodulatedSignal.size()/audio_decim), 0);
    for (int i = 0; i < 5; i++){

    //  std::cout <<"mono: " << audio_block[i] << std::endl;

    }
  //  std::vector<float> audio_block_df((demodulatedSignal.size()/audio_decim), 0);
  //  blockConvolve(allPass_coeff, audio_block, allPass_State, audio_taps, audio_block_filtered);
  //  blockConvolve(allPass_coeff, audio_block_filtered, allPass_State, audio_taps, audio_block_df);
    audio_data_final.insert(audio_data_final.end(), audio_block.begin(), audio_block.end() );
    //audio_data_final.insert(audio_data_final.end(), audio_block.begin(), audio_block.end() );

    blockCount += 1;
  }

}

int main()
{
  int mode = 0;

  const std::string out_fname = "fmMonoBlock(cpp).wav";
  std::ofstream fdout(out_fname, std::ios::out | std::ios::binary);
  	std::vector <short int> play (stereo_data_final.size(), 0);
  // Read file

	// assume the wavio.py script was run beforehand to produce a binary file
	const std::string in_fname = "stereo.raw";
	// declare vector where the audio data will be stored
	std::vector<uint8_t> iq_data;

	// note: we allocate memory for audio_data from within this read function
	read_audio_data(in_fname, iq_data);

	std::vector<float> audio_data(iq_data.size(),0);

  for (int i = 0; i < iq_data.size(); i++){
    audio_data[i] = ((float)iq_data[i] - 128.0)/128.0;
  }



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
  int blockCount = 0;

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
  impulseResponseLPF(audio_Fs, audio_Fc, audio_taps, mono_extraction_coeff);

  std::vector<float> allPass_coeff;
	impulseResponseLPF(48000,20000, allPass_taps, allPass_coeff);

  unsigned short int startIndex_mono = 0;

  std::vector<float> allPass_State(allPass_taps - 1,0);
  // End Mono Specific Variables


  while ((blockCount + 1) * blockSize < audio_data.size()){

    //std::cout <<"Processing block: " << blockCount << std::endl;

    // Processing IQ samples to acquire them within 100khz range
    std::fill (filtered_i.begin(),filtered_i.end(),0);
    makeOddEvenSubList(block,audio_data,blockCount*blockSize,(blockCount + 1)*blockSize);

    blockProcessing(rf_coeff, block, state_i, rf_taps, filtered_i,startIndex_i,rf_decim);

    std::fill (filtered_q.begin(),filtered_q.end(),0);
    makeOddEvenSubList(block,audio_data,blockCount*blockSize + 1,(blockCount + 1)*blockSize);
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
      blockProcessing(mono_extraction_coeff, demodulatedSignal, mono_extraction_state, audio_taps, mono_extraction_block,startIndex_mono, audio_decim);

		}else{
      blockResample(mono_extraction_coeff, demodulatedSignal, mono_extraction_state, audio_taps, mono_extraction_block,audio_decim,audio_upSample);
		}
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
  //  std::cout <<"sdfghj block: " << blockCount << std::endl;
    for (int i = 0; i < stereoLeft.size(); i++){
      //std::cout <<"lll block: " << blockCount << std::endl;
      stereoLeft [i] = (filtered_stereo[i] + mono_block_filtered[i])/2;
      stereoRight [i] = (mono_block_filtered[i] - filtered_stereo[i])/2;

  		stereo[2*i] = stereoLeft[i];
  		stereo[2*i+1] = stereoRight[i];
  	}
/*
    for (int i = 0; i < 15; i++){


      std::cout <<"left: " << stereoLeft[i] << std::endl;

    }

    for (int i = 0; i < 15; i++){


      std::cout <<"right: " << stereoRight[i] << std::endl;

    }*/
    //std::cout <<"poiuytr block: " << blockCount << std::endl;
//    stereo_data_final.insert(stereo_data_final.end(), stereo.begin(), stereo.end() );


  	write_stereo_data(stereo_data_final, audio_Fs/2, play);

  	//write to file and play in terminal

    fwrite(&play[0], sizeof(short int), play.size(), stdout);

		blockCount += 1;
  }
/*
  for (int i = 0; i < 15; i++){


    std::cout <<"stereo final: " << stereo_data_final[i] << std::endl;

  }*/
	//g++ (filename).cpp -o (file)
	//     ./(file) | aplay -c 2 -f S16_LE -r 48000

//	std::vector<float> monoSignal;
//	mono(monoSignal);

  fdout.close();

	return 0;
}
