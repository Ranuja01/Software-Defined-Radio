/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Copyright by Nicola Nicolici
Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

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
#include <thread>
#include <queue>
#include <chrono>
//#include "MonoBlock.cpp"


#define PI 3.14159265358979323846

void readStdinBlockData(unsigned int num_samples, unsigned int block_id, std::vector<float> &block_data){
	std::vector<char> raw_data(num_samples);
	std::cin.read(reinterpret_cast<char*>(&raw_data[0]), num_samples*sizeof(char));
	for (int k=0; k<(int)num_samples; k++){
		block_data[k] = float(((unsigned char)raw_data[k]-128)/128.0);
	}
}

//std::ref(rf_coeff),std::ref(rf_taps), std::ref(rf_decim),std::ref(blockSize), std::ref(demod_Q)
void rfThread(std::vector<float> &rf_coeff, unsigned short int &rf_taps, unsigned short int &rf_decim, int &blockSize,std::queue<float> &demod_Q,std::queue<float> &rds_demod_Q){


  int blockCount = 0;
  std::vector<float> state_i(rf_taps - 1,0);
  std::vector<float> state_q(rf_taps - 1,0);

  float prevI = 0;
  float prevQ = 0;

  std::vector<float> filtered_i((blockSize/(2 * rf_decim)), 0);
	std::vector<float> filtered_q((blockSize/(2 * rf_decim)), 0);
  std::vector<float> block;

  unsigned short int startIndex_i = 0;
	unsigned short int startIndex_q = 0;
  std::vector<float> audio_data(blockSize);
  std::vector<float> demodulatedSignal ((int)(filtered_i.size()),0);

  for (unsigned int blockCount =0; ;blockCount++){
    // Processing IQ samples to acquire them within 100khz range
    //std::cout <<"RF block: " << blockCount << std::endl;

    readStdinBlockData(blockSize,blockCount,audio_data);
    if (std::cin.rdstate()!=0){
      break;
    }

    std::fill (filtered_i.begin(),filtered_i.end(),0);
    makeOddEvenSubList(block,audio_data,0,blockSize);
    blockProcessing(rf_coeff, block, state_i, rf_taps, filtered_i,startIndex_i,rf_decim);
    std::fill (filtered_q.begin(),filtered_q.end(),0);
    makeOddEvenSubList(block,audio_data,1,blockSize);
    blockProcessing(rf_coeff, block, state_q, rf_taps, filtered_q,startIndex_q,rf_decim);

    fmDemod (demodulatedSignal, filtered_i, filtered_q,prevI,prevQ);


    for (int i = 0; i < demodulatedSignal.size(); i++){
    //  std::cout <<"dmSignal " << i << " " << demodulatedSignal[i] << std::endl;
      demod_Q.push(demodulatedSignal[i]);
			rds_demod_Q.push(demodulatedSignal[i]);
    }

    //blockCount++;
  }

}

//std::thread audio_thread(audioThread,std::ref(stereo_extraction_coeff),std::ref(carrier_coeff), std::ref(mono_extraction_coeff), std::ref(allPass_coeff),std::ref(stereo_coeff),std::ref(audio_taps), std::ref(rf_decim), std::ref(audio_decim), std::ref(audio_upSample),blockSize,std::ref(audio_data), std::ref(demod_Q),std::ref(stereo_data_final));
void audioThread(std::vector<float> &stereo_extraction_coeff, std::vector<float> &carrier_coeff, std::vector<float> &mono_extraction_coeff, std::vector<float> &allPass_coeff, std::vector<float> &stereo_coeff , unsigned short int &audio_taps, unsigned short int &rf_decim, unsigned short int &audio_decim, unsigned short int &audio_upSample, int &blockSize, std::queue<float> &demod_Q, std::vector<float> &stereo_data_final, float &audio_Fs){

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
  int blockCount = 0;

  std::vector<float> stereo_state(audio_taps - 1, 0);
  std::vector<float> state_pll(audio_taps - 1,0);

	unsigned short int startIndex_stereo = 0;
  unsigned short int startIndex_pll = 0;

  std::vector<float> stereo_extraction_state(audio_taps - 1,0);

  std::vector<float> mono_extraction_state(audio_taps - 1,0);

  unsigned short int allPass_taps = 74;
  unsigned short int startIndex_mono = 0;

  std::vector<float> allPass_State(allPass_taps - 1,0);

  std::vector<float> demodulatedSignal (blockSize/(2 * rf_decim),0);
  std::vector<float> stereo_extraction_block((demodulatedSignal.size()/audio_decim)*audio_upSample, 0);
  std::vector<float> pll_block(stereo_extraction_block.size(), 0);
  std::vector<float> mono_extraction_block((demodulatedSignal.size()/audio_decim)*audio_upSample, 0);
  std::vector<float> mono_block_filtered((demodulatedSignal.size()/audio_decim), 0);
  std::vector<float> pll_processed ((pll_block.size()),0);
  std::vector<float> stereo_block ((pll_processed.size()),0);
  std::vector<float> filtered_stereo ((pll_processed.size()),0);
  std::vector<float> stereoLeft ((filtered_stereo.size()),0);
  std::vector<float> stereoRight ((filtered_stereo.size()),0);
  std::vector<float> stereo ((filtered_stereo.size()*2),0);

  for (unsigned int blockCount =0; ;blockCount++){
		//std::cout <<"Audio block: " << blockCount << std::endl;
    if (demod_Q.empty()){
			//std::cout << "kjhygtf: " << std::endl;
      exit(1);
		}
std::cout << "AA: " << std::endl;
    for (int i = 0; i < demodulatedSignal.size(); i++){

      demodulatedSignal[i] = demod_Q.front();
      demod_Q.pop();

    }
//std::cout << "BB: " << std::endl;
    std::fill (stereo_extraction_block.begin(),stereo_extraction_block.end(),0);
    std::fill (pll_block.begin(),pll_block.end(),0);

    if (audio_upSample == 1){
      blockProcessing(stereo_extraction_coeff, demodulatedSignal, stereo_extraction_state, audio_taps, stereo_extraction_block,startIndex_stereo, audio_decim);
      blockProcessing(carrier_coeff, demodulatedSignal, state_pll, audio_taps, pll_block,startIndex_pll, audio_decim);
    }else{
      blockResample(stereo_extraction_coeff, demodulatedSignal, stereo_extraction_state, audio_taps, stereo_extraction_block,audio_decim,audio_upSample);
      blockResample(carrier_coeff, demodulatedSignal, state_pll, audio_taps, pll_block, audio_decim,audio_upSample);
    }

    std::fill (mono_extraction_block.begin(),mono_extraction_block.end(),0);

    if (audio_upSample == 1){
      blockProcessing(mono_extraction_coeff, demodulatedSignal, mono_extraction_state, audio_taps, mono_extraction_block,startIndex_mono, audio_decim);
		}else{
      blockResample(mono_extraction_coeff, demodulatedSignal, mono_extraction_state, audio_taps, mono_extraction_block,audio_decim,audio_upSample);
    }

    std::fill (mono_block_filtered.begin(),mono_block_filtered.end(),0);

    blockConvolve(allPass_coeff, mono_extraction_block, allPass_State, audio_taps, mono_block_filtered);

    std::fill (pll_processed.begin(),pll_processed.end(),0);

    fmPLL(pll_processed, pll_block, freq, fs, ncoScale, phaseadjust, normBandwidth,pll_variables);

    std::fill (stereo_block.begin(),stereo_block.end(),0);

    std::fill (filtered_stereo.begin(),filtered_stereo.end(),0);

		for (int i = 0; i < pll_block.size(); i++){
	    stereo_block [i] = stereo_extraction_block[i] * pll_processed[i];
	  }

		blockConvolve(stereo_coeff, stereo_block, stereo_state, audio_taps, filtered_stereo);

    std::fill (stereoLeft.begin(),stereoLeft.end(),0);

    std::fill (stereoRight.begin(),stereoRight.end(),0);

    std::fill (stereo.begin(),stereo.end(),0);

    for (int i = 0; i < stereoLeft.size(); i++){
      stereoLeft [i] = (filtered_stereo[i] + mono_block_filtered[i])/2;
      stereoRight [i] = (mono_block_filtered[i] - filtered_stereo[i])/2;

  	//	stereo[2*i] = (stereoLeft[i] + stereoRight[i] - filtered_stereo[i])/2 + (stereoLeft[i] + stereoRight[i] - mono_block_filtered[i]);
      stereo[2*i] = stereoLeft[i];
      stereo[2*i+1] = stereoRight[i];
  	}


    //write to file and play in terminal
    std::vector<short int> final_data(stereo.size());
     for (unsigned int k=0;k<stereo.size();k++){
       if(std::isnan(stereo[k])) final_data[k]=0;
       else final_data[k] = static_cast<short int>(stereo[k]*16384);

     }

     fwrite(&final_data[0], sizeof(short int), final_data.size(), stdout);

		//blockCount += 1;
  }

}
//std::ref(rds_extraction_coeff),std::ref(rds_pll_coeff), std::ref(RDS_resample_coeff), std::ref(rds_allPass_coeff),std::ref(audio_taps), std::ref(rf_decim),std::ref(blockSize), std::ref(rds_demod_Q), std::ref(audio_Fs)
void rdsThread(std::vector<float> &rds_extraction_coeff, std::vector<float> &rds_pll_coeff, std::vector<float> &RDS_resample_coeff, std::vector<float> &rds_allPass_coeff, unsigned short int &audio_taps, unsigned short int &rf_decim, int &blockSize, std::queue<float> &rds_demod_Q, float &audio_Fs, unsigned short int &audio_decim, unsigned short int &audio_upSample,int &mode,std::vector<float> &cosFilt_coeff){

	if (mode == 0 || mode == 2){

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
	  int blockCount = 0;

		unsigned short int allPass_taps = 74;
		std::vector<float> rds_state(audio_taps - 1, 0);
		std::vector<float> state_pll(audio_taps - 1,0);
		std::vector<float> allPass_State(allPass_taps - 1,0);
		std::vector<float> rds_extraction_state(audio_taps - 1,0);
		std::vector<float> RDS_resample_state(audio_taps - 1,0);
		std::vector<float> cosFilt_state(audio_taps - 1,0);


		unsigned short int startIndex_rds = 0;

		//unsigned short int allPass_taps = 74;
		unsigned short int startIndex_mono = 0;

		std::vector<float> demodulatedSignal (blockSize/(2 * rf_decim),0);
		std::vector<float> rds_extraction_block((demodulatedSignal.size()), 0);
		std::vector<float> pll_block(rds_extraction_block.size(), 0);
		std::vector<float> pll_filtered (pll_block.size(), 0);
		std::vector<float> filtered_RDS ((pll_filtered.size()),0);
		std::vector<float> pll_processed (pll_block.size(), 0);
		std::vector<float> RDS_block ((pll_processed.size()),0);
		std::vector<float> resampled_RDS ((pll_processed.size()/audio_decim)*audio_upSample,0);
		std::vector<float> RDS_Cos_Filtered ((resampled_RDS.size()),0);


		float rds_upSample = 171;
		float rds_downSample = 640;
		if (mode == 2){
			rds_upSample = 19;
			rds_downSample = 48;
		}

		for (unsigned int blockCount =0; ;blockCount++){
	  //  std::cout <<"RDS block: " << blockCount << std::endl;

	/*
	    if (demod_Q.empty())
	      exit(1);*/

			//std::cout << "AA: " << std::endl;

	    for (int i = 0; i < demodulatedSignal.size(); i++){
			//	std::cout << "BB: " << i << std::endl;
	      demodulatedSignal[i] = rds_demod_Q.front();
	      rds_demod_Q.pop();

	    }
	    //std::cout << "BB: " << std::endl;

			std::fill (rds_extraction_block.begin(),rds_extraction_block.end(),0);
			std::fill (pll_block.begin(),pll_block.end(),0);

			blockConvolve(rds_extraction_coeff, demodulatedSignal, rds_extraction_state, audio_taps, rds_extraction_block);
			for (int i = 0; i < pll_block.size(); i++){
				pll_block[i] *= pll_block[i];
			}

			std::fill (pll_filtered.begin(),pll_filtered.end(),0);
			std::fill (filtered_RDS.begin(),filtered_RDS.end(),0);

	// ALL pass filter
	    blockConvolve(rds_allPass_coeff, rds_extraction_block, allPass_State, audio_taps, filtered_RDS);
	// Bandpass filter
	    blockConvolve(rds_pll_coeff, pll_block, state_pll, audio_taps, pll_filtered);

			std::fill (pll_processed.begin(),pll_processed.end(),0);

	    fmPLL(pll_processed,pll_filtered, freq, fs, ncoScale, phaseadjust, normBandwidth, pll_variables);

			std::fill (RDS_block.begin(),RDS_block.end(),0);

			for (int i = 0; i < pll_block.size(); i++){
		    RDS_block [i] = rds_extraction_block[i] * pll_processed[i];
		  }

			std::fill (resampled_RDS.begin(),resampled_RDS.end(),0);

			blockResample(RDS_resample_coeff, RDS_block, RDS_resample_state, audio_taps, resampled_RDS,rds_downSample,rds_upSample);


			std::fill (RDS_Cos_Filtered.begin(),RDS_Cos_Filtered.end(),0);
			blockConvolve(cosFilt_coeff, resampled_RDS, cosFilt_state, audio_taps, RDS_Cos_Filtered);



	  }
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
  float stereo_Fb = 23000.0;
  float stereo_Fe = 53000.0;

  float carrier_Fb = 18500.0;
  float carrier_Fe = 19500.0;

	float rds_Fb = 54000.0;
  float rds_Fe = 60000.0;

  //variables for PLL block process

  unsigned short int rf_upSample = 1;
  unsigned short int rf_taps = 51;
  unsigned short int rf_decim = 10;

  unsigned short int audio_upSample = 1;
  unsigned short int audio_taps = 51;
  unsigned short int audio_decim = 5;
  float audio_Fc = 16000.0;
  int blockSize = 1024 * rf_decim * audio_decim * 2;

	// RDS

	float rds_sampleFreq = 64125;

  if(mode == 0){
    rf_Fs = 2400000.0;
  //	rf_taps = 151;
    rf_decim = 10;

  	// audio variables
  	audio_Fs = 48000.0;
  //	audio_taps = 151;
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

		std::cerr << "RDS IS NOT SUPPORTED!"<< std::endl;

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
		rds_sampleFreq = 95000;

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
		std::cerr << "RDS IS NOT SUPPORTED!"<< std::endl;
  }
  int blockCount = 0;

  const std::string out_fname = "fmMonoBlock(cpp).wav";
    std::ofstream fdout(out_fname, std::ios::out | std::ios::binary);
//  std::vector <short int> play (stereo_data_final.size(), 0);

	// impulse response (reuse code from the previous experiment)
	std::vector<float> rf_coeff;
  impulseResponseLPF(rf_Fs, rf_Fc, rf_taps, rf_coeff,1);
  // Stereo Specific Variables

  std::vector<float> stereo_extraction_coeff;
	impulseResponseBPF(audio_Fs, stereo_Fb, stereo_Fe, audio_taps, stereo_extraction_coeff,audio_upSample);
  std::vector<float> carrier_coeff;
	impulseResponseBPF(audio_Fs, carrier_Fb, carrier_Fe, audio_taps, carrier_coeff,audio_upSample);
	std::vector<float> stereo_coeff;
	impulseResponseLPF(audio_Fs, 19000, audio_taps, stereo_coeff,audio_upSample);
  // End stereo specific

  // Start Mono Specific variables
  //std::vector<float> mono_extraction_state(audio_taps - 1,0);
  unsigned short int allPass_taps = 74;
  std::vector<float> mono_extraction_coeff;
  impulseResponseLPF(audio_Fs, audio_Fc, audio_taps, mono_extraction_coeff,audio_upSample);
  std::vector<float> allPass_coeff;
	impulseResponseLPF(audio_Fs,16500, allPass_taps, allPass_coeff,audio_upSample);

  std::vector<float> stereo_data_final;
  // End Mono Specific Variables

	// Start RDS Specific variables

	std::vector<float> rds_extraction_coeff;
	impulseResponseBPF(240000, rds_Fb, rds_Fe, audio_taps, rds_extraction_coeff,audio_upSample);

	std::vector<float> rds_allPass_coeff;
	impulseResponseLPF(240000,57000, allPass_taps, rds_allPass_coeff,audio_upSample);

	std::vector<float> rds_pll_coeff;
	impulseResponseBPF(240000,113500, 114500, audio_taps, rds_pll_coeff,audio_upSample);

	std::vector<float> RDS_resample_coeff;
	impulseResponseLPF(240000,57000, audio_taps, RDS_resample_coeff,audio_upSample);

	std::vector<float> cosFilt_coeff;
	impulseResponseRootRaisedCosine(rds_sampleFreq, audio_taps, cosFilt_coeff);

	// End RDS Specific Variables

  std::queue<float> demod_Q;
	std::queue<float> rds_demod_Q;

  std::thread rf_thread(rfThread,std::ref(rf_coeff),std::ref(rf_taps), std::ref(rf_decim),std::ref(blockSize), std::ref(demod_Q), std::ref(rds_demod_Q));
  std::this_thread::sleep_for(std::chrono::milliseconds(10000));

  std::thread audio_thread(audioThread,std::ref(stereo_extraction_coeff),std::ref(carrier_coeff), std::ref(mono_extraction_coeff), std::ref(allPass_coeff),std::ref(stereo_coeff),std::ref(audio_taps), std::ref(rf_decim), std::ref(audio_decim), std::ref(audio_upSample),std::ref(blockSize), std::ref(demod_Q),std::ref(stereo_data_final), std::ref(audio_Fs));

	std::thread rds_thread(rdsThread,std::ref(rds_extraction_coeff),std::ref(rds_pll_coeff), std::ref(RDS_resample_coeff), std::ref(rds_allPass_coeff),std::ref(audio_taps), std::ref(rf_decim),std::ref(blockSize), std::ref(rds_demod_Q), std::ref(audio_Fs),std::ref(audio_decim),std::ref(audio_upSample),std::ref(mode),std::ref(cosFilt_coeff));
//	void rdsThread(std::vector<float> &rds_extraction_coeff, std::vector<float> &rds_pll_coeff, std::vector<float> &RDS_resample_coeff, std::vector<float> &rds_allPass_coeff, unsigned short int &audio_taps, unsigned short int &rf_decim, int &blockSize, std::queue<float> &rds_demod_Q, std::vector<float> &stereo_data_final, float &audio_Fs)

  rf_thread.join();
//  std::cout <<"EEE "<< std::endl;
  // Pause for 1 second

  audio_thread.join();

	rds_thread.join();
//  std::cout <<"FFFF "<< std::endl;


  for (int i = 0; i < 150; i++){


  //  std::cout <<"stereo final: " << stereo_data_final[i] << std::endl;

  }
	//g++ (filename).cpp -o (file)
	//     ./test | aplay -c 2 -f S16_LE -r 48000

//	std::vector<float> monoSignal;
//	mono(monoSignal);


  fdout.close();

	return 0;
}
