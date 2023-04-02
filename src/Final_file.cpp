/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Copyright by Nicola Nicolici
Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/
//Mono
// mode 0:works
// mode 1: works
// mode 2: works
// mode 3: works

//stereo
// mode 0:works
// mode 1: works
// mode 2: no works
// mode 3: no works

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

//rtl_sdr -f 99.9M -s 2.4M -| ./project 0 1 | aplay -c 1  -f S16_LE -r 48000

#define PI 3.14159265358979323846
//Reads Block Data
void readStdinBlockData(unsigned int num_samples, unsigned int block_id, std::vector<float> &block_data){
	std::vector<char> raw_data(num_samples); //vector of chars to store raw data read
	std::cin.read(reinterpret_cast<char*>(&raw_data[0]), num_samples*sizeof(char));
	for (int k=0; k<(int)num_samples; k++){ 
		block_data[k] = float(((unsigned char)raw_data[k]-128)/128.0); //converts raw data to float
	}
	}

//Main function for RF thread
void rfThread(std::vector<float> &rf_coeff, unsigned short int &rf_taps, unsigned short int &rf_decim, int &blockSize,std::queue<float> &demod_Q, std::queue<float> &rds_demod_Q){


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
   
    readStdinBlockData(blockSize,blockCount,audio_data); //reads a block of audio data into the vector
    if (std::cin.rdstate()!=0){ //if there was error while readin then break
      break;
    }

    std::fill (filtered_i.begin(),filtered_i.end(),0); //intialize to all 0
    makeOddEvenSubList(block,audio_data,0,blockSize); //extracts even samples from current block of audio data
    blockProcessing(rf_coeff, block, state_i, rf_taps, filtered_i,startIndex_i,rf_decim); //RF filtering to even samples stored in filtered_i - I samples
	  
    std::fill (filtered_q.begin(),filtered_q.end(),0); //intialize to all 0
    makeOddEvenSubList(block,audio_data,1,blockSize); //extracts odd samples from current block of audio data
    blockProcessing(rf_coeff, block, state_q, rf_taps, filtered_q,startIndex_q,rf_decim); //RF filtering to odd samples stored in filtered_q - Q samples
 
// Demodulate the FM signal using IQ compponents of filtered signal
    fmDemod (demodulatedSignal, filtered_i, filtered_q,prevI,prevQ);

    for (int i = 0; i < demodulatedSignal.size(); i++){ 
      demod_Q.push(demodulatedSignal[i]); //pushes samples of demod signal to demod_Q
      rds_demod_Q.push(demodulatedSignal[i]); //pushes samples of demod signal to rds_demod_Q
    }

  }

}

//Main Function for audio thread
void audioThread(std::vector<float> &stereo_extraction_coeff, std::vector<float> &carrier_coeff, std::vector<float> &mono_extraction_coeff, std::vector<float> &allPass_coeff, std::vector<float> &stereo_coeff , unsigned short int &audio_taps, unsigned short int &rf_decim, unsigned short int &audio_decim, unsigned short int &audio_upSample, int &blockSize, std::queue<float> &demod_Q, std::vector<float> &stereo_data_final, float &audio_Fs, bool &exitFlag, int type){

  
  std::vector<float> pll_variables(5,0);
	pll_variables[0] = 0;
	pll_variables[1] = 0;
	pll_variables[2] = 1;
	pll_variables[3] = 0;
	pll_variables[4] = 0;
  

	float freq = 19000;
	float fs = 44100;
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
  std::vector<float> mono_block_filtered((demodulatedSignal.size()/audio_decim)*audio_upSample, 0);
  std::vector<float> pll_processed ((pll_block.size()),0);
  std::vector<float> stereo_block ((pll_processed.size()),0);
  std::vector<float> filtered_stereo ((pll_processed.size()),0);
  std::vector<float> stereoLeft ((filtered_stereo.size()),0);
  std::vector<float> stereoRight ((filtered_stereo.size()),0);
  std::vector<float> stereo ((filtered_stereo.size()*2),0);
  std::vector<float>filtMonoz(audio_taps,0); // Added

// Loop for Audio Processing 
  for (unsigned int blockCount =0; ;blockCount++){

    if (demod_Q.empty()) //checks if there is data in queue
      exit(1);

    for (int i = 0; i < demodulatedSignal.size(); i++){

      demodulatedSignal[i] = demod_Q.front(); // extract samples from queue
      demod_Q.pop();
    }

    std::fill (stereo_extraction_block.begin(),stereo_extraction_block.end(),0); //intialize to all 0
    std::fill (pll_block.begin(),pll_block.end(),0); //intialize to all 0

 //if no upsampling is required
    if (audio_upSample == 1){
      blockProcessing(stereo_extraction_coeff, demodulatedSignal, stereo_extraction_state, audio_taps, stereo_extraction_block,startIndex_stereo, audio_decim);
      blockProcessing(carrier_coeff, demodulatedSignal, state_pll, audio_taps, pll_block,startIndex_pll, audio_decim);

    }else{
// start upsampling the demod signal and process filters
      blockResample(stereo_extraction_coeff, demodulatedSignal, stereo_extraction_state, audio_taps, stereo_extraction_block,audio_decim,audio_upSample);
      blockResample(carrier_coeff, demodulatedSignal, state_pll, audio_taps, pll_block, audio_decim,audio_upSample);
     
    }

   
    std::fill (mono_extraction_block.begin(),mono_extraction_block.end(),0); //intialize to all 0

    

    if (audio_upSample == 1){ //use block processing if value is 1 to extract mono signal
      blockProcessing(mono_extraction_coeff, demodulatedSignal, mono_extraction_state, audio_taps, mono_extraction_block,startIndex_mono, audio_decim);
    }else{ //use block resample otherwise to extract mono signal
      blockResample(mono_extraction_coeff, demodulatedSignal, mono_extraction_state, audio_taps, mono_extraction_block,audio_decim,audio_upSample);
    
    }

//Stereo type
if(type==2){
    std::fill (mono_block_filtered.begin(),mono_block_filtered.end(),0); //intialize to all 0
// Convolvution of mono_extraction_block with allPass_coeff stored into mono_block_filtered
    blockConvolve(allPass_coeff, mono_extraction_block, allPass_State, audio_taps, mono_block_filtered); 

    std::fill (pll_processed.begin(),pll_processed.end(),0); //intialize to all 0
// Apply FM PLL to pll_block and store in pll_processed
    fmPLL(pll_processed, pll_block, freq, fs, ncoScale, phaseadjust, normBandwidth,pll_variables);

    std::fill (stereo_block.begin(),stereo_block.end(),0); //intialize to all 0
    std::fill (filtered_stereo.begin(),filtered_stereo.end(),0); //intialize to all 0

    for (int i = 0; i < pll_block.size(); i++){
	    stereo_block [i] = stereo_extraction_block[i] * pll_processed[i]; //multiply stereo_extraction_block and pll_processed to form stereo_block
	  }
// Convolution of stereo_block with stereo_coefficients to get filtered stereo signal
    blockConvolve(stereo_coeff, stereo_block, stereo_state, audio_taps, filtered_stereo);

    std::fill (stereoLeft.begin(),stereoLeft.end(),0); //intialize to all 0
    std::fill (stereoRight.begin(),stereoRight.end(),0); //intialize to all 0
    std::fill (stereo.begin(),stereo.end(),0); //intialize to all 0
	
// Combine mono_block_filtered and filtered_stereo signals to form stereo signal
    for (int i = 0; i < stereoLeft.size(); i++){
      stereoLeft [i] = (filtered_stereo[i] + mono_block_filtered[i])/2; //average of samples is left channel value
      stereoRight [i] = (mono_block_filtered[i] - filtered_stereo[i])/2;//difference of samples is right channel value
      stereo[2*i] = stereoLeft[i]; // left channel samples stored in stereo buffer
      stereo[2*i+1] = stereoRight[i];} //right channel samples stored in stereo buffer
	
    std::vector<short int> final_data(stereo.size()); //vector to store the final stereo audio data
    for (unsigned int k=0;k<stereo.size();k++){ // loops thru each element in stereo vector and convert float value to int value
       if(std::isnan(stereo[k])) final_data[k]=0; // if sample is not a number set to 0
       else final_data[k] = static_cast<short int>(stereo[k]*16384); // if sample is not not a number, convert to short int
     }

     fwrite(&final_data[0], sizeof(short int), final_data.size(), stdout);} // output final processed audio data to standard output
     
     //if audio input is mono
    else{
		std::vector<short int> final_data(mono_extraction_block.size()); //vector to store the final mono audio data
	for (unsigned int k=0;k<mono_extraction_block.size();k++){// iterates thru mono extraction block and processes each sample
       if(std::isnan(mono_extraction_block[k])) final_data[k]=0; // if sample is not a number set to 0
       else final_data[k] = static_cast<short int>(mono_extraction_block[k]*16384); // if sample is not not a number, convert to short int

     }

     fwrite(&final_data[0], sizeof(short int), final_data.size(), stdout);}  // output final processed audio data to standard output

  }

}
//RDS thread for RDS data
void rdsThread(std::vector<float> &rds_extraction_coeff, std::vector<float> &rds_pll_coeff, std::vector<float> &RDS_resample_coeff, std::vector<float> &rds_allPass_coeff, unsigned short int &audio_taps, unsigned short int &rf_decim, int &blockSize, std::queue<float> &rds_demod_Q, float &audio_Fs, unsigned short int &audio_decim, unsigned short int &audio_upSample,int &mode,std::vector<float> &cosFilt_coeff){

	if (mode == 0 || mode == 2){ // checks if mode is 0 or 2

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


		float rds_upSample = 171; // set the upsample and downsample ratio for RDS signal based on mode
		float rds_downSample = 640;
		if (mode == 2){
			rds_upSample = 19;
			rds_downSample = 48;
		}

		for (unsigned int blockCount =0; ;blockCount++){ //processing demodulated signal in blocks
	  

	    for (int i = 0; i < demodulatedSignal.size(); i++){ 
	    
	      demodulatedSignal[i] = rds_demod_Q.front(); //reads demodulated signal from queue
	      rds_demod_Q.pop();

	    }
	
			std::fill (rds_extraction_block.begin(),rds_extraction_block.end(),0);//intialize to all 0
			std::fill (pll_block.begin(),pll_block.end(),0);//intialize to all 0
			
// Convolution of the RDS extraction filter coefficients with demodulated RDS signal
			blockConvolve(rds_extraction_coeff, demodulatedSignal, rds_extraction_state, audio_taps, rds_extraction_block);
			for (int i = 0; i < pll_block.size(); i++){
				pll_block[i] *= pll_block[i]; // square the values in PLL block
			}

			std::fill (pll_filtered.begin(),pll_filtered.end(),0);//intialize to all 0
			std::fill (filtered_RDS.begin(),filtered_RDS.end(),0);//intialize to all 0

	// ALL pass filter
	    blockConvolve(rds_allPass_coeff, rds_extraction_block, allPass_State, audio_taps, filtered_RDS); //convolution of all- pass filter coeff with fitered RDS signal
	// Bandpass filter
	    blockConvolve(rds_pll_coeff, pll_block, state_pll, audio_taps, pll_filtered); // convolution of PLL filter coeffs with PLL block

			std::fill (pll_processed.begin(),pll_processed.end(),0);//intialize to all 0
			
// Apply phase-locked loop (PLL) to PLL filtered signal
	    fmPLL(pll_processed,pll_filtered, freq, fs, ncoScale, phaseadjust, normBandwidth, pll_variables);

			std::fill (RDS_block.begin(),RDS_block.end(),0);//intialize to all 0

			for (int i = 0; i < pll_block.size(); i++){
		    RDS_block [i] = rds_extraction_block[i] * pll_processed[i]; // Multiply RDS extraction block with PLL processed block to obtain RDS block

		  }

			std::fill (resampled_RDS.begin(),resampled_RDS.end(),0);//intialize to all 0
			
// Resampling RDS block with resample filter coefficients
			blockResample(RDS_resample_coeff, RDS_block, RDS_resample_state, audio_taps, resampled_RDS,rds_downSample,rds_upSample);
			
			std::fill (RDS_Cos_Filtered.begin(),RDS_Cos_Filtered.end(),0);//intialize to all 0

			blockConvolve(cosFilt_coeff, resampled_RDS, cosFilt_state, audio_taps, RDS_Cos_Filtered);// convolution of cosine filter coeffs with resampled RDS signal

	  }
	}
}

// Main Function
int main(int argc, char* argv[])
{
  int mode; //can be 0,1,2,3
  int type; //can be 1:mono or 2:stereo
  
  // Mode Selection
  if (argc<3){ //Default mode and type
    std::cerr << "Operating in default mode 0 and mono" << std::endl;
    mode = 0;
    type = 1;
    
  } else if (argc==3){
    mode=atoi(argv[1]);
    type=atoi(argv[2]);
    if (mode>3){ // Invalid mode
      std::cerr << "Wrong mode " << mode << std::endl;
      exit(1);}
  } else { // Prints usage help and exit
    std::cerr << "Usage: " << argv[0] << std::endl;
    std::cerr << "or " << std::endl;
    std::cerr << "Usage: " << argv[0] << " <mode>" << std::endl;
    std::cerr << "\t\t <mode> is a value from 0-3" << argv[0] << std::endl;
    exit(1);}

  std::cerr << "Operating in mode " << mode << std::endl; // prints selected mode
  std::cerr << "Operating in Type " << mode << std::endl; // prints selected type
  
  
  // RF variables
  float rf_Fs = 2400000.0;
  float rf_Fc = 100000.0;
  float audio_Fs = 48000.0;	// sample rate for our "assumed" audio (change as needed for 48 ksamples/sec audio files)
  float audio_Fb = 23000.0;
  float audio_Fe = 53000.0;
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
  float rds_sampleFreq = 64125;
	
// Selects RF and audio parameters based on the mode selected
	// MODE 0
  if(mode == 0){ 
    rf_Fs = 2400000.0;
    rf_decim = 10;
  	audio_Fs = 48000.0;
    audio_decim = 5;
	// MODE 1
  }else if (mode == 1){
    rf_Fs = 1440000.0;
  	rf_taps = 51;
    rf_decim = 5;
    blockSize = 1024 * rf_decim * audio_decim * 2;
  	audio_Fs = 48000.0;
  	audio_taps = 51;
    audio_decim = 6;
 	// MODE 2
  }else if (mode == 2){
    rf_Fs = 2400000.0;
    rf_taps = 51; //changed
    rf_decim = 10;
    blockSize = 1024 * rf_decim * audio_decim * 2; //changed
  	// audio variables
  	audio_Fs = 44100.0;
  	audio_taps = 51;
	audio_upSample = 147;
	audio_decim = 800;
	rds_sampleFreq = 95000;
 	// MODE 3
  }else if (mode == 3){
    rf_Fs = 1920000.0;
    rf_taps = 51;
    rf_decim = 5;
    blockSize = 1024 * rf_decim * audio_decim * 2;
    // audio variables
    audio_Fs = 44100.0;
    audio_taps = 51;
    audio_upSample = 147;
    audio_decim = 1280;
  }
  int blockCount = 0;


// generates low pass filter coefficients for RF filtering
	std::vector<float> rf_coeff; 
	impulseResponseLPF(rf_Fs, rf_Fc, rf_taps, rf_coeff,1);
// generates band pass filter coefficients for stereo extraction
	std::vector<float> stereo_extraction_coeff;
	impulseResponseBPF(audio_Fs*audio_upSample, audio_Fb, audio_Fe, audio_taps*audio_upSample, stereo_extraction_coeff,1);
// generates band pass filter coefficients for carrier extraction
	std::vector<float> carrier_coeff;
	impulseResponseBPF(audio_Fs*audio_upSample, carrier_Fb, carrier_Fe, audio_taps*audio_upSample, carrier_coeff,1);
// generates low pass filter coefficients for stereo filtering
	std::vector<float> stereo_coeff;
	impulseResponseLPF(48000*audio_upSample, 19000, audio_taps*audio_upSample, stereo_coeff,audio_upSample);	
// generates low pass filter coefficients for mono extraction	
	std::vector<float> mono_extraction_coeff;
	impulseResponseLPF(audio_Fs*audio_upSample, audio_Fc, audio_taps*audio_upSample, mono_extraction_coeff,audio_upSample);

// all pass filter for stereo audio	
	unsigned short int allPass_taps = 74;
	std::vector<float> allPass_coeff;
	impulseResponseLPF(48000*audio_upSample,16500, allPass_taps*audio_upSample, allPass_coeff, audio_upSample);

	// Start Mono Specific variables
	std::vector<float> mono_extraction_state(audio_taps - 1,0);
	
	// Start RDS Specific variables
	
// band pass filter for extracting RDS signal
	std::vector<float> rds_extraction_coeff;
	impulseResponseBPF(240000, rds_Fb, rds_Fe, audio_taps, rds_extraction_coeff,audio_upSample);
// all pass filter for RDS signal
	std::vector<float> rds_allPass_coeff;
	impulseResponseLPF(240000,57000, allPass_taps, rds_allPass_coeff,audio_upSample);
// band pass filter for RDS signal PLL
	std::vector<float> rds_pll_coeff;
	impulseResponseBPF(240000,113500, 114500, audio_taps, rds_pll_coeff,audio_upSample);
// low pass filter for resampling RDS signal
	std::vector<float> RDS_resample_coeff;
	impulseResponseLPF(240000,57000, audio_taps, RDS_resample_coeff,audio_upSample);
// root raised cosine filter for RDS signal
	std::vector<float> cosFilt_coeff;
	impulseResponseRootRaisedCosine(rds_sampleFreq, audio_taps, cosFilt_coeff);

	// End RDS Specific Variables
	unsigned short int startIndex_mono = 0;
	std::vector<float> stereo_data_final; /
    	std::queue<float> demod_Q;
	std::queue<float> rds_demod_Q;
	
	bool exitFlag = false;
 
	std::thread rf_thread(rfThread,std::ref(rf_coeff),std::ref(rf_taps), std::ref(rf_decim),std::ref(blockSize), std::ref(demod_Q), std::ref(rds_demod_Q));
	std::this_thread::sleep_for(std::chrono::milliseconds(5000));
	
	std::thread audio_thread(audioThread,std::ref(stereo_extraction_coeff),std::ref(carrier_coeff), std::ref(mono_extraction_coeff), std::ref(allPass_coeff),std::ref(stereo_coeff),std::ref(audio_taps), std::ref(rf_decim), std::ref(audio_decim), std::ref(audio_upSample),std::ref(blockSize), std::ref(demod_Q),std::ref(stereo_data_final), std::ref(audio_Fs), std::ref(exitFlag),std::ref(type));
	//std::thread rds_thread(rdsThread,std::ref(rds_extraction_coeff),std::ref(rds_pll_coeff), std::ref(RDS_resample_coeff), std::ref(rds_allPass_coeff),std::ref(audio_taps), std::ref(rf_decim),std::ref(blockSize), std::ref(rds_demod_Q), std::ref(audio_Fs),std::ref(audio_decim),std::ref(audio_upSample),std::ref(mode),std::ref(cosFilt_coeff));
	
	rf_thread.join();
	audio_thread.join();
	//rds_thread.join();


	return 0;
}
