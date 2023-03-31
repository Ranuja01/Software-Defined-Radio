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

//std::thread rf_thread(rfThread,std::ref(rf_taps), std::ref(rf_decim),std::ref(blockSize),std::ref(audio_data), std::ref(demod_Q));
void rfThread(std::vector<float> &rf_coeff, unsigned short int &rf_taps, unsigned short int &rf_decim, int &blockSize,std::queue<float> &demod_Q, bool &exitFlag){



//  std::cout <<"AAA "<< std::endl;
  int blockCount = 0;
  std::vector<float> state_i(rf_taps - 1,0);
  std::vector<float> state_q(rf_taps - 1,0);


  float prevI = 0;
  float prevQ = 0;

  std::vector<float> filtered_i((blockSize/(2 * rf_decim)), 0);
	std::vector<float> filtered_q((blockSize/(2 * rf_decim)), 0);
  std::vector<float> block;
  //std::cout <<"BBB "<< std::endl;
  unsigned short int startIndex_i = 0;
	unsigned short int startIndex_q = 0;
  std::vector<float> audio_data(blockSize);
  std::vector<float> demodulatedSignal ((int)(filtered_i.size()),0);

  for (unsigned int blockCount =0; ;blockCount++){
    // Processing IQ samples to acquire them within 100khz range
    //std::cout <<"RF block: " << blockCount << std::endl;

    readStdinBlockData(blockSize,blockCount,audio_data);
    if (std::cin.rdstate()!=0){
      //exitFlag == true;
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
    }

    //blockCount++;
  }

}

//std::thread audio_thread(audioThread,std::ref(stereo_extraction_coeff),std::ref(carrier_coeff), std::ref(mono_extraction_coeff), std::ref(allPass_coeff),std::ref(stereo_coeff),std::ref(audio_taps), std::ref(rf_decim), std::ref(audio_decim), std::ref(audio_upSample),blockSize,std::ref(audio_data), std::ref(demod_Q),std::ref(stereo_data_final));
void audioThread(std::vector<float> &stereo_extraction_coeff, std::vector<float> &carrier_coeff, std::vector<float> &mono_extraction_coeff, std::vector<float> &allPass_coeff, std::vector<float> &stereo_coeff , unsigned short int &audio_taps, unsigned short int &rf_decim, unsigned short int &audio_decim, unsigned short int &audio_upSample, int &blockSize, std::queue<float> &demod_Q, std::vector<float> &stereo_data_final, float &audio_Fs, bool &exitFlag){
//std::cout <<"zzzzz "<< std::endl;
  std::vector<float> pll_variables(5,0);
//  std::cout <<"CLKIJHRWCCC "<< std::endl;
	pll_variables[0] = 0;
	pll_variables[1] = 0;
	pll_variables[2] = 1;
	pll_variables[3] = 0;
	pll_variables[4] = 0;
  //std::cout <<"DDDD "<< std::endl;

	float freq = 19000;
	float fs = 48000;
	float phaseadjust = 0.0;

//  std::cout <<"aaaa "<< std::endl;
	float normBandwidth = 0.01;
  float trigOffset = 0.0;
  float ncoScale = 2.0;
  int blockCount = 0;

  std::vector<float> stereo_state(audio_taps - 1, 0);
//  std::cout <<"GGGG "<< std::endl;
  std::vector<float> state_pll(audio_taps - 1,0);
//  std::cout <<"HHHH "<< std::endl;

	unsigned short int startIndex_stereo = 0;
  unsigned short int startIndex_pll = 0;

  std::vector<float> stereo_extraction_state(audio_taps - 1,0);
//  std::cout <<"IIII "<< std::endl;

  std::vector<float> mono_extraction_state(audio_taps - 1,0);

  unsigned short int allPass_taps = 74;
  unsigned short int startIndex_mono = 0;
//  std::cout <<"JJJJ "<< std::endl;

  std::vector<float> allPass_State(allPass_taps - 1,0);
//  std::cout <<"KKKK "<< std::endl;


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

//  std::cout <<"ppppppppppppppppppppppp "<< demodulatedSignal.size() <<  std::endl;
//  std::cout <<"LLLL "<< std::endl;

  for (unsigned int blockCount =0; ;blockCount++){
    //std::cout <<"Audio block: " << blockCount << std::endl;

    if (demod_Q.empty())
      exit(1);

    for (int i = 0; i < demodulatedSignal.size(); i++){

      demodulatedSignal[i] = demod_Q.front();
      //std::cout <<"dddd123 "<< demodulatedSignal[i] <<  std::endl;
      demod_Q.pop();
    //  std::cout << i<< std::endl;

    }
    //std::cout <<"gfd "<< std::endl;

    std::fill (stereo_extraction_block.begin(),stereo_extraction_block.end(),0);
    //std::cout <<"NNNN "<< std::endl;


    std::fill (pll_block.begin(),pll_block.end(),0);
//    std::cout <<"OOOO "<< stereo_extraction_block.size() <<  std::endl;
//    std::cout <<"OOOO1 "<< pll_block.size() <<  std::endl;
//    std::cout <<"OOOO2 "<< stereo_extraction_state.size() <<  std::endl;
//    std::cout <<"OOOO3 "<< demodulatedSignal.size() <<  std::endl;
//    std::cout <<"OOOO4 "<<blockSize/(2 * rf_decim) <<  std::endl;


    if (audio_upSample == 1){
      blockProcessing(stereo_extraction_coeff, demodulatedSignal, stereo_extraction_state, audio_taps, stereo_extraction_block,startIndex_stereo, audio_decim);
//      std::cout <<"POIJHNJH "<< std::endl;
      blockProcessing(carrier_coeff, demodulatedSignal, state_pll, audio_taps, pll_block,startIndex_pll, audio_decim);
//      std::cout <<"PPPP "<< std::endl;
    }else{
      blockResample(stereo_extraction_coeff, demodulatedSignal, stereo_extraction_state, audio_taps, stereo_extraction_block,audio_decim,audio_upSample);
      blockResample(carrier_coeff, demodulatedSignal, state_pll, audio_taps, pll_block, audio_decim,audio_upSample);
//      std::cout <<"QQQQ "<< std::endl;
    }


    std::fill (mono_extraction_block.begin(),mono_extraction_block.end(),0);

    if (audio_upSample == 1){
      blockProcessing(mono_extraction_coeff, demodulatedSignal, mono_extraction_state, audio_taps, mono_extraction_block,startIndex_mono, audio_decim);
//      std::cout <<"RRRR "<< std::endl;
		}else{
      blockResample(mono_extraction_coeff, demodulatedSignal, mono_extraction_state, audio_taps, mono_extraction_block,audio_decim,audio_upSample);
//      std::cout <<"SSSS "<< std::endl;
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

  //  std::cout <<"sdfghj block: " << blockCount << std::endl;
    for (int i = 0; i < stereoLeft.size(); i++){
      //std::cout <<"lll block: " << blockCount << std::endl;
      stereoLeft [i] = (filtered_stereo[i] + mono_block_filtered[i])/2;
      stereoRight [i] = (mono_block_filtered[i] - filtered_stereo[i])/2;

  	//	stereo[2*i] = (stereoLeft[i] + stereoRight[i] - filtered_stereo[i])/2 + (stereoLeft[i] + stereoRight[i] - mono_block_filtered[i]);
      stereo[2*i] = stereoLeft[i];
      stereo[2*i+1] = stereoRight[i];
  	}

    for (int i = 0; i < 15; i++){


    //  std::cout <<"stereo: " << pll_processed[i] << std::endl;

    }

    //stereo_data_final.insert(stereo_data_final.end(), stereo.begin(), stereo.end() );

    //write_stereo_data(stereo, audio_Fs/2);

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

  unsigned short int rf_upSample = 1;
  unsigned short int rf_taps = 51;
  unsigned short int rf_decim = 10;

  unsigned short int audio_upSample = 1;
  unsigned short int audio_taps = 51;
  unsigned short int audio_decim = 5;
  float audio_Fc = 16000.0;
  int blockSize = 1024 * rf_decim * audio_decim * 2;

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

  const std::string out_fname = "fmMonoBlock(cpp).wav";
    std::ofstream fdout(out_fname, std::ios::out | std::ios::binary);
//  std::vector <short int> play (stereo_data_final.size(), 0);




	// impulse response (reuse code from the previous experiment)
	std::vector<float> rf_coeff;
  impulseResponseLPF(rf_Fs, rf_Fc, rf_taps, rf_coeff);
  // Stereo Specific Variables


  std::vector<float> stereo_extraction_coeff;
	impulseResponseBPF(audio_Fs, audio_Fb, audio_Fe, audio_taps, stereo_extraction_coeff);

  std::vector<float> carrier_coeff;
	impulseResponseBPF(audio_Fs, carrier_Fb, carrier_Fe, audio_taps, carrier_coeff);

	std::vector<float> stereo_coeff;
	impulseResponseLPF(48000, 19000, audio_taps, stereo_coeff);

  // End stereo specific

  // Start Mono Specific variables

  std::vector<float> mono_extraction_state(audio_taps - 1,0);
  unsigned short int allPass_taps = 74;
  std::vector<float> mono_extraction_coeff;
  impulseResponseLPF(audio_Fs, audio_Fc, audio_taps, mono_extraction_coeff);

  std::vector<float> allPass_coeff;
	impulseResponseLPF(48000,16500, allPass_taps, allPass_coeff);

  unsigned short int startIndex_mono = 0;

//  std::vector<float> allPass_State(allPass_taps - 1,0);
  std::vector<float> stereo_data_final;
  // End Mono Specific Variables

  std::queue<float> demod_Q;
  bool exitFlag = false;
  //std::cout <<"CCCC "<< std::endl;
  std::thread rf_thread(rfThread,std::ref(rf_coeff),std::ref(rf_taps), std::ref(rf_decim),std::ref(blockSize), std::ref(demod_Q), std::ref(exitFlag));
  std::this_thread::sleep_for(std::chrono::milliseconds(10000));
  //std::cout <<"DDD "<< std::endl;
  std::thread audio_thread(audioThread,std::ref(stereo_extraction_coeff),std::ref(carrier_coeff), std::ref(mono_extraction_coeff), std::ref(allPass_coeff),std::ref(stereo_coeff),std::ref(audio_taps), std::ref(rf_decim), std::ref(audio_decim), std::ref(audio_upSample),std::ref(blockSize), std::ref(demod_Q),std::ref(stereo_data_final), std::ref(audio_Fs), std::ref(exitFlag));

  rf_thread.join();
//  std::cout <<"EEE "<< std::endl;
  // Pause for 1 second

  audio_thread.join();
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
