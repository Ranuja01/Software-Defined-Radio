#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>

#include "PLL.cpp"
#include "filter.cpp"
#include "iofunc.cpp"




void readStdinBlockData(unsigned int num_samples, unsigned int block_id, std::vector<float> &block_data){
	std::vector<char> raw_data(num_samples);
	std::cin.read(reinterpret_cast<char*>(&raw_data[0]), num_samples*sizeof(char));
	for (int k=0; k<(int)num_samples; k++){
		block_data[k] = float(((unsigned char)raw_data[k]-128)/128.0);
	}}


  int main(int argc, char* argv[])
  {

  	// Default mode 0
  	int mode = 0;
    int type=1; //1 means mono, 2 means stereo
  	// Mode Selection
  	if (argc<3){
  		std::cerr << "Operating in default mode 0 and mono" << std::endl;
  	} else if (argc==3){
  		mode=atoi(argv[1]);
      type = atoi(argv[2]);
  		if (mode>3){
  			std::cerr << "Wrong mode " << mode << std::endl;
  			exit(1);
  		}
  	} else {
  		std::cerr << "Usage: " << argv[0] << std::endl;
  		std::cerr << "or " << std::endl;
  		std::cerr << "Usage: " << argv[0] << " <mode>" << std::endl;
  		std::cerr << "\t\t <mode> is a value from 0-3" << argv[0] << std::endl;
  		exit(1);
  	}

  	std::cerr << "Operating in mode " << mode << std::endl;

//Defining variables
      float rf_Fs;
      float rf_Fc;
      unsigned short int rf_taps;
      unsigned short int rf_decim;
      float audio_Fs;// sample rate for our "assumed" audio (change as needed for 48 ksamples/sec audio files)
      float audio_Fc ;// cutoff frequency (explore ... but up-to Nyquist only!)
      unsigned short int rf_upSample;
      unsigned short int audio_upSample;
      unsigned short int audio_taps;
      unsigned short int audio_decim;
      int blockSize;
      int blockCount=0;


    //RF variable
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

      //Setting parameters based on mode of operation
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
       	 rf_taps = 101;
         rf_decim = 10;
       	// audio variables
       	 audio_Fs = 44100.0;	// sample rate for our "assumed" audio (change as needed for 48 ksamples/sec audio files)
       	 audio_Fc = 16000.0;	// cutoff frequency (explore ... but up-to Nyquist only!)
       	 audio_taps = 101;

         audio_upSample = 147;
         audio_decim = 800;
     		 blockSize = 1029 * rf_decim  * audio_decim/audio_upSample*2;
        //blockSize = audio_data.size();

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
         blockSize = 8 * rf_decim * audio_decim * 2;
       }

       //low pass filter impulse response
       std::vector<float> rf_coeff;
       impulseResponseLPF(rf_Fs, rf_Fc, rf_taps, rf_coeff);
       std::vector<float> audio_coeff;
       impulseResponseLPF(audio_Fs, audio_Fc, audio_taps, audio_coeff);
       std::vector<float> allPass_coeff;
     	 impulseResponseLPF(48000,100000, audio_taps, allPass_coeff);
       std::vector<float> carrier_coeff;
     	impulseResponseBPF(audio_Fs, carrier_Fb, carrier_Fe, audio_taps, carrier_coeff);
     	std::vector<float> stereo_coeff;
     	impulseResponseLPF(48000, 19000, audio_taps, stereo_coeff);

       std::vector<float> allPass_State(audio_taps - 1,0);
       std::vector<float> state_i(rf_taps - 1,0);
       std::vector<float> state_q(rf_taps - 1,0);
       std::vector<float> state_audio(audio_taps - 1,0);
       std::vector<float> filtered_block((blockSize/rf_decim), 0);

       float prevI = 0;
       float prevQ = 0;

       std::vector<float> mono_audio;
       std::vector<float> audio_data_final_stereo;

       std::vector<float> filtered_i((blockSize/(2 * rf_decim)), 0);
       std::vector<float> filtered_q((blockSize/(2 * rf_decim)), 0);
       std::vector<float> block;

       unsigned short int startIndex_i = 0;
       unsigned short int startIndex_q = 0;
       unsigned short int startIndex_audio = 0;
       unsigned short int startIndex_pll = 0;
       std::vector<float> state_pll(audio_taps - 1,0);
       std::vector<float> block_data(blockSize);

       std::vector<float> stereo_state(audio_taps - 1, 0);
       std::vector<float> stereo_data_final;

       std::vector<float> pll_block(filtered_i.size(), 0);
       std::vector<float> pll;


       for (unsigned int block_id =0; ;block_id++){


         //reading
         readStdinBlockData(blockSize,block_id,block_data);
         if (std::cin.rdstate()!=0){
           exit(1);
           }

         //RF front-end processing
         std::fill (filtered_i.begin(),filtered_i.end(),0);
         makeOddEvenSubList(block,block_data,0,blockSize);
         blockProcessing(rf_coeff, block, state_i, rf_taps, filtered_i,startIndex_i,rf_decim);
         std::fill (filtered_q.begin(),filtered_q.end(),0);
         makeOddEvenSubList(block,block_data,1,blockSize);
         blockProcessing(rf_coeff, block, state_q, rf_taps, filtered_q,startIndex_q,rf_decim);
         std::vector<float> demodulatedSignal ((int)(filtered_i.size()),0);
         fmDemod (demodulatedSignal, filtered_i, filtered_q,prevI,prevQ);


         //Mono processing //calculate mono for both cases
         std::vector<float> audio_block((demodulatedSignal.size()/audio_decim)*audio_upSample, 0);
         //check if resample or down sample is needed
         if (audio_upSample == 1){
           blockProcessing(audio_coeff, demodulatedSignal, state_audio, audio_taps, audio_block,startIndex_audio, audio_decim);
     			std::cout << "FFF: " << audio_block.size() << std::endl;
     		}else{
           blockResample(audio_coeff, demodulatedSignal, state_audio, audio_taps, audio_block,audio_decim,audio_upSample);
           //impulseResponseLPF2((rf_Fs/rf_decim)*audio_upSample, 16000, audio_upSample*151, audio_upSample, h);
     			//convolveFIR2(audio_block, interF, h, filtMonoz, audio_decim, audio_upSample);
     		}

        if (type ==2){

          std::vector<float> audio_block_filtered((demodulatedSignal.size()/audio_decim), 0);
          blockConvolve(allPass_coeff, audio_block, allPass_State, audio_taps, audio_block_filtered);
          mono_audio.insert(mono_audio.end(), audio_block_filtered.begin(), audio_block_filtered.end() ); //mono data after passed through an all pass filter


          std::vector<float> pll_block(audio_block.size(), 0);
      		//std::fill (filtered_block.begin(),filtered_block.end(),0);
      		// blockProcessing(rf_coeff, block, state_q, rf_taps, filtered_q,startIndex_q,rf_decim);
          if (audio_upSample == 1){
            blockProcessing(carrier_coeff, demodulatedSignal, state_pll, audio_taps, pll_block,startIndex_pll, audio_decim);
          }else{
            blockResample(carrier_coeff, demodulatedSignal, state_pll, audio_taps, pll_block, audio_decim,audio_upSample);
          }

          std::vector<float> pll_processed ((pll_block.size()),0);
          fmPLL(pll_processed, pll_block, freq, fs, ncoScale, phaseadjust, normBandwidth,pll_variables);

          std::vector<float> stereo_block ((pll_processed.size()),0);
          std::vector<float> filtered_stereo ((pll_processed.size()),0);

      		for (int i = 0; i < pll_block.size(); i++){
      	    stereo_block [i] = audio_block[i] * pll_processed[i];
      	  }

          audio_data_final_stereo.insert(audio_data_final_stereo.end(), audio_block.begin(), audio_block.end() );
          pll.insert(pll.end(), pll_block.begin(), pll_block.end());
          //std::cout<< "asdbdoooooooooo" << std::endl;
          blockConvolve(stereo_coeff, stereo_block, stereo_state, audio_taps, filtered_stereo);
        //  std::cout<< "1111bd" << std::endl;
          stereo_data_final.insert(stereo_data_final.end(), filtered_stereo.begin(), filtered_stereo.end() );

          //Extracting channels

          std::vector<float> stereoLeft ((stereo_data_final.size()),0);
        	std::vector<float> stereoRight ((stereo_data_final.size()),0);
        	std::vector<float> stereo ((stereo_data_final.size()*2),0);

          for (int i = 0; i < stereoLeft.size(); i++){

            stereoLeft [i] = (stereo_data_final[i] + mono_audio[i])/2;
            stereoRight [i] = (stereo_data_final[i] - mono_audio[i])/2;

        		stereo[2*i] = stereoLeft[i];
        		stereo[2*i+1] = stereoRight[i];
        	}

        	std::vector <short int> play (stereo.size(), 0);
        	write_stereo_data(stereo, audio_Fs/2, play);

        	//write to file and play in terminal
          const std::string out_fname = "fmMonoBlock(cpp).wav";
          std::ofstream fdout(out_fname, std::ios::out | std::ios::binary); //not sure how to write for stereo
          fwrite(&play[0], sizeof(short int), play.size(), stdout);
          fdout.close();



        }





       //write to file and play in terminal
       //const std::string out_fname = "fmMonoBlock(cpp).bin";
       else{
       std::vector<short int> audio_data(audio_block.size());
       for (unsigned int k=0;k<audio_block.size();k++){
         if(std::isnan(audio_block[k])) audio_data[k]=0;
         else audio_data[k] = static_cast<short int>(audio_block[k]*16384);

       }
       fwrite(&audio_data[0], sizeof(short int), audio_data.size(), stdout);}




       }

       return 0;
     }
