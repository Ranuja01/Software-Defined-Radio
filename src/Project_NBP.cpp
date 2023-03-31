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
#include <thread>
#include <queue>
#include <chrono>




#include "PLL.cpp"
//#include "MonoBlock.cpp"


#define PI 3.14159265358979323846

void genIndexVector(std::vector<float> &x, const int size) {
  x.clear(); x.resize(size, static_cast<float>(0));
  for (int i=0; i<size; i++) {
    x[i] = static_cast<float>(i);
  }
}

// function to be used for logging a float vector in a .dat file (for .gnuplot)
// can be reused for different types of vectors with 32-bit floating point vals
void logVector(const std::string filename, \
  const std::vector<float> &x, \
  const std::vector<float> &y)
{
  // write data in text format to be parsed by gnuplot (change as needed)
  const std::string dat_filename = filename + ".dat";
  std::fstream fd;
  fd.open(dat_filename, std::ios::out);
  fd << "#\tx_axis\ty_axis\n";

  for (int i = 0; i < (int)x.size(); i++) {
    fd << "\t " << x[i] << "\t";
    // if the number of values on the Y axis is less than on the X tx_axis
    // then we just do not write anything on the Y axis
    if (i < (int)y.size())
      fd << y[i];
    fd << "\n";
  }
  std::cout << "Generated " << dat_filename << " to be used by gnuplot\n";
  fd.close();
}

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


void write_audio_data(std::vector<float> &audio, float audio_Fs, std::vector<short int> &play)
{

  //block process writing to file
  std::vector<short int>sample(audio.size());

  for(int k=0; k<audio.size(); k++){
    if(std::isnan(audio[k]))sample[k] = 0;
    else sample[k] = static_cast<short int>(audio[k]*16384);
  }
  //append each block
  play.insert(play.end(), sample.begin(), sample.end());
}

void blockConvolve(std::vector<float> &h, const std::vector<float> &block, std::vector<float> &state, int num_taps, std::vector<float> &filtered_block){

  //std::cout << filtered_block.size() << " " << h.size() << "  " << block.size() << std::endl;
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
	makeSubList(state,block,block.size() - num_taps + 1, block.size());
}

void blockProcessing(std::vector<float> &h, const std::vector<float> &block, std::vector<float> &state, int num_taps, std::vector<float> &filtered_block, unsigned short int &startIndex, int dRate){
  int endIndex = 0;
  for (int n = 0; n < block.size(); n+=dRate){
    for (int k = 0; k < h.size(); k++){
      if (n - k >= 0){
        if (n - k < block.size()){
          filtered_block[(n - startIndex)/dRate] += h[k] * block[n-k];
          if(k == 0 && n-k == 0){
          }
          if(k == 0 && n-k == 10){
          }
        }
      }else {
        if (n - k + num_taps - 1 < state.size()){
          filtered_block[(n - startIndex)/dRate] += h[k] * state[(n - k) + num_taps - 1];
        }
      }
    }
    endIndex = n;
  }
  startIndex = (endIndex - block.size()) + dRate;

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

	makeSubList(state,block,block.size() - num_taps + 1, block.size());
}


void makeSubList (std::vector<float> &subList, const std::vector<float> &list, int first, int last){
  subList.clear();
  for (int i = first; i < last; i++){


//			std::cout << "uytrew: " << list[i] << std::endl;


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
  /*
  std::cout <<"I "<< I[0] <<std::endl;
  std::cout <<"Q "<<Q[0] <<std::endl;
  std::cout <<"prevI "<<prevI <<std::endl;
  std::cout <<"prevQ "<<prevQ <<std::endl;
*/
  for (int k = 0; k < I.size(); k++){
    if (!( (I[k] * I[k]) + (Q[k] * Q[k]) == 0)){
      demodulatedSignal[k] = (1.0/( (I[k] * I[k]) + (Q[k] * Q[k]) ) ) * (I[k] * (Q[k] - prevQ) - Q[k] * (I[k] - prevI));


    }else {
      demodulatedSignal[k] = 0;
    }
    prevI = I[k];
    prevQ = Q[k];
  }

//	std::cout <<demodulatedSignal[0] <<std::endl;
}

//void write_stereo_data(std::vector<float> &audio, float audio_Fs, std::vector<short int> &play)
void write_stereo_data(std::vector<float> &audio, float audio_Fs)
{

  //block process writing to file
  std::vector<short int>sample(audio.size());

  for(int k=0; k<audio.size(); k++){
    if(std::isnan(audio[k]))sample[k] = 0;
    else sample[k] = static_cast<short int>(audio[k]*16384);
  }
  //append each block
  fwrite(&sample[0], sizeof(short int), sample.size(), stdout);
//  play.insert(play.end(), sample.begin(), sample.end());
}

//std::thread rf_thread(rfThread,std::ref(rf_taps), std::ref(rf_decim),std::ref(blockSize),std::ref(audio_data), std::ref(demod_Q));
void rfThread(std::vector<float> &rf_coeff, unsigned short int &rf_taps, unsigned short int &rf_decim, int &blockSize, std::vector<float> &audio_data,std::queue<float> &demod_Q){


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
  std::vector<float> demodulatedSignal ((int)(filtered_i.size()),0);

  while ((blockCount + 1) * blockSize < audio_data.size()){
    // Processing IQ samples to acquire them within 100khz range
    //std::cout <<"RF block: " << blockCount << std::endl;

    std::fill (filtered_i.begin(),filtered_i.end(),0);
    makeOddEvenSubList(block,audio_data,blockCount*blockSize,(blockCount + 1)*blockSize);
    blockProcessing(rf_coeff, block, state_i, rf_taps, filtered_i,startIndex_i,rf_decim);

    std::fill (filtered_q.begin(),filtered_q.end(),0);
    makeOddEvenSubList(block,audio_data,blockCount*blockSize + 1,(blockCount + 1)*blockSize);
    blockProcessing(rf_coeff, block, state_q, rf_taps, filtered_q,startIndex_q,rf_decim);

    // Demodulating to merge IQ samples
    std::fill (demodulatedSignal.begin(),demodulatedSignal.end(),0);
    fmDemod (demodulatedSignal, filtered_i, filtered_q,prevI,prevQ);


    for (int i = 0; i < demodulatedSignal.size(); i++){
    //  std::cout <<"dmSignal " << i << " " << demodulatedSignal[i] << std::endl;
      demod_Q.push(demodulatedSignal[i]);
    }

    blockCount++;
  }

}

//std::thread audio_thread(audioThread,std::ref(stereo_extraction_coeff),std::ref(carrier_coeff), std::ref(mono_extraction_coeff), std::ref(allPass_coeff),std::ref(stereo_coeff),std::ref(audio_taps), std::ref(rf_decim), std::ref(audio_decim), std::ref(audio_upSample),blockSize,std::ref(audio_data), std::ref(demod_Q),std::ref(stereo_data_final));
void audioThread(std::vector<float> &stereo_extraction_coeff, std::vector<float> &carrier_coeff, std::vector<float> &mono_extraction_coeff, std::vector<float> &allPass_coeff, std::vector<float> &stereo_coeff , unsigned short int &audio_taps, unsigned short int &rf_decim, unsigned short int &audio_decim, unsigned short int &audio_upSample, int &blockSize, std::vector<float> &audio_data, std::queue<float> &demod_Q, std::vector<float> &stereo_data_final, float &audio_Fs){
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

  while ((blockCount + 1) * blockSize < audio_data.size()){
    //std::cout <<"Audio block: " << blockCount << std::endl;

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

    stereo_data_final.insert(stereo_data_final.end(), stereo.begin(), stereo.end() );

    //write_stereo_data(stereo, audio_Fs/2);

    //write to file and play in terminal





		blockCount += 1;
  }

}





int main()
{
  int mode = 0;

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

  //const std::string out_fname = "fmMonoBlock(cpp).wav";
  //  std::ofstream fdout(out_fname, std::ios::out | std::ios::binary);
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
  //std::cout <<"CCCC "<< std::endl;
  std::thread rf_thread(rfThread,std::ref(rf_coeff),std::ref(rf_taps), std::ref(rf_decim),std::ref(blockSize),std::ref(audio_data), std::ref(demod_Q));
  std::this_thread::sleep_for(std::chrono::milliseconds(7000));
  //std::cout <<"DDD "<< std::endl;
  std::thread audio_thread(audioThread,std::ref(stereo_extraction_coeff),std::ref(carrier_coeff), std::ref(mono_extraction_coeff), std::ref(allPass_coeff),std::ref(stereo_coeff),std::ref(audio_taps), std::ref(rf_decim), std::ref(audio_decim), std::ref(audio_upSample),std::ref(blockSize),std::ref(audio_data), std::ref(demod_Q),std::ref(stereo_data_final), std::ref(audio_Fs));

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
std::vector <short int> play (stereo_data_final.size(), 0);
write_stereo_data(stereo_data_final, audio_Fs/2);

//write to file and play in terminal
const std::string out_fname = "fmMonoBlock(cpp).wav";
std::ofstream fdout(out_fname, std::ios::out | std::ios::binary);
fwrite(&play[0], sizeof(short int), play.size(), stdout);


  fdout.close();

	return 0;
}
