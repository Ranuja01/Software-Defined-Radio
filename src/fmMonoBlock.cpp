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
void read_audio_data(const std::string in_fname, std::vector<float> &audio_data)
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
						num_samples*sizeof(float));

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
void write_audio_data(const std::string out_fname, const std::vector<float> &audio_left, const std::vector<float> &audio_right)
{
	// file descriptor for the output to be written
	if (audio_left.size() != audio_right.size()) {
		std::cout << "Something got messed up with audio channels\n";
		std::cout << "They must have the same size ... exiting\n";
		exit(1);
	} else {
		std::cout << "Writing raw audio to \"" << out_fname << "\"\n";
	}
	std::ofstream fdout(out_fname, std::ios::binary);
	for (int i=0; i<(int)audio_left.size(); i++) {
		// we assume we have handled a stereo audio file
		// hence, we must interleave the two channels
		// (change as needed if testing with mono files)
		fdout.write(reinterpret_cast<const char*>(&audio_left[i]),\
								sizeof(audio_left[i]));
		fdout.write(reinterpret_cast<const char*>(&audio_right[i]),\
								sizeof(audio_right[i]));
	}
	fdout.close();
}

void blockProcessing(std::vector<float> &h, const std::vector<float> &block, std::vector<float> &state, int num_taps, std::vector<float> &filtered_block, int &startIndex, int dRate){
	int endIndex = 0;
	for (int n = startIndex; n < block.size(); n++){
		for (int k = 0; k < h.size(); k++){
			if (n - k >= 0){
				if (n - k < block.size()){
					filtered_block[(n - startIndex)/dRate] += h[k] * block[n-k];
				}
			}else {
		  	if (n - k + num_taps < state.size()){
					filtered_block[(n - startIndex)/dRate] += h[k] * state[(n - k) + num_taps];
				}
			}
		}
		endIndex = n;
	}
	startIndex = 1 + (endIndex + block.size() % dRate) - blockSize;
	makeSubList(state,block,block.size() - num_taps + 1, block.size());
}

void makeSubList (std::vector<float> &subList, const std::vector<float> &list, int first, int last){
	subList.clear();
	for (int i = first; i < last; i++){
		subList.push_back(list[i]);
	}

}

void makeOddEvenSubList (std::vector<float> &subList, const std::vector<float> &list, int first, int last){
	subList.clear();
	for (int i = first; i < last; i+=2){
		subList.push_back(list[i]);
	}

}

void downSample (){

}


int main()
{
	// assume the wavio.py script was run beforehand to produce a binary file
	const std::string in_fname = "../data/float32samples.bin";
	// declare vector where the audio data will be stored
	std::vector<float> iq_data;
  std::vector<float> audio_data;
	// note: we allocate memory for audio_data from within this read function
	read_audio_data(in_fname, iq_data);

  for (int i = 0; i < iq_data.size(); i++){
    audio_data[i] = ((float)audio_data[i] - 128.0)/128.0;
  }

  // RF variables
  float rf_Fs = 240000.0;
	float rf_Fc = 100000.0;
	unsigned short int rf_taps = 151;
  unsigned short int rf_decim = 10;

	// audio variables
	float Fs = 48000.0;	// sample rate for our "assumed" audio (change as needed for 48 ksamples/sec audio files)
	float Fc = 16000.0;	// cutoff frequency (explore ... but up-to Nyquist only!)
	unsigned short int audio_taps = 151;
  unsigned short int audio_decim = 10;

	int blockSize = 1024 * rf_decim * audio_decim * 2;
  int blockCount = 0;
  int input = 2;

	// impulse response (reuse code from the previous experiment)
	std::vector<float> rf_coeff;
	impulseResponseLPF(rf_Fs, rf_Fc, rf_taps, rf_coeff);

  std::vector<float> audio_coeff;
	impulseResponseLPF(audio_Fs, audio_Fc, audio_taps, audio_coeff);
// MAYBE WE NEED HANN WINDOW

  std::vector<float> state_i(rf_taps - 1,0);
  std::vector<float> state_q(rf_taps - 1,0);
  std::vector<float> state_audio(audio_taps - 1,0);
  std::vector<float> filtered_block(blockSize, 0);

  float prevI = 0;
  float prevQ = 0;

  std::vector<float> filtered_i, filtered_q;
  std::vector<float> block;

	unsigned short int startIndex_i = 0;
	unsigned short int startIndex_q = 0;

  while (blockCount + 1) * blockSize < iq_data.size()){
    std::cout <<"Processing block: " << blockCount << std::endl;

    std::fill (filtered_block.begin(),filtered_block.end(),0);
    makeOddEvenSubList(block,iq_data,blockCount*blockSize,(blockCount + 1)*blockSize);
    blockProcessing(rf_coeff, block, state_i, rf_taps, filtered_block,startIndex_i,rf_decim);
    filtered_i.insert( filtered_i.end(), filtered_block.begin(), filtered_block.end() );

    std::vector<float> filtered_block(blockSize, 0);
    makeOddEvenSubList(block,iq_data,blockCount*blockSize + 1,(blockCount + 1)*blockSize);
    blockProcessing(rf_coeff, block, state_q, rf_taps, filtered_block,startIndex_q,rf_decim);
    filtered_q.insert(filtered_q.end(), filtered_block.begin(), filtered_block.end() );







  }




	// declare vectors where the audio left/right channels will be stored
	std::vector<float> audio_left, audio_right;
	// note: we allocate the memory for the left/right channels
	// from within the split function that is called in the code below
	split_audio_into_channels(audio_data, audio_left, audio_right);

	// convolution code for filtering (reuse from the previous experiment)
	std::vector<float> single_pass_left, single_pass_right;

	if (input == 1){
		std::cout << "BBBB";
		convolveFIR(single_pass_left, audio_left, h);
		convolveFIR(single_pass_right, audio_right, h);
		// note: by default the above convolution produces zero on the output stream
		// YOU will need to update the convolveFIR and impulseResponseLPF functions

		// create a binary file to be read by wavio.py script to produce a .wav file
		// note: small adjustments will need to be made to wavio.py, i.e., you should
		// match the filenames, no need for self-checks as default Python code, ...
		const std::string out_fname = "../data/float32filtered_single.bin";
		write_audio_data(out_fname, single_pass_left,	single_pass_right);
  } else {

		std::cout << "AAAAA";

		std::vector<float> filtered_left, filtered_right;

		audio_left.resize(audio_left.size() + blockSize - audio_left.size() % blockSize,0);
		audio_right.resize(audio_right.size() + blockSize - audio_right.size() % blockSize,0);

		std::vector<float> block;
    std::vector<float> block_right;

		std::vector<float> state_left(h.size() - 1, 0);
		std::vector<float> state_right(h.size() - 1, 0);
		std::vector<float> filtered_block(blockSize, 0);
		int blockNum = (int)(audio_left.size() / blockSize);

		for (int i = 0; i < blockNum; i++){

			std::fill (filtered_block.begin(),filtered_block.end(),0);
			makeSubList(block,audio_left,i*blockSize,(i + 1)*blockSize);
			blockProcessing(h, block, state_left, num_taps, filtered_block);
		  filtered_left.insert( filtered_left.end(), filtered_block.begin(), filtered_block.end() );

			//std::cout << "CCCCCCCCCCCC"<< std::endl;
			//std::fill (filtered_block.begin(),filtered_block.end(),0);
			std::vector<float> filtered_block_right(blockSize, 0);
			makeSubList(block_right,audio_right,i*blockSize,(i + 1)*blockSize);
		//
			blockProcessing(h, block_right, state_right, num_taps, filtered_block_right);
			//std::cout << "DDDDDD"<< std::endl;
			filtered_right.insert(filtered_right.end(), filtered_block_right.begin(), filtered_block_right.end() );

		}

		const std::string out_fname = "../data/float32filtered_block.bin";


		write_audio_data(out_fname, filtered_left,filtered_right);

	}
	return 0;
}
