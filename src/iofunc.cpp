
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <vector>
#include <complex>
#include <cmath>
void write_stereo_data(std::vector<float> &audio, float audio_Fs, std::vector<short int> &play)
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
