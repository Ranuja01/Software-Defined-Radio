#ifndef DY4_IOFUNC_H
#define DY4_IOFUNC_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>

//void write_stereo_data(std::vector<float>& , float, std::vector<short int>&);
void write_stereo_data(std::vector<float> &audio, float audio_Fs);
void read_audio_data(const std::string in_fname, std::vector<uint8_t> &audio_data);
void write_audio_data(std::vector<float> &audio, float audio_Fs);

#endif // DY4_IOFUNC_H
