/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#ifndef DY4_FILTER_H
#define DY4_FILTER_H

// add headers as needed
#include <iostream>
#include <vector>

// declaration of a function prototypes
void makeSubList (std::vector<float> &subList, const std::vector<float> &list, int first, int last);
void makeOddEvenSubList (std::vector<float> &subList, const std::vector<float> &list, int first, int last);
void blockConvolve(std::vector<float> &h, const std::vector<float> &block, std::vector<float> &state, int num_taps, std::vector<float> &filtered_block);
void blockProcessing(std::vector<float> &h, const std::vector<float> &block, std::vector<float> &state, int num_taps, std::vector<float> &filtered_block, unsigned short int &startIndex, int dRate);
void blockResample(std::vector<float> &h, const std::vector<float> &block, std::vector<float> &state, int num_taps, std::vector<float> &filtered_block, int dRate, int uRate);
void fmDemod (std::vector<float> &demodulatedSignal, const std::vector<float> &I, const std::vector<float> &Q, float &prevI, float &prevQ);
void impulseResponseLPF(float Fs, float Fc, unsigned short int num_taps, std::vector<float> &h, const int &gain);
void impulseResponseBPF(float Fs, float Fb, float Fe, unsigned short int num_taps, std::vector<float> &h, const int &gain);
void convolveFIR(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h);
void impulseResponseRootRaisedCosine(float Fs, unsigned short int num_taps, std::vector<float> &h);


#endif // DY4_FILTER_H
