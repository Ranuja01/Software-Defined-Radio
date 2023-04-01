/*
Comp Eng 3DY4 (Computer Systems Integration Project)

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

// Makes sublist
void makeSubList (std::vector<float> &subList, const std::vector<float> &list, int first, int last){
  subList.clear();
  for (int i = first; i < last; i++){
    subList.push_back(list[i]);
  }

}

// Convolution with no up/down samples
void blockConvolve(std::vector<float> &h, const std::vector<float> &block, std::vector<float> &state, int num_taps, std::vector<float> &filtered_block){

  for (int n = 0; n < block.size(); n++){
		for (int k = 0; k < h.size(); k++){
			if (n - k >=0){
				if (n - k < block.size()){
					filtered_block[n] += h[k] * block[n-k]; //convolution over block samples
				}
			}else {
		  	if (n - k + num_taps - 1 < state.size()){
					filtered_block[n] += h[k] * state[(n - k) + num_taps - 1]; //convolution over states
				}
			}
		}
	}
  // Makes the state as a sublist
	makeSubList(state,block,block.size() - num_taps + 1, block.size());
}

// Block convolution with downsample
void blockProcessing(std::vector<float> &h, const std::vector<float> &block, std::vector<float> &state, int num_taps, std::vector<float> &filtered_block, unsigned short int &startIndex, int dRate){
  // Holds the end index
  int endIndex = 0;
  // Begins at the startIndex since each block would not start at 0
  // Increments the loop variable by the downsample rate
  for (int n = startIndex; n < block.size(); n+=dRate){
    for (int k = 0; k < h.size(); k++){
      if (n - k >= 0){
        if (n - k < block.size()){
          // Subtract the index by the start and divide by the downsample rate to normalize
          filtered_block[(n - startIndex)/dRate] += h[k] * block[n-k];
        }
      }else {
        if (n - k + num_taps - 1 < state.size()){
          filtered_block[(n - startIndex)/dRate] += h[k] * state[(n - k) + num_taps - 1];
        }
      }

    }
    endIndex = n;
  }
  // Define the start index for the next block
  startIndex = (endIndex - block.size()) + dRate;
  // Creates the state
  makeSubList(state,block,block.size() - num_taps + 1, block.size());
}

// Block convolution with resampling
void blockResample(std::vector<float> &h, const std::vector<float> &block, std::vector<float> &state, int num_taps, std::vector<float> &filtered_block, int dRate, int uRate){
	int p = 0;
  int j = 0;

  // Start at the state size, such that any variable under the size would represent state values
  // This is so that negative indices won't be acquired
	for (int n = state.size(); n < block.size() + state.size(); n++){
		// Acquire the phase
    // This defines the start position of the base vector
    p = (n * dRate)%uRate;
		for (int k = p; k < h.size(); k+= uRate){

      // By multiplying by the base index, subtracting the offset k and divide by the upsample uRate
      // This normalizes the index to where it should be
      // This is the fast resample
      j = (int)((n*dRate - k)/uRate);

      // Subtract by state size to normalize the offset done at the beginning
			if (j >= 0) {    
        if (j < state.size()){
          filtered_block[n - state.size()] += h[k] * state[j - state.size()]; // if index within state size, use value of state
        } else if (j < block.size()){ 
          filtered_block[n - state.size()] += h[k] * block[j - state.size()]; //if index is within block size, use value of block
        }
      }
		}
	}
	makeSubList(state,block,block.size() - num_taps + 1, block.size());
}


void makeOddEvenSubList (std::vector<float> &subList, const std::vector<float> &list, int first, int last){
  subList.clear();
  for (int i = first; i < last; i+=2){
    subList.push_back(list[i]); //adds current element to sub list
  }

}

void fmDemod (std::vector<float> &demodulatedSignal, const std::vector<float> &I, const std::vector<float> &Q, float &prevI, float &prevQ){
  std::fill (demodulatedSignal.begin(),demodulatedSignal.end(),0); //initialized demod signal with zeros
  demodulatedSignal[0] = 0;

  for (int k = 0; k < I.size(); k++){
    if (!( (I[k] * I[k]) + (Q[k] * Q[k]) == 0)){ //checks if denominator is not 0
      demodulatedSignal[k] = (1.0/( (I[k] * I[k]) + (Q[k] * Q[k]) ) ) * (I[k] * (Q[k] - prevQ) - Q[k] * (I[k] - prevI)); //calculates demod signal using prev IQ values

    }else {
      demodulatedSignal[k] = 0; //if denominator is 0, signal is 0
    }
    prevI = I[k]; //set prev I value to current val
    prevQ = Q[k]; //set prev Q value to current val
  }

}

// function for computing the impulse response (reuse from previous experiment)
void impulseResponseLPF(float Fs, float Fc, unsigned short int num_taps, std::vector<float> &h,const int &gain)
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
    h[i] *= gain; //apply gain
  }
}

void impulseResponseBPF(float Fs, float Fb, float Fe, unsigned short int num_taps, std::vector<float> &h,const int &gain)
{
	// allocate memory for the impulse response
	h.clear(); h.resize(num_taps, 0.0);

	// the rest of the code in this function is to be completed by you
	// based on your understanding and the Python code from the first lab

	float norm_center = ((Fe + Fb)/2)/(Fs/2);
  float norm_pass = (Fe - Fb)/(Fs/2);

	for (int i = 0; i < num_taps; i++){
    //  ****TRY SPLITTING IF STATEMENT INTO TWO FOR LOOPS**** - from previous exp
		if (i == (num_taps - 1)/2){
			h[i] = norm_pass;
		} else {
			h[i] = norm_pass * sin (PI * norm_pass * (i - (num_taps - 1)/2)) / (PI * norm_pass * (i - (num_taps - 1)/2));
		}
    h[i] *= cos(i*PI*norm_center);
		h[i] *= pow(sin(PI * i / num_taps),2) * gain;

	}
}

void impulseResponseRootRaisedCosine(float Fs, unsigned short int num_taps, std::vector<float> &h)
{
    // allocate memory for the impulse response
    h.clear(); h.resize(num_taps, 0.0);
  /*
  """
  Root raised cosine (RRC) filter

  Fs  		sampling rate at the output of the resampler in the RDS path
        sampling rate must be an integer multipler of 2375
        this integer multiple is the number of samples per symbol

  N_taps  	number of filter taps

  """
  */

  // duration for each symbol - should NOT be changed for RDS!
  float T_symbol = 1/2375.0;
  float t = 0;
  // roll-off factor (must be greater than 0 and smaller than 1)
  // for RDS a value in the range of 0.9 is a good trade-off between
  // the excess bandwidth and the size/duration of ripples in the time-domain
  float beta = 0.90;

  // the RRC inpulse response that will be computed in this function
  //impulseResponseRRC = np.empty(N_taps)
//  RRCimpulse


  //for k in range(N_taps):
    for(int i = 0; i<num_taps; i++){ 
    t = ((i-num_taps/2))/Fs; //calculates time offset
    // we ignore the 1/T_symbol scale factor
    
    if (t == 0.0){ //checks for case where t =0
      h[i] = 1.0 + beta*((4/PI)-1); //calcukates filter coeff
    }

    else if (t == -T_symbol/(4*beta) or t == T_symbol/(4*beta)){ 
      h[i] = (beta/sqrt(2))*(((1+2/PI)* \
          (sin(PI/(4*beta)))) + ((1-2/PI)*(cos(PI/(4*beta)))));
        }

    else{ h[i] = (sin(PI*t*(1-beta)/T_symbol) + 4*beta*(t/T_symbol)*cos(PI*t*(1+beta)/T_symbol))/(PI*t*(1-(4*beta*t/T_symbol)*(4*beta*t/T_symbol))/T_symbol);
        }
  }
  // returns the RRC impulse response to be used by convolution
  //return impulseResponseRRC
}

// function for computing the impulse response - reused from previous experiment
void convolveFIR(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h)
{

  // allocates memory for the output (filtered) data
  y.clear(); y.resize(x.size()+h.size()-1, 0.0); // resizes output vector
  for (int n = 0; n < y.size(); n++){
    for (int k = 0; k < h.size(); k++){
      if (n - k >= 0 && n - k < x.size()){ // if index is within bounds
        y[n] += h[k] * x[n-k]; //convolve filter coeff with input
      }
    }
  }
}
