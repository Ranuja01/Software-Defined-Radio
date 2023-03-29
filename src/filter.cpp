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


void makeSubList (std::vector<float> &subList, const std::vector<float> &list, int first, int last){
  subList.clear();
  for (int i = first; i < last; i++){
    subList.push_back(list[i]);
  }

}

//sublist for state saving for demodulation
void makeOddEvenSubList (std::vector<float> &subList, const std::vector<float> &list, int first, int last){
  subList.clear();
  for (int i = first; i < last; i+=2){
    subList.push_back(list[i]);
  }

}

void impulseResponseLPF(float Fs, float Fc, unsigned short int num_taps, std::vector<float> &h)
{
  // allocate memory for the impulse response
  h.clear(); h.resize(num_taps, 0.0);

  float norm_cutoff = Fc/(Fs/2);

  //low pass filter computation
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


void blockProcessing(std::vector<float> &h, const std::vector<float> &block, std::vector<float> &state, int num_taps, std::vector<float> &filtered_block, unsigned short int &startIndex, int dRate){
  int endIndex = 0;

  //convolution
  for (int n = startIndex; n < block.size(); n+=dRate){

    for (int k = 0; k < h.size(); k++){
      if (n - k >= 0){
        if (n - k < block.size()){
          filtered_block[(n - startIndex)/dRate] += h[k] * block[n-k];
        }
      }
      //use the saved state
      else {
        if (n - k + num_taps - 1 < state.size()){
          filtered_block[(n - startIndex)/dRate] += h[k] * state[(n - k) + num_taps - 1];
        }
      }
    }
    endIndex = n;
  }

  //save the start index for the next block
  startIndex = (endIndex - block.size()) + dRate;
  makeSubList(state,block,block.size() - num_taps + 1, block.size());
}

void fmDemod (std::vector<float> &demodulatedSignal, const std::vector<float> &I, const std::vector<float> &Q, float &prevI, float &prevQ){
  std::fill (demodulatedSignal.begin(),demodulatedSignal.end(),0);
  demodulatedSignal[0] = 0;

  //difference equation for demodulation
  for (int k = 0; k < I.size(); k++){
    if (!( (I[k] * I[k]) + (Q[k] * Q[k]) == 0)){
      demodulatedSignal[k] = (1.0/( (I[k] * I[k]) + (Q[k] * Q[k]) ) ) * (I[k] * (Q[k] - prevQ) - Q[k] * (I[k] - prevI));
    }else {
      demodulatedSignal[k] = 0;
    }

    //state saving
    prevI = I[k];
    prevQ = Q[k];
  }

}

void blockResample(std::vector<float> &h, const std::vector<float> &block, std::vector<float> &state, int num_taps, std::vector<float> &filtered_block, int dRate, int uRate){
	int p = 0;
  int j = 0;

  //convolution
	for (int n = state.size(); n < block.size() + state.size(); n++){
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


void blockConvolve(std::vector<float> &h, const std::vector<float> &block, std::vector<float> &state, int num_taps, std::vector<float> &filtered_block){

  std::cout << filtered_block.size() << " " << h.size() << "  " << block.size() << std::endl;
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

  std::cout << "blockconv" << std::endl;
	makeSubList(state,block,block.size() - num_taps + 1, block.size());
}




// function to compute the filtered output "y" by doing the convolution
// of the input data "x" with the impulse response "h
