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

//  std::cout << "blockconv" << std::endl;
	makeSubList(state,block,block.size() - num_taps + 1, block.size());
}

//blockProcessing(rf_coeff, block, state_i, rf_taps, filtered_i,startIndex_i,rf_decim);
void blockProcessing(std::vector<float> &h, const std::vector<float> &block, std::vector<float> &state, int num_taps, std::vector<float> &filtered_block, unsigned short int &startIndex, int dRate){
  int endIndex = 0;
  for (int n = 0; n < block.size(); n+=dRate){
    //std::cout << "bbbb: " << n << std::endl;
    for (int k = 0; k < h.size(); k++){
      if (n - k >= 0){
        if (n - k < block.size()){
          filtered_block[(n - startIndex)/dRate] += h[k] * block[n-k];
          if(k == 0 && n-k == 0){
          //	std::cout << "bbbb: " <<  h[k] << std::endl;
          //	std::cout << "cccc: " <<  block[n-k] << std::endl;
          }
          if(k == 0 && n-k == 10){
          //	std::cout << "dddd: " <<  h[k] << std::endl;
          //	std::cout << "eeee: " <<  block[n-k] << std::endl;
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
  //startIndex = 1 + (endIndex + block.size() % dRate) - block.size();
  startIndex = (endIndex - block.size()) + dRate;
/*
  for (int i = filtered_block.size() - 5140; i < filtered_block.size() - 5100;i++){
    //std::cout << "aaaa: " << filtered_block[i] << std::endl;
    //std::cout << "bbbb: " <<  block[i] << std::endl;
  //	std::cout << "cccc: " <<  h[i] << std::endl;

}
  /*
  for (int i = 51000; i < block.size();i++){
  //	std::cout << "aaaa: " << filtered_block[i] << std::endl;
    //std::cout << "bbbb: " <<  block[i] << std::endl;

      std::cout << i << " dddd: " <<  block[i] << std::endl;
      std::cout << (i)/dRate << " eeee: " <<  filtered_block[(i)/dRate] << std::endl;
  }*/
  //std::cout << "asdfghjkl: " << filtered_block[filtered_block.size() - 2] << std::endl;
  //std::cout << startIndex << std::endl;
  //startIndex = 0;
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
    h[i] *= gain;
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
