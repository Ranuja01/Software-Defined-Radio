#
# Comp Eng 3DY4 (Computer Systems Integration Project)
#
# Copyright by Nicola Nicolici
# Department of Electrical and Computer Engineering
# McMaster University
# Ontario, Canada
#

import numpy as np
import math

void impulseResponseRootRaisedCosine(float Fs, float Fc, unsigned short int num_taps, std::vector<float> &h)
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

  for(int i = 0; i<N_taps; i++){
    t = ((i-N_taps/2))/Fs;
    // we ignore the 1/T_symbol scale factor
    if (t == 0.0){
      h[i] = 1.0 + beta*((4/PI)-1);
    }

    else if (t == -T_symbol/(4*beta) or t == T_symbol/(4*beta)){
      h[i] = (beta/sqrt(2))*(((1+2/PI)* \
          (sin(PI/(4*beta)))) + ((1-2/PI)*(cos(PI/(4*beta)))));
        }
    else{ h[i] = (sin(PI*t*(1-beta)/T_symbol) +  \
          4*beta*(t/T_symbol)*cos(PI*t*(1+beta)/T_symbol))/ \
          (PI*t*(1-(4*beta*t/T_symbol)*(4*beta*t/T_symbol))/T_symbol);
        }
  }
  // returns the RRC impulse response to be used by convolution
  //return impulseResponseRRC
}
