
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstring>
#include <algorithm>
#include <vector>
#include "PLL.h"

void fmPLL(std::vector <float>&pllOut, std::vector <float> &pllIn, float freq, float Fs, float ncoScale, float phaseAdjust, float normBandwidth,std::vector <float> &pllVariables ){

    /*
    pllIn 	 		array of floats
                                 input signal to the PLL (assume known frequency)

 freq 			float
                                 reference frequency to which the PLL locks

 Fs  			float
                                 sampling rate for the input/output signals

 ncoScale		float
                                 frequency scale factor for the NCO output
                                 = 2 for stereo pll
                                 = 0.5 for rds pll

 phaseAdjust		float
                                 phase adjust to be added to the NCO output only

 normBandwidth	float
                                 normalized bandwidth for the loop filter
                                 (relative to the sampling rate)

 state 			to be added*/
/*
 scale factors for proportional/integrator terms
 these scale factors were derived assuming the following:
 damping factor of 0.707 (1 over square root of 2)
 there is no oscillator gain and no phase detector gain
*/

 float Cp = 2.666;
 float Ci = 3.555;

 float pi = 4.0*atan(1);
 	// gain for the proportional term
 	float Kp = (normBandwidth)*Cp;
 	// gain for the integrator term
 	float Ki = (normBandwidth*normBandwidth)*Ci;


  std::vector <float> ncoOut(pllIn.size()+1, 0);

  	// initialize internal state
    //before block process, the pointers below were initialized to 0
    //float integrator = 0.0
  	//float phaseEst = 0.0
  	//float feedbackI = 1.0
  	//float feedbackQ = 0.0
    //float trigOffset = 0

    ncoOut[0] = 1.0;
  	// note: state saving will be needed for block processing
    float errorI = 0 ;
    float errorQ = 0;

    float errorD = 0;
    float trigArg = 0;

    float kerr = 0;

    //float integ = 0;
    /*
    pllVariables[i] to state save for next block
    i = 0, float integrator = 0;
    i = 1, float phaseEst = 0;
    i = 2, float feedbackI = 1.0;
    i = 3, float feedbackQ = 0.0;
    i = 4, float trigOffset = 0.0;
    */

    for(int i=0; i<pllIn.size(); i++){
      // phase detector
  		errorI = pllIn[i] * (pllVariables[2]);  // complex conjugate of the
  		errorQ = pllIn[i] * (-1.0*(pllVariables[3]));  // feedback complex exponential

      // four-quadrant arctangent discriminator for phase error detection
      //account for divide by 0
      if (errorI == 0){
        errorD = 2*atan(1);
      }else{
  		    errorD = atan(errorQ/errorI);
      }

      pllVariables[0] += Ki*errorD;

  		// update phase estimate
      pllVariables[1] +=  Kp*errorD + pllVariables[0];

  		// internal oscillator

  		pllVariables[4] += 1;
  		trigArg = 2*pi*(freq/Fs)*(pllVariables[4]);
  		pllVariables[2] = cos(trigArg); //in phase
  		pllVariables[3] = sin(trigArg); //quartature
  		ncoOut[i+1] = cos(trigArg*ncoScale + pllVariables[1] + phaseAdjust);

    }
    pllOut = ncoOut;
  }
