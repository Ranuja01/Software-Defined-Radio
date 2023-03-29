
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

    //before block process, initialize the pointers below to 0
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
    float integrator = 0;
    float phaseEst = 0;
    float feedbackI = 1.0;
    float feedbackQ = 0.0;
    float phaseadjust = 0.0;
    */

    for(int i=0; i<pllVariables.size(); i++){
      //  std::cout << i << "asdasd: " << pllVariables[i] << std::endl;
    }

    for(int i=0; i<pllIn.size(); i++){
      //std::cout << "3: " << integrator << std::endl;
      // phase detector
  		errorI = pllIn[i] * (pllVariables[2]);  // complex conjugate of the
  		errorQ = pllIn[i] * (-1.0*(pllVariables[3]));  // feedback complex exponential

      // four-quadrant arctangent discriminator for phase error detection
      if (errorI == 0){
        errorD = 2*atan(1);
      }else{
  		    errorD = atan(errorQ/errorI);
      }
  		//loop filter
      //integrator = integ;


      pllVariables[0] += Ki*errorD;

      //pllIn[0] += errorD;
  		// update phase estimate
      pllVariables[1] +=  Kp*errorD + pllVariables[0];
  		//*phaseEst = *phaseEst + Kp*errorD + *integrator;
/*
      std::cout << "1: " << pllIn[0] << std::endl;
      std::cout << "2: " << errorQ << std::endl;
      std::cout << "3: " << pllVariables[0] << std::endl;
      std::cout << "4: " << pllVariables[1]<< std::endl;
      std::cout << "5: " << Ki*errorD << std::endl;*/
  		// internal oscillator
      //std::cout << "fdes: " << pllIn[1024] << std::endl;
  		pllVariables[4] += 1;
  		trigArg = 2*pi*(freq/Fs)*(pllVariables[4]);
  		pllVariables[2] = cos(trigArg);
  		pllVariables[3] = sin(trigArg);
  		ncoOut[i+1] = cos(trigArg*ncoScale + pllVariables[1] + phaseAdjust);

  	// for stereo only the in-phase NCO component should be returned
  	// for block processing you should also return the state
  	//return ncoOut
  	//for RDS add also the quadrature NCO component to the output
    //std::cout<< "d:   " << ncoOut[i+1] << std::endl;
    }
    pllOut = ncoOut;
    //PLLwave.insert(PLLwave.end(), ncoOut.begin(), ncoOut.end());
}
