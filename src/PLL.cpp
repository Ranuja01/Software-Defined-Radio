
#include <vector>


void fmPLL(std::vector <float>pllIn, float freq, float Fs, float ncoScale = 2.0, float phaseAdjust = 0.0, float normBandwidth = 0.01, std::vector <float> PLLwave, float &integrator, float &feedbackI, float &feedbackQ, float &trigOffset, float &phaseEst){

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

    ncoOut[0] = 1.0
  	// note: state saving will be needed for block processing
    float errorI = 0 ;
    float errorQ = 0;

    for(int i=0; i<pllIn.size() i++){

      // phase detector
  		errorI = pllIn[i] * (feedbackI);  # complex conjugate of the
  		errorQ = pllIn[i] * (-1.0*feedbackQ);  # feedback complex exponential

      // four-quadrant arctangent discriminator for phase error detection
  		errorD = atan(errorQ/errorI);

  		//loop filter
  		integrator = integrator + Ki*errorD;

  		// update phase estimate
  		phaseEst = phaseEst + Kp*errorD + integrator;

  		// internal oscillator
  		trigOffset += 1;
  		trigArg = 2*pi*(freq/Fs)*(trigOffset);
  		feedbackI = cos(trigArg);
  		feedbackQ = sin(trigArg);
  		ncoOut[i+1] = cos(trigArg*ncoScale + phaseEst + phaseAdjust);

  	// for stereo only the in-phase NCO component should be returned
  	// for block processing you should also return the state
  	//return ncoOut
  	//for RDS add also the quadrature NCO component to the output

    }
    return ncoOut;
    //PLLwave.insert(PLLwave.end(), ncoOut.begin(), ncoOut.end());
}
