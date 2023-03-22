#ifndef PLL
#define PLL

#include <vector>

std::vector <float> fmPLL(std::vector <float>pllIn, float freq, float Fs, float ncoScale, float phaseAdjust, float normBandwidth, float &integrator, float &feedbackI, float &feedbackQ, float &trigOffset, float &phaseEst);

#endif
