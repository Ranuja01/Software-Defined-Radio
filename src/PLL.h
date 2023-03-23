#ifndef PLL
#define PLL
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstring>
#include <algorithm>
#include <vector>

void fmPLL(std::vector <float>&pllOut, std::vector <float> &pllIn, float freq, float Fs, float ncoScale, float phaseAdjust, float normBandwidth, float &integrator, float* feedbackI, float &feedbackQ, float &trigOffset, float &phaseEst);

#endif
