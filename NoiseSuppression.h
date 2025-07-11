//
// Created by prana on 08-07-2025.
//

#ifndef NOISESUPPRESSION_H
#define NOISESUPPRESSION_H
#include <complex>
#include <vector>


class NoiseSuppression {
public:
    std::vector<double> generateAverageNoise(std::vector<std::vector<std::complex<double> > > &frames,
                                             int initialFrames, int frame_size);

    void subtractSpectralNoise(std::vector<std::complex<double> > &frames, std::vector<double> &spectralNoise, double freqBinWidth);

};



#endif //NOISESUPPRESSION_H
