//
// Created by prana on 07-07-2025.
//


#ifndef AUDIOPROCESSING_H
#define AUDIOPROCESSING_H

#include <complex>
#include <iostream>
#include <ostream>

#include "dr_wav.h"
#include <vector>


class AudioProcessing {

    public:
    std::vector<double> loadWav(const char* filename,  uint32_t& sampleRate);

    void saveWav(const char* filename, const std::vector<double> &samples, uint32_t sampleRate);

    std::vector<double> generateWindow(int N);

    std::vector<std::vector<std::complex<double>>> generateFrames(std::vector<double> &audio, std::vector<double> &window, int frame_size, int hop_size);

    void FFT(std::vector<std::complex<double>> &frame, bool invert);

private:

    uint32_t nearestPowerOfTwo(uint32_t x);


};



#endif //AUDIOPROCESSING_H
