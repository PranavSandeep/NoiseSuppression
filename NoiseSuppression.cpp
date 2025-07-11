//
// Created by prana on 08-07-2025.
//

#include "NoiseSuppression.h"
#include "AudioProcessing.h"

std::vector<double> NoiseSuppression::generateAverageNoise(
    std::vector<std::vector<std::complex<double> > > &frames, int initialFrames, int frame_size)
{
    std::vector<double> noiseFrames(frame_size, 0.0);
    for (int i = 0; i < initialFrames; i++)
    {
        auto& frame = frames[i];
        AudioProcessing().FFT(frame, false);

        for (int j = 0; j < frame.size(); j++)
        {
            noiseFrames[j] += std::abs(frame[j]);
        }
    }

    for (int k = 0; k < frame_size; ++k)
    {
        noiseFrames[k] /= 70.0;
    }

    return noiseFrames;

}

void NoiseSuppression::subtractSpectralNoise(std::vector<std::complex<double>> &frames,
                                             std::vector<double> &spectralNoise, double freqBinWidth)
{

    for (int i = 0; i < frames.size(); i++)
    {
        double noiseMag = spectralNoise[i];
        double mag = std::abs(frames[i]);
        double phase = std::arg(frames[i]);




        double cleanedMag = std::max(0.0, mag - noiseMag);
        double ratio = mag / (noiseMag + 1e-8);
        if (ratio < 1.5) {  // You can tune this threshold
            cleanedMag = 0.0;
        }
        frames[i] = std::polar(cleanedMag, phase);
    }
}
