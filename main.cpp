#include <iostream>
#include "AudioProcessing.h"
#include <complex>


uint32_t nearestPowerOfTwo(uint32_t x)
{
    if (x == 0) return 0;
    --x;

    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    return x+1;
}

int main()
{
    AudioProcessing audioProcesser;

    uint32_t sampleRate;
    auto samples = audioProcesser.loadWav(R"(C:\Users\prana\Downloads\Test.wav)", sampleRate);

    std::vector<double> window = audioProcesser.generateWindow(512);

    auto frames = audioProcesser.generateFrames(samples, window, 512, 256);

    size_t padded_size = frames[0].size();
    std::vector<double> output(samples.size(), 0.0);

    for (size_t i = 0; i < frames.size(); ++i)
    {
        auto& frame = frames[i];

        audioProcesser.FFT(frame, false);

        //TODO Do spectral noise subtraction here.

        audioProcesser.FFT(frame, true);

        for (int j = 0; j < 512; ++j)
        {
            output[i*256 + j] += frame[j].real() * window[j];
        }
    }

    audioProcesser.saveWav(R"(C:\Users\prana\Downloads\Output.wav)", output, sampleRate);





    return 0;

}