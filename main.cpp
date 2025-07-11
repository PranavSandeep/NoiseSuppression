#include <iostream>
#include "AudioProcessing.h"
#include "NoiseSuppression.h"
#include <complex>



int main()
{
    AudioProcessing audioProcesser;

    uint32_t sampleRate;
    auto samples = audioProcesser.loadWav(R"(C:\Users\prana\Downloads\Test3JG.wav)", sampleRate);

    std::vector<double> window = audioProcesser.generateWindow(512);

    auto frames = audioProcesser.generateFrames(samples, window, 512, 256);

    size_t padded_size = frames[0].size();
    std::vector<double> output(samples.size(), 0.0);

    std::vector<double> noise = NoiseSuppression().generateAverageNoise(frames, 150, 512);

    for (size_t i = 0; i < frames.size(); ++i)
    {
        auto& frame = frames[i];

        audioProcesser.FFT(frame, false);
        double freqBinWidth = sampleRate / frame.size();



        NoiseSuppression().subtractSpectralNoise(frame, noise, freqBinWidth);


        audioProcesser.FFT(frame, true);


        for (int j = 0; j < 512; ++j)
        {
            output[i*256 + j] += frame[j].real() * window[j];
        }
    }

    output.resize(output.size()/2);

    for (auto& sample : output)
        sample *= 1.5;

    audioProcesser.saveWav(R"(C:\Users\prana\Downloads\Output10.wav)", output, sampleRate);





    return 0;

}