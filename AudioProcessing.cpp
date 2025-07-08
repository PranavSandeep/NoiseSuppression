//
// Created by prana on 07-07-2025.
//
#define DR_WAV_IMPLEMENTATION
#include "AudioProcessing.h"

#include <cmath>
#include <complex>


std::vector<double> AudioProcessing::loadWav(const char *filename, uint32_t &sampleRate)
{
    drwav wav;
    if (!drwav_init_file(&wav, filename))
    {
        std::cerr << "drwav_init_file failed" << std::endl;
        exit(1);
    }

    sampleRate = wav.sampleRate;
    size_t totalSamples = wav.totalSampleCount * wav.channels;

    std::vector<float> tempData(totalSamples);
    drwav_read_f32(&wav, wav.totalSampleCount, tempData.data());
    drwav_uninit(&wav);


    std::vector<double> samples(wav.totalSampleCount);
    for (size_t i = 0; i < wav.totalSampleCount; i++)
    {
        float sum = 0.f;
        for (uint32_t j = 0; j < wav.channels; j++)
        {
            sum+= tempData[i * wav.channels + j];
        }
        samples[i] = static_cast<double>(sum/wav.channels);
    }

    return samples;

}

void AudioProcessing::saveWav(const char *filename, const std::vector<double> &samples, uint32_t sampleRate)
{
    std::vector<float> floatSamples(samples.size());
    for (size_t i = 0; i < samples.size(); i++)
    {
        double s = samples[i];
        if (s < -1) { s = -1; }
        if (s > 1) { s = 1; }
        floatSamples[i] = static_cast<float>(s);
    }

    //Format for wav files:
    drwav_data_format format;
    format.container = drwav_container_riff;
    format.format = DR_WAVE_FORMAT_IEEE_FLOAT;
    format.channels = 1;
    format.sampleRate = sampleRate;
    format.bitsPerSample = 32;

    drwav wav;
    if (!drwav_init_file_write(&wav, filename, &format))
    {
        std::cerr << "drwav_init_file_write failed" << std::endl;
    }

    //Write the audio
    drwav_uint64 framesWritten = drwav_write(&wav, floatSamples.size(), floatSamples.data());
    if (framesWritten != floatSamples.size())
    {
        std::cout << "Only written " << framesWritten << " frames of data" <<std::endl;
    }
    drwav_uninit(&wav);

}

std::vector<double> AudioProcessing::generateWindow(int N)
{
    std::vector<double> window(N);
    for (size_t i = 0; i < N; i++)
    {
        window[i] = 0.54 - 0.46 * cos((2 * M_PI * i) / (N-1));
    }
    return window;
}

std::vector<std::vector<std::complex<double>>> AudioProcessing::generateFrames(
    std::vector<double> &audio, std::vector<double> &window, int frame_size, int hop_size)
{
    std::vector<std::vector<std::complex<double>>> frames;

    for (size_t i = 0; i + frame_size < audio.size(); i+= hop_size)
    {
        std::vector<std::complex<double>> frame(frame_size);
        for (int j = 0; j < frame_size; j++)
        {
            frame[j] = audio[i + j] * window[j];
        }

        size_t padded = nearestPowerOfTwo(frame.size());
        frame.resize(padded, std::complex<double>(0.f, 0.f));
        frames.push_back(frame);
    }
    return frames;
}

void AudioProcessing::FFT(std::vector<std::complex<double>> &frame, bool invert)
{

    int size = frame.size();
    if (size <=1) return;

    std::vector<std::complex<double>> even(size/2), odd(size/2);

    for (int i = 0; i < size/2; i++)
    {
        even[i] = frame[i*2];
        odd[i] = frame[i*2 + 1];
    }

    FFT(even, invert);
    FFT(odd, invert);

    std::complex<double> w;
    w = 1.0;

    std::complex<double> wn = std::polar(1.0, -2*M_PI*(invert ? -1 : 1)/size);
    for (int i = 0; i < size/2; i++)
    {
        frame[i] = even[i] + w * odd[i];
        frame[i + size/2] = even[i] - w * odd[i];
        w *= wn;
        if (invert)
        {
            frame[i] /= 2;
            frame[i + size/2] /= 2;
        }


    }
}

uint32_t AudioProcessing::nearestPowerOfTwo(uint32_t x)
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
