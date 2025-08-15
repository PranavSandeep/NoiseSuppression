//
// Created by prana on 29-07-2025.
//
#define MINIAUDIO_IMPLEMENTATION

#include "AudioInput.h"

#include <chrono>
#include "AudioProcessing.h"
#include "Model.h"
#include <deque>

std::vector<float> inputBuffer;
std::vector<float> outputBuffer;
std::vector<double> window = AudioProcessing().generateWindow(512);
int FRAME_SIZE = 512;
int HOP_SIZE = 256;
size_t writeIndex = 0;
std::vector<float> phase;

auto model = Model(L"model-QuickRCED.onnx");

int frame_c = 0;

std::vector<float> buffer(257*251, 0);

void AudioInput::data_callback(ma_device *device, void *output, const void *input, ma_uint32 frame_count)
{

    auto start = std::chrono::high_resolution_clock::now();
    const float* in = (const float*)input;
    float* out = (float*)output;


    inputBuffer.insert(inputBuffer.end(), in, in + frame_count);
    AudioProcessing audio_processor;
    while (inputBuffer.size() >= FRAME_SIZE)
    {


        std::vector<std::complex<double>> frame(FRAME_SIZE);

        //Apply window
        for (int i = 0; i < FRAME_SIZE; i++)
            frame[i] = inputBuffer[i] * window[i];





        audio_processor.FFT(frame, false);

        buffer.erase(buffer.begin(), buffer.begin() + FRAME_SIZE);
        for (auto i : frame)
        {
            buffer.push_back(std::abs(i));
            phase.push_back(std::arg(i));
        }

        frame_c++;
        if(frame_c >=64)
        {
            auto start1 = std::chrono::high_resolution_clock::now();
            model.predict(buffer);
            auto end1 = std::chrono::high_resolution_clock::now();
            auto duration1 = std::chrono::duration_cast<std::chrono::microseconds>(end1 - start1);
            //std::cout << duration1.count() << " microseconds\n";

            frame_c = 0;
        }

        audio_processor.FFT(frame, true);


        //Overlap add
        for (int i = 0; i < FRAME_SIZE; i++)
        {
            double val = frame[i].real() * window[i];
            size_t  outIndex = writeIndex + i;
            if (outputBuffer.size()  <= outIndex)
                outputBuffer.resize(outIndex + 1, 0.f);
            outputBuffer[outIndex] += static_cast<float>(val);
        }
        writeIndex += HOP_SIZE;

        inputBuffer.erase(inputBuffer.begin(), inputBuffer.begin() + HOP_SIZE);

    }

    if (outputBuffer.size() >= frame_count) {
        for (uint32_t i = 0; i < frame_count; ++i) {
            out[i] = outputBuffer[i];
        }

    } else {
        for (uint32_t i = 0; i < frame_count; ++i)
            out[i] = 0.0f;
    }


    //Remove written bins from output buffer
    if (outputBuffer.size() >= frame_count)
        outputBuffer.erase(outputBuffer.begin(), outputBuffer.begin() + frame_count);

    if (writeIndex >= frame_count)
        writeIndex -= frame_count;
    else
        writeIndex = 0;



    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    // std::cout << duration.count() << " microseconds\n";


}

