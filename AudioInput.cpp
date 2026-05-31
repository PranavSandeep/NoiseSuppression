//
// Created by prana on 29-07-2025.
//


#include "AudioInput.h"

#include "AudioProcessing.h"


#include "concurrentqueue.h"
#include "librosa/librosa.h"

#include <atomic>
#include <chrono>
#include <thread>
#include <vector>

constexpr int FRAME_SIZE  = 512;
constexpr int HOP_SIZE    = 256;
constexpr int MODEL_INPUT = 1024;

std::vector<double> window = AudioProcessing().generateWindow(MODEL_INPUT);

std::vector<float> in(FRAME_SIZE);
std::vector<float> out(FRAME_SIZE);

moodycamel::ConcurrentQueue<std::vector<float>> outputQueue;
moodycamel::ConcurrentQueue<std::vector<float>> inputQueue;



std::atomic<bool> AudioInput::running = true;



bool playbackStarted = false;

void AudioInput::data_callback(
    ma_device* device,
    void* output,
    const void* input,
    ma_uint32 frame_count)
{
    std::memcpy(
        in.data(),
        input,
        frame_count *
        ma_get_bytes_per_frame(
            device->capture.format,
            device->capture.channels));

    inputQueue.enqueue(in);

    if (outputQueue.size_approx() > 0 || playbackStarted)
    {
        playbackStarted = true;

        if (outputQueue.try_dequeue(out))
        {
            std::memcpy(
                output,
                out.data(),
                frame_count *
                ma_get_bytes_per_frame(
                    device->playback.format,
                    device->playback.channels));
        }


    }
}

void AudioInput::proccessAudio()
{

    constexpr int CHUNK   = 512;   // receive 512 at a time
    constexpr int MODEL_INPUT = 1024;
    constexpr int FFT     = 512;
    constexpr int HOP     = 256;

    std::vector<float> contextBuffer(MODEL_INPUT * 2, 0.0f);
    std::vector<float> accumBuffer;
    accumBuffer.reserve(MODEL_INPUT);



    std::vector<float> ola_buffer(MODEL_INPUT * 2, 0.0f);

    //Take larger context, and save history to try and fix the crackle. Gave some success.
    Eigen::MatrixXcf stftHistory = Eigen::MatrixXcf::Zero(257, 9);
    while (running.load(std::memory_order_relaxed))
    {
        std::vector<float> chunk512;

        if (!inputQueue.try_dequeue(chunk512))
        {
            std::this_thread::yield();
            continue;
        }

        if (chunk512.size() != (size_t)CHUNK)
            continue;


        accumBuffer.insert(
            accumBuffer.end(),
            chunk512.begin(),
            chunk512.end());

        if (accumBuffer.size() < (size_t)MODEL_INPUT)
            continue;


        std::memmove(
            contextBuffer.data(),
            contextBuffer.data() + MODEL_INPUT,
            MODEL_INPUT * sizeof(float));

        std::memcpy(
            contextBuffer.data() + MODEL_INPUT,
            accumBuffer.data(),
            MODEL_INPUT * sizeof(float));

        accumBuffer.clear();


        auto STFT = librosa::Feature::stft(
            contextBuffer, FFT, HOP, "hann", true, "constant");

        auto STFT_Matrix = librosa::internal::toEigen(STFT);
        Eigen::MatrixXcf S = STFT_Matrix.transpose();


        int halfCols = S.cols() / 2 + 1;
        Eigen::MatrixXcf S_curr = S.rightCols(halfCols);

        Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> mag
            = S_curr.cwiseAbs();
        Eigen::ArrayXXcf phase
            = S_curr.array() / (S_curr.array().abs() + 1e-8f);


        Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> output;
        md->fast_predict(mag.data(), mag.size(), output);

        const float current_threshold = noiseThreshold.load(std::memory_order_relaxed);
        const float current_attentuation = floorAttenuation.load(std::memory_order_relaxed);


        output = output.unaryExpr([current_threshold, current_attentuation](float maskValue) {

            if (maskValue < current_threshold) {
                return maskValue * current_attentuation;
            }
            return maskValue;
        });

        Eigen::MatrixXcf S_reconstructed = output.array() * phase;

        stftHistory.leftCols(halfCols) = stftHistory.rightCols(halfCols).eval();

        stftHistory.rightCols(halfCols) = S_reconstructed;

        std::vector<float> reconstructed = librosa::Feature::istft(
            stftHistory, FFT, HOP, "hann", true, "constant");


        //OLA

        for (int  i = 0; i < reconstructed.size(); i++)
        {
            ola_buffer[i] += reconstructed[i];
        }


        outputQueue.enqueue(std::vector<float>(
            ola_buffer.begin(),
            ola_buffer.begin() + CHUNK));

        outputQueue.enqueue(std::vector<float>(
            ola_buffer.begin() + CHUNK,
            ola_buffer.begin() + MODEL_INPUT));



        std::memmove(
            ola_buffer.data(),
            ola_buffer.data() + MODEL_INPUT,
            MODEL_INPUT * sizeof(float)

            );

        std::fill(ola_buffer.begin() + MODEL_INPUT, ola_buffer.end(), 0.0f);


    }
}


