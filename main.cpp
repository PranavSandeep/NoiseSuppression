#include <iostream>
#include "AudioProcessing.h"
#include "NoiseSuppression.h"
#include <complex>
#include "Model.h"
#include "AudioInput.h"

#define SAMPLE_RATE 0
#define CHANNELS 1
#define FRAMES_PER_BUFF 1024


void denoiseAudio();

int main()
{
    ma_device device;
    ma_device_config device_config;

    device_config = ma_device_config_init(ma_device_type_duplex); //Input and output

    device_config.capture.format = ma_format_f32;
    device_config.capture.channels = CHANNELS;

    device_config.sampleRate = SAMPLE_RATE;

    device_config.playback.channels = CHANNELS;
    device_config.playback.format = ma_format_f32;

    device_config.dataCallback = AudioInput::data_callback;
    device_config.pUserData = nullptr;

    if (ma_device_init(nullptr, &device_config, &device) != MA_SUCCESS)
    {
        std::cerr << "Failed to initialize device \n";
        return -1;
    }

    std::cout << "Starting duplex stream \n";
    ma_device_start(&device);

    std::cout << "Press Enter to stop\n";
    std::cin.get();

    ma_device_uninit(&device);




    return 0;


}

void denoiseAudio()
{
    AudioProcessing audioProcesser;

    uint32_t sampleRate;
    auto samples = audioProcesser.loadWav(R"(audioPath.wav)", sampleRate);
    auto Airnoise = audioProcesser.loadWav(R"(noisePath.wav)", sampleRate);

    std::vector<double> window = audioProcesser.generateWindow(512);



    auto frames = audioProcesser.generateFrames(samples, window, 512, 256);
    auto noiseFrames = audioProcesser.generateFrames(Airnoise, window, 512, 256);

    size_t padded_size = frames[0].size();
    std::vector<double> output(samples.size(), 0.0);

    std::vector<double> noise = NoiseSuppression().generateAverageNoise(noiseFrames, 150, 512);

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

    audioProcesser.saveWav(R"(Output.wav)", output, sampleRate);




}

void PredictionAudio()
{
    // auto model = Model();
    //
    // Ort::Session s = model.GetSession(L"model.onnx");
    //
    // Ort::AllocatorWithDefaultOptions allocator;
    // auto inputs = s.GetInputNames();
    // auto outputs = s.GetOutputNames();
    //
    // const char* input_name = inputs[0].c_str();
    // const char* output_name = outputs[0].c_str();
    //
    // std::vector<float> input(1*257*251*1, 1.0);
    // input.at(1) = 2.0;
    //
    // std::array<int64_t, 4> shape = {1, 257, 251, 1};
    //
    // Ort::MemoryInfo memory_info = Ort::MemoryInfo("Cpu", OrtAllocatorType::OrtArenaAllocator, 0, OrtMemTypeDefault);
    //Ort::MemoryInfo memory_info_cuda("Cuda", OrtArenaAllocator, /*device_id*/0,
                                     //OrtMemTypeDefault);
    // auto type_info = s.GetInputTypeInfo(0);
    // auto tensor_info = type_info.GetTensorTypeAndShapeInfo();
    // auto expected_shape = tensor_info.GetShape();
    //
    // std::cout << "Model expects input shape: ";
    // for (auto dim : expected_shape) std::cout << dim << " ";
    // std::cout << std::endl;
    //
    // auto type = s.GetInputTypeInfo(0)
    //             .GetTensorTypeAndShapeInfo()
    //             .GetElementType();
    // std::cout << "Model expects type: " << type << std::endl;
    // auto input_tensor = Ort::Value::CreateTensor<float>(
    //     memory_info,
    //     input.data(),
    //     input.size(),
    //     shape.data(),
    //     shape.size());
    //
    // auto output_tensor = s.Run(Ort::RunOptions{nullptr},
    //     &input_name,
    //     &input_tensor,
    //     1,
    //     &output_name,
    //     1);
    //
    // float* output = output_tensor[0].GetTensorMutableData<float>();
    //
    // auto& output_value = output_tensor[0];
    // auto shape_info = output_value.GetTensorTypeAndShapeInfo();
    // std::vector<int64_t> output_shape = shape_info.GetShape();
    //
    // //Total elements
    // size_t total_length = 1;
    // for (auto dim: shape)
    //     total_length *= dim;
    //
    // float* float_array = output_value.GetTensorMutableData<float>();
    //
    //
    // for (size_t i = 0; i < total_length; i++)
    // {
    //     std::cout << float_array[i] << std::endl;
    //     if ((i+1) % output_shape.back() == 0) std::cout << "\n";
    // }
    //



}