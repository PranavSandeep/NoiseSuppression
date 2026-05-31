#define MINIAUDIO_IMPLEMENTATION
#include <iostream>
#include <complex>
#include <thread>
#include "AudioInput.h"
#include "AppConfig.hpp"

#define SAMPLE_RATE 0
#define CHANNELS 1
#define FRAMES_PER_BUFF 1024




int main(int argc, char* argv[])
{

    AppConfig cfg = load_config("config.cfg");



    AudioInput::md = new Model(L"model23kgru.onnx", cfg.threads);
    AudioInput::noiseThreshold = cfg.noise_floor;
    AudioInput::floorAttenuation = cfg.attenuation;
    ma_context context;
    ma_context_init(NULL, 0, NULL, &context);

    ma_device_info* pPlaybackInfos;
    ma_uint32 playbackCount;
    ma_device_info* pCaptureInfos;
    ma_uint32 captureCount;

    ma_context_get_devices(
        &context,
        &pPlaybackInfos, &playbackCount,
        &pCaptureInfos,  &captureCount);

    std::cerr << "=== PLAYBACK DEVICES ===\n";
    for (ma_uint32 i = 0; i < playbackCount; i++)
        std::cerr << i << ": " << pPlaybackInfos[i].name << "\n";

    std::cerr << "\n=== CAPTURE DEVICES ===\n";
    for (ma_uint32 i = 0; i < captureCount; i++)
        std::cerr << i << ": " << pCaptureInfos[i].name << "\n";


    int capture_index = -1;
    int playback_index = -1;

    for (int i = 1; i < argc; i++)
    {
        std::string arg = argv[i];
        if (arg == "--capture" && i + 1 < argc) capture_index = std::stoi(argv[++i]);
        else if (arg == "--playback" && i + 1 < argc) playback_index = std::stoi(argv[++i]);
        else if (arg == "--help")
        {
            std::cerr << "Usage: NoiseSuppression.exe [--capture N] [--playback N] [options]\n";
            return 0;
        }
    }

    if (capture_index > (int)captureCount) {std::cerr << "Invalid capture index!\n";}
    if (playback_index > (int)playbackCount) {std::cerr << "Invalid playback index!\n";}






    std::thread audioThread(AudioInput::proccessAudio);

    SetThreadPriority(audioThread.native_handle(), THREAD_PRIORITY_HIGHEST);
    ma_device device;
    ma_device_config device_config;

    device_config = ma_device_config_init(ma_device_type_duplex); //Input and output

    device_config.capture.format = ma_format_f32;
    device_config.capture.channels = CHANNELS;

    device_config.sampleRate = SAMPLE_RATE;

    device_config.noPreSilencedOutputBuffer = false;


    device_config.periodSizeInFrames = 512;

    device_config.playback.channels = CHANNELS;
    device_config.playback.format = ma_format_f32;

    device_config.dataCallback = AudioInput::data_callback;
    device_config.pUserData = nullptr;


    if (capture_index >= 0) device_config.capture.pDeviceID = &pCaptureInfos[capture_index].id;
    if (playback_index >= 0) device_config.playback.pDeviceID = &pPlaybackInfos[playback_index].id;

    if (ma_device_init(nullptr, &device_config, &device) != MA_SUCCESS)
    {
        std::cerr << "Failed to initialize device \n";
        return -1;
    }

    std::cout << "Starting duplex stream \n";
    ma_device_start(&device);

    std::cout << "Press Enter to stop\n";
    std::cin.get();
    ma_context_uninit(&context);
    ma_device_uninit(&device);
    AudioInput::running = false;
    audioThread.join();




    return 0;


}
