//
// Created by prana on 29-07-2025.
//

#include "miniaudio.h"
#include <atomic>
#include "Model.h"



#ifndef AUDIOINPUT_H
#define AUDIOINPUT_H


class AudioInput {
public:
    static void data_callback(ma_device* device, void* output, const void* input,ma_uint32 frame_count);

    static void proccessAudio();

    static std::atomic<bool> running;
    static inline std::atomic<float> noiseThreshold;
    static inline std::atomic<float> floorAttenuation;

    static inline Model* md = nullptr;




};



#endif //AUDIOINPUT_H
