//
// Created by prana on 29-07-2025.
//

#include "miniaudio.h"

#ifndef AUDIOINPUT_H
#define AUDIOINPUT_H



class AudioInput {
public:
    static void data_callback(ma_device* device, void* output, const void* input,ma_uint32 frame_count);

};



#endif //AUDIOINPUT_H
