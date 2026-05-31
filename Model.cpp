//
// Created by prana on 01-08-2025.
//

#include "Model.h"

#include <iostream>
#include <thread>
#include <__msvc_ostream.hpp>

Model::Model(const ORTCHAR_T* model_name, int threads)
    : memory_info(Ort::MemoryInfo("Cpu", OrtAllocatorType::OrtArenaAllocator, 0, OrtMemTypeDefault))
{

    sessionOptions.SetIntraOpNumThreads(threads);
    sessionOptions.SetInterOpNumThreads(threads);
    sessionOptions.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_ALL);

    sessionOptions.SetExecutionMode(ExecutionMode::ORT_SEQUENTIAL);


    //Don't really need CUDA here. But just leaving it here incase it's needed.
    // OrtCUDAProviderOptions cuda_options;
    //
    // cuda_options.device_id = 0;
    //
    //
    // sessionOptions.AppendExecutionProvider_CUDA(cuda_options);


    session = Ort::Session(env, model_name, sessionOptions);



    inputs = session.GetInputNames();
    outputs = session.GetOutputNames();

    input_name = inputs[0].c_str();
    output_name = outputs[0].c_str();

    shape = {B, F, T};

    output_buffer.resize(shape[0] * shape[1] * shape[2]);
    output_tensor = Ort::Value::CreateTensor<float>(
            memory_info,
            output_buffer.data(),
    output_buffer.size(),
        shape.data(),
    shape.size());


}

std::vector<std::string> Model::getProviders() const
{
    std::vector<std::string> providers = Ort::GetAvailableProviders();
    return providers;
}

Ort::Session Model::GetSession(const ORTCHAR_T* model_name) const {
    return Ort::Session(env, model_name, sessionOptions);
}




void Model::fast_predict(float *data, int size, Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> &output)
{
    if (size != B*F*T)
    {
        std::cerr << "Input size is not " << B*F*T << "but " << size << "\n";
        return;
    }

    input_tensor = Ort::Value::CreateTensor<float>(
        memory_info,
        data,
        B*F*T,
        shape.data(),
        shape.size());



    auto output_tensors = session.Run(Ort::RunOptions{nullptr},
    &input_name, &input_tensor, 1,
    &output_name, 1);

    float* output_data = output_tensors[0].GetTensorMutableData<float>();

    output = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(
        output_data, F, T);
}


