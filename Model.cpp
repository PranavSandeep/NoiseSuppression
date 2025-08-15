//
// Created by prana on 01-08-2025.
//

#include "Model.h"

#include <iostream>
#include <thread>
#include <__msvc_ostream.hpp>

Model::Model(const ORTCHAR_T* model_name)
    : memory_info(Ort::MemoryInfo("Cpu", OrtAllocatorType::OrtArenaAllocator, 0, OrtMemTypeDefault))
{
    sessionOptions.EnableProfiling(L"profiledata.json");  // ✅ must be before session is created
    sessionOptions.SetIntraOpNumThreads(12);
    sessionOptions.SetInterOpNumThreads(12);
    sessionOptions.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_ALL);

    sessionOptions.SetExecutionMode(ExecutionMode::ORT_PARALLEL);

    OrtCUDAProviderOptions cuda_options;

    cuda_options.device_id = 0;


    sessionOptions.AppendExecutionProvider_CUDA(cuda_options);

    // ✅ Create session *after* all options are set
    session = Ort::Session(env, model_name, sessionOptions);

    inputs = session.GetInputNames();
    outputs = session.GetOutputNames();

    input_name = inputs[0].c_str();
    output_name = outputs[0].c_str();

    shape = {1, 257, 251, 1};
}

std::vector<std::string> Model::getProviders() const
{
    std::vector<std::string> providers = Ort::GetAvailableProviders();
    return providers;
}

Ort::Session Model::GetSession(const ORTCHAR_T* model_name) const {
    return Ort::Session(env, model_name, sessionOptions);
}

std::vector<float> Model::predict(std::vector<float>& input)
{

    if (input.size() != 1*257*251*1)
    {
        std::cerr << "Input size is not " << 1*257*251*1 << "but " << input.size() << "\n";
        return {};

    }


    auto input_tensor = Ort::Value::CreateTensor<float>(
        memory_info,
        input.data(),
        input.size(),
        shape.data(),
        shape.size());


    auto output_tensors = session.Run(Ort::RunOptions{nullptr},
    &input_name, &input_tensor, 1,
    &output_name, 1);

    Ort::AllocatorWithDefaultOptions allocator;
    Ort::AllocatedStringPtr profile_path = session.EndProfilingAllocated(allocator);

    // Use it like a normal char*
    std::cout << "Profiling file saved to: " << profile_path.get() << "\n";


    Ort::Value& output_tensor = output_tensors[0];
    float* output_data = output_tensor.GetTensorMutableData<float>();



    constexpr size_t tensor_size = 1*257*251*1;

    std::vector<float> result(output_data, output_data + tensor_size);

    return result;

}

void Model::setSessionOptions(Ort::Session &session)
{


}
