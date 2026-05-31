//
// Created by prana on 01-08-2025.
//

#include <librosa/eigen3/Eigen/Dense>
#include <onnxruntime_cxx_api.h>

#include <string>
#include <vector>


#ifndef MODEL_H
#define MODEL_H



class Model{
    public:
    Model(const ORTCHAR_T* model_name, int threads=1);

    [[nodiscard]] std::vector<std::string> getProviders() const;

    Ort::Session GetSession(const ORTCHAR_T* model_name) const;

    void fast_predict(float* data, int size, Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& output);


    private:
    Ort::Env env;
    Ort::RunOptions runOptions;
    Ort::SessionOptions sessionOptions;

    Ort::Session session = Ort::Session(nullptr);

    std::vector<std::string> inputs;
    std::vector<std::string> outputs;

    const char* input_name;
    const char* output_name;
    Ort::MemoryInfo memory_info;
    std::array<int64_t, 3> shape = {B, F, T};

    int B = 1;
    int F = 257;
    int T = 5;

    std::vector<float> output_buffer;
    Ort::Value output_tensor{nullptr};
    Ort::Value input_tensor{nullptr};



};



#endif //MODEL_H
