//
// Created by prana on 01-08-2025.
//
#include <onnxruntime_cxx_api.h>
#include <string>
#include <vector>

#ifndef MODEL_H
#define MODEL_H



class Model{
    public:
    Model(const ORTCHAR_T* model_name);

    [[nodiscard]] std::vector<std::string> getProviders() const;

    Ort::Session GetSession(const ORTCHAR_T* model_name) const;

    std::vector<float> predict(std::vector<float>& input);

    void setSessionOptions(Ort::Session& session);


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
    std::array<int64_t, 4> shape = {1, 257, 251, 1};




};



#endif //MODEL_H
