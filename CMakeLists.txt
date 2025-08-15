cmake_minimum_required(VERSION 3.31)
project(NoiseSuppression)

set(CMAKE_CXX_STANDARD 20)

#Set ONNX runtime root
set(ORT_ROOT "C:/Dev/Dudes making out/NoiseSuppression/onnxruntime")
if(NOT DEFINED ORT_ROOT)
    message(FATAL_ERROR "Please set ORT_ROOT: cmake -DORT_ROOT=...</onnxruntime> ..")
endif()

#Include ONNX headers
include_directories("${ORT_ROOT}/include")

#Link libraries
link_directories("${ORT_ROOT}/lib")
message(STATUS ${ORT_ROOT})


add_executable(NoiseSuppression main.cpp
        AudioProcessing.cpp
        AudioProcessing.h
        NoiseSuppression.cpp
        NoiseSuppression.h
        AudioInput.cpp
        AudioInput.h
        Model.cpp
        Model.h
)
target_include_directories(NoiseSuppression PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

#Link ONNX core runtime
target_link_libraries(NoiseSuppression PRIVATE onnxruntime)

#Link CUDA compilier
target_link_libraries(NoiseSuppression PRIVATE onnxruntime_providers_cuda)
