cmake_minimum_required(VERSION 3.31)
project(NoiseSuppression)

set(CMAKE_CXX_STANDARD 20)

add_executable(NoiseSuppression main.cpp
        AudioProcessing.cpp
        AudioProcessing.h
        NoiseSuppression.cpp
        NoiseSuppression.h
)
target_include_directories(NoiseSuppression PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
