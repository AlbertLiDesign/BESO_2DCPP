cmake_minimum_required(VERSION 3.21)
project(BESO)
find_package(Eigen3 3.4 REQUIRED)

INCLUDE_DIRECTORIES("/opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3/")

file(GLOB SRC_FILES
        "${PROJECT_SOURCE_DIR} /src/*.h"
        "${PROJECT_SOURCE_DIR} /src/*.cpp"
        )

set(CMAKE_CXX_STANDARD 17)

add_executable(${CMAKE_PROJECT_NAME} src/main.cpp src/BESO.cpp include/BESO.h)

