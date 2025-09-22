# Specify the minimum version of CMake required
cmake_minimum_required(VERSION 3.10)

# Define the project name
project(DAISY_SEED_CLAPPER_RMS_FIRMWARE CXX)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Define the source files
set(SOURCES src/main_clapper_rms.cpp)

# Define the include directories
include_directories(
    "./src/"
    "../DaisyExamples/libDaisy/"
    "../DaisyExamples/DaisySP/"
)

# Add subdirectories for external libraries
add_subdirectory(../DaisyExamples/libDaisy/ libDaisy_build)
add_subdirectory(../DaisyExamples/DaisySP/ DaisySP_build)

# Add the executable target
add_executable(${PROJECT_NAME} ${SOURCES})

# Link the executable with the libraries
target_link_libraries(${PROJECT_NAME} libDaisy DaisySP)