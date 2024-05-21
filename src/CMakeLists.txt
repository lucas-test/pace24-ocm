cmake_minimum_required(VERSION 3.10) # Minimum version of CMake required
project(main) # Name of your project

# Set the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

file(GLOB SOURCES "*.cpp" "solver/*.cpp")

# Add the executable with all.cpp files
add_executable(${PROJECT_NAME} ${SOURCES})