# CMAKE TERMINAL COMMAND:
# cmake -S <top_level_src_dir> -B <build_dir> -G "MSYS Makefiles"
cmake_minimum_required(VERSION 3.14.0)
project(ising)
# https://github.com/alandefreitas/matplotplusplus/blob/master/README.md#cmake
add_subdirectory(matplot++)

add_executable(${PROJECT_NAME} ising.cpp)
target_link_libraries(${PROJECT_NAME} PUBLIC matplot)