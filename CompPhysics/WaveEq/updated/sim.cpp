// std library headers
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
// project headers
#include "matrix.hpp"
#include "rk4.hpp"
#include "manifold.hpp"

template <typename UNIT, size_t NX, size_t NY, size_t NT>
void WriteData(matrix::Matrix3<UNIT,NX,NY,NT>& pos, const UNIT dx, const UNIT dy){
    for(int t = 0; t < NT; ++t){
        // std::cout << t << "\n";
        if( t == 1){
            std::stringstream ss;
            ss << t;
            std::string str;
            ss >> str;
            std::string filename = "new_time" + str + ".dat";
            // had to replace ofstream with fopen, ofstream was crashing due to size
            FILE* output = fopen(filename.c_str(),"w+");
    
            for(int i = 0; i < NX; ++i) {
                std::string line{};
                for(int j = 0; j < NY; ++j) {
                    line = std::to_string(i*dx) + " " + std::to_string(j*dy) + " " + std::to_string(pos[t][i][j]) + "\n";
                    fputs(line.c_str(),output);
                }
            }
            fclose(output);
        }
    }
}

int main() 
{
    // Initializing sim values
    constexpr size_t nx { 101 };                        // size of mesh
    constexpr size_t ny { 101 };
    constexpr size_t nsteps { 51 };                     // num of time steps to traverse
    constexpr double dx { 1.0 / (nx - 1.0) };           // mesh step size
    constexpr double dy { 1.0 / (ny - 1.0) };
    constexpr double dt { 1.0e-2 / (nsteps - 1.0) };    // time step size

    // Initialize position & velocity matrix pointers
    auto z_ptr = std::make_unique<matrix::Matrix3<double,nx,ny,nsteps>>();
    auto vz_ptr = std::make_unique<matrix::Matrix3<double,nx,ny,nsteps>>();
    // Initialize position matrix for manifold size of radius R
    manifold::init_pos((*z_ptr),dx,dy);

    // RK4 time step loop 
    for(size_t t = 0; t < (nsteps-1); ++t) {
        std::cout << "time step: " << t << ", processing...";
        rk4::runge_kutta_4((*z_ptr),(*vz_ptr),t,dt,dx,dy);
        std::cout << "...complete.\n";
    }

    WriteData((*z_ptr),dx,dy);
    
    return 0;
}





