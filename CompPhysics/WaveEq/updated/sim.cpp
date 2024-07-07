// std library headers
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
// project headers
#include "matrix.hpp"
#include "rk4.hpp"
#include "manifold.hpp"

int main() 
{
    // Initializing sim values
    constexpr size_t nx { 51 };                         // size of mesh
    constexpr size_t ny { 51 };
    constexpr size_t n_steps { 4001 };                  // num of time steps to traverse
    constexpr double dx { 1.0 / (nx - 1.0) };           // mesh step size
    constexpr double dy { 1.0 / (ny - 1.0) };
    // time step size needs to be tweaked for various n x n grids for some reason?
    // constexpr double dt { 1e-2 / (n_steps - 1.0) };  // time step size
    // constexpr double dt { 1e-5 / (1.0 - dx*dy) };    
    constexpr double dt { 1e-5 };                       // stable-ish time step for 25x25 to 101-101 mesh
    // std::cout << dt << '\n';

    // Initialize position & velocity matrix pointers
    auto z_ptr = std::make_unique<matrix::Matrix3<double,nx,ny,n_steps>>();
    auto vz_ptr = std::make_unique<matrix::Matrix3<double,nx,ny,n_steps>>();

    // Initialize position matrix for manifold size of radius R
    manifold::init_pos(z_ptr,dx,dy);

    // RK4 time step loop 
    for(size_t t_step = 0; t_step < (n_steps-1); ++t_step) {
        std::cout << "time step: " << t_step << ", processing...";
        rk4::runge_kutta_4((*z_ptr),(*vz_ptr),t_step,dt,dx,dy);
        std::cout << "...complete.\n";
    }

    for(size_t t_step = 0; t_step < n_steps; ++t_step) {
        // plot a meshgrid for every xth time step
        if ((t_step % 10) == 0) { manifold::plot(z_ptr, dx, t_step); }
    }
    return 0;
}







