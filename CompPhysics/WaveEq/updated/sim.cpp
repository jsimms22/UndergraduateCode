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
    constexpr size_t nx {51};                   // size of mesh
    constexpr size_t ny {51};
    constexpr size_t n_steps {1001};            // num of time steps to traverse
    constexpr double x_side {10.0};              // length of our 2d mesh sides
    constexpr double y_side {10.0};
    constexpr double dx {x_side / (nx - 1.0)};  // mesh step size
    constexpr double dy {y_side / (ny - 1.0)};
    constexpr double cfl {0.2};                 // courant-friedrichs-lewy condition
    constexpr double dt {cfl*(dx+dy) / 343.0};  // time step size, where speed = 343 mph
    // std::cout << dt << '\n';

    // Initialize position & velocity matrix pointers
    auto z_ptr = std::make_unique<matrix::Matrix3<double,nx,ny,n_steps>>();
    auto vz_ptr = std::make_unique<matrix::Matrix3<double,nx,ny,n_steps>>();

    // Initialize position matrix for manifold size of radius R
    manifold::init_pos(z_ptr,dx,dy);

    // RK4 time step loop 
    for(size_t t_step = 0; t_step < (n_steps-1); ++t_step) {
        std::cout << "time step: " << t_step << ", processing...";
        rk4::runge_kutta_4(z_ptr,vz_ptr,t_step,dt,dx,dy);
        std::cout << "...complete.\n";
    }

    for(size_t t_step = 0; t_step < n_steps; ++t_step) {
        // plot a meshgrid for every xth time step
        if ((t_step % 5) == 0) { manifold::plot(z_ptr, dx, t_step, x_side); }
    }
    return 0;
}







