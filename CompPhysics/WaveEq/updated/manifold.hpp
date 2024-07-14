#pragma once
// std library headers
#include <cmath>
#include <string>
#include <utility>
// project headers
#include "matrix.hpp"   // uses std::array as container primative
// 3rd-party libraries
#include <matplot/matplot.h>

namespace manifold
{
    constexpr double rho {1.0};                                     // density of manifold medium
    constexpr double speed {343.0};                                 // speed of wave through medium (speed of sound in air)
    constexpr double c {speed / rho};                               // speed - density ratio
    constexpr double radius {1.5};                                  // radius of manifold

    template <typename UNIT>
    UNIT wave_eq(const size_t t,                                    // time step index
                 const size_t i, const size_t j,                    // spacial step indices
                 const UNIT dx, const UNIT dy)                      // step size through space
    {
        return exp(-((i * dx) * (i * dx) + (j * dy) * (j * dy)) / (radius * radius)) - exp(-1);
    }

    // initialize position matrix for manifold
    template <typename UNIT, size_t NX, size_t NY, size_t NT>
    void init_pos(const std::unique_ptr<matrix::Matrix3<UNIT,NX,NY,NT>>& pos,// position matrix
                  const UNIT dx, const UNIT dy)                     // step size through space
    {
        int cx { (NX - 1) / 2 };                                    // center of manifold circle
        int cy { (NY - 1) / 2 };

        size_t i1{}, i2{}, j1{}, j2{};                              // index elements of our quadrants around our center
        double r{};                                                 // Use a comparison radius to initialize a circular manifold
        // initialize our z matrix
        for(int i = 0; i <= (cx/2); ++i) {
            for(int j = 0; j <= (cy/2); ++j) {
                i1 = i + cx; i2 = -i + cx;                          // index offset from the center in +/- x dir (+/- row index)
                j1 = j + cy; j2 = -j + cy;                          // index offset from the center in +/- y dir (+/- col index)
                r = sqrt((i*dx)*(i*dx) + (j*dy)*(j*dy));            // calculate distance from the center
                if(radius >= r) {   // if the point is within our desired circle size
                    (*pos)[0][i1][j1] = wave_eq<UNIT>(0,i,j,dx,dy); // top right quad element
                    (*pos)[0][i1][j2] = wave_eq<UNIT>(0,i,j,dx,dy); // top left quad element
                    (*pos)[0][i2][j1] = wave_eq<UNIT>(0,i,j,dx,dy); // bottom right quad element
                    (*pos)[0][i2][j2] = wave_eq<UNIT>(0,i,j,dx,dy); // bottom left quad element
                }
            }
        }
    }

    // simply dampen scaling operation
    template <typename UNIT>
    UNIT damp_func(UNIT a) { return 0.9*a; }

    // calculate acceleration for position1
    // accerlation depends on spatial coordinate, not time
    template <typename UNIT>
    UNIT accel_func(const std::pair<UNIT,UNIT> pos2, 
                    const std::pair<UNIT,UNIT> pos1, 
                    const std::pair<UNIT,UNIT> pos0, 
                    const UNIT dx, const UNIT dy) 
    {
        UNIT a{}, dx2{}, dy2{};
        dx2 = (pos2.first - 2*pos1.first + pos0.first) / (dx*dx);
        dy2 = (pos2.second - 2*pos1.second + pos0.second) / (dy*dy);
        a = c*c*(dx2 + dy2);
        // a = damp_func(a);
        return a;
    }

    // handles stepping for velocity and position since they depend on time
    template <typename UNIT>
    UNIT ode_func(const UNIT term1, const UNIT term2, const UNIT dt)
    {
        return term1 + term2 * dt;
    }

    template <typename UNIT, size_t NX, size_t NY, size_t NT>
    void plot(const std::unique_ptr<matrix::Matrix3<UNIT,NX,NY,NT>>& z_matrix, 
              const double dx, const double dy, 
              const std::pair<double,double> length, const size_t t_step)
    {
        // initialize our mesh grid matrices, matplot++ restricts us to doubles
        auto [X,Y] = matplot::meshgrid(matplot::iota(0.0,dx,(length.first+dx)), 
                                       matplot::iota(0.0,dy,(length.second+dy)));
        auto Z = matplot::transform(X, Y, [](double x, double y){ return 0.0; });
        // std::cout << X[0].size() << ", " << Y[0].size() << ", " << Z[0].size() << ", " << (*z_matrix)[t_step].size() << std::endl;
        for (size_t i = 0; i < NX; ++i) {
            for (size_t j = 0; j < NY; ++j) {
                // matplot++ has the reverse order indices I am using
                Z[j][i] = (*z_matrix)[t_step][i][j];
            }
        }
        // std::cout << X.size() << ", " << Y.size() << ", " << Z.size() << std::endl;
        // build meshgrid surface plot
        matplot::mesh(X, Y, Z);
        matplot::zlim({-0.6,0.6});
        // save graph
        matplot::save("img/time"+std::to_string(t_step)+".gif");
        // show last plot configured
        // matplot::show();
    }
} // namespace manifold