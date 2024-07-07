#pragma once
// std library headers
#include <memory>
// project headers
#include "matrix.hpp"   // uses std::array as container primative
#include "manifold.hpp" // contains global variables for sim

namespace rk4
{
    template <typename UNIT, size_t NX, size_t NY, size_t NT>
    void runge_kutta_4(matrix::Matrix3<UNIT,NX,NY,NT>& pos,     // position matrix[time][x][y]
                       matrix::Matrix3<UNIT,NX,NY,NT>& vel,     // velocity matrix[time][x][y]
                       const size_t t, const UNIT dt,           // time index & size of step
                       const UNIT dx, const UNIT dy)            // size of spacial step
    {
        // Initialize our accerlation sub-components for rk4 loop
        auto a0_ptr = std::make_unique<matrix::Matrix2<UNIT,NX,NY>>();
        auto a1_ptr = std::make_unique<matrix::Matrix2<UNIT,NX,NY>>();
        auto a2_ptr = std::make_unique<matrix::Matrix2<UNIT,NX,NY>>();
        auto a3_ptr = std::make_unique<matrix::Matrix2<UNIT,NX,NY>>();
        auto a4_ptr = std::make_unique<matrix::Matrix2<UNIT,NX,NY>>();
        // Initialize our k-value velocity sub-components for rk4 loop
        auto k1_ptr = std::make_unique<matrix::Matrix2<UNIT,NX,NY>>();
        auto k2_ptr = std::make_unique<matrix::Matrix2<UNIT,NX,NY>>();
        auto k3_ptr = std::make_unique<matrix::Matrix2<UNIT,NX,NY>>();
        auto k4_ptr = std::make_unique<matrix::Matrix2<UNIT,NX,NY>>();

        /*------- RK4 Term 1 -------*/
        manifold::accel(t,pos,(*a0_ptr),dx,dy);
        for(int i = 0; i < NX - 1; ++i) {
            for(int j = 0; j < NY - 1; ++j) {
                vel[t+1][i][j] = vel[t][i][j] + (*a0_ptr)[i][j] * dt;
                (*k1_ptr)[i][j] = vel[t][i][j] + (*a0_ptr)[i][j] * dt;
                pos[t+1][i][j] = pos[t][i][j] + (*k1_ptr)[i][j] * dt;
            }
        }
        /*------- RK4 Term 2 -------*/
        manifold::accel(t,pos,(*a1_ptr),dx,dy);
        for(int i = 1; i < NX - 1; ++i) {
            for(int j = 1; j < NY - 1; ++j) {
                (*k2_ptr)[i][j] = vel[t+1][i][j] + (*a1_ptr)[i][j] * (dt / 2.0);
                pos[t+1][i][j] = pos[t][i][j] + (*k2_ptr)[i][j] * (dt / 2.0);
            }
        }
        /*------- RK4 Term 3 -------*/
        manifold::accel(t,pos,(*a2_ptr),dx,dy);
        for(int i = 1; i < NX - 1; ++i) {
            for(int j = 1; j < NY - 1; ++j) {
                (*k3_ptr)[i][j] = vel[t+1][i][j] + (*a2_ptr)[i][j] * (dt / 2.0);
                pos[t+1][i][j] = pos[t][i][j] + (*k3_ptr)[i][j] * (dt / 2.0);
            }
        }
        /*------- RK4 Term 4 -------*/
        manifold::accel(t,pos,(*a3_ptr),dx,dy);
        for(int i = 1; i < NX - 1; ++i) {
            for(int j = 1; j < NY - 1; ++j) {
                (*k4_ptr)[i][j] = vel[t+1][i][j] + (*a3_ptr)[i][j] * dt;
                pos[t+1][i][j] = pos[t][i][j] + (*k4_ptr)[i][j] * dt;
            }
        }
        /*------- RK4 Term Sum  -------*/
        manifold::accel(t,pos,(*a4_ptr),dx,dy);
        for(int i = 1; i < NX - 1; ++i) {
            for(int j = 1; j < NY - 1; ++j) {
                vel[t+1][i][j] = vel[t][i][j] + 
                                       (1.0/6.0) * ((*a1_ptr)[i][j] + 
                                       2 * (*a2_ptr)[i][j] + 
                                       2 * (*a3_ptr)[i][j] + 
                                       (*a4_ptr)[i][j]) * dt;

                pos[t+1][i][j] = pos[t][i][j] + 
                                      1.0/6.0 * ((*k1_ptr)[i][j] + 
                                      2 * (*k2_ptr)[i][j] +
                                      2 * (*k3_ptr)[i][j] + 
                                      (*k4_ptr)[i][j]) * dt;
            }
        }
    }    
} // namespace rk4
