#pragma once
// std library headers
#include <memory>
// project headers
#include "matrix.hpp"   // uses std::array as container primative
#include "manifold.hpp" // contains global variables for sim

namespace rk4
{
    template <typename UNIT, size_t NX, size_t NY, size_t NT>
    void runge_kutta_4(const std::unique_ptr<matrix::Matrix3<UNIT,NX,NY,NT>>& pos,  // position matrix[time][x][y]
                       const std::unique_ptr<matrix::Matrix3<UNIT,NX,NY,NT>>& vel,  // velocity matrix[time][x][y]
                       const size_t t, const UNIT dt,   // time index & size of step
                       const UNIT dx, const UNIT dy)    // size of spacial step
    {
        // Initialize our accerlation sub-components for rk4 loop
        auto a0_ptr = std::make_unique<matrix::Matrix2<UNIT,NX,NY>>();
        auto a1_ptr = std::make_unique<matrix::Matrix2<UNIT,NX,NY>>();
        auto a2_ptr = std::make_unique<matrix::Matrix2<UNIT,NX,NY>>();
        auto a3_ptr = std::make_unique<matrix::Matrix2<UNIT,NX,NY>>();
        auto a4_ptr = std::make_unique<matrix::Matrix2<UNIT,NX,NY>>();
        // Initialize our k-value velocity sub-components for rk4 loop
        auto vk1_ptr = std::make_unique<matrix::Matrix2<UNIT,NX,NY>>();
        auto vk2_ptr = std::make_unique<matrix::Matrix2<UNIT,NX,NY>>();
        auto vk3_ptr = std::make_unique<matrix::Matrix2<UNIT,NX,NY>>();
        auto vk4_ptr = std::make_unique<matrix::Matrix2<UNIT,NX,NY>>();
        // Initialize our k-value position sub-components for rk4 loop
        auto pk1_ptr = std::make_unique<matrix::Matrix2<UNIT,NX,NY>>();
        auto pk2_ptr = std::make_unique<matrix::Matrix2<UNIT,NX,NY>>();
        auto pk3_ptr = std::make_unique<matrix::Matrix2<UNIT,NX,NY>>();
        auto pk4_ptr = std::make_unique<matrix::Matrix2<UNIT,NX,NY>>();

        /*------- RK4 Term 1 -------*/
        manifold::accel(t,pos,a0_ptr,dx,dy);
        for(int i = 0; i < NX - 1; ++i) {
            for(int j = 0; j < NY - 1; ++j) {
                (*vk1_ptr)[i][j] = (*vel)[t][i][j] + (*a0_ptr)[i][j] * dt;
                (*pk1_ptr)[i][j] = (*pos)[t][i][j] + (*vel)[t][i][j] * dt;
            }
        }
        /*------- RK4 Term 2 -------*/
        manifold::accel(pk1_ptr,a1_ptr,dx,dy);
        for(int i = 1; i < NX - 1; ++i) {
            for(int j = 1; j < NY - 1; ++j) {
                (*vk2_ptr)[i][j] = ((*vel)[t][i][j] + (0.5*(*vk1_ptr)[i][j])) + (*a1_ptr)[i][j] * (dt / 2.0);
                (*pk2_ptr)[i][j] = ((*pos)[t][i][j] + (0.5*(*pk1_ptr)[i][j])) + (*vk2_ptr)[i][j] * (dt / 2.0);
            }
        }
        /*------- RK4 Term 3 -------*/
        manifold::accel(pk2_ptr,a2_ptr,dx,dy);
        for(int i = 1; i < NX - 1; ++i) {
            for(int j = 1; j < NY - 1; ++j) {
                (*vk3_ptr)[i][j] = ((*vel)[t][i][j] + (0.5*(*vk2_ptr)[i][j])) + (*a2_ptr)[i][j] * (dt / 2.0);
                (*pk3_ptr)[i][j] = ((*pos)[t][i][j] + (0.5*(*pk2_ptr)[i][j])) + (*vk3_ptr)[i][j] * (dt / 2.0);
            }
        }
        /*------- RK4 Term 4 -------*/
        manifold::accel(pk3_ptr,a3_ptr,dx,dy);
        for(int i = 1; i < NX - 1; ++i) {
            for(int j = 1; j < NY - 1; ++j) {
                (*vk4_ptr)[i][j] = ((*vel)[t][i][j] + (*vk3_ptr)[i][j]) + (*a3_ptr)[i][j] * dt;
                (*pk4_ptr)[i][j] = ((*pos)[t][i][j] +(*pk3_ptr)[i][j]) + (*vk4_ptr)[i][j] * dt;
            }
        }
        /*------- RK4 Term Sum  -------*/
        manifold::accel(pk4_ptr,a4_ptr,dx,dy);
        for(int i = 1; i < NX - 1; ++i) {
            for(int j = 1; j < NY - 1; ++j) {
                // calculate velocity for next time step
                (*vel)[t+1][i][j] = (*vel)[t][i][j] + (dt/6.0) * 
                                    ((*a1_ptr)[i][j] + 2 * (*a2_ptr)[i][j] + 2 * (*a3_ptr)[i][j] + (*a4_ptr)[i][j]);
                // calculate position for next time step
                (*pos)[t+1][i][j] = (*pos)[t][i][j] + (dt/6.0) * 
                                    ((*vk1_ptr)[i][j] + 2*(*vk2_ptr)[i][j] + 2*(*vk3_ptr)[i][j] + (*vk4_ptr)[i][j]);
            }
        }
    }    
} // namespace rk4
