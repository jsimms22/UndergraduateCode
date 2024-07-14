#pragma once
// std library headers
#include <functional>
#include <memory>
// project headers
#include "matrix.hpp"   // uses std::array as container primative
#include "manifold.hpp" // contains global variables for sim

namespace rk4
{
    template <typename UNIT, size_t NX, size_t NY, size_t NT>
    using Matrix3 = matrix::Matrix3<UNIT,NX,NY,NT>;

    template <typename UNIT, size_t NX, size_t NY, size_t NT>
    void runge_kutta_4(const std::unique_ptr<Matrix3<UNIT,NX,NY,NT>>& pos,  // position matrix[time][x][y]
                       const std::unique_ptr<Matrix3<UNIT,NX,NY,NT>>& vel,  // velocity matrix[time][x][y]
                       const size_t t, const UNIT dt,                               // time index & size of step
                       const UNIT dx, const UNIT dy)                                // size of spacial step
    {
        // Initialize k-value acceleration sub-components
        auto ak_ptr = std::make_unique<Matrix3<UNIT,NX,NY,5>>();    // a0, a1, a2, a3, a4
        // Initialize k-value velocity sub-components
        auto vk_ptr = std::make_unique<Matrix3<UNIT,NX,NY,4>>();    // v0, v1, v2, v3
        // Initialize k-value position sub-components
        auto pk_ptr = std::make_unique<Matrix3<UNIT,NX,NY,4>>();    // p0, p1, p2, p3

        /*------- RK4 Term 1 -------*/
        for(int i = 1; i < NX - 1; ++i) {
            for(int j = 1; j < NY - 1; ++j) {
                (*ak_ptr)[0][i][j] = manifold::accel_func({(*pos)[t][i+1][j], (*pos)[t][i][j+1]},
                                                          {(*pos)[t][i+0][j], (*pos)[t][i][j+0]},
                                                          {(*pos)[t][i-1][j], (*pos)[t][i][j-1]},
                                                          dx, dy);
                (*vk_ptr)[0][i][j] = manifold::ode_func((*vel)[t][i][j],(*ak_ptr)[0][i][j],dt);
                (*pk_ptr)[0][i][j] = manifold::ode_func((*pos)[t][i][j],(*vel)[t][i][j],dt);
            }
        } /*------- RK4 Term 2 -------*/
        for(int i = 1; i < NX - 1; ++i) {
            for(int j = 1; j < NY - 1; ++j) {
                (*ak_ptr)[1][i][j] = manifold::accel_func({(*pk_ptr)[0][i+1][j], (*pk_ptr)[0][i][j+1]},
                                                          {(*pk_ptr)[0][i+0][j], (*pk_ptr)[0][i][j+0]},
                                                          {(*pk_ptr)[0][i-1][j], (*pk_ptr)[0][i][j-1]},
                                                          dx, dy);
                (*vk_ptr)[1][i][j] = manifold::ode_func(((*vel)[t][i][j] + (0.5*(*vk_ptr)[0][i][j])), (*ak_ptr)[1][i][j], (dt/2.0));
                (*pk_ptr)[1][i][j] = manifold::ode_func(((*pos)[t][i][j] + (0.5*(*pk_ptr)[0][i][j])), (*vk_ptr)[1][i][j], (dt/2.0));
            }
        }/*------- RK4 Term 3 -------*/
        for(int i = 1; i < NX - 1; ++i) {
            for(int j = 1; j < NY - 1; ++j) {
                (*ak_ptr)[2][i][j] = manifold::accel_func({(*pk_ptr)[1][i+1][j], (*pk_ptr)[1][i][j+1]},
                                                          {(*pk_ptr)[1][i+0][j], (*pk_ptr)[1][i][j+0]},
                                                          {(*pk_ptr)[1][i-1][j], (*pk_ptr)[1][i][j-1]},
                                                          dx, dy);
                (*vk_ptr)[2][i][j] = manifold::ode_func(((*vel)[t][i][j] + (0.5*(*vk_ptr)[1][i][j])), (*ak_ptr)[2][i][j], (dt/2.0));
                (*pk_ptr)[2][i][j] = manifold::ode_func(((*pos)[t][i][j] + (0.5*(*pk_ptr)[1][i][j])), (*vk_ptr)[2][i][j], (dt/2.0));
            }
        } /*------- RK4 Term 4 -------*/
        for(int i = 1; i < NX - 1; ++i) {
            for(int j = 1; j < NY - 1; ++j) {
                (*ak_ptr)[3][i][j] = manifold::accel_func({(*pk_ptr)[2][i+1][j], (*pk_ptr)[2][i][j+1]},
                                                          {(*pk_ptr)[2][i+0][j], (*pk_ptr)[2][i][j+0]},
                                                          {(*pk_ptr)[2][i-1][j], (*pk_ptr)[2][i][j-1]},
                                                          dx, dy);
                (*vk_ptr)[3][i][j] = manifold::ode_func(((*vel)[t][i][j] + (*vk_ptr)[2][i][j]), (*ak_ptr)[3][i][j], dt);
                (*pk_ptr)[3][i][j] = manifold::ode_func(((*pos)[t][i][j] + (*pk_ptr)[2][i][j]), (*vk_ptr)[3][i][j], dt);
            }
        }/*------- RK4 Term Sum  -------*/
        for(int i = 1; i < NX - 1; ++i) {
            for(int j = 1; j < NY - 1; ++j) {
                (*ak_ptr)[4][i][j] = manifold::accel_func({(*pk_ptr)[3][i+1][j], (*pk_ptr)[3][i][j+1]},
                                                          {(*pk_ptr)[3][i+0][j], (*pk_ptr)[3][i][j+0]},
                                                          {(*pk_ptr)[3][i-1][j], (*pk_ptr)[3][i][j-1]},
                                                          dx, dy);
                // calculate velocity for next time step
                (*vel)[t+1][i][j] = (*vel)[t][i][j] + (dt/6.0) * 
                                    ((*ak_ptr)[1][i][j] + 2 * (*ak_ptr)[2][i][j] + 2 * (*ak_ptr)[3][i][j] + (*ak_ptr)[4][i][j]);
                // calculate position for next time step
                (*pos)[t+1][i][j] = (*pos)[t][i][j] + (dt/6.0) * 
                                    ((*vk_ptr)[0][i][j] + 2*(*vk_ptr)[1][i][j] + 2*(*vk_ptr)[2][i][j] + (*vk_ptr)[3][i][j]);
            }
        }
    }
} // namespace rk4
