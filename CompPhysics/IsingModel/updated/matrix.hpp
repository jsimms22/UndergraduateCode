#pragma once
// std library headers
#include <fstream>
#include <array>

namespace matrix
{
    // Matrix layout [depth][row][column]
    template <typename UNIT, size_t NX, size_t NY>
    using Matrix2 = std::array<std::array<UNIT,NY>,NX>;

    // Matrix layout [depth][row][column]
    template <typename UNIT, size_t NX, size_t NY, size_t NT>
    using Matrix3 = std::array<Matrix2<UNIT,NX,NY>,NT>;
} // namespace matrix

// operator<< overload for Matrix3
template <typename UNIT, size_t NX, size_t NY, size_t NT>
std::ostream& operator<<(std::ostream& os, const matrix::Matrix3<UNIT,NX,NY,NT>& matrix)
{
    for(size_t t = 0; t < NT; ++t) {
        for(size_t i = 0; i < NX; ++i) {
            for(size_t j = 0; j < NY; ++j) {
                // [time depth][row][column]
                os << matrix[t][i][j] << ' ';
            }
            os << '\n';
        }
        os << '\n';
    }
    return os;
}

// operator<< overload for Matrix2
template <typename UNIT, size_t NX, size_t NY>
std::ostream& operator<<(std::ostream& os, const matrix::Matrix2<UNIT,NX,NY>& matrix)
{
    for(size_t i = 0; i < NX; ++i) {
        for(size_t j = 0; j < NY; ++j) {
            // [row][column]
            os << matrix[i][j] << ' ';
        }
        os << '\n';
    }
    return os;
}





