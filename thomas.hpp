#pragma once
#include <cstddef>
#include <vector>

// Solve A x = rhs for a tridiagonal matrix A.
// lower[i] = A(i,i-1)  (lower[0] unused)
// diag[i]  = A(i,i)
// upper[i] = A(i,i+1)  (upper[n-1] unused)
// diag and rhs are taken by value because the algorithm overwrites them.
std::vector<double> thomas_solve(const std::vector<double>& lower,
                                 const std::vector<double>& upper,
                                 std::vector<double> diag,
                                 std::vector<double> rhs);