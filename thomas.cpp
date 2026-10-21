#include "thomas.hpp"

std::vector<double> thomas_solve(const std::vector<double>& lower,
                                 const std::vector<double>& upper,
                                 std::vector<double> diag,
                                 std::vector<double> rhs)
{
    const std::size_t n = diag.size();

    // Forward sweep: eliminate the sub-diagonal, top to bottom.
    for (std::size_t i = 1; i < n; ++i) {
        double w = lower[i] / diag[i - 1];
        diag[i] -= w * upper[i - 1];
        rhs[i]  -= w * rhs[i - 1];
    }

    // Back substitution: solve the last unknown, then walk upward.
    std::vector<double> x(n);
    x[n - 1] = rhs[n - 1] / diag[n - 1];
    for (std::size_t i = n - 1; i-- > 0; ) {
        x[i] = (rhs[i] - upper[i] * x[i + 1]) / diag[i];
    }
    return x;
}