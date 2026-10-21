#include <cmath>
#include <iostream>
#include <vector>
#include "thomas.hpp"

// True if every entry of got matches expected within tol.
bool all_close(const std::vector<double>& got,
               const std::vector<double>& expected, double tol)
{
    if (got.size() != expected.size()) return false;
    for (std::size_t i = 0; i < got.size(); ++i) {
        if (std::fabs(got[i] - expected[i]) > tol) return false;
    }
    return true;
}

int main()
{
    int failures = 0;

    // Test 1: the 3x3 system from the hand calculation. Answer: (1,1,1).
    {
        std::vector<double> lower = {0.0, -1.0, -1.0};
        std::vector<double> upper = {-1.0, -1.0, 0.0};
        std::vector<double> diag  = {2.0, 2.0, 2.0};
        std::vector<double> rhs   = {1.0, 0.0, 1.0};
        std::vector<double> expected = {1.0, 1.0, 1.0};

        std::vector<double> x = thomas_solve(lower, upper, diag, rhs);
        bool ok = all_close(x, expected, 1e-12);
        std::cout << "Test 1 (3x3 hand check): " << (ok ? "PASS" : "FAIL") << "\n";
        if (!ok) ++failures;
    }

    // Test 2: the N=4 diffusion matrix. rhs = A * (1,2,3,4).
    {
        const double c = 0.0016;
        std::vector<double> lower = {0.0, -c, -c, -c};
        std::vector<double> upper = {-c, -c, -c, 0.0};
        std::vector<double> diag  = {0.0748, 0.0732, 0.0732, 0.0748};
        std::vector<double> rhs   = {0.0716, 0.14, 0.21, 0.2944};
        std::vector<double> expected = {1.0, 2.0, 3.0, 4.0};

        std::vector<double> x = thomas_solve(lower, upper, diag, rhs);
        bool ok = all_close(x, expected, 1e-12);
        std::cout << "Test 2 (N=4 diffusion matrix): " << (ok ? "PASS" : "FAIL") << "\n";
        std::cout << "  solution: ";
        for (double value : x) std::cout << value << " ";
        std::cout << "\n";
        if (!ok) ++failures;
    }

    return failures;   // 0 = all passed (this becomes the exit code)
}