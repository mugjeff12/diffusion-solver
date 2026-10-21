#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>
#include "thomas.hpp"

int main()
{
    // Test-case parameters (cm, 1/cm)
    const double D = 1.0;
    const double sigma_a = 0.07;
    const double nu_sigma_f = 0.08;
    const double a = 100.0;
    const std::size_t N = 100;

    const double h = a / N;
    const double c = D / (h * h);

    // Build the tridiagonal loss matrix A (Milestone 2).
    std::vector<double> lower(N, -c), upper(N, -c), diag(N, 2.0 * c + sigma_a);
    diag[0]     = 3.0 * c + sigma_a;   // ghost cell at left wall
    diag[N - 1] = 3.0 * c + sigma_a;   // ghost cell at right wall
    lower[0]    = 0.0;                 // unused
    upper[N - 1] = 0.0;                // unused

    // Initial guess: flat flux, k = 1.
    std::vector<double> phi(N, 1.0);
    double k = 1.0;

    const double tol = 1e-10;
    const int max_iter = 100000;

    std::cout << std::setprecision(12);
    std::cout << "iter  k_eff           dk\n";

    for (int iter = 1; iter <= max_iter; ++iter) {
        // Old fission source.
        double source_old = 0.0;
        for (double value : phi) source_old += nu_sigma_f * value;

        // Right-hand side = old fission source / k.
        std::vector<double> rhs(N);
        for (std::size_t i = 0; i < N; ++i) rhs[i] = nu_sigma_f * phi[i] / k;

        // Solve A phi_new = rhs.
        std::vector<double> phi_new = thomas_solve(lower, upper, diag, rhs);

        // New fission source and k update.
        double source_new = 0.0;
        for (double value : phi_new) source_new += nu_sigma_f * value;
        const double k_new = k * source_new / source_old;

        // Rescale so the flux sums to N (arbitrary scale, fixed for convenience).
        double total = 0.0;
        for (double value : phi_new) total += value;
        for (double& value : phi_new) value *= static_cast<double>(N) / total;

        // Convergence measures.
        double max_dphi = 0.0;
        for (std::size_t i = 0; i < N; ++i) {
            max_dphi = std::fmax(max_dphi, std::fabs(phi_new[i] - phi[i]));
        }
        const double dk = std::fabs(k_new - k);

        std::cout << iter << "  " << k_new << "  " << dk << "\n";

        phi = phi_new;
        k = k_new;

        if (dk < tol && max_dphi < tol) {
            std::cout << "Converged after " << iter << " iterations, k_eff = " << k << "\n";
            return 0;
        }
    }

    std::cout << "Did NOT converge in " << max_iter << " iterations\n";
    return 1;
}