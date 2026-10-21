#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include "thomas.hpp"

int main(int argc, char* argv[])
{
    // Usage: ./diffusion [N]   (default N = 100)
    const std::size_t N = (argc > 1) ? std::strtoul(argv[1], nullptr, 10) : 100UL;

    const double D = 1.0, sigma_a = 0.07, nu_sigma_f = 0.08, a = 100.0;
    const double h = a / static_cast<double>(N);
    const double c = D / (h * h);

    // Loss matrix A (Milestone 2)
    std::vector<double> lower(N, -c), upper(N, -c), diag(N, 2.0 * c + sigma_a);
    diag[0] = 3.0 * c + sigma_a;
    diag[N - 1] = 3.0 * c + sigma_a;
    lower[0] = 0.0;
    upper[N - 1] = 0.0;

    // Power iteration (Milestone 4)
    std::vector<double> phi(N, 1.0);
    double k = 1.0;
    const double tol = 1e-12;
    int iterations = 0;
    bool converged = false;

    for (int iter = 1; iter <= 200000; ++iter) {
        double source_old = 0.0;
        for (double v : phi) source_old += nu_sigma_f * v;

        std::vector<double> rhs(N);
        for (std::size_t i = 0; i < N; ++i) rhs[i] = nu_sigma_f * phi[i] / k;

        std::vector<double> phi_new = thomas_solve(lower, upper, diag, rhs);

        double source_new = 0.0;
        for (double v : phi_new) source_new += nu_sigma_f * v;
        const double k_new = k * source_new / source_old;

        double total = 0.0;
        for (double v : phi_new) total += v;
        for (double& v : phi_new) v *= static_cast<double>(N) / total;

        double max_dphi = 0.0;
        for (std::size_t i = 0; i < N; ++i)
            max_dphi = std::fmax(max_dphi, std::fabs(phi_new[i] - phi[i]));
        const double dk = std::fabs(k_new - k);

        phi = phi_new;
        k = k_new;
        iterations = iter;
        if (dk < tol && max_dphi < tol) { converged = true; break; }
    }

    // Verification (Milestone 5)
    const double pi = std::acos(-1.0);
    const double B2 = (pi / a) * (pi / a);
    const double k_analytic = nu_sigma_f / (sigma_a + D * B2);
    const double k_error = std::fabs(k - k_analytic);

    // Analytic flux at the cell centres, then normalize both to max = 1.
    std::vector<double> x(N), phi_analytic(N);
    for (std::size_t i = 0; i < N; ++i) {
        x[i] = -a / 2.0 + (static_cast<double>(i) + 0.5) * h;   // i is 0-based
        phi_analytic[i] = std::cos(pi * x[i] / a);
    }
    double max_num = 0.0, max_ana = 0.0;
    for (std::size_t i = 0; i < N; ++i) {
        max_num = std::fmax(max_num, phi[i]);
        max_ana = std::fmax(max_ana, phi_analytic[i]);
    }
    double max_flux_error = 0.0;
    std::ofstream out("results/flux.csv");
    out << std::setprecision(12) << "x,numerical,analytic\n";
    for (std::size_t i = 0; i < N; ++i) {
        const double num = phi[i] / max_num;
        const double ana = phi_analytic[i] / max_ana;
        max_flux_error = std::fmax(max_flux_error, std::fabs(num - ana));
        out << x[i] << "," << num << "," << ana << "\n";
    }

    // Convergence study (Milestone 6): append one row per run.
    std::ofstream conv("results/convergence.csv", std::ios::app);
    conv << std::setprecision(12) << N << "," << k << "," << k_analytic << ","
         << k_error << "," << max_flux_error << "," << iterations << "\n";

    std::cout << std::setprecision(12);
    std::cout << "N              = " << N << "\n"
              << "converged      = " << (converged ? "yes" : "NO") << " (" << iterations << " iterations)\n"
              << "k_numerical    = " << k << "\n"
              << "k_analytic     = " << k_analytic << "\n"
              << "k_error        = " << k_error << "\n"
              << "max_flux_error = " << max_flux_error << "\n";
    return converged ? 0 : 1;
}