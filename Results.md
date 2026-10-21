# Results (Milestone 5)
- k_analytic = 1.12696753747
- N=100: k_numerical = 1.12696882613, k_error = 1.289e-06, max flux error = 6.58e-12
- N=200: k_numerical = 1.12696785965, k_error = 3.222e-07, max flux error = 6.51e-12
- Error ratio (N=100 / N=200) = 4.000
- Power iteration: 238 iterations at tol 1e-12 (both N)
- NumPy eigvalsh cross-check (N=100): k = 1.126968826131

# Results (Milestone 6)
- Fitted log-log slope of k error vs N (N = 10..1000): -2.000
- Error ratio per doubling of N: 3.99 to 4.00
- k error: N=10 -> 1.285e-04, N=100 -> 1.289e-06, N=1000 -> 1.289e-08
- Max flux error vs cos(pi x/a): ~6.6e-12 (N=200)
- Plots: results/convergence.png, results/flux_overlay.png