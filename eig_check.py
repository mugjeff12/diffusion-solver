import numpy as np

D, sa, nsf, a, N = 1.0, 0.07, 0.08, 100.0, 100
h = a / N
c = D / h**2

A = np.diag(np.full(N, 2*c + sa)) \
  + np.diag(np.full(N-1, -c), 1) + np.diag(np.full(N-1, -c), -1)
A[0, 0] = A[-1, -1] = 3*c + sa

lam = np.linalg.eigvalsh(A)          # sorted ascending (A is symmetric)
k = nsf / lam[0]
print(f"numpy k_eff  = {k:.12f}")