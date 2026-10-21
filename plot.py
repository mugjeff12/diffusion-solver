import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

d = np.genfromtxt("results/convergence.csv", delimiter=",", names=True)
order = np.argsort(d["N"])
N, err = d["N"][order], d["k_error"][order]

slope, intercept = np.polyfit(np.log10(N), np.log10(err), 1)
print(f"fitted slope = {slope:.4f}")

plt.figure()
plt.loglog(N, err, "o-", label="numerical")
plt.loglog(N, 10**intercept * N**slope, "--", label=f"fit, slope = {slope:.3f}")
plt.xlabel("N (cells)"); plt.ylabel("|k_eff error|")
plt.title("Convergence of k_eff"); plt.legend(); plt.grid(True, which="both")
plt.savefig("results/convergence.png", dpi=150)

f = np.genfromtxt("results/flux_N50.csv", delimiter=",", names=True)
plt.figure()
plt.plot(f["x"], f["analytic"], "-", label="analytic cos(pi x/a)")
plt.plot(f["x"], f["numerical"], "o", markersize=4, label="numerical (N=50)")
plt.xlabel("x (cm)"); plt.ylabel("normalized flux"); plt.legend(); plt.grid(True)
plt.savefig("results/flux_overlay.png", dpi=150)