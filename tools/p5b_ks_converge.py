#!/usr/bin/env python3
"""Koster-Slater truncation convergence: ||T(Rcut)|| vs the truncation radius,
with the mismatched-gauge filukk (old) vs the gauge-consistent filukk_150 (new).
Values are the verified outputs of tools/tmatrix_p5b.py on the two wann_data files
(also tabulated on the results page)."""
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt

Rcut = np.array([0, 1, 2, 3, 4, 5, 6])
dim  = np.array([11, 99, 275, 539, 891, 1331, 1584])
old  = np.array([0.006, 0.043, 0.097, 0.189, 0.594, 2.061, 2.063])   # 17-band-run filukk (gauge mismatch)
new  = np.array([0.003, 0.007, 0.057, 1.717, 2.059, 2.064, 2.059])   # filukk_150 (consistent gauge)
full = 2.06

fig, ax = plt.subplots(figsize=(6.6, 4.6))
ax.axhline(full, ls=":", c="gray", lw=1, label="full block ($R_{\\rm cut}=6$, dim 1584)")
ax.plot(Rcut, old, "s--", c="#d62728", ms=6, label="old filukk (17-band run, gauge mismatch)")
ax.plot(Rcut, new, "o-",  c="#1f77b4", ms=6, label="filukk_150 (re-Wannierized, consistent gauge)")
ax.axvline(4, ls="--", c="#1f77b4", lw=0.8, alpha=0.5)
ax.annotate("converged at $R_{\\rm cut}=4$\n(dim 891, 56% of full)", xy=(4, 2.059),
            xytext=(1.4, 1.55), fontsize=9, color="#1f77b4",
            arrowprops=dict(arrowstyle="->", color="#1f77b4"))
for r, d in zip(Rcut, dim):
    ax.annotate(f"{d}", (r, -0.02), ha="center", va="top", fontsize=7, color="0.4")
ax.set_xlabel("Koster-Slater cutoff $R_{\\rm cut}$  (subspace dim below each tick)")
ax.set_ylabel(r"$\|T(R_{\rm cut})\|_F$  (Ry, Wannier convention)")
ax.set_title("Koster-Slater truncation converges once the gauge is consistent")
ax.set_ylim(-0.15, 2.3); ax.set_xlim(-0.3, 6.3); ax.legend(fontsize=8, loc="center right"); ax.grid(alpha=0.3)
plt.tight_layout(); plt.savefig("p5b_ks_converge.png", dpi=120); print("wrote p5b_ks_converge.png")
