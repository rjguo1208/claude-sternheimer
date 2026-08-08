#!/usr/bin/env python3
"""Vtilde(omega) on the 12x12 grid: spatial decay across the gap, and how well it
interpolates.  omega0 sits below the gap, so freezing Sigma there is the largest
single error in the chain; sampling it at Chebyshev nodes costs 6 s each."""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib import cm
RY = 13.605693122994; NK = 144
z = np.load("/home/rjguo/edi_tmatrix/sternA/vsrlx/vtilde_omega.npz")
E, d12 = z["E"], z["d12"]
fig, ax = plt.subplots(1, 2, figsize=(11.2, 4.3))
cols = cm.viridis(np.linspace(.08, .88, len(E)))

ed = np.arange(-0.5, 7.6, 1.0)
def env(v):
    x, y = [], []
    for lo, hi in zip(ed[:-1], ed[1:]):
        m = (d12 >= lo) & (d12 < hi)
        if m.any(): x.append(max(0.0, .5*(lo+hi))); y.append(v[m].max())
    return np.array(x), np.array(y)
xb, yb = env(z["env_bare"])
ax[0].semilogy(xb, yb, "--", color="#8a8f98", lw=2.2, label="bare $M_{AA}$ ($\\omega$-independent)")
for i, e in enumerate(E):
    x, y = env(z["env_%d" % i])
    ax[0].semilogy(x, y, "-o", color=cols[i], lw=1.5, ms=4, label="$\\omega=%+.2f$ eV" % e)
ax[0].axvline(3.0, color="#c1121f", lw=.9, ls=":")
ax[0].text(3.05, 2e-5, " 6$\\times$6 cell edge", color="#c1121f", fontsize=8)
ax[0].set_xlim(-0.4, 7.2); ax[0].set_ylim(5e-6, 6)
ax[0].set_xlabel("$|R_e - R_\\mathrm{def}|$   (lattice constants)")
ax[0].set_ylabel("$\\max_{R_p}|\\tilde{\\mathcal{V}}(R_e,R_p;\\omega)|$")
ax[0].legend(frameon=False, fontsize=7.6, ncol=2, loc="upper right")
ax[0].set_title("(a)  the dressed vertex stays local across the gap", fontsize=10, loc="left")

med = np.array([1e3*RY*np.median(z["err_%d" % i])/NK for i in range(len(E))])
p90 = np.array([1e3*RY*np.percentile(z["err_%d" % i], 90)/NK for i in range(len(E))])
mx = np.array([1e3*RY*z["err_%d" % i].max()/NK for i in range(len(E))])
ax[1].fill_between(E, med, p90, color="#1f3b73", alpha=.18, lw=0, label="median to p90")
ax[1].plot(E, med, "-o", color="#1f3b73", lw=2.0, ms=5, label="median")
ax[1].plot(E, mx, "--^", color="#1f3b73", lw=1.3, ms=5, alpha=.75, label="worst element")
ax[1].axhline(0.74, color="#8a8f98", lw=1.4, ls="--")
ax[1].text(0.10, 0.78, "bare $M_{AA}$ median", color="#666", fontsize=8.5)
for x, lab, c in ((0.0058, "VBM", "#c1121f"), (1.6679, "CBM", "#c1121f")):
    ax[1].axvline(x, color=c, lw=1.0, ls=":")
    ax[1].text(x, 3.75, lab, color=c, fontsize=8.5, ha="center")
ax[1].set_xlim(-0.06, 1.74); ax[1].set_ylim(0, 4.0)
ax[1].set_xlabel("$\\omega$   (eV above the pristine VBM)")
ax[1].set_ylabel("interpolation error of $\\tilde V(\\omega)$   (meV)")
ax[1].legend(frameon=False, fontsize=8.5, loc="upper left")
ax[1].set_title("(b)  leave-one-out, $\\omega$ by $\\omega$", fontsize=10, loc="left")
for a in ax: a.tick_params(labelsize=9); a.grid(alpha=.18, lw=.6)
fig.tight_layout()
out = "/home/rjguo/edi_tmatrix/sternA/claude-sternheimer/docs/assets/vtilde_omega.png"
fig.savefig(out, dpi=170); print("wrote", out)
print("  held-out node %.4f eV" % z["held"])
print("  loo median  %.3f -> %.3f meV across the gap (bare vertex 0.74)" % (med[0], med[-1]))
print("  on-site |V| %.3f -> %.3f ;  at 6a  %.2e -> %.2e"
      % (z["env_0"][np.argmin(d12)], z["env_%d" % (len(E)-1)][np.argmin(d12)],
         z["env_0"][np.argmax(d12)], z["env_%d" % (len(E)-1)][np.argmax(d12)]))
