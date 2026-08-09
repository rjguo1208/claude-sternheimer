#!/usr/bin/env python3
"""Electron-defect spectral function on Gamma-M-K-Gamma.

  A(k,w) = -(1/pi) Im Tr [ w + i.eta - H_W(k) - n_d T(k,w) ]^{-1}

T(k,w) comes from the real-space Koster-Slater cluster built on the omega-resolved
downfolded vertex, and is independent of n_d -- so the same cache serves any
concentration.  No fit in omega: the continued fraction is evaluated at every
point, so the dressed rest state near +2.7 eV is carried exactly."""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import PowerNorm
RY = 13.605693122994; GATE_VBM = -5.9359
z = np.load("/home/rjguo/edi_tmatrix/sternA/vsrlx/kpath_T.npz")
E, T, Hk, ekp, klab = z["E"], z["T"], z["Hk"], z["ekp"], z["klab"]
eta = float(z["eta"]); nw = Hk.shape[1]; nkp = Hk.shape[0]
print("cache: %d omega x %d k, RCUT=%d, fine grid %d, eta=%.0f meV"
      % (len(E), nkp, int(z["rcut"]), int(z["nfine"]), eta*RY*1e3))

def spectral(nd):
    A = np.zeros((len(E), nkp))
    I = np.eye(nw)
    for i, e in enumerate(E):
        zz = (e + GATE_VBM)/RY
        M = (zz + 1j*eta)*I[None] - Hk - nd*T[i]
        A[i] = -np.trace(np.linalg.inv(M), axis1=1, axis2=2).imag/np.pi
    return A

NDS = [(1/144, "$n_d=1/144$  (this calculation's dilution)"),
       (1/36,  "$n_d=1/36$  (the 6$\\times$6 supercell)")]
fig, ax = plt.subplots(1, 2, figsize=(11.6, 4.6), sharey=True)
x = np.arange(nkp)
vmax = None
for p, (nd, lab) in enumerate(NDS):
    A = spectral(nd)
    if vmax is None: vmax = np.percentile(A, 99.5)
    im = ax[p].pcolormesh(x, E, A, cmap="magma", norm=PowerNorm(0.45, vmin=0, vmax=vmax),
                          shading="gouraud", rasterized=True)
    for b in range(nw):
        ax[p].plot(x, ekp[:, b], color="#5fb3ff", lw=.6, alpha=.55)
    for xl in klab[1:-1]: ax[p].axvline(xl, color="w", lw=.6, alpha=.5)
    ax[p].set_xticks(klab); ax[p].set_xticklabels(["$\\Gamma$", "M", "K", "$\\Gamma$"])
    ax[p].set_xlim(0, nkp-1); ax[p].set_ylim(E.min(), E.max())
    ax[p].set_title("(%s)  %s" % ("ab"[p], lab), fontsize=10, loc="left")
    ax[p].tick_params(labelsize=9)
ax[0].set_ylabel("$E-E_\\mathrm{VBM}$   (eV)")
cb = fig.colorbar(im, ax=ax, pad=.015, fraction=.035)
cb.set_label("$A(k,\\omega)$   (states / eV)", fontsize=9); cb.ax.tick_params(labelsize=8)
out = "/home/rjguo/edi_tmatrix/sternA/claude-sternheimer/docs/assets/kpath_spectral.png"
fig.savefig(out, dpi=170, bbox_inches="tight"); print("wrote", out)
for nd, lab in NDS:
    A = spectral(nd)
    g = (E > 0.10) & (E < 1.60)
    print("  n_d=%.4f : in-gap weight %.3f ; max A in gap %.2f at E=%+.3f eV"
          % (nd, np.trapezoid(A[g].sum(1), E[g]) if hasattr(np, "trapezoid") else 0,
             A[g].max(), E[g][A[g].max(axis=1).argmax()]))
