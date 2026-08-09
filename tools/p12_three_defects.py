#!/usr/bin/env python3
"""Electron-defect spectral function of three point defects in monolayer MoS2,
from the same pipeline: omega-resolved downfolded vertex -> pair-Wannier ->
real-space Koster-Slater cluster -> Sigma^ed = n_d T_kk -> A(k,omega).

The vacancy binds two in-gap states; neither isovalent substitution binds any --
max|T| never approaches the singularity of [1 - V G^A]."""
import sys, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import PowerNorm
RY = 13.605693122994; GATE_VBM = -5.9359
B = "/home/rjguo/edi_tmatrix/sternA/vsrlx"
ND = float(sys.argv[1]) if len(sys.argv) > 1 else 1/36
ZOOM = (float(sys.argv[2]), float(sys.argv[3])) if len(sys.argv) > 3 else (-0.35, 2.05)
TAG = sys.argv[4] if len(sys.argv) > 4 else ""
CASES = [("kpath_T450.npz", "V$_S$ vacancy"), ("kpath_T_OS.npz", "O$_S$ substitution"),
         ("kpath_T_SES.npz", "Se$_S$ substitution")]

fig, ax = plt.subplots(1, 3, figsize=(13.4, 4.4), sharey=True)
vmax = None
for p, (fn, lab) in enumerate(CASES):
    z = np.load(f"{B}/{fn}")
    E, T, Hk, ekp, klab = z["E"], z["T"], z["Hk"], z["ekp"], z["klab"]
    eta = float(z["eta"]); nw = Hk.shape[1]; nkp = Hk.shape[0]; I = np.eye(nw)
    A = np.zeros((len(E), nkp))
    for i, e in enumerate(E):
        M = ((e + GATE_VBM)/RY + 1j*eta)*I[None] - Hk - ND*T[i]
        A[i] = -np.trace(np.linalg.inv(M), axis1=1, axis2=2).imag/np.pi
    m = (E >= ZOOM[0]) & (E <= ZOOM[1])
    if vmax is None: vmax = np.percentile(A[m], 99.5)
    x = np.arange(nkp)
    im = ax[p].pcolormesh(x, E, A, cmap="magma", norm=PowerNorm(0.45, vmin=0, vmax=vmax),
                          shading="gouraud", rasterized=True)
    for b in range(nw): ax[p].plot(x, ekp[:, b], color="#5fb3ff", lw=.6, alpha=.5)
    for xl in klab[1:-1]: ax[p].axvline(xl, color="w", lw=.6, alpha=.45)
    ax[p].set_xticks(klab); ax[p].set_xticklabels(["$\\Gamma$", "M", "K", "$\\Gamma$"])
    ax[p].set_xlim(0, nkp-1); ax[p].set_ylim(*ZOOM); ax[p].tick_params(labelsize=9)
    tk = np.abs(T).max(axis=(1, 2, 3))
    ax[p].set_title("(%s)  %s     $\\max|T|=%.2f$" % ("abc"[p], lab, tk.max()),
                    fontsize=10, loc="left")
ax[0].set_ylabel("$E-E_\\mathrm{VBM}$   (eV)")
cb = fig.colorbar(im, ax=ax, pad=.012, fraction=.028)
cb.set_label("$A(k,\\omega)$   (states / eV)", fontsize=9); cb.ax.tick_params(labelsize=8)
fig.suptitle("$n_d = %s$" % ("1/36" if abs(ND-1/36) < 1e-9 else "1/144" if abs(ND-1/144) < 1e-9
                             else "%.4f" % ND), fontsize=10, x=.13, y=.98)
out = "/home/rjguo/edi_tmatrix/sternA/claude-sternheimer/docs/assets/three_defects%s.png" % TAG
fig.savefig(out, dpi=165, bbox_inches="tight"); print("wrote", out)
for fn, lab in CASES:
    z = np.load(f"{B}/{fn}"); tk = np.abs(z["T"]).max(axis=(1, 2, 3))
    print("  %-22s max|T| %7.2f   (a bound state needs [1-V G] singular)" % (lab, tk.max()))
