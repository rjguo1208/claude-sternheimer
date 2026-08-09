#!/usr/bin/env python3
"""Zoom on the band edges at an experimentally relevant defect density.

n_d is given in cm^-2 and converted with the primitive-cell area
A = (sqrt(3)/2) a^2 = 8.786 A^2, so 1e12 cm^-2 = 8.79e-4 defects per cell
(one per ~1140 cells).  At that dilution the defect contribution is small and a
log colour scale is required to see it -- the same convention as the earlier
1%/5% figures."""
import sys, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
RY = 13.605693122994; GATE_VBM = -5.9359
A_CELL = np.sqrt(3)/2 * 3.18517668**2 * 1e-16          # cm^2 per primitive cell
B = "/home/rjguo/edi_tmatrix/sternA/vsrlx"
NCM = float(sys.argv[1]) if len(sys.argv) > 1 else 1e12
ND = NCM * A_CELL
CASES = [("kpath_T_VSe15.npz", "V$_S$ vacancy"), ("kpath_T_OSe15.npz", "O$_S$"),
         ("kpath_T_SESe15.npz", "Se$_S$")]
WIN = [(-0.32, 0.32, "valence edge"), (1.08, 1.88, "conduction edge")]
print("n_d = %.1e cm^-2  ->  %.3e per cell  (1 per %.0f cells)" % (NCM, ND, 1/ND))

fig, ax = plt.subplots(2, 3, figsize=(13.2, 7.4), sharex="col")
dat = []
for fn, lab in CASES:
    z = np.load(f"{B}/{fn}")
    E, T, Hk, ekp, klab = z["E"], z["T"], z["Hk"], z["ekp"], z["klab"]
    eta = float(z["eta"]); nw = Hk.shape[1]; nkp = Hk.shape[0]; I = np.eye(nw)
    A = np.zeros((len(E), nkp))
    for i, e in enumerate(E):
        M = ((e + GATE_VBM)/RY + 1j*eta)*I[None] - Hk - ND*T[i]
        A[i] = -np.trace(np.linalg.inv(M), axis1=1, axis2=2).imag/np.pi
    dat.append((E, A, ekp, klab, nkp, lab))
vmax = max(a[1].max() for a in dat)
norm = LogNorm(vmin=vmax*3e-4, vmax=vmax)
for r, (lo, hi, wl) in enumerate(WIN):
    for c, (E, A, ekp, klab, nkp, lab) in enumerate(dat):
        x = np.arange(nkp)
        im = ax[r, c].pcolormesh(x, E, np.maximum(A, vmax*3e-4), cmap="magma",
                                 norm=norm, shading="gouraud", rasterized=True)
        for b in range(ekp.shape[1]):
            ax[r, c].plot(x, ekp[:, b], color="#5fb3ff", lw=.7, alpha=.6)
        for xl in klab[1:-1]: ax[r, c].axvline(xl, color="w", lw=.6, alpha=.4)
        ax[r, c].set_xticks(klab); ax[r, c].set_xticklabels(["$\\Gamma$", "M", "K", "$\\Gamma$"])
        ax[r, c].set_xlim(0, nkp-1); ax[r, c].set_ylim(lo, hi); ax[r, c].tick_params(labelsize=9)
        if r == 0: ax[r, c].set_title(lab, fontsize=11, loc="left")
        if c == 0: ax[r, c].set_ylabel("$E-E_\\mathrm{VBM}$   (eV)")
        ax[r, c].text(.02, .04, wl, transform=ax[r, c].transAxes, fontsize=8.5,
                      color="w", alpha=.85)
cb = fig.colorbar(im, ax=ax, pad=.012, fraction=.022)
cb.set_label("$A(k,\\omega)$   (states / eV, log)", fontsize=9); cb.ax.tick_params(labelsize=8)
eta_mev = float(np.load(f"{B}/{CASES[0][0]}")["eta"])*RY*1e3
nf = int(np.load(f"{B}/{CASES[0][0]}")["nfine"])
ttl = ("$n_d = %.0f\\times10^{%d}$ cm$^{-2}$   ($%.2f\\times10^{-3}$ per cell, one per %.0f cells)"
       "     $\\eta = %.0f$ meV,  %d$\\times$%d mesh"
       % (NCM/10**int(np.log10(NCM)), int(np.log10(NCM)), ND*1e3, 1/ND, eta_mev, nf, nf))
fig.suptitle(ttl, fontsize=10.5, x=.42, y=.965)
out = "/home/rjguo/edi_tmatrix/sternA/claude-sternheimer/docs/assets/zoom_edges_1e12.png"
fig.savefig(out, dpi=165, bbox_inches="tight"); print("wrote", out)
for (E, A, ekp, klab, nkp, lab) in dat:
    g = (E > 0.02) & (E < 1.64)
    print("  %-16s peak A in gap %8.3f   at E=%+.4f eV   band-edge A %.1f"
          % (lab, A[g].max(), E[g][A[g].max(axis=1).argmax()], A[np.argmin(np.abs(E))].max()))
