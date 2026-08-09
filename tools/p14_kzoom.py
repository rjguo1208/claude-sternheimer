#!/usr/bin/env python3
"""Conduction-edge zoom about K, in the ranges of the earlier Banff T-matrix
figure: k from K-0.3*GK to K+0.45*KM, E-E_VBM in [1.45,2.15] eV, eta = 5 meV.

Top: A(k,omega) for the three defects, log scale, bare bands overlaid.
Bottom left: line cuts at K.  Bottom right: the flat-band detector, median_k A,
which picks out k-independent (localized) weight and suppresses the dispersive
host band."""
import sys, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
RY = 13.605693122994; GATE_VBM = -5.9359
A_CELL = np.sqrt(3)/2 * 3.18517668**2 * 1e-16
B = "/home/rjguo/edi_tmatrix/sternA/vsrlx"
ND = float(sys.argv[1]) if len(sys.argv) > 1 else 1/36
CASES = [("kpath_T_OSkz.npz", "O$_S$", "#c1121f"), ("kpath_T_SESkz.npz", "Se$_S$", "#1f3b73"),
         ("kpath_T_VSkz.npz", "V$_S$", "#2a7a2a")]

dat = []
for fn, lab, col in CASES:
    z = np.load(f"{B}/{fn}")
    E, T, Hk, ekp, klab = z["E"], z["T"], z["Hk"], z["ekp"], z["klab"]
    eta = float(z["eta"]); nw = Hk.shape[1]; nkp = Hk.shape[0]; I = np.eye(nw)
    A = np.zeros((len(E), nkp))
    for i, e in enumerate(E):
        M = ((e + GATE_VBM)/RY + 1j*eta)*I[None] - Hk - ND*T[i]
        A[i] = -np.trace(np.linalg.inv(M), axis1=1, axis2=2).imag/np.pi
    dat.append((E, A, ekp, klab, nkp, lab, col, eta*RY*1e3, int(z["nfine"])))

fig = plt.figure(figsize=(13.6, 8.6))
gs = fig.add_gridspec(2, 3, height_ratios=[1.45, 1], hspace=.30, wspace=.13)
vmax = max(d[1].max() for d in dat)
norm = LogNorm(vmin=vmax*2e-4, vmax=vmax)
for c, (E, A, ekp, klab, nkp, lab, col, etam, nf) in enumerate(dat):
    ax = fig.add_subplot(gs[0, c]); x = np.arange(nkp)
    im = ax.pcolormesh(x, E, np.maximum(A, vmax*2e-4), cmap="magma", norm=norm,
                       shading="gouraud", rasterized=True)
    for b in range(ekp.shape[1]):
        ax.plot(x, ekp[:, b], color="#5fb3ff", lw=.55, alpha=.55)
    ax.axvline(klab[1], color="w", lw=.7, alpha=.55)
    ax.set_xticks(klab); ax.set_xticklabels(["$-0.3$", "K", "$+0.45$"], fontsize=9)
    ax.set_xlabel("$K + f\\,\\overline{\\Gamma K}$   |   $K + f\\,\\overline{KM}$", fontsize=8.5)
    ax.set_xlim(0, nkp-1); ax.set_ylim(E.min(), E.max()); ax.tick_params(labelsize=9)
    ax.set_title(lab, fontsize=11)
    if c == 0: ax.set_ylabel("$E-E_\\mathrm{VBM}$   (eV)")
    else: ax.set_yticklabels([])
cb = fig.colorbar(im, ax=fig.axes[:3], pad=.012, fraction=.022)
cb.set_label("$A(k,\\omega)$   (states/eV, log)", fontsize=9); cb.ax.tick_params(labelsize=8)

axL = fig.add_subplot(gs[1, 0:2])
for (E, A, ekp, klab, nkp, lab, col, etam, nf) in dat:
    axL.semilogy(E, A[:, klab[1]], color=col, lw=1.6, label="%s @ K" % lab)
axL.set_xlim(E.min(), E.max()); axL.set_xlabel("$E-E_\\mathrm{VBM}$   (eV)")
axL.set_ylabel("$A$   (states/eV)"); axL.legend(frameon=False, fontsize=9)
axL.tick_params(labelsize=9); axL.grid(alpha=.18, lw=.5)
axL.set_title("line cuts at K   ($\\eta=%.0f$ meV, $\\Delta\\omega=%.1f$ meV, %d $k$, $n_d=%s$)"
              % (dat[0][7], 1e3*(dat[0][0][1]-dat[0][0][0]), dat[0][4],
                 "1/36" if abs(ND-1/36) < 1e-9 else "%.2e" % ND), fontsize=9.5, loc="left")

axR = fig.add_subplot(gs[1, 2])
for (E, A, ekp, klab, nkp, lab, col, etam, nf) in dat:
    axR.plot(E, np.median(A, axis=1), color=col, lw=1.6, label=lab)
axR.set_xlim(1.75, E.max()); axR.set_xlabel("$E-E_\\mathrm{VBM}$   (eV)")
axR.set_ylabel("median$_k\\,A(\\omega)$"); axR.legend(frameon=False, fontsize=9)
axR.tick_params(labelsize=9); axR.grid(alpha=.18, lw=.5)
axR.set_title("flat-band detector", fontsize=9.5, loc="left")
out = "/home/rjguo/edi_tmatrix/sternA/claude-sternheimer/docs/assets/kzoom_cbm.png"
fig.savefig(out, dpi=160, bbox_inches="tight"); print("wrote", out)
for (E, A, ekp, klab, nkp, lab, col, etam, nf) in dat:
    m = np.median(A, axis=1); s = (E > 1.75)
    i = np.where((m[1:-1] > m[:-2]) & (m[1:-1] >= m[2:]) & (m[1:-1] > 2*np.median(m[s])))[0]+1
    print("  %-6s eta %.0f meV, Nf %d;  flat features above 1.75 eV: %s"
          % (lab.replace("$_S$", "_S"), etam, nf,
             " ".join("%+.4f" % E[j] for j in i if E[j] > 1.75) or "(none)"))
