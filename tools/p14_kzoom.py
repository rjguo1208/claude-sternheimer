#!/usr/bin/env python3
"""Conduction-edge zoom about K, in the ranges of the earlier Banff T-matrix
figure: k from K-0.3*GK to K+0.45*KM, E-E_VBM in [1.45,2.15] eV, eta = 5 meV.

Top: A(k,omega) for the three defects, log scale, bare bands overlaid.
Bottom left: line cuts at K.  Bottom right: the flat-band detector, median_k A,
which picks out k-independent (localized) weight and suppresses the dispersive
host band."""
import os, sys, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
RY = 13.605693122994; GATE_VBM = -5.9359
A_CELL = np.sqrt(3)/2 * 3.18517668**2 * 1e-16
B = "/home/rjguo/edi_tmatrix/sternA/vsrlx"
ND = float(sys.argv[1]) if len(sys.argv) > 1 else 1/36
SUF = sys.argv[2] if len(sys.argv) > 2 else "kz"          # kz = conduction edge, vz = valence edge
CASES = [("kpath_T_OS%s.npz" % SUF, "O$_S$", "#c1121f"), ("kpath_T_SES%s.npz" % SUF, "Se$_S$", "#1f3b73"),
         ("kpath_T_VS%s.npz" % SUF, "V$_S$", "#2a7a2a")]
# supercell reference levels for V_S at n_d = 1/36, gate frame (section 4)
TRUTH = [-0.2981, -0.2895, -0.2558, -0.0899, -0.0248, +0.0596]
B_ = "/home/rjguo/edi_tmatrix/sternA/vsrlx"
CASES = [c for c in CASES if os.path.exists(f"{B_}/{c[0]}")]        # plot what exists
print("panels:", ", ".join(c[1].replace("$_S$", "_S") for c in CASES))

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

fig = plt.figure(figsize=(4.6*len(CASES)+0.5, 8.6))
NP = len(CASES)
gs = fig.add_gridspec(2, NP, height_ratios=[1.45, 1], hspace=.30, wspace=.13)
vmax = max(d[1].max() for d in dat)
norm = LogNorm(vmin=vmax*2e-4, vmax=vmax)
for c, (E, A, ekp, klab, nkp, lab, col, etam, nf) in enumerate(dat):
    ax = fig.add_subplot(gs[0, c]); x = np.arange(nkp)
    im = ax.pcolormesh(x, E, np.maximum(A, vmax*2e-4), cmap="magma", norm=norm,
                       shading="gouraud", rasterized=True)
    for b in range(ekp.shape[1]):
        ax.plot(x, ekp[:, b], color="#5fb3ff", lw=.55, alpha=.55)
    ax.axvline(klab[1], color="w", lw=.7, alpha=.55)
    if SUF == "vz" and lab.startswith("V$_S$"):
        for tv in TRUTH:
            if E.min() < tv < E.max():
                ax.axhline(tv, color="#7CFC00", lw=.8, ls=":", alpha=.85)
    ax.set_xticks(klab); ax.set_xticklabels(["$-0.3$", "K", "$+0.45$"], fontsize=9)
    ax.set_xlabel("$K + f\\,\\overline{\\Gamma K}$   |   $K + f\\,\\overline{KM}$", fontsize=8.5)
    ax.set_xlim(0, nkp-1); ax.set_ylim(E.min(), E.max()); ax.tick_params(labelsize=9)
    ax.set_title(lab, fontsize=11)
    if c == 0: ax.set_ylabel("$E-E_\\mathrm{VBM}$   (eV)")
    else: ax.set_yticklabels([])
cb = fig.colorbar(im, ax=fig.axes[:NP], pad=.012, fraction=.022)
cb.set_label("$A(k,\\omega)$   (states/eV, log)", fontsize=9); cb.ax.tick_params(labelsize=8)

axL = fig.add_subplot(gs[1, 0:max(1, NP-1)])
for (E, A, ekp, klab, nkp, lab, col, etam, nf) in dat:
    axL.semilogy(E, A[:, klab[1]], color=col, lw=1.6, label="%s @ K" % lab)
axL.set_xlim(E.min(), E.max()); axL.set_xlabel("$E-E_\\mathrm{VBM}$   (eV)")
axL.set_ylabel("$A$   (states/eV)"); axL.legend(frameon=False, fontsize=9)
axL.tick_params(labelsize=9); axL.grid(alpha=.18, lw=.5)
axL.set_title("line cuts at K", fontsize=9.5, loc="left")
fig.suptitle("$\\eta=%.0f$ meV,  $\\Delta\\omega=%.1f$ meV,  %d path points,  $n_d=%s$"
             % (dat[0][7], 1e3*(dat[0][0][1]-dat[0][0][0]), dat[0][4],
                "1/36" if abs(ND-1/36) < 1e-9 else "%.2e" % ND), fontsize=10, y=.965)

axR = fig.add_subplot(gs[1, NP-1])
for (E, A, ekp, klab, nkp, lab, col, etam, nf) in dat:
    axR.plot(E, np.median(A, axis=1), color=col, lw=1.6, label=lab)
axR.set_xlim(*( (1.75, E.max()) if SUF == "kz" else (E.min(), E.max()) )); axR.set_xlabel("$E-E_\\mathrm{VBM}$   (eV)")
axR.set_ylabel("median$_k\\,A(\\omega)$"); axR.legend(frameon=False, fontsize=9)
axR.tick_params(labelsize=9); axR.grid(alpha=.18, lw=.5)
axR.set_title("flat-band detector", fontsize=9.5, loc="left")
out = "/home/rjguo/edi_tmatrix/sternA/claude-sternheimer/docs/assets/kzoom_%s.png" % ("cbm" if SUF == "kz" else "vbm")
fig.savefig(out, dpi=160, bbox_inches="tight"); print("wrote", out)
for (E, A, ekp, klab, nkp, lab, col, etam, nf) in dat:
    m = np.median(A, axis=1)
    i = np.where((m[1:-1] > m[:-2]) & (m[1:-1] >= m[2:]) & (m[1:-1] > 2*np.median(m)))[0]+1
    print("  %-6s eta %.0f meV, Nf %d;  flat features: %s"
          % (lab.replace("$_S$", "_S"), etam, nf,
             " ".join("%+.4f" % E[j] for j in i) or "(none)"))
