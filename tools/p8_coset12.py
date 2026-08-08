#!/usr/bin/env python3
"""12x12 by coset decomposition: density of states and defect levels.

DeltaV of a 6x6-periodic defect has Fourier weight only on supercell G-vectors, so a
12x12 primitive k-grid splits into 4 mutually uncoupled 36-k cosets, each of which is
one shifted 6x6 run.  (a) DOS of the downfolded active manifold, 6x6 (Gamma coset)
vs 12x12 (4 cosets), against the supercell truth on the same 4 supercell k-points.
(b) the omega-resolved defect levels at the supercell Gamma and M points.  Both are
plotted with the single calibrated vacuum-alignment constant removed from the
downfold; nothing else is fitted."""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

BASE = "/home/rjguo/edi_tmatrix/sternA/vsrlx"
GATE_VBM = -5.9359   # zero of the published gate (defect supercell VBM), eV
REF = 0.0596         # sc_eigs_4k zero -> gate zero
PRED = 0.0179        # parameter-free vacuum-alignment prediction for V_S, eV
ETA = 0.030
C_DF, C_TR, C_66 = "#1f3b73", "#8a8f98", "#c1121f"

d = np.load(f"{BASE}/statdos_VS_4coset.npz"); s = np.load(f"{BASE}/sc_eigs_4k.npz")
vb = float(s["vbm"])

# --- omega-resolved defect levels (gate_omega QP fixed points) ----------------
QP = {"$\\Gamma$": [-0.0652, -0.0010, +0.0840, +1.1985, +1.7003],
      "M":         [-0.0758, -0.0226, +0.0703, +1.2002, +1.2073]}
TR = {"$\\Gamma$": sorted(set(np.round(s["gamma"] - vb + REF, 4))),
      "M":         sorted(set(np.round(s["k0"] - vb + REF, 4)))}
for k in TR:
    TR[k] = [x for x in TR[k] if -0.12 < x < 1.80]

print("omega-resolved downfold minus supercell truth (meV), order-preserving match")
res = {}
for k in QP:
    assert len(QP[k]) == len(TR[k]), (k, len(QP[k]), len(TR[k]))
    dd = np.array([1000*(a - b) for a, b in zip(sorted(QP[k]), sorted(TR[k]))])
    res[k] = dd
    print("  %-9s %s   mean %+.1f  rms %.1f   -> after reference %+.1f +- %.1f meV"
          % (k.replace("$\\Gamma$", "Gamma"), np.round(dd, 1), dd.mean(), dd.std(),
             dd.mean() - 1000*PRED, dd.std()))

# --- figure ------------------------------------------------------------------
fig, ax = plt.subplots(1, 2, figsize=(11.2, 4.2), gridspec_kw=dict(width_ratios=[1.5, 1]))

w = np.arange(-0.60, 2.30, 0.004)
def dos(ev, shift):
    g = np.zeros_like(w)
    for e in np.asarray(ev) - shift:
        if w[0]-0.3 < e < w[-1]+0.3:
            g += np.exp(-0.5*((w - e)/ETA)**2)
    return g / (ETA*np.sqrt(2*np.pi))

dos66 = dos(d["e0"], GATE_VBM + PRED)
dos12 = np.mean([dos(d["e%d" % i], GATE_VBM + PRED) for i in range(4)], axis=0)
tr = np.mean([dos(s[k], vb - REF) for k in ("gamma", "k0", "k1", "k2")], axis=0)

ax[0].fill_between(w, tr, color=C_TR, alpha=.40, lw=0, label="supercell truth (4 $k_{sc}$)")
ax[0].plot(w, dos66, color=C_66, lw=1.5, ls="--", label="downfold 6$\\times$6 ($\\Gamma$ coset)")
ax[0].plot(w, dos12, color=C_DF, lw=2.0, label="downfold 12$\\times$12 (4 cosets)")
ax[0].set_xlim(-0.60, 2.15)
m = (w > -0.60) & (w < 2.15)
ax[0].set_ylim(0, 1.14*max(tr[m].max(), dos12[m].max(), dos66[m].max()))
ax[0].set_xlabel("$E-E_{\\rm VBM}$ (eV)"); ax[0].set_ylabel("states / eV / coset")
ax[0].legend(frameon=False, fontsize=9, loc="upper left")
ax[0].set_title("(a)  DOS of the downfolded active manifold", fontsize=10, loc="left")
ax[0].annotate("defect levels", xy=(0.02, 0.20*ax[0].get_ylim()[1]),
               xytext=(0.45, 0.55*ax[0].get_ylim()[1]), fontsize=9, color="#444",
               arrowprops=dict(arrowstyle="->", color="#444", lw=.9))

for i, k in enumerate(["$\\Gamma$", "M"]):
    x0 = i*1.0
    for b in TR[k]:
        ax[1].plot([x0-0.34, x0-0.06], [b, b], color=C_TR, lw=3.0, solid_capstyle="butt")
    for a in QP[k]:
        ax[1].plot([x0+0.06, x0+0.34], [a-PRED, a-PRED], color=C_DF, lw=3.0, solid_capstyle="butt")
    ax[1].text(x0, 1.99, k, ha="center", fontsize=11)
    ax[1].text(x0, 1.86, "%+.1f $\\pm$ %.1f meV" % (res[k].mean()-1000*PRED, res[k].std()),
               ha="center", fontsize=8.5, color="#444")
ax[1].plot([], [], color=C_TR, lw=3.0, label="supercell truth")
ax[1].plot([], [], color=C_DF, lw=3.0, label="downfold, reference removed")
ax[1].set_xlim(-0.6, 1.6); ax[1].set_ylim(-0.16, 2.14)
ax[1].set_xticks([]); ax[1].set_ylabel("$E-E_{\\rm VBM}$ (eV)")
ax[1].legend(frameon=False, fontsize=9, loc="center right")
ax[1].set_title("(b)  $\\omega$-resolved defect levels", fontsize=10, loc="left")
for a in ax: a.tick_params(labelsize=9)
fig.tight_layout()
out = "/home/rjguo/edi_tmatrix/sternA/claude-sternheimer/docs/assets/coset12_VS.png"
fig.savefig(out, dpi=170); print("wrote", out)
