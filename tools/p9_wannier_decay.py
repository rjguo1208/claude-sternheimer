#!/usr/bin/env python3
"""Does the defect vertex survive Wannier interpolation?

Both indices of the pair kernel are ELECTRON Wannier positions,
    M(R_e,R_p) = <w_{R_e}| dV |w_{R_p}> ,
exactly Lu-Bernardi Eq. (5) and the same convention as this project's earlier
tmatrix_p6_wannier.py.  The defect is not a third index: it is where dV sits, so
it is the point the kernel decays about.  In EDI's own supercells the defect is
at the origin and never appears; ours sits at the cell centre, R_def = (3,3).

(a) every one of the 2.98M matrix elements against the radius at which a
    spherical truncation would drop it.
(b) the diagonal element of one Wannier orbital, M_mm(R,R) -- the S p orbitals
    on the vacancy sublattice see the full potential, the far S barely anything.
(c) leave-one-out: build the kernel from the Gamma coset (6.2% of pairs),
    predict the rest."""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

RY = 13.605693122994; NK = 144
z = np.load("/home/rjguo/edi_tmatrix/sternA/vsrlx/loo.npz")
C12, C66, CG = "#1f3b73", "#c1121f", "#8a8f98"
CMO, CSV, CSF = "#7a5c00", "#c1121f", "#2a7a2a"
fig, ax = plt.subplots(1, 3, figsize=(14.6, 4.3))
rng = np.random.default_rng(0)

# ---- (a) every matrix element ------------------------------------------------
A, d = z["absM_1212"], z["draw_1212"]
rad = np.maximum(d[:, None], d[None, :])                       # truncation radius
flat_r = np.repeat(rad.ravel(), A.shape[2]*A.shape[3])
flat_v = A.reshape(-1)
keep = rng.choice(flat_v.size, 120000, replace=False)
ax[0].semilogy(flat_r[keep], np.maximum(flat_v[keep], 1e-12), ".", color=C12,
               ms=1.0, alpha=.10, rasterized=True)
ed = np.arange(-0.5, 7.6, 1.0)
xs, env, med = [], [], []
for lo, hi in zip(ed[:-1], ed[1:]):
    m = (flat_r >= lo) & (flat_r < hi)
    if m.any():
        xs.append(max(0.0, 0.5*(lo+hi))); env.append(flat_v[m].max())
        med.append(np.median(flat_v[m][flat_v[m] > 0]))
ax[0].semilogy(xs, env, "-o", color="k", lw=1.8, ms=5, label="largest element in shell")
ax[0].semilogy(xs, med, "--s", color=CG, lw=1.5, ms=4.5, label="median element")
ax[0].axvline(3.0, color=C66, lw=.9, ls=":")
ax[0].text(3.05, 2e-9, " 6$\\times$6 cell edge", color=C66, fontsize=8, va="bottom")
ax[0].set_xlim(-0.4, 7.2); ax[0].set_ylim(1e-9, 12)
ax[0].set_xlabel("$\\max(|R_e-R_\\mathrm{def}|,\\,|R_p-R_\\mathrm{def}|)$   (lattice constants)")
ax[0].set_ylabel("$|\\mathcal{M}_{mn}(R_e,R_p)|$")
ax[0].legend(frameon=False, fontsize=8.5, loc="upper right")
ax[0].set_title("(a)  all 2.98M matrix elements", fontsize=10, loc="left")

# ---- (b) one orbital's diagonal ---------------------------------------------
o = np.argsort(d)
grp = ((range(0, 5), CMO, "Mo $d$ (5)"), (range(5, 8), CSV, "S $p$, vacancy sublattice (3)"),
       (range(8, 11), CSF, "S $p$, far S (3)"))
edb = np.arange(-0.5, 7.6, 1.0)
for orbs, c, lab in grp:
    v = np.concatenate([A[np.arange(len(d)), np.arange(len(d)), m, m] for m in orbs])
    x = np.tile(d, len(list(orbs)))
    ax[1].semilogy(x, np.maximum(v, 1e-12), "o", color=c, ms=2.6, alpha=.30, mew=0)
    xs2, ev2 = [], []
    for lo, hi in zip(edb[:-1], edb[1:]):
        m = (x >= lo) & (x < hi)
        if m.any(): xs2.append(max(0.0, 0.5*(lo+hi))); ev2.append(v[m].max())
    ax[1].semilogy(xs2, ev2, "-", color=c, lw=1.9, label=lab)
ax[1].axvline(3.0, color=C66, lw=.9, ls=":")
ax[1].set_xlim(-0.3, 7.2); ax[1].set_ylim(1e-8, 12)
ax[1].set_xlabel("$|R-R_\\mathrm{def}|$   (lattice constants)")
ax[1].set_ylabel("$|\\mathcal{M}_{mm}(R,R)|$")
ax[1].legend(frameon=False, fontsize=8.5, loc="upper right")
ax[1].set_title("(b)  diagonal element, per Wannier orbital", fontsize=10, loc="left")

# ---- (c) leave-one-out -------------------------------------------------------
err, on = z["err"], z["on"]
mev = 1e3*RY*err/NK
for lab, m, c, ls in (("one off-grid", np.outer(on, ~on) | np.outer(~on, on), C66, "--"),
                      ("both off-grid (93.8%)", np.outer(~on, ~on), C12, "-")):
    v = np.sort(mev[m].ravel())
    ax[2].semilogx(np.maximum(v, 1e-9), np.linspace(0, 100, len(v)), ls, color=c, lw=2.0, label=lab)
for x, lab, c in ((7.6, "downfold residual", CMO), (11.0, "6$\\times$6 image error", CSF)):
    ax[2].axvline(x, color=c, lw=1.1, ls="-.", alpha=.85)
    ax[2].text(x*1.09, 30, lab, color=c, fontsize=8.5, rotation=90, va="center")
ax[2].text(0.135, 88, "training pairs return\nat $4\\times10^{-15}$ (off scale)",
           fontsize=8.5, color="#666", va="top")
ax[2].set_xlim(0.12, 40); ax[2].set_ylim(0, 100)
ax[2].set_xlabel("interpolation error of one element   (meV in $H_\\mathrm{eff}$)")
ax[2].set_ylabel("percentile of $(k,k')$ pairs")
ax[2].legend(frameon=False, fontsize=8.5, loc="center left", bbox_to_anchor=(0.02, 0.45))
ax[2].set_title("(c)  leave-one-out prediction", fontsize=10, loc="left")
for a in ax: a.tick_params(labelsize=9); a.grid(alpha=.18, lw=.6)
fig.tight_layout()
out = "/home/rjguo/edi_tmatrix/sternA/claude-sternheimer/docs/assets/wannier_decay.png"
fig.savefig(out, dpi=165); print("wrote", out)
print("  on-site diagonal by orbital:", " ".join("%.3f" % A[o[0], o[0], m, m] for m in range(11)))
for orbs, c, lab in grp:
    m0 = list(orbs)[0]
    v = A[o, o, m0, m0]
    print("  %-32s on site %.3f  ->  at 3a %.2e  ->  at 6a %.2e"
          % (lab, v[0], v[np.argmin(np.abs(d[o]-3.0))], v[np.argmin(np.abs(d[o]-6.0))]))
